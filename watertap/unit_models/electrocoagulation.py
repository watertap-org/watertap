#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

from copy import deepcopy

# Import Pyomo libraries
from pyomo.environ import (
    Constraint,
    Var,
    check_optimal_termination,
    Param,
    Suffix,
    log,
    log10,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.initialization import ModularInitializerBase
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.constants import Constants
from idaes.core.util.exceptions import ConfigurationError, InitializationError
from idaes.core.util.misc import StrEnum
from idaes.core.scaling import CustomScalerBase, ConstraintScalingScheme
from idaes.core.util.tables import create_stream_table_dataframe
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.core import InitializationMixin
from watertap.core.solvers import get_solver
from watertap.costing.unit_models.electrocoagulation import cost_electrocoagulation

__author__ = "Kurban Sitterley, Mukta Hardikar"


class ElectrodeMaterial(StrEnum):
    aluminum = "aluminum"
    iron = "iron"


class ReactorMaterial(StrEnum):
    pvc = "pvc"
    stainless_steel = "stainless_steel"


class OverpotentialCalculation(StrEnum):
    nernst = "nernst"
    fixed = "fixed"
    regression = "regression"


"""
REFERENCES:

K. L. Dubrawski, C. Du and M. Mohseni
Electrochimica Acta 2014 Vol. 129 Pages 187-195
DOI: 10.1016/j.electacta.2014.02.089

Z. Gu, Z. Liao, M. Schulz, J. R. Davis, J. C. Baygents and J. Farrell
Industrial & Engineering Chemistry Research 2009 Vol. 48 Issue 6 Pages 3112-3117
DOI: 10.1021/ie801086c

"""


class ElectrocoagulationScaler(CustomScalerBase):
    """
    Scaler class for Electrocoagulation models
    """

    DEFAULT_SCALING_FACTORS = {
        "electrode_thickness": 1e3,
        "electrode_gap": 10,
        "electrolysis_time": 0.1,
        "applied_current": 1,
        "current_efficiency": 1,
        "cell_volume": 0.1,
        "ohmic_resistance": 1e3,
        "charge_loading_rate": 1e-2,
        "current_density": 1e-2,
        "cell_voltage": 0.1,
        "overpotential": 0.1,
    }

    def variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        """
        Apply variable scaling to model.

        Args:
            model: model to be scaled
            overwrite: whether to overwrite existing scaling factors
            submodel_scalers: ComponentMap of Scaler objects to be applied to sub-models

        Returns:
            None

        """
        # Scale inlet state
        self.call_submodel_scaler_method(
            submodel=model.properties_in,
            method="variable_scaling_routine",
            submodel_scalers=submodel_scalers,
            overwrite=overwrite,
        )

        # Propagate scaling from inlet state
        self.propagate_state_scaling(
            target_state=model.properties_out,
            source_state=model.properties_in,
            overwrite=overwrite,
        )
        self.propagate_state_scaling(
            target_state=model.properties_waste,
            source_state=model.properties_in,
            overwrite=overwrite,
        )

        # Scale remaining states
        self.call_submodel_scaler_method(
            submodel=model.properties_retentate,
            method="variable_scaling_routine",
            submodel_scalers=submodel_scalers,
            overwrite=overwrite,
        )
        self.call_submodel_scaler_method(
            submodel=model.properties_permeate,
            method="variable_scaling_routine",
            submodel_scalers=submodel_scalers,
            overwrite=overwrite,
        )

        # Scale remaining variables
        # Variables with default scaling factors
        for n in self.DEFAULT_SCALING_FACTORS.keys():
            c = model.find_component(n)
            for v in c.values():
                self.scale_variable_by_default(v, overwrite=overwrite)

    def constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        """
        Apply constraint scaling to model.

        Args:
            model: model to be scaled
            overwrite: whether to overwrite existing scaling factors
            submodel_scalers: ComponentMap of Scaler objects to be applied to sub-models

        Returns:
            None

        """
        # Call scaling methods for sub-models
        self.call_submodel_scaler_method(
            submodel=model.properties_in,
            method="constraint_scaling_routine",
            submodel_scalers=submodel_scalers,
            overwrite=overwrite,
        )
        self.call_submodel_scaler_method(
            submodel=model.properties_out,
            method="constraint_scaling_routine",
            submodel_scalers=submodel_scalers,
            overwrite=overwrite,
        )
        self.call_submodel_scaler_method(
            submodel=model.properties_waste,
            method="constraint_scaling_routine",
            submodel_scalers=submodel_scalers,
            overwrite=overwrite,
        )

        # Scale remaining constraints
        for c in model.component_data_objects(Constraint, descend_into=False):
            self.scale_constraint_by_nominal_value(
                c,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )


@declare_process_block_class("Electrocoagulation")
class ElectrocoagulationData(InitializationMixin, UnitModelBlockData):
    """
    Electrocoagulation model
    """

    CONFIG = ConfigBlock()

    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
    **default** = False.""",
        ),
    )

    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Indicates whether holdup terms should be constructed or not.
    **default** - False.""",
        ),
    )

    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property calculations,
    **default** - useDefault.
    **Valid values:** {
    **useDefault** - use default package from parent model or flowsheet,
    **PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )

    CONFIG.declare(
        "property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property packages",
            doc="""A ConfigBlock with arguments to be passed to a property block(s)
    and used when constructing these,
    **default** - None.
    **Valid values:** {
    see property package for documentation.}""",
        ),
    )

    CONFIG.declare(
        "electrode_material",
        ConfigValue(
            default="aluminum",
            domain=In(ElectrodeMaterial),
            description="Electrode material",
        ),
    )

    CONFIG.declare(
        "reactor_material",
        ConfigValue(
            default="pvc",
            domain=In(ReactorMaterial),
            description="Reactor material",
        ),
    )

    CONFIG.declare(
        "overpotential_calculation",
        ConfigValue(
            default="fixed",
            domain=In(OverpotentialCalculation),
            description="Determination of overpotential",
        ),
    )

    def build(self):

        super().build()

        solutes = self.config.property_package.solute_set

        if "TDS" not in solutes:
            raise ConfigurationError(
                "TDS must be in feed stream for solution conductivity estimation."
            )

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = True  # inlet block is an inlet
        self.properties_in = self.config.property_package.state_block_class(
            self.flowsheet().config.time, doc="Material properties of inlet", **tmp_dict
        )

        # Add outlet and waste block
        tmp_dict["defined_state"] = False
        self.properties_out = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of outlet",
            **tmp_dict,
        )

        self.properties_waste = self.config.property_package.state_block_class(
            self.flowsheet().config.time, doc="Material properties of waste", **tmp_dict
        )

        prop_in = self.properties_in[0]
        prop_out = self.properties_out[0]
        prop_waste = self.properties_waste[0]

        self.add_port(name="inlet", block=self.properties_in)
        self.add_port(name="outlet", block=self.properties_out)
        self.add_port(name="waste", block=self.properties_waste)

        removal_frac_dict = dict(
            zip(
                solutes,
                [0.7 for _ in solutes],
            )
        )

        self.removal_frac_mass_comp = Param(
            solutes,
            initialize=removal_frac_dict,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Component removal efficiency on mass basis",
        )

        self.recovery_frac_mass_water = Param(
            initialize=0.99,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Water recovery on mass basis",
        )

        self.tds_to_cond_conversion = Param(
            initialize=5e3,
            mutable=True,
            units=(pyunits.mg * pyunits.m) / (pyunits.liter * pyunits.S),
            doc="Conersion factor for mg/L TDS to S/m",
        )

        self.standard_temperature = Param(
            initialize=298.15,  # 25C
            mutable=True,
            units=pyunits.degK,
            doc="Standard temperature",
        )

        self.ec_ion_mw = Param(
            initialize=100e-3,
            units=pyunits.kg / pyunits.mol,
            doc="Molecular weight of aluminum",
        )

        self.ec_ion_z = Param(
            initialize=1,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Charge of reactive species",
        )

        self.stoic_coeff = Param(
            initialize=1,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Stoichiometric coefficient for electrode material",
        )

        self.density_electrode_material = Param(
            initialize=2000,
            units=pyunits.kg / pyunits.m**3,
            doc="Density of electrode material",
        )

        if self.config.overpotential_calculation == OverpotentialCalculation.nernst:

            self.anode_cell_potential_std = Param(
                initialize=-0.5,
                mutable=True,
                units=pyunits.volt,
                doc="Anodic non-equilibrium cell potential, standard @ 25C",
            )

            self.anode_entropy_change_std = Param(
                initialize=1e-4,
                mutable=True,
                units=pyunits.volt / pyunits.K,
                doc="Entropy change",
            )

            self.anodic_exchange_current_density = Param(
                initialize=2e-5,
                mutable=True,
                units=pyunits.ampere / pyunits.m**2,
                doc="Anodic exchange current density",
            )

            self.cathodic_exchange_current_density = Param(
                initialize=1.0e-4,
                mutable=True,
                units=pyunits.ampere / pyunits.m**2,
                doc="Cathodic exchange current density",
            )

            self.cathode_surface_pH = Param(
                initialize=11,
                mutable=True,
                units=pyunits.dimensionless,
                doc="pH at cathode surface",
            )

            self.cathode_cell_potential_std = Param(
                initialize=-0.83,
                mutable=True,
                units=pyunits.volt,
                doc="Cathode equilibrium cell potential, standard H2 @ 25C",
            )

            self.cathode_entropy_change_std = Param(
                initialize=-0.000836,  # = S / (Z * F)
                mutable=True,
                units=pyunits.volt / pyunits.K,
                doc="Entropy change for H2",
            )

            self.partial_pressure_H2 = Param(
                initialize=1,
                mutable=True,
                units=pyunits.atm,
                doc="Partial pressure of hydrogen gas",
            )

            self.frac_increase_temperature = Param(
                initialize=1.05,
                mutable=True,
                units=pyunits.dimensionless,
                doc="Fractional increase in water temperature from inlet to outlet",
            )

            self.tafel_slope_cathode = Var(
                initialize=0.146,
                bounds=(0, None),
                units=pyunits.volt,  # volt per 10x change (decade)
                doc="Tafel slope for cathode",
            )

            self.tafel_slope_anode = Var(
                initialize=0.093,
                bounds=(0, None),
                units=pyunits.volt,  # volt per 10x change (decade)
                doc="Tafel slope for anode",
            )

        if self.config.electrode_material == ElectrodeMaterial.aluminum:

            """
            Reaction at an anode is:
            Al_3+ + 3e_- <--> Al(s);  Ea0 = -1.66 V
            Reaction in bulk solution is:
            Al_3+ + 3OH_- <--> Al(OH)3(s)
            """

            self.ec_ion_mw.set_value(26.98e-3)
            self.ec_ion_z.set_value(3)
            self.stoic_coeff.set_value(1)
            self.density_electrode_material.set_value(2710)

            if self.config.overpotential_calculation == OverpotentialCalculation.nernst:
                self.anode_cell_potential_std.set_value(
                    -1.66
                )  # Dubrawski et al., 2014; Zhang et al., 2020
                self.anode_entropy_change_std.set_value(
                    0.000533
                )  # = S / (Z * F); Bratsch, 1989; Dubrawski et al., 2014
                self.anodic_exchange_current_density.set_value(
                    2.602e-5
                )  # Electrochemical Thermodynamics and Kinetics, pg 276
                self.cathodic_exchange_current_density.set_value(
                    1.0e-4
                )  # Electrochemical Thermodynamics and Kinetics, pg 276

        if self.config.electrode_material == ElectrodeMaterial.iron:

            """
            Reaction at an anode is:
            Fe_2+ + 2e_- <--> Fe(s);  Ea0 = -0.41 V
            Reaction in bulk solution is:
            Fe_2+ + 2OH_- <--> Fe(OH)2(s)
            """

            self.ec_ion_mw.set_value(55.845e-3)
            self.ec_ion_z.set_value(2)
            self.stoic_coeff.set_value(1)
            self.density_electrode_material.set_value(7860)

            if self.config.overpotential_calculation == OverpotentialCalculation.nernst:
                self.anode_cell_potential_std.set_value(
                    -0.41
                )  # Dubrawski et al., 2014; Zhang et al., 2020
                self.anode_entropy_change_std.set_value(
                    7e-5
                )  # = S / (Z * F); Bratsch, 1989; Dubrawski et al., 2014
                self.anodic_exchange_current_density.set_value(
                    2.5e-4
                )  # Electrochemical Thermodynamics and Kinetics, pg 276
                self.cathodic_exchange_current_density.set_value(
                    1e-3
                )  # Electrochemical Thermodynamics and Kinetics, pg 276

        self.floc_basin_vol = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.m**3,
            doc="Reactor volume total (flotation + sedimentation)",
        )

        self.floc_retention_time = Var(
            initialize=30,
            bounds=(2, 200),
            units=pyunits.minute,
            doc="Electrolysis time",
        )

        self.coagulant_dose = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.g / pyunits.liter,
            doc="Metal dose to the feed in g/L",
        )

        self.electrode_thickness = Var(
            initialize=0.001,
            bounds=(0, 0.1),
            units=pyunits.m,
            doc="Electrode thickness",
        )

        self.electrode_mass = Var(
            initialize=10,
            bounds=(0, None),
            units=pyunits.kg,
            doc="Electrode mass",
        )

        self.electrode_volume = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.m**3,
            doc="Electrode volume",
        )

        self.cathode_area = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.m**2,
            doc="Area of cathode",
        )

        self.anode_area = Var(
            initialize=1,
            bounds=(1e-6, None),
            units=pyunits.m**2,
            doc="Area of anode",
        )

        self.electrode_gap = Var(
            initialize=0.005,
            bounds=(0.00001, 0.5),
            units=pyunits.m,
            doc="Electrode gap",
        )

        self.electrolysis_time = Var(
            initialize=30,
            bounds=(0.1, 200),
            units=pyunits.minute,
            doc="Electrolysis time",
        )

        self.cell_volume = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.m**3,
            doc="Volume of electrolytic cell",
        )

        self.current_density = Var(
            initialize=1,
            bounds=(0, 2000),
            units=pyunits.ampere / pyunits.m**2,
            doc="Current density",
        )

        self.applied_current = Var(
            initialize=1e4,
            bounds=(0, None),
            units=pyunits.ampere,
            doc="Applied current",
        )

        self.ohmic_resistance = Var(
            initialize=1e-5,
            bounds=(0, None),
            units=pyunits.ohm * pyunits.m**2,
            doc="Ohmic resistance of solution",
        )

        self.charge_loading_rate = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.coulomb / pyunits.liter,
            doc="Charge loading rate",
        )

        self.current_efficiency = Var(
            initialize=1,
            bounds=(0, 5),
            units=pyunits.kg / pyunits.kg,
            doc="Current efficiency",
        )

        self.cell_voltage = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.volt,
            doc="Cell voltage",
        )

        self.overpotential = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.volt,
            doc="Overpotential",
        )

        @self.Expression(doc="Conductivity")
        def conductivity(b):
            return pyunits.convert(
                prop_in.conc_mass_phase_comp["Liq", "TDS"] / b.tds_to_cond_conversion,
                to_units=pyunits.S / pyunits.m,
            )

        @self.Expression(doc="Theoretical metal loading")
        def theoretical_coagulant_dose(b):
            return pyunits.convert(
                (b.applied_current * b.ec_ion_mw)
                / (
                    Constants.faraday_constant
                    * b.ec_ion_z
                    * prop_in.flow_vol_phase["Liq"]
                ),
                to_units=pyunits.g / pyunits.liter,
            )

        @self.Expression(doc="Total electrode area")
        def electrode_area_total(b):
            return b.anode_area

        @self.Expression(doc="Power required")
        def power_required(b):
            return pyunits.convert(
                b.cell_voltage * b.applied_current, to_units=pyunits.kW
            )

        @self.Expression(doc="Power density Faradaic")
        def power_density_faradaic(b):
            return pyunits.convert(
                (b.overpotential * b.applied_current) / b.electrode_area_total,
                to_units=pyunits.microwatts / pyunits.cm**2,
            )

        @self.Expression(doc="Power density total")
        def power_density_total(b):
            return pyunits.convert(
                b.power_required / (b.electrode_area_total * 0.5),
                to_units=pyunits.microwatts / pyunits.cm**2,
            )

        if self.config.overpotential_calculation == OverpotentialCalculation.regression:

            self.overpotential_k1 = Var(
                units=pyunits.millivolt,
                doc="Constant k1 in overpotential equation from Gu (2009)",
            )

            self.overpotential_k2 = Var(
                units=pyunits.millivolt,
                doc="Constant k2 in overpotential equation from Gu (2009)",
            )

            @self.Constraint(doc="Overpotential calculation")
            def eq_overpotential(b):
                cd_dimensionless = pyunits.convert(
                    b.current_density * pyunits.cm**2 / pyunits.milliampere,
                    to_units=pyunits.dimensionless,
                )
                return b.overpotential == pyunits.convert(
                    (b.overpotential_k1 * log(cd_dimensionless) + b.overpotential_k2),
                    to_units=pyunits.volt,
                )

        if self.config.overpotential_calculation == OverpotentialCalculation.nernst:

            @self.Expression(doc="Hydroxide concentration at cathode surface")
            def cathode_conc_mol_hydroxide(b):
                pOH = 14 - b.cathode_surface_pH
                return 10 ** (-pOH) * pyunits.mol / pyunits.liter

            @self.Constraint(doc="Temperature change")
            def eq_temperature_change(b):
                return (
                    prop_out.temperature
                    == b.frac_increase_temperature * prop_in.temperature
                )

            @self.Expression(doc="Change in effluent temperature relative to standard")
            def temp_diff_std(b):
                return prop_out.temperature - b.standard_temperature

            @self.Expression(
                doc="Anode equilibrium potential adjusted for outlet temperature"
            )
            def anode_cell_potential_temp_adj(b):
                return pyunits.convert(
                    b.anode_cell_potential_std
                    + b.anode_entropy_change_std * b.temp_diff_std,
                    to_units=pyunits.volt,
                )

            @self.Expression(doc="Anode equilibrium cell potential")
            def anode_cell_potential(b):
                ec_dose_dimensionless = pyunits.convert(
                    (b.coagulant_dose / b.ec_ion_mw) * pyunits.liter * pyunits.mol**-1,
                    to_units=pyunits.dimensionless,
                )
                return b.anode_cell_potential_temp_adj - (
                    pyunits.convert(
                        (
                            (2.303 * Constants.gas_constant * prop_out.temperature)
                            / (Constants.faraday_constant * b.ec_ion_z)
                        ),
                        to_units=pyunits.volt,
                    )
                    * log10(ec_dose_dimensionless**-b.stoic_coeff)
                )

            @self.Expression(
                doc="Cathode equilibrium potential adjusted for outlet temperature"
            )
            def cathode_cell_potential_temp_adj(b):
                return pyunits.convert(
                    b.cathode_cell_potential_std
                    + b.cathode_entropy_change_std * b.temp_diff_std,
                    to_units=pyunits.volt,
                )

            @self.Expression(doc="Cathode equilibrium cell potential")
            def cathode_cell_potential(b):
                ccmh_dimensionless = pyunits.convert(
                    b.cathode_conc_mol_hydroxide * pyunits.liter / pyunits.mol,
                    to_units=pyunits.dimensionless,
                )
                ph2_dimensionless = pyunits.convert(
                    b.partial_pressure_H2 * pyunits.atm**-1,
                    to_units=pyunits.dimensionless,
                )
                return b.cathode_cell_potential_temp_adj - pyunits.convert(
                    (
                        (2.303 * Constants.gas_constant * prop_out.temperature)
                        / (Constants.faraday_constant * b.ec_ion_z)
                    ),
                    to_units=pyunits.volt,
                ) * log10(ccmh_dimensionless**2 * ph2_dimensionless)

            @self.Expression(doc="Anode overpotential")
            def anode_overpotential(b):
                return b.tafel_slope_anode * log10(
                    pyunits.convert(
                        b.current_density / b.anodic_exchange_current_density,
                        to_units=pyunits.dimensionless,
                    )
                )

            @self.Expression(doc="Cathode overpotential")
            def cathode_overpotential(b):
                return -1 * (
                    b.tafel_slope_cathode
                    * log10(
                        pyunits.convert(
                            b.current_density / b.cathodic_exchange_current_density,
                            to_units=pyunits.dimensionless,
                        )
                    )
                )

            @self.Constraint(doc="Overpotential calculation")
            def eq_nernst_overpotential(b):
                return b.overpotential == abs(
                    b.cathode_cell_potential - b.anode_cell_potential
                ) + b.anode_overpotential + abs(b.cathode_overpotential)

        ### MASS BALANCE ###

        @self.Constraint(doc="Water recovery")
        def eq_water_recovery(b):
            return (
                prop_out.flow_mass_phase_comp["Liq", "H2O"]
                == prop_in.flow_mass_phase_comp["Liq", "H2O"]
                * b.recovery_frac_mass_water
            )

        @self.Constraint(doc="Water mass balance")
        def eq_water_mass_balance(b):
            return (
                prop_in.flow_mass_phase_comp["Liq", "H2O"]
                == prop_out.flow_mass_phase_comp["Liq", "H2O"]
                + prop_waste.flow_mass_phase_comp["Liq", "H2O"]
            )

        @self.Constraint(solutes, doc="Component removal")
        def eq_component_removal(b, j):
            return (
                b.removal_frac_mass_comp[j] * prop_in.flow_mass_phase_comp["Liq", j]
                == prop_waste.flow_mass_phase_comp["Liq", j]
            )

        @self.Constraint(solutes, doc="Component mass balance")
        def eq_comp_mass_balance(b, j):
            return (
                prop_in.flow_mass_phase_comp["Liq", j]
                == prop_out.flow_mass_phase_comp["Liq", j]
                + prop_waste.flow_mass_phase_comp["Liq", j]
            )

        @self.Constraint(doc="Charge loading rate equation")
        def eq_charge_loading_rate(b):
            return (
                pyunits.convert(
                    b.charge_loading_rate * prop_in.flow_vol_phase["Liq"],
                    to_units=pyunits.ampere,
                )
                == b.applied_current
            )

        @self.Constraint(doc="Total flocculation tank volume")
        def eq_floc_reactor_volume(b):
            return b.floc_basin_vol == pyunits.convert(
                prop_in.flow_vol * b.floc_retention_time,
                to_units=pyunits.m**3,
            )

        @self.Constraint(doc="Faraday's Law")
        def eq_faraday(b):
            return (
                b.applied_current * b.current_efficiency * b.ec_ion_mw
                == pyunits.convert(
                    prop_in.flow_vol_phase["Liq"]
                    * b.coagulant_dose
                    * b.ec_ion_z
                    * Constants.faraday_constant,
                    to_units=(pyunits.kg * pyunits.amp) / pyunits.mol,
                )
            )

        @self.Constraint(doc="Total anode area required")
        def eq_anode_area(b):
            return b.anode_area == pyunits.convert(
                b.applied_current / b.current_density, to_units=pyunits.m**2
            )

        @self.Constraint(doc="Cell voltage")
        def eq_cell_voltage(b):
            return b.cell_voltage == pyunits.convert(
                b.overpotential
                + b.applied_current * (b.ohmic_resistance / b.anode_area),
                to_units=pyunits.volt,
            )

        @self.Constraint(doc="Electrode volume")
        def eq_electrode_volume(b):
            return (
                b.electrode_volume
                == (b.anode_area + b.cathode_area) * b.electrode_thickness
            )

        @self.Constraint(doc="Cathode and anode areas are equal")
        def eq_cathode_anode(b):
            return b.cathode_area == b.anode_area

        @self.Constraint(doc="Total reactor volume")
        def eq_reactor_volume(b):
            return b.cell_volume == pyunits.convert(
                prop_in.flow_vol_phase["Liq"] * b.electrolysis_time,
                to_units=pyunits.m**3,
            )

        @self.Constraint(doc="Ohmic resistance")
        def eq_ohmic_resistance(b):
            return (
                pyunits.convert(b.ohmic_resistance * b.conductivity, to_units=pyunits.m)
                == b.electrode_gap
            )

        @self.Constraint(doc="Electrode mass")
        def eq_electrode_mass(b):
            return b.electrode_mass == b.electrode_volume * b.density_electrode_material

    def initialize_build(
        self,
        state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        General wrapper for initialization routines

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = {}).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None)
            solver : str indicating which solver to use during
                     initialization (default = None)

        Returns: None
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        flags = self.properties_in.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
            hold_state=True,
        )
        init_log.info("Initialization Step 1a Complete.")
        # ---------------------------------------------------------------------
        # Initialize other state blocks
        # Set state_args from inlet state
        if state_args is None:
            self.state_args = state_args = {}
            state_dict = self.properties_in[
                self.flowsheet().config.time.first()
            ].define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        self.state_args_out = state_args_out = deepcopy(state_args)

        for j in self.properties_out.component_list:
            if j == "H2O":
                state_args_out["flow_mass_phase_comp"][("Liq", j)] = (
                    state_args["flow_mass_phase_comp"][("Liq", j)]
                    * self.recovery_frac_mass_water.value
                )
                continue
            else:
                state_args_out["flow_mass_phase_comp"][("Liq", j)] = state_args[
                    "flow_mass_phase_comp"
                ][("Liq", j)] * (1 - self.removal_frac_mass_comp[j].value)

        self.properties_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_out,
        )
        init_log.info("Initialization Step 1b Complete.")

        state_args_waste = deepcopy(state_args)
        for j in self.properties_waste.component_list:
            if j == "H2O":
                state_args_waste["flow_mass_phase_comp"][("Liq", j)] = state_args[
                    "flow_mass_phase_comp"
                ][("Liq", j)] * (1 - self.recovery_frac_mass_water.value)
            else:
                state_args_waste["flow_mass_phase_comp"][("Liq", j)] = (
                    state_args["flow_mass_phase_comp"][("Liq", j)]
                    * self.removal_frac_mass_comp[j].value
                )

        self.properties_waste.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_waste,
        )

        init_log.info("Initialization Step 1c Complete.")

        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
            if not check_optimal_termination(res):
                init_log.warning(
                    f"Trouble solving unit model {self.name}, trying one more time"
                )
                res = opt.solve(self, tee=slc.tee)

        init_log.info("Initialization Step 2 {}.".format(idaeslog.condition(res)))
        # ---------------------------------------------------------------------
        # Release Inlet state
        self.properties_in.release_state(flags, outlvl=outlvl)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {self.name} failed to initialize")

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        if iscale.get_scaling_factor(self.electrode_thickness) is None:
            iscale.set_scaling_factor(self.electrode_thickness, 1e3)

        if iscale.get_scaling_factor(self.electrode_mass) is None:
            iscale.set_scaling_factor(self.electrode_mass, 1)

        if iscale.get_scaling_factor(self.electrode_gap) is None:
            iscale.set_scaling_factor(self.electrode_gap, 10)

        if iscale.get_scaling_factor(self.electrolysis_time) is None:
            iscale.set_scaling_factor(self.electrolysis_time, 0.1)

        if iscale.get_scaling_factor(self.applied_current) is None:
            iscale.set_scaling_factor(self.applied_current, 1)

        if iscale.get_scaling_factor(self.current_efficiency) is None:
            iscale.set_scaling_factor(self.current_efficiency, 1)

        if iscale.get_scaling_factor(self.cell_volume) is None:
            iscale.set_scaling_factor(self.cell_volume, 0.1)

        if iscale.get_scaling_factor(self.ohmic_resistance) is None:
            iscale.set_scaling_factor(self.ohmic_resistance, 1e3)

        if iscale.get_scaling_factor(self.charge_loading_rate) is None:
            iscale.set_scaling_factor(self.charge_loading_rate, 1e-2)

        if iscale.get_scaling_factor(self.current_density) is None:
            iscale.set_scaling_factor(self.current_density, 1e-2)

        if iscale.get_scaling_factor(self.cell_voltage) is None:
            iscale.set_scaling_factor(self.cell_voltage, 0.1)

        if iscale.get_scaling_factor(self.overpotential) is None:
            iscale.set_scaling_factor(self.overpotential, 0.1)

        sf = iscale.get_scaling_factor(self.charge_loading_rate)
        iscale.constraint_scaling_transform(self.eq_charge_loading_rate, sf)

        sf = iscale.get_scaling_factor(self.ohmic_resistance)
        iscale.constraint_scaling_transform(self.eq_ohmic_resistance, sf)

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Feed Inlet": self.inlet,
                "Liquid Outlet": self.outlet,
                "Waste Outlet": self.waste,
            },
            time_point=time_point,
        )

    def _get_performance_contents(self, time_point=0):

        # TODO
        var_dict = {}

        return {"vars": var_dict}

    @property
    def default_costing_method(self):
        return cost_electrocoagulation
