#################################################################################
# WaterTAP Copyright (c) 2020-2025, The Regents of the University of California,
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
    Var,
    check_optimal_termination,
    Param,
    Suffix,
    log,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.constants import Constants
from idaes.core.util.exceptions import ConfigurationError, InitializationError
from idaes.core.util.misc import StrEnum
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
    carbon_steel = "carbon_steel"
    stainless_steel = "stainless_steel"


class OverpotentialCalculation(StrEnum):
    detailed = "detailed"
    fixed = "fixed"
    regression = "regression"


"""
References:

K. L. Dubrawski, C. Du and M. Mohseni (2014)
Electrochimica Acta 2014 Vol. 129 Pages 187-195
DOI: 10.1016/j.electacta.2014.02.089

Z. Gu, Z. Liao, M. Schulz, J. R. Davis, J. C. Baygents and J. Farrell (2009)
Industrial & Engineering Chemistry Research 2009 Vol. 48 Issue 6 Pages 3112-3117
DOI: 10.1021/ie801086c

Bratsch, S. G. (1989). 
Standard Electrode Potentials and Temperature Coefficients in Water at 298.15 K. 
Journal of Physical and Chemical Reference Data, 18(1), 1-21. 
DOI: 10.1063/1.555839 

Zhang, F., Yang, C., Zhu, H., Li, Y., & Gui, W. (2020). 
An integrated prediction model of heavy metal ion concentration for 
  iron electrocoagulation process. 
Chemical Engineering Journal, 391, 123628. 
DOI: 10.1016/j.cej.2019.123628 

"""


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
            default="carbon_steel",
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
        tmp_dict["defined_state"] = True
        self.properties_in = self.config.property_package.state_block_class(
            self.flowsheet().config.time, doc="Material properties of inlet", **tmp_dict
        )

        tmp_dict["defined_state"] = False
        self.properties_out = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of outlet",
            **tmp_dict,
        )

        self.properties_byproduct = self.config.property_package.state_block_class(
            self.flowsheet().config.time, doc="Material properties of waste", **tmp_dict
        )

        prop_in = self.properties_in[0]
        prop_out = self.properties_out[0]
        prop_byproduct = self.properties_byproduct[0]

        self.add_port(name="inlet", block=self.properties_in)
        self.add_port(name="outlet", block=self.properties_out)
        self.add_port(name="byproduct", block=self.properties_byproduct)

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

        self.recovery_frac_mass_H2O = Param(
            initialize=0.99,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Water recovery on mass basis",
        )

        self.tds_to_cond_conversion = Param(
            initialize=5e3,
            mutable=True,
            units=(pyunits.mg * pyunits.m) / (pyunits.liter * pyunits.S),
            doc="Conversion factor for mg/L TDS to S/m",
        )

        self.standard_temperature = Param(
            initialize=298.15,  # 25C
            mutable=True,
            units=pyunits.degK,
            doc="Standard temperature",
        )

        self.mw_electrode_material = Param(
            initialize=100e-3,
            units=pyunits.kg / pyunits.mol,
            doc="Molecular weight of coagulant species",
        )

        self.charge_transfer_number = Param(
            initialize=1,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Charge transfer number of coagulant species",
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

        self.frac_increase_temperature = Param(
            initialize=1.05,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Fractional increase in water temperature from inlet to outlet",
        )

        if self.config.overpotential_calculation == OverpotentialCalculation.detailed:

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
                initialize=1e-4,
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

            self.tafel_slope_cathode = Var(
                initialize=0.146 / 2.303,  # default value used log10
                bounds=(0, None),
                units=pyunits.volt,  # volt per 10x change (decade)
                doc="Tafel slope for cathode",
            )

            self.tafel_slope_anode = Var(
                initialize=0.093 / 2.303,  # default value used log10
                bounds=(0, None),
                units=pyunits.volt,  # volt per 10x change (decade)
                doc="Tafel slope for anode",
            )

        if self.config.electrode_material == ElectrodeMaterial.aluminum:

            """
            Reaction at anode is:
            Al_3+ + 3e_- <--> Al(s);  Ea0 = -1.66 V
            Reaction in bulk solution is:
            Al_3+ + 3OH_- <--> Al(OH)3(s)
            """

            self.mw_electrode_material.set_value(26.98e-3)
            self.charge_transfer_number.set_value(3)
            self.stoic_coeff.set_value(1)
            self.density_electrode_material.set_value(2710)

            if (
                self.config.overpotential_calculation
                == OverpotentialCalculation.detailed
            ):
                # Dubrawski et al., 2014; Zhang et al., 2020
                self.anode_cell_potential_std.set_value(-1.66)
                # = S / (Z * F); Bratsch, 1989; Dubrawski et al., 2014
                self.anode_entropy_change_std.set_value(0.000533)
                # Electrochemical Thermodynamics and Kinetics, pg 276
                self.anodic_exchange_current_density.set_value(2.602e-5)
                # Electrochemical Thermodynamics and Kinetics, pg 276
                self.cathodic_exchange_current_density.set_value(1.0e-4)

        if self.config.electrode_material == ElectrodeMaterial.iron:

            """
            Reaction at anode is:
            Fe_2+ + 2e_- <--> Fe(s);  Ea0 = -0.41 V
            Reaction in bulk solution is:
            Fe_2+ + 2OH_- <--> Fe(OH)2(s)
            """

            self.mw_electrode_material.set_value(55.845e-3)
            self.charge_transfer_number.set_value(2)
            self.stoic_coeff.set_value(1)
            self.density_electrode_material.set_value(7860)

            if (
                self.config.overpotential_calculation
                == OverpotentialCalculation.detailed
            ):
                # Dubrawski et al., 2014; Zhang et al., 2020
                self.anode_cell_potential_std.set_value(-0.41)
                # = S / (Z * F); Bratsch, 1989; Dubrawski et al., 2014
                self.anode_entropy_change_std.set_value(7e-5)
                # Dubrawski et al., 2014
                self.anodic_exchange_current_density.set_value(2.5e-4)
                # Dubrawski et al., 2014
                self.cathodic_exchange_current_density.set_value(1e-3)

        self.floc_basin_vol = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.m**3,
            doc="Floc basin volume total (flotation + sedimentation)",
        )

        self.floc_retention_time = Var(
            initialize=30,
            bounds=(2, 200),
            units=pyunits.minute,
            doc="Floc basin retention time",
        )

        self.coagulant_dose = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.g / pyunits.liter,
            doc="Assumed coagulant dose to the feed in g/L",
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

        self.reactor_volume = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.m**3,
            doc="Volume of electrocoagulation reactor",
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
            initialize=0.1,
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
                (b.applied_current * b.mw_electrode_material)
                / (
                    Constants.faraday_constant
                    * b.charge_transfer_number
                    * prop_in.flow_vol_phase["Liq"]
                ),
                to_units=pyunits.g / pyunits.liter,
            )

        @self.Expression(doc="Total electrode area")
        def electrode_area_total(b):
            return b.anode_area + b.cathode_area

        @self.Expression(doc="Power required")
        def power_required(b):
            return pyunits.convert(
                b.cell_voltage * b.applied_current, to_units=pyunits.kW
            )

        @self.Expression(doc="Power density Faradaic")
        def power_density_faradaic(b):
            return pyunits.convert(
                (b.overpotential * b.applied_current) / b.anode_area,
                to_units=pyunits.microwatts / pyunits.cm**2,
            )

        @self.Expression(doc="Power density total")
        def power_density_total(b):
            return pyunits.convert(
                b.power_required / b.anode_area,
                to_units=pyunits.microwatts / pyunits.cm**2,
            )

        if self.config.overpotential_calculation == OverpotentialCalculation.regression:

            self.overpotential_k1 = Var(
                units=pyunits.millivolt,
                doc="Constant k1 in overpotential equation from Gu et al. (2009)",
            )

            self.overpotential_k2 = Var(
                units=pyunits.millivolt,
                doc="Constant k2 in overpotential equation from Gu et al. (2009)",
            )

            @self.Constraint(
                doc="Overpotential calculation - adapted from Eq. 18 in Gu et al. (2009)"
            )
            def eq_overpotential(b):
                cd_dimensionless = pyunits.convert(
                    b.current_density * pyunits.cm**2 / pyunits.milliampere,
                    to_units=pyunits.dimensionless,
                )
                return b.overpotential == pyunits.convert(
                    (b.overpotential_k1 * log(cd_dimensionless) + b.overpotential_k2),
                    to_units=pyunits.volt,
                )

        if self.config.overpotential_calculation == OverpotentialCalculation.detailed:

            @self.Expression(doc="Hydroxide concentration at cathode surface")
            def cathode_conc_mol_hydroxide(b):
                pOH = 14 - b.cathode_surface_pH
                return 10 ** (-pOH) * pyunits.mol / pyunits.liter

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

            @self.Expression(
                doc="Anode non-equilibrium cell potential via Nernst equation"
            )
            def anode_cell_potential(b):
                ec_dose_dimensionless = pyunits.convert(
                    (b.coagulant_dose / b.mw_electrode_material)
                    * pyunits.liter
                    * pyunits.mol**-1,
                    to_units=pyunits.dimensionless,
                )
                return b.anode_cell_potential_temp_adj - (
                    pyunits.convert(
                        (
                            (Constants.gas_constant * prop_out.temperature)
                            / (Constants.faraday_constant * b.charge_transfer_number)
                        ),
                        to_units=pyunits.volt,
                    )
                    * log(ec_dose_dimensionless**-b.stoic_coeff)
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

            @self.Expression(doc="Cathode cell potential via Nernst equation")
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
                        (Constants.gas_constant * prop_out.temperature)
                        / (Constants.faraday_constant * b.charge_transfer_number)
                    ),
                    to_units=pyunits.volt,
                ) * log(ccmh_dimensionless**2 * ph2_dimensionless)

            # Eq. 6 in Zhang et al. (2020); Eq. 12 in Dubrawski et al. (2014)
            @self.Expression(doc="Anode activation overpotential")
            def anode_overpotential(b):
                return b.tafel_slope_anode * log(
                    pyunits.convert(
                        b.current_density / b.anodic_exchange_current_density,
                        to_units=pyunits.dimensionless,
                    )
                )

            # Eq. 7 in Zhang et al. (2020); Eq. 11 in Dubrawski et al. (2014)
            @self.Expression(doc="Cathode activation overpotential")
            def cathode_overpotential(b):
                return -1 * (
                    b.tafel_slope_cathode
                    * log(
                        pyunits.convert(
                            b.current_density / b.cathodic_exchange_current_density,
                            to_units=pyunits.dimensionless,
                        )
                    )
                )

            # See Eq. 3 in Zhang et al. (2020), Eq. 10 in Dubrawski et al. (2014)
            # neglecting concentration/mass-transfer overpotential
            @self.Constraint(doc="Overpotential calculation")
            def eq_total_overpotential(b):
                return b.overpotential == abs(
                    b.cathode_cell_potential - b.anode_cell_potential
                ) + b.anode_overpotential + abs(b.cathode_overpotential)

        @self.Constraint(doc="Temperature change")
        def eq_temperature_change(b):
            return (
                prop_out.temperature
                == b.frac_increase_temperature * prop_in.temperature
            )

        @self.Constraint(doc="Water recovery")
        def eq_water_recovery(b):
            return (
                prop_out.flow_mass_phase_comp["Liq", "H2O"]
                == prop_in.flow_mass_phase_comp["Liq", "H2O"] * b.recovery_frac_mass_H2O
            )

        @self.Constraint(doc="Water mass balance")
        def eq_water_mass_balance(b):
            return (
                prop_in.flow_mass_phase_comp["Liq", "H2O"]
                == prop_out.flow_mass_phase_comp["Liq", "H2O"]
                + prop_byproduct.flow_mass_phase_comp["Liq", "H2O"]
            )

        @self.Constraint(solutes, doc="Component removal")
        def eq_component_removal(b, j):
            return (
                b.removal_frac_mass_comp[j] * prop_in.flow_mass_phase_comp["Liq", j]
                == prop_byproduct.flow_mass_phase_comp["Liq", j]
            )

        @self.Constraint(solutes, doc="Component mass balance")
        def eq_comp_mass_balance(b, j):
            return (
                prop_in.flow_mass_phase_comp["Liq", j]
                == prop_out.flow_mass_phase_comp["Liq", j]
                + prop_byproduct.flow_mass_phase_comp["Liq", j]
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
                prop_in.flow_vol_phase["Liq"] * b.floc_retention_time,
                to_units=pyunits.m**3,
            )

        @self.Constraint(doc="Faraday's Law")
        def eq_faraday(b):
            return (
                b.applied_current * b.current_efficiency * b.mw_electrode_material
                == pyunits.convert(
                    prop_in.flow_vol_phase["Liq"]
                    * b.coagulant_dose
                    * b.charge_transfer_number
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
            return b.reactor_volume == pyunits.convert(
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
                    * self.recovery_frac_mass_H2O.value
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
        for j in self.properties_byproduct.component_list:
            if j == "H2O":
                state_args_waste["flow_mass_phase_comp"][("Liq", j)] = state_args[
                    "flow_mass_phase_comp"
                ][("Liq", j)] * (1 - self.recovery_frac_mass_H2O.value)
            else:
                state_args_waste["flow_mass_phase_comp"][("Liq", j)] = (
                    state_args["flow_mass_phase_comp"][("Liq", j)]
                    * self.removal_frac_mass_comp[j].value
                )

        self.properties_byproduct.initialize(
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

        if iscale.get_scaling_factor(self.reactor_volume) is None:
            iscale.set_scaling_factor(self.reactor_volume, 0.1)

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
                "Byproduct Outlet": self.byproduct,
            },
            time_point=time_point,
        )

    def _get_performance_contents(self, time_point=0):

        var_dict = {}
        var_dict["Voltage Required"] = self.cell_voltage
        var_dict["Overpotential"] = self.overpotential
        var_dict["Applied Current"] = self.applied_current
        var_dict["Current Density"] = self.current_density
        var_dict["Charge Loading Rate"] = self.charge_loading_rate
        var_dict["Ohmic Resistance"] = self.ohmic_resistance
        var_dict["Coagulant Dose"] = self.coagulant_dose
        var_dict["Electrode Mass"] = self.electrode_mass
        var_dict["Electrode Volume"] = self.electrode_volume
        var_dict["Electrode Thickness"] = self.electrode_thickness
        var_dict["Electrode Gap"] = self.electrode_gap
        var_dict["Anode Area"] = self.anode_area
        var_dict["Cathode Area"] = self.cathode_area
        var_dict["Electrolysis Time"] = self.electrolysis_time
        var_dict["Current Efficiency"] = self.current_efficiency
        var_dict["Flocculation Basin Volume"] = self.floc_basin_vol
        var_dict["Flocculation Retention Time"] = self.floc_retention_time
        var_dict["EC Reactor Volume"] = self.reactor_volume

        if self.config.overpotential_calculation == OverpotentialCalculation.regression:
            var_dict["Overpotential Equation k1 Parameter"] = self.overpotential_k1
            var_dict["Overpotential Equation k2 Parameter"] = self.overpotential_k2

        if self.config.overpotential_calculation == OverpotentialCalculation.detailed:
            var_dict["Electrode Material Tafel Slope Anode"] = self.tafel_slope_anode
            var_dict["Electrode Material Tafel Slope Cathode"] = (
                self.tafel_slope_cathode
            )

        expr_dict = {}
        expr_dict["Conductivity"] = self.conductivity
        expr_dict["Total Power Required"] = self.power_required
        expr_dict["Power Density Faradaic"] = self.power_density_faradaic
        expr_dict["Power Density Total"] = self.power_density_total
        expr_dict["Total Electrode Area"] = self.electrode_area_total

        if self.config.overpotential_calculation == OverpotentialCalculation.detailed:
            expr_dict["Temp. Adjusted Std Anode Potential"] = (
                self.anode_cell_potential_temp_adj
            )
            expr_dict["Temp. Adjusted Std Cathode Potential"] = (
                self.cathode_cell_potential_temp_adj
            )
            expr_dict["Anode Cell Potential"] = self.anode_cell_potential
            expr_dict["Cathode Cell Potential"] = self.cathode_cell_potential
            expr_dict["Anode Activation Overpotential"] = self.anode_overpotential
            expr_dict["Cathode Activation Overpotential"] = self.cathode_overpotential

        param_dict = {}
        for s in self.config.property_package.solute_set:
            param_dict["Removal Fraction " + s] = self.removal_frac_mass_comp[s]
        param_dict["Water Recovery Fraction"] = self.recovery_frac_mass_H2O
        param_dict["TDS to Conductivity Conversion"] = self.tds_to_cond_conversion
        param_dict["Standard Temperature"] = self.standard_temperature
        param_dict["Electrode Material Density"] = self.density_electrode_material
        param_dict["Electrode Material MW"] = self.mw_electrode_material
        param_dict["Electrode Material Stoichiometric Coefficient"] = self.stoic_coeff
        param_dict["Electrode Material Charge Transfer Number"] = (
            self.charge_transfer_number
        )
        if self.config.overpotential_calculation == OverpotentialCalculation.detailed:
            param_dict["Electrode Material Standard Anodic Potential"] = (
                self.anode_cell_potential_std
            )
            param_dict["Electrode Material Anodic Exchange Current Density"] = (
                self.anodic_exchange_current_density
            )
            param_dict["Electrode Material Anodic Entropy Change"] = (
                self.anode_entropy_change_std
            )
            param_dict["Electrode Material Standard Cathodic Potential"] = (
                self.cathode_cell_potential_std
            )
            param_dict["Electrode Material Cathodic Exchange Current Density"] = (
                self.cathodic_exchange_current_density
            )
            param_dict["Electrode Material Cathodic Entropy Change"] = (
                self.cathode_entropy_change_std
            )

        return {"vars": var_dict, "params": param_dict, "exprs": expr_dict}

    @property
    def default_costing_method(self):
        return cost_electrocoagulation
