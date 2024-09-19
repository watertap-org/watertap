###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################

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
from idaes.core.solvers.get_solver import get_solver
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.constants import Constants
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.misc import StrEnum
from idaes.core.util.exceptions import ConfigurationError, InitializationError
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.core import InitializationMixin

__author__ = "Kurban Sitterley, Abdiel Lugo"

_log = idaeslog.getLogger(__name__)


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
        cats = self.config.property_package.cation_set

        if "TDS" not in solutes:
            raise ConfigurationError(
                "TDS must be in feed stream for solution conductivity estimation."
            )

        if self.config.electrode_material == ElectrodeMaterial.aluminum:
            if "Al_3+" not in cats:
                raise ConfigurationError(
                    "Electrode material ion must be in feed stream with concentration set to target electrocoagulation dose."
                )

            self.ec_ion = "Al_3+"
            self.density_electrode_material = Param(
                initialize=2710,
                units=pyunits.kg / pyunits.m**3,
                doc="Density of electrode material",
            )

            if self.config.overpotential_calculation == OverpotentialCalculation.nernst:

                self.anode_cell_potential_std = Param(
                    initialize=-1.66,
                    mutable=True,
                    units=pyunits.volt,
                    doc="Anodic non-equilibrium cell potential, standard @ 25C",
                )

                self.anode_entropy_change_std = Param(
                    initialize=0.000533,  # = S / (Z * F)
                    mutable=True,
                    units=pyunits.volt / pyunits.K,
                    doc="Entropy change",
                )

                self.anodic_exchange_current_density = Param(
                    initialize=2.602e-5,
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

        if self.config.electrode_material == ElectrodeMaterial.iron:
            if "Fe_2+" not in cats:
                raise ConfigurationError(
                    "Electrode material ion must be in feed stream with concentration set to target electrocoagulation dose."
                )

            self.ec_ion = "Fe_2+"
            self.density_electrode_material = Param(
                initialize=7860,
                units=pyunits.kg / pyunits.m**3,
                doc="Density of electrode material",
            )

            if self.config.overpotential_calculation == OverpotentialCalculation.nernst:

                self.anode_cell_potential_std = Param(
                    initialize=-0.41,
                    mutable=True,
                    units=pyunits.volt,
                    doc="Anodic non-equilibrium cell potential, standard @ 25C",
                )

                self.anode_entropy_change_std = Param(
                    initialize=7e-5,  # = S / (Z * F)
                    mutable=True,
                    units=pyunits.volt / pyunits.K,
                    doc="Entropy change",
                )

                self.anodic_exchange_current_density = Param(
                    initialize=0.00025,
                    mutable=True,
                    units=pyunits.ampere / pyunits.m**2,
                    doc="Anodic exchange current density",
                )

                self.cathodic_exchange_current_density = Param(
                    initialize=0.001,
                    mutable=True,
                    units=pyunits.ampere / pyunits.m**2,
                    doc="Cathodic exchange current density",
                )

        if self.config.overpotential_calculation == OverpotentialCalculation.nernst:

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

            self.cathode_conc_mol_hydroxide = Param(
                initialize=1e-3,
                mutable=True,
                units=pyunits.mol / pyunits.liter,
                doc="Hydroxide concentration at cathode surface- assume pH = 11",
            )

            self.cathode_conc_mol_metal = Param(
                initialize=0.08,
                mutable=True,
                units=pyunits.mol / pyunits.liter,
                doc="Metal concentration at cathode surface- assume pH = 11",
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

        self.add_port(name="inlet", block=self.properties_in)
        self.add_port(name="outlet", block=self.properties_out)
        self.add_port(name="waste", block=self.properties_waste)

        self.metal_loading = pyunits.convert(
            self.properties_in[0].conc_mass_phase_comp["Liq", self.ec_ion],
            to_units=pyunits.kg / pyunits.liter,
        )

        self.ec_ion_mw = self.config.property_package.mw_comp[self.ec_ion]

        self.ec_ion_z = self.config.property_package.config.charge[self.ec_ion]

        self.current_per_reactor = Param(
            initialize=3000,
            units=pyunits.ampere,
            mutable=True,
            doc="Current required per reactor",
        )

        self.tds_to_cond_conversion = Param(
            initialize=5e3,
            mutable=True,
            units=(pyunits.mg * pyunits.m) / (pyunits.liter * pyunits.S),
            doc="Conersion factor for mg/L TDS to S/m",
        )

        @self.Expression(doc="Conductivity")
        def conductivity(b):
            tds = pyunits.convert(
                b.properties_in[0].conc_mass_phase_comp["Liq", "TDS"],
                to_units=pyunits.mg / pyunits.L,
            )
            return tds / b.tds_to_cond_conversion

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
            doc="Component removal efficiency on mass basis",
        )

        self.recovery_frac_mass_water = Param(
            initialize=0.99, mutable=True, doc="Water recovery on mass basis"
        )

        self.electrode_width = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.m,
            doc="Electrode width",
        )

        self.electrode_height = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.m,
            doc="Electrode height",
        )

        self.electrode_thick = Var(
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

        self.electrode_area_total = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.m**2,
            doc="Total electrode area",
        )

        self.electrode_area_per = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.m**2,
            doc="Electrode area",
        )

        self.electrode_volume_per = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.m**3,
            doc="Electrode volume",
        )

        self.electrode_gap = Var(
            initialize=0.005,
            bounds=(0.0001, 0.2),
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
            doc="Reactor volume total (electrochemical + flotation + sedimentation)",
        )

        self.number_electrode_pairs = Var(
            initialize=5,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Number of electrode pairs",
        )

        self.number_cells = Var(
            initialize=5,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Number of cells",
        )

        self.current_density = Var(
            initialize=1,
            bounds=(1, 2000),
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
            units=pyunits.ohm,
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
            bounds=(0.9, 2.5),
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

        if self.config.overpotential_calculation == OverpotentialCalculation.regression:

            self.overpotential_k1 = Var(
                units=pyunits.millivolt,
                doc="Constant k1 in overpotential equation",
            )

            self.overpotential_k2 = Var(
                units=pyunits.millivolt,
                doc="Constant k2 in overpotential equation",
            )

            @self.Constraint(doc="Overpotential calculation")
            def eq_overpotential(b):
                cd = pyunits.convert(
                    b.current_density, to_units=pyunits.milliampere / pyunits.cm**2
                )
                cd_dimensionless = pyunits.convert(
                    cd * pyunits.cm**2 / pyunits.milliampere,
                    to_units=pyunits.dimensionless,
                )
                ea_tot = pyunits.convert(
                    b.electrode_area_total, to_units=pyunits.cm**2
                )
                return b.overpotential == pyunits.convert(
                    (
                        (
                            cd
                            * (
                                b.overpotential_k1 * log(cd_dimensionless)
                                + b.overpotential_k2
                            )
                        )
                        * ea_tot
                    )
                    / b.applied_current,
                    to_units=pyunits.volt,
                )

        if self.config.overpotential_calculation == OverpotentialCalculation.nernst:

            @self.Constraint(doc="Temperature change")
            def eq_temperature_change(b):
                return (
                    b.properties_out[0].temperature
                    == b.frac_increase_temperature * b.properties_in[0].temperature
                )

            @self.Expression(doc="Change in temperature")
            def delta_temp(b):
                return b.properties_out[0].temperature - b.properties_in[0].temperature

            @self.Expression(doc="Anode equilibrium cell potential")
            def anode_cell_potential(b):
                ccmm_dimensionless = pyunits.convert(
                    b.cathode_conc_mol_metal * pyunits.liter / pyunits.mol,
                    to_units=pyunits.dimensionless,
                )
                return (
                    b.anode_cell_potential_std
                    + b.anode_entropy_change_std * b.delta_temp
                    + (
                        (Constants.gas_constant * b.properties_out[0].temperature)
                        / (Constants.faraday_constant * b.ec_ion_z)
                        * log(ccmm_dimensionless)
                    )
                )

            @self.Expression(doc="Cathode equilibrium cell potential")
            def cathode_cell_potential(b):
                ccmh_dimensionless = pyunits.convert(
                    b.cathode_conc_mol_hydroxide * pyunits.liter / pyunits.mol,
                    to_units=pyunits.dimensionless,
                )
                ph2_dimesionless = pyunits.convert(
                    b.partial_pressure_H2 * pyunits.atm**-1,
                    to_units=pyunits.dimensionless,
                )
                return (
                    b.cathode_cell_potential_std
                    + b.cathode_entropy_change_std * b.delta_temp
                    - (
                        (Constants.gas_constant * b.properties_out[0].temperature)
                        / (Constants.faraday_constant * b.ec_ion_z)
                        * log(ccmh_dimensionless * ph2_dimesionless**2)
                    )
                )

            @self.Expression(doc="Anode overpotential")
            def anode_overpotential(b):
                aecd_dimensionless = pyunits.convert(
                    b.anodic_exchange_current_density * pyunits.m**2 / pyunits.ampere,
                    to_units=pyunits.dimensionless,
                )
                cd_dimensionless = pyunits.convert(
                    b.current_density * pyunits.m**2 / pyunits.ampere,
                    to_units=pyunits.dimensionless,
                )
                return -1 * (
                    (
                        Constants.gas_constant
                        * b.properties_out[0].temperature
                        * log(aecd_dimensionless)
                    )
                    / (Constants.faraday_constant * b.ec_ion_z * 0.5)
                ) + (
                    (
                        Constants.gas_constant
                        * b.properties_out[0].temperature
                        * log(cd_dimensionless)
                    )
                    / (Constants.faraday_constant * b.ec_ion_z * 0.5)
                )

            @self.Expression(doc="Cathode overpotential")
            def cathode_overpotential(b):
                cecd_dimensionless = pyunits.convert(
                    b.cathodic_exchange_current_density
                    * pyunits.m**2
                    / pyunits.ampere,
                    to_units=pyunits.dimensionless,
                )
                cd_dimensionless = pyunits.convert(
                    b.current_density * pyunits.m**2 / pyunits.ampere,
                    to_units=pyunits.dimensionless,
                )
                return (
                    (
                        Constants.gas_constant
                        * b.properties_out[0].temperature
                        * log(cecd_dimensionless)
                    )
                    / (Constants.faraday_constant * b.ec_ion_z * 0.5)
                ) - (
                    (
                        Constants.gas_constant
                        * b.properties_out[0].temperature
                        * log(cd_dimensionless)
                    )
                    / (Constants.faraday_constant * b.ec_ion_z * 0.5)
                )

            @self.Constraint(doc="Overpotential calculation")
            def eq_nernst_overpotential(b):
                return b.overpotential == abs(
                    b.cathode_cell_potential - b.anode_cell_potential
                ) + b.anode_overpotential + abs(b.cathode_overpotential)

        @self.Constraint(doc="Charge loading rate equation")
        def eq_charge_loading_rate(b):
            return b.charge_loading_rate == (
                b.applied_current
                * pyunits.convert(b.electrolysis_time, to_units=pyunits.second)
            ) / pyunits.convert(b.reactor_volume, to_units=pyunits.liter)

        @self.Constraint(doc="Metal loading rate equation")
        def eq_metal_loading_rate(b):
            return b.metal_loading == (
                b.current_efficiency * b.charge_loading_rate * b.ec_ion_mw
            ) / (b.ec_ion_z * Constants.faraday_constant)

        @self.Constraint(doc="Total current required")
        def eq_applied_current(b):
            flow_in = pyunits.convert(
                b.properties_in[0].flow_vol_phase["Liq"],
                to_units=pyunits.liter / pyunits.second,
            )
            return b.applied_current == (
                flow_in * b.metal_loading * b.ec_ion_z * Constants.faraday_constant
            ) / (b.current_efficiency * b.ec_ion_mw)

        @self.Constraint(doc="Total electrode area required")
        def eq_electrode_area_total(b):
            return b.electrode_area_total == b.applied_current / b.current_density

        @self.Constraint(doc="Cell voltage")
        def eq_cell_voltage(b):
            return (
                b.cell_voltage
                == b.overpotential + b.applied_current * b.ohmic_resistance
            )

        @self.Constraint(doc="Area per electrode")
        def eq_electrode_area_per(b):
            return b.electrode_area_per == b.electrode_area_total / (
                b.number_electrode_pairs * 2
            )

        @self.Constraint(doc="Electrode width")
        def eq_electrode_width(b):
            return b.electrode_width == (2 * b.electrode_area_per) ** 0.5

        @self.Constraint(doc="Electrode height")
        def eq_electrode_height(b):
            return b.electrode_height == b.electrode_area_per / b.electrode_width

        @self.Constraint(doc="Electrode volume")
        def eq_electrode_volume_per(b):
            return (
                b.electrode_volume_per
                == b.electrode_width * b.electrode_height * b.electrode_thick
            )

        @self.Constraint(doc="Total reactor volume")
        def eq_reactor_volume(b):
            flow_vol = b.properties_in[0].flow_vol_phase["Liq"]
            return (
                b.reactor_volume
                == pyunits.convert(
                    flow_vol * b.electrolysis_time,
                    to_units=pyunits.m**3,
                )
                / b.number_cells
            )

        @self.Constraint(doc="Ohmic resistance")
        def eq_ohmic_resistance(b):
            return b.ohmic_resistance == b.electrode_gap / (
                b.conductivity * b.electrode_area_per * b.number_cells
            )

        @self.Constraint(doc="Electrode mass")
        def eq_electrode_mass(b):
            return (
                b.electrode_mass
                == b.electrode_volume_per * b.density_electrode_material
            )

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

        ### MASS BALANCE ###

        @self.Constraint(doc="Water recovery")
        def eq_water_recovery(b):
            prop_in = b.properties_in[0]
            prop_out = b.properties_out[0]
            return (
                prop_out.flow_mass_phase_comp["Liq", "H2O"]
                == prop_in.flow_mass_phase_comp["Liq", "H2O"]
                * b.recovery_frac_mass_water
            )

        @self.Constraint(doc="Water mass balance")
        def eq_water_mass_balance(b):
            prop_in = b.properties_in[0]
            prop_out = b.properties_out[0]
            prop_waste = b.properties_waste[0]
            return (
                prop_in.flow_mass_phase_comp["Liq", "H2O"]
                == prop_out.flow_mass_phase_comp["Liq", "H2O"]
                + prop_waste.flow_mass_phase_comp["Liq", "H2O"]
            )

        @self.Constraint(solutes, doc="Component removal")
        def eq_component_removal(b, j):
            prop_in = b.properties_in[0]
            prop_waste = b.properties_waste[0]
            return (
                b.removal_frac_mass_comp[j] * prop_in.flow_mass_phase_comp["Liq", j]
                == prop_waste.flow_mass_phase_comp["Liq", j]
            )

        @self.Constraint(solutes, doc="Component mass balance")
        def eq_comp_mass_balance(b, j):
            prop_in = b.properties_in[0]
            prop_out = b.properties_out[0]
            prop_waste = b.properties_waste[0]
            return (
                prop_in.flow_mass_phase_comp["Liq", j]
                == prop_out.flow_mass_phase_comp["Liq", j]
                + prop_waste.flow_mass_phase_comp["Liq", j]
            )

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
                state_args_out["flow_mol_phase_comp"][("Liq", j)] = (
                    state_args["flow_mol_phase_comp"][("Liq", j)]
                    * self.recovery_frac_mass_water.value
                )
                continue
            else:
                state_args_out["flow_mol_phase_comp"][("Liq", j)] = state_args[
                    "flow_mol_phase_comp"
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
                state_args_waste["flow_mol_phase_comp"][("Liq", j)] = state_args[
                    "flow_mol_phase_comp"
                ][("Liq", j)] * (1 - self.recovery_frac_mass_water.value)
            else:
                state_args_waste["flow_mol_phase_comp"][("Liq", j)] = (
                    state_args["flow_mol_phase_comp"][("Liq", j)]
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

        iscale.set_scaling_factor(self.electrode_width, 10)

        iscale.set_scaling_factor(self.electrode_height, 10)

        iscale.set_scaling_factor(self.electrode_thick, 1e3)

        iscale.set_scaling_factor(self.electrode_mass, 1)

        iscale.set_scaling_factor(self.electrode_area_total, 1e-2)

        iscale.set_scaling_factor(self.electrode_area_per, 10)

        iscale.set_scaling_factor(self.electrode_volume_per, 1)

        iscale.set_scaling_factor(self.electrode_gap, 10)

        iscale.set_scaling_factor(self.electrolysis_time, 0.1)

        iscale.set_scaling_factor(self.number_electrode_pairs, 0.1)

        iscale.set_scaling_factor(self.number_cells, 1)

        iscale.set_scaling_factor(self.applied_current, 1e-4)

        iscale.set_scaling_factor(self.current_efficiency, 1)

        iscale.set_scaling_factor(self.reactor_volume, 0.1)

        iscale.set_scaling_factor(self.ohmic_resistance, 1e5)

        iscale.set_scaling_factor(self.charge_loading_rate, 1e-2)

        iscale.set_scaling_factor(self.current_density, 1e-2)

        iscale.set_scaling_factor(self.cell_voltage, 0.1)

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
