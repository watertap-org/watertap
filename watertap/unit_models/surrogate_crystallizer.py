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

# Import Pyomo libraries
from pyomo.environ import (
    Set,
    Var,
    check_optimal_termination,
    Param,
    Constraint,
    Suffix,
    units as pyunits,
    NonNegativeReals,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.solvers import get_solver
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.constants import Constants
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import InitializationError
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.core.surrogate.surrogate_block import SurrogateBlock
from idaes.core.surrogate.keras_surrogate import KerasSurrogate
from idaes.core.surrogate.sampling.scaling import OffsetScaler
from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.core.util.initialization import assert_degrees_of_freedom, check_solve

# import watertap.property_models.cryst_prop_pack as props2
import watertap.property_models.water_prop_pack as props3
from watertap.costing.unit_models.surrogate_crystallizer import (
    cost_surrogate_crystallizer,
)


_log = idaeslog.getLogger(__name__)

__author__ = "Oluwamayowa Amusat, Adam Atia"


@declare_process_block_class("SurrogateCrystallizer")
class SurrogateCrystallizerData(UnitModelBlockData):
    """
    ML-based crystallizer model
    """

    CONFIG = ConfigBlock()

    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
    **default** = False. NF units do not support dynamic
    behavior.""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Indicates whether holdup terms should be constructed or not.
    **default** - False. NF units do not have defined volume, thus
    this must be False.""",
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
        "vapor_property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for vapor phase",
            doc="""Property parameter object used to define property calculations
    for the vapor phase,
    **default** - useDefault.
    **Valid values:** {
    **useDefault** - use default package from parent model or flowsheet,
    **PropertyParameterObject** - a PropertyParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "vapor_property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing vapor phase properties",
            doc="""A ConfigBlock with arguments to be passed to vapor phase
property block(s) and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
        ),
    )
    CONFIG.declare(
        "input_ion_list",
        ConfigValue(
            default={},
            domain=list,
            description="Input ions",
            doc="""list of ions in feed""",
        ),
    )
    CONFIG.declare(
        "solids_ions_dict",
        ConfigValue(
            default={},
            domain=dict,
            description="Specification of ion makeup of each solid in system",
            doc="""Solids makeup""",
        ),
    )

    def build(self):


        # Call UnitModel.build to setup dynamics
        super().build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)
        units_meta = self.config.property_package.get_metadata().get_derived_units

        self.ion_list = Set(
            initialize=self.config.input_ion_list
        )  # Could extract from property package
        self.solids_list = Set(initialize=self.config.solids_ions_dict.keys())

        # Add ion in solid ratios as parameters
        dict1 = dict()
        for m in self.solids_list:
            for p in self.ion_list:
                dict1[m, p] = float(
                    self.config.solids_ions_dict[m][p]
                    if p in self.config.solids_ions_dict[m]
                    else 0.0
                )
        self.mwc = Param(self.solids_list, self.ion_list, initialize=dict1)

        # Add pther variables
        self.temperature_operating = Var(
            initialize=298.15,
            domain=NonNegativeReals,
            units=pyunits.K,
            bounds=(273, 1000),
            doc="Crystallizer operating temperature in K",
        )
        self.pressure_operating = Var(
            initialize=101325,
            domain=NonNegativeReals,
            units=pyunits.Pa,
            bounds=(0.001, 1e6),
            doc="Crystallizer pressure in Pa",
        )
        self.evaporation_percent = Var(
            initialize=50,
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            bounds=(10, 95),
            doc="Crystallizer percentage of water evaporation",
        )
        self.S = Var(
            self.ion_list,
            initialize=0,
            domain=NonNegativeReals,
            units=pyunits.kg / pyunits.s,
            doc="True species solid mass flow (kg)",
        )
        self.mixed_solids = Var(
            self.solids_list,
            initialize=0,
            domain=NonNegativeReals,
            units=pyunits.kg / pyunits.s,
            doc="Apparent species solid mass flow (kg)",
        )
        self.S_total = Var(
            initialize=0,
            domain=NonNegativeReals,
            units=pyunits.kg / pyunits.s,
            doc="Total solid flow",
        )
        self.L_total = Var(
            initialize=0,
            domain=NonNegativeReals,
            units=pyunits.kg / pyunits.s,
            doc="Total liquid flow",
        )
        self.V_total = Var(
            initialize=0,
            domain=NonNegativeReals,
            units=pyunits.kg / pyunits.s,
            doc="Total vapor flow",
        )

        self.Q = Var(
            initialize=1e5,
            domain=NonNegativeReals,
            units=pyunits.kW,
            bounds=(-5e6, 5e6),
            doc="Heat requirement for crystallizer (kW)",
        )

        # # Add state blocks for inlet and three outlets
        # Add inlet block
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = True  # inlet block is an inlet
        self.properties_in = self.config.property_package.state_block_class(
            self.flowsheet().config.time, doc="Material properties of inlet", **tmp_dict
        )

        # Add liquid outlet
        tmp_dict["defined_state"] = False  # outlet and waste block is not an inlet
        self.properties_out_liq = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of liquid outlet",
            **tmp_dict,
        )

        # Add vapor outlet
        tmp_dict2 = dict()
        tmp_dict2["has_phase_equilibrium"] = False
        tmp_dict2["parameters"] = self.config.vapor_property_package
        tmp_dict["defined_state"] = False
        self.properties_out_vapor = self.config.vapor_property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of vapor outlet",
            **tmp_dict2,
        )

        # Add ports - oftentimes users interact with these rather than the state blocks
        self.add_port(name="inlet", block=self.properties_in)
        self.add_port(name="outlet", block=self.properties_out_liq)
        self.add_port(name="vapor", block=self.properties_out_vapor)

        # Add constraints
        # 1. Fix empty flows:
        
        for p in self.config.vapor_property_package.phase_list:
            if p == "Liq":
                self.properties_out_vapor[0].flow_mass_phase_comp[p, "H2O"].fix(0.0)

        # 2. Material balances
        @self.Constraint(
            self.config.property_package.component_list,
            doc="Mass balance for components",
        )
        def eq_mass_balance_constraints(b, j):
            if j != "H2O":
                return (
                    sum(
                        b.properties_in[0].flow_mass_phase_comp[p, j]
                        for p in self.config.property_package.phase_list
                    )
                    == sum(
                        b.properties_out_liq[0].flow_mass_phase_comp[p, j]
                        for p in self.config.property_package.phase_list
                    )
                    + b.S[j]
                )
            else:
                return (
                    sum(
                        b.properties_in[0].flow_mass_phase_comp[p, j]
                        for p in self.config.property_package.phase_list
                    )
                    == sum(
                        b.properties_out_liq[0].flow_mass_phase_comp[p, j]
                        for p in self.config.property_package.phase_list
                    )
                    + sum(
                        b.properties_out_vapor[0].flow_mass_phase_comp[p, j]
                        for p in self.config.vapor_property_package.phase_list
                    )
                    # + b.S[j]
                )

        # 3. Temperature and pressure balances:
        @self.Constraint()
        def eq_T_con1(b):
            return self.temperature_operating == b.properties_out_liq[0].temperature

        @self.Constraint()
        def eq_T_con2(b):
            return self.temperature_operating == b.properties_out_vapor[0].temperature

        @self.Constraint()
        def eq_P_con1(b):
            return self.pressure_operating == b.properties_out_liq[0].pressure

        @self.Constraint()
        def eq_P_con2(b):
            return self.pressure_operating == b.properties_out_vapor[0].pressure

        @self.Constraint(
            self.config.property_package.component_list,
            doc="Liquid-Vapour phase water split",
        )
        def eq_water_balance_constraints(b, j):
            if j == "H2O":
                return (
                    sum(
                        b.properties_out_vapor[0].flow_mass_phase_comp[p, j]
                        for p in self.config.vapor_property_package.phase_list
                    )
                    == b.properties_in[0].flow_mass_phase_comp["Liq", j]
                    * b.evaporation_percent
                    / 100
                )
            else:
                return Constraint.Skip

        @self.Constraint(
            self.config.property_package.component_list,
            doc="Mass balance for components",
        )
        def eq_solid_amount(b, j):
            if j == "H2O":
                return Constraint.Skip
            else:
                return b.S[j] == sum(
                    (b.mwc[q, j] * b.properties_in[0].params.mw_comp[j])
                    / sum(
                        b.mwc[q, j] * b.properties_in[0].params.mw_comp[j]
                        for j in b.ion_list
                    )
                    * b.mixed_solids[q]
                    for q in b.solids_list
                )

        # 5. Constraints calculating total flows for each outlet stream
        @self.Constraint(
            doc="Total solid flow constraint",
        )
        def eq_total_solids_constraint(b):
            return b.S_total == sum(b.mixed_solids[k] for k in b.solids_list)

        @self.Constraint(
            doc="Total liquid flow constraint",
        )
        def eq_total_liquid_constraint(b):
            return b.L_total == sum(
                sum(
                    b.properties_out_liq[0].flow_mass_phase_comp[p, j]
                    for p in self.config.property_package.phase_list
                )
                for j in self.config.property_package.component_list
            )

        @self.Constraint(
            doc="Total vapor flow constraint",
        )
        def eq_total_vapor_constraint(b):
            return b.V_total == sum(
                b.properties_out_vapor[0].flow_mass_phase_comp[p, "H2O"]
                for p in self.config.vapor_property_package.phase_list
            )

        # 6. Constraints for computing heat requirement
        @self.Constraint(
            doc="Energy balance",
        )
        def eq_energy_balance_constraint(b):
            return pyunits.convert(b.Q, to_units=pyunits.W) + b.properties_in[
                0
            ].enth_flow == b.properties_out_liq[0].enth_flow + sum(
                b.properties_out_vapor[0].enth_flow_phase[p]
                for p in self.config.vapor_property_package.phase_list
            )

    def initialize_build(
        self,
        state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        General wrapper for pressure changer initialization routines

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
        # Initialize holdup block
        flags = self.properties_in.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
            hold_state=True,
        )
        init_log.info_high("Initialization Step 1 Complete.")

        # ---------------------------------------------------------------------
        # Initialize other state blocks
        # Set state_args from inlet state
        if state_args is None:
            state_args = {}
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

        self.properties_out_liq.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )

        state_args_vapor = {}
        state_dict_vapor = self.properties_out_vapor[
            self.flowsheet().config.time.first()
        ].define_port_members()

        for k in state_dict_vapor.keys():
            if state_dict_vapor[k].is_indexed():
                state_args_vapor[k] = {}
                for m in state_dict_vapor[k].keys():
                    state_args_vapor[k][m] = state_dict_vapor[k][m].value
            else:
                state_args_vapor[k] = state_dict_vapor[k].value
        for p in self.config.vapor_property_package.phase_list:
            # for j in self.ion_list:
            if p == "Vap":
                state_args_vapor["flow_mass_phase_comp"][p, "H2O"] = state_args[
                    "flow_mass_phase_comp"
                ]["Liq", "H2O"]
            else:
                state_args_vapor["flow_mass_phase_comp"][p, "H2O"] = 0
        self.properties_out_vapor.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_vapor,
        )
        init_log.info_high("Initialization Step 2 Complete.")
        # ---------------------------------------------------------------------
        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=True)  # slc.tee)
        # res = fsTools.standard_solve(self, tee=True, check_close_to_bounds=True)
        # self.display()
        # self.report()
        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))
        # ---------------------------------------------------------------------
        # Release Inlet state
        self.properties_in.release_state(flags, outlvl=outlvl)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

        # from pyomo.util.infeasible import (
        #     log_active_constraints,
        #     log_close_to_bounds,
        #     log_infeasible_bounds,
        #     log_infeasible_constraints,
        # )
        # from pyomo.common.log import LoggingIntercept
        # import logging
        # from io import StringIO

        # output = StringIO()
        # with LoggingIntercept(output, "pyomo.util.infeasible", logging.INFO):
        #     log_infeasible_constraints(self)
        # print(output.getvalue().splitlines())

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {self.name} failed to initialize")

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Feed Inlet": self.inlet,
                "Liquid Outlet": self.outlet,
                "Vapor Outlet": self.vapor,
            },
            time_point=time_point,
        )

    def _get_performance_contents(self, time_point=0):
        var_dict = {}
        var_dict["Operating Temperature (K)"] = self.temperature_operating
        var_dict["Vapor Pressure (Pa)"] = self.pressure_operating
        var_dict["Total solids at outlet (Kg)"] = self.S_total
        var_dict["Total liquid flow at outlet (Kg)"] = self.L_total
        var_dict["Total water vapor flow at outlet (Kg)"] = self.V_total
        var_dict["Heat requirement"] = self.Q
        return {"vars": var_dict}

    def calculate_scaling_factors(self):

        super().calculate_scaling_factors()

        iscale.set_scaling_factor(
            self.temperature_operating,
            iscale.get_scaling_factor(self.properties_in[0].temperature),
        )
        iscale.set_scaling_factor(self.pressure_operating, 1e-3)
        iscale.set_scaling_factor(
            self.Q,
            iscale.get_scaling_factor(
                self.properties_in[0].flow_mass_phase_comp["Liq", "H2O"]
            )
            * iscale.get_scaling_factor(self.properties_in[0].enth_mass_phase["Liq"])
            * 1e3,
        )

        # transforming constraints
        for ind, c in self.eq_T_con1.items():
            sf = iscale.get_scaling_factor(self.properties_in[0].temperature)
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_T_con2.items():
            sf = iscale.get_scaling_factor(self.properties_in[0].temperature)
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_P_con1.items():
            sf = iscale.get_scaling_factor(self.properties_in[0].pressure)
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_P_con2.items():
            sf = iscale.get_scaling_factor(self.properties_in[0].pressure)
            iscale.constraint_scaling_transform(c, sf)

        for j, c in self.eq_mass_balance_constraints.items():
            sf = iscale.get_scaling_factor(
                self.properties_in[0].flow_mass_phase_comp["Liq", j]
            )
            iscale.constraint_scaling_transform(c, sf)

        for j, c in self.eq_water_balance_constraints.items():
            sf = iscale.get_scaling_factor(
                self.properties_in[0].flow_mass_phase_comp["Liq", "H2O"]
            )
            iscale.constraint_scaling_transform(c, sf)

        for j, c in self.eq_energy_balance_constraint.items():
            sf = iscale.get_scaling_factor(
                self.properties_in[0].flow_mass_phase_comp["Liq", "H2O"]
            ) * iscale.get_scaling_factor(self.properties_in[0].enth_mass_phase["Liq"])
            iscale.constraint_scaling_transform(c, sf)

    @property
    def default_costing_method(self):
        return cost_surrogate_crystallizer
