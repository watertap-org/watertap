#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
from pyomo.environ import (
    Var,
    Reals,
    Constraint,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    useDefault,
    MaterialFlowBasis,
)
from idaes.core.solvers import get_solver
from idaes.core.util.config import is_physical_parameter_block
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog


__author__ = "Alexander V. Dudchenko"

_log = idaeslog.getLogger(__name__)


@declare_process_block_class("GenericDesalter")
class GenericDesalterData(UnitModelBlockData):
    """
    GenericDesalter - users must provide water recovery
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
        "material_balance_type",
        ConfigValue(
            default=MaterialBalanceType.useDefault,
            domain=In(MaterialBalanceType),
            description="Material balance construction flag",
            doc="""Indicates what type of mass balance should be constructed,
    **default** - MaterialBalanceType.useDefault.
    **Valid values:** {
    **MaterialBalanceType.useDefault - refer to property package for default
    balance type
    **MaterialBalanceType.none** - exclude material balances,
    **MaterialBalanceType.componentPhase** - use phase component balances,
    **MaterialBalanceType.componentTotal** - use total component balances,
    **MaterialBalanceType.elementTotal** - use total element balances,
    **MaterialBalanceType.total** - use total material balance.}""",
        ),
    )
    CONFIG.declare(
        "energy_balance_type",
        ConfigValue(
            default=EnergyBalanceType.useDefault,
            domain=In(EnergyBalanceType),
            description="Energy balance construction flag",
            doc="""Indicates what type of energy balance should be constructed,
    **default** - EnergyBalanceType.useDefault.
    **Valid values:** {
    **EnergyBalanceType.useDefault - refer to property package for default
    balance type
    **EnergyBalanceType.none** - exclude energy balances,
    **EnergyBalanceType.enthalpyTotal** - single enthalpy balance for material,
    **EnergyBalanceType.enthalpyPhase** - enthalpy balances for each phase,
    **EnergyBalanceType.energyTotal** - single energy balance for material,
    **EnergyBalanceType.energyPhase** - energy balances for each phase.}""",
        ),
    )
    CONFIG.declare(
        "momentum_balance_type",
        ConfigValue(
            default=MomentumBalanceType.pressureTotal,
            domain=In(MomentumBalanceType),
            description="Momentum balance construction flag",
            doc="""Indicates what type of momentum balance should be constructed,
    **default** - MomentumBalanceType.pressureTotal.
    **Valid values:** {
    **MomentumBalanceType.none** - exclude momentum balances,
    **MomentumBalanceType.pressureTotal** - single pressure balance for material,
    **MomentumBalanceType.pressurePhase** - pressure balances for each phase,
    **MomentumBalanceType.momentumTotal** - single momentum balance for material,
    **MomentumBalanceType.momentumPhase** - momentum balances for each phase.}""",
        ),
    )
    CONFIG.declare(
        "has_pressure_change",
        ConfigValue(
            default=False,
            domain=In([True, False]),
            description="Pressure change term construction flag",
            doc="""Indicates whether terms for pressure change should be
    constructed,
    **default** - False.
    **Valid values:** {
    **True** - include pressure change terms,
    **False** - exclude pressure change terms.}""",
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
        "tracked_solids_list",
        ConfigValue(
            default=None,
            domain=list,
            description="List of solids that can form in desalter unit",
            doc="""A list of solids that will should be tracked in the waste of desalter unit, 
            if None, then model will not include calculation of water content in solids or 
            solids concentration""",
        ),
    )

    def build(self):
        # Call UnitModel.build to setup dynamics
        super().build()

        # self.scaling_factor = Suffix(direction=Suffix.EXPO)
        self._get_property_package()
        units_meta = self.config.property_package.get_metadata().get_derived_units

        self.water_recovery = Var(
            initialize=80,
            bounds=(None, 99.9999),
            domain=Reals,
            units=pyunits.dimensionless,
            doc="water recovery",
        )
        self.water_recovery.fix()
        if self.config.tracked_solids_list is not None:
            self.brine_solids_concentration = Var(
                initialize=80,
                bounds=(None, None),
                domain=Reals,
                units=pyunits.kg / pyunits.m**3,
                doc="water recovery",
            )
            self.brine_water_mass_percent = Var(
                initialize=80,
                bounds=(None, None),
                domain=Reals,
                units=pyunits.dimensionless,
                doc="water recovery",
            )
        self.brine_unit = ControlVolume0DBlock(
            dynamic=False,
            has_holdup=False,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )
        self.brine_unit.add_state_blocks(has_phase_equilibrium=False)

        self.brine_unit.add_material_balances(
            balance_type=self.config.material_balance_type,
            has_mass_transfer=True,
        )
        self.product_properties = self.config.property_package.build_state_block(
            self.flowsheet().time,
            doc="Material properties in feed",
            defined_state=True,
            has_phase_equilibrium=False,
            **self.config.property_package_args
        )
        self.product_properties[0].define_state_vars()
        self.add_inlet_port(
            name="inlet",
            block=self.brine_unit.properties_in,
        )
        self.add_port(name="brine", block=self.brine_unit.properties_out)
        self.add_port(name="product", block=self.product_properties)

        @self.brine_unit.Constraint(
            self.flowsheet().config.time,
            doc="isothermal energy balance for reactor",
        )
        def eq_isothermal(b, t):
            return b.properties_in[t].temperature == b.properties_out[t].temperature

        @self.brine_unit.Constraint(
            self.flowsheet().config.time,
            doc="isothermal energy balance for reactor",
        )
        def eq_isobaric(b, t):
            return b.properties_in[t].pressure == b.properties_out[t].pressure

        @self.Constraint(
            self.flowsheet().config.time,
            doc="isothermal energy balance for reactor",
        )
        def eq_isothermal(b, t):
            return (
                b.brine_unit.properties_in[t].temperature
                == b.product_properties[t].temperature
            )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="isothermal energy balance for reactor",
        )
        def eq_isobaric(b, t):
            return (
                b.brine_unit.properties_in[t].pressure
                == b.product_properties[t].pressure
            )

        @self.brine_unit.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Mass balance",
        )
        def eq_ion_transfer(b, t, phase, ion):
            if ion == "H2O":
                return Constraint.Skip
            else:
                if (
                    b.config.property_package.config.material_flow_basis
                    == MaterialFlowBasis.mass
                ):
                    return (
                        b.properties_in[t].flow_mass_phase_comp[phase, ion]
                        == b.properties_out[t].flow_mass_phase_comp[phase, ion]
                    )
                elif (
                    b.config.property_package.config.material_flow_basis
                    == MaterialFlowBasis.molar
                ):
                    return (
                        b.properties_in[t].flow_mol_phase_comp[phase, ion]
                        == b.properties_out[t].flow_mol_phase_comp[phase, ion]
                    )

        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Mass balance",
        )
        def eq_product_ion_transfer(b, t, phase, ion):
            if ion == "H2O":
                return Constraint.Skip
            else:
                if (
                    b.config.property_package.config.material_flow_basis
                    == MaterialFlowBasis.mass
                ):
                    return (
                        b.product_properties[t].flow_mass_phase_comp[phase, ion] == 0.0
                    )
                elif (
                    b.config.property_package.config.material_flow_basis
                    == MaterialFlowBasis.molar
                ):
                    return (
                        b.product_properties[t].flow_mass_phase_comp[phase, ion] == 0.0
                    )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Mass balance with reaction terms",
        )
        def eq_water_recovery(b, t):
            return (
                b.brine_unit.properties_in[t].flow_vol_phase["Liq"]
                * (self.water_recovery / 100)
                == b.product_properties[t].flow_vol_phase["Liq"]
            )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Mass balance with reaction terms",
        )
        def eq_brine_flow(b, t):
            return b.brine_unit.properties_out[t].flow_vol_phase[
                "Liq"
            ] == b.brine_unit.properties_in[t].flow_vol_phase["Liq"] * (
                1 - self.water_recovery / 100
            )

        if self.config.tracked_solids_list is not None:

            @self.Constraint(
                self.flowsheet().config.time,
                doc="Mass balance with reaction terms",
            )
            def solids_concentration_eq(b, t):
                return b.brine_solids_concentration * b.brine_unit.properties_out[
                    t
                ].flow_vol_phase["Liq"] == sum(
                    [
                        b.brine_unit.properties_out[t].flow_mass_phase_comp["Liq", tds]
                        for tds in self.config.tracked_solids_list
                    ]
                )

            @self.Constraint(
                self.flowsheet().config.time,
                doc="Mass balance with reaction terms",
            )
            def water_mass_frac_eq(b, t):
                return (
                    b.brine_water_mass_percent
                    * (
                        sum(
                            [
                                b.brine_unit.properties_out[0].flow_mass_phase_comp[
                                    "Liq", tds
                                ]
                                for tds in self.config.tracked_solids_list
                            ]
                        )
                        + b.brine_unit.properties_out[0].flow_mass_phase_comp[
                            "Liq", "H2O"
                        ]
                    )
                    == b.brine_unit.properties_out[0].flow_mass_phase_comp["Liq", "H2O"]
                    * 100
                )

    def initialize_build(
        self, state_args=None, outlvl=idaeslog.NOTSET, solver=None, optarg=None
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
        control_vol_flags = self.brine_unit.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )
        init_log.info_high("Dissolution reactor Step 1 Complete.")
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))
        self.brine_unit.release_state(control_vol_flags, outlvl)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
        for (t, p, j), con in self.brine_unit.eq_ion_transfer.items():
            sf = iscale.get_scaling_factor(
                self.brine_unit.properties_in[t].get_material_flow_terms(p, j)
            )
            iscale.constraint_scaling_transform(con, sf)
        sf = iscale.get_scaling_factor(
            self.brine_unit.properties_in[0].pressure, 1 / 1e5
        )
        iscale.constraint_scaling_transform(self.eq_isobaric[0], sf)
        iscale.constraint_scaling_transform(self.brine_unit.eq_isobaric[0], sf)
        sf = iscale.get_scaling_factor(
            self.brine_unit.properties_in[0].temperature, 1 / 100
        )
        iscale.constraint_scaling_transform(self.brine_unit.eq_isothermal[0], sf)
        iscale.constraint_scaling_transform(self.eq_isothermal[0], sf)
        iscale.constraint_scaling_transform(self.eq_water_recovery[0], 1)
        iscale.set_scaling_factor(self.water_recovery, 1)
