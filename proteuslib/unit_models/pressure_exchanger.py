##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################

# Import Pyomo libraries
from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.environ import Block, Constraint, Var, Suffix, NonNegativeReals, Reals, \
    SolverFactory, units as pyunits
from pyomo.network import Port

# Import IDAES cores
import idaes.logger as idaeslog
from idaes.core import (ControlVolume0DBlock,
                        declare_process_block_class,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        UnitModelBlockData,
                        useDefault)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.testing import get_default_solver
import idaes.core.util.scaling as iscale

_log = idaeslog.getLogger(__name__)

@declare_process_block_class("PressureExchanger")
class PressureExchangerData(UnitModelBlockData):
    """
    Standard Pressure Exchanger Unit Model Class:
    - steady state only
    """
    CONFIG = ConfigBlock()

    CONFIG.declare("dynamic", ConfigValue(
        domain=In([False]),
        default=False,
        description="Dynamic model flag - must be False",
        doc="""Indicates whether this model will be dynamic or not,
    **default** = False. Pressure exchangers do not support dynamic behavior."""))
    CONFIG.declare("has_holdup", ConfigValue(
        default=False,
        domain=In([False]),
        description="Holdup construction flag - must be False",
        doc="""Indicates whether holdup terms should be constructed or not.
    **default** - False. Pressure exchangers do not have defined volume, thus
    this must be False."""))
    CONFIG.declare("material_balance_type", ConfigValue(
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
    **MaterialBalanceType.total** - use total material balance.}"""))
    CONFIG.declare("energy_balance_type", ConfigValue(
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
    **EnergyBalanceType.energyPhase** - energy balances for each phase.}"""))
    CONFIG.declare("momentum_balance_type", ConfigValue(
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
    **MomentumBalanceType.momentumPhase** - momentum balances for each phase.}"""))
    CONFIG.declare("is_isothermal", ConfigValue(
        default=True,
        domain=In([True, False]),
        description="Pressure change term construction flag",
        doc="""Indicates whether terms for pressure change should be
        constructed,
        **default** - False.
        **Valid values:** {
        **True** - include pressure change terms,
        **False** - exclude pressure change terms.}"""))
    CONFIG.declare("property_package", ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for control volume",
        doc="""Property parameter object used to define property calculations,
        **default** - useDefault.
        **Valid values:** {
        **useDefault** - use default package from parent model or flowsheet,
        **PhysicalParameterObject** - a PhysicalParameterBlock object.}"""))
    CONFIG.declare("property_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing property packages",
        doc="""A ConfigBlock with arguments to be passed to a property block(s)
        and used when constructing these,
        **default** - None.
        **Valid values:** {
        see property package for documentation.}"""))

    def build(self):
        super().build()

        # Pressure exchanger supports only liquid phase
        if (len(self.config.property_package.phase_list) > 1
                or 'Liq' not in [p for p in self.config.property_package.phase_list]):
            raise ConfigurationError(
                "RO model only supports one liquid phase ['Liq'],"
                "the property package has specified the following phases {}"
                    .format([p for p in self.config.property_package.phase_list]))

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        units_meta = self.config.property_package.get_metadata().get_derived_units

        self.efficiency_pressure_exchanger = Var(
            self.flowsheet().config.time,
            initialize=0.9,
            bounds=(1e-6, 1),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc='Pressure exchanger efficiency')
        self.work_transfer = Var(
            self.flowsheet().config.time,
            initialize=1,
            bounds=(-1e6, 1e6),
            domain=Reals,
            units=pyunits.dimensionless,
            doc='Work transferred to low pressure side')

        # Build control volume for high pressure side
        self.high_pressure_side = ControlVolume0DBlock(default={
            "dynamic": False,
            "has_holdup": False,
            "property_package": self.config.property_package,
            "property_package_args": self.config.property_package_args})

        self.high_pressure_side.add_state_blocks(
            has_phase_equilibrium=False)

        self.high_pressure_side.add_material_balances(
            balance_type=self.config.material_balance_type)

        if self.config.is_isothermal:
            @self.high_pressure_side.Constraint(
                self.flowsheet().config.time,
                doc="Isothermal constraint")
            def isothermal_temperature(b, t):
                return b.properties_in[t].temperature == b.properties_out[t].temperature
        else:
            self.high_pressure_side.add_energy_balances(
                balance_type=self.config.energy_balance_type)

        if not self.config.is_isothermal:
            self.high_pressure_side.add_energy_balances(
                balance_type=self.config.energy_balance_type)

        self.high_pressure_side.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=True)

        self.high_pressure_side.work_fluid = Var(
            self.flowsheet().config.time,
            initialize=1.0,
            bounds=(-1e6, 1e6),
            domain=Reals,
            doc='Work required to increase the pressure of the liquid',
            units=units_meta('power'))

        # Build control volume for high pressure side
        self.low_pressure_side = ControlVolume0DBlock(default={
            "dynamic": False,
            "has_holdup": False,
            "property_package": self.config.property_package,
            "property_package_args": self.config.property_package_args})

        self.low_pressure_side.add_state_blocks(
            has_phase_equilibrium=False)

        self.low_pressure_side.add_material_balances(
            balance_type=self.config.material_balance_type)

        if self.config.is_isothermal:
            @self.low_pressure_side.Constraint(
                self.flowsheet().config.time,
                doc="Isothermal constraint")
            def isothermal_temperature(b, t):
                return b.properties_in[t].temperature == b.properties_out[t].temperature
        else:
            self.low_pressure_side.add_energy_balances(
                balance_type=self.config.energy_balance_type)

        self.low_pressure_side.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=True)

        self.low_pressure_side.work_fluid = Var(
            self.flowsheet().config.time,
            initialize=1.0,
            bounds=(-1e6, 1e6),
            domain=Reals,
            doc="Work required to increase the pressure of the liquid",
            units=units_meta('power'))

        # Add Ports
        self.add_inlet_port(name='high_pressure_inlet', block=self.high_pressure_side)
        self.add_outlet_port(name='high_pressure_outlet', block=self.high_pressure_side)
        self.add_inlet_port(name='low_pressure_inlet', block=self.low_pressure_side)
        self.add_outlet_port(name='low_pressure_outlet', block=self.low_pressure_side)

        # Performance equation
        @self.high_pressure_side.Constraint(
            self.flowsheet().config.time,
            doc="Fluid work term")
        def eq_work_fluid(b, t):
            return b.work_fluid[t] == b.properties_out[t].flow_vol * b.deltaP[t]

        @self.low_pressure_side.Constraint(
            self.flowsheet().config.time,
            doc="Fluid work term")
        def eq_work_fluid(b, t):
            return b.work_fluid[t] == b.properties_out[t].flow_vol * b.deltaP[t]

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Work transfer from high pressure side")
        def eq_work_transfer_high_pressure(b, t):
            return (b.work_transfer[t] ==
                    b.efficiency_pressure_exchanger[t]
                    * -b.high_pressure_side.work_fluid[t])

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Work transfer to low pressure side")
        def eq_work_transfer_low_pressure(b, t):
            return b.work_transfer[t] == b.low_pressure_side.work_fluid[t]

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Equal volumetric flow rate constraint")
        def eq_equal_flow_vol(b, t):
            return (b.high_pressure_side.properties_in[t].flow_vol ==
                    b.low_pressure_side.properties_in[t].flow_vol)
    #
    # def initialize(
    #         self,
    #         state_args=None,
    #         routine=None,
    #         outlvl=idaeslog.NOTSET,
    #         solver=None,
    #         optarg={"tol": 1e-6}):
    #     """
    #     General wrapper for pressure changer initialization routines
    #     Keyword Arguments:
    #         routine : str stating which initialization routine to execute
    #                     * None - currently no specialized routine for RO unit
    #         state_args : a dict of arguments to be passed to the property
    #                      package(s) to provide an initial state for
    #                      initialization (see documentation of the specific
    #                      property package) (default = {}).
    #         outlvl : sets output level of initialization routine
    #         optarg : solver options dictionary object (default={'tol': 1e-6})
    #         solver : str indicating whcih solver to use during
    #                  initialization (default = 'ipopt')
    #     Returns:
    #         None
    #     """
    #     # Get loggers
    #     init_log = idaeslog.getInitLogger(self.name, outlvl, tag="properties")
    #     solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="properties")
    #
    #     # Set solver and options
    #     # TODO: clean up once IDAES new API for initialize solvers is released
    #     if isinstance(solver, str):
    #         opt = SolverFactory(solver)
    #         opt.options = optarg
    #     else:
    #         if solver is None:
    #             opt = get_default_solver()
    #         else:
    #             opt = solver
    #             opt.options = optarg
    #
    #     # Solve simple unit
    #     with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
    #         res = opt.solve(self, tee=slc.tee)
    #     init_log.info_high(
    #         "Initialized zero order separator.".format(idaeslog.condition(res)))
    #
    # def calculate_scaling_factors(self):
    #     super().calculate_scaling_factors()
    #
    #     # scale variables
    #     for (t, p, j), v in self.recovery_frac_phase_comp.items():
    #         recovery_var = v
    #         removal_var = self.removal_frac_phase_comp[t, p, j]
    #         if (iscale.get_scaling_factor(recovery_var) is None
    #                 and iscale.get_scaling_factor(removal_var) is None):
    #             sf = iscale.get_scaling_factor(removal_var, default=1, warning=True)
    #             iscale.set_scaling_factor(recovery_var, sf)
    #             iscale.set_scaling_factor(removal_var, sf)
    #         elif iscale.get_scaling_factor(recovery_var) is None:
    #             sf = 1 / iscale.get_scaling_factor(removal_var)
    #             iscale.set_scaling_factor(recovery_var, sf)
    #         elif iscale.get_scaling_factor(removal_var) is None:
    #             sf = 1 / iscale.get_scaling_factor(recovery_var)
    #             iscale.set_scaling_factor(removal_var, sf)
    #
    #     for t, v in self.deltaP_outlet.items():
    #         if iscale.get_scaling_factor(v) is None:
    #             sf = iscale.get_scaling_factor(self.inlet_state[t].pressure)
    #             sf = sf / 10
    #             iscale.set_scaling_factor(v, sf)
    #
    #     for t, v in self.deltaP_waste.items():
    #         if iscale.get_scaling_factor(v) is None:
    #             sf = iscale.get_scaling_factor(self.inlet_state[t].pressure)
    #             sf = sf / 10
    #             iscale.set_scaling_factor(v, sf)
    #
    #     # transform constraints
    #     for (t, p, j), c in self.eq_component_mass_balance.items():
    #         sf = iscale.get_scaling_factor(self.inlet_state[t].get_material_flow_terms(p, j))
    #         iscale.constraint_scaling_transform(c, sf)
    #
    #     for (t, p, j), c in self.eq_component_removal.items():
    #         sf = iscale.get_scaling_factor(self.inlet_state[t].get_material_flow_terms(p, j))
    #         iscale.constraint_scaling_transform(c, sf)
    #
    #     for (t, p, j), c in self.eq_removal_to_recovery.items():
    #         sf = iscale.get_scaling_factor(self.removal_frac_phase_comp[t, p, j])
    #         iscale.constraint_scaling_transform(c, sf)
    #
    #     for t, c in self.eq_outlet_temperature.items():
    #         sf = iscale.get_scaling_factor(self.inlet_state[t].temperature)
    #         iscale.constraint_scaling_transform(c, sf)
    #
    #     for t, c in self.eq_waste_temperature.items():
    #         sf = iscale.get_scaling_factor(self.inlet_state[t].temperature)
    #         iscale.constraint_scaling_transform(c, sf)
    #
    #     for t, c in self.eq_outlet_pressure.items():
    #         sf = iscale.get_scaling_factor(self.inlet_state[t].pressure)
    #         iscale.constraint_scaling_transform(c, sf)
    #
    #     for t, c in self.eq_waste_pressure.items():
    #         sf = iscale.get_scaling_factor(self.inlet_state[t].pressure)
    #         iscale.constraint_scaling_transform(c, sf)
