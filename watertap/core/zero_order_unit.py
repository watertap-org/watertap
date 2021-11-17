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
"""
This module contains the base class for all zero order single inlet-two outlet
(SITO) unit models.
"""
from idaes.core import UnitModelBlockData, useDefault
from idaes.core.util.config import is_physical_parameter_block
import idaes.logger as idaeslog
from idaes.core.util import get_solver
import idaes.core.util.scaling as iscale
from idaes.core.util.exceptions import ConfigurationError

from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.environ import NonNegativeReals, Var, units as pyunits

# Some more inforation about this module
__author__ = "Andrew Lee"


class SITOBaseData(UnitModelBlockData):
    """
    Standard base class for single inlet-two outlet unit models.

    This class is intended to be used for creating derived model classes and
    cannot be instantiated as by itself. When creating derived classes,
    developers must set the '_has_deltaP_outlet` and `_has_deltaP_waste`
    attributes on the model before calling `super().build()`. These attributes
    determine whetehr the deltaP terms for the outlet and waste streams are
    added to the model and included in the pressure constraints.
    """

    CONFIG = ConfigBlock()
    CONFIG.declare('dynamic', ConfigValue(
        domain=In([False]),
        default=False,
        description='Dynamic model flag - must be False',
        doc='''All zero-order models are steady-state only'''))
    CONFIG.declare('has_holdup', ConfigValue(
        default=False,
        domain=In([False]),
        description='Holdup construction flag - must be False',
        doc='''Zero order models do not include holdup'''))
    CONFIG.declare('property_package', ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description='Property package to use for control volume',
        doc='''Property parameter object used to define property  calculations,
        **default** - useDefault.
        **Valid values:** {
        **useDefault** - use default package from parent model or flowsheet,
        **PhysicalParameterObject** - a PhysicalParameterBlock object.}'''))
    CONFIG.declare('property_package_args', ConfigBlock(
        implicit=True,
        description='Arguments to use for constructing property packages',
        doc='''A ConfigBlock with arguments to be passed to a property block(s)
        and used when constructing these, **default** - None.
        **Valid values:** {see property package for documentation.}'''))

    def build(self):
        super().build()

        # Check that derived class has implemented flags for pressure change
        if not hasattr(self, "_has_deltaP_outlet"):
            raise NotImplementedError(
                f"{self.name} derived class class has not implemented "
                f"_has_deltaP_outlet.")
        if not hasattr(self, "_has_deltaP_waste"):
            raise NotImplementedError(
                f"{self.name} derived class class has not implemented "
                f"_has_deltaP_waste.")

        # Check that property package meets requirements
        if self.config.property_package.phase_list != ["Liq"]:
            raise ConfigurationError(
                f"{self.name} configured with invalid property package. "
                f"Zero-order models only support proeprty packages with a "
                f"single phase named 'Liq'.")
        if (not hasattr(self.config.property_package, "solvent_set") or
                self.config.property_package.solvent_set != ["H2O"]):
            raise ConfigurationError(
                f"{self.name} configured with invalid property package. "
                f"Zero-order models only support property packages which "
                f"include 'H2O' as the only Solvent.")
        if not hasattr(self.config.property_package, "solute_set"):
            raise ConfigurationError(
                f"{self.name} configured with invalid property package. "
                f"Zero-order models require property packages to declare all "
                f"dissolved species as Solutes.")
        if (len(self.config.property_package.solute_set) !=
                len(self.config.property_package.component_list)-1):
            raise ConfigurationError(
                f"{self.name} configured with invalid property package. "
                f"Zero-order models only support `H2O` as a solvent and all "
                f"other species as Solutes.")

        # Get units metadata
        units_meta = \
            self.config.property_package.get_metadata().get_derived_units

        # Create state blocks for inlet and outlets
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["defined_state"] = True

        self.properties_in = self.config.property_package.build_state_block(
            self.flowsheet().time,
            doc="Material properties at inlet",
            default=tmp_dict)

        tmp_dict_2 = dict(**tmp_dict)
        tmp_dict_2["defined_state"] = False

        self.properties_out = self.config.property_package.build_state_block(
            self.flowsheet().time,
            doc="Material properties at outlet",
            default=tmp_dict_2)
        self.properties_waste = self.config.property_package.build_state_block(
            self.flowsheet().time,
            doc="Material properties at waste outlet",
            default=tmp_dict_2)

        # Create Ports
        self.add_port("inlet", self.properties_in, doc="Inlet port")
        self.add_port("outlet", self.properties_out, doc="Outlet port")
        self.add_port("waste", self.properties_waste, doc="Waste port")

        # Add performance variables
        self.recovery_vol = Var(
            self.flowsheet().time,
            initialize=0.8,
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            bounds=(1E-8, 1.0000001),
            doc='Volumetric recovery fraction of water in the outlet stream')
        self.removal_mass_solute = Var(
            self.flowsheet().time,
            self.config.property_package.solute_set,
            domain=NonNegativeReals,
            initialize=0.01,
            units=pyunits.dimensionless,
            doc='Solute removal fraction on a mass basis')

        if self._has_deltaP_outlet:
            self.deltaP_outlet = Var(
                self.flowsheet().time,
                initialize=0,
                units=units_meta('pressure'),
                doc='Pressure change between inlet and outlet')
        if self._has_deltaP_waste:
            self.deltaP_waste = Var(
                self.flowsheet().time,
                initialize=0,
                units=units_meta('pressure'),
                doc='Pressure change between inlet and waste')

        # Add performance constraints
        # Water recovery
        @self.Constraint(self.flowsheet().time, doc='Water recovery equation')
        def water_recovery_equation(b, t):
            return (b.recovery_vol[t] * b.properties_in[t].flow_vol ==
                    b.properties_out[t].flow_vol)

        # Flow balance
        @self.Constraint(self.flowsheet().time, doc='Overall flow balance')
        def flow_balance(b, t):
            return (b.properties_in[t].flow_vol ==
                    b.properties_out[t].flow_vol +
                    b.properties_waste[t].flow_vol)

        # Solute removal
        @self.Constraint(self.flowsheet().time,
                         self.config.property_package.solute_set,
                         doc='Solute removal equations')
        def solute_removal_equation(b, t, j):
            return (b.removal_mass_solute[t, j] *
                    b.properties_in[t].get_material_flow_terms("Liq", j) ==
                    b.properties_waste[t].get_material_flow_terms("Liq", j))

        # Solute balances
        @self.Constraint(self.flowsheet().time,
                         self.config.property_package.solute_set,
                         doc='Solute mass balance equations')
        def solute_mass_balances(b, t, j):
            return (b.properties_in[t].get_material_flow_terms("Liq", j) ==
                    b.properties_out[t].get_material_flow_terms("Liq", j) +
                    b.properties_waste[t].get_material_flow_terms("Liq", j))

        # Pressure drop
        @self.Constraint(self.flowsheet().time, doc='Outlet pressure equation')
        def outlet_pressure_constraint(b, t):
            if self._has_deltaP_outlet:
                dp = b.deltaP_outlet[t]
            else:
                dp = 0
            return (b.properties_in[t].pressure + dp ==
                    b.properties_out[t].pressure)

        @self.Constraint(self.flowsheet().time, doc='Waste pressure equation')
        def waste_pressure_constraint(b, t):
            if self._has_deltaP_waste:
                dp = b.deltaP_waste[t]
            else:
                dp = 0
            return (b.properties_in[t].pressure + dp ==
                    b.properties_waste[t].pressure)

        # Temperature equality
        @self.Constraint(self.flowsheet().time,
                         doc='Outlet temperature equality')
        def outlet_temperature_equality(b, t):
            return (b.properties_in[t].temperature ==
                    b.properties_out[t].temperature)

        @self.Constraint(self.flowsheet().time,
                         doc='Waste temperature equality')
        def waste_temperature_equality(b, t):
            return (b.properties_in[t].temperature ==
                    b.properties_waste[t].temperature)

    def initialize(blk, state_args=None, outlvl=idaeslog.NOTSET,
                   solver=None, optarg=None):
        '''
        Initialization routine for single inlet-two outlet unit models.

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                           package(s) to provide an initial state for
                           initialization (see documentation of the specific
                           property package) (default = {}).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None, use default IDAES solver)

        Returns:
            None
        '''
        if optarg is None:
            optarg = {}

        # Set solver options
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

        solver_obj = get_solver(solver, optarg)

        # Get initial guesses for inlet if none provided
        if state_args is None:
            state_args = {}
            state_dict = (
                blk.properties_in[
                    blk.flowsheet().time.first()]
                .define_port_members())

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        # ---------------------------------------------------------------------
        # Initialize control volume block
        flags = blk.properties_in.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
            hold_state=True
        )
        blk.properties_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
            hold_state=False
        )
        blk.properties_waste.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
            hold_state=False
        )

        init_log.info_high('Initialization Step 1 Complete.')

        # ---------------------------------------------------------------------
        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = solver_obj.solve(blk, tee=slc.tee)

        init_log.info_high(
            "Initialization Step 2 {}.".format(idaeslog.condition(results))
        )

        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.properties_in.release_state(flags, outlvl)

        init_log.info('Initialization Complete: {}'
                      .format(idaeslog.condition(results)))

    def calculate_scaling_factors(self):
        # Get default scale factors and do calculations from base classes
        super().calculate_scaling_factors()

        for t, v in self.water_recovery_equation.items():
            iscale.constraint_scaling_transform(
                v, iscale.get_scaling_factor(
                    self.properties_in[t].flow_vol,
                    default=1,
                    warning=True,
                    hint=" for water recovery"))

        for t, v in self.flow_balance.items():
            iscale.constraint_scaling_transform(
                v, iscale.get_scaling_factor(
                    self.properties_in[t].flow_vol,
                    default=1,
                    warning=False))  # would just be a duplicate of above

        for (t, j), v in self.solute_removal_equation.items():
            iscale.constraint_scaling_transform(
                v, iscale.get_scaling_factor(
                    self.properties_in[t].flow_mass_comp[j],
                    default=1,
                    warning=True,
                    hint=" for solute removal"))

        for (t, j), v in self.solute_mass_balances.items():
            iscale.constraint_scaling_transform(
                v, iscale.get_scaling_factor(
                    self.properties_in[t].flow_mass_comp[j],
                    default=1,
                    warning=False))  # would just be a duplicate of above

        for t, v in self.outlet_pressure_constraint.items():
            iscale.constraint_scaling_transform(
                v, iscale.get_scaling_factor(
                    self.properties_in[t].pressure,
                    default=1e-5,
                    warning=True,
                    hint=" for outlet pressure constraint"))

        for t, v in self.waste_pressure_constraint.items():
            iscale.constraint_scaling_transform(
                v, iscale.get_scaling_factor(
                    self.properties_in[t].pressure,
                    default=1e-5,
                    warning=False))  # would just be a duplicate of above

        for t, v in self.outlet_temperature_equality.items():
            iscale.constraint_scaling_transform(
                v, iscale.get_scaling_factor(
                    self.properties_in[t].temperature,
                    default=1e-2,
                    warning=True,
                    hint=" for outlet temperature equality"))

        for t, v in self.waste_temperature_equality.items():
            iscale.constraint_scaling_transform(
                v, iscale.get_scaling_factor(
                    self.properties_in[t].temperature,
                    default=1e-2,
                    warning=False))  # would just be a duplicate of above
