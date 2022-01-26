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
This module contains the base class for all zero order single inlet-double
outlet (SIDO) unit models.
"""

import idaes.logger as idaeslog
from idaes.core.util import get_solver
import idaes.core.util.scaling as iscale
from idaes.core.util.tables import create_stream_table_dataframe

from pyomo.environ import NonNegativeReals, Var, units as pyunits

# Some more inforation about this module
__author__ = "Andrew Lee"

# Set up logger
_log = idaeslog.getLogger(__name__)


def build_sido(self):
    self._has_recovery_removal = True
    self._initialize = initialize_sido
    self._scaling = calculate_scaling_factors_sido

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

    self.properties_treated = \
        self.config.property_package.build_state_block(
            self.flowsheet().time,
            doc="Material properties of treated water",
            default=tmp_dict_2)
    self.properties_byproduct = \
        self.config.property_package.build_state_block(
            self.flowsheet().time,
            doc="Material properties of byproduct stream",
            default=tmp_dict_2)

    # Create Ports
    self.add_port("inlet", self.properties_in, doc="Inlet port")
    self.add_port("treated",
                  self.properties_treated,
                  doc="Treated water outlet port")
    self.add_port("byproduct",
                  self.properties_byproduct,
                  doc="Byproduct outlet port")

    # Add performance variables
    self.recovery_vol = Var(
        self.flowsheet().time,
        initialize=0.8,
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        bounds=(1E-8, 1.0000001),
        doc='Volumetric recovery fraction of water in the treated stream')
    self.removal_frac_mass_solute = Var(
        self.flowsheet().time,
        self.config.property_package.solute_set,
        domain=NonNegativeReals,
        initialize=0.01,
        units=pyunits.dimensionless,
        doc='Solute removal fraction on a mass basis')

    # Add performance constraints
    # Water recovery
    @self.Constraint(self.flowsheet().time, doc='Water recovery equation')
    def water_recovery_equation(b, t):
        return (b.recovery_vol[t] * b.properties_in[t].flow_vol ==
                b.properties_treated[t].flow_vol)

    # Flow balance
    @self.Constraint(self.flowsheet().time, doc='Overall flow balance')
    def flow_balance(b, t):
        return (b.properties_in[t].flow_vol ==
                b.properties_treated[t].flow_vol +
                b.properties_byproduct[t].flow_vol)

    # Solute removal
    @self.Constraint(self.flowsheet().time,
                     self.config.property_package.solute_set,
                     doc='Solute removal equations')
    def solute_removal_equation(b, t, j):
        return (b.removal_frac_mass_solute[t, j] *
                b.properties_in[t].conc_mass_comp[j] ==
                (1 - b.recovery_vol[t]) *
                b.properties_byproduct[t].conc_mass_comp[j])

    # Solute concentration of treated stream
    @self.Constraint(self.flowsheet().time,
                     self.config.property_package.solute_set,
                     doc='Constraint for solute concentration in treated '
                     'stream.')
    def solute_treated_equation(b, t, j):
        return ((1 - b.removal_frac_mass_solute[t, j]) *
                b.properties_in[t].conc_mass_comp[j] ==
                b.recovery_vol[t] *
                b.properties_treated[t].conc_mass_comp[j])

    # Add electricity consumption to model
    self.electricity = Var(self.flowsheet().time,
                           units=pyunits.kW,
                           doc="Electricity consumption of unit")
    self.energy_electric_flow_vol_inlet = Var(
        units=pyunits.kWh/pyunits.m**3,
        doc="Electricity intensity with respect to inlet flowrate of unit")

    @self.Constraint(self.flowsheet().time,
                     doc='Constraint for electricity consumption base on '
                     'feed flowrate.')
    def electricity_consumption(b, t):
        return b.electricity[t] == (
            b.energy_electric_flow_vol_inlet *
            pyunits.convert(b.properties_in[t].flow_vol,
                            to_units=pyunits.m**3/pyunits.hour))

    self._stream_table_dict = {"Inlet": self.inlet,
                               "Treated": self.treated,
                               "Byproduct": self.byproduct}

    self._perf_var_dict = {
        "Water Recovery": self.recovery_vol,
        "Solute Removal": self.removal_frac_mass_solute,
        "Electricity Demand": self.electricity,
        "Electricity Intensity": self.energy_electric_flow_vol_inlet}


def initialize_sido(blk, state_args=None, outlvl=idaeslog.NOTSET,
                    solver=None, optarg=None):
    '''
    Initialization routine for single inlet-double outlet unit models.

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
    blk.properties_treated.initialize(
        outlvl=outlvl,
        optarg=optarg,
        solver=solver,
        state_args=state_args,
        hold_state=False
    )
    blk.properties_byproduct.initialize(
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


def calculate_scaling_factors_sido(self):
    # Get default scale factors and do calculations from base classes
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

    for (t, j), v in self.solute_treated_equation.items():
        iscale.constraint_scaling_transform(
            v, iscale.get_scaling_factor(
                self.properties_in[t].flow_mass_comp[j],
                default=1,
                warning=False))  # would just be a duplicate of above
