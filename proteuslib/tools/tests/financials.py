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
from pyomo.environ import (
    Block, Constraint, Expression, Var, Param, NonNegativeReals, units as pyunits)
from idaes.core.util.exceptions import ConfigurationError


def add_costing_param_block(self):
    self.costing_param = Block()
    b = self.costing_param

    b.load_factor = Var(
        initialize=0.9,
        doc='Load factor [fraction of uptime]')
    b.factor_total_investment = Var(
        initialize=2,
        doc='Total investment factor [investment cost/equipment cost]')
    b.factor_MLC = Var(
        initialize=0.03,
        doc='Maintenance-labor-chemical factor [fraction of investment cost/year]')
    b.factor_capital_annualization = Var(
        initialize=0.1,
        doc='Capital annualization factor [fraction of investment cost/year]')
    b.factor_membrane_replacement = Var(
        initialize=0.2,
        doc='Membrane replacement factor [fraction of membrane replaced/year]')
    b.electricity_cost = Var(
        initialize=0.07,
        doc='Electricity cost [$/kWh]')
    b.mem_cost = Var(
        initialize=30,
        doc='Membrane cost [$/m2]')
    b.hp_pump_cost = Var(
        initialize=53 / 1e5 * 3600,
        doc='High pressure pump cost [$/W]')
    b.erd_cost = Var(
        ['A', 'B'],
        initialize={'A': 3134.7, 'B': 0.58},
        doc='Energy recovery device cost parameters')

    # traditional parameters are the only Vars on the block and should be fixed
    for v in b.component_objects(Var, descend_into=True):
        for i in v:
            if v[i].value is None:
                raise ConfigurationError(
                    "{} parameter {} was not assigned"
                    " a value. Please check your configuration "
                    "arguments.".format(b.name, v.local_name))
            v[i].fix()


def get_system_costing(self):
    if not hasattr(self, 'costing'):
        self.costing = Block()
    b = self.costing

    b.capital_cost_total = Var(
        initialize=1e6,
        domain=NonNegativeReals,
        doc='Total capital cost [$]')
    b.investment_cost_total = Var(
        initialize=1e6,
        domain=NonNegativeReals,
        doc='Total investment cost [$]')
    b.operating_cost_MLC = Var(
        initialize=1e5,
        domain=NonNegativeReals,
        doc='Maintenance-labor-chemical operating cost [$/year]')
    b.operating_cost_total = Var(
        initialize=1e5,
        domain=NonNegativeReals,
        doc='Total operating cost [$/year]')
    b.LCOW = Var(
        initialize=1e5,
        domain=NonNegativeReals,
        doc='Levelized cost of water [$/m3]')

    capital_cost_var_lst = []
    operating_cost_var_lst = []
    for b_unit in self.component_objects(Block, descend_into=True):
        if hasattr(b_unit, 'costing'):
            capital_cost_var_lst.append(b_unit.costing.capital_cost)
            operating_cost_var_lst.append(b_unit.costing.operating_cost)
    operating_cost_var_lst.append(b.operating_cost_MLC)

    b.eq_capital_cost_total = Constraint(
        expr=b.capital_cost_total == sum(capital_cost_var_lst))
    b.eq_investment_cost_total = Constraint(
        expr=(b.investment_cost_total ==
              b.capital_cost_total * self.costing_param.factor_total_investment))
    b.eq_operating_cost_MLC = Constraint(
        expr=(b.operating_cost_MLC ==
              b.investment_cost_total * self.costing_param.factor_MLC))
    b.eq_operating_cost_total = Constraint(
        expr=b.operating_cost_total == sum(operating_cost_var_lst))
    b.eq_LCOW = Constraint(
        expr=b.LCOW == (b.investment_cost_total * self.costing_param.factor_capital_annualization
                        + b.operating_cost_total) / (self.AWP / (pyunits.m**3/pyunits.year)))


def _make_vars(self):
    # build generic costing variables (all costing models need these vars)
    # self.base_cost = Var(initialize=1e5,
    #                      domain=NonNegativeReals,
    #                      doc='Unit base cost [$]')
    self.capital_cost = Var(initialize=1e5,
                            domain=NonNegativeReals,
                            doc='Unit capital cost [$]')
    self.operating_cost = Var(initialize=1e5,
                              domain=NonNegativeReals,
                              doc='Operating cost [$/year]')


def RO_costing(self):
    _make_vars(self)

    b_RO = self.parent_block()
    b_fs = b_RO.parent_block()

    # capital cost
    self.eq_capital_cost = Constraint(
        expr=self.capital_cost == b_fs.costing_param.mem_cost * b_RO.area/pyunits.m**2)

    # operating cost
    self.eq_operating_cost = Constraint(
        expr=self.operating_cost == b_fs.costing_param.factor_membrane_replacement
             * b_fs.costing_param.mem_cost * b_RO.area / pyunits.m ** 2)

def pressure_changer_costing(self, Mat_factor="stain_steel",
                             # applies for all (pump, compressor, fan, blower)
                             mover_type="compressor",
                             # fan, blower, compressor
                             compressor_type="centrifugal",
                             # only for compressor
                             driver_mover_type="electrical_motor",
                             # only for compressors
                             pump_type="centrifugal",
                             # centrifugal, external_gear, reciprocating
                             pump_type_factor='1.4',
                             # 1.1 to 1.4, 2.1 and 2.2
                             # (needs to be wise-selected by user see table)
                             pump_motor_type_factor='open',
                             # centrifugal_backward, centrifugal_straight
                             # vane_axial, tube_axial
                             fan_type='centrifugal_backward',
                             # select from table depends on fan's head
                             fan_head_factor=1.45,
                             # centrifugal and rotary
                             blower_type='centrifugal'):
    _make_vars(self)

    b_PC = self.parent_block()
    b_fs = b_PC.parent_block()

    if pump_type == 'High pressure':
        # capital cost
        self.eq_capital_cost = Constraint(
            expr=self.capital_cost == b_fs.costing_param.hp_pump_cost * b_PC.work_mechanical[0] / pyunits.W)

        # operating cost
        self.eq_operating_cost = Constraint(
            expr=self.operating_cost == (b_PC.work_mechanical[0] / pyunits.W
                                         * 3600 * 24 * 365 * b_fs.costing_param.load_factor)
            * b_fs.costing_param.electricity_cost / 3600 / 1000)

    elif pump_type == 'Pressure exchanger':
        # capital cost
        b_cv_in = b_PC.control_volume.properties_in[0]
        self.eq_capital_cost = Constraint(
            expr=(self.capital_cost == b_fs.costing_param.erd_cost['A']
                 * (sum(b_cv_in.flow_mass_comp[j] / (pyunits.kg/pyunits.s)
                        for j in b_PC.config.property_package.component_list)
                 / (b_cv_in.dens_mass / (pyunits.kg/pyunits.m**3)) * 3600) ** 0.58))

        # operating cost
        self.operating_cost.fix(0)

