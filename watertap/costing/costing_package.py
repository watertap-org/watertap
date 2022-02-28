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
from pyomo.environ import (
    Block, Constraint, Expression, Var, Param, Reals, NonNegativeReals, units as pyunits)
from idaes.core.util.exceptions import ConfigurationError


class CostingPackage():
    # Define standard material flows and costs
    defined_flows = {
        "electricity": (0.07, pyo.units.dollars / (pyo.units.kW * pyo.units.h))}

    # # Define default mapping of costing methods to unit models
    # default_mapping = {
    #     "Pump": pump_costing,  # TODO: this may be carrying over a limitation of IDAES
    #     "ReverseOsmosis0D": ro_costing}

    # Define units of measurement and conversions
    # TBD, based on pint's API

    @staticmethod
    def build_global_params(blk):
        # Build flowsheet level costing components
        # This is package specific

        # TODO: add units to parameters
        blk.load_factor = Var(
            initialize=0.9,
            doc='Load factor [fraction of uptime]')
        blk.factor_total_investment = Var(
            initialize=2,
            doc='Total investment factor [investment cost/equipment cost]')
        blk.factor_MLC = Var(
            initialize=0.03,
            doc='Maintenance-labor-chemical factor [fraction of investment cost/year]')
        blk.factor_capital_annualization = Var(
            initialize=0.1,
            doc='Capital annualization factor [fraction of investment cost/year]')
        blk.factor_membrane_replacement = Var(
            initialize=0.2,
            doc='Membrane replacement factor [fraction of membrane replaced/year]')
        blk.mem_cost = Var(
            initialize=30,
            doc='Membrane cost [$/m2]')
        blk.hp_pump_cost = Var(
            initialize=53 / 1e5 * 3600,
            doc='High pressure pump cost [$/W]')
        blk.pxr_cost = Var(
            initialize=535,
            doc='Pressure exchanger cost [$/(m3/h)]')

    @staticmethod
    def build_process_cost(blk):
        # blk = m.fs.costing
        # Build flowsheet level total costs
        blk.total_capital_cost = Var()
        blk.total_capital_cost_constraint = Constraint(
            expr=blk.total_capital_cost ==
            blk.aggregate_capital_cost * blk.factor_total_investment)

        blk.total_fixed_operating_cost = Var()
        blk.total_fixed_operating_cost_constraint = Constraint(
            expr=blk.total_fixed_operating_cost ==
                 blk.aggregate_fixed_operating_cost
                 + blk.total_capital_cost * blk.factor_MLC)

        blk.total_variable_operating_cost = Var()  # Unnecessary
        blk.total_variable_operating_cost_constraint = Constraint(
            expr=blk.total_variable_operating_cost ==
                 blk.aggregate_variable_operating_cost)


    # Define costing methods supported by package
    @staticmethod
    def cost_ro(blk):
        # blk = m.fs.unit_cost['fs.RO1']
        # Add code to create Expressions for desired Unit Level Costs
        blk.capital_cost = Var(initialize=1e5,
                               bounds=(0, 1e8),
                               domain=NonNegativeReals,
                               doc='Unit capital cost [$]')
        blk.TBD_capital_cost = Expression(
            expr=blk.params.mem_cost * blk.unit_model.area)

        blk.TBD_fixed_operating_cost = Expression(
            expr=blk.params.factor_membrane_replacement * blk.TBD_capital_cost)

    @staticmethod
    def cost_pump(blk):
        # Add code to create Expressions for desired Unit Level Costs
        blk.TBD_capital_cost = Expression(
            expr=blk.params.hp_pump_cost * blk.unit_model.work_mechanical[0])

# def add_costing_param_block(self):
#     self.costing_param = Block()
#     b = self.costing_param
#
#     b.load_factor = Var(
#         initialize=0.9,
#         doc='Load factor [fraction of uptime]')
#     b.factor_total_investment = Var(
#         initialize=2,
#         doc='Total investment factor [investment cost/equipment cost]')
#     b.factor_MLC = Var(
#         initialize=0.03,
#         doc='Maintenance-labor-chemical factor [fraction of investment cost/year]')
#     b.factor_capital_annualization = Var(
#         initialize=0.1,
#         doc='Capital annualization factor [fraction of investment cost/year]')
#     b.factor_membrane_replacement = Var(
#         initialize=0.2,
#         doc='Membrane replacement factor [fraction of membrane replaced/year]')
#     b.electricity_cost = Var(
#         initialize=0.07,
#         doc='Electricity cost [$/kWh]')
#     b.mem_cost = Var(
#         initialize=30,
#         doc='Membrane cost [$/m2]')
#     b.hp_pump_cost = Var(
#         initialize=53 / 1e5 * 3600,
#         doc='High pressure pump cost [$/W]')
#     b.pxr_cost = Var(
#         initialize=535,
#         doc='Pressure exchanger cost [$/(m3/h)]')
#
#     # traditional parameters are the only Vars on the block and should be fixed
#     for v in b.component_objects(Var, descend_into=True):
#         for i in v:
#             if v[i].value is None:
#                 raise ConfigurationError(
#                     "{} parameter {} was not assigned"
#                     " a value. Please check your configuration "
#                     "arguments.".format(b.name, v.local_name))
#             v[i].fix()
#
#
# def get_system_costing(self):
#     if not hasattr(self, 'costing'):
#         self.costing = Block()
#     b = self.costing
#
#     b.capital_cost_total = Var(
#         initialize=1e3,
#         domain=NonNegativeReals,
#         doc='Total capital cost [$]')
#     b.investment_cost_total = Var(
#         initialize=1e3,
#         domain=NonNegativeReals,
#         doc='Total investment cost [$]')
#     b.operating_cost_MLC = Var(
#         initialize=1e3,
#         domain=NonNegativeReals,
#         doc='Maintenance-labor-chemical operating cost [$/year]')
#     b.operating_cost_total = Var(
#         initialize=1e3,
#         domain=NonNegativeReals,
#         doc='Total operating cost [$/year]')
#     b.LCOW = Var(
#         initialize=1,
#         domain=NonNegativeReals,
#         doc='Levelized cost of water [$/m3]')
#
#     capital_cost_var_lst = []
#     operating_cost_var_lst = []
#     for b_unit in self.component_objects(Block, descend_into=True):
#         if hasattr(b_unit, 'costing'):
#             capital_cost_var_lst.append(b_unit.costing.capital_cost)
#             operating_cost_var_lst.append(b_unit.costing.operating_cost)
#     operating_cost_var_lst.append(b.operating_cost_MLC)
#
#     b.eq_capital_cost_total = Constraint(
#         expr=b.capital_cost_total == sum(capital_cost_var_lst))
#     b.eq_investment_cost_total = Constraint(
#         expr=(b.investment_cost_total ==
#               b.capital_cost_total * self.costing_param.factor_total_investment))
#     b.eq_operating_cost_MLC = Constraint(
#         expr=(b.operating_cost_MLC ==
#               b.investment_cost_total * self.costing_param.factor_MLC))
#     b.eq_operating_cost_total = Constraint(
#         expr=b.operating_cost_total == sum(operating_cost_var_lst))
#     b.eq_LCOW = Constraint(
#         expr=b.LCOW == (b.investment_cost_total * self.costing_param.factor_capital_annualization
#                         + b.operating_cost_total) / (self.annual_water_production / (pyunits.m**3/pyunits.year)))
#
#
# def _make_vars(self):
#     # build generic costing variables (all costing models need these vars)
#     self.capital_cost = Var(initialize=1e5,
#                             domain=NonNegativeReals,
#                             doc='Unit capital cost [$]')
#     self.operating_cost = Var(initialize=1e5,
#                               domain=Reals,
#                               bounds=(0, 1e6),
#                               doc='Unit operating cost [$/year]')
#
#
# def ReverseOsmosis_costing(self):
#     _make_vars(self)
#
#     b_RO = self.parent_block()
#     b_fs = b_RO.parent_block()
#
#     # capital cost
#     self.eq_capital_cost = Constraint(
#         expr=self.capital_cost == b_fs.costing_param.mem_cost * b_RO.area/pyunits.m**2)
#
#     # operating cost
#     self.eq_operating_cost = Constraint(
#         expr=self.operating_cost == b_fs.costing_param.factor_membrane_replacement
#              * b_fs.costing_param.mem_cost * b_RO.area / pyunits.m ** 2)
#
#
# def PressureExchanger_costing(self):
#     _make_vars(self)
#
#     b_PXR = self.parent_block()
#     b_fs = b_PXR.parent_block()
#
#     # capital cost
#     self.eq_capital_cost = Constraint(
#         expr=self.capital_cost == b_fs.costing_param.pxr_cost
#              * b_PXR.low_pressure_side.properties_in[0].flow_vol*3600 / (pyunits.m**3/pyunits.s))
#
#     # operating cost
#     self.operating_cost.fix(0)
#
# def pressure_changer_costing(self, pump_type="centrifugal"):
#     _make_vars(self)
#
#     b_PC = self.parent_block()
#     b_fs = b_PC.parent_block()
#
#     self.purchase_cost = Var()
#     self.cp_cost_eq = Constraint(expr=self.purchase_cost == 0)
#
#     if pump_type == 'High pressure':
#         # capital cost
#         self.eq_capital_cost = Constraint(
#             expr=self.capital_cost == b_fs.costing_param.hp_pump_cost * b_PC.work_mechanical[0] / pyunits.W)
#
#         # operating cost
#         self.eq_operating_cost = Constraint(
#             expr=self.operating_cost == (b_PC.work_mechanical[0] / pyunits.W
#                                          * 3600 * 24 * 365 * b_fs.costing_param.load_factor)
#             * b_fs.costing_param.electricity_cost / 3600 / 1000)
#
#     elif pump_type == 'Pressure exchanger':
#         # capital cost
#         b_cv_in = b_PC.control_volume.properties_in[0]
#         self.eq_capital_cost = Constraint(
#             expr=(self.capital_cost == b_fs.costing_param.erd_cost['A']
#                  * (sum(b_cv_in.flow_mass_comp[j] / (pyunits.kg/pyunits.s)
#                         for j in b_PC.config.property_package.component_list)
#                  / (b_cv_in.dens_mass / (pyunits.kg/pyunits.m**3)) * 3600) ** 0.58))
#
#         # operating cost
#         self.operating_cost.setlb(-1e6)
#         self.eq_operating_cost = Constraint(
#             expr=self.operating_cost == (b_PC.work_mechanical[0] / pyunits.W
#                                          * 3600 * 24 * 365 * b_fs.costing_param.load_factor)
#                  * b_fs.costing_param.electricity_cost / 3600 / 1000)
