###############################################################################
# ProteusLib Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/nawi-hub/proteuslib/"
#
###############################################################################

from pyomo.environ import (
    Block, ConcreteModel, Constraint, Expression, Var, Param, value, TransformationFactory, units as pyunits)
import proteuslib.flowsheets.full_treatment_train.flowsheet_components.financials as financials
from proteuslib.flowsheets.full_treatment_train.flowsheet_components import feed_block
from proteuslib.flowsheets.full_treatment_train.model_components import unit_separator, unit_0DRO, unit_1DRO, property_models

from proteuslib.flowsheets.full_treatment_train.flowsheet_components.desalination import (build_desalination,
                                                                                          solve_desalination,
                                                                                          scale_desalination,
                                                                                          initialize_desalination,
                                                                                          display_desalination)
# from proteuslib.flowsheets.full_treatment_train.analysis.flowsheet_softening_two_stage import build, scale, initialize
from proteuslib.flowsheets.full_treatment_train.util import (solve_with_user_scaling,
                                                             check_dof)
import idaes.core.util.scaling as iscale
from idaes.core import FlowsheetBlock


cost_capacity_flag = False

def build_costing(m, module=financials, **kwargs):
    '''
    Add costing to a given flowsheet. This function will
        1) call the get_costing method for each unit model (note: unit model must have a get_costing method
        to be detected), and
        2) call get_system_costing which will tally up all capex and opex for each process
    m : model
    module : financials module
    '''
    crf = m.fs.costing_param.factor_capital_annualization

    # call get_costing for each unit model
    #TODO: add in other components as they become available
    m.fs.costing = Block()

    # Nanofiltration
    if hasattr(m.fs, 'NF'):
        if kwargs['NF_type'] == 'ZO':
            m.fs.NF.get_costing(module=module, section='pretreatment')
        elif kwargs['NF_type'] == 'Sep':
            raise NotImplementedError("get_costing will not be implemented for the NF separator model.")
    if hasattr(m.fs, 'pump_NF'):
        m.fs.pump_NF.get_costing(module=module, section='pretreatment', pump_type="low pressure", cost_capacity=cost_capacity_flag)

    # Reverse Osmosis
    if hasattr(m.fs, 'RO'):
        if kwargs['RO_type'] == '0D' or kwargs['RO_type'] == '1D':
            m.fs.RO.get_costing(module=module, section='primary')
        elif kwargs['RO_type'] == 'Sep':
            raise NotImplementedError("get_costing will not be implemented for the RO separator model.")

    # Stage 2 RO
    if hasattr(m.fs, 'RO2'):
        m.fs.RO2.get_costing(module=module, RO_type='high pressure', section='primary')

    # Pump
    if hasattr(m.fs, 'pump_RO'):
        m.fs.pump_RO.get_costing(module=module, section='primary', pump_type="high pressure")

    # Stage 2 pump
    if hasattr(m.fs, 'pump_RO2'):
        m.fs.pump_RO2.get_costing(module=module, section='primary', pump_type="high pressure")

    # ERD
    if hasattr(m.fs, 'ERD'):
        m.fs.ERD.get_costing(module=module, section='primary', pump_type='pressure exchanger')

    # Pretreatment
    if hasattr(m.fs,'stoich_softening_mixer_unit'):
        m.fs.stoich_softening_mixer_unit.get_costing(module=module, section='pretreatment', mixer_type="lime_softening",cost_capacity=cost_capacity_flag)
        m.fs.lime_softening_unit_capex = Expression(expr=m.fs.stoich_softening_mixer_unit.costing.capital_cost/m.fs.annual_water_production *crf)
        m.fs.lime_softening_unit_opex = Expression(expr=m.fs.stoich_softening_mixer_unit.costing.operating_cost / m.fs.annual_water_production)
    else:
        m.fs.lime_softening_unit_capex = Expression(expr=0)
        m.fs.lime_softening_unit_opex = Expression(expr=0)

    # Post-treatment
    if hasattr(m.fs,'ideal_naocl_mixer_unit'):
        # print('FOUND CHLORINATION UNIT')
        m.fs.ideal_naocl_mixer_unit.get_costing(module=module, section='post_treatment', mixer_type='naocl_mixer', cost_capacity=cost_capacity_flag)
        m.fs.chlorination_unit_capex = Expression(expr=m.fs.ideal_naocl_mixer_unit.costing.capital_cost/m.fs.annual_water_production *crf)
        m.fs.chlorination_unit_opex = Expression(expr=m.fs.ideal_naocl_mixer_unit.costing.operating_cost / m.fs.annual_water_production)
    else:
        m.fs.chlorination_unit_capex = Expression(expr=0)
        m.fs.chlorination_unit_opex = Expression(expr=0)

    if hasattr(m.fs,'mixer_permeate'):
        m.fs.mixer_permeate.get_costing(module=module, section='primary', cost_capacity=cost_capacity_flag)

    # call get_system_costing for whole flowsheet
    module.get_system_costing(m.fs)

    # apply scaling to cost variables and constraints
    scale_costing(m)

def scale_costing(self):
    for b_unit in self.component_objects(Block, descend_into=True):
        if hasattr(b_unit, 'costing'):
            base = b_unit.costing
            for var in base.component_objects(Var):
                if iscale.get_scaling_factor(var) is None:
                    iscale.set_scaling_factor(var, 1e-3)
                for con in base.component_objects(Constraint):
                    sf = iscale.get_scaling_factor(var)
                    iscale.constraint_scaling_transform(con, sf)

def display_costing(m):
    crf = m.fs.costing_param.factor_capital_annualization
    if not hasattr(m.fs, 'pump_RO2'):
        m.fs.pump_RO2 = Block()
        m.fs.pump_RO2.costing = Block()
        m.fs.pump_RO2.costing.operating_cost = Param(initialize=0)
    if not hasattr(m.fs, 'NF'):
        m.fs.NF = Block()
        m.fs.NF.costing = Block()
        m.fs.NF.costing.operating_cost = Param(initialize=0)
    if not hasattr(m.fs, 'RO2'):
        m.fs.RO2 = Block()
        m.fs.RO2.costing = Block()
        m.fs.RO2.costing.operating_cost = Param(initialize=0)

    # UNITS FOR ALL COST COMPONENTS [=] $/m3 of permeate water produced
    cost_dict={
        'LCOW': m.fs.costing.LCOW, # Total LCOW
        'Total CAPEX': m.fs.costing.investment_cost_total * crf
                       / m.fs.annual_water_production,  # Direct + Indirect CAPEX
        'Direct CAPEX': m.fs.costing.capital_cost_total * crf
                        / m.fs.annual_water_production,  # Direct CAPEX for all system components
        'Indirect CAPEX': (m.fs.costing.investment_cost_total - m.fs.costing.capital_cost_total) * crf
                        / m.fs.annual_water_production,  # Indirect CAPEX for miscellaneous items
        'Total OPEX': m.fs.costing.operating_cost_total / m.fs.annual_water_production,  # Total OPEX
        'Labor & Maintenance Costs': m.fs.costing.operating_cost_labor_maintenance
                                     / m.fs.annual_water_production,
        'Total Electricity Cost': m.fs.costing.electricity_cost_total / m.fs.annual_water_production,
        'Total Membrane Replacement Cost': (m.fs.NF.costing.operating_cost
                                            + m.fs.RO.costing.operating_cost
                                            + m.fs.RO2.costing.operating_cost) / m.fs.annual_water_production,
        'Total Pretreatment Cost': m.fs.costing.pretreatment_cost_total / m.fs.annual_water_production,
        'Total Primary Cost': m.fs.costing.primary_cost_total / m.fs.annual_water_production,
        'Total Post-treatment Cost': m.fs.costing.post_treatment_cost_total / m.fs.annual_water_production,
        'Lime softener CAPEX': m.fs.lime_softening_unit_capex,
        'Lime softener OPEX': m.fs.lime_softening_unit_opex,
        'Chlorination CAPEX': m.fs.chlorination_unit_capex,
        'Chlorination OPEX': m.fs.chlorination_unit_opex,

    }

    for item, val in cost_dict.items():
        print(f"{item} = {value(val)}")

    return cost_dict
