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
    Block, Constraint, Expression, Var, Param, units as pyunits, value)
import proteuslib.flowsheets.full_treatment_train.example_flowsheets.financials as financials
import idaes.core.util.scaling as iscale

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
    #get_costing_sweep(m.fs, module=financials)
    #TODO: add in other components as they become available
    m.fs.costing = Block()
    # Nanofiltration
    if hasattr(m.fs, 'NF'):
        if kwargs['NF_type'] == 'ZO':
            m.fs.NF.get_costing(module=module, section='pretreatment')
        elif kwargs['NF_type'] == 'Sep':
            raise NotImplementedError("get_costing will not be implemented for the NF separator model.")
    if hasattr(m.fs, 'pump_NF'):
        m.fs.pump_NF.get_costing(module=module, section='pretreatment', pump_type="low pressure", cost_capacity=False)
    # Reverse Osmosis
    if hasattr(m.fs, 'RO'):
        if kwargs['RO_type'] == '0D':
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
        print('FOUND LIME SOFTENER')
        m.fs.stoich_softening_mixer_unit.get_costing(module=module, section='pretreatment', mixer_type="lime_softening",cost_capacity=False)
        m.fs.lime_softening_unit_capex = Expression(expr=m.fs.stoich_softening_mixer_unit.costing.capital_cost/m.fs.annual_water_production *crf)
        m.fs.lime_softening_unit_opex = Expression(expr=m.fs.stoich_softening_mixer_unit.costing.operating_cost / m.fs.annual_water_production)
    else:
        m.fs.lime_softening_unit_capex = Expression(expr=0)
        m.fs.lime_softening_unit_opex = Expression(expr=0)
    if hasattr(m.fs,'ideal_naocl_mixer_unit'): #TODO: check how posttreatment (chlorination) was implemented on flowsheet (once added in)
        # print('FOUND CHLORINATION UNIT')
        m.fs.ideal_naocl_mixer_unit.get_costing(module=module, section='post_treatment', mixer_type='naocl_mixer')
        m.fs.chlorination_unit_capex = Expression(expr=m.fs.ideal_naocl_mixer_unit.costing.capital_cost/m.fs.annual_water_production *crf)
        m.fs.chlorination_unit_opex = Expression(expr=m.fs.ideal_naocl_mixer_unit.costing.operating_cost / m.fs.annual_water_production)
    else:
        m.fs.chlorination_unit_capex = Expression(expr=0)
        m.fs.chlorination_unit_opex = Expression(expr=0)

    # TODO: Solver currently fails to converge when trying to account for permeate mixer costing (two-stage not converging)
    if hasattr(m.fs,'mixer_permeate'):
        m.fs.mixer_permeate.get_costing(module=module, section='primary')
    # call get_system_costing for whole flowsheet
    module.get_system_costing(m.fs)

    scale_costing(m)


#def get_costing_sweep(self, **kwargs):
    ## Initial attempt to do a general sweep across unit models and call get_costing
    # for b_unit in self.component_objects(Block, descend_into=True):
    #     # print(b_unit)
    #     if hasattr(b_unit, 'get_costing') and callable(b_unit.get_costing):
    #         name = getattr(b_unit, 'local_name')
    #         # if getattr(b_unit, '__class__') == 'idaes.core.process_block._ScalarPump':
    #         if isinstance(b_unit, PumpData):
    #             print(f"We got ourselves a pump called {name}!")
    #         else:
    #             print(f"We got ourselves a {name}!")
    #             # b_unit.get_costing(module=module)

def scale_costing(self):
    for b_unit in self.component_objects(Block, descend_into=True):
        # print(b_unit)
        if hasattr(b_unit, 'costing'):
            base = b_unit.costing
            for var in base.component_objects(Var):
                if iscale.get_scaling_factor(var) is None:
                    iscale.set_scaling_factor(var, 1e-3)
                for con in base.component_objects(Constraint):
                    sf = iscale.get_scaling_factor(var)
                    iscale.constraint_scaling_transform(con, sf)

def display_costing(m, **kwargs):
    crf = m.fs.costing_param.factor_capital_annualization
    # #TODO: add expressions for all cost components that we may want in LCOW breakdown bar charts
    dummy_val=0.05
    if not hasattr(m.fs, 'pump_RO2'): # TODO: remove this temporary fix meant for adding to cost_dict without error
        m.fs.pump_RO2 = Block()
        m.fs.pump_RO2.costing = Block()
        m.fs.pump_RO2.costing.operating_cost = Param(initialize=0)
    if not hasattr(m.fs, 'NF'): # TODO: remove this temporary fix meant for adding to cost_dict without error
        m.fs.NF = Block()
        m.fs.NF.costing = Block()
        m.fs.NF.costing.operating_cost = Param(initialize=0)
    if not hasattr(m.fs, 'RO2'): # TODO: remove this temporary fix meant for adding to cost_dict without error
        m.fs.RO2 = Block()
        m.fs.RO2.costing = Block()
        m.fs.RO2.costing.operating_cost = Param(initialize=0)

    # UNITS FOR ALL COST COMPONENTS [=] $/m3 of permeate water produced
    # TODO: still need to add separate chemical costs and other cost components
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
                                     / m.fs.annual_water_production,  # TODO: Presumably for RO Plant - may revise; separate chemical costs and add membrane cleaning
        # 'Total Chemical Costs': dummy_val,
        'Total Electricity Cost': m.fs.costing.electricity_cost_total / m.fs.annual_water_production,  # TODO: should other energy costs be accounted for, i.e., pretreatment? probably
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

    # cost_dict_extras = {
    #     'NF Pump Electricity Cost': m.fs.pump_NF.costing.operating_cost / m.fs.annual_water_production,
    #     'Stage 1 HP Pump Electricity Cost': m.fs.pump_RO.costing.operating_cost / m.fs.annual_water_production,
    #     'Stage 2 HP Pump Electricity Cost': m.fs.pump_RO2.costing.operating_cost / m.fs.annual_water_production,
    #     'NF Membrane Replacement Cost': m.fs.NF.costing.operating_cost / m.fs.annual_water_production,
    #     'Stage 1 RO Membrane Replacement Cost': m.fs.RO.costing.operating_cost / m.fs.annual_water_production,
    #     'Stage 2 RO Membrane Replacement Cost': m.fs.RO2.costing.operating_cost / m.fs.annual_water_production,
    #
    # }
    for item, val in cost_dict.items():
        print(f"{item} = {value(val)}")

    # print(f'LCOW = ${round(m.fs.costing.LCOW.value, 5)}/m3')
    #
    # if hasattr(m.fs,'pump_RO'):
    #     pump_RO_spec_opex= m.fs.pump_RO.costing.operating_cost.value/value(m.fs.annual_water_production)
    #     print(f'RO Pump 1 specific Opex = ${round(pump_RO_spec_opex,3)}/m3')
    # if hasattr(m.fs,'pump_RO2'):
    #     pump_RO2_spec_opex= m.fs.pump_RO2.costing.operating_cost.value/value(m.fs.annual_water_production)
    #     print(f'RO Pump 2 specific Opex = ${round(pump_RO2_spec_opex,3)}/m3')
    #
    # if hasattr(m.fs,'stoich_softening_mixer_unit'): #TODO: check if pretreatment by lime softening was implemented on flowsheet (once added in)
    #     lime_softener_spec_capex= m.fs.stoich_softening_mixer_unit.costing.capital_cost.value/value(m.fs.annual_water_production) *crf.value
    #     lime_softener_spec_opex= m.fs.stoich_softening_mixer_unit.costing.operating_cost.value/value(m.fs.annual_water_production)
    #
    #     print(f'Lime Softening specific CAPEX = ${round(lime_softener_spec_capex,5)}/m3')
    #     print(f'Lime Softening specific OPEX = ${round(lime_softener_spec_opex,5)}/m3')
    #
    # if hasattr(m.fs,'ideal_naocl_mixer_unit'):
    #     chlorination_spec_capex= m.fs.ideal_naocl_mixer_unit.costing.capital_cost.value/value(m.fs.annual_water_production) *crf.value
    #     chlorination_spec_opex= m.fs.ideal_naocl_mixer_unit.costing.operating_cost.value/value(m.fs.annual_water_production)
    #
    #     print(f'Chlorination specific CAPEX = ${round(chlorination_spec_capex,5)}/m3')
    #     print(f'Chlorination specific OPEX = ${round(chlorination_spec_opex,5)}/m3')
    # if hasattr(m.fs,'ERD'):
    #     pump_ERD_spec_opex = m.fs.ERD.costing.operating_cost.value / value(m.fs.annual_water_production)
    #     print(f'ERD specific Opex = ${round(pump_ERD_spec_opex, 3)}/m3')
    return cost_dict

if __name__ == "__main__":
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    kwargs_flowsheet= {'has_bypass': True,
                       'has_desal_feed': False,
                       'is_twostage': True,
                       'NF_type': 'ZO',
                       'NF_base': 'ion',
                       'RO_type': '0D',
                       'RO_base': 'TDS',
                       'RO_level': 'simple'
                       }
    m = solve_optimization(system_recovery=0.78, max_conc_factor=3, **kwargs_flowsheet)
