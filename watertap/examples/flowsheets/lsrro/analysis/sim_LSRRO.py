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

from watertap.examples.flowsheets.lsrro.lsrro import build_lsrro_case
from pyomo.environ import value

def main():
    num_stages = []
    total_area = []
    final_lcow = []
    final_sec = []
    final_perm = []
    lcow_breakdown = {}
    for stage in range(4, 5):
        m = build_lsrro_case(number_of_stages=stage,
                             # water_recovery=0.7,
                             Cin=70,
                             Cbrine=100,  # mg/L
                             A_case="optimize",
                             B_case="optimize",
                             AB_tradeoff="inequality constraint",
                             nacl_solubility_limit=True,
                             permeate_quality_limit=1000e-6,
                             has_CP=True,
                             has_Pdrop=True,
                             A_fixed=1.5 / 3.6e11  # 2.78e-12
                             )

        num_stages.append(value(m.fs.NumberOfStages))
        total_area.append(value(sum(m.fs.ROUnits[a].area for a in range(1, m.fs.NumberOfStages + 1))))
        final_lcow.append(value(m.fs.costing.LCOW))
        final_sec.append(value(m.fs.specific_energy_consumption))
        final_perm.append(value(m.fs.product.flow_mass_phase_comp[0, 'Liq', 'NaCl'] /
                                sum(m.fs.product.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl'])))
        lcow_breakdown[f'primary_pump_capex {stage}'] = (
            value(m.fs.costing_param.factor_capital_annualization * sum(
                m.fs.PrimaryPumps[stage].costing.capital_cost for stage in m.fs.StageSet) / m.fs.annual_water_production))
        lcow_breakdown[f'booster_pump_capex {stage}'] = (
            value(m.fs.costing_param.factor_capital_annualization * sum(
                m.fs.BoosterPumps[stage].costing.capital_cost for stage in
                m.fs.LSRRO_StageSet) / m.fs.annual_water_production))
        lcow_breakdown[f'erd_capex {stage}'] = (
            value(
                m.fs.costing_param.factor_capital_annualization * m.fs.EnergyRecoveryDevice.costing.capital_cost / m.fs.annual_water_production))
        lcow_breakdown[f'membrane_capex {stage}'] = (
            value(m.fs.costing_param.factor_capital_annualization * sum(
                m.fs.ROUnits[stage].costing.capital_cost for stage in m.fs.StageSet) / m.fs.annual_water_production))
        lcow_breakdown[f'indirect_capex {stage}'] = (
            value(m.fs.costing_param.factor_capital_annualization * (
                        m.fs.costing.investment_cost_total - m.fs.costing.capital_cost_total) / m.fs.annual_water_production))
        lcow_breakdown[f'electricity {stage}'] = (
            value((sum(m.fs.PrimaryPumps[stage].costing.operating_cost for stage in m.fs.StageSet)
                   + sum(m.fs.BoosterPumps[stage].costing.operating_cost for stage in m.fs.LSRRO_StageSet)
                   + m.fs.EnergyRecoveryDevice.costing.operating_cost) / m.fs.annual_water_production))

    return m, num_stages, total_area, final_lcow, final_sec, final_perm, lcow_breakdown

if __name__ == '__main__':
    m, num_stages, total_area, final_lcow, final_sec, final_perm, lcow_breakdown = main()