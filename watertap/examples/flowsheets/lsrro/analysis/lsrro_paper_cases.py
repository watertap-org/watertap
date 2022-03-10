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

import sys
import os
import time

from watertap.tools.parameter_sweep import _init_mpi, LinearSample, parameter_sweep
from watertap.examples.flowsheets.lsrro.analysis import lsrro_paper_analysis as lsrro_case
from pyomo.environ import units as pyunits, check_optimal_termination


def run_case(number_of_stages, Cin, water_recovery, A_fixed, permeate_quality_limit, has_CP, nx):
    sweep_params = {}
    outputs = {}

    m = lsrro_case.build(number_of_stages, has_CP=has_CP)
    lsrro_case.set_operating_conditions(m, Cin=Cin)
    lsrro_case.initialize(m)
    lsrro_case.solve(m)
    lsrro_case.optimize_set_up(m, water_recovery=water_recovery, A_fixed=A_fixed,
                               permeate_quality_limit=permeate_quality_limit)
    m, res = lsrro_case.solve(m, raise_on_failure=False, tee=False)

    # Sweep parameters ------------------------------------------------------------------------
    sweep_params['Volumetric Recovery Rate'] = LinearSample(m.fs.water_recovery, 0.3, 0.9, nx)
    sweep_params['Feed Concentration'] = LinearSample(
        m.fs.feed.properties[0].conc_mass_phase_comp['Liq', 'NaCl'], 5, 160, nx)

    # Outputs  -------------------------------------------------------------------------------
    outputs['LCOW'] = m.fs.costing.LCOW
    outputs['SEC'] = m.fs.specific_energy_consumption
    outputs['SEC of Feed'] = m.fs.specific_energy_consumption_feed
    outputs['Number of Stages'] = m.fs.NumberOfStages
    outputs['Final Brine Concentration'] = m.fs.disposal.properties[0].conc_mass_phase_comp['Liq','NaCl']
    outputs['Final Permeate Concentration (ppm)'] = \
        m.fs.product.properties[0].mass_frac_phase_comp['Liq', 'NaCl'] * 1e6

    outputs['Total Membrane Area'] = m.fs.total_membrane_area
    outputs['Total Capex LCOW'] = (m.fs.costing.investment_cost_total * m.fs.costing_param.factor_capital_annualization
                                   / m.fs.annual_water_production)
    outputs['Total Opex LCOW'] = m.fs.costing.operating_cost_total / m.fs.annual_water_production

    outputs['Primary Pump Capex LCOW'] = m.fs.costing.primary_pump_capex_lcow
    outputs['Booster Pump Capex LCOW'] = m.fs.costing.booster_pump_capex_lcow
    outputs['ERD Capex LCOW'] = m.fs.costing.erd_capex_lcow
    outputs['Membrane Capex LCOW'] = m.fs.costing.membrane_capex_lcow
    outputs['Indirect Capex LCOW'] = m.fs.costing.indirect_capex_lcow
    outputs['Electricity LCOW'] = m.fs.costing.electricity_lcow
    outputs['Membrane Replacement LCOW'] = m.fs.costing.membrane_replacement_lcow
    outputs['Chem,labor,maintenance LCOW'] = m.fs.costing.chemical_labor_maintenance_lcow

    outputs['Pumping Energy Agg LCOW'] = m.fs.costing.pumping_energy_aggregate_lcow
    outputs['Membrane Agg LCOW'] = m.fs.costing.membrane_aggregate_lcow

    outputs.update({f'Feed Pressure (bar), Stage {idx}':
                    pyunits.convert(pump.control_volume.properties_out[0].pressure, to_units=pyunits.bar)
                    for idx, pump in m.fs.PrimaryPumps.items()})

    outputs.update({f'Membrane Area, Stage {idx}':
                    stage.area
                    for idx, stage in m.fs.ROUnits.items()})

    outputs.update({f'Observed Rejection (%), Stage {idx}':
                    stage.rejection_phase_comp[0, 'Liq', 'NaCl']
                    for idx, stage in m.fs.ROUnits.items()})

    outputs.update({f'A-Value (LMH/bar), Stage {idx}':
                    stage.A_comp[0, 'H2O'] * 3.6e11
                    for idx, stage in m.fs.ROUnits.items()})

    outputs.update({f'B-Value (LMH), Stage {idx}':
                    stage.B_comp[0, 'NaCl'] * 3.6e6
                    for idx, stage in m.fs.ROUnits.items()})

    outputs.update({f'Average Water Flux (LMH), Stage {idx}':
                    stage.flux_mass_phase_comp_avg[0, 'Liq', 'H2O'] * 3.6e3
                    for idx, stage in m.fs.ROUnits.items()})

    outputs.update({f'Average NaCl Flux (GMH), Stage {idx}':
                    stage.flux_mass_phase_comp_avg[0, 'Liq', 'NaCl'] * 3.6e6
                    for idx, stage in m.fs.ROUnits.items()})

    #TODO:
    #  - consider pressure drop
    #  - stage-wise volumetric and mass water recovery
    #  - system-level mass water recovery
    #  - stage-wise MASS FLOW salt passage/salt recovery (mass flow salt perm/mass flow salt feed to stage)
    #  - system-wide rejection ( 1- cp,product/cf,feed)
    #  - levelized cost of water wrt to feed
    #  - Fraction of energy recovered (-ERD total energy recovered/Pump total energy)
    #  - stage-wise Inlet Re number
    #  - Stage-wise pressure drop
    #

    return m, res, sweep_params, outputs


if __name__ == "__main__":
    m, res, sweep_params, outputs = run_case(number_of_stages=2,
                                             Cin=35,
                                             water_recovery=0.75,
                                             A_fixed=5 / 3.6e11,
                                             permeate_quality_limit=1000e-6,
                                             has_CP=True,
                                             nx=4
                                             )

    from pyomo.environ import value

    for i, v in outputs.items():
        print(i, value(v))
    if not check_optimal_termination(res):
        print("solve failed")