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
from pyomo.environ import units as pyunits, check_optimal_termination, value


def lsrro_optimize(number_of_stages=2, water_recovery=0.5, A_fixed=5 / 3.6e11, permeate_quality_limit=1000e-6, has_CP=True):
    m = lsrro_case.build(number_of_stages=number_of_stages, nacl_solubility_limit=True, has_CP =has_CP, has_Pdrop=True)
    lsrro_case.set_operating_conditions(m)
    lsrro_case.initialize(m)
    lsrro_case.solve(m)
    lsrro_case.optimize_set_up(m, A_fixed=A_fixed, permeate_quality_limit=permeate_quality_limit)
    m= lsrro_case.solve(m, raise_on_failure=False, tee=False)

    return m

def run_case(number_of_stages, nx):
    sweep_params = {}
    outputs = {}


    # lsrro_kwargs = {
    #                 # 'number_of_stages': 2,
    #                 # 'Cin': 35,
    #                 # 'water_recovery': 0.75,
    #                 'A_fixed': 5 / 3.6e11,
    #                 'permeate_quality_limit': 1000e-6,
    #                 'has_CP': True}
    # m, _ = lsrro_optimize(number_of_stages)

    m = lsrro_case.build(number_of_stages=number_of_stages, nacl_solubility_limit=True, has_CP =True, has_Pdrop=True)
    lsrro_case.set_operating_conditions(m)
    lsrro_case.initialize(m)
    lsrro_case.solve(m)
    lsrro_case.optimize_set_up(m, A_fixed=5 / 3.6e11, permeate_quality_limit=1000e-6)


    opt_function = lsrro_case.solve

    # Sweep parameters ------------------------------------------------------------------------
    # sweep_params['Number of Stages'] = LinearSample(m.fs.NumberOfStages, 2, 8) #
    #TODO: parameter_sweep is producing all nans in results, even when just sweeping reovery. Find out what the issue is.
    sweep_params['Volumetric Recovery Rate'] = LinearSample(m.fs.water_recovery, 0.5, 0.55, nx)
    #TODO: want to sweep feed concentration but don't see how state vars of feed will get updated/unfixed since they're fixed
    # sweep_params['Feed Concentration'] = LinearSample(
    #     m.fs.feed.properties[0].conc_mass_phase_comp['Liq', 'NaCl'], 35, 70, nx)

    output_filename = f'param_sweep_output/{number_of_stages}_stage/results_LSRRO.csv'

    # Outputs  -------------------------------------------------------------------------------
    outputs['LCOW'] = m.fs.costing.LCOW
    outputs['LCOW wrt Feed Flow'] = m.fs.LCOW_feed
    outputs['SEC'] = m.fs.specific_energy_consumption
    outputs['SEC wrt Feed'] = m.fs.specific_energy_consumption_feed
    outputs['Number of Stages'] = m.fs.NumberOfStages
    outputs['Final Brine Concentration'] = m.fs.disposal.properties[0].conc_mass_phase_comp['Liq','NaCl']
    outputs['Final Permeate Concentration (ppm)'] = \
        m.fs.product.properties[0].mass_frac_phase_comp['Liq', 'NaCl'] * 1e6
    outputs['Annual Feed Flow'] = m.fs.annual_feed
    outputs['Annual Water Production'] = m.fs.annual_water_production

    outputs['Pump Work In (kW)'] = m.fs.total_work_in / 1000
    outputs['Pump Work Recovered (kW)'] = m.fs.total_work_recovered / 1000
    outputs['Net Pump Work In (kW)'] = m.fs.net_pump_work / 1000
    outputs['Energy Recovery (%)'] = -m.fs.total_work_recovered / m.fs.total_work_in * 100

    outputs['Mass Water Recovery Rate (%)'] = m.fs.mass_water_recovery * 100
    outputs['System Salt Rejection (%)'] = m.fs.system_salt_rejection * 100

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
                    stage.rejection_phase_comp[0, 'Liq', 'NaCl'] * 100
                    for idx, stage in m.fs.ROUnits.items()})

    outputs.update({f'Observed Salt Passage (%), Stage {idx}':
                    (1 - stage.rejection_phase_comp[0, 'Liq', 'NaCl']) * 100
                    for idx, stage in m.fs.ROUnits.items()})

    outputs.update({f'Mass Salt Passage (%), Stage {idx}':
                    (stage.mixed_permeate[0].flow_mass_phase_comp['Liq', 'NaCl']
                    / stage.feed_side.properties[0, 0].flow_mass_phase_comp['Liq', 'NaCl']) * 100
                    for idx, stage in m.fs.ROUnits.items()})

    outputs.update({f'Volumetric Module Recovery Rate (%), Stage {idx}':
                    stage.recovery_vol_phase[0, 'Liq'] * 100
                    for idx, stage in m.fs.ROUnits.items()})

    outputs.update({f'Mass Water Module Recovery Rate (%), Stage {idx}':
                    stage.recovery_mass_phase_comp[0, 'Liq', 'H2O'] * 100
                    for idx, stage in m.fs.ROUnits.items()})

    outputs.update({f'Volumetric Stage Recovery Rate (%), Stage {idx}':
                    getattr(m.fs, f'stage{idx}_recovery_vol') * 100
                    for idx in m.fs.StageSet})

    outputs.update({f'Mass Water Stage Recovery Rate (%), Stage {idx}':
                    getattr(m.fs, f'stage{idx}_recovery_mass_H2O') * 100
                    for idx in m.fs.StageSet})

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

    outputs.update({f'Pressure Drop (bar), Stage {idx}':
                    -pyunits.convert(stage.deltaP[0], to_units=pyunits.bar)
                    for idx, stage in m.fs.ROUnits.items()})
    outputs.update({f'Inlet Reynolds Number, Stage {idx}':
                    stage.N_Re[0, 0]
                    for idx, stage in m.fs.ROUnits.items()})
    outputs.update({f'Outlet Reynolds Number, Stage {idx}':
                    stage.N_Re[0, 1]
                    for idx, stage in m.fs.ROUnits.items()})
    outputs.update({f'Inlet Crossflow Velocity, Stage {idx}':
                    stage.velocity[0, 0]
                    for idx, stage in m.fs.ROUnits.items()})
    outputs.update({f'Outlet Crossflow Velocity, Stage {idx}':
                    stage.velocity[0, 1]
                    for idx, stage in m.fs.ROUnits.items()})

    global_results = parameter_sweep(m, sweep_params, outputs, csv_results_file=output_filename,
                                     optimize_function=opt_function,
                                     optimize_kwargs={'solver': None, 'tee': True, 'raise_on_failure': True},
                                     # reinitialize_function=lsrro_case.initialize,
                                     # reinitialize_kwargs={'verbose': True, 'solver': None},
                                     # reinitialize_before_sweep=True,
                                     debugging_data_dir=os.path.split(output_filename)[0]+'/local',
                                     interpolate_nan_outputs=True)

    return global_results, sweep_params, m


if __name__ == "__main__":
    for n in range(3, 4):
        # m = lsrro_optimize(number_of_stages=n)
        global_results, sweep_params, m = run_case(number_of_stages=n, nx=4)


