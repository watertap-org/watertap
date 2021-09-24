import numpy as np
import sys
import os
from idaes.core.util import get_solver

from proteuslib.tools.parameter_sweep import _init_mpi, LinearSample, parameter_sweep

# Import all analysis flowsheets
import proteuslib.flowsheets.full_treatment_train.analysis.flowsheet_NF as fs_NF
import proteuslib.flowsheets.full_treatment_train.analysis.flowsheet_NF_no_bypass as fs_NF_no_bypass
import proteuslib.flowsheets.full_treatment_train.analysis.flowsheet_NF_two_stage as fs_NF_two_stage
import proteuslib.flowsheets.full_treatment_train.analysis.flowsheet_single_stage as fs_single_stage
import proteuslib.flowsheets.full_treatment_train.analysis.flowsheet_two_stage as fs_two_stage 
import proteuslib.flowsheets.full_treatment_train.analysis.flowsheet_softening as fs_softening
import proteuslib.flowsheets.full_treatment_train.analysis.flowsheet_softening_two_stage as fs_softening_two_stage

# Start MPI communicator
comm, rank, num_procs = _init_mpi()

# Get the case number to run
try:
    case_num = int(sys.argv[1])
except:
    # Default to running all cases
    case_num = 6

# Get the default number of discretization points
try:
    nx = int(sys.argv[2])
except:
    # Default to a 4-point discretization
    nx = 4

def add_costing_outputs_to_dict(m, units_to_cost, outputs):
    for unit_name in units_to_cost:
        for cost_type in ['capital_cost', 'operating_cost']:
            unit = getattr(m.fs, unit_name)
            cost = getattr(unit.costing, cost_type)

            outputs['%s_%s' % (cost_type, unit_name)] = cost

    return outputs


if case_num == 1:
    # ================================================================
    # Single Stage RO
    # ================================================================
    m = fs_single_stage.optimize_flowsheet(RO_type='0D')

    sweep_params = {}
    # sweep_params['System Recovery'] = LinearSample(m.fs.system_recovery_target, 0.3, 0.95, nx)
    sweep_params['RO Recovery'] = LinearSample(m.fs.RO_recovery, 0.3, 0.95, nx)

    outputs = {}
    outputs['LCOW'] = m.fs.costing.LCOW
    outputs['Saturation Index'] = m.fs.desal_saturation.saturation_index
    outputs['Pump Pressure'] = m.fs.pump_RO.control_volume.properties_out[0].pressure
    outputs = add_costing_outputs_to_dict(m, ['RO', 'pump_RO', 'ERD'], outputs)

    output_filename = 'output/fs_single_stage/results_%d.csv' % case_num

    global_results = parameter_sweep(m, sweep_params, outputs, output_filename, 
                                     optimize_function=fs_single_stage.optimize,
                                     debugging_data_dir=os.path.split(output_filename)[0]+'/local',
                                     interpolate_nan_outputs=True)

elif case_num == 2:
    # ================================================================
    # Two Stage RO
    # ================================================================
    m = fs_two_stage.optimize_flowsheet(RO_type='0D')

    sweep_params = {}
    # sweep_params['System Recovery'] = LinearSample(m.fs.system_recovery_target, 0.3, 0.95, nx)
    sweep_params['RO Recovery'] = LinearSample(m.fs.RO_recovery, 0.3, 0.95, nx)
    sweep_params['Max Pressure'] = LinearSample(m.fs.max_allowable_pressure, 300e5, 75e5, nx)

    outputs = {}
    outputs['LCOW'] = m.fs.costing.LCOW
    outputs['Saturation Index'] = m.fs.desal_saturation.saturation_index
    outputs['RO-1 Pump Pressure'] = m.fs.pump_RO.control_volume.properties_out[0].pressure
    outputs['RO-2 Pump Pressure'] = m.fs.pump_RO2.control_volume.properties_out[0].pressure
    outputs = add_costing_outputs_to_dict(m, ['RO', 'pump_RO', 'RO2', 'pump_RO2', 'ERD'], outputs)

    output_filename = 'output/fs_two_stage/results_%d.csv' % case_num

    global_results = parameter_sweep(m, sweep_params, outputs, output_filename, 
                                     optimize_function=fs_two_stage.optimize,
                                     debugging_data_dir=os.path.split(output_filename)[0]+'/local',
                                     interpolate_nan_outputs=True)

elif case_num == 3:
    # ================================================================
    # NF No Bypass
    # ================================================================
    m = fs_NF_no_bypass.solve_flowsheet()

    sweep_params = {}
    sweep_params['NF Recovery'] = LinearSample(m.fs.NF.recovery_mass_phase_comp[0, 'Liq', 'H2O'], .3, .95, nx)

    outputs = {}
    outputs['LCOW'] = m.fs.costing.LCOW
    outputs['Saturation Index'] = m.fs.pretrt_saturation.saturation_index

    output_filename = 'output/fs_NF_no_bypass/results_%d.csv' % case_num

    global_results = parameter_sweep(m, sweep_params, outputs, output_filename, 
                                     optimize_function=fs_NF_no_bypass.simulate,
                                     debugging_data_dir=os.path.split(output_filename)[0]+'/local',
                                     interpolate_nan_outputs=True)

elif case_num == 4:
    # ================================================================
    # NF Two Stage
    # ================================================================
    m = fs_NF_two_stage.optimize_flowsheet(RO_type='0D')

    sweep_params = {}
    sweep_params['NF Split Fraction'] = LinearSample(m.fs.splitter.split_fraction[0, 'bypass'], .01, .99, nx)

    outputs = {}
    outputs['Ca Removal'] = m.fs.removal_Ca

    output_filename = 'output/fs_NF_two_stage/results_%d.csv' % case_num

    global_results = parameter_sweep(m, sweep_params, outputs, output_filename, 
                                     optimize_function=fs_NF_two_stage.optimize,
                                     debugging_data_dir=os.path.split(output_filename)[0]+'/local',
                                     interpolate_nan_outputs=True)

elif case_num == 5:
    # ================================================================
    # NF Two Stage
    # ================================================================
    m = fs_NF_two_stage.optimize_flowsheet(RO_type='0D')

    sweep_params = {}
    sweep_params['NF Split Fraction'] = LinearSample(m.fs.splitter.split_fraction[0, 'bypass'], 0.25, 1.0, 4)
    sweep_params['RO Recovery'] = LinearSample(m.fs.RO_recovery, 0.5, 0.95, nx)

    outputs = {}
    outputs['LCOW'] = m.fs.costing.LCOW

    output_filename = 'output/fs_NF_two_stage/results_%d.csv' % case_num

    global_results = parameter_sweep(m, sweep_params, outputs, output_filename, 
                                     optimize_function=fs_NF_two_stage.optimize,
                                     debugging_data_dir=os.path.split(output_filename)[0]+'/local',
                                     interpolate_nan_outputs=True)

elif case_num == 6:
    # ================================================================
    # NF Two Stage
    # ================================================================
    m = fs_NF_two_stage.optimize_flowsheet(RO_type='0D')

    sweep_params = {}
    sweep_params['System Recovery'] = LinearSample(m.fs.system_recovery_target, 0.5, 0.95, nx)

    outputs = {}
    outputs['NF Split Fraction'] = m.fs.splitter.split_fraction[0, 'bypass']

    output_filename = 'output/fs_NF_two_stage/results_%d.csv' % case_num

    global_results = parameter_sweep(m, sweep_params, outputs, output_filename, 
                                     optimize_function=fs_NF_two_stage.optimize,
                                     debugging_data_dir=os.path.split(output_filename)[0]+'/local',
                                     interpolate_nan_outputs=True)

elif case_num == 7:
    # ================================================================
    # NF Two Stage
    # ================================================================
    m = fs_NF_two_stage.optimize_flowsheet(RO_type='0D')

    sweep_params = {}
    sweep_params['NF Split Fraction'] = LinearSample(m.fs.splitter.split_fraction[0, 'bypass'], 0.25, 1.0, 4)
    sweep_params['System Recovery'] = LinearSample(m.fs.system_recovery_target, 0.5, 0.95, nx)

    outputs = {}
    outputs['LCOW'] = m.fs.costing.LCOW

    output_filename = 'output/fs_NF_two_stage/results_%d.csv' % case_num

    global_results = parameter_sweep(m, sweep_params, outputs, output_filename, 
                                     optimize_function=fs_NF_two_stage.optimize,
                                     debugging_data_dir=os.path.split(output_filename)[0]+'/local',
                                     interpolate_nan_outputs=True)

elif case_num == 8:
    # ================================================================
    # Softening
    # ================================================================
    m = fs_softening.solve_flowsheet()

    sweep_params = {}
    sweep_params['Lime Dosing'] = LinearSample(m.fs.stoich_softening_mixer_unit.lime_stream_state[0].flow_mol, 0.0, 0.128, nx)

    outputs = {}
    outputs['LCOW'] = m.fs.costing.LCOW
    outputs['Ca Removal'] = m.fs.removal_Ca
    outputs['Mg Removal'] = m.fs.removal_Mg

    output_filename = 'output/fs_softening/results_%d.csv' % case_num

    global_results = parameter_sweep(m, sweep_params, outputs, output_filename, 
                                     optimize_function=fs_softening.simulate,
                                     debugging_data_dir=os.path.split(output_filename)[0]+'/local',
                                     interpolate_nan_outputs=True)

elif case_num == 9:
    # ================================================================
    # Softening Two Stage
    # ================================================================
    m = fs_softening_two_stage.optimize_flowsheet()

    sweep_params = {}
    sweep_params['System Recovery'] = LinearSample(m.fs.system_recovery_target, 0.5, 0.95, nx)

    outputs = {}
    outputs['LCOW'] = m.fs.costing.LCOW
    outputs['Lime Dosing'] = m.fs.stoich_softening_mixer_unit.lime_stream_state[0].flow_mol

    output_filename = 'output/fs_softening_two_stage/results_%d.csv' % case_num

    global_results = parameter_sweep(m, sweep_params, outputs, output_filename, 
                                     optimize_function=fs_softening_two_stage.optimize,
                                     debugging_data_dir=os.path.split(output_filename)[0]+'/local',
                                     interpolate_nan_outputs=True)

else:
    raise ValueError('case_num = %d not recognized.' % (case_num))
