import sys
import os
import time

from pyomo.environ import Expression
from watertap.tools.parameter_sweep import _init_mpi, LinearSample, parameter_sweep

def run_analysis(case_num, nx, RO_type, interp_nan_outputs=True):

    desal_kwargs = {'has_desal_feed': False, 'is_twostage': True, 'has_ERD': True,
                    'RO_type': RO_type, 'RO_base': 'TDS', 'RO_level': 'detailed'}

    sweep_params = {}
    outputs = {}
    optimize_kwargs = {'check_termination': False}  # None

    if case_num == 1:
        # ================================================================
        # Two Stage RO
        # ================================================================
        import watertap.examples.flowsheets.full_treatment_train.analysis.flowsheet_two_stage as fs_two_stage

        m = fs_two_stage.optimize_flowsheet(**desal_kwargs)
        # deactivate eNRTL constraint
        m.fs.eq_max_saturation_index_desal.deactivate()

        mw_comp = {'H2O': 18.015e-3,
                   'Na_+': 22.990e-3,
                   'Ca_2+': 40.078e-3,
                   'Mg_2+': 24.305e-3,
                   'SO4_2-': 96.06e-3,
                   'Cl_-': 35.453e-3}

        def mass_frac_calc(blk, j):
            return (blk.desal_saturation.properties[0].flow_mol_phase_comp['Liq', j] * mw_comp[j]
                    / sum(blk.desal_saturation.properties[0].flow_mol_phase_comp['Liq', jj] * mw_comp[jj]
                          for jj in mw_comp.keys()))
        m.fs.mass_frac_comp = Expression(mw_comp.keys(), rule=mass_frac_calc)

        sweep_params['System Recovery'] = LinearSample(m.fs.system_recovery_target, 0.3, 0.85, nx)
        # sweep_params['RO Recovery'] = LinearSample(m.fs.RO_recovery, 0.3, 0.95, nx)
        # sweep_params['Max Pressure'] = LinearSample(m.fs.max_allowable_pressure, 300e5, 75e5, nx)

        for j in mw_comp.keys():
            outputs[j] = m.fs.mass_frac_comp[j]
        outputs['LCOW'] = m.fs.costing.LCOW

        output_filename = 'output/oli_cases/results_two_stage_RO.csv'

        opt_function = fs_two_stage.optimize

    # elif case_num == 2:
        # ================================================================
        # NF
        # ================================================================

    else:
        raise ValueError('case_num = %d not recognized.' % (case_num))


    global_results = parameter_sweep(m, sweep_params, outputs, csv_results_file=output_filename,
                                     optimize_function=opt_function,
                                     optimize_kwargs=optimize_kwargs,
                                     debugging_data_dir=os.path.split(output_filename)[0]+'/local',
                                     interpolate_nan_outputs=interp_nan_outputs)

    return global_results, sweep_params


if __name__ == "__main__":

    # Start MPI communicator
    comm, rank, num_procs = _init_mpi()

    # Get the case number to run
    try:
        case_num = int(sys.argv[1])
    except:
        # Default to running case 1
        case_num = 1

    # Get the default number of discretization points
    try:
        nx = int(sys.argv[2])
    except:
        # Default to a 4-point discretization
        nx = 4

    tic = time.time()
    global_results, sweep_params = run_analysis(case_num, nx, RO_type='1D')
    print(global_results)
    toc = time.time()

    if rank == 0:
        total_samples = 1

        for k, v in sweep_params.items():
            total_samples *= v.num_samples

        print('Finished case_num = %d.' % (case_num))
        print('Processed %d swept parameters comprising %d total points.' % (len(sweep_params), total_samples))
        print('Elapsed time = %.1f s.' % (toc-tic))