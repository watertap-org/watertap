import sys
import os
import time

from idaes.core.util import get_solver
from pyomo.environ import Expression, Param, Constraint
from watertap.tools.parameter_sweep import _init_mpi, LinearSample, parameter_sweep
import watertap.examples.flowsheets.high_pressure_RO.high_pressure_RO as high_pressure_RO


def run_analysis(case, nx, interp_nan_outputs=False):
    sweep_params = {}
    outputs = {}
    optimize_kwargs = {'check_termination': False}

    if case == 'seawater':
        # build, set, and initialize
        m = high_pressure_RO.build(case=case)
        high_pressure_RO.specify_model(m)
        high_pressure_RO.initialize_model(m)
        # simulate
        high_pressure_RO.solve(m)
        # set up optimize
        high_pressure_RO.set_up_optimization(m)
        high_pressure_RO.optimize(m)

        # set up parameter sweep
        sweep_params['System Recovery'] = LinearSample(m.fs.product_recovery, 0.3, 0.85, nx)

        for j in m.fs.disposal.properties[0].conc_mass_comp:
            outputs[j] = m.fs.disposal.properties[0].conc_mass_comp[j]

        outputs['LCOW'] = m.fs.costing.LCOW

        output_filename = 'output/oli_cases/seawater.csv'

        opt_function = high_pressure_RO.optimize
    else:
        raise ValueError('case %s not recognized.' % case)


    global_results = parameter_sweep(m, sweep_params, outputs, csv_results_file=output_filename,
                                     optimize_function=opt_function,
                                     optimize_kwargs=optimize_kwargs,
                                     debugging_data_dir=os.path.split(output_filename)[0]+'/local',
                                     interpolate_nan_outputs=interp_nan_outputs)

    return global_results, sweep_params


if __name__ == "__main__":

    # Start MPI communicator
    comm, rank, num_procs = _init_mpi()

    # analysis
    case = 'seawater'
    nx = 30

    tic = time.time()
    global_results, sweep_params = run_analysis(case, nx)
    print(global_results)
    toc = time.time()

    if rank == 0:
        total_samples = 1

        for k, v in sweep_params.items():
            total_samples *= v.num_samples

        print('Finished case = %s.' % case)
        print('Processed %d swept parameters comprising %d total points.' % (len(sweep_params), total_samples))
        print('Elapsed time = %.1f s.' % (toc-tic))
