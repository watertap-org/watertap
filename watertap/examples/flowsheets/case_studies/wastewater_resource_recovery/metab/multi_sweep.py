import sys
import os
import time

from pyomo.environ import Constraint
from watertap.tools.parameter_sweep import _init_mpi, LinearSample, parameter_sweep
import watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.metab.metab as metab


def set_up_sensitivity(m):
    outputs = {}
    optimize_kwargs = {"check_termination": False}  # None
    opt_function = metab.solve

    # tie parameters together
    m.fs.costing.metab.eq_bead_cost = Constraint(
        expr=m.fs.costing.metab.bead_cost['hydrogen']
             == m.fs.costing.metab.bead_cost['methane']
    )
    m.fs.costing.metab.bead_cost['methane'].unfix()

    # new basline parameters
    m.fs.costing.metab.bead_cost['hydrogen'].fix(14.4)
    # m.fs.costing.metab.hydraulic_retention_time['hydrogen'].fix()
    # m.fs.costing.metab.hydraulic_retention_time['methane'].fix()

    # create outputs
    outputs["LCOW"] = m.fs.costing.LCOW
    outputs["LCOH"] = m.fs.costing.LCOH
    outputs["LCOM"] = m.fs.costing.LCOM

    return outputs, optimize_kwargs, opt_function


def run_analysis(case_num, nx, interp_nan_outputs=True):
    m, _ = metab.main()

    outputs, optimize_kwargs, opt_function = set_up_sensitivity(m)

    sweep_params = {}
    if case_num == 1:
        # bead cost
        m.fs.costing.metab.bead_cost.display()
        sweep_params["bead_cost"] = LinearSample(
            m.fs.costing.metab.bead_cost['hydrogen'], 5, 50, nx
        )
        output_filename = 'sensitivity_'+str(case_num)+'.csv'

    elif case_num == 2:
        pass
        # # ================================================================
        # # Two Stage RO
        # # ================================================================
        # import watertap.examples.flowsheets.full_treatment_train.analysis.flowsheet_two_stage as fs_two_stage
        #
        # m = fs_two_stage.optimize_flowsheet(**desal_kwargs)
        #
        # sweep_params["System Recovery"] = LinearSample(
        #     m.fs.system_recovery_target, 0.3, 0.95, nx
        # )
        # # sweep_params['RO Recovery'] = LinearSample(m.fs.RO_recovery, 0.3, 0.95, nx)
        # sweep_params["Max Pressure"] = LinearSample(
        #     m.fs.max_allowable_pressure, 300e5, 75e5, nx
        # )
        #
        # outputs["LCOW"] = m.fs.costing.LCOW
        # outputs["Saturation Index"] = m.fs.desal_saturation.saturation_index
        # outputs["RO-1 Pump Pressure"] = m.fs.pump_RO.control_volume.properties_out[
        #     0
        # ].pressure
        # outputs["RO-2 Pump Pressure"] = m.fs.pump_RO2.control_volume.properties_out[
        #     0
        # ].pressure
        # outputs["RO Recovery"] = m.fs.RO_recovery
        # outputs["Annual Water Production"] = m.fs.costing.annual_water_production
        # outputs = append_costing_outputs(
        #     m, outputs, ["RO", "pump_RO", "RO2", "pump_RO2", "ERD"]
        # )
        #
        # output_filename = (
        #     base_path
        #     + f"output/fs_two_stage/results_{case_num}_{desal_kwargs['RO_type']}RO.csv"
        # )
        #
        # opt_function = fs_two_stage.optimize

    else:
        raise ValueError("case_num = %d not recognized." % (case_num))

    global_results = parameter_sweep(
        m,
        sweep_params,
        outputs,
        results_file_name=output_filename,
        write_csv=True,
        optimize_function=opt_function,
        optimize_kwargs=optimize_kwargs,
        # debugging_data_dir=os.path.split(output_filename)[0] + '/local',
        interpolate_nan_outputs=interp_nan_outputs,
    )

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
        nx = 11

    tic = time.time()
    global_results, sweep_params = run_analysis(case_num, nx)
    print(global_results)
    toc = time.time()

    if rank == 0:
        total_samples = 1

        for k, v in sweep_params.items():
            total_samples *= v.num_samples

        print("Finished case_num = %d." % (case_num))
        print(
            "Processed %d swept parameters comprising %d total points."
            % (len(sweep_params), total_samples)
        )
        print("Elapsed time = %.1f s." % (toc - tic))
