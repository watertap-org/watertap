import os
from datetime import datetime

from pyomo.environ import check_optimal_termination, units as pyunits

from idaes.core.util.model_statistics import degrees_of_freedom

from parameter_sweep.loop_tool.loop_tool import loopTool, get_working_dir

import watertap.flowsheets.ccro.ro_w_erd.RO_w_ERD_stage as ro
from watertap.flowsheets.ccro.utils.utils import report_n_stage_system
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.core.util.initialization import *
from watertap.core.solvers import get_solver


here = os.path.dirname(os.path.abspath(__file__))


def solve_ro_w_erd(m, solver=None, max_iter=3000, tee=True, raise_on_failure=True):

    if solver is None:
        solver = get_solver()
    # m.fs.system_recovery.unfix()
    # m.fs.system_recovery.fix(0.5)
    for n, stage in m.fs.stage.items():
        stage.RO.area.unfix()
        if stage.add_pump:
            stage.pump.control_volume.properties_out[0].pressure.unfix()

    print("\n--------- SOLVING ---------\n")
    print(f"Degrees of Freedom: {degrees_of_freedom(m)}")

    results = solver.solve(m, tee=tee, options={"max_iter": max_iter})

    if check_optimal_termination(results):
        print("\n--------- OPTIMAL SOLVE!!! ---------\n")
        return results
    print_close_to_bounds(m)
    if raise_on_failure:
        print("\n--------- INFEASIBLE SOLVE!!! ---------\n")

        print("\n--------- CLOSE TO BOUNDS ---------\n")
        print_close_to_bounds(m)

        print("\n--------- INFEASIBLE BOUNDS ---------\n")
        print_infeasible_bounds(m)

        print("\n--------- INFEASIBLE CONSTRAINTS ---------\n")
        print_infeasible_constraints(m)

    return results


def main():

    timestamp = datetime.now().strftime("%Y-%m-%d_%H:%M")
    save_file = "ro_w_erd"

    loopTool(
        f"{here}/salinity_recovery_sweep.yaml",
        build_function=ro.run_n_stage_system,
        optimize_function=solve_ro_w_erd,
        save_name=save_file,
        saving_dir=here,
        number_of_subprocesses=1,
        num_loop_workers=8,
    )

    # m = ro.run_n_stage_system(
    #     flow_vol=1e-3,
    #     salt_mass_frac=35e-3,
    #     water_recovery=0.5,
    #     pump_dict={1: True, 2: False},
    # )

    # results = solve_ro_w_erd(m)

    # report_n_stage_system(m)


if __name__ == "__main__":
    main()
