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
from watertap.flowsheets.ccro.utils.utils import (
    calculate_operating_pressure,
    report_pump,
    report_ro,
    report_costing,
    report_n_stage_system,
    relax_bounds_for_low_salinity_waters,
)

here = os.path.dirname(os.path.abspath(__file__))


def optimize_for_recovery(m, **kwargs):

    if m.salt_mass_frac < 15e-3:
        for n, stage in m.fs.stage.items():
            relax_bounds_for_low_salinity_waters(stage.RO)

    for n, stage in m.fs.stage.items():
        stage.RO.area.unfix()
        stage.RO.width.unfix()
        if stage.add_pump:
            stage.pump.control_volume.properties_out[0].pressure.unfix()

    results = ro.solve_model(m, **kwargs)

    return results


def main():

    timestamp = datetime.now().strftime("%Y-%m-%d_%H:%M")
    save_file = "ro_w_erd"

    loopTool(
        f"{here}/salinity_recovery_sweep.yaml",
        build_function=ro.run_n_stage_system,
        optimize_function=optimize_for_recovery,
        save_name=save_file,
        saving_dir=here,
        number_of_subprocesses=1,
        h5_backup=False, 
        num_loop_workers=8,
    )

    # m = ro.run_n_stage_system(
    #     flow_vol=1e-3,
    #     salt_mass_frac=35e-3,
    #     water_recovery=0.5,
    #     pump_dict={1: True, 2: False},
    # )

    # results = optimize_for_recovery(m)

    # report_n_stage_system(m)


if __name__ == "__main__":
    main()
