import os
from parameter_sweep.loop_tool.loop_tool import loopTool, get_working_dir
import watertap.flowsheets.ccro.analysis_scripts.ccro_setup as ccro
from pathlib import Path

here = os.path.dirname(os.path.abspath(__file__))

def main():

    # loopTool(
    #     get_working_dir() + "/brine_sweep.yaml",
    #     build_function=ccro.build,
    #     optimize_function=ccro.solve_model,
    #     save_name="ccro_brine_sweep",
    #     saving_dir=get_working_dir(),
    #     number_of_subprocesses=1,
    #     num_loop_workers=1,
    # )

    # loopTool(
    #     get_working_dir() + "/map_sweep.yaml",
    #     build_function=ccro.build,
    #     optimize_function=ccro.solve_model,
    #     save_name="ccro_map_sweep",
    #     saving_dir=get_working_dir(),
    #     number_of_subprocesses=1,
    #     num_loop_workers=1,
    # )
    # loopTool(
    #     get_working_dir() + "/flow_sweep.yaml",
    #     build_function=ccro.build,
    #     optimize_function=ccro.solve_model,
    #     save_name="ccro_flow_sweep",
    #     saving_dir=get_working_dir(),
    #     number_of_subprocesses=1,
    #     num_loop_workers=1,
    # )
    # loopTool(
    #     get_working_dir() + "/mesh_sweep.yaml",
    #     build_function=ccro.build,
    #     optimize_function=ccro.solve_model,
    #     save_name="ccro_bw_mesh",
    #     saving_dir=get_working_dir(),
    #     number_of_subprocesses=1,
    #     num_loop_workers=2,
    # )

    loopTool(
        here + "/recovery_sweep.yaml",
        build_function=ccro.build,
        optimize_function=ccro.solve_model,
        save_name="ccro_recovery_sweep",
        saving_dir=here,
        number_of_subprocesses=1,
        num_loop_workers=1,
    )

    # loopTool(
    #     here + "/flush_eff_sweep.yaml",
    #     build_function=ccro.build_for_flush_eff,
    #     optimize_function=ccro.solve_model,
    #     save_name="ccro_flush_eff_recovery",
    #     saving_dir=here,
    #     number_of_subprocesses=1,
    #     num_loop_workers=1,
    # )


if __name__ == "__main__":
    main()
