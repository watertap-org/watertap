import os
from parameter_sweep.loop_tool.loop_tool import loopTool, get_working_dir
import watertap.flowsheets.ccro.analysis_scripts.flush_setup as psro
from pathlib import Path


def main():

    loopTool(
        get_working_dir() + "/flush_sweep.yaml",
        build_function=psro.build,
        optimize_function=psro.solve_model,
        save_name="flush_sweep",
        saving_dir=get_working_dir(),
        number_of_subprocesses=1,
        num_loop_workers=1,
    )


if __name__ == "__main__":
    main()
