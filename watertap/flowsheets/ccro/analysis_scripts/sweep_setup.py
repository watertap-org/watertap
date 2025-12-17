import os
from parameter_sweep.loop_tool.loop_tool import loopTool, get_working_dir
import watertap.flowsheets.ccro.analysis_scripts.ccro_setup as ccro
from pathlib import Path


def main():

    loopTool(
        get_working_dir() + "/brine_sweep.yaml",
        build_function=ccro.build,
        optimize_function=ccro.solve_model,
        save_name="ccro_brine_sweep",
        saving_dir=get_working_dir(),
        number_of_subprocesses=1,
        num_loop_workers=1,
    )


if __name__ == "__main__":
    main()
