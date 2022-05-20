import sys
import os
import time

from pyomo.environ import Constraint
from watertap.tools.parameter_sweep import _init_mpi, LinearSample, parameter_sweep
import watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.dye_desalination.dye_desalination as dye_desalination


def set_up_sensitivity(m):
    outputs = {}
    optimize_kwargs = {"check_termination": False}
    opt_function = dye_desalination.solve

    # create outputs
    outputs["LCOW"] = m.fs.costing.LCOW_dye_recovered
    # outputs["Total_Cost"] = m.fs.costing.total_annualized_cost
    # outputs["LCODS"] = m.fs.costing.LCODS

    return outputs, optimize_kwargs, opt_function


def run_analysis(case_num, nx, interpolate_nan_outputs=True):

    m = dye_desalination.main()[0]

    outputs, optimize_kwargs, opt_function = set_up_sensitivity(m)

    sweep_params = {}

    if case_num == 1:
        m.fs.costing.dye_cost.unfix()
        sweep_params["dye_cost"] = LinearSample(m.fs.costing.dye_cost, 0.1, 5, nx)

    elif case_num == 2:
        m.fs.costing.waste_disposal_cost.unfix()
        sweep_params["disposal_cost"] = LinearSample(
            m.fs.costing.waste_disposal_cost, 1, 10, nx
        )
    elif case_num == 3:
        m.fs.nanofiltration.removal_frac_mass_solute[0, "dye"].unfix()
        sweep_params["dye_removal"] = LinearSample(
            m.fs.nanofiltration.removal_frac_mass_solute[0, "dye"], 0.2, 0.999, nx
        )
    elif case_num == 4:
        # sweep across membrane properties
        m.fs.nanofiltration.water_permeability_coefficient[0].unfix()
        m.fs.nanofiltration.removal_frac_mass_solute[0, "dye"].unfix()
        sweep_params["water_permeability"] = LinearSample(
            m.fs.nanofiltration.water_permeability_coefficient[0], 1, 100, nx
        )
        sweep_params["dye_removal"] = LinearSample(
            m.fs.nanofiltration.removal_frac_mass_solute[0, "dye"], 0.2, 0.999, nx
        )
    elif case_num == 5:
        # sweep across membrane cost parameters
        m.fs.costing.nanofiltration.membrane_cost.unfix()

        sweep_params["membrane_cost"] = LinearSample(
            m.fs.costing.nanofiltration.membrane_cost, 1, 100, nx
        )
    elif case_num == 6:
        m.fs.costing.dye_cost.unfix()
        m.fs.costing.waste_disposal_cost.unfix()

        sweep_params["dye_cost"] = LinearSample(m.fs.costing.dye_cost, 0, 10, nx)
        sweep_params["waste_disposal"] = LinearSample(
            m.fs.costing.waste_disposal_cost, 0, 10, nx
        )
    elif case_num == 7:
        m.fs.feed.flow_vol.unfix()
        sweep_params["feed_flow"] = LinearSample(
            m.fs.feed.flow_vol, 0.01 / 3600, 120 / 3600, nx
        )

    else:
        raise ValueError("case_num = %d not recognized." % (case_num))

    output_filename = "sensitivity_" + str(case_num)
    global_results = parameter_sweep(
        m,
        sweep_params,
        outputs,
        csv_results_file_name=output_filename,
        optimize_function=opt_function,
        optimize_kwargs=optimize_kwargs,
        interpolate_nan_outputs=interpolate_nan_outputs,
    )

    return global_results, sweep_params, m


def main(case_num=7, nx=10, interpolate_nan_outputs=True):
    # when from the command line
    case_num = int(case_num)
    nx = int(nx)
    interpolate_nan_outputs = bool(interpolate_nan_outputs)

    # comm, rank, num_procs = _init_mpi()

    global_results, sweep_params, m = run_analysis(
        case_num, nx, interpolate_nan_outputs
    )
    print(global_results)

    return global_results, m


if __name__ == "__main__":
    results, model = main()
