#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
from watertap.tools.parameter_sweep import (
    LinearSample,
    parameter_sweep,
)
import watertap.examples.flowsheets.case_studies.electroNP.electroNP_flowsheet as electroNP_flowsheet
from pyomo.environ import units as pyunits


def set_up_sensitivity(m):
    outputs = {}
    optimize_kwargs = {"fail_flag": False}
    opt_function = electroNP_flowsheet.solve

    # create outputs
    # outputs["ElectroNP Capital Cost"] = m.fs.electroNP.costing.capital_cost
    outputs["Electricity"] = pyunits.convert(
        m.fs.costing.aggregate_flow_costs["electricity"] / m.fs.AD.inlet.flow_vol[0],
        to_units=pyunits.USD_2018 / pyunits.m**3,
    )

    return outputs, optimize_kwargs, opt_function


def run_analysis(case_num=3, nx=11, interpolate_nan_outputs=True, results_path=None):

    if results_path is None:
        results_path = f"sensitivity_analysis{case_num}.csv"

    m = electroNP_flowsheet.build_flowsheet()[0]

    outputs, optimize_kwargs, opt_function = set_up_sensitivity(m)

    sweep_params = {}
    if case_num == 1:
        # sensitivity analysis
        m.fs.costing.electroNP.sizing_cost.unfix()
        sweep_params["sizing_cost"] = LinearSample(
            m.fs.costing.electroNP.sizing_cost, 1, 30, nx
        )
    elif case_num == 2:
        m.fs.costing.electroNP.sizing_cost.unfix()
        sweep_params["sizing_cost"] = LinearSample(
            m.fs.costing.electroNP.sizing_cost, 1, 30, nx
        )
        m.fs.costing.electroNP.HRT.unfix()
        sweep_params["HRT"] = LinearSample(m.fs.costing.electroNP.HRT, 1, 5, nx)
    elif case_num == 3:
        sweep_params["energy_consumption"] = LinearSample(
            m.fs.electroNP.energy_electric_flow_mass, 0, 3, nx
        )
    elif case_num == 4:
        sweep_params["energy_consumption"] = LinearSample(
            m.fs.electroNP.energy_electric_flow_mass, 0, 3, nx
        )
        m.fs.costing.electricity_cost.unfix()
        sweep_params["electricity_cost"] = LinearSample(
            m.fs.costing.electricity_cost, 0.05, 0.1, nx
        )
    elif case_num == 5:
        sweep_params["electricity_cost"] = LinearSample(
            m.fs.costing.electricity_cost, 0.05, 0.1, nx
        )
        sweep_params["MgCl2_cost"] = LinearSample(
            m.fs.costing.magnesium_chloride_cost, 0.07, 0.1, nx
        )
    else:
        raise ValueError(f"{case_num} is not yet implemented")

    global_results = parameter_sweep(
        m,
        sweep_params,
        outputs,
        csv_results_file_name=results_path,
        optimize_function=opt_function,
        optimize_kwargs=optimize_kwargs,
        interpolate_nan_outputs=interpolate_nan_outputs,
    )

    return global_results, sweep_params, m


if __name__ == "__main__":
    results, sweep_params, m = run_analysis()
