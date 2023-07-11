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
    outputs["LCOW"] = m.fs.costing.LCOW
    outputs["ElectroNP Capital Cost"] = m.fs.electroNP.costing.capital_cost
    outputs["Electricity"] = pyunits.convert(
        m.fs.costing.aggregate_flow_costs["electricity"] / m.fs.AD.inlet.flow_vol[0],
        to_units=pyunits.USD_2018 / pyunits.m**3,
    )
    outputs[
        "S_PO4 Concentration, treated effluent stream"
    ] = m.fs.electroNP.treated.conc_mass_comp[0, "S_PO4"]

    return outputs, optimize_kwargs, opt_function


def run_analysis(case_num=1, nx=11, interpolate_nan_outputs=True, output_filename=None):

    if output_filename is None:
        output_filename = "sensitivity_" + str(case_num) + ".csv"

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
            m.fs.costing.electroNP.sizing_cost, 0, 30, nx
        )
        m.fs.costing.electroNP.HRT.unfix()
        sweep_params["HRT"] = LinearSample(m.fs.costing.electroNP.HRT, 0.1, 5, nx)
    elif case_num == 3:
        sweep_params["energy_consumption"] = LinearSample(
            m.fs.electroNP.energy_electric_flow_mass, 0.04, 2.5, nx
        )
    elif case_num == 4:
        sweep_params["energy_consumption"] = LinearSample(
            m.fs.electroNP.energy_electric_flow_mass, 0, 2.5, nx
        )
        sweep_params["electricity_cost"] = LinearSample(
            m.fs.costing.electricity_cost, 0.02, 0.3, nx
        )
    elif case_num == 5:
        sweep_params["electricity_cost"] = LinearSample(
            m.fs.costing.electricity_cost, 0.02, 0.3, nx
        )
        sweep_params["MgCl2_cost"] = LinearSample(
            m.fs.costing.magnesium_chloride_cost, 0.06, 0.1, nx
        )
    elif case_num == 6:
        sweep_params["energy_consumption"] = LinearSample(
            m.fs.electroNP.energy_electric_flow_mass, 0, 2.5, nx
        )
        sweep_params["MgCl2_dosage"] = LinearSample(
            m.fs.electroNP.magnesium_chloride_dosage, 0.1, 0.5, nx
        )
    elif case_num == 7:
        sweep_params["S_IP_concentration"] = LinearSample(
            m.fs.AD.inlet.conc_mass_comp[0, "S_IP"], 0.01, 0.5, nx
        )
    elif case_num == 8:
        sweep_params["S_IP_concentration"] = LinearSample(
            m.fs.AD.inlet.conc_mass_comp[0, "S_IP"], 0.01, 0.5, nx
        )
        sweep_params["X_PHA_concentration"] = LinearSample(
            m.fs.AD.inlet.conc_mass_comp[0, "X_PHA"], 0, 0.5, nx
        )
    elif case_num == 9:
        sweep_params["X_PAO_concentration"] = LinearSample(
            m.fs.AD.inlet.conc_mass_comp[0, "X_PAO"], 3.2, 3.8, nx
        )
        sweep_params["X_PHA_concentration"] = LinearSample(
            m.fs.AD.inlet.conc_mass_comp[0, "X_PHA"], 1e-6, 0.5, nx
        )
    elif case_num == 10:
        sweep_params["X_PAO_concentration"] = LinearSample(
            m.fs.AD.inlet.conc_mass_comp[0, "X_PHA"], 0, 0.5, nx
        )
        sweep_params["X_PP_concentration"] = LinearSample(
            m.fs.AD.inlet.conc_mass_comp[0, "X_PP"], 0.8, 1.2, nx
        )
    else:
        raise ValueError(f"{case_num} is not yet implemented")

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


if __name__ == "__main__":
    results, sweep_params, m = run_analysis()
