#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
import logging

logging.getLogger("idaes.core.util.scaling").disabled = True
from parameter_sweep import LinearSample, ParameterSweep
import pyomo.environ as pyo
import watertap.flowsheets.electroNP.BSM2_genericNP_no_bioP as genericNP_flowsheet
from pyomo.environ import units as pyunits


def build_and_initialize_model():
    """Construct flowsheet with costing-first pattern used by the working test."""

    m = genericNP_flowsheet.build_flowsheet(has_genericNP=True, basis="mass")

    # Costing before conditions (critical) + deactivate capital cost constraints
    genericNP_flowsheet.add_costing(m)
    for c in m.fs.component_objects(pyo.Constraint, descend_into=True):
        if "capital_cost" in c.name:
            c.deactivate()

    # Set default operating conditions
    genericNP_flowsheet.set_operating_conditions(m)

    # Deactivate pressure equality constraints before init
    for mx in m.fs.mixers:
        mx.pressure_equality_constraints[0.0, 2].deactivate()
    m.fs.MX3.pressure_equality_constraints[0.0, 2].deactivate()
    m.fs.MX3.pressure_equality_constraints[0.0, 3].deactivate()

    # Initialize with costing present
    genericNP_flowsheet.initialize_system(m, has_genericNP=True)

    # Re-deactivate after init
    for mx in m.fs.mixers:
        mx.pressure_equality_constraints[0.0, 2].deactivate()
    m.fs.MX3.pressure_equality_constraints[0.0, 2].deactivate()
    m.fs.MX3.pressure_equality_constraints[0.0, 3].deactivate()

    # Solve once to establish a consistent starting point
    solver = genericNP_flowsheet.get_solver()
    solver.options["max_iter"] = 5000
    results = solver.solve(m, tee=False)
    pyo.assert_optimal_termination(results)

    return m


def build_outputs(model, **kwargs):
    """Build output dictionary from model. Must always return valid Pyomo components."""
    outputs = {}

    # parameter_sweep needs Pyomo components for outputs
    outputs["LCOW"] = model.fs.costing.LCOW
    outputs["GenericNP Capital Cost"] = model.fs.genericNP.costing.capital_cost
    outputs["Electricity Cost"] = model.fs.costing.aggregate_flow_costs["electricity"]
    outputs["S_PO4 Concentration"] = model.fs.genericNP.treated.conc_mass_comp[
        0, "S_PO4"
    ]
    outputs["S_NH4 Concentration"] = model.fs.genericNP.treated.conc_mass_comp[
        0, "S_NH4"
    ]
    outputs["NH4_removal"] = model.fs.genericNP.removal_factors["S_NH4"]
    outputs["P_removal"] = model.fs.genericNP.removal_factors["S_PO4"]

    return outputs


def build_sweep_params(model, case_num=1, nx=11, **kwargs):
    """Build sweep parameters that directly vary model variables."""
    sweep_params = {}

    if case_num == 1:
        # 1D: NH4 removal fraction sensitivity
        sweep_params["NH4_removal_fraction"] = LinearSample(
            model.fs.genericNP.removal_factors["S_NH4"], 0.1, 0.95, nx
        )
    elif case_num == 2:
        # 2D: NH4 removal fraction and energy intensity (heatmap case)
        sweep_params["NH4_removal_fraction"] = LinearSample(
            model.fs.genericNP.removal_factors["S_NH4"], 0.1, 0.95, nx
        )
        sweep_params["NH4_energy_intensity"] = LinearSample(
            model.fs.genericNP.energy_electric_flow["S_NH4"], 0.04, 2.5, nx
        )
    else:
        raise ValueError(f"Unsupported case_num: {case_num}")

    return sweep_params


def get_optimize_function(solver):
    """Solve wrapper that handles failures gracefully."""

    def optimize_function(m_inner):
        try:
            results = solver.solve(m_inner, tee=True)
            # Check if solve was successful
            if not pyo.check_optimal_termination(results):
                print(
                    f"Warning: Solver returned {results.solver.termination_condition}"
                )
            return results
        except Exception as e:
            print(f"Solve failed with exception: {e}")
            # Return a failed result object
            from pyomo.opt import SolverResults

            results = SolverResults()
            results.solver.termination_condition = "error"
            return results

    return optimize_function


def run_analysis(case_num=1, nx=11, interpolate_nan_outputs=True, output_filename=None):

    if output_filename is None:
        output_filename = f"genericnp_sensitivity_{case_num}.csv"

    solver = genericNP_flowsheet.get_solver()
    solver.options["max_iter"] = 5000
    opt_function = get_optimize_function(solver)

    ps = ParameterSweep(
        csv_results_file_name=output_filename,
        interpolate_nan_outputs=interpolate_nan_outputs,
        optimize_function=opt_function,
        optimize_kwargs={},
        initialize_before_sweep=False,
        reinitialize_before_sweep=False,
        number_of_subprocesses=1,
        parallel_back_end=None,
    )

    results_array, results_dict = ps.parameter_sweep(
        build_model=build_and_initialize_model,
        build_sweep_params=build_sweep_params,
        build_sweep_params_kwargs=dict(case_num=case_num, nx=nx),
        build_outputs=build_outputs,
        num_samples=nx,
    )

    return results_array, results_dict, None


if __name__ == "__main__":
    case_num = 1
    nx = 3

    print(f"Running GenericNP sensitivity case {case_num} with nx={nx}")
    results_array, results_dict, _ = run_analysis(case_num=case_num, nx=nx)
    print(f"Sweep complete! Results saved to genericnp_sensitivity_{case_num}.csv")
