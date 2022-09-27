import sys
import os
import time

from pyomo.environ import Constraint
from watertap.tools.parameter_sweep import LinearSample, parameter_sweep
import watertap.tools.MPI as MPI
import watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.metab.metab as metab


def set_up_sensitivity(m):
    outputs = {}
    optimize_kwargs = {"check_termination": False}  # None
    opt_function = metab.solve

    # tie parameters together
    m.fs.costing.metab.eq_bead_cost = Constraint(
        expr=m.fs.costing.metab.bead_cost["hydrogen"]
        == m.fs.costing.metab.bead_cost["methane"]
    )
    m.fs.costing.metab.bead_cost["methane"].unfix()

    m.fs.costing.metab.eq_bead_replacement_factor = Constraint(
        expr=m.fs.costing.metab.bead_replacement_factor["hydrogen"]
        == m.fs.costing.metab.bead_replacement_factor["methane"]
    )
    m.fs.costing.metab.bead_replacement_factor["methane"].unfix()

    # new baseline parameters
    m.fs.costing.metab.bead_cost["hydrogen"].fix(14.4)
    # m.fs.costing.metab.hydraulic_retention_time['hydrogen'].fix()
    # m.fs.costing.metab.hydraulic_retention_time['methane'].fix()

    # create outputs
    outputs["LCOW"] = m.fs.costing.LCOW
    outputs["LCOH"] = m.fs.costing.LCOH
    outputs["LCOM"] = m.fs.costing.LCOM
    # outputs["Effluent Volumetric Flow Rate"] = m.fs.product_H2O.properties[0].flow_vol

    return outputs, optimize_kwargs, opt_function


def run_analysis(case_num, nx, interpolate_nan_outputs=True):
    m, _ = metab.main()

    outputs, optimize_kwargs, opt_function = set_up_sensitivity(m)

    sweep_params = {}
    if case_num == 1:
        # bead cost
        sweep_params["bead_cost"] = LinearSample(
            m.fs.costing.metab.bead_cost["hydrogen"], 1, 50, nx
        )

    elif case_num == 2:
        # bead replacement rate
        # baseline corresponds to replacement rate of 0.3 years; sensitivity on replacement_factor corresponding
        # to 0.3 yr to 5 yr
        import numpy as np
        from pyomo.environ import value

        upper_replacement_rate = 5  # replace every x years
        replacement_intervals = np.arange(
            upper_replacement_rate,
            value(m.fs.costing.plant_lifetime),
            upper_replacement_rate,
        )
        rep_factor_upper_lim = sum(
            1 / (1 + value(m.fs.costing.wacc)) ** x for x in replacement_intervals
        ) * value(m.fs.costing.capital_recovery_factor)

        sweep_params["bead_cost"] = LinearSample(
            m.fs.costing.metab.bead_replacement_factor["hydrogen"],
            3.376,
            rep_factor_upper_lim,
            nx,
        )

    elif case_num == 3:
        # Hydrogen METAB HRT
        sweep_params["hydrogen_hrt"] = LinearSample(
            m.fs.metab_hydrogen.hydraulic_retention_time, 0.75, 24, nx
        )

    elif case_num == 4:
        # Methane METAB HRT
        sweep_params["methane_hrt"] = LinearSample(
            m.fs.metab_methane.hydraulic_retention_time, 47.25, 360, nx
        )

    elif case_num == 5:
        # Hydrogen Conversion Rate: sweep from 0.06 to 0.6 L H2/g-COD removed
        sweep_params["hydrogen_conversion_rate"] = LinearSample(
            m.fs.metab_hydrogen.generation_ratio["cod_to_hydrogen", "hydrogen"],
            5.03e-3,
            5e-2,
            nx,
        )

    elif case_num == 6:
        # Methane Conversion Rate: sweep from ~ 0.1 to 0.3 L H2/g-COD removed
        sweep_params["methane_conversion_rate"] = LinearSample(
            m.fs.metab_methane.generation_ratio["cod_to_methane", "methane"],
            6.68e-2,
            2e-1,
            nx,
        )
    elif case_num == 7:
        # Recovered resource selling price:
        sweep_params["Hydrogen selling price"] = LinearSample(
            m.fs.costing.hydrogen_product_cost, -2, -10, nx
        )
        sweep_params["Methane selling price"] = LinearSample(
            m.fs.costing.methane_product_cost, -0.305, -1, nx
        )

    else:
        raise ValueError("case_num = %d not recognized." % (case_num))

    output_filename = "sensitivity_" + str(case_num) + ".csv"

    global_results = parameter_sweep(
        m,
        sweep_params,
        outputs,
        csv_results_file_name=output_filename,
        optimize_function=opt_function,
        optimize_kwargs=optimize_kwargs,
        # debugging_data_dir=os.path.split(output_filename)[0] + '/local',
        interpolate_nan_outputs=interpolate_nan_outputs,
    )

    return global_results, sweep_params


def main(case_num=1, nx=11, interpolate_nan_outputs=True):

    # when from the command line
    case_num = int(case_num)
    nx = int(nx)
    interpolate_nan_outputs = bool(interpolate_nan_outputs)

    tic = time.time()
    global_results, sweep_params = run_analysis(case_num, nx, interpolate_nan_outputs)
    print(global_results)
    toc = time.time()

    if MPI.COMM_WORLD.rank == 0:

        total_samples = 1

        for k, v in sweep_params.items():
            total_samples *= v.num_samples

        print("Finished case_num = %d." % (case_num))
        print(
            "Processed %d swept parameters comprising %d total points."
            % (len(sweep_params), total_samples)
        )
        print("Elapsed time = %.1f s." % (toc - tic))

    return global_results, sweep_params


if __name__ == "__main__":
    global_results, sweep_params = main(*sys.argv[1:])
