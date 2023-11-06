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
from watertap.tools.parameter_sweep import LinearSample, parameter_sweep
import watertap.examples.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery as RO
from watertap.examples.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery import (
    ERDtype,
)


def set_up_sensitivity():
    outputs = {}

    m = RO.build(erd_type=ERDtype.no_ERD)
    RO.set_operating_conditions(
        m,
        water_recovery=0.7,
        over_pressure=0.3,
        flow_vol=1e-3,
        salt_mass_conc=5e-3,
    )
    RO.initialize_system(m)
    RO.solve(m)
    m.fs.feed.properties[0].flow_mass_phase_comp.unfix()
    m.fs.feed.properties[0].flow_vol_phase["Liq"].fix()
    m.fs.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].fix()
    RO.optimize_set_up(m)
    RO.solve(m)
    print("\n***---Optimization results---***")
    RO.display_system(m)
    RO.display_design(m)
    # create outputs
    outputs["LCOW"] = m.fs.costing.LCOW
    # Not used in Excel
    outputs["RO Energy Consumption"] = m.fs.costing.specific_energy_consumption
    # No conversion in Excel
    outputs["System Capital Cost"] = m.fs.costing.aggregate_capital_cost
    # Multiplied by flow rate scaling factor in Excel
    outputs["RO Operating Cost"] = m.fs.RO.costing.fixed_operating_cost
    # Multiplied by flow rate scaling factor in Excel
    outputs[
        "MLC Operating Cost"
    ] = m.fs.costing.maintenance_labor_chemical_operating_cost
    # Multiplied by flow rate scaling factor in Excel
    outputs["Feed Flow Rate"] = m.fs.feed.properties[0].flow_vol_phase["Liq"]
    # Fixed based on desired flow rates given by the dye desal team
    outputs["Permeate Flow Rate"] = m.fs.product.properties[0].flow_vol_phase["Liq"]
    # Multiplied by flow rate scaling factor in Excel
    outputs["Retentate Flow Rate"] = m.fs.disposal.properties[0].flow_vol_phase["Liq"]
    # Multiplied by flow rate scaling factor in Excel
    outputs["RO Operating Pressure"] = m.fs.RO.inlet.pressure[0]
    # No conversion in Excel
    outputs["RO Permeate H2O Mass Flow"] = m.fs.RO.permeate.flow_mass_phase_comp[
        0, "Liq", "H2O"
    ]
    # Multiplied by flow rate scaling factor in Excel
    outputs["RO Permeate Salt Mass Flow"] = m.fs.RO.permeate.flow_mass_phase_comp[
        0, "Liq", "NaCl"
    ]
    # Multiplied by flow rate scaling factor in Excel
    outputs["RO Retentate H2O Mass Flow"] = m.fs.RO.retentate.flow_mass_phase_comp[
        0, "Liq", "H2O"
    ]
    # Multiplied by flow rate scaling factor in Excel
    outputs["RO Retentate Salt Mass Flow"] = m.fs.RO.retentate.flow_mass_phase_comp[
        0, "Liq", "NaCl"
    ]
    # Multiplied by flow rate scaling factor in Excel

    return outputs, m


def run_analysis(
    case_num=1, nx=3, interpolate_nan_outputs=True, withERD=True, output_filename=None
):

    if output_filename is None:
        output_filename = "sensitivity_" + str(case_num) + ".csv"

    # when from the command line
    case_num = int(case_num)
    nx = int(nx)
    interpolate_nan_outputs = bool(interpolate_nan_outputs)
    # withERD = bool(withERD)

    # # select flowsheet configuration
    # if withERD:
    #     m = RO.main(erd_type=ERDtype.pump_as_turbine)[0]
    # else:
    #     m = RO.main(erd_type=ERDtype.pressure_exchanger)[0]

    outputs, m = set_up_sensitivity()

    # choose parameter sweep from case structure
    sweep_params = {}

    if case_num == 1:
        # Need to unfix mass recovery of water (or simply sweep across it instead of recovery_vol)
        m.fs.RO.recovery_mass_phase_comp.unfix()

        sweep_params["mass_concentration"] = LinearSample(
            m.fs.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"],
            0.963,
            4.816,
            nx,
        )
        sweep_params["volumetric_recovery"] = LinearSample(
            m.fs.RO.recovery_vol_phase[0, "Liq"], 0.7, 0.9, nx
        )
    else:
        raise ValueError(f"{case_num} is not yet implemented")

    global_results = parameter_sweep(
        m,
        sweep_params,
        outputs,
        csv_results_file_name=output_filename,
        optimize_function=RO.solve,
        # optimize_kwargs=optimize_kwargs,
        interpolate_nan_outputs=interpolate_nan_outputs,
    )

    return global_results, sweep_params, m


if __name__ == "__main__":
    results, sweep_params, m = run_analysis()
    print(results)
