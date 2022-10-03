###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################

import sys
import os
import time

from pyomo.environ import Constraint
from watertap.tools.parameter_sweep import LinearSample, parameter_sweep
import watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.dye_desalination.dye_desalination as dye_desalination
import watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.dye_desalination.dye_desalination_withRO as dye_desalination_withRO


def set_up_sensitivity(m, withRO):
    outputs = {}

    # LCOT is an output for both flowsheets
    outputs["LCOT"] = m.fs.LCOT

    # choose the right flowsheet and if ro is enabled add lcow
    if withRO:
        outputs["LCOW"] = m.fs.LCOW
        opt_function = dye_desalination_withRO.solve
    else:
        opt_function = dye_desalination.solve

    optimize_kwargs = {"check_termination": False}
    return outputs, optimize_kwargs, opt_function


def run_analysis(
    case_num=4, nx=11, interpolate_nan_outputs=True, withRO=True, save_path=None
):
    # when from the command line
    case_num = int(case_num)
    nx = int(nx)
    interpolate_nan_outputs = bool(interpolate_nan_outputs)
    withRO = bool(withRO)

    if not withRO and case_num > 7:
        raise ValueError(
            "Case numbers 8 and above are only for dye_desalination_withRO. Please set 'withRO=True'"
        )

    # select flowsheet
    if withRO:
        m = dye_desalination_withRO.main()[0]
    else:
        m = dye_desalination.main()[0]

    # set up sensitivities
    outputs, optimize_kwargs, opt_function = set_up_sensitivity(m, withRO)

    # choose parameter sweep from case structure
    sweep_params = {}

    if case_num == 1:
        m.fs.zo_costing.dye_mass_cost.unfix()
        sweep_params["dye_mass_cost"] = LinearSample(
            m.fs.zo_costing.dye_mass_cost, 0.1, 5, nx
        )

    elif case_num == 2:
        m.fs.zo_costing.waste_disposal_cost.unfix()
        sweep_params["disposal_cost"] = LinearSample(
            m.fs.zo_costing.waste_disposal_cost, 1, 10, nx
        )
    elif case_num == 3:
        m.fs.dye_separation.nanofiltration.removal_frac_mass_comp[0, "dye"].unfix()
        sweep_params["dye_removal"] = LinearSample(
            m.fs.dye_separation.nanofiltration.removal_frac_mass_comp[0, "dye"],
            0.2,
            0.999,
            nx,
        )
    elif case_num == 4:
        # sweep across membrane properties
        m.fs.dye_separation.nanofiltration.water_permeability_coefficient[0].unfix()
        m.fs.dye_separation.nanofiltration.removal_frac_mass_comp[0, "dye"].unfix()
        sweep_params["water_permeability"] = LinearSample(
            m.fs.dye_separation.nanofiltration.water_permeability_coefficient[0],
            1,
            100,
            nx,
        )
        sweep_params["dye_removal"] = LinearSample(
            m.fs.dye_separation.nanofiltration.removal_frac_mass_comp[0, "dye"],
            0.2,
            0.999,
            nx,
        )
    elif case_num == 5:
        # sweep across membrane cost parameters
        m.fs.zo_costing.nanofiltration.membrane_cost.unfix()

        sweep_params["membrane_cost"] = LinearSample(
            m.fs.zo_costing.nanofiltration.membrane_cost, 1, 100, nx
        )
    elif case_num == 6:
        m.fs.zo_costing.dye_mass_cost.unfix()
        m.fs.zo_costing.waste_disposal_cost.unfix()

        sweep_params["dye_cost"] = LinearSample(m.fs.zo_costing.dye_mass_cost, 0, 1, nx)
        sweep_params["waste_disposal"] = LinearSample(
            m.fs.zo_costing.waste_disposal_cost, 0, 10, nx
        )
    elif case_num == 7:
        m.fs.zo_costing.electricity_cost.unfix()
        sweep_params["electricity_cost"] = LinearSample(
            m.fs.zo_costing.electricity_cost, 0.0, 0.25, nx
        )
    elif case_num == 8:
        m.fs.zo_costing.recovered_water_cost.unfix()
        sweep_params["recovered_water_value"] = LinearSample(
            m.fs.zo_costing.recovered_water_cost, 0.0, 1.5, nx
        )
        m.fs.zo_costing.dye_mass_cost.unfix()
        sweep_params["recovered_dye_value"] = LinearSample(
            m.fs.zo_costing.dye_mass_cost, 0.0, 1, nx
        )
    elif case_num == 9:
        m.fs.zo_costing.recovered_water_cost.unfix()
        sweep_params["recovered_water_value"] = LinearSample(
            m.fs.zo_costing.recovered_water_cost, 0.0, 1.5, nx
        )
        m.fs.zo_costing.waste_disposal_cost.unfix()
        sweep_params["waste_disposal_cost"] = LinearSample(
            m.fs.zo_costing.waste_disposal_cost, 0.0, 10, nx
        )
    elif case_num == 10:
        m.fs.zo_costing.electricity_cost.unfix()
        m.fs.zo_costing.recovered_water_cost.unfix()
        sweep_params["recovered_water_value"] = LinearSample(
            m.fs.zo_costing.recovered_water_cost, 0.0, 1.5, nx
        )
        sweep_params["electricity_cost"] = LinearSample(
            m.fs.zo_costing.electricity_cost, 0.0, 0.25, nx
        )
    elif case_num == 11:
        m.fs.dye_separation.nanofiltration.removal_frac_mass_comp[0, "dye"].unfix()
        sweep_params["NF_dye_removal"] = LinearSample(
            m.fs.dye_separation.nanofiltration.removal_frac_mass_comp[0, "dye"],
            0.2,
            1,
            nx,
        )
        m.fs.desalination.RO.A_comp.unfix()
        sweep_params["RO_permeability"] = LinearSample(
            m.fs.desalination.RO.A_comp, 1e-12, 1e-11, nx
        )
    elif case_num == 12:
        desal = m.fs.desalination
        desal.RO.recovery_vol_phase[0, "Liq"].unfix()
        desal.RO.velocity[0, 0].unfix()

        desal.P2.control_volume.properties_out[0].pressure.unfix()
        desal.P2.control_volume.properties_out[0].pressure.setub(8300000)
        desal.P2.control_volume.properties_out[0].pressure.setlb(100000)

        desal.RO.area.unfix()
        desal.RO.area.setub(5000)
        desal.RO.area.setlb(50)

        sweep_params["RO_recovery"] = LinearSample(
            m.fs.desalination.RO.recovery_vol_phase[0, "Liq"], 0.1, 0.75, nx
        )
    else:
        raise ValueError("case_num = %d not recognized." % (case_num))

    # run sweep
    global_results = parameter_sweep(
        m,
        sweep_params,
        outputs,
        csv_results_file_name=save_path,
        optimize_function=opt_function,
        optimize_kwargs=optimize_kwargs,
        interpolate_nan_outputs=interpolate_nan_outputs,
    )

    return global_results, sweep_params, m


if __name__ == "__main__":
    results, params, model = run_analysis(*sys.argv[1:])
