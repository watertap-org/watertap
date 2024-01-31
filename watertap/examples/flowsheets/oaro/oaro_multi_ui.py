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
from idaes.core.solvers import get_solver
from watertap.ui.fsapi import FlowsheetInterface
from watertap.examples.flowsheets.oaro.oaro_multi import (
    build,
    set_operating_conditions,
    initialize_system,
    optimize_set_up,
    solve,
    ERDtype,
)
from pyomo.environ import units as pyunits


def export_to_ui():
    return FlowsheetInterface(
        name="OARO",
        do_export=export_variables,
        do_build=build_flowsheet,
        do_solve=solve_flowsheet,
    )


def export_variables(flowsheet=None, exports=None, build_options=None, **kwargs):
    fs = flowsheet
    # --- Input data ---
    # Feed conditions
    exports.add(
        obj=fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"],
        name="Water mass flowrate",
        ui_units=pyunits.kg / pyunits.s,
        display_units="kg/s",
        rounding=3,
        description="Inlet water mass flowrate",
        is_input=True,
        input_category="Feed",
        is_output=False,
    )


def build_flowsheet(
    number_of_stages=3, system_recovery=0.5, build_options=None, **kwargs
):
    # get solver
    solver = get_solver()

    # build, set, and initialize
    m = build(number_of_stages=number_of_stages, erd_type=ERDtype.pump_as_turbine)
    set_operating_conditions(m)
    initialize_system(
        m,
        number_of_stages,
        solvent_multiplier=0.5,
        solute_multiplier=0.7,
        solver=solver,
    )

    optimize_set_up(
        m, number_of_stages=number_of_stages, water_recovery=system_recovery
    )

    # display
    solve(m, solver=solver)
    return m


def solve_flowsheet(flowsheet=None):
    fs = flowsheet
    results = solve(fs)
    return results
