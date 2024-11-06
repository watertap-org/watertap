#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
from idaes_flowsheet_processor.api import FlowsheetInterface
from watertap.flowsheets.mvc.mvc_single_stage import (
    build,
    set_operating_conditions,
    add_Q_ext,
    initialize_system,
    scale_costs,
    fix_outlet_pressures,
    solve,
    set_up_optimization,
)
from pyomo.environ import units as pyunits, assert_optimal_termination, Objective
from watertap.core.solvers import get_solver


def export_to_ui():
    return FlowsheetInterface(
        name="Mechanical vapor compression",
        do_export=export_variables,
        do_build=build_flowsheet,
        do_solve=solve_flowsheet,
    )


def export_variables(flowsheet, exports=None, build_options=None, **kwargs):
    fs = flowsheet
    # --- Input data ---
    exports.add(
        obj=fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"],
        name="Feed mass flow",
        ui_units=pyunits.kg / pyunits.s,
        display_units="kg/s",
        rounding=2,
        description="Liquid feed mass flow",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.properties[0].mass_frac_phase_comp["Liq", "TDS"],
        name="TDS mass fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Salt mass fraction",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.properties[0].temperature,
        name="Feed temperature",
        ui_units=pyunits.K,
        display_units="K",
        rounding=2,
        description="Feed temperature",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.properties[0].pressure,
        name="Feed pressure",
        ui_units=pyunits.Pa,
        display_units="Pa",
        rounding=2,
        description="Feed pressure",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.recovery[0],
        name="Volumetric recovery",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Volumetric recovery of the flowsheet",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )

    exports.add(
        obj=fs.pump_feed.efficiency_pump[0],
        name="Feed pump efficiency",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Feed pump efficiency",
        is_input=True,
        input_category="Feed pump",
        is_output=False,
    )
    exports.add(
        obj=fs.pump_feed.control_volume.deltaP[0],
        name="Feed pump pressure change",
        ui_units=pyunits.Pa,
        display_units="Pa",
        rounding=2,
        description="Feed pump pressure change",
        is_input=True,
        input_category="Feed pump",
        is_output=True,
        output_category="Feed pump",
    )

    exports.add(
        obj=fs.hx_distillate.overall_heat_transfer_coefficient[0],
        name="Distillate HEX heat transfer coefficient",
        ui_units=pyunits.J * pyunits.s**-1 * pyunits.m**-2 * pyunits.K**-1,
        display_units="W/K-m2",
        rounding=2,
        description="Distillate heat exchanger heat transfer coefficient",
        is_input=True,
        input_category="Distillate heat exchanger",
        is_output=True,
        output_category="Distillate heat exchanger",
    )
    exports.add(
        obj=fs.hx_distillate.area,
        name="Distillate HEX area",
        ui_units=pyunits.m**2,
        display_units="m2",
        rounding=2,
        description="Distillate heat exchanger area",
        is_input=True,
        input_category="Distillate heat exchanger",
        is_output=True,
        output_category="Distillate heat exchanger",
    )
    exports.add(
        obj=fs.hx_distillate.cold.deltaP[0],
        name="Distillate HEX cold side pressure change",
        ui_units=pyunits.Pa,
        display_units="Pa",
        rounding=2,
        description="Distillate heat exchanger cold side pressure change",
        is_input=True,
        input_category="Distillate heat exchanger",
        is_output=True,
        output_category="Distillate heat exchanger",
    )
    exports.add(
        obj=fs.hx_distillate.hot.deltaP[0],
        name="Distillate HEX hot side pressure change",
        ui_units=pyunits.Pa,
        display_units="Pa",
        rounding=2,
        description="Distillate heat exchanger hot side pressure change",
        is_input=True,
        input_category="Distillate heat exchanger",
        is_output=True,
        output_category="Distillate heat exchanger",
    )

    exports.add(
        obj=fs.hx_brine.overall_heat_transfer_coefficient[0],
        name="Brine HEX heat transfer coefficient",
        ui_units=pyunits.J * pyunits.s**-1 * pyunits.m**-2 * pyunits.K**-1,
        display_units="W/K-m2",
        rounding=2,
        description="Brine heat exchanger heat transfer coefficient",
        is_input=True,
        input_category="Brine heat exchanger",
        is_output=True,
        output_category="Brine heat exchanger",
    )
    exports.add(
        obj=fs.hx_brine.area,
        name="Brine HEX area",
        ui_units=pyunits.m**2,
        display_units="m2",
        rounding=2,
        description="Brine heat exchanger area",
        is_input=True,
        input_category="Brine heat exchanger",
        is_output=True,
        output_category="Brine heat exchanger",
    )
    exports.add(
        obj=fs.hx_brine.cold.deltaP[0],
        name="Brine HEX cold side pressure change",
        ui_units=pyunits.Pa,
        display_units="Pa",
        rounding=2,
        description="Brine heat exchanger cold side pressure change",
        is_input=True,
        input_category="Brine heat exchanger",
        is_output=True,
        output_category="Brine heat exchanger",
    )
    exports.add(
        obj=fs.hx_brine.hot.deltaP[0],
        name="Brine HEX hot side pressure change",
        ui_units=pyunits.Pa,
        display_units="Pa",
        rounding=2,
        description="Brine heat exchanger hot side pressure change",
        is_input=True,
        input_category="Brine heat exchanger",
        is_output=True,
        output_category="Brine heat exchanger",
    )
    exports.add(
        obj=fs.evaporator.area,
        name="Evaporator area",
        ui_units=pyunits.m**2,
        display_units="m2",
        rounding=2,
        description="Evaporator area",
        is_input=True,
        input_category="Evaporator",
        is_output=True,
        output_category="Evaporator",
    )
    exports.add(
        obj=fs.evaporator.inlet_feed.temperature[0],
        name="Evaporator inlet temperature",
        ui_units=pyunits.K,
        display_units="K",
        rounding=2,
        description="Evaporator inlet feed temperature",
        is_input=True,
        input_category="Evaporator",
        is_output=True,
        output_category="Evaporator",
    )
    exports.add(
        obj=fs.evaporator.outlet_brine.temperature[0],
        name="Evaporator brine temperature",
        ui_units=pyunits.K,
        display_units="K",
        rounding=2,
        description="Evaporator outlet brine temperature",
        is_input=True,
        input_category="Evaporator",
        is_output=True,
        output_category="Evaporator",
    )
    exports.add(
        obj=fs.evaporator.U,
        name="Evaporator heat transfer coefficient",
        ui_units=pyunits.J * pyunits.s**-1 * pyunits.m**-2 * pyunits.K**-1,
        display_units="W/K-m2",
        rounding=2,
        description="Evaporator heat transfer coefficient",
        is_input=True,
        input_category="Evaporator",
        is_output=False,
    )
    exports.add(
        obj=fs.compressor.pressure_ratio,
        name="Compressor pressure ratio",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Compressor pressure ratio",
        is_input=True,
        input_category="Compressor",
        is_output=True,
        output_category="Compressor",
    )
    exports.add(
        obj=fs.compressor.efficiency,
        name="Compressor efficiency",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Compressor efficiency",
        is_input=True,
        input_category="Compressor",
        is_output=False,
    )
    exports.add(
        obj=fs.pump_brine.efficiency_pump[0],
        name="Brine pump efficiency",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Brine pump efficiency",
        is_input=True,
        input_category="Brine pump",
        is_output=False,
    )
    exports.add(
        obj=fs.pump_brine.control_volume.deltaP[0],
        name="Brine pump pressure change",
        ui_units=pyunits.Pa,
        display_units="Pa",
        rounding=2,
        description="Brine pump pressure change",
        is_input=True,
        input_category="Brine pump",
        is_output=True,
        output_category="Brine pump",
    )
    exports.add(
        obj=fs.pump_distillate.efficiency_pump[0],
        name="Distillate pump efficiency",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Distillate pump efficiency",
        is_input=True,
        input_category="Distillate pump",
        is_output=False,
    )
    exports.add(
        obj=fs.pump_distillate.control_volume.deltaP[0],
        name="Distillate pump pressure change",
        ui_units=pyunits.Pa,
        display_units="Pa",
        rounding=2,
        description="Distillate pump pressure change",
        is_input=True,
        input_category="Distillate pump",
        is_output=True,
        output_category="Distillate pump",
    )

    # System costing
    exports.add(
        obj=fs.costing.utilization_factor,
        name="Utilization factor",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Utilization factor - [annual use hours/total hours in year]",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.TIC,
        name="Practical investment factor",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=1,
        description="Practical investment factor - [total investment cost/direct "
        "capital costs]",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.electricity_cost,
        name="Electricity cost",
        ui_units=fs.costing.base_currency / pyunits.kWh,
        display_units="$/kWh",
        rounding=3,
        description="Electricity cost",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.total_investment_factor,
        name="Total investment factor",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Utilization factor - [annual use hours/total hours in year]",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.heat_exchanger.material_factor_cost,
        name="Heat exchanger material cost factor",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Heat exchanger material cost factor",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.evaporator.material_factor_cost,
        name="Evaporator material cost factor",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Evaporator material cost factor",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.compressor.unit_cost,
        name="Compressor unit cost",
        ui_units=pyunits.USD_2001,
        display_units="$(2001)",
        rounding=2,
        description="Compressor unit cost",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )

    # Outlets

    feed_salinity = pyunits.convert(
        fs.feed.properties[0].mass_frac_phase_comp["Liq", "TDS"],
        to_units=pyunits.g / pyunits.kg,
    )
    exports.add(
        obj=feed_salinity,
        name="Feed salinity",
        ui_units=pyunits.g / pyunits.kg,
        display_units="g/kg",
        rounding=2,
        description="Feed salinity",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    brine_salinity = pyunits.convert(
        fs.evaporator.properties_brine[0].mass_frac_phase_comp["Liq", "TDS"],
        to_units=pyunits.g / pyunits.kg,
    )
    exports.add(
        obj=brine_salinity,
        name="Brine salinity",
        ui_units=pyunits.g / pyunits.kg,
        display_units="g/kg",
        rounding=2,
        description="Brine salinity",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.separator_feed.split_fraction[0, "hx_distillate_cold"],
        name="Separator split fraction for distillate HEX",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Separator split fraction going into distillate heat exchanger",
        is_input=False,
        is_output=True,
        output_category="Separator",
    )
    exports.add(
        obj=fs.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"],
        name="Product flow rate",
        ui_units=pyunits.kg / pyunits.s,
        display_units="kg/s",
        rounding=2,
        description="Evaporator steam flow rate",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.compressor.control_volume.properties_out[0].temperature,
        name="Compressed vapor temperature",
        ui_units=pyunits.K,
        display_units="K",
        rounding=2,
        description="Compressed vapor temperature",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.compressor.control_volume.properties_out[0].pressure,
        name="Compressed vapor pressure",
        ui_units=pyunits.Pa,
        display_units="Pa",
        rounding=2,
        description="Compressed vapor pressure",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.condenser.control_volume.properties_out[0].temperature,
        name="Condensed vapor temperature",
        ui_units=pyunits.K,
        display_units="K",
        rounding=2,
        description="Condensed vapor temperature",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.evaporator.lmtd,
        name="Evaporator log-mean temperature difference",
        ui_units=pyunits.K,
        display_units="K",
        rounding=2,
        description="Evaporator LMTD",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )

    # System metrics
    exports.add(
        obj=fs.costing.LCOW,
        name="LCOW",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of centrate",
        rounding=3,
        description="Levelized cost of water including operating and capital costs",
        is_input=False,
        is_output=True,
        output_category="Levelized cost metrics",
    )
    # Normalized metrics
    total_capital_norm = fs.costing.total_capital_cost / fs.feed.properties[0].flow_vol
    exports.add(
        obj=total_capital_norm,
        name="Total capital",
        ui_units=fs.costing.base_currency / (pyunits.m**3 / pyunits.day),
        display_units="$/(m3/day)",
        rounding=1,
        description="Normalized total capital costs accounting for indirect "
        "capital and installation - [total capital costs/feed flow rate]",
        is_input=False,
        is_output=True,
        output_category="Normalized metrics",
    )

    # performance metrics
    exports.add(
        obj=fs.costing.specific_energy_consumption,
        name="Specific energy consumption",
        ui_units=pyunits.kWh / pyunits.m**3,
        display_units="kWh/m3",
        rounding=2,
        description="Specific energy consumption",
        is_input=False,
        is_output=True,
        output_category="Normalized metrics",
    )

    # Capital costs
    exports.add(
        obj=fs.costing.total_capital_cost,
        name="Total",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="Total capital costs - including investment factor to account "
        "for indirect capital",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )
    exports.add(
        obj=fs.hx_distillate.costing.capital_cost,
        name="Distillate heat exchanger",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="Distillate heat exchanger capital cost",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )
    exports.add(
        obj=fs.hx_brine.costing.capital_cost,
        name="Brine heat exchanger",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="Brine heat exchanger capital cost",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )
    exports.add(
        obj=fs.mixer_feed.costing.capital_cost,
        name="Feed mixer",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="Feed mixer capital cost",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )
    exports.add(
        obj=fs.evaporator.costing.capital_cost,
        name="Evaporator",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="Evaporator capital cost",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )
    exports.add(
        obj=fs.compressor.costing.capital_cost,
        name="Compressor",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="Compressor capital cost",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )
    # Operating costs
    exports.add(
        obj=fs.costing.total_operating_cost,
        name="Total",
        ui_units=fs.costing.base_currency / pyunits.year,
        display_units="$/year",
        rounding=0,
        description="Total annual operating costs - including electricity, heating, "
        "and fixed",
        is_input=False,
        is_output=True,
        output_category="Operating costs",
    )
    exports.add(
        obj=fs.costing.aggregate_flow_costs["electricity"],
        name="Electricity",
        ui_units=fs.costing.base_currency / pyunits.year,
        display_units="$/year",
        rounding=0,
        description="Annual electricity costs ",
        is_input=False,
        is_output=True,
        output_category="Operating costs",
    )


def build_flowsheet(build_options=None, **kwargs):
    # build and solve initial flowsheet
    m = build()
    set_operating_conditions(m)
    add_Q_ext(m, time_point=m.fs.config.time)
    initialize_system(m)
    # rescale costs after initialization because scaling depends on flow rates
    scale_costs(m)
    fix_outlet_pressures(m)  # outlet pressure are initially unfixed for initialization

    m.fs.objective = Objective(expr=m.fs.Q_ext[0])

    solver = get_solver()
    results = solve(m, solver=solver, tee=False)
    assert_optimal_termination(results)

    m.fs.Q_ext[0].fix(0)  # no longer want external heating in evaporator
    del m.fs.objective
    set_up_optimization(m)
    results = solve(m, solver=solver, tee=False)

    return m


def solve_flowsheet(flowsheet=None):
    fs = flowsheet
    results = solve(fs)
    return results
