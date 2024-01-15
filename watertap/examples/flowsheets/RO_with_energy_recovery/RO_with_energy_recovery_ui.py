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
from watertap.examples.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery import (
    build,
    set_operating_conditions,
    initialize_system,
    solve,
    ERDtype,
)
from pyomo.environ import units as pyunits


def export_to_ui():
    return FlowsheetInterface(
        name="RO with energy recovery",
        do_export=export_variables,
        do_build=build_flowsheet,
        do_solve=solve_flowsheet,
    )


def export_variables(model=None, exports=None, build_options=None, **kwargs):
    # --- Input data ---
    # Feed conditions
    exports.add(
        obj=model.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"],
        name="Water mass flowrate",
        ui_units=pyunits.kg / pyunits.s,
        display_units="kg/s",
        rounding=3,
        description="Inlet water mass flowrate",
        is_input=True,
        input_category="Feed",
        is_output=False,
    )
    exports.add(
        obj=model.fs.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"],
        name="NaCl mass flowrate",
        ui_units=pyunits.kg / pyunits.s,
        display_units="kg/s",
        rounding=3,
        description="Inlet NaCl mass flowrate",
        is_input=True,
        input_category="Feed",
        is_output=False,
    )

    # Unit model data, feed pump
    exports.add(
        obj=model.fs.P1.efficiency_pump[0],
        name="Pump efficiency",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Efficiency of feed pump",
        is_input=True,
        input_category="Feed Pump",
        is_output=False,
    )
    exports.add(
        obj=model.fs.P1.control_volume.properties_out[0].pressure,
        name="Operating pressure",
        ui_units=pyunits.bar,
        display_units="bar",
        rounding=2,
        description="Operating pressure of feed pump",
        is_input=True,
        input_category="Feed Pump",
        is_output=False,
    )

    # Unit model data, RO
    exports.add(
        obj=model.fs.RO.A_comp[0, "H2O"],
        name="Water permeability coefficient",
        ui_units=pyunits.L / pyunits.hr / pyunits.m**2 / pyunits.bar,
        display_units="LMH/bar",
        rounding=2,
        description="Water permeability coefficient",
        is_input=True,
        input_category="Reverse Osmosis",
        is_output=False,
    )

    exports.add(
        obj=model.fs.RO.B_comp[0, "NaCl"],
        name="Salt permeability coefficient",
        ui_units=pyunits.L / pyunits.hr / pyunits.m**2,
        display_units="LMH",
        rounding=2,
        description="Salt permeability coefficient",
        is_input=True,
        input_category="Reverse Osmosis",
        is_output=False,
    )
    exports.add(
        obj=model.fs.RO.feed_side.channel_height,
        name="Feed-side channel height",
        ui_units=pyunits.mm,
        display_units="mm",
        rounding=2,
        description="Feed-side channel height",
        is_input=True,
        input_category="Reverse Osmosis",
        is_output=False,
    )
    exports.add(
        obj=model.fs.RO.feed_side.spacer_porosity,
        name="Feed-side spacer porosity",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Feed-side spacer porosity",
        is_input=True,
        input_category="Reverse Osmosis",
        is_output=False,
    )
    exports.add(
        obj=model.fs.RO.permeate.pressure[0],
        name="Permeate-side pressure",
        ui_units=pyunits.bar,
        display_units="bar",
        rounding=2,
        description="Permeate-side pressure",
        is_input=True,
        input_category="Reverse Osmosis",
        is_output=False,
    )
    exports.add(
        obj=model.fs.RO.width,
        name="Width",
        ui_units=pyunits.m,
        display_units="m",
        rounding=2,
        description="Stage width",
        is_input=True,
        input_category="Reverse Osmosis",
        is_output=False,
    )
    exports.add(
        obj=model.fs.RO.recovery_mass_phase_comp[0, "Liq", "H2O"],
        name="Water mass recovery",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Water mass recovery of RO unit",
        is_input=True,
        input_category="Reverse Osmosis",
        is_output=True,
        output_category="System metrics",
    )

    # Unit model data, ERD
    exports.add(
        obj=model.fs.ERD.efficiency_pump[0],
        name="Pump efficiency",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Efficiency of energy recovery device",
        is_input=True,
        input_category="Energy Recovery Device",
        is_output=False,
    )
    exports.add(
        obj=model.fs.ERD.control_volume.properties_out[0].pressure,
        name="Operating pressure",
        ui_units=pyunits.bar,
        display_units="bar",
        rounding=2,
        description="Operating pressure of energy recovery device",
        is_input=True,
        input_category="Energy Recovery Device",
        is_output=False,
    )

    # System costing
    exports.add(
        obj=model.fs.costing.utilization_factor,
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
        obj=model.fs.costing.TIC,
        name="Total Installed Cost",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=1,
        description="Total Installed Cost (TIC)",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=model.fs.costing.TPEC,
        name="Total Purchased Equipment Cost",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=1,
        description="Total Purchased Equipment Cost (TPEC)",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=model.fs.costing.factor_total_investment,
        name="Total investment factor",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Total investment factor [investment cost/equipment cost]",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=model.fs.costing.factor_maintenance_labor_chemical,
        name="Maintenance-labor-chemical factor",
        ui_units=1 / pyunits.year,
        display_units="fraction/year",
        rounding=2,
        description="Maintenance-labor-chemical factor [fraction of investment cost/year]",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=model.fs.costing.factor_capital_annualization,
        name="Capital annualization factor",
        ui_units=1 / pyunits.year,
        display_units="fraction/year",
        rounding=2,
        description="Capital annualization factor [fraction of investment cost/year]",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=model.fs.costing.electricity_cost,
        name="Electricity cost",
        ui_units=model.fs.costing.base_currency / pyunits.kWh,
        display_units="$/kWh",
        rounding=3,
        description="Electricity cost",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )

    # Feed
    exports.add(
        obj=model.fs.feed.properties[0].flow_vol_phase["Liq"],
        name="Volumetric flow rate",
        ui_units=pyunits.m**3 / pyunits.hour,
        display_units="m3/hr",
        rounding=2,
        description="Inlet volumetric flow rate",
        is_input=False,
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=model.fs.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"],
        name="NaCl concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=2,
        description="Inlet NaCl concentration",
        is_input=False,
        is_output=True,
        output_category="Feed",
    )

    # Product
    exports.add(
        obj=model.fs.product.properties[0].flow_vol,
        name="Volumetric flow rate",
        ui_units=pyunits.m**3 / pyunits.hr,
        display_units="m3/h",
        rounding=2,
        description="Outlet product water volumetric flow rate",
        is_input=False,
        is_output=True,
        output_category="Product",
    )
    exports.add(
        obj=model.fs.product.properties[0].conc_mass_phase_comp["Liq", "NaCl"],
        name="NaCl concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=3,
        description="Outlet product water NaCl concentration",
        is_input=False,
        is_output=True,
        output_category="Product",
    )

    # Disposal
    exports.add(
        obj=model.fs.disposal.properties[0].flow_vol,
        name="Volumetric flow rate",
        ui_units=pyunits.m**3 / pyunits.hr,
        display_units="m3/h",
        rounding=2,
        description="Outlet product water volumetric flow rate",
        is_input=False,
        is_output=True,
        output_category="Disposal",
    )
    exports.add(
        obj=model.fs.disposal.properties[0].conc_mass_phase_comp["Liq", "NaCl"],
        name="NaCl concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=3,
        description="Outlet product water NaCl concentration",
        is_input=False,
        is_output=True,
        output_category="Disposal",
    )

    # System metrics
    exports.add(
        obj=model.fs.RO.area,
        name="Area",
        ui_units=pyunits.m**2,
        display_units="m^2",
        rounding=2,
        description="Membrane area",
        is_input=False,
        is_output=True,
        output_category="System metrics",
    )
    exports.add(
        obj=model.fs.costing.specific_energy_consumption,
        name="Specific energy consumption",
        ui_units=pyunits.kWh / pyunits.m**3,
        display_units="kWh/m3 of product water",
        rounding=3,
        description="Specific energy consumption (SEC)",
        is_input=False,
        is_output=True,
        output_category="System metrics",
    )
    exports.add(
        obj=model.fs.costing.LCOW,
        name="Levelized cost of water",
        ui_units=model.fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of product water",
        rounding=3,
        description="Levelized cost of water (LCOW)",
        is_input=False,
        is_output=True,
        output_category="System metrics",
    )


def build_flowsheet(erd_type=ERDtype.pump_as_turbine, build_options=None, **kwargs):
    # build and solve initial flowsheet
    m = build()

    solver = get_solver()

    # build, set, and initialize
    m = build(erd_type=erd_type)
    set_operating_conditions(m)
    initialize_system(m, solver=solver)

    # NOTE: GUI aims to solve simulation-based flowsheet
    # optimize_set_up(m)

    # display
    solve(m, solver=solver)
    return m


def solve_flowsheet(model=None):
    results = solve(model)
    return results
