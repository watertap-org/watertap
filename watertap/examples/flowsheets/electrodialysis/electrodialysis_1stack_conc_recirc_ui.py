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
from watertap.core.solvers import get_solver
from watertap.ui.fsapi import FlowsheetInterface
from watertap.examples.flowsheets.electrodialysis.electrodialysis_1stack_conc_recirc import (
    build,
    _condition_base,
    initialize_system,
    solve,
)
from pyomo.environ import units as pyunits


def export_to_ui():
    return FlowsheetInterface(
        name="Electrodialysis with concentrate recirculation",
        do_export=export_variables,
        do_build=build_flowsheet,
        do_solve=solve_flowsheet,
    )


def export_variables(flowsheet=None, exports=None, build_options=None, **kwargs):
    fs = flowsheet
    # --- Input data ---
    # Feed conditions
    exports.add(
        obj=fs.feed.properties[0].flow_vol_phase["Liq"],
        name="Feed solution volume flowrate",
        ui_units=pyunits.m ** 3/ pyunits.s,
        display_units="m^3/s",
        rounding=3,
        description="Inlet water mass flowrate",
        is_input=True,
        input_category="Feed",
        is_output=False,
    )
    exports.add(
        obj=fs.feed.properties[0].conc_mol_phase_comp["Liq", "Na_+"],
        name="Na_+ molar concentration",
        ui_units=pyunits.mol / pyunits.m **3,
        display_units="mol/m^3",
        rounding=3,
        description="Molcar concentration of Na_+ ions in the feed solution",
        is_input=True,
        input_category="Feed",
        is_output=False,
    )

    exports.add(
        obj=fs.feed.properties[0].conc_mol_phase_comp["Liq", "Cl_-"],
        name="Cl_- molar concentration",
        ui_units=pyunits.mol / pyunits.m **3,
        display_units="mol/m^3",
        rounding=3,
        description="Molcar concentration of Cl_- ions in the feed solution",
        is_input=True,
        input_category="Feed",
        is_output=False,
    )

    exports.add(
        obj=fs.EDstack.voltage_applied[0],
        name="Stack voltage",
        ui_units=pyunits.mol / pyunits.m **3,
        display_units="V",
        rounding=3,
        description="Applied constant voltage on the ED stack",
        is_input=True,
        input_category="ED stack",
        is_output=True,
    )

    exports.add(
        obj=fs.recovery_vol_H2O,
        name="Water recovery by volume",
        ui_units=pyunits.dimensionless,
        display_units="1",
        rounding=2,
        description="Product water recovery by volume",
        is_input=True,
        input_category="ED stack",
        is_output=True,
    )
    
    m.fs.feed.properties[0].pressure.fix(101325)
    m.fs.feed.properties[0].temperature.fix(298.15)
    m.fs.pump0.control_volume.properties_in[0].pressure.fix(101325)
    m.fs.pump1.efficiency_pump.fix(0.8)
    m.fs.pump0.efficiency_pump.fix(0.8)

    m.fs.prod.properties[0].pressure.fix(101325)
    m.fs.disp.properties[0].pressure.fix(101325)
    m.fs.disp.properties[0].temperature.fix(298.15)

    # Set ED unit vars
    # membrane properties
    m.fs.EDstack.water_trans_number_membrane["cem"].fix(5.8)
    m.fs.EDstack.water_trans_number_membrane["aem"].fix(4.3)
    m.fs.EDstack.water_permeability_membrane["cem"].fix(2.16e-14)
    m.fs.EDstack.water_permeability_membrane["aem"].fix(1.75e-14)
    m.fs.EDstack.membrane_areal_resistance["cem"].fix(1.89e-4)
    m.fs.EDstack.membrane_areal_resistance["aem"].fix(1.77e-4)
    m.fs.EDstack.solute_diffusivity_membrane["cem", "Na_+"].fix(3.28e-11)
    m.fs.EDstack.solute_diffusivity_membrane["aem", "Na_+"].fix(3.28e-11)
    m.fs.EDstack.solute_diffusivity_membrane["cem", "Cl_-"].fix(3.28e-11)
    m.fs.EDstack.solute_diffusivity_membrane["aem", "Cl_-"].fix(3.28e-11)
    m.fs.EDstack.ion_trans_number_membrane["cem", "Na_+"].fix(1)
    m.fs.EDstack.ion_trans_number_membrane["aem", "Na_+"].fix(0)
    m.fs.EDstack.ion_trans_number_membrane["cem", "Cl_-"].fix(0)
    m.fs.EDstack.ion_trans_number_membrane["aem", "Cl_-"].fix(1)
    m.fs.EDstack.membrane_thickness["aem"].fix(1.3e-4)
    m.fs.EDstack.membrane_thickness["cem"].fix(1.3e-4)

    # Stack properties
    m.fs.EDstack.cell_pair_num.fix(56)
    m.fs.EDstack.channel_height.fix(7.1e-4)
    m.fs.EDstack.cell_width.fix(0.197)
    m.fs.EDstack.cell_length.fix(1.68)

    # Spacer properties
    m.fs.EDstack.spacer_porosity.fix(0.83)
    m.fs.EDstack.spacer_specific_area.fix(10400)

    # Electrochemical properties
    m.fs.EDstack.electrodes_resistance.fix(0)
    m.fs.EDstack.current_utilization.fix(1)
    m.fs.EDstack.diffus_mass.fix(1.6e-9)
    # Unit model data, feed pump
    exports.add(
        obj=fs.pump0.efficiency_pump[0],
        name="Concentrate Pump efficiency",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Efficiency of concentrate feed pump",
        is_input=True,
        input_category="Feed Pump",
        is_output=False,
    )

    exports.add(
        obj=fs.pump1.efficiency_pump[0],
        name="Diluate Pump efficiency",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Efficiency of diluate feed pump",
        is_input=True,
        input_category="Feed Pump",
        is_output=False,
    )

    # exports.add(
    #     obj=fs.pump0.control_volume.properties_out[0].pressure,
    #     name="Concentrate Pump Operating pressure",
    #     ui_units=pyunits.bar,
    #     display_units="bar",
    #     rounding=2,
    #     description="Operating pressure of feed pump",
    #     is_input=True,
    #     input_category="Feed Pump",
    #     is_output=False,
    # )

    # Unit model data, RO
    exports.add(
        obj=fs.RO.A_comp[0, "H2O"],
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
        obj=fs.RO.B_comp[0, "NaCl"],
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
        obj=fs.RO.feed_side.channel_height,
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
        obj=fs.RO.feed_side.spacer_porosity,
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
        obj=fs.RO.permeate.pressure[0],
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
        obj=fs.RO.width,
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
        obj=fs.RO.recovery_mass_phase_comp[0, "Liq", "H2O"],
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
        obj=fs.ERD.efficiency_pump[0],
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
        obj=fs.ERD.control_volume.properties_out[0].pressure,
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
        obj=fs.costing.TPEC,
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
        obj=fs.costing.total_investment_factor,
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
        obj=fs.costing.maintenance_labor_chemical_factor,
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
        obj=fs.costing.capital_recovery_factor,
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

    # Feed
    exports.add(
        obj=fs.feed.properties[0].flow_vol_phase["Liq"],
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
        obj=fs.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"],
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
        obj=fs.product.properties[0].flow_vol,
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
        obj=fs.product.properties[0].conc_mass_phase_comp["Liq", "NaCl"],
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
        obj=fs.disposal.properties[0].flow_vol,
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
        obj=fs.disposal.properties[0].conc_mass_phase_comp["Liq", "NaCl"],
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
        obj=fs.RO.area,
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
        obj=fs.costing.specific_energy_consumption,
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
        obj=fs.costing.LCOW,
        name="Levelized cost of water",
        ui_units=fs.costing.base_currency / pyunits.m**3,
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

    # the UI sets `capital_recovery_factor`, so unfix `wacc`
    m.fs.costing.wacc.unfix()
    m.fs.costing.capital_recovery_factor.fix()

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


def solve_flowsheet(flowsheet=None):
    fs = flowsheet
    results = solve(fs)
    return results
