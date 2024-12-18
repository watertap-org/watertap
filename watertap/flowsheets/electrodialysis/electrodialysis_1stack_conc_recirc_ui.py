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
from idaes_flowsheet_processor.api import FlowsheetInterface
from watertap.flowsheets.electrodialysis.electrodialysis_1stack_conc_recirc import (
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
        ui_units=pyunits.m**3 / pyunits.s,
        display_units="m^3 s^-1",
        rounding=3,
        description="Inlet water mass flowrate",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.properties[0].conc_mol_phase_comp["Liq", "Na_+"],
        name="Na_+ molar concentration",
        ui_units=pyunits.mol / pyunits.m**3,
        display_units="mol m^-3",
        rounding=3,
        description="Molar concentration of Na_+ ions in the feed solution",
        is_input=True,
        input_category="Feed",
        is_output=False,
    )

    exports.add(
        obj=fs.feed.properties[0].conc_mol_phase_comp["Liq", "Cl_-"],
        name="Cl_- molar concentration",
        ui_units=pyunits.mol / pyunits.m**3,
        display_units="mol m^-3",
        rounding=3,
        description="Molar concentration of Cl_- ions in the feed solution",
        is_input=True,
        input_category="Feed",
        is_output=False,
    )

    # ED stack conditions
    exports.add(
        obj=fs.EDstack.voltage_applied[0],
        name="Stack voltage",
        ui_units=pyunits.volt,
        display_units="V",
        rounding=3,
        description="Applied constant voltage on the ED stack",
        is_input=True,
        input_category="ED stack",
        is_output=True,
        output_category="ED stack",
    )
    exports.add(
        obj=fs.EDstack.cell_pair_num,
        name="Cell pair number",
        ui_units=pyunits.dimensionless,
        display_units="",
        rounding=0,
        description="Cell pair number in a single ED stack",
        is_input=True,
        input_category="ED stack",
        is_output=True,
        output_category="ED stack",
    )
    exports.add(
        obj=fs.EDstack.channel_height,
        name="Channel height",
        ui_units=pyunits.meter,
        display_units="m",
        rounding=4,
        description="Channel height of the ED flow channels",
        is_input=True,
        input_category="ED stack",
        is_output=True,
        output_category="ED stack",
    )
    exports.add(
        obj=fs.EDstack.cell_width,
        name="Cell width",
        ui_units=pyunits.meter,
        display_units="m",
        rounding=3,
        description="The width of ED cell or stack",
        is_input=True,
        input_category="ED stack",
        is_output=True,
        output_category="ED stack",
    )
    exports.add(
        obj=fs.EDstack.cell_length,
        name="Channel length",
        ui_units=pyunits.meter,
        display_units="m",
        rounding=3,
        description="The length of ED cell or stack",
        is_input=True,
        input_category="ED stack",
        is_output=True,
        output_category="ED stack",
    )

    exports.add(
        obj=fs.EDstack.spacer_porosity,
        name="Spacer porosity",
        ui_units=pyunits.dimensionless,
        display_units="",
        rounding=2,
        description="Porosity of the flow spacer",
        is_input=True,
        input_category="ED stack",
        is_output=False,
    )

    exports.add(
        obj=fs.EDstack.spacer_specific_area,
        name="Spacer specific surface area",
        ui_units=pyunits.meter**-1,
        display_units="m^-1",
        rounding=0,
        description="Specific surface area of the flow spacer",
        is_input=True,
        input_category="ED stack",
        is_output=False,
    )

    exports.add(
        obj=fs.EDstack.electrodes_resistance,
        name="Electrode resistance",
        ui_units=pyunits.ohm * pyunits.meter**2,
        display_units="ohm m^2",
        rounding=2,
        description="Areal resistance of the two electrodes",
        is_input=True,
        input_category="ED stack",
        is_output=False,
    )

    exports.add(
        obj=fs.EDstack.current_utilization,
        name="Current utilization coefficient",
        ui_units=pyunits.dimensionless,
        display_units="",
        rounding=2,
        description="Current utilization coefficient",
        is_input=True,
        input_category="ED stack",
        is_output=False,
    )

    exports.add(
        obj=fs.recovery_vol_H2O,
        name="Water recovery by volume",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Product water recovery by volume",
        is_input=True,
        input_category="ED stack",
        is_output=True,
        output_category="ED stack",
    )

    # Membrane-related properties
    exports.add(
        obj=fs.EDstack.water_trans_number_membrane["cem"],
        name="Water electroosmosis transport number of CEM",
        ui_units=pyunits.dimensionless,
        display_units="",
        rounding=2,
        description="Water electroosmosis transport number of the cation exchange membrane",
        is_input=True,
        input_category="Membrane properties",
        is_output=False,
    )
    exports.add(
        obj=fs.EDstack.water_trans_number_membrane["aem"],
        name="Water electroosmosis transport number of AEM",
        ui_units=pyunits.dimensionless,
        display_units="",
        rounding=2,
        description="Water electroosmosis transport number of the anion exchange membrane",
        is_input=True,
        input_category="Membrane properties",
        is_output=False,
    )

    exports.add(
        obj=fs.EDstack.water_permeability_membrane["cem"],
        name="Water osmosis permeability of CEM",
        ui_units=pyunits.meter * pyunits.second**-1 * pyunits.pascal**-1,
        display_units="m s^-1 Pa^-1",
        rounding=2,
        description="Water osmosis permeability of the cation exchange membrane",
        is_input=True,
        input_category="Membrane properties",
        is_output=False,
    )
    exports.add(
        obj=fs.EDstack.water_permeability_membrane["aem"],
        name="Water osmosis permeability of AEM",
        ui_units=pyunits.meter * pyunits.second**-1 * pyunits.pascal**-1,
        display_units="m s^-1 Pa^-1",
        rounding=2,
        description="Water osmosis permeability of the anion exchange membrane",
        is_input=True,
        input_category="Membrane properties",
        is_output=False,
    )

    exports.add(
        obj=fs.EDstack.membrane_areal_resistance["cem"],
        name="Areal resistnace of CEM",
        ui_units=pyunits.ohm * pyunits.meter**2,
        display_units="ohm m^2",
        rounding=2,
        description="Constant areal resistance of the cation exchange membrane measured in concentrated electrolyte.",
        is_input=True,
        input_category="Membrane properties",
        is_output=False,
    )
    exports.add(
        obj=fs.EDstack.membrane_areal_resistance["aem"],
        name="Areal resistnace of AEM",
        ui_units=pyunits.ohm * pyunits.meter**2,
        display_units="ohm m^2",
        rounding=2,
        description="Constant areal resistance of the anion exchange membrane measured in concentrated electrolyte.",
        is_input=True,
        input_category="Membrane properties",
        is_output=False,
    )

    exports.add(
        obj=fs.EDstack.membrane_thickness["cem"],
        name="Thickness of CEM",
        ui_units=pyunits.meter,
        display_units="m",
        rounding=2,
        description="Thickness of the cation exchange membrane",
        is_input=True,
        input_category="Membrane properties",
        is_output=False,
    )
    exports.add(
        obj=fs.EDstack.membrane_thickness["aem"],
        name="Thickness of AEM",
        ui_units=pyunits.meter,
        display_units="m",
        rounding=2,
        description="Thickness of the anion exchange membrane",
        is_input=True,
        input_category="Membrane properties",
        is_output=False,
    )

    exports.add(
        obj=fs.EDstack.solute_diffusivity_membrane["cem", "Na_+"],
        name="Na_+ diffusivity as solute in CEM",
        ui_units=pyunits.meter**2 * pyunits.second**-1,
        display_units="m^2 s^-1",
        rounding=2,
        description="Na_+ diffusivity as solute in the cation exchange membrane (the mass diffusivity of the corresponding neutral solute should be used.)",
        is_input=True,
        input_category="Membrane properties",
        is_output=False,
    )

    exports.add(
        obj=fs.EDstack.solute_diffusivity_membrane["aem", "Na_+"],
        name="Na_+ diffusivity as solute in AEM",
        ui_units=pyunits.meter**2 * pyunits.second**-1,
        display_units="m^2 s^-1",
        rounding=2,
        description="Na_+ diffusivity as solute in the anion exchange membrane (the mass diffusivity of the corresponding neutral solute should be used.)",
        is_input=True,
        input_category="Membrane properties",
        is_output=False,
    )

    exports.add(
        obj=fs.EDstack.solute_diffusivity_membrane["cem", "Cl_-"],
        name="Cl_- diffusivity as solute in CEM",
        ui_units=pyunits.meter**2 * pyunits.second**-1,
        display_units="m^2 s^-1",
        rounding=2,
        description="Cl_- diffusivity as solute in the cation exchange membrane (the mass diffusivity of the corresponding neutral solute should be used.)",
        is_input=True,
        input_category="Membrane properties",
        is_output=False,
    )

    exports.add(
        obj=fs.EDstack.solute_diffusivity_membrane["aem", "Cl_-"],
        name="Cl_- diffusivity as solute in AEM",
        ui_units=pyunits.meter**2 * pyunits.second**-1,
        display_units="m^2 s^-1",
        rounding=2,
        description="Cl_- diffusivity as solute in the anion exchange membrane (the mass diffusivity of the corresponding neutral solute should be used.)",
        is_input=True,
        input_category="Membrane properties",
        is_output=False,
    )

    exports.add(
        obj=fs.EDstack.ion_trans_number_membrane["cem", "Na_+"],
        name="Na_+ transport number in CEM",
        ui_units=pyunits.dimensionless,
        display_units="",
        rounding=1,
        description="Na_+ transport number in the cation exchange membrane",
        is_input=True,
        input_category="Membrane properties",
        is_output=False,
    )
    exports.add(
        obj=fs.EDstack.ion_trans_number_membrane["aem", "Na_+"],
        name="Na_+ transport number in AEM",
        ui_units=pyunits.dimensionless,
        display_units="",
        rounding=1,
        description="Na_+ transport number in the anion exchange membrane",
        is_input=True,
        input_category="Membrane properties",
        is_output=False,
    )
    exports.add(
        obj=fs.EDstack.ion_trans_number_membrane["cem", "Cl_-"],
        name="Cl_- transport number in CEM",
        ui_units=pyunits.dimensionless,
        display_units="",
        rounding=1,
        description="Cl_- transport number in the cation exchange membrane",
        is_input=True,
        input_category="Membrane properties",
        is_output=False,
    )
    exports.add(
        obj=fs.EDstack.ion_trans_number_membrane["aem", "Cl_-"],
        name="Cl_- transport number in AEM",
        ui_units=pyunits.dimensionless,
        display_units="",
        rounding=1,
        description="Cl_- transport number in the anion exchange membrane",
        is_input=True,
        input_category="Membrane properties",
        is_output=False,
    )

    # Feed pump properties
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
        obj=fs.feed_salinity,
        name="NaCl mass concentration",
        ui_units=pyunits.kg / pyunits.m**3,
        display_units="kg m^-3",
        rounding=2,
        description="Feed NaCl mass concentration",
        is_input=False,
        is_output=True,
        output_category="Feed",
    )

    # Product
    exports.add(
        obj=fs.prod.properties[0].flow_vol_phase["Liq"],
        name="Volumetric flow rate",
        ui_units=pyunits.m**3 / pyunits.hr,
        display_units="m3/h",
        rounding=2,
        description="Product water volumetric flow rate",
        is_input=False,
        is_output=True,
        output_category="Product",
    )
    exports.add(
        obj=fs.prod.properties[0].conc_mol_phase_comp["Liq", "Na_+"],
        name="NaCl molar concentration",
        ui_units=pyunits.mol / pyunits.meter**3,
        display_units="mol m^-3",
        rounding=3,
        description="Product water NaCl molar concentration",
        is_input=False,
        is_output=True,
        output_category="Product",
    )
    exports.add(
        obj=fs.product_salinity,
        name="NaCl mass concentration",
        ui_units=pyunits.kg / pyunits.meter**3,
        display_units="kg m^-3",
        rounding=3,
        description="Product water NaCl mass concentration",
        is_input=False,
        is_output=True,
        output_category="Product",
    )

    # Disposal
    exports.add(
        obj=fs.disp.properties[0].flow_vol_phase["Liq"],
        name="Volumetric flow rate",
        ui_units=pyunits.m**3 / pyunits.hr,
        display_units="m3/h",
        rounding=2,
        description="Disposal water volumetric flow rate",
        is_input=False,
        is_output=True,
        output_category="Disposal",
    )
    exports.add(
        obj=fs.disp.properties[0].conc_mol_phase_comp["Liq", "Na_+"],
        name="NaCl molar concentration",
        ui_units=pyunits.mol / pyunits.meter**3,
        display_units="mol m^-3",
        rounding=3,
        description="Disposal water NaCl molar concentration",
        is_input=False,
        is_output=True,
        output_category="Disposal",
    )
    exports.add(
        obj=fs.disposal_salinity,
        name="NaCl mass concentration",
        ui_units=pyunits.kg / pyunits.meter**3,
        display_units="kg m^-3",
        rounding=3,
        description="Disposal water NaCl mass concentration",
        is_input=False,
        is_output=True,
        output_category="Product",
    )

    # System metrics

    exports.add(
        obj=fs.mem_area,
        name="Membrane area",
        ui_units=pyunits.m**2,
        display_units="m^2",
        rounding=2,
        description="Membrane area of the cation or anion exchange membranes",
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


def build_flowsheet(build_options=None, **kwargs):

    # the UI sets `capital_recovery_factor`, so unfix `wacc`
    m = build()
    m.fs.costing.wacc.unfix()
    m.fs.costing.capital_recovery_factor.fix()

    solver = get_solver()

    # build, set, and initialize
    m = build()
    init_arg = {
        ("flow_vol_phase", ("Liq")): 5.2e-4,
        ("conc_mol_phase_comp", ("Liq", "Na_+")): 34.188,
        ("conc_mol_phase_comp", ("Liq", "Cl_-")): 34.188,
    }  # Corresponding to C_feed = 2g/L
    m.fs.feed.properties.calculate_state(
        init_arg,
        hold_state=True,
    )
    m.fs.EDstack.voltage_applied[0].fix(10)
    m.fs.recovery_vol_H2O.fix(0.7)
    _condition_base(m)
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
