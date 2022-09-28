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
from watertap.ui.fsapi import FlowsheetInterface
from watertap.core.util.initialization import assert_degrees_of_freedom
from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.dye_desalination.dye_desalination_withRO import (
    build,
    set_operating_conditions,
    initialize_system,
    solve,
    add_costing,
)
from idaes.core.solvers import get_solver
from pyomo.environ import units as pyunits, assert_optimal_termination
from pyomo.util.check_units import assert_units_consistent


def export_to_ui():
    return FlowsheetInterface(
        name="Dye Desalination",
        do_export=export_variables,
        do_build=build_flowsheet,
        do_solve=solve_flowsheet,
    )


def export_variables(flowsheet=None, exports=None):
    fs = flowsheet
    # --- Input data ---
    # Feed conditions
    exports.add(
        obj=fs.feed.flow_vol[0],
        name="Volumetric flow rate",
        ui_units=pyunits.m**3 / pyunits.hr,
        display_units="m3/h",
        rounding=0,
        description="Inlet volumetric flow rate",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.conc_mass_comp[0, "dye"],
        name="Dye concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",  # can this be done by default?
        rounding=2,
        description="Inlet dye concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.conc_mass_comp[0, "tds"],
        name="TDS concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",  # can this be done by default?
        rounding=2,
        description="Inlet total dissolved solids (TDS) concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.tb_nf_ro.properties_out[0].temperature,
        name="Solution temperature",
        ui_units=pyunits.K,
        display_units="K",
        rounding=2,
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    # Unit model data, NF Pump
    exports.add(
        obj=fs.dye_separation.P1.eta_motor,
        name="NF pump- motor efficiency",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="NF pump- motor efficiency",
        is_input=True,
        input_category="rHGO Nanofiltration",
        is_output=False,
    )
    exports.add(
        obj=fs.dye_separation.P1.eta_pump,
        name="NF pump efficiency",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="NF pump efficiency",
        is_input=True,
        input_category="rHGO Nanofiltration",
        is_output=False,
    )
    # Unit model data, rHGO Nanofiltration
    exports.add(
        obj=fs.dye_separation.nanofiltration.recovery_frac_mass_H2O[0],
        name="Water recovery",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Water recovery [g-H2O treated/g-H2O inlet]",
        is_input=True,
        input_category="rHGO Nanofiltration",
        is_output=False,
    )
    exports.add(
        obj=fs.dye_separation.nanofiltration.removal_frac_mass_comp[0, "dye"],
        name="Mass removal fraction, dye",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Dye mass removal [kg-Dye removed/kg-Dye inlet]",
        is_input=True,
        input_category="rHGO Nanofiltration",
        is_output=False,
    )
    exports.add(
        obj=fs.dye_separation.nanofiltration.removal_frac_mass_comp[0, "tds"],
        name="Mass removal fraction, TDS",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="TDS mass removal [kg-Dye removed/kg-Dye inlet]",
        is_input=True,
        input_category="rHGO Nanofiltration",
        is_output=False,
    )
    exports.add(
        obj=fs.dye_separation.nanofiltration.water_permeability_coefficient[0],
        name="NF Water permeability Coefficient, A",
        ui_units=pyunits.L / pyunits.m**2 / pyunits.hour / pyunits.bar,
        display_units="LMH/bar",
        rounding=2,
        description="Membrane water permeability Coefficient, A",
        is_input=True,
        input_category="rHGO Nanofiltration",
        is_output=False,
    )
    v = fs.dye_separation.nanofiltration.applied_pressure
    exports.add(
        obj=v[0],
        name=v.doc,
        ui_units=getattr(pyunits, str(v._units)),
        display_units=str(v._units),
        rounding=1,
        description=v.doc,
        is_input=True,
        input_category="rHGO Nanofiltration",
        is_output=False,
    )

    # Unit model, secondary WWTP
    exports.add(
        obj=fs.pretreatment.wwtp.energy_electric_flow_vol_inlet,
        name="Specific energy consumption per inlet flow rate",
        ui_units=pyunits.kWh / pyunits.m**3,
        display_units="kWh/m3",
        rounding=2,
        description="Electrical energy consumption with respect to influent volumetric flow rate",
        is_input=True,
        input_category="Secondary Wastewater Treatment",
        is_output=False,
    )

    # Unit model, RO
    v = fs.desalination.RO.A_comp
    exports.add(
        obj=v[0, "H2O"],
        name="RO Water permeability coefficient, A",
        ui_units=getattr(pyunits, str(v._units)),
        display_units=str(v._units),
        rounding=12,
        description="Membrane water permeability Coefficient, A",
        is_input=True,
        input_category="Reverse Osmosis",
        is_output=False,
    )
    v = fs.desalination.RO.B_comp
    exports.add(
        obj=v[0, "TDS"],
        name="RO Salt permeability coefficient, B",
        ui_units=getattr(pyunits, str(v._units)),
        display_units=str(v._units),
        rounding=8,
        description="Membrane salt permeability Coefficient, B",
        is_input=True,
        input_category="Reverse Osmosis",
        is_output=False,
    )
    exports.add(
        obj=fs.desalination.P2.efficiency_pump[0],
        name="RO high-pressure pump efficiency",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="RO high-pressure pump efficiency",
        is_input=True,
        input_category="Reverse Osmosis",
        is_output=False,
    )
    exports.add(
        obj=fs.desalination.P3.efficiency_pump[0],
        name="RO booster pump efficiency",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="RO booster pump efficiency",
        is_input=True,
        input_category="Reverse Osmosis",
        is_output=False,
    )
    exports.add(
        obj=fs.desalination.PXR.efficiency_pressure_exchanger[0],
        name="Isobaric pressure exchanger efficiency",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Isobaric pressure exchanger efficiency",
        is_input=True,
        input_category="Reverse Osmosis",
        is_output=False,
    )
    # Unit cost data, NF pump
    exports.add(
        obj=fs.zo_costing.pump_electricity.pump_cost["default"],
        name="Pump cost",
        ui_units=fs.zo_costing.base_currency / (pyunits.m**3 / pyunits.hr),
        display_units="$/(m^3/hr)",
        rounding=0,
        description="Pump capital cost parameter",
        is_input=True,
        input_category="rHGO NF costing",
        is_output=False,
    )
    # Unit cost data, NF
    v = fs.zo_costing.nanofiltration.membrane_cost
    exports.add(
        obj=v["rHGO_dye_rejection"],
        name="NF membrane cost",
        ui_units=getattr(pyunits, str(v._units)),
        display_units=str(v._units),
        rounding=2,
        is_input=True,
        input_category="rHGO NF costing",
        is_output=False,
    )
    v = fs.zo_costing.nanofiltration.membrane_replacement_rate
    exports.add(
        obj=v["rHGO_dye_rejection"],
        name="NF membrane replacement rate",
        ui_units=getattr(pyunits, str(v._units)),
        display_units="fraction",
        rounding=2,
        is_input=True,
        input_category="rHGO NF costing",
        is_output=False,
    )
    # Unit cost data, RO
    v = fs.ro_costing.reverse_osmosis.membrane_cost
    exports.add(
        obj=v,
        name="RO membrane cost",
        ui_units=getattr(pyunits, str(v._units)),
        display_units=str(v._units),
        rounding=2,
        is_input=True,
        input_category="RO costing",
        is_output=False,
    )
    v = fs.ro_costing.reverse_osmosis.factor_membrane_replacement
    exports.add(
        obj=v,
        name="RO membrane replacement rate",
        ui_units=getattr(pyunits, str(v._units)),
        display_units="fraction",
        rounding=2,
        is_input=True,
        input_category="RO costing",
        is_output=False,
    )
    v = fs.ro_costing.high_pressure_pump.cost
    exports.add(
        obj=v,
        name="RO unit pump cost",
        ui_units=getattr(pyunits, str(v._units)),
        display_units=str(v._units),
        rounding=2,
        description="Unit pump cost for high-pressure and booster pump",
        is_input=True,
        input_category="RO costing",
        is_output=False,
    )
    v = fs.ro_costing.energy_recovery_device.pressure_exchanger_cost
    exports.add(
        obj=v,
        name=v.doc,
        ui_units=getattr(pyunits, str(v._units)),
        display_units=str(v._units),
        rounding=2,
        is_input=True,
        input_category="RO costing",
        is_output=False,
    )
    # System costs
    v = fs.ro_costing.electricity_base_cost
    exports.add(
        obj=v,
        name="Electricity cost",
        ui_units=getattr(pyunits, str(v._units)),
        display_units=str(v._units),
        rounding=2,
        is_input=True,
        input_category="System Costs",
        is_output=False,
    )
    v = fs.zo_costing.waste_disposal_cost
    exports.add(
        obj=v,
        name="Waste disposal cost per volume",
        ui_units=getattr(pyunits, str(v._units)),
        display_units=str(v._units),
        rounding=2,
        is_input=True,
        input_category="System Costs",
        is_output=False,
    )
    v = fs.zo_costing.dye_mass_cost
    exports.add(
        obj=v,
        name="Value of recovered dye",
        ui_units=getattr(pyunits, str(v._units)),
        display_units=str(v._units),
        rounding=2,
        is_input=True,
        input_category="System Costs",
        is_output=False,
    )
    v = fs.zo_costing.recovered_water_cost
    exports.add(
        obj=v,
        name="Value of recovered water",
        ui_units=getattr(pyunits, str(v._units)),
        display_units=str(v._units),
        rounding=2,
        is_input=True,
        input_category="System Costs",
        is_output=False,
    )
    # NF Results for rejection and membrane area
    exports.add(
        obj=fs.dye_separation.nanofiltration.rejection_comp[0, "dye"],
        name="Solute Rejection- dye",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Rejection rate, dye",
        is_input=False,
        output_category="rHGO Nanofiltration",
        is_output=True,
    )
    exports.add(
        obj=fs.dye_separation.nanofiltration.rejection_comp[0, "tds"],
        name="Solute Rejection- tds",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Rejection rate, tds",
        is_input=False,
        output_category="rHGO Nanofiltration",
        is_output=True,
    )
    v = fs.dye_separation.nanofiltration.area
    exports.add(
        obj=fs.dye_separation.nanofiltration.area,
        name=v.doc,
        ui_units=getattr(pyunits, str(v._units)),
        display_units=str(v._units),
        rounding=2,
        description="NF Membrane area",
        is_input=False,
        output_category="rHGO Nanofiltration",
        is_output=True,
    )
    # Outlets
    exports.add(
        obj=fs.dye_retentate.properties[0].flow_vol,
        name="Volumetric NF retentate flow rate",
        ui_units=pyunits.m**3 / pyunits.hr,
        display_units="m3/h",
        rounding=2,
        description="NF retentate flow rate",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.dye_retentate.properties[0].conc_mass_comp["dye"],
        name="NF retentate concentration, dye",
        ui_units=pyunits.kg / pyunits.m**3,
        display_units="kg/m3",
        rounding=2,
        description="NF retentate concentration, total dissolved solids (TDS)",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.dye_retentate.properties[0].conc_mass_comp["tds"],
        name="NF retentate concentration, tds",
        ui_units=pyunits.kg / pyunits.m**3,
        display_units="kg/m3",
        rounding=2,
        description="NF retentate concentration, total dissolved solids (TDS)",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.permeate.properties[0].flow_vol,
        name="Volumetric RO permeate flow rate",
        ui_units=pyunits.m**3 / pyunits.hr,
        display_units="m3/h",
        rounding=2,
        description="RO permeate flow rate",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.permeate.properties[0].conc_mass_phase_comp["Liq", "TDS"],
        name="RO permeate concentration, TDS",
        ui_units=pyunits.kg / pyunits.m**3,
        display_units="kg/m3",
        rounding=2,
        description="RO permeate concentration, total dissolved solids (TDS), including residual dye,",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.brine.properties[0].flow_vol,
        name="Volumetric RO brine flow rate",
        ui_units=pyunits.m**3 / pyunits.hr,
        display_units="m3/h",
        rounding=2,
        description="RO brine flow rate",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.brine.properties[0].conc_mass_phase_comp["Liq", "TDS"],
        name="RO brine concentration, TDS",
        ui_units=pyunits.kg / pyunits.m**3,
        display_units="kg/m3",
        rounding=2,
        description="RO brine concentration, total dissolved solids (TDS), including residual dye,",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )


def build_flowsheet():
    # build and solve initial flowsheet
    m = build()

    set_operating_conditions(m)
    assert_degrees_of_freedom(m, 0)
    assert_units_consistent(m)

    initialize_system(m)

    results = solve(m)
    assert_optimal_termination(results)

    add_costing(m)
    assert_degrees_of_freedom(m, 0)

    results = solve(m)
    assert_optimal_termination(results)
    return m.fs


def solve_flowsheet(flowsheet=None):
    fs = flowsheet
    results = solve(fs)
    return results
