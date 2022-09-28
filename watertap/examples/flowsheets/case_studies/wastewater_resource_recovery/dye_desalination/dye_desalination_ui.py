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
    # Unit model data, NF Pump
    exports.add(
        obj=fs.dye_separation.P1.eta_motor,
        name="NF pump- motor efficiency",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="NF pump- motor efficiency",
        is_input=True,
        input_category="NF Pump",
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
        input_category="NF Pump",
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
        name="Water permeability Coefficient, A",
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
    # Unit cost data, NF pump
    exports.add(
        obj=fs.zo_costing.pump_electricity.pump_cost["default"],
        name="Pump cost",
        ui_units=fs.zo_costing.base_currency / (pyunits.m**3 / pyunits.hr),
        display_units="$/(m^3/hr)",
        rounding=0,
        description="Pump capital cost parameter",
        is_input=True,
        input_category="Pump costing",
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
        is_input=True,
        input_category="rHGO NF costing",
        is_output=False,
    )
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
