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
from watertap.ui.fsapi import FlowsheetInterface
from watertap.core.util.initialization import assert_degrees_of_freedom
from watertap.examples.flowsheets.case_studies.full_water_resource_recovery_facility.BSM2 import (
    build,
    set_operating_conditions,
    initialize_system,
    solve,
    add_costing,
)
from pyomo.environ import units as pyunits, assert_optimal_termination
from pyomo.util.check_units import assert_units_consistent


def export_to_ui():
    return FlowsheetInterface(
        name="BSM2",
        do_export=export_variables,
        do_build=build_flowsheet,
        do_solve=solve_flowsheet,
    )


def export_variables(flowsheet=None, exports=None):
    fs = flowsheet
    # --- Input data ---
    # Feed conditions
    exports.add(
        obj=fs.FeedWater.flow_vol[0],
        name="Volumetric flow rate",
        ui_units=pyunits.m**3 / pyunits.day,
        display_units="m3/day",
        rounding=2,
        description="Inlet volumetric flow rate",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.FeedWater.conc_mass_comp[0, "S_I"],
        name="S_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=2,
        description="Inlet soluble inert organic matter concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.FeedWater.conc_mass_comp[0, "S_S"],
        name="S_S concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=2,
        description="Inlet readily biodegradable substrate concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.FeedWater.conc_mass_comp[0, "X_I"],
        name="X_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=2,
        description="Inlet particulate inert organic matter concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.FeedWater.conc_mass_comp[0, "X_S"],
        name="X_S concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=2,
        description="Inlet slowly biodegradable substrate concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.FeedWater.conc_mass_comp[0, "X_BH"],
        name="X_BH concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=2,
        description="Inlet active heterotrophic biomass concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.FeedWater.conc_mass_comp[0, "X_BA"],
        name="X_BA concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=2,
        description="Inlet active autotrophic biomass concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.FeedWater.conc_mass_comp[0, "X_P"],
        name="X_P concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=2,
        description="Inlet particulate products arising from biomass decay concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.FeedWater.conc_mass_comp[0, "S_O"],
        name="S_O concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=2,
        description="Inlet oxygen concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.FeedWater.conc_mass_comp[0, "S_NO"],
        name="S_NO concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=2,
        description="Inlet nitrate and nitrite nitrogen concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.FeedWater.conc_mass_comp[0, "S_NH"],
        name="S_NH concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=2,
        description="Inlet ammonium and ammonia nitrogen concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.FeedWater.conc_mass_comp[0, "S_ND"],
        name="S_ND concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=2,
        description="Inlet biodegradable organic nitrogen concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.FeedWater.conc_mass_comp[0, "X_ND"],
        name="X_ND concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=2,
        description="Inlet particulate biodegradable organic nitrogen concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.FeedWater.properties[0].alkalinity,
        name="S_ALK",
        ui_units=pyunits.mol / pyunits.m**3,
        display_units="mol/m3",
        rounding=2,
        description="Alkalinity",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )

    # Unit model data, secondary clarifier
    exports.add(
        obj=fs.CL1.split_fraction[0, "effluent", "H2O"],
        name="H2O split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Water split fraction",
        is_input=True,
        input_category="Secondary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL1.split_fraction[0, "effluent", "S_I"],
        name="S_I split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Soluble inert organic matter split fraction",
        is_input=True,
        input_category="Secondary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL1.split_fraction[0, "effluent", "S_S"],
        name="S_S split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Readily biodegradable substrate split fraction",
        is_input=True,
        input_category="Secondary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL1.split_fraction[0, "effluent", "X_I"],
        name="X_I split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Particulate inert organic matter split fraction",
        is_input=True,
        input_category="Secondary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL1.split_fraction[0, "effluent", "X_S"],
        name="X_S split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Slowly biodegradable substrate split fraction",
        is_input=True,
        input_category="Secondary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL1.split_fraction[0, "effluent", "X_BH"],
        name="X_BH split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Active heterotrophic biomass split fraction",
        is_input=True,
        input_category="Secondary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL1.split_fraction[0, "effluent", "X_BA"],
        name="X_BA split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Active autotrophic biomass split fraction",
        is_input=True,
        input_category="Secondary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL1.split_fraction[0, "effluent", "X_P"],
        name="X_P split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Particulate products arising from biomass decay split fraction",
        is_input=True,
        input_category="Secondary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL1.split_fraction[0, "effluent", "S_O"],
        name="S_O split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Oxygen split fraction",
        is_input=True,
        input_category="Secondary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL1.split_fraction[0, "effluent", "S_NO"],
        name="S_NO split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Nitrate and nitrite nitrogen split fraction",
        is_input=True,
        input_category="Secondary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL1.split_fraction[0, "effluent", "S_NH"],
        name="S_NH split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Ammonium and ammonia nitrogen split fraction",
        is_input=True,
        input_category="Secondary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL1.split_fraction[0, "effluent", "S_ND"],
        name="S_ND split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Soluble biodegradable organic nitrogen split fraction",
        is_input=True,
        input_category="Secondary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL1.split_fraction[0, "effluent", "X_ND"],
        name="X_ND split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Particulate biodegradable organic nitrogen split fraction",
        is_input=True,
        input_category="Secondary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL1.split_fraction[0, "effluent", "S_ALK"],
        name="S_ALK split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Alkalinity split fraction",
        is_input=True,
        input_category="Secondary clarifier",
        is_output=False,
    )
    # Unit model data, primary clarifier
    exports.add(
        obj=fs.CL.split_fraction[0, "effluent", "H2O"],
        name="H2O split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Water split fraction",
        is_input=True,
        input_category="Primary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL.split_fraction[0, "effluent", "S_I"],
        name="S_I split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Soluble inert organic matter split fraction",
        is_input=True,
        input_category="Primary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL.split_fraction[0, "effluent", "S_S"],
        name="S_S split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Readily biodegradable substrate split fraction",
        is_input=True,
        input_category="Primary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL.split_fraction[0, "effluent", "X_I"],
        name="X_I split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Particulate inert organic matter split fraction",
        is_input=True,
        input_category="Primary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL.split_fraction[0, "effluent", "X_S"],
        name="X_S split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Slowly biodegradable substrate split fraction",
        is_input=True,
        input_category="Primary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL.split_fraction[0, "effluent", "X_BH"],
        name="X_BH split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Active heterotrophic biomass split fraction",
        is_input=True,
        input_category="Primary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL.split_fraction[0, "effluent", "X_BA"],
        name="X_BA split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Active autotrophic biomass split fraction",
        is_input=True,
        input_category="Primary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL.split_fraction[0, "effluent", "X_P"],
        name="X_P split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Particulate products arising from biomass decay split fraction",
        is_input=True,
        input_category="Primary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL.split_fraction[0, "effluent", "S_O"],
        name="S_O split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Oxygen split fraction",
        is_input=True,
        input_category="Primary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL.split_fraction[0, "effluent", "S_NO"],
        name="S_NO split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Nitrate and nitrite nitrogen split fraction",
        is_input=True,
        input_category="Primary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL.split_fraction[0, "effluent", "S_NH"],
        name="S_NH split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Ammonium and ammonia nitrogen split fraction",
        is_input=True,
        input_category="Primary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL.split_fraction[0, "effluent", "S_ND"],
        name="S_ND split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Soluble biodegradable organic nitrogen split fraction",
        is_input=True,
        input_category="Primary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL.split_fraction[0, "effluent", "X_ND"],
        name="X_ND split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Particulate biodegradable organic nitrogen split fraction",
        is_input=True,
        input_category="Primary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL.split_fraction[0, "effluent", "S_ALK"],
        name="S_ALK split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Alkalinity split fraction",
        is_input=True,
        input_category="Primary clarifier",
        is_output=False,
    )
    # Unit model data, anaerobic digester
    exports.add(
        obj=fs.RADM.volume_liquid[0],
        name="Liquid volume",
        ui_units=pyunits.m**3,
        display_units="m3",
        rounding=1,
        description="Liquid volume",
        is_input=True,
        input_category="Anaerobic digester",
        is_output=False,
    )
    exports.add(
        obj=fs.RADM.volume_vapor[0],
        name="Vapor volume",
        ui_units=pyunits.m**3,
        display_units="m3",
        rounding=1,
        description="Vapor volume",
        is_input=True,
        input_category="Anaerobic digester",
        is_output=False,
    )

    # System costing
    # exports.add(
    #     obj=fs.costing.utilization_factor,
    #     name="Utilization factor",
    #     ui_units=pyunits.dimensionless,
    #     display_units="fraction",
    #     rounding=2,
    #     description="Utilization factor - [annual use hours/total hours in year]",
    #     is_input=True,
    #     input_category="System costing",
    #     is_output=False,
    # )
    # exports.add(
    #     obj=fs.costing.TIC,
    #     name="Practical investment factor",
    #     ui_units=pyunits.dimensionless,
    #     display_units="fraction",
    #     rounding=1,
    #     description="Practical investment factor - [total investment cost/direct "
    #     "capital costs]",
    #     is_input=True,
    #     input_category="System costing",
    #     is_output=False,
    # )
    # exports.add(
    #     obj=fs.costing.plant_lifetime,
    #     name="Plant lifetime",
    #     ui_units=pyunits.year,
    #     display_units="years",
    #     rounding=1,
    #     description="Plant lifetime",
    #     is_input=True,
    #     input_category="System costing",
    #     is_output=False,
    # )
    # exports.add(
    #     obj=fs.costing.wacc,
    #     name="Discount rate",
    #     ui_units=pyunits.dimensionless,
    #     display_units="fraction",
    #     rounding=2,
    #     description="Discount rate used in calculating the capital annualization",
    #     is_input=True,
    #     input_category="System costing",
    #     is_output=False,
    # )
    # exports.add(
    #     obj=fs.costing.maintenance_costs_percent_FCI,
    #     name="Fixed operating cost factor",
    #     ui_units=1 / pyunits.year,
    #     display_units="fraction/year",
    #     rounding=2,
    #     description="Fixed operating cost factor - [annual fixed operating cost/total "
    #     "investment cost]",
    #     is_input=True,
    #     input_category="System costing",
    #     is_output=False,
    # )
    # exports.add(
    #     obj=fs.costing.electricity_cost,
    #     name="Electricity cost",
    #     ui_units=fs.costing.base_currency / pyunits.kWh,
    #     display_units="$/kWh",
    #     rounding=3,
    #     description="Electricity cost",
    #     is_input=True,
    #     input_category="System costing",
    #     is_output=False,
    # )

    # Outlets
    exports.add(
        obj=fs.Treated.properties[0].flow_vol,
        name="Flow rate",
        ui_units=pyunits.m**3 / pyunits.day,
        display_units="m3/day",
        rounding=2,
        description="Outlet treated stream flow rate",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Effluent",
    )
    exports.add(
        obj=fs.Treated.properties[0].conc_mass_comp["S_I"],
        name="S_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet soluble inert organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Effluent",
    )
    exports.add(
        obj=fs.Treated.properties[0].conc_mass_comp["S_S"],
        name="S_S concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet readily biodegradable substrate concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Effluent",
    )
    exports.add(
        obj=fs.Treated.properties[0].conc_mass_comp["X_I"],
        name="X_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet particulate inert organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Effluent",
    )
    exports.add(
        obj=fs.Treated.properties[0].conc_mass_comp["X_S"],
        name="X_S concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet slowly biodegradable substrate concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Effluent",
    )
    exports.add(
        obj=fs.Treated.properties[0].conc_mass_comp["X_BH"],
        name="X_BH concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet active heterotrophic biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Effluent",
    )
    exports.add(
        obj=fs.Treated.properties[0].conc_mass_comp["X_BA"],
        name="X_BA concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet active autotrophic biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Effluent",
    )
    exports.add(
        obj=fs.Treated.properties[0].conc_mass_comp["X_P"],
        name="X_P concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet particulate products arising from biomass decay concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Effluent",
    )
    exports.add(
        obj=fs.Treated.properties[0].conc_mass_comp["S_O"],
        name="S_O concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet oxygen concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Effluent",
    )
    exports.add(
        obj=fs.Treated.properties[0].conc_mass_comp["S_NO"],
        name="S_NO concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet nitrate and nitrite nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Effluent",
    )
    exports.add(
        obj=fs.Treated.properties[0].conc_mass_comp["S_NH"],
        name="S_NH concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet ammonium and ammonia nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Effluent",
    )
    exports.add(
        obj=fs.Treated.properties[0].conc_mass_comp["S_ND"],
        name="S_ND concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet soluble biodegradable organic nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Effluent",
    )
    exports.add(
        obj=fs.Treated.properties[0].conc_mass_comp["X_ND"],
        name="X_ND concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet particulate biodegradable organic nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Effluent",
    )
    exports.add(
        obj=fs.Treated.properties[0].alkalinity,
        name="S_ALK concentration",
        ui_units=pyunits.mol / pyunits.m**3,
        display_units="mol/m3",
        rounding=5,
        description="Outlet alkalinity concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Effluent",
    )

    exports.add(
        obj=fs.DU.underflow.flow_vol[0],
        name="Flow rate",
        ui_units=pyunits.m**3 / pyunits.day,
        display_units="m3/day",
        rounding=2,
        description="Outlet sludge disposal flow rate",
        is_input=False,
        is_output=True,
        output_category="Dewatered Sludge",
    )
    exports.add(
        obj=fs.DU.underflow.conc_mass_comp[0, "S_I"],
        name="S_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet soluble inert organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatered Sludge",
    )
    exports.add(
        obj=fs.DU.underflow.conc_mass_comp[0, "S_S"],
        name="S_S concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet readily biodegradable substrate concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatered Sludge",
    )
    exports.add(
        obj=fs.DU.underflow.conc_mass_comp[0, "X_I"],
        name="X_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet particulate inert organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatered Sludge",
    )
    exports.add(
        obj=fs.DU.underflow.conc_mass_comp[0, "X_S"],
        name="X_S concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet slowly biodegradable substrate concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatered Sludge",
    )
    exports.add(
        obj=fs.DU.underflow.conc_mass_comp[0, "X_BH"],
        name="X_BH concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet active heterotrophic biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatered Sludge",
    )
    exports.add(
        obj=fs.DU.underflow.conc_mass_comp[0, "X_BA"],
        name="X_BA concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet active autotrophic biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatered Sludge",
    )
    exports.add(
        obj=fs.DU.underflow.conc_mass_comp[0, "X_P"],
        name="X_P concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet particulate products arising from biomass decay concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatered Sludge",
    )
    exports.add(
        obj=fs.DU.underflow.conc_mass_comp[0, "S_O"],
        name="S_O concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet oxygen concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatered Sludge",
    )
    exports.add(
        obj=fs.DU.underflow.conc_mass_comp[0, "S_NO"],
        name="S_NO concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet nitrate and nitrite nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatered Sludge",
    )
    exports.add(
        obj=fs.DU.underflow.conc_mass_comp[0, "S_NH"],
        name="S_NH concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet ammonium and ammonia nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatered Sludge",
    )
    exports.add(
        obj=fs.DU.underflow.conc_mass_comp[0, "S_ND"],
        name="Sludge S_ND concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet soluble biodegradable organic nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatered Sludge",
    )
    exports.add(
        obj=fs.DU.underflow.conc_mass_comp[0, "X_ND"],
        name="X_ND concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet particulate biodegradable organic nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatered Sludge",
    )
    exports.add(
        obj=fs.DU.underflow.alkalinity[0],
        name="S_ALK concentration",
        ui_units=pyunits.mol / pyunits.m**3,
        display_units="mol/m3",
        rounding=5,
        description="Outlet alkalinity concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatered Sludge",
    )

    exports.add(
        obj=fs.DU.overflow.flow_vol[0],
        name="Flow rate",
        ui_units=pyunits.m**3 / pyunits.day,
        display_units="m3/day",
        rounding=2,
        description="Dewatering unit liquid outlet flow rate",
        is_input=False,
        is_output=True,
        output_category="Dewatering Unit Liquid Outlet",
    )
    exports.add(
        obj=fs.DU.overflow.conc_mass_comp[0, "S_I"],
        name="S_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet soluble inert organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatering Unit Liquid Outlet",
    )
    exports.add(
        obj=fs.DU.overflow.conc_mass_comp[0, "S_S"],
        name="S_S concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet readily biodegradable substrate concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatering Unit Liquid Outlet",
    )
    exports.add(
        obj=fs.DU.overflow.conc_mass_comp[0, "X_I"],
        name="X_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet particulate inert organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatering Unit Liquid Outlet",
    )
    exports.add(
        obj=fs.DU.overflow.conc_mass_comp[0, "X_S"],
        name="X_S concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet slowly biodegradable substrate concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatering Unit Liquid Outlet",
    )
    exports.add(
        obj=fs.DU.overflow.conc_mass_comp[0, "X_BH"],
        name="X_BH concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet active heterotrophic biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatering Unit Liquid Outlet",
    )
    exports.add(
        obj=fs.DU.overflow.conc_mass_comp[0, "X_BA"],
        name="X_BA concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet active autotrophic biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatering Unit Liquid Outlet",
    )
    exports.add(
        obj=fs.DU.overflow.conc_mass_comp[0, "X_P"],
        name="X_P concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet particulate products arising from biomass decay concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatering Unit Liquid Outlet",
    )
    exports.add(
        obj=fs.DU.overflow.conc_mass_comp[0, "S_O"],
        name="S_O concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet oxygen concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatering Unit Liquid Outlet",
    )
    exports.add(
        obj=fs.DU.overflow.conc_mass_comp[0, "S_NO"],
        name="S_NO concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet nitrate and nitrite nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatering Unit Liquid Outlet",
    )
    exports.add(
        obj=fs.DU.overflow.conc_mass_comp[0, "S_NH"],
        name="S_NH concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet ammonium and ammonia nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatering Unit Liquid Outlet",
    )
    exports.add(
        obj=fs.DU.overflow.conc_mass_comp[0, "S_ND"],
        name="S_ND concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet soluble biodegradable organic nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatering Unit Liquid Outlet",
    )
    exports.add(
        obj=fs.DU.overflow.conc_mass_comp[0, "X_ND"],
        name="X_ND concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet particulate biodegradable organic nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatering Unit Liquid Outlet",
    )
    exports.add(
        obj=fs.DU.overflow.alkalinity[0],
        name="S_ALK concentration",
        ui_units=pyunits.mol / pyunits.m**3,
        display_units="mol/m3",
        rounding=5,
        description="Outlet alkalinity concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatering Unit Liquid Outlet",
    )

    exports.add(
        obj=fs.asm_adm.inlet.flow_vol[0],
        name="Flow rate",
        ui_units=pyunits.m**3 / pyunits.day,
        display_units="m3/day",
        rounding=2,
        description="ASM1/ADM1 translator inlet flow rate",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Inlet",
    )
    exports.add(
        obj=fs.asm_adm.inlet.conc_mass_comp[0, "S_I"],
        name="S_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet soluble inert organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Inlet",
    )
    exports.add(
        obj=fs.asm_adm.inlet.conc_mass_comp[0, "S_S"],
        name="S_S concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet readily biodegradable substrate concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Inlet",
    )
    exports.add(
        obj=fs.asm_adm.inlet.conc_mass_comp[0, "X_I"],
        name="X_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet particulate inert organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Inlet",
    )
    exports.add(
        obj=fs.asm_adm.inlet.conc_mass_comp[0, "X_S"],
        name="X_S concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet slowly biodegradable substrate concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Inlet",
    )
    exports.add(
        obj=fs.asm_adm.inlet.conc_mass_comp[0, "X_BH"],
        name="X_BH concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet active heterotrophic biomass concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Inlet",
    )
    exports.add(
        obj=fs.asm_adm.inlet.conc_mass_comp[0, "X_BA"],
        name="X_BA concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet active autotrophic biomass concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Inlet",
    )
    exports.add(
        obj=fs.asm_adm.inlet.conc_mass_comp[0, "X_P"],
        name="X_P concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet particulate products arising from biomass decay concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Inlet",
    )
    exports.add(
        obj=fs.asm_adm.inlet.conc_mass_comp[0, "S_O"],
        name="S_O concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet oxygen concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Inlet",
    )
    exports.add(
        obj=fs.asm_adm.inlet.conc_mass_comp[0, "S_NO"],
        name="S_NO concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet nitrate and nitrite nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Inlet",
    )
    exports.add(
        obj=fs.asm_adm.inlet.conc_mass_comp[0, "S_NH"],
        name="S_NH concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet ammonium and ammonia nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Inlet",
    )
    exports.add(
        obj=fs.asm_adm.inlet.conc_mass_comp[0, "S_ND"],
        name="S_ND concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet soluble biodegradable organic nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Inlet",
    )
    exports.add(
        obj=fs.asm_adm.inlet.conc_mass_comp[0, "X_ND"],
        name="X_ND concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet particulate biodegradable organic nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Inlet",
    )
    exports.add(
        obj=fs.asm_adm.inlet.alkalinity[0],
        name="S_ALK concentration",
        ui_units=pyunits.mol / pyunits.m**3,
        display_units="mol/m3",
        rounding=5,
        description="Inlet alkalinity concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Inlet",
    )

    exports.add(
        obj=fs.asm_adm.outlet.flow_vol[0],
        name="Flow rate",
        ui_units=pyunits.m**3 / pyunits.day,
        display_units="m3/day",
        rounding=2,
        description="ASM1/ADM1 translator outlet flow rate",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Outlet",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "S_su"],
        name="S_su concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Monosaccharides concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Outlet",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "S_aa"],
        name="S_aa concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Amino acids concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Outlet",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "S_fa"],
        name="S_fa concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Long chain fatty acids concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Outlet",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "S_va"],
        name="S_va concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Total valerate concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Outlet",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "S_bu"],
        name="S_bu concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Total butyrate concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Outlet",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "S_pro"],
        name="S_pro concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Total propionate concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Outlet",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "S_ac"],
        name="S_ac concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Total acetate concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Outlet",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "S_h2"],
        name="S_h2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Hydrogen gas concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Outlet",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "S_ch4"],
        name="S_ch4 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Methane gas concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Outlet",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "S_IC"],
        name="S_IC concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inorganic carbon concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Outlet",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "S_IN"],
        name="S_IN concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inorganic nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Outlet",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "S_I"],
        name="S_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Soluble inerts concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Outlet",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "X_c"],
        name="X_c concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Composites concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Outlet",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "X_ch"],
        name="X_ch concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Carbohydrates concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Outlet",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "X_pr"],
        name="X_pr concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Proteins concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Outlet",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "X_li"],
        name="X_li concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Lipids concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Outlet",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "X_su"],
        name="X_su concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Sugar degraders concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Outlet",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "X_aa"],
        name="X_aa concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Amino acid degraders concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Outlet",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "X_fa"],
        name="X_fa concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Long chain fatty acid (LCFA) degraders concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Outlet",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "X_c4"],
        name="X_c4 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Valerate and butyrate degraders concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Outlet",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "X_pro"],
        name="X_pro concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Propionate degraders concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Outlet",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "X_ac"],
        name="X_ac concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Acetate degraders concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Outlet",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "X_h2"],
        name="X_h2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Hydrogen degraders concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Outlet",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "X_I"],
        name="X_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Particulate inerts concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Outlet",
    )

    exports.add(
        obj=fs.adm_asm.inlet.flow_vol[0],
        name="Flow rate",
        ui_units=pyunits.m**3 / pyunits.day,
        display_units="m3/day",
        rounding=2,
        description="ADM1/ASM1 translator inlet flow rate",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Inlet",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "S_su"],
        name="S_su concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Monosaccharides concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Translator Outlet",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "S_aa"],
        name="S_aa concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Amino acids concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Inlet",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "S_fa"],
        name="S_fa concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Long chain fatty acids concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Inlet",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "S_va"],
        name="S_va concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Total valerate concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Inlet",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "S_bu"],
        name="S_bu concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Total butyrate concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Inlet",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "S_pro"],
        name="S_pro concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Total propionate concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Inlet",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "S_ac"],
        name="S_ac concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Total acetate concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Inlet",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "S_h2"],
        name="S_h2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Hydrogen gas concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Inlet",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "S_ch4"],
        name="S_ch4 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Methane gas concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Inlet",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "S_IC"],
        name="S_IC concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inorganic carbon concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Inlet",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "S_IN"],
        name="S_IN concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inorganic nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Inlet",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "S_I"],
        name="S_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Soluble inerts concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Inlet",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "X_c"],
        name="X_c concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Composites concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Inlet",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "X_ch"],
        name="X_ch concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Carbohydrates concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Inlet",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "X_pr"],
        name="X_pr concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Proteins concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Inlet",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "X_li"],
        name="X_li concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Lipids concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Inlet",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "X_su"],
        name="X_su concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Sugar degraders concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Inlet",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "X_aa"],
        name="X_aa concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Amino acid degraders concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Inlet",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "X_fa"],
        name="X_fa concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Long chain fatty acid (LCFA) degraders concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Inlet",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "X_c4"],
        name="X_c4 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Valerate and butyrate degraders concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Inlet",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "X_pro"],
        name="X_pro concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Propionate degraders concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Inlet",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "X_ac"],
        name="X_ac concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Acetate degraders concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Inlet",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "X_h2"],
        name="X_h2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Hydrogen degraders concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Inlet",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "X_I"],
        name="X_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Particulate inerts concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Inlet",
    )

    exports.add(
        obj=fs.adm_asm.outlet.flow_vol[0],
        name="Flow rate",
        ui_units=pyunits.m**3 / pyunits.day,
        display_units="m3/day",
        rounding=2,
        description="ADM1/ASM1 translator outlet flow rate",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Outlet",
    )
    exports.add(
        obj=fs.adm_asm.outlet.conc_mass_comp[0, "S_I"],
        name="S_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet soluble inert organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Outlet",
    )
    exports.add(
        obj=fs.adm_asm.outlet.conc_mass_comp[0, "S_S"],
        name="S_S concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet readily biodegradable substrate concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Outlet",
    )
    exports.add(
        obj=fs.adm_asm.outlet.conc_mass_comp[0, "X_I"],
        name="X_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet particulate inert organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Outlet",
    )
    exports.add(
        obj=fs.adm_asm.outlet.conc_mass_comp[0, "X_S"],
        name="X_S concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet slowly biodegradable substrate concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Outlet",
    )
    exports.add(
        obj=fs.adm_asm.outlet.conc_mass_comp[0, "X_BH"],
        name="X_BH concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet active heterotrophic biomass concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Outlet",
    )
    exports.add(
        obj=fs.adm_asm.outlet.conc_mass_comp[0, "X_BA"],
        name="X_BA concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet active autotrophic biomass concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Outlet",
    )
    exports.add(
        obj=fs.adm_asm.outlet.conc_mass_comp[0, "X_P"],
        name="X_P concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet particulate products arising from biomass decay concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Outlet",
    )
    exports.add(
        obj=fs.adm_asm.outlet.conc_mass_comp[0, "S_O"],
        name="S_O concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet oxygen concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Outlet",
    )
    exports.add(
        obj=fs.adm_asm.outlet.conc_mass_comp[0, "S_NO"],
        name="S_NO concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet nitrate and nitrite nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Outlet",
    )
    exports.add(
        obj=fs.adm_asm.outlet.conc_mass_comp[0, "S_NH"],
        name="S_NH concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet ammonium and ammonia nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Outlet",
    )
    exports.add(
        obj=fs.adm_asm.outlet.conc_mass_comp[0, "S_ND"],
        name="S_ND concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet soluble biodegradable organic nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Outlet",
    )
    exports.add(
        obj=fs.adm_asm.outlet.conc_mass_comp[0, "X_ND"],
        name="X_ND concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet particulate biodegradable organic nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Outlet",
    )
    exports.add(
        obj=fs.adm_asm.outlet.alkalinity[0],
        name="S_ALK concentration",
        ui_units=pyunits.mol / pyunits.m**3,
        display_units="mol/m3",
        rounding=5,
        description="Outlet alkalinity concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Translator Outlet",
    )

    exports.add(
        obj=fs.RADM.vapor_outlet.flow_vol[0],
        name="Flow rate",
        ui_units=pyunits.m**3 / pyunits.day,
        display_units="m3/day",
        rounding=2,
        description="Anaerobic digestor vapor outlet flow rate",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digestor Vapor Outlet",
    )
    exports.add(
        obj=fs.RADM.vapor_outlet.conc_mass_comp[0, "S_h2"],
        name="S_h2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Hydrogen gas concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digestor Vapor Outlet",
    )
    exports.add(
        obj=fs.RADM.vapor_outlet.conc_mass_comp[0, "S_ch4"],
        name="S_ch4 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Methane gas concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digestor Vapor Outlet",
    )
    exports.add(
        obj=fs.RADM.vapor_outlet.conc_mass_comp[0, "S_co2"],
        name="S_co2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Carbon dioxide gas concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digestor Vapor Outlet",
    )

    # performance metrics
    recovery_vol = (
        fs.Treated.properties[0].flow_vol / fs.FeedWater.properties[0].flow_vol
    )
    exports.add(
        obj=recovery_vol,
        name="Volumetric recovery",
        ui_units=pyunits.dimensionless,
        display_units="m3 of product/m3 of feed",
        rounding=5,
        description="Normalized volumetric recovery",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
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

    # add_costing(m)
    # assert_degrees_of_freedom(m, 0)
    # m.fs.costing.initialize()
    #
    # results = solve(m)
    # assert_optimal_termination(results)
    return m


def solve_flowsheet(flowsheet=None):
    fs = flowsheet
    results = solve(fs)
    return results
