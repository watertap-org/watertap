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
"""
GUI configuration for the base BSM2 flowsheet.
"""

from pyomo.environ import units as pyunits, assert_optimal_termination
from pyomo.util.check_units import assert_units_consistent

from watertap.ui.fsapi import FlowsheetInterface

from watertap.core.util.initialization import assert_degrees_of_freedom

from watertap.examples.flowsheets.case_studies.full_water_resource_recovery_facility.BSM2 import (
    build,
    set_operating_conditions,
    initialize_system,
    solve,
)


def export_to_ui():
    """
    Exports the variables, flowsheet build, and solver results to the GUI.
    """
    return FlowsheetInterface(
        name="BSM2",
        do_export=export_variables,
        do_build=build_flowsheet,
        do_solve=solve_flowsheet,
    )


def export_variables(flowsheet=None, exports=None, build_options=None, **kwargs):
    """
    Exports the variables to the GUI.
    """
    fs = flowsheet
    # --- Input data ---
    # Feed conditions
    exports.add(
        obj=fs.FeedWater.flow_vol[0],
        name="Feed volumetric flow rate",
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
        name="Feed S_I concentration",
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
        name="Feed S_S concentration",
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
        name="Feed X_I concentration",
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
        name="Feed X_S concentration",
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
        name="Feed X_BH concentration",
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
        name="Feed X_BA concentration",
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
        name="Feed X_P concentration",
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
        name="Feed S_O concentration",
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
        name="Feed S_NO concentration",
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
        name="Feed S_NH concentration",
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
        name="Feed S_ND concentration",
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
        name="Feed X_ND concentration",
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
        name="Feed S_ALK concentration",
        ui_units=pyunits.mol / pyunits.m**3,
        display_units="mol/m3",
        rounding=2,
        description="Alkalinity",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )

    # Unit model data, activated sludge process
    exports.add(
        obj=fs.R1.volume[0],
        name="First anoxic reactor volume",
        ui_units=pyunits.m**3,
        display_units="m3",
        rounding=1,
        description="CSTR volume",
        is_input=True,
        input_category="Activated sludge process",
        is_output=False,
    )
    exports.add(
        obj=fs.R2.volume[0],
        name="Second anoxic reactor volume",
        ui_units=pyunits.m**3,
        display_units="m3",
        rounding=1,
        description="CSTR volume",
        is_input=True,
        input_category="Activated sludge process",
        is_output=False,
    )
    exports.add(
        obj=fs.R3.volume[0],
        name="First aerobic reactor volume",
        ui_units=pyunits.m**3,
        display_units="m3",
        rounding=1,
        description="CSTR volume",
        is_input=True,
        input_category="Activated sludge process",
        is_output=False,
    )
    exports.add(
        obj=fs.R4.volume[0],
        name="Second aerobic reactor volume",
        ui_units=pyunits.m**3,
        display_units="m3",
        rounding=1,
        description="CSTR volume",
        is_input=True,
        input_category="Activated sludge process",
        is_output=False,
    )
    exports.add(
        obj=fs.R5.volume[0],
        name="Third aerobic reactor volume",
        ui_units=pyunits.m**3,
        display_units="m3",
        rounding=1,
        description="CSTR volume",
        is_input=True,
        input_category="Activated sludge process",
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

    # TODO: uncomment and revise below once costing is merged
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
        name="Treated flow rate",
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
        name="Treated S_I concentration",
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
        name="Treated S_S concentration",
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
        name="Treated X_I concentration",
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
        name="Treated X_S concentration",
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
        name="Treated X_BH concentration",
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
        name="Treated X_BA concentration",
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
        name="Treated X_P concentration",
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
        name="Treated S_O concentration",
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
        name="Treated S_NO concentration",
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
        name="Treated S_NH concentration",
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
        name="Treated S_ND concentration",
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
        name="Treated X_ND concentration",
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
        name="Treated S_ALK concentration",
        ui_units=pyunits.mol / pyunits.m**3,
        display_units="mol/m3",
        rounding=5,
        description="Outlet alkalinity concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Effluent",
    )

    exports.add(
        obj=fs.CL.effluent.flow_vol[0],
        name="Primary clarifier effluent flow rate",
        ui_units=pyunits.m**3 / pyunits.day,
        display_units="m3/day",
        rounding=2,
        description="Outlet primary clarifier flow rate",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Effluent",
    )
    exports.add(
        obj=fs.CL.effluent.conc_mass_comp[0, "S_I"],
        name="Primary clarifier effluent S_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet soluble inert organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Effluent",
    )
    exports.add(
        obj=fs.CL.effluent.conc_mass_comp[0, "S_S"],
        name="Primary clarifier effluent S_S concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet readily biodegradable substrate concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Effluent",
    )
    exports.add(
        obj=fs.CL.effluent.conc_mass_comp[0, "X_I"],
        name="Primary clarifier effluent X_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet particulate inert organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Effluent",
    )
    exports.add(
        obj=fs.CL.effluent.conc_mass_comp[0, "X_S"],
        name="Primary clarifier effluent X_S concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet slowly biodegradable substrate concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Effluent",
    )
    exports.add(
        obj=fs.CL.effluent.conc_mass_comp[0, "X_BH"],
        name="Primary clarifier effluent X_BH concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet active heterotrophic biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Effluent",
    )
    exports.add(
        obj=fs.CL.effluent.conc_mass_comp[0, "X_BA"],
        name="Primary clarifier effluent X_BA concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet active autotrophic biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Effluent",
    )
    exports.add(
        obj=fs.CL.effluent.conc_mass_comp[0, "X_P"],
        name="Primary clarifier effluent X_P concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet particulate products arising from biomass decay concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Effluent",
    )
    exports.add(
        obj=fs.CL.effluent.conc_mass_comp[0, "S_O"],
        name="Primary clarifier effluent S_O concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet oxygen concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Effluent",
    )
    exports.add(
        obj=fs.CL.effluent.conc_mass_comp[0, "S_NO"],
        name="Primary clarifier effluent S_NO concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet nitrate and nitrite nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Effluent",
    )
    exports.add(
        obj=fs.CL.effluent.conc_mass_comp[0, "S_NH"],
        name="Primary clarifier effluent S_NH concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet ammonium and ammonia nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Effluent",
    )
    exports.add(
        obj=fs.CL.effluent.conc_mass_comp[0, "S_ND"],
        name="Primary clarifier effluent S_ND concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet soluble biodegradable organic nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Effluent",
    )
    exports.add(
        obj=fs.CL.effluent.conc_mass_comp[0, "X_ND"],
        name="Primary clarifier effluent X_ND concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet particulate biodegradable organic nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Effluent",
    )
    exports.add(
        obj=fs.CL.effluent.alkalinity[0],
        name="Primary clarifier effluent S_ALK concentration",
        ui_units=pyunits.mol / pyunits.m**3,
        display_units="mol/m3",
        rounding=5,
        description="Outlet alkalinity concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Effluent",
    )

    exports.add(
        obj=fs.CL.underflow.flow_vol[0],
        name="Primary clarifier underflow flow rate",
        ui_units=pyunits.m**3 / pyunits.day,
        display_units="m3/day",
        rounding=2,
        description="Outlet primary clarifier flow rate",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Underflow",
    )
    exports.add(
        obj=fs.CL.underflow.conc_mass_comp[0, "S_I"],
        name="Primary clarifier underflow S_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet soluble inert organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Underflow",
    )
    exports.add(
        obj=fs.CL.underflow.conc_mass_comp[0, "S_S"],
        name="Primary clarifier underflow S_S concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet readily biodegradable substrate concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Underflow",
    )
    exports.add(
        obj=fs.CL.underflow.conc_mass_comp[0, "X_I"],
        name="8Primary clarifier underflow X_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet particulate inert organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Underflow",
    )
    exports.add(
        obj=fs.CL.underflow.conc_mass_comp[0, "X_S"],
        name="Primary clarifier underflow X_S concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet slowly biodegradable substrate concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Underflow",
    )
    exports.add(
        obj=fs.CL.underflow.conc_mass_comp[0, "X_BH"],
        name="Primary clarifier underflow X_BH concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet active heterotrophic biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Underflow",
    )
    exports.add(
        obj=fs.CL.underflow.conc_mass_comp[0, "X_BA"],
        name="Primary clarifier underflow X_BA concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet active autotrophic biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Underflow",
    )
    exports.add(
        obj=fs.CL.underflow.conc_mass_comp[0, "X_P"],
        name="Primary clarifier underflow X_P concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet particulate products arising from biomass decay concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Underflow",
    )
    exports.add(
        obj=fs.CL.underflow.conc_mass_comp[0, "S_O"],
        name="Primary clarifier underflow S_O concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet oxygen concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Underflow",
    )
    exports.add(
        obj=fs.CL.underflow.conc_mass_comp[0, "S_NO"],
        name="Primary clarifier underflow S_NO concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet nitrate and nitrite nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Underflow",
    )
    exports.add(
        obj=fs.CL.underflow.conc_mass_comp[0, "S_NH"],
        name="Primary clarifier underflow S_NH concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet ammonium and ammonia nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Underflow",
    )
    exports.add(
        obj=fs.CL.underflow.conc_mass_comp[0, "S_ND"],
        name="Primary clarifier underflow S_ND concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet soluble biodegradable organic nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Underflow",
    )
    exports.add(
        obj=fs.CL.underflow.conc_mass_comp[0, "X_ND"],
        name="Primary clarifier underflow X_ND concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet particulate biodegradable organic nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Underflow",
    )
    exports.add(
        obj=fs.CL.underflow.alkalinity[0],
        name="Primary clarifier underflow S_ALK concentration",
        ui_units=pyunits.mol / pyunits.m**3,
        display_units="mol/m3",
        rounding=5,
        description="Outlet alkalinity concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Underflow",
    )

    exports.add(
        obj=fs.SP6.recycle.flow_vol[0],
        name="ASP recycle inlet flow rate",
        ui_units=pyunits.m**3 / pyunits.day,
        display_units="m3/day",
        rounding=2,
        description="ASP recycle flow rate",
        is_input=False,
        is_output=True,
        output_category="ASP Recycle Inlet",
    )
    exports.add(
        obj=fs.SP6.recycle.conc_mass_comp[0, "S_I"],
        name="ASP recycle inlet S_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet soluble inert organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="ASP Recycle Inlet",
    )
    exports.add(
        obj=fs.SP6.recycle.conc_mass_comp[0, "S_S"],
        name="ASP recycle inlet S_S concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet readily biodegradable substrate concentration",
        is_input=False,
        is_output=True,
        output_category="ASP Recycle Inlet",
    )
    exports.add(
        obj=fs.SP6.recycle.conc_mass_comp[0, "X_I"],
        name="ASP recycle inlet X_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet particulate inert organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="ASP Recycle Inlet",
    )
    exports.add(
        obj=fs.SP6.recycle.conc_mass_comp[0, "X_S"],
        name="ASP recycle inlet X_S concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet slowly biodegradable substrate concentration",
        is_input=False,
        is_output=True,
        output_category="ASP Recycle Inlet",
    )
    exports.add(
        obj=fs.SP6.recycle.conc_mass_comp[0, "X_BH"],
        name="ASP recycle inlet X_BH concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet active heterotrophic biomass concentration",
        is_input=False,
        is_output=True,
        output_category="ASP Recycle Inlet",
    )
    exports.add(
        obj=fs.SP6.recycle.conc_mass_comp[0, "X_BA"],
        name="ASP recycle inlet X_BA concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet active autotrophic biomass concentration",
        is_input=False,
        is_output=True,
        output_category="ASP Recycle Inlet",
    )
    exports.add(
        obj=fs.SP6.recycle.conc_mass_comp[0, "X_P"],
        name="ASP recycle inlet X_P concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet particulate products arising from biomass decay concentration",
        is_input=False,
        is_output=True,
        output_category="ASP Recycle Inlet",
    )
    exports.add(
        obj=fs.SP6.recycle.conc_mass_comp[0, "S_O"],
        name="ASP recycle inlet S_O concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet oxygen concentration",
        is_input=False,
        is_output=True,
        output_category="ASP Recycle Inlet",
    )
    exports.add(
        obj=fs.SP6.recycle.conc_mass_comp[0, "S_NO"],
        name="ASP recycle inlet S_NO concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet nitrate and nitrite nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="ASP Recycle Inlet",
    )
    exports.add(
        obj=fs.SP6.recycle.conc_mass_comp[0, "S_NH"],
        name="ASP recycle inlet S_NH concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet ammonium and ammonia nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="ASP Recycle Inlet",
    )
    exports.add(
        obj=fs.SP6.recycle.conc_mass_comp[0, "S_ND"],
        name="ASP recycle inlet S_ND concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet soluble biodegradable organic nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="ASP Recycle Inlet",
    )
    exports.add(
        obj=fs.SP6.recycle.conc_mass_comp[0, "X_ND"],
        name="ASP recycle inlet X_ND concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet particulate biodegradable organic nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="ASP Recycle Inlet",
    )
    exports.add(
        obj=fs.SP6.recycle.alkalinity[0],
        name="ASP recycle inlet S_ALK concentration",
        ui_units=pyunits.mol / pyunits.m**3,
        display_units="mol/m3",
        rounding=5,
        description="Inlet alkalinity concentration",
        is_input=False,
        is_output=True,
        output_category="ASP Recycle Inlet",
    )

    exports.add(
        obj=fs.SP5.overflow.flow_vol[0],
        name="Secondary clarifier inlet flow rate",
        ui_units=pyunits.m**3 / pyunits.day,
        display_units="m3/day",
        rounding=2,
        description="Secondary clarifier inlet flow rate",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Inlet",
    )
    exports.add(
        obj=fs.SP5.overflow.conc_mass_comp[0, "S_I"],
        name="Secondary clarifier inlet S_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet soluble inert organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Inlet",
    )
    exports.add(
        obj=fs.SP5.overflow.conc_mass_comp[0, "S_S"],
        name="Secondary clarifier inlet S_S concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet readily biodegradable substrate concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Inlet",
    )
    exports.add(
        obj=fs.SP5.overflow.conc_mass_comp[0, "X_I"],
        name="Secondary clarifier inlet X_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet particulate inert organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Inlet",
    )
    exports.add(
        obj=fs.SP5.overflow.conc_mass_comp[0, "X_S"],
        name="Secondary clarifier inlet X_S concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet slowly biodegradable substrate concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Inlet",
    )
    exports.add(
        obj=fs.SP5.overflow.conc_mass_comp[0, "X_BH"],
        name="Secondary clarifier inlet X_BH concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet active heterotrophic biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Inlet",
    )
    exports.add(
        obj=fs.SP5.overflow.conc_mass_comp[0, "X_BA"],
        name="Secondary clarifier inlet X_BA concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet active autotrophic biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Inlet",
    )
    exports.add(
        obj=fs.SP5.overflow.conc_mass_comp[0, "X_P"],
        name="Secondary clarifier inlet X_P concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet particulate products arising from biomass decay concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Inlet",
    )
    exports.add(
        obj=fs.SP5.overflow.conc_mass_comp[0, "S_O"],
        name="Secondary clarifier inlet S_O concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet oxygen concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Inlet",
    )
    exports.add(
        obj=fs.SP5.overflow.conc_mass_comp[0, "S_NO"],
        name="Secondary clarifier inlet S_NO concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet nitrate and nitrite nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Inlet",
    )
    exports.add(
        obj=fs.SP5.overflow.conc_mass_comp[0, "S_NH"],
        name="Secondary clarifier inlet S_NH concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet ammonium and ammonia nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Inlet",
    )
    exports.add(
        obj=fs.SP5.overflow.conc_mass_comp[0, "S_ND"],
        name="Secondary clarifier inlet S_ND concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet soluble biodegradable organic nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Inlet",
    )
    exports.add(
        obj=fs.SP5.overflow.conc_mass_comp[0, "X_ND"],
        name="Secondary clarifier inlet X_ND concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet particulate biodegradable organic nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Inlet",
    )
    exports.add(
        obj=fs.SP5.overflow.alkalinity[0],
        name="Secondary clarifier inlet S_ALK concentration",
        ui_units=pyunits.mol / pyunits.m**3,
        display_units="mol/m3",
        rounding=5,
        description="Inlet alkalinity concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Inlet",
    )

    exports.add(
        obj=fs.TU.inlet.flow_vol[0],
        name="Thickener inlet flow rate",
        ui_units=pyunits.m**3 / pyunits.day,
        display_units="m3/day",
        rounding=2,
        description="Thickener inlet flow rate",
        is_input=False,
        is_output=True,
        output_category="Thickener Inlet",
    )
    exports.add(
        obj=fs.TU.inlet.conc_mass_comp[0, "S_I"],
        name="Thickener inlet S_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet soluble inert organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="Thickener Inlet",
    )
    exports.add(
        obj=fs.TU.inlet.conc_mass_comp[0, "S_S"],
        name="Thickener inlet S_S concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet readily biodegradable substrate concentration",
        is_input=False,
        is_output=True,
        output_category="Thickener Inlet",
    )
    exports.add(
        obj=fs.TU.inlet.conc_mass_comp[0, "X_I"],
        name="Thickener inlet X_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet particulate inert organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="Thickener Inlet",
    )
    exports.add(
        obj=fs.TU.inlet.conc_mass_comp[0, "X_S"],
        name="Thickener inlet X_S concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet slowly biodegradable substrate concentration",
        is_input=False,
        is_output=True,
        output_category="Thickener Inlet",
    )
    exports.add(
        obj=fs.TU.inlet.conc_mass_comp[0, "X_BH"],
        name="Thickener inlet X_BH concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet active heterotrophic biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Thickener Inlet",
    )
    exports.add(
        obj=fs.TU.inlet.conc_mass_comp[0, "X_BA"],
        name="Thickener inlet X_BA concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet active autotrophic biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Thickener Inlet",
    )
    exports.add(
        obj=fs.TU.inlet.conc_mass_comp[0, "X_P"],
        name="Thickener inlet X_P concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet particulate products arising from biomass decay concentration",
        is_input=False,
        is_output=True,
        output_category="Thickener Inlet",
    )
    exports.add(
        obj=fs.TU.inlet.conc_mass_comp[0, "S_O"],
        name="Thickener inlet S_O concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet oxygen concentration",
        is_input=False,
        is_output=True,
        output_category="Thickener Inlet",
    )
    exports.add(
        obj=fs.TU.inlet.conc_mass_comp[0, "S_NO"],
        name="Thickener inlet S_NO concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet nitrate and nitrite nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Thickener Inlet",
    )
    exports.add(
        obj=fs.TU.inlet.conc_mass_comp[0, "S_NH"],
        name="Thickener inlet S_NH concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet ammonium and ammonia nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Thickener Inlet",
    )
    exports.add(
        obj=fs.TU.inlet.conc_mass_comp[0, "S_ND"],
        name="Thickener inlet S_ND concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet soluble biodegradable organic nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Thickener Inlet",
    )
    exports.add(
        obj=fs.TU.inlet.conc_mass_comp[0, "X_ND"],
        name="Thickener inlet X_ND concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet particulate biodegradable organic nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Thickener Inlet",
    )
    exports.add(
        obj=fs.TU.inlet.alkalinity[0],
        name="Thickener inlet S_ALK concentration",
        ui_units=pyunits.mol / pyunits.m**3,
        display_units="mol/m3",
        rounding=5,
        description="Inlet alkalinity concentration",
        is_input=False,
        is_output=True,
        output_category="Thickener Inlet",
    )

    exports.add(
        obj=fs.asm_adm.inlet.flow_vol[0],
        name="ASM1/ADM1 inlet flow rate",
        ui_units=pyunits.m**3 / pyunits.day,
        display_units="m3/day",
        rounding=2,
        description="ASM1/ADM1 interface inlet flow rate",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Interface Inlet",
    )
    exports.add(
        obj=fs.asm_adm.inlet.conc_mass_comp[0, "S_I"],
        name="ASM1/ADM1 inlet S_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet soluble inert organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Interface Inlet",
    )
    exports.add(
        obj=fs.asm_adm.inlet.conc_mass_comp[0, "S_S"],
        name="ASM1/ADM1 inlet S_S concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet readily biodegradable substrate concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Interface Inlet",
    )
    exports.add(
        obj=fs.asm_adm.inlet.conc_mass_comp[0, "X_I"],
        name="ASM1/ADM1 inlet X_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet particulate inert organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Interface Inlet",
    )
    exports.add(
        obj=fs.asm_adm.inlet.conc_mass_comp[0, "X_S"],
        name="ASM1/ADM1 inlet X_S concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet slowly biodegradable substrate concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Interface Inlet",
    )
    exports.add(
        obj=fs.asm_adm.inlet.conc_mass_comp[0, "X_BH"],
        name="ASM1/ADM1 inlet X_BH concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet active heterotrophic biomass concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Interface Inlet",
    )
    exports.add(
        obj=fs.asm_adm.inlet.conc_mass_comp[0, "X_BA"],
        name="ASM1/ADM1 inlet X_BA concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet active autotrophic biomass concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Interface Inlet",
    )
    exports.add(
        obj=fs.asm_adm.inlet.conc_mass_comp[0, "X_P"],
        name="ASM1/ADM1 inlet X_P concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet particulate products arising from biomass decay concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Interface Inlet",
    )
    exports.add(
        obj=fs.asm_adm.inlet.conc_mass_comp[0, "S_O"],
        name="ASM1/ADM1 inlet S_O concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet oxygen concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Interface Inlet",
    )
    exports.add(
        obj=fs.asm_adm.inlet.conc_mass_comp[0, "S_NO"],
        name="ASM1/ADM1 inlet S_NO concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet nitrate and nitrite nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Interface Inlet",
    )
    exports.add(
        obj=fs.asm_adm.inlet.conc_mass_comp[0, "S_NH"],
        name="ASM1/ADM1 inlet S_NH concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet ammonium and ammonia nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Interface Inlet",
    )
    exports.add(
        obj=fs.asm_adm.inlet.conc_mass_comp[0, "S_ND"],
        name="ASM1/ADM1 inlet S_ND concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet soluble biodegradable organic nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Interface Inlet",
    )
    exports.add(
        obj=fs.asm_adm.inlet.conc_mass_comp[0, "X_ND"],
        name="ASM1/ADM1 inlet X_ND concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inlet particulate biodegradable organic nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Interface Inlet",
    )
    exports.add(
        obj=fs.asm_adm.inlet.alkalinity[0],
        name="ASM1/ADM1 inlet S_ALK concentration",
        ui_units=pyunits.mol / pyunits.m**3,
        display_units="mol/m3",
        rounding=5,
        description="Inlet alkalinity concentration",
        is_input=False,
        is_output=True,
        output_category="ASM1/ADM1 Interface Inlet",
    )

    exports.add(
        obj=fs.asm_adm.outlet.flow_vol[0],
        name="Anaerobic digester inlet flow rate",
        ui_units=pyunits.m**3 / pyunits.day,
        display_units="m3/day",
        rounding=2,
        description="Anaerobic digester inlet flow rate",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "S_su"],
        name="Anaerobic digester inlet S_su concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Monosaccharides concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "S_aa"],
        name="Anaerobic digester inlet S_aa concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Amino acids concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "S_fa"],
        name="Anaerobic digester inlet S_fa concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Long chain fatty acids concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "S_va"],
        name="Anaerobic digester inlet S_va concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Total valerate concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "S_bu"],
        name="Anaerobic digester inlet S_bu concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Total butyrate concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "S_pro"],
        name="Anaerobic digester inlet S_pro concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Total propionate concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "S_ac"],
        name="Anaerobic digester inlet S_ac concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Total acetate concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "S_h2"],
        name="Anaerobic digester inlet S_h2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Hydrogen gas concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "S_ch4"],
        name="Anaerobic digester inlet S_ch4 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Methane gas concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "S_IC"],
        name="Anaerobic digester inlet S_IC concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inorganic carbon concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "S_IN"],
        name="Anaerobic digester inlet S_IN concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inorganic nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "S_I"],
        name="Anaerobic digester inlet S_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Soluble inerts concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "X_c"],
        name="Anaerobic digester inlet X_c concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Composites concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "X_ch"],
        name="Anaerobic digester inlet X_ch concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Carbohydrates concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "X_pr"],
        name="Anaerobic digester inlet X_pr concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Proteins concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "X_li"],
        name="Anaerobic digester inlet X_li concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Lipids concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "X_su"],
        name="Anaerobic digester inlet X_su concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Sugar degraders concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "X_aa"],
        name="Anaerobic digester inlet X_aa concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Amino acid degraders concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "X_fa"],
        name="Anaerobic digester inlet X_fa concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Long chain fatty acid (LCFA) degraders concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "X_c4"],
        name="Anaerobic digester inlet X_c4 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Valerate and butyrate degraders concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "X_pro"],
        name="Anaerobic digester inlet X_pro concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Propionate degraders concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "X_ac"],
        name="Anaerobic digester inlet X_ac concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Acetate degraders concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "X_h2"],
        name="Anaerobic digester inlet X_h2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Hydrogen degraders concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.asm_adm.outlet.conc_mass_comp[0, "X_I"],
        name="Anaerobic digester inlet X_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Particulate inerts concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )

    exports.add(
        obj=fs.adm_asm.inlet.flow_vol[0],
        name="Anaerobic digester liquid outlet flow rate",
        ui_units=pyunits.m**3 / pyunits.day,
        display_units="m3/day",
        rounding=2,
        description="Anaerobic digester liquid outlet flow rate",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Liquid Outlet (pre-interface)",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "S_su"],
        name="Anaerobic digester liquid outlet S_su concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Monosaccharides concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Liquid Outlet (pre-interface)",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "S_aa"],
        name="Anaerobic digester liquid outlet S_aa concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Amino acids concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Liquid Outlet (pre-interface)",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "S_fa"],
        name="Anaerobic digester liquid outlet S_fa concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Long chain fatty acids concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Liquid Outlet (pre-interface)",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "S_va"],
        name="Anaerobic digester liquid outlet S_va concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Total valerate concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Liquid Outlet (pre-interface)",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "S_bu"],
        name="Anaerobic digester liquid outlet S_bu concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Total butyrate concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Liquid Outlet (pre-interface)",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "S_pro"],
        name="Anaerobic digester liquid outlet S_pro concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Total propionate concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Liquid Outlet (pre-interface)",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "S_ac"],
        name="Anaerobic digester liquid outlet S_ac concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Total acetate concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Liquid Outlet (pre-interface)",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "S_h2"],
        name="Anaerobic digester liquid outlet S_h2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Hydrogen gas concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Liquid Outlet (pre-interface)",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "S_ch4"],
        name="Anaerobic digester liquid outlet S_ch4 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Methane gas concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Liquid Outlet (pre-interface)",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "S_IC"],
        name="Anaerobic digester liquid outlet S_IC concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inorganic carbon concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Liquid Outlet (pre-interface)",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "S_IN"],
        name="Anaerobic digester liquid outlet S_IN concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inorganic nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Liquid Outlet (pre-interface)",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "S_I"],
        name="Anaerobic digester liquid outlet S_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Soluble inerts concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Liquid Outlet (pre-interface)",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "X_c"],
        name="Anaerobic digester liquid outlet X_c concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Composites concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Liquid Outlet (pre-interface)",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "X_ch"],
        name="Anaerobic digester liquid outlet X_ch concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Carbohydrates concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Liquid Outlet (pre-interface)",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "X_pr"],
        name="Anaerobic digester liquid outlet X_pr concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Proteins concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Liquid Outlet (pre-interface)",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "X_li"],
        name="Anaerobic digester liquid outlet X_li concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Lipids concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Liquid Outlet (pre-interface)",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "X_su"],
        name="Anaerobic digester liquid outlet X_su concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Sugar degraders concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Liquid Outlet (pre-interface)",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "X_aa"],
        name="Anaerobic digester liquid outlet X_aa concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Amino acid degraders concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Liquid Outlet (pre-interface)",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "X_fa"],
        name="Anaerobic digester liquid outlet X_fa concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Long chain fatty acid (LCFA) degraders concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Liquid Outlet (pre-interface)",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "X_c4"],
        name="Anaerobic digester liquid outlet X_c4 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Valerate and butyrate degraders concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Liquid Outlet (pre-interface)",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "X_pro"],
        name="Anaerobic digester liquid outlet X_pro concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Propionate degraders concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Liquid Outlet (pre-interface)",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "X_ac"],
        name="Anaerobic digester liquid outlet X_ac concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Acetate degraders concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Liquid Outlet (pre-interface)",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "X_h2"],
        name="Anaerobic digester liquid outlet X_h2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Hydrogen degraders concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Liquid Outlet (pre-interface)",
    )
    exports.add(
        obj=fs.adm_asm.inlet.conc_mass_comp[0, "X_I"],
        name="Anaerobic digester liquid outlet X_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Particulate inerts concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Liquid Outlet (pre-interface)",
    )

    exports.add(
        obj=fs.adm_asm.outlet.flow_vol[0],
        name="ADM1/ASM1 outlet flow rate",
        ui_units=pyunits.m**3 / pyunits.day,
        display_units="m3/day",
        rounding=2,
        description="ADM1/ASM1 interface outlet flow rate",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Interface Outlet",
    )
    exports.add(
        obj=fs.adm_asm.outlet.conc_mass_comp[0, "S_I"],
        name="ADM1/ASM1 outlet S_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet soluble inert organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Interface Outlet",
    )
    exports.add(
        obj=fs.adm_asm.outlet.conc_mass_comp[0, "S_S"],
        name="ADM1/ASM1 outlet S_S concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet readily biodegradable substrate concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Interface Outlet",
    )
    exports.add(
        obj=fs.adm_asm.outlet.conc_mass_comp[0, "X_I"],
        name="ADM1/ASM1 outlet X_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet particulate inert organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Interface Outlet",
    )
    exports.add(
        obj=fs.adm_asm.outlet.conc_mass_comp[0, "X_S"],
        name="ADM1/ASM1 outlet X_S concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet slowly biodegradable substrate concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Interface Outlet",
    )
    exports.add(
        obj=fs.adm_asm.outlet.conc_mass_comp[0, "X_BH"],
        name="ADM1/ASM1 outlet X_BH concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet active heterotrophic biomass concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Interface Outlet",
    )
    exports.add(
        obj=fs.adm_asm.outlet.conc_mass_comp[0, "X_BA"],
        name="ADM1/ASM1 outlet X_BA concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet active autotrophic biomass concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Interface Outlet",
    )
    exports.add(
        obj=fs.adm_asm.outlet.conc_mass_comp[0, "X_P"],
        name="ADM1/ASM1 outlet X_P concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet particulate products arising from biomass decay concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Interface Outlet",
    )
    exports.add(
        obj=fs.adm_asm.outlet.conc_mass_comp[0, "S_O"],
        name="ADM1/ASM1 outlet S_O concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet oxygen concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Interface Outlet",
    )
    exports.add(
        obj=fs.adm_asm.outlet.conc_mass_comp[0, "S_NO"],
        name="ADM1/ASM1 outlet S_NO concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet nitrate and nitrite nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Interface Outlet",
    )
    exports.add(
        obj=fs.adm_asm.outlet.conc_mass_comp[0, "S_NH"],
        name="ADM1/ASM1 outlet S_NH concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet ammonium and ammonia nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Interface Outlet",
    )
    exports.add(
        obj=fs.adm_asm.outlet.conc_mass_comp[0, "S_ND"],
        name="ADM1/ASM1 outlet S_ND concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet soluble biodegradable organic nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Interface Outlet",
    )
    exports.add(
        obj=fs.adm_asm.outlet.conc_mass_comp[0, "X_ND"],
        name="ADM1/ASM1 outlet X_ND concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet particulate biodegradable organic nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Interface Outlet",
    )
    exports.add(
        obj=fs.adm_asm.outlet.alkalinity[0],
        name="ADM1/ASM1 outlet S_ALK concentration",
        ui_units=pyunits.mol / pyunits.m**3,
        display_units="mol/m3",
        rounding=5,
        description="Outlet alkalinity concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM1 Interface Outlet",
    )

    exports.add(
        obj=fs.RADM.vapor_outlet.flow_vol[0],
        name="Anaerobic digester vapor flow rate",
        ui_units=pyunits.m**3 / pyunits.day,
        display_units="m3/day",
        rounding=2,
        description="Anaerobic digester vapor outlet flow rate",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Vapor Outlet",
    )
    exports.add(
        obj=fs.RADM.vapor_outlet.conc_mass_comp[0, "S_h2"],
        name="Anaerobic digester vapor S_h2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Hydrogen gas concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Vapor Outlet",
    )
    exports.add(
        obj=fs.RADM.vapor_outlet.conc_mass_comp[0, "S_ch4"],
        name="Anaerobic digester vapor S_ch4 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Methane gas concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Vapor Outlet",
    )
    exports.add(
        obj=fs.RADM.vapor_outlet.conc_mass_comp[0, "S_co2"],
        name="Anaerobic digester vapor S_co2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Carbon dioxide gas concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Vapor Outlet",
    )

    exports.add(
        obj=fs.DU.underflow.flow_vol[0],
        name="Dewatered sludge flow rate",
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
        name="Dewatered sludge S_I concentration",
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
        name="Dewatered sludge S_S concentration",
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
        name="Dewatered sludge X_I concentration",
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
        name="Dewatered sludge X_S concentration",
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
        name="Dewatered sludge X_BH concentration",
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
        name="Dewatered sludge X_BA concentration",
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
        name="Dewatered sludge X_P concentration",
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
        name="Dewatered sludge S_O concentration",
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
        name="Dewatered sludge S_NO concentration",
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
        name="Dewatered sludge S_NH concentration",
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
        name="Dewatered sludge S_ND concentration",
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
        name="Dewatered sludge X_ND concentration",
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
        name="Dewatered sludge S_ALK concentration",
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
        name="Dewatering unit liquid outlet flow rate",
        ui_units=pyunits.m**3 / pyunits.day,
        display_units="m3/day",
        rounding=2,
        description="Dewatering unit overflow flow rate",
        is_input=False,
        is_output=True,
        output_category="Dewatering Unit Liquid Outlet",
    )
    exports.add(
        obj=fs.DU.overflow.conc_mass_comp[0, "S_I"],
        name="Dewatering unit liquid outlet S_I concentration",
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
        name="Dewatering unit liquid outlet S_S concentration",
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
        name="Dewatering unit liquid outlet X_I concentration",
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
        name="Dewatering unit liquid outlet X_S concentration",
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
        name="Dewatering unit liquid outlet X_BH concentration",
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
        name="Dewatering unit liquid outlet X_BA concentration",
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
        name="Dewatering unit liquid outlet X_P concentration",
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
        name="Dewatering unit liquid outlet S_O concentration",
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
        name="Dewatering unit liquid outlet S_NO concentration",
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
        name="Dewatering unit liquid outlet S_NH concentration",
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
        name="Dewatering unit liquid outlet S_ND concentration",
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
        name="Dewatering unit liquid outlet X_ND concentration",
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
        name="Dewatering unit liquid outlet S_ALK concentration",
        ui_units=pyunits.mol / pyunits.m**3,
        display_units="mol/m3",
        rounding=5,
        description="Outlet alkalinity concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatering Unit Liquid Outlet",
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


def build_flowsheet(build_options=None, **kwargs):
    """
    Builds the initial flowsheet.
    """
    m = build()

    set_operating_conditions(m)
    assert_degrees_of_freedom(m, 0)
    assert_units_consistent(m)

    initialize_system(m)

    results = solve(m)
    assert_optimal_termination(results)

    # TODO: incorporate costing when merged
    # add_costing(m)
    # assert_degrees_of_freedom(m, 0)
    # m.fs.costing.initialize()
    #
    # results = solve(m)
    # assert_optimal_termination(results)
    return m


def solve_flowsheet(flowsheet=None):
    """
    Solves the initial flowsheet.
    """
    fs = flowsheet
    results = solve(fs)
    return results
