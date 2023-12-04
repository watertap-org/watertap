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
GUI configuration for the base GAC model.
"""

from pyomo.environ import units as pyunits
from watertap.examples.flowsheets.case_studies.gac_basic import (
    gac_basic_flowsheet as gac_fs,
)
from watertap.ui.fsapi import FlowsheetInterface


def export_to_ui():
    """
    Exports the variables, flowsheet build, and solver results to the GUI.
    """
    return FlowsheetInterface(
        name="GAC",
        do_export=export_variables,
        do_build=build_flowsheet,
        do_solve=solve_flowsheet,
    )


def export_variables(flowsheet=None, exports=None, build_options=None, **kwargs):
    """
    Exports the variables to the GUI.
    """
    fs = flowsheet
    rounding = 3
    solute_name = "solute"

    # input data
    # ---------------------------------------------------------------------
    # feed conditions
    category = "Feed"
    exports.add(
        obj=fs.feed.properties[0].temperature,
        name="Feed temperature",
        ui_units=pyunits.kelvin,
        display_units="K",
        rounding=rounding,
        description="Inlet temperature",
        is_input=True,
        input_category=category,
        is_output=False,
    )
    exports.add(
        obj=fs.feed.properties[0].pressure,
        name="Feed pressure",
        ui_units=pyunits.bar,
        display_units="bar",
        rounding=rounding,
        description="Inlet pressure",
        is_input=True,
        input_category=category,
        is_output=False,
    )
    exports.add(
        obj=fs.feed.properties[0].flow_vol_phase["Liq"],
        name="Feed volumetric flow rate",
        ui_units=pyunits.Mgallon / pyunits.day,
        display_units="MGD",
        rounding=rounding,
        description="Inlet volumetric flow rate",
        is_input=False,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.feed.properties[0].conc_mass_phase_comp["Liq", solute_name],
        name="Feed solute concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=rounding,
        description="Inlet solute concentration",
        is_input=False,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.feed.properties[0].flow_mol_phase_comp["Liq", "H2O"],
        name="Feed water molar flow rate",
        ui_units=pyunits.mol / pyunits.s,
        display_units="mol/s",
        rounding=rounding,
        description="Inlet water molar flow rate",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.feed.properties[0].flow_mol_phase_comp["Liq", solute_name],
        name="Feed solute molar flow rate",
        ui_units=pyunits.mol / pyunits.s,
        display_units="mol/s",
        rounding=rounding,
        description="Inlet solute molar flow rate",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    # ---------------------------------------------------------------------
    # product conditions
    category = "Treated"
    exports.add(
        obj=fs.product.properties[0].flow_vol_phase["Liq"],
        name="Treated water volumetric flow rate",
        ui_units=pyunits.Mgallon / pyunits.day,
        display_units="MGD",
        rounding=rounding,
        description="Outlet volumetric flow rate",
        is_input=False,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.product.properties[0].conc_mass_phase_comp["Liq", solute_name],
        name="Treated water solute concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=rounding,
        description="Outlet solute concentration",
        is_input=False,
        is_output=True,
        output_category=category,
    )
    # ---------------------------------------------------------------------
    # gac media
    category = "GAC media"
    exports.add(
        obj=fs.gac.particle_dens_app,
        name="Media apparent density",
        ui_units=pyunits.kg * pyunits.m**-3,
        display_units="kg/m3",
        rounding=rounding,
        description="Apparent density of the GAC media",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.gac.particle_dens_bulk,
        name="Media bulk density in the bed",
        ui_units=pyunits.kg * pyunits.m**-3,
        display_units="kg/m3",
        rounding=rounding,
        description="Media bulk density in the bed",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.gac.particle_dia,
        name="Media particle diameter",
        ui_units=pyunits.mm,
        display_units="mm",
        rounding=rounding,
        description="Particle diameter of the GAC media",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    # ---------------------------------------------------------------------
    # adsorption parameters
    category = "Adsorption parameters"
    exports.add(
        obj=fs.gac.freund_k,
        name="Freundlich isotherm parameter k",
        ui_units=pyunits.dimensionless,
        display_units="(m3/kg)**(1/n)",
        rounding=rounding,
        description="Freundlich isotherm parameter k",
        is_input=True,
        input_category=category,
        is_output=False,
    )
    exports.add(
        obj=fs.gac.freund_ninv,
        name="Freundlich isotherm parameter 1/n",
        ui_units=pyunits.dimensionless,
        display_units="",
        rounding=rounding,
        description="Freundlich isotherm parameter 1/n",
        is_input=True,
        input_category=category,
        is_output=False,
    )
    exports.add(
        obj=fs.gac.ds,
        name="Surface diffusion coefficient",
        ui_units=pyunits.m**2 * pyunits.s**-1,
        display_units="m2/s",
        rounding=rounding,
        description="Surface diffusion coefficient",
        is_input=True,
        input_category=category,
        is_output=False,
    )
    exports.add(
        obj=fs.gac.kf,
        name="Film transfer coefficient",
        ui_units=pyunits.m * pyunits.s**-1,
        display_units="m/s",
        rounding=rounding,
        description="Liquid phase film transfer coefficient",
        is_input=True,
        input_category=category,
        is_output=False,
    )
    # ---------------------------------------------------------------------
    # bed design
    category = "Bed design"
    exports.add(
        obj=fs.gac.velocity_sup,
        name="Superficial velocity",
        ui_units=pyunits.m * pyunits.hour**-1,
        display_units="m/h",
        rounding=rounding,
        description="Superficial velocity",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.gac.velocity_int,
        name="Interstitial velocity",
        ui_units=pyunits.m * pyunits.hour**-1,
        display_units="m/h",
        rounding=rounding,
        description="Interstitial velocity",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.gac.bed_length,
        name="Bed length",
        ui_units=pyunits.m,
        display_units="m",
        rounding=rounding,
        description="Length of the GAC bed",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.gac.bed_diameter,
        name="Effective bed diameter",
        ui_units=pyunits.m,
        display_units="m",
        rounding=rounding,
        description="Effective bed diameter as if only once bed was present",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.gac.bed_area,
        name="Effective total bed area",
        ui_units=pyunits.m**2,
        display_units="m2",
        rounding=rounding,
        description="Effective total bed area for flow",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.gac.bed_volume,
        name="Effective total bed volume",
        ui_units=pyunits.m**3,
        display_units="m3",
        rounding=rounding,
        description="Effective total bed volume",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.gac.ebct,
        name="Empty bed contact time",
        ui_units=pyunits.min,
        display_units="min",
        rounding=rounding,
        description="Empty bed contact time",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.gac.residence_time,
        name="Residence time",
        ui_units=pyunits.min,
        display_units="min",
        rounding=rounding,
        description="Residence time",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.gac.bed_voidage,
        name="Bed voidage",
        ui_units=pyunits.dimensionless,
        display_units="",
        rounding=rounding,
        description="Voidage of the GAC bed",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.gac.bed_mass_gac,
        name="Mass of media in the bed",
        ui_units=pyunits.kg,
        display_units="kg",
        rounding=rounding,
        description="Mass of media in the bed",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    # ---------------------------------------------------------------------
    # system performance
    category = "System performance"
    exports.add(
        obj=fs.gac.conc_ratio_replace,
        name="Concentration ratio at breakthrough",
        ui_units=pyunits.dimensionless,
        display_units="",
        rounding=rounding,
        description="Concentration ratio at breakthrough, corresponding to the time of replacing the media",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.gac.conc_ratio_avg,
        name="Average concentration ratio in operational time",
        ui_units=pyunits.dimensionless,
        display_units="",
        rounding=rounding,
        description="Average concentration ratio in operational time",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.gac.operational_time,
        name="Time operating until breakthrough",
        ui_units=pyunits.day,
        display_units="days",
        rounding=rounding,
        description="Time operating until breakthrough",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.gac.bed_volumes_treated,
        name="Bed volumes treated",
        ui_units=pyunits.dimensionless,
        display_units="",
        rounding=rounding,
        description="Bed volumes treated",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.gac.mass_adsorbed,
        name="Total mass of solute adsorbed in operational time",
        ui_units=pyunits.kg,
        display_units="kg",
        rounding=rounding,
        description="Total mass of solute adsorbed in operational time",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.gac.gac_usage_rate,
        name="Mass usage rate of media in operational time",
        ui_units=pyunits.kg * pyunits.day**-1,
        display_units="kg/d",
        rounding=rounding,
        description="Mass usage rate of media in operational time",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    # ---------------------------------------------------------------------
    # empirical parameters
    category = "Empirical parameters"
    exports.add(
        obj=fs.gac.a0,
        name="Stanton equation parameter 0",
        ui_units=pyunits.dimensionless,
        display_units="",
        rounding=rounding,
        description="Stanton equation parameter 0",
        is_input=True,
        input_category=category,
        is_output=False,
    )
    exports.add(
        obj=fs.gac.a1,
        name="Stanton equation parameter 1",
        ui_units=pyunits.dimensionless,
        display_units="",
        rounding=rounding,
        description="Stanton equation parameter 1",
        is_input=True,
        input_category=category,
        is_output=False,
    )
    exports.add(
        obj=fs.gac.b0,
        name="Throughput equation parameter 0",
        ui_units=pyunits.dimensionless,
        display_units="",
        rounding=rounding,
        description="Throughput equation parameter 0",
        is_input=True,
        input_category=category,
        is_output=False,
    )
    exports.add(
        obj=fs.gac.b1,
        name="Throughput equation parameter 1",
        ui_units=pyunits.dimensionless,
        display_units="",
        rounding=rounding,
        description="Throughput equation parameter 1",
        is_input=True,
        input_category=category,
        is_output=False,
    )
    exports.add(
        obj=fs.gac.b2,
        name="Throughput equation parameter 2",
        ui_units=pyunits.dimensionless,
        display_units="",
        rounding=rounding,
        description="Throughput equation parameter 2",
        is_input=True,
        input_category=category,
        is_output=False,
    )
    exports.add(
        obj=fs.gac.b3,
        name="Throughput equation parameter 3",
        ui_units=pyunits.dimensionless,
        display_units="",
        rounding=rounding,
        description="Throughput equation parameter 3",
        is_input=True,
        input_category=category,
        is_output=False,
    )
    exports.add(
        obj=fs.gac.b4,
        name="Throughput equation parameter 4",
        ui_units=pyunits.dimensionless,
        display_units="",
        rounding=rounding,
        description="Throughput equation parameter 4",
        is_input=True,
        input_category=category,
        is_output=False,
    )
    # ---------------------------------------------------------------------
    # costing options
    category = "Costing options"
    exports.add(
        obj=fs.costing.gac_pressure.regen_frac,
        name="Fraction of media regenerated after breakthrough",
        ui_units=pyunits.dimensionless,
        display_units="",
        rounding=rounding,
        description="Fraction of media regenerated after breakthrough",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.costing.gac_pressure.num_contactors_op,
        name="Number of operating beds",
        ui_units=pyunits.dimensionless,
        display_units="",
        rounding=rounding,
        description="Number of operating beds",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.costing.gac_pressure.num_contactors_redundant,
        name="Number of redundant beds",
        ui_units=pyunits.dimensionless,
        display_units="",
        rounding=rounding,
        description="Number of redundant beds",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
    )


def build_flowsheet(build_options=None, **kwargs):
    """
    Builds the initial flowsheet.
    """

    m = gac_fs.build()
    gac_fs.initialize_model(m)

    return m


def solve_flowsheet(flowsheet=None):
    """
    Solves the initial flowsheet.
    """

    res = gac_fs.solve_model(flowsheet)

    return res
