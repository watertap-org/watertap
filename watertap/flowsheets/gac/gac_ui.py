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
"""
GUI configuration for the GAC model.
"""

from pyomo.environ import units as pyunits
from idaes.core import MaterialFlowBasis
from watertap.core.solvers import get_solver
from watertap.property_models.multicomp_aq_sol_prop_pack import DiffusivityCalculation
from watertap.unit_models.gac import (
    FilmTransferCoefficientType,
    SurfaceDiffusionCoefficientType,
)
from watertap.costing.unit_models.gac import ContactorType
from watertap.flowsheets.gac import gac as gac_fs
from idaes_flowsheet_processor.api import FlowsheetInterface, FlowsheetCategory

__author__ = "Hunter Barber"


def export_to_ui():
    """
    Exports the variables, flowsheet build, and solver results to the GUI.
    """
    return FlowsheetInterface(
        name="Granular Activated Carbon (GAC)",
        do_export=export_variables,
        do_build=build_flowsheet,
        do_solve=solve_flowsheet,
        requires_idaes_solver=True,
        category=FlowsheetCategory.wastewater,
        build_options={
            "MaterialFlowBasis": {
                "name": "MaterialFlowBasis",
                "display_name": "Material Flow Basis",
                "values_allowed": [
                    x.name for x in [MaterialFlowBasis.molar, MaterialFlowBasis.mass]
                ],
                "value": MaterialFlowBasis.molar.name,
            },
            "FilmTransferCoefficientType": {
                "name": "FilmTransferCoefficientType",
                "display_name": "Film Transfer Coefficient Type",
                "values_allowed": [x.name for x in FilmTransferCoefficientType],
                "value": FilmTransferCoefficientType.fixed.name,
            },
            "SurfaceDiffusionCoefficientType": {
                "name": "SurfaceDiffusionCoefficientType",
                "display_name": "Surface Diffusion Coefficient Type",
                "values_allowed": [x.name for x in SurfaceDiffusionCoefficientType],
                "value": SurfaceDiffusionCoefficientType.fixed.name,
            },
            "DiffusivityCalculation": {
                "name": "DiffusivityCalculation",
                "display_name": "Diffusivity Calculation",
                "values_allowed": [x.name for x in DiffusivityCalculation],
                "value": DiffusivityCalculation.none.name,
            },
            "ContactorType": {
                "name": "ContactorType",
                "display_name": "Cost Contactor Type",
                "values_allowed": [x.name for x in ContactorType],
                "value": ContactorType.pressure.name,
            },
        },
    )


def export_variables(flowsheet=None, exports=None, build_options=None, **kwargs):
    """
    Exports the variables to the GUI.
    """
    fs = flowsheet
    rounding = 3
    sci_not_rounding = 20
    solute_name = "solute"
    fs.costing.display()
    # input data
    # ---------------------------------------------------------------------
    # solute properties
    category = "Solute properties"
    exports.add(
        obj=fs.properties.mw_comp["solute"],
        name="MW",
        ui_units=pyunits.gram / pyunits.mol,
        display_units="g/mol",
        rounding=rounding,
        description="Molecular weight",
        is_input=True,
        input_category=category,
        is_output=False,
    )
    if (
        build_options["FilmTransferCoefficientType"].value == "calculated"
        or build_options["SurfaceDiffusionCoefficientType"].value == "calculated"
    ):
        if build_options["DiffusivityCalculation"].value == "none":
            exports.add(
                obj=fs.properties.diffus_phase_comp["Liq", "solute"],
                name="Diffusivity",
                ui_units=pyunits.m**2 / pyunits.s,
                display_units="m2/s",
                rounding=sci_not_rounding,
                description="Diffusivity",
                is_input=True,
                input_category=category,
                is_output=False,
            )
        else:
            exports.add(
                obj=fs.properties.molar_volume_phase_comp["Liq", "solute"],
                name="Molar volume",
                ui_units=pyunits.m**3 / pyunits.mol,
                display_units="m3/mol",
                rounding=sci_not_rounding,
                description="Molar volume",
                is_input=True,
                input_category=category,
                is_output=False,
            )
    # ---------------------------------------------------------------------
    # feed conditions
    category = "Feed"
    exports.add(
        obj=fs.feed.properties[0].temperature,
        name="Temperature",
        ui_units=pyunits.kelvin,
        display_units="K",
        rounding=rounding,
        description="Feed temperature",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.feed.properties[0].pressure,
        name="Pressure",
        ui_units=pyunits.bar,
        display_units="bar",
        rounding=rounding,
        description="Feed pressure",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.feed.properties[0].flow_vol_phase["Liq"],
        name="Volumetric flow rate",
        ui_units=pyunits.Mgallon / pyunits.day,
        display_units="MGD",
        rounding=rounding,
        description="Feed volumetric flow rate",
        is_input=False,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.feed.properties[0].conc_mass_phase_comp["Liq", solute_name],
        name="Solute concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=rounding,
        description="Feed solute concentration",
        is_input=False,
        is_output=True,
        output_category=category,
    )
    if build_options["MaterialFlowBasis"].value == "molar":
        fix_molar = True
        fix_mass = False
    else:
        fix_molar = False
        fix_mass = True
    exports.add(
        obj=fs.feed.properties[0].flow_mol_phase_comp["Liq", "H2O"],
        name="Molar flow rate water",
        ui_units=pyunits.mol / pyunits.s,
        display_units="mol/s",
        rounding=rounding,
        description="Feed molar flow rate of water",
        is_input=fix_molar,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.feed.properties[0].flow_mol_phase_comp["Liq", solute_name],
        name="Molar flow rate solute",
        ui_units=pyunits.mol / pyunits.s,
        display_units="mol/s",
        rounding=rounding,
        description="Feed molar flow rate of solute",
        is_input=fix_molar,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"],
        name="Mass flow rate water",
        ui_units=pyunits.kg / pyunits.s,
        display_units="kg/s",
        rounding=rounding,
        description="Feed mass flow rate of water",
        is_input=fix_mass,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.feed.properties[0].flow_mass_phase_comp["Liq", solute_name],
        name="Mass flow rate solute",
        ui_units=pyunits.kg / pyunits.s,
        display_units="kg/s",
        rounding=rounding,
        description="Feed mass flow rate of solute",
        is_input=fix_mass,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    # ---------------------------------------------------------------------
    # product conditions
    category = "Product"
    exports.add(
        obj=fs.product.properties[0].temperature,
        name="Temperature",
        ui_units=pyunits.kelvin,
        display_units="K",
        rounding=rounding,
        description="Product temperature",
        is_input=False,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.product.properties[0].pressure,
        name="Pressure",
        ui_units=pyunits.bar,
        display_units="bar",
        rounding=rounding,
        description="Product pressure",
        is_input=False,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.product.properties[0].flow_vol_phase["Liq"],
        name="Volumetric flow rate",
        ui_units=pyunits.Mgallon / pyunits.day,
        display_units="MGD",
        rounding=rounding,
        description="Product volumetric flow rate",
        is_input=False,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.product.properties[0].conc_mass_phase_comp["Liq", solute_name],
        name="Solute concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=rounding,
        description="Product solute concentration",
        is_input=False,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.product.properties[0].flow_mol_phase_comp["Liq", "H2O"],
        name="Molar flow rate water",
        ui_units=pyunits.mol / pyunits.s,
        display_units="mol/s",
        rounding=rounding,
        description="Product molar flow rate of water",
        is_input=False,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.product.properties[0].flow_mol_phase_comp["Liq", solute_name],
        name="Molar flow rate solute",
        ui_units=pyunits.mol / pyunits.s,
        display_units="mol/s",
        rounding=rounding,
        description="Product molar flow rate of solute",
        is_input=False,
        is_output=True,
        output_category=category,
    )
    # ---------------------------------------------------------------------
    # gac media
    category = "GAC media properties"
    exports.add(
        obj=fs.gac.particle_dens_app,
        name="Apparent density",
        ui_units=pyunits.kg * pyunits.m**-3,
        display_units="kg/m3",
        rounding=rounding,
        description="Apparent density of the GAC media",
        is_input=True,
        input_category=category,
        is_output=False,
    )
    exports.add(
        obj=fs.gac.particle_dens_bulk,
        name="Bulk density",
        ui_units=pyunits.kg * pyunits.m**-3,
        display_units="kg/m3",
        rounding=rounding,
        description="Bulk density of the GAC media in the bed",
        is_input=False,
        is_output=False,
    )
    exports.add(
        obj=fs.gac.particle_dia,
        name="Diameter",
        ui_units=pyunits.mm,
        display_units="mm",
        rounding=rounding,
        description="Particle diameter of the GAC media",
        is_input=True,
        input_category=category,
        is_output=False,
    )
    if build_options["SurfaceDiffusionCoefficientType"].value == "calculated":
        exports.add(
            obj=fs.gac.particle_porosity,
            name="Porosity",
            ui_units=pyunits.dimensionless,
            display_units="-",
            rounding=rounding,
            description="Particle porosity of the GAC media",
            is_input=True,
            input_category=category,
            is_output=False,
        )
    # ---------------------------------------------------------------------
    # adsorption parameters
    category = "Adsorption parameters"
    exports.add(
        obj=fs.gac.freund_k,
        name="Freundlich k",
        ui_units=pyunits.dimensionless,
        display_units="(m3/kg)(1/n)",
        rounding=rounding,
        description="Freundlich isotherm parameter k",
        is_input=True,
        input_category=category,
        is_output=False,
    )
    exports.add(
        obj=fs.gac.freund_ninv,
        name="Freundlich 1/n",
        ui_units=pyunits.dimensionless,
        display_units="-",
        rounding=rounding,
        description="Freundlich isotherm parameter 1/n",
        is_input=True,
        input_category=category,
        is_output=False,
    )
    if build_options["SurfaceDiffusionCoefficientType"].value == "fixed":
        ds_is_input = True
    else:
        ds_is_input = False
        exports.add(
            obj=fs.gac.tort,
            name="Tortuosity",
            ui_units=pyunits.dimensionless,
            display_units="-",
            rounding=rounding,
            description="Tortuosity of the path that the adsorbate must take as compared to the radius",
            is_input=True,
            input_category=category,
            is_output=False,
        )
        exports.add(
            obj=fs.gac.spdfr,
            name="SPDFR",
            ui_units=pyunits.dimensionless,
            display_units="-",
            rounding=rounding,
            description="Surface-to-pore diffusion flux ratio",
            is_input=True,
            input_category=category,
            is_output=False,
        )
    exports.add(
        obj=fs.gac.ds,
        name="Surface diffusion coefficient",
        ui_units=pyunits.m**2 * pyunits.s**-1,
        display_units="m2/s",
        rounding=sci_not_rounding,
        description="Surface diffusion coefficient",
        is_input=ds_is_input,
        input_category=category,
        is_output=False,
    )
    if build_options["FilmTransferCoefficientType"].value == "fixed":
        kf_is_input = True
    else:
        kf_is_input = False
        exports.add(
            obj=fs.gac.shape_correction_factor,
            name="SCF",
            ui_units=pyunits.dimensionless,
            display_units="-",
            rounding=rounding,
            description="Film transfer coefficient shape correction factor",
            is_input=True,
            input_category=category,
            is_output=False,
        )
    exports.add(
        obj=fs.gac.kf,
        name="Film transfer coefficient",
        ui_units=pyunits.m * pyunits.s**-1,
        display_units="m/s",
        rounding=sci_not_rounding,
        description="Liquid phase film transfer coefficient",
        is_input=kf_is_input,
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
        name="Bed diameter",
        ui_units=pyunits.m,
        display_units="m",
        rounding=rounding,
        description="Effective bed diameter as if only once bed was present",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
        is_readonly=True,
    )
    exports.add(
        obj=fs.gac.bed_area,
        name="Bed area",
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
        name="Bed volume",
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
        name="EBCT",
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
        display_units="-",
        rounding=rounding,
        description="Voidage of the GAC bed",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.gac.bed_mass_gac,
        name="Media mass",
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
        name="Breakthrough concentration ratio",
        ui_units=pyunits.dimensionless,
        display_units="-",
        rounding=rounding,
        description="Concentration ratio at breakthrough, corresponding to the time of replacing the media",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.gac.conc_ratio_avg,
        name="Average concentration ratio",
        ui_units=pyunits.dimensionless,
        display_units="-",
        rounding=rounding,
        description="Average concentration ratio during operational time",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.gac.operational_time,
        name="Breakthrough time",
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
        name="BVT",
        ui_units=pyunits.dimensionless,
        display_units="-",
        rounding=rounding,
        description="Bed volumes treated",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.gac.mass_adsorbed,
        name="Solute adsorption rate",
        ui_units=pyunits.kg,
        display_units="kg",
        rounding=rounding,
        description="Mass adsorption rate of solute during operational time",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.gac.gac_usage_rate,
        name="Media usage rate",
        ui_units=pyunits.kg * pyunits.day**-1,
        display_units="kg/day",
        rounding=rounding,
        description="Mass usage rate of media during operational time",
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
        name="Stanton parameter 0",
        ui_units=pyunits.dimensionless,
        display_units="-",
        rounding=6,
        description="Stanton parameter 0",
        is_input=True,
        input_category=category,
        is_output=False,
    )
    exports.add(
        obj=fs.gac.a1,
        name="Stanton parameter 1",
        ui_units=pyunits.dimensionless,
        display_units="-",
        rounding=6,
        description="Stanton parameter 1",
        is_input=True,
        input_category=category,
        is_output=False,
    )
    exports.add(
        obj=fs.gac.b0,
        name="Throughput parameter 0",
        ui_units=pyunits.dimensionless,
        display_units="-",
        rounding=6,
        description="Throughput parameter 0",
        is_input=True,
        input_category=category,
        is_output=False,
    )
    exports.add(
        obj=fs.gac.b1,
        name="Throughput parameter 1",
        ui_units=pyunits.dimensionless,
        display_units="-",
        rounding=6,
        description="Throughput parameter 1",
        is_input=True,
        input_category=category,
        is_output=False,
    )
    exports.add(
        obj=fs.gac.b2,
        name="Throughput parameter 2",
        ui_units=pyunits.dimensionless,
        display_units="-",
        rounding=6,
        description="Throughput parameter 2",
        is_input=True,
        input_category=category,
        is_output=False,
    )
    exports.add(
        obj=fs.gac.b3,
        name="Throughput parameter 3",
        ui_units=pyunits.dimensionless,
        display_units="-",
        rounding=6,
        description="Throughput parameter 3",
        is_input=True,
        input_category=category,
        is_output=False,
    )
    exports.add(
        obj=fs.gac.b4,
        name="Throughput parameter 4",
        ui_units=pyunits.dimensionless,
        display_units="-",
        rounding=6,
        description="Throughput parameter 4",
        is_input=True,
        input_category=category,
        is_output=False,
    )
    # ---------------------------------------------------------------------
    # costing options
    category = "Costing options"
    exports.add(
        obj=fs.costing.electricity_cost,
        name="Electricity cost",
        ui_units=fs.costing.base_currency / pyunits.kWh,
        display_units="$/kW",
        rounding=rounding,
        description="Electricity cost",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    if build_options["ContactorType"].value == "pressure":
        unit_model_costing = fs.costing.gac_pressure
    else:
        unit_model_costing = fs.costing.gac_gravity
    exports.add(
        obj=unit_model_costing.regen_frac,
        name="Regeneration fraction",
        ui_units=pyunits.dimensionless,
        display_units="-",
        rounding=rounding,
        description="Fraction of media that can be regenerated after breakthrough (remainder is replaced)",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=unit_model_costing.num_contactors_op,
        name="Number of operating beds",
        ui_units=pyunits.dimensionless,
        display_units="-",
        rounding=rounding,
        description="Number of beds in operation at any steady state time",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=unit_model_costing.num_contactors_redundant,
        name="Number of redundant beds",
        ui_units=pyunits.dimensionless,
        display_units="-",
        rounding=rounding,
        description="Number of beds offline (redundant) at any steady state time",
        is_input=True,
        input_category=category,
        is_output=True,
        output_category=category,
    )
    # ---------------------------------------------------------------------
    # cost results
    category = "Cost Results"
    exports.add(
        obj=fs.costing.aggregate_capital_cost,
        name="Capital costs",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="Capital costs",
        is_input=False,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.costing.aggregate_fixed_operating_cost,
        name="Fixed operating costs",
        ui_units=fs.costing.base_currency / pyunits.year,
        display_units="$/yr",
        rounding=0,
        description="Fixed operating costs",
        is_input=False,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.costing.aggregate_flow_costs["electricity"],
        name="Electricity costs",
        ui_units=fs.costing.base_currency / pyunits.year,
        display_units="$/yr",
        rounding=0,
        description="Electricity costs",
        is_input=False,
        is_output=True,
        output_category=category,
    )
    exports.add(
        obj=fs.costing.LCOW,
        name="LCOW",
        ui_units=fs.costing.base_currency / pyunits.meter**3,
        display_units="$/yr",
        rounding=rounding,
        description="Levelized cost of water",
        is_input=False,
        is_output=True,
        output_category=category,
    )


def build_flowsheet(build_options=None, **kwargs):
    """
    Build and solve the initial flowsheet.
    """

    if build_options is not None:
        m = gac_fs.build(
            material_flow_basis=build_options["MaterialFlowBasis"].value,
            film_transfer_coefficient_type=build_options[
                "FilmTransferCoefficientType"
            ].value,
            surface_diffusion_coefficient_type=build_options[
                "SurfaceDiffusionCoefficientType"
            ].value,
            diffusivity_calculation=build_options["DiffusivityCalculation"].value,
            cost_contactor_type=build_options["ContactorType"].value,
        )
    else:
        m = gac_fs.build()
    gac_fs.initialize(m)

    return m


def solve_flowsheet(flowsheet=None):
    """
    Solves the flowsheet.
    """

    fs = flowsheet
    solver = get_solver()
    res = gac_fs.optimize(fs, solver)

    return res
