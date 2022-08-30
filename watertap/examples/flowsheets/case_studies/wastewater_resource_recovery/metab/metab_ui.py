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
from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.metab.metab import (
    main,
    solve,
)
from idaes.core.solvers import get_solver
from pyomo.environ import units as pyunits


def export_to_ui():
    return FlowsheetInterface(
        name="Metab",
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
        ui_units=pyunits.L / pyunits.hr,
        display_units="L/h",  # can this be done by default?
        rounding=0,
        description="Inlet volumetric flow rate",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.conc_mass_comp[0, "cod"],
        name="COD concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",  # can this be done by default?
        rounding=2,
        description="Inlet chemical oxygen demand (COD) concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )

    # Unit model data, hydrogen reactor
    exports.add(
        obj=fs.metab_hydrogen.recovery_frac_mass_H2O[0],
        name="Water recovery",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Water recovery [g-H2O treated/g-H2O inlet]",
        is_input=True,
        input_category="Hydrogen reactor",
        is_output=False,
    )
    exports.add(
        obj=fs.metab_hydrogen.reaction_conversion[0, "cod_to_hydrogen"],
        name="COD conversion",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="COD conversion [g-COD reacted/g-COD inlet]",
        is_input=True,
        input_category="Hydrogen reactor",
        is_output=False,
    )
    exports.add(
        obj=fs.metab_hydrogen.generation_ratio["cod_to_hydrogen", "hydrogen"],
        name="H2 conversion ratio",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="H2 mass conversion ratio with respect to COD [g-H2 produced/g-COD reacted]",
        is_input=True,
        input_category="Hydrogen reactor",
        is_output=False,
    )
    exports.add(
        obj=fs.metab_hydrogen.generation_ratio["cod_to_hydrogen", "hydrogen"],
        name="H2 conversion ratio",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="H2 mass conversion ratio with respect to COD [g-H2 produced/g-COD reacted]",
        is_input=True,
        input_category="Hydrogen reactor",
        is_output=False,
    )
    exports.add(
        obj=fs.metab_hydrogen.hydraulic_retention_time,
        name="HRT",
        ui_units=pyunits.hr,
        display_units="h",
        rounding=1,
        description="Hydraulic retention time",
        is_input=True,
        input_category="Hydrogen reactor",
        is_output=False,
    )
    exports.add(
        obj=fs.metab_hydrogen.energy_electric_mixer_vol,
        name="Mixer specific power",
        ui_units=pyunits.kW / pyunits.m**3,
        display_units="kW/m3 of reactor",
        rounding=3,
        description="Mixer specific power relating the power to the volume of the reactor",
        is_input=True,
        input_category="Hydrogen reactor",
        is_output=False,
    )
    exports.add(
        obj=fs.metab_hydrogen.energy_electric_vacuum_flow_vol_byproduct,
        name="Vacuum specific power",
        ui_units=pyunits.kW / (pyunits.kg / pyunits.hr),
        display_units="kW/(kg-H2/h)",
        rounding=0,
        description="Vacuum specific power relating the power to the production of H2",
        is_input=True,
        input_category="Hydrogen reactor",
        is_output=False,
    )
    exports.add(
        obj=fs.metab_hydrogen.energy_thermal_flow_vol_inlet,
        name="Specific heating",
        ui_units=pyunits.MJ / pyunits.m**3,
        display_units="MJ/m3 of reactor",
        rounding=0,
        description="Specific heating relating the thermal energy input to reactor volume",
        is_input=True,
        input_category="Hydrogen reactor",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.metab.bead_bulk_density["hydrogen"],
        name="Bead bulk density",
        ui_units=pyunits.kg / pyunits.m**3,
        display_units="kg/m3 of reactor",
        rounding=1,
        description="Bead bulk density [kg-beads/m3 of reactor]",
        is_input=True,
        input_category="Hydrogen reactor",
        is_output=False,
    )

    # Unit cost data, hydrogen reactor
    exports.add(
        obj=fs.costing.metab.reactor_cost["hydrogen"],
        name="Reactor cost",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of reactor",
        rounding=0,
        description="Reactor capital cost parameter",
        is_input=True,
        input_category="Hydrogen reactor costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.metab.mixer_cost["hydrogen"],
        name="Mixer cost",
        ui_units=fs.costing.base_currency / pyunits.kW,
        display_units="$/kW of mixer",
        rounding=0,
        description="Mixer capital cost parameter",
        is_input=True,
        input_category="Hydrogen reactor costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.metab.bead_cost["hydrogen"],
        name="Bead cost",
        ui_units=fs.costing.base_currency / pyunits.kg,
        display_units="$/kg",
        rounding=0,
        description="Bead cost parameter",
        is_input=True,
        input_category="Hydrogen reactor costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.metab.bead_replacement_factor["hydrogen"],
        name="Bead replacement factor",
        ui_units=1 / pyunits.year,
        display_units="1/year",
        rounding=2,
        description="Bead replacement factor - amount of initial beads replaced per year",
        is_input=True,
        input_category="Hydrogen reactor costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.metab.membrane_sidestream_fraction["hydrogen"],
        name="Membrane side stream",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Fraction of reactor volumetric flow that recirculates in membrane side stream",
        is_input=True,
        input_category="Hydrogen reactor costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.metab.membrane_specific_size["hydrogen"],
        name="Membrane specific size",
        ui_units=pyunits.m**2 / (pyunits.m**3 / pyunits.hr),
        display_units="m2/(m3/h of side stream)",
        rounding=2,
        description="Membrane specific size relating membrane area to side stream flow rate",
        is_input=True,
        input_category="Hydrogen reactor costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.metab.membrane_cost["hydrogen"],
        name="Membrane cost",
        ui_units=fs.costing.base_currency / pyunits.m**2,
        display_units="$/m2",
        rounding=0,
        description="Membrane cost parameter",
        is_input=True,
        input_category="Hydrogen reactor costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.metab.vacuum_cost["hydrogen"],
        name="Vacuum cost",
        ui_units=fs.costing.base_currency / (pyunits.kg / pyunits.hr),
        display_units="$/(kg-H2/h)",
        rounding=0,
        description="Vacuum cost parameter",
        is_input=True,
        input_category="Hydrogen reactor costing",
        is_output=False,
    )

    # Unit model data, methane reactor
    exports.add(
        obj=fs.metab_methane.recovery_frac_mass_H2O[0],
        name="Water recovery",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Water recovery [g-H2O treated/g-H2O inlet]",
        is_input=True,
        input_category="Methane reactor",
        is_output=False,
    )
    exports.add(
        obj=fs.metab_methane.reaction_conversion[0, "cod_to_methane"],
        name="COD conversion",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="COD conversion [g-COD reacted/g-COD inlet]",
        is_input=True,
        input_category="Methane reactor",
        is_output=False,
    )
    exports.add(
        obj=fs.metab_methane.generation_ratio["cod_to_methane", "methane"],
        name="CH4 conversion ratio",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="CH4 mass conversion ratio with respect to COD [g-CH4 produced/g-COD reacted]",
        is_input=True,
        input_category="Methane reactor",
        is_output=False,
    )
    exports.add(
        obj=fs.metab_methane.hydraulic_retention_time,
        name="HRT",
        ui_units=pyunits.hr,
        display_units="h",
        rounding=1,
        description="Hydraulic retention time",
        is_input=True,
        input_category="Methane reactor",
        is_output=False,
    )
    exports.add(
        obj=fs.metab_methane.energy_electric_mixer_vol,
        name="Mixer specific power",
        ui_units=pyunits.kW / pyunits.m**3,
        display_units="kW/m3 of reactor",
        rounding=3,
        description="Mixer specific power relating the power to the volume of the reactor",
        is_input=True,
        input_category="Methane reactor",
        is_output=False,
    )
    exports.add(
        obj=fs.metab_methane.energy_electric_vacuum_flow_vol_byproduct,
        name="Vacuum specific power",
        ui_units=pyunits.kW / (pyunits.kg / pyunits.hr),
        display_units="kW/(kg-H2/h)",
        rounding=0,
        description="Vacuum specific power relating the power to the production of H2",
        is_input=True,
        input_category="Methane reactor",
        is_output=False,
    )
    exports.add(
        obj=fs.metab_methane.energy_thermal_flow_vol_inlet,
        name="Specific heating",
        ui_units=pyunits.MJ / pyunits.m**3,
        display_units="MJ/m3 of reactor",
        rounding=0,
        description="Specific heating relating the thermal energy input to reactor volume",
        is_input=True,
        input_category="Methane reactor",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.metab.bead_bulk_density["methane"],
        name="Bead bulk density",
        ui_units=pyunits.kg / pyunits.m**3,
        display_units="kg/m3 of reactor",
        rounding=1,
        description="Bead bulk density [kg-beads/m3 of reactor]",
        is_input=True,
        input_category="Methane reactor",
        is_output=False,
    )

    # Unit cost data, methane reactor
    exports.add(
        obj=fs.costing.metab.reactor_cost["methane"],
        name="Reactor cost",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of reactor",
        rounding=0,
        description="Reactor capital cost parameter",
        is_input=True,
        input_category="Methane reactor costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.metab.mixer_cost["methane"],
        name="Mixer cost",
        ui_units=fs.costing.base_currency / pyunits.kW,
        display_units="$/kW of mixer",
        rounding=0,
        description="Mixer capital cost parameter",
        is_input=True,
        input_category="Methane reactor costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.metab.bead_cost["methane"],
        name="Bead cost",
        ui_units=fs.costing.base_currency / pyunits.kg,
        display_units="$/kg",
        rounding=0,
        description="Bead cost parameter",
        is_input=True,
        input_category="Methane reactor costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.metab.bead_replacement_factor["methane"],
        name="Bead replacement factor",
        ui_units=1 / pyunits.year,
        display_units="1/year",
        rounding=2,
        description="Bead replacement factor - amount of initial beads replaced per year",
        is_input=True,
        input_category="Methane reactor costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.metab.membrane_sidestream_fraction["methane"],
        name="Membrane side stream",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Fraction of reactor volumetric flow that recirculates in membrane side stream",
        is_input=True,
        input_category="Methane reactor costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.metab.membrane_specific_size["methane"],
        name="Membrane specific size",
        ui_units=pyunits.m**2 / (pyunits.m**3 / pyunits.hr),
        display_units="m2/(m3/h of side stream)",
        rounding=2,
        description="Membrane specific size relating membrane area to side stream flow rate",
        is_input=True,
        input_category="Methane reactor costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.metab.membrane_cost["methane"],
        name="Membrane cost",
        ui_units=fs.costing.base_currency / pyunits.m**2,
        display_units="$/m2",
        rounding=0,
        description="Membrane cost parameter",
        is_input=True,
        input_category="Methane reactor costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.metab.vacuum_cost["methane"],
        name="Vacuum cost",
        ui_units=fs.costing.base_currency / (pyunits.kg / pyunits.hr),
        display_units="$/(kg-CH4/h)",
        rounding=0,
        description="Vacuum cost parameter",
        is_input=True,
        input_category="Methane reactor costing",
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
        is_input=False,
        is_output=True,
        output_category="System costing",
    )
    exports.add(
        obj=fs.costing.TIC,
        name="Practical investment factor",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=1,
        description="Practical investment factor - [total investment cost/direct capital costs]",
        is_input=False,
        is_output=True,
        output_category="System costing",
    )
    exports.add(
        obj=fs.costing.plant_lifetime,
        name="Plant lifetime",
        ui_units=pyunits.year,
        display_units="years",
        rounding=1,
        description="Plant lifetime",
        is_input=False,
        is_output=True,
        output_category="System costing",
    )
    exports.add(
        obj=fs.costing.wacc,
        name="Discount rate",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=1,
        description="Discount rate used in calculating the capital annualization",
        is_input=False,
        is_output=True,
        output_category="System costing",
    )
    exports.add(
        obj=fs.costing.maintenance_costs_percent_FCI,
        name="Fixed operating cost factor",
        ui_units=1 / pyunits.year,
        display_units="fraction/year",
        rounding=1,
        description="Fixed operating cost factor - [annual fixed operating cost/total investment cost]",
        is_input=False,
        is_output=True,
        output_category="System costing",
    )
    exports.add(
        obj=fs.costing.electricity_cost,
        name="Electricity cost",
        ui_units=fs.costing.base_currency / pyunits.kWh,
        display_units="$/kWh",
        rounding=3,
        description="Electricity cost",
        is_input=False,
        is_output=True,
        output_category="System costing",
    )
    exports.add(
        obj=fs.costing.heat_cost,
        name="Heating cost",
        ui_units=fs.costing.base_currency / pyunits.MJ,
        display_units="$/MJ",
        rounding=3,
        description="Heating cost",
        is_input=False,
        is_output=True,
        output_category="System costing",
    )
    exports.add(
        obj=fs.costing.hydrogen_product_cost,
        name="Hydrogen cost",
        ui_units=fs.costing.base_currency / pyunits.kg,
        display_units="$/kg",
        rounding=3,
        description="Hydrogen cost is negative because it is sold",
        is_input=False,
        is_output=True,
        output_category="System costing",
    )
    exports.add(
        obj=fs.costing.methane_product_cost,
        name="Methane cost",
        ui_units=fs.costing.base_currency / pyunits.kg,
        display_units="$/kg",
        rounding=3,
        description="Methane cost is negative because it is sold",
        is_input=False,
        is_output=True,
        output_category="System costing",
    )

    # System metrics
    exports.add(
        obj=fs.costing.LCOT,
        name="Levelized cost of treatment",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of feed",
        rounding=1,
        description="Levelized cost of treatment including operating and capital costs",
        is_input=False,
        is_output=True,
        output_category="System metrics",
    )
    exports.add(
        obj=fs.costing.LCOW,
        name="Levelized cost of water",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of product",
        rounding=1,
        description="Levelized cost of water including operating and capital costs",
        is_input=False,
        is_output=True,
        output_category="System metrics",
    )
    exports.add(
        obj=fs.costing.LCOH,
        name="Levelized cost of hydrogen",
        ui_units=fs.costing.base_currency / pyunits.kg,
        display_units="$/kg-H2",
        rounding=1,
        description="Levelized cost of hydrogen including operating and capital costs",
        is_input=False,
        is_output=True,
        output_category="System metrics",
    )
    exports.add(
        obj=fs.costing.LCOM,
        name="Levelized cost of methane",
        ui_units=fs.costing.base_currency / pyunits.kg,
        display_units="$/kg-CH4",
        rounding=1,
        description="Levelized cost of methane including operating and capital costs",
        is_input=False,
        is_output=True,
        output_category="System metrics",
    )
    # TODO: add all metric results


def build_flowsheet():
    # build and solve initial flowsheet
    (m, results) = main()
    return m.fs


def solve_flowsheet(flowsheet=None):
    fs = flowsheet
    results = solve(fs)
    return results
