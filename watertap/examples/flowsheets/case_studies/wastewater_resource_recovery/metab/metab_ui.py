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
    # TODO: add all unit costing data
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
    # TODO: add all unit costing data
    # TODO: add system costing data
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
