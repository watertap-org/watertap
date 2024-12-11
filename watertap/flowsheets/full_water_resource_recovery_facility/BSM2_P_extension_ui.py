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
GUI configuration for the extended BSM2 flowsheet.
"""

from pyomo.environ import units as pyunits

import idaes.logger as idaeslog

from idaes_flowsheet_processor.api import FlowsheetInterface

from watertap.flowsheets.full_water_resource_recovery_facility.BSM2_P_extension import (
    build,
    set_operating_conditions,
    initialize_system,
    solve,
    add_costing,
)

from watertap.core.util.initialization import (
    assert_degrees_of_freedom,
    interval_initializer,
)

# Set up logger
_log = idaeslog.getLogger(__name__)


def export_to_ui():
    """
    Exports the variables, flowsheet build, and solver results to the GUI.
    """
    return FlowsheetInterface(
        name="BSM2_P_extension",
        do_export=export_variables,
        do_build=build_flowsheet,
        do_solve=solve_flowsheet,
        requires_idaes_solver=True,
        build_options={
            "BioP": {
                "name": "BioP",
                "display_name": "Phosphorus Biomass Transformation",
                "values_allowed": ["False", "True"],
                "value": "False",  # default value
            },
        },
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
        obj=fs.FeedWater.conc_mass_comp[0, "S_A"],
        name="Feed S_A concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=2,
        description="Fermentation products (acetate) concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.FeedWater.conc_mass_comp[0, "S_F"],
        name="Feed S_F concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=2,
        description="Fermentable organic substrates concentration",
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
        obj=fs.FeedWater.conc_mass_comp[0, "S_IC"],
        name="Feed S_IC concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=2,
        description="Inlet inorganic carbon concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.FeedWater.conc_mass_comp[0, "S_K"],
        name="Feed S_K concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=2,
        description="Inlet potassium concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.FeedWater.conc_mass_comp[0, "S_Mg"],
        name="Feed S_Mg concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=2,
        description="Inlet magnesium concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.FeedWater.conc_mass_comp[0, "S_N2"],
        name="Feed S_N2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=2,
        description="Inlet dinitrogen product (from denitrification) concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.FeedWater.conc_mass_comp[0, "S_NH4"],
        name="Feed S_NH4 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=2,
        description="Inlet ammonium plus ammonia nitrogen concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.FeedWater.conc_mass_comp[0, "S_NO3"],
        name="Feed S_NO3 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=2,
        description="Inlet nitrate plus nitrite nitrogen concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.FeedWater.conc_mass_comp[0, "S_O2"],
        name="Feed S_O2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=2,
        description="Inlet dissolved oxygen concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.FeedWater.conc_mass_comp[0, "S_PO4"],
        name="Feed S_PO4 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=2,
        description="Inlet inorganic soluble phosphorus concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.FeedWater.conc_mass_comp[0, "X_AUT"],
        name="Feed X_AUT concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=2,
        description="Inlet autotrophic nitrifying biomass concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.FeedWater.conc_mass_comp[0, "X_H"],
        name="Feed X_H concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=2,
        description="Inlet heterotrophic biomass concentration",
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
        description="Inlet inert particulate organic matter concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.FeedWater.conc_mass_comp[0, "X_PAO"],
        name="Feed X_PAO concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=2,
        description="Inlet phosphate-accumulating biomass concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.FeedWater.conc_mass_comp[0, "X_PHA"],
        name="Feed X_PHA concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=2,
        description="Inlet polyhydroxyalkanoates concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.FeedWater.conc_mass_comp[0, "X_PP"],
        name="Feed X_PP concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=2,
        description="Inlet poly-phosphate concentration",
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
        name="Third anoxic reactor volume",
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
        name="Fourth anoxic reactor volume",
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
        obj=fs.R6.volume[0],
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
        obj=fs.R7.volume[0],
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
        obj=fs.AD.volume_liquid[0],
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
        obj=fs.AD.volume_vapor[0],
        name="Vapor volume",
        ui_units=pyunits.m**3,
        display_units="m3",
        rounding=1,
        description="Vapor volume",
        is_input=True,
        input_category="Anaerobic digester",
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
        obj=fs.CL.split_fraction[0, "effluent", "S_A"],
        name="S_A split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Fermentation products (acetate) split fraction",
        is_input=True,
        input_category="Primary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL.split_fraction[0, "effluent", "S_F"],
        name="S_F split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Fermentable organic substrates split fraction",
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
        description="Inlet soluble inert organic matter split fraction",
        is_input=True,
        input_category="Primary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL.split_fraction[0, "effluent", "S_IC"],
        name="S_IC split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Inlet inorganic carbon split fraction",
        is_input=True,
        input_category="Primary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL.split_fraction[0, "effluent", "S_K"],
        name="S_K split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Inlet potassium split fraction",
        is_input=True,
        input_category="Primary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL.split_fraction[0, "effluent", "S_Mg"],
        name="S_Mg split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Inlet magnesium split fraction",
        is_input=True,
        input_category="Primary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL.split_fraction[0, "effluent", "S_N2"],
        name="S_N2 split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Inlet dinitrogen product (from denitrification) split fraction",
        is_input=True,
        input_category="Primary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL.split_fraction[0, "effluent", "S_NH4"],
        name="S_NH4 split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Inlet ammonium plus ammonia nitrogen split fraction",
        is_input=True,
        input_category="Primary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL.split_fraction[0, "effluent", "S_NO3"],
        name="S_NO3 split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Inlet nitrate plus nitrite nitrogen split fraction",
        is_input=True,
        input_category="Primary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL.split_fraction[0, "effluent", "S_O2"],
        name="S_O2 split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Inlet dissolved oxygen split fraction",
        is_input=True,
        input_category="Primary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL.split_fraction[0, "effluent", "S_PO4"],
        name="S_PO4 split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Inlet inorganic soluble phosphorus split fraction",
        is_input=True,
        input_category="Primary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL.split_fraction[0, "effluent", "X_AUT"],
        name="X_AUT split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Inlet autotrophic nitrifying biomass split fraction",
        is_input=True,
        input_category="Primary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL.split_fraction[0, "effluent", "X_H"],
        name="X_H split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Inlet heterotrophic biomass split fraction",
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
        description="Inlet inert particulate organic matter split fraction",
        is_input=True,
        input_category="Primary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL.split_fraction[0, "effluent", "X_PAO"],
        name="X_PAO split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Inlet phosphate-accumulating biomass split fraction",
        is_input=True,
        input_category="Primary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL.split_fraction[0, "effluent", "X_PHA"],
        name="X_PHA split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Inlet polyhydroxyalkanoates split fraction",
        is_input=True,
        input_category="Primary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL.split_fraction[0, "effluent", "X_PP"],
        name="X_PP split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Inlet poly-phosphate split fraction",
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
        description="Inlet slowly biodegradable substrate split fraction",
        is_input=True,
        input_category="Primary clarifier",
        is_output=False,
    )

    # Unit model data, secondary clarifier
    exports.add(
        obj=fs.CL2.split_fraction[0, "effluent", "H2O"],
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
        obj=fs.CL2.split_fraction[0, "effluent", "S_A"],
        name="S_A split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Fermentation products (acetate) split fraction",
        is_input=True,
        input_category="Secondary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL2.split_fraction[0, "effluent", "S_F"],
        name="S_F split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Fermentable organic substrates split fraction",
        is_input=True,
        input_category="Secondary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL2.split_fraction[0, "effluent", "S_I"],
        name="S_I split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Inlet soluble inert organic matter split fraction",
        is_input=True,
        input_category="Secondary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL2.split_fraction[0, "effluent", "S_IC"],
        name="S_IC split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Inlet inorganic carbon split fraction",
        is_input=True,
        input_category="Secondary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL2.split_fraction[0, "effluent", "S_K"],
        name="S_K split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Inlet potassium split fraction",
        is_input=True,
        input_category="Secondary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL2.split_fraction[0, "effluent", "S_Mg"],
        name="S_Mg split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Inlet magnesium split fraction",
        is_input=True,
        input_category="Secondary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL2.split_fraction[0, "effluent", "S_N2"],
        name="S_N2 split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Inlet dinitrogen product (from denitrification) split fraction",
        is_input=True,
        input_category="Secondary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL2.split_fraction[0, "effluent", "S_NH4"],
        name="S_NH4 split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Inlet ammonium plus ammonia nitrogen split fraction",
        is_input=True,
        input_category="Secondary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL2.split_fraction[0, "effluent", "S_NO3"],
        name="S_NO3 split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Inlet nitrate plus nitrite nitrogen split fraction",
        is_input=True,
        input_category="Secondary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL2.split_fraction[0, "effluent", "S_O2"],
        name="S_O2 split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Inlet dissolved oxygen split fraction",
        is_input=True,
        input_category="Secondary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL2.split_fraction[0, "effluent", "S_PO4"],
        name="S_PO4 split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Inlet inorganic soluble phosphorus split fraction",
        is_input=True,
        input_category="Secondary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL2.split_fraction[0, "effluent", "X_AUT"],
        name="X_AUT split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Inlet autotrophic nitrifying biomass split fraction",
        is_input=True,
        input_category="Secondary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL2.split_fraction[0, "effluent", "X_H"],
        name="X_H split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Inlet heterotrophic biomass split fraction",
        is_input=True,
        input_category="Secondary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL2.split_fraction[0, "effluent", "X_I"],
        name="X_I split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Inlet inert particulate organic matter split fraction",
        is_input=True,
        input_category="Secondary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL2.split_fraction[0, "effluent", "X_PAO"],
        name="X_PAO split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Inlet phosphate-accumulating biomass split fraction",
        is_input=True,
        input_category="Secondary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL2.split_fraction[0, "effluent", "X_PHA"],
        name="X_PHA split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Inlet polyhydroxyalkanoates split fraction",
        is_input=True,
        input_category="Secondary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL2.split_fraction[0, "effluent", "X_PP"],
        name="X_PP split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Inlet poly-phosphate split fraction",
        is_input=True,
        input_category="Secondary clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.CL2.split_fraction[0, "effluent", "X_S"],
        name="X_S split fraction",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=5,
        description="Inlet slowly biodegradable substrate split fraction",
        is_input=True,
        input_category="Secondary clarifier",
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
        name="Practical investment factor",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=1,
        description="Practical investment factor - [total investment cost/direct "
        "capital costs]",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.plant_lifetime,
        name="Plant lifetime",
        ui_units=pyunits.year,
        display_units="years",
        rounding=1,
        description="Plant lifetime",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.wacc,
        name="Discount rate",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Discount rate used in calculating the capital annualization",
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

    # Cost metrics
    exports.add(
        obj=fs.costing.LCOW,
        name="Levelized cost of water",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3",
        rounding=3,
        description="Levelized cost of water with respect to product water",
        is_input=False,
        is_output=True,
        output_category="Cost metrics",
    )
    exports.add(
        obj=fs.costing.total_operating_cost,
        name="Total operating cost",
        ui_units=fs.costing.base_currency / pyunits.yr,
        display_units="$/yr",
        rounding=3,
        description="Total operating cost",
        is_input=False,
        is_output=True,
        output_category="Cost metrics",
    )
    exports.add(
        obj=fs.costing.total_capital_cost,
        name="Total capital cost",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=3,
        description="Total capital cost",
        is_input=False,
        is_output=True,
        output_category="Cost metrics",
    )
    exports.add(
        obj=fs.costing.total_annualized_cost,
        name="Total annualized cost",
        ui_units=fs.costing.base_currency / pyunits.yr,
        display_units="$/yr",
        rounding=3,
        description="Total annualized cost",
        is_input=False,
        is_output=True,
        output_category="Cost metrics",
    )
    exports.add(
        obj=fs.costing.specific_energy_consumption,
        name="Specific energy consumption",
        ui_units=pyunits.kWh / pyunits.m**3,
        display_units="kWh/m3",
        rounding=3,
        description="Specific energy consumption with respect to influent flowrate",
        is_input=False,
        is_output=True,
        output_category="Cost metrics",
    )

    # Capital costs
    exports.add(
        obj=fs.R1.costing.capital_cost,
        name="Reactor 1 capital cost",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=3,
        description="Capital cost of first reactor in activated sludge process",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )
    exports.add(
        obj=fs.R2.costing.capital_cost,
        name="Reactor 2 capital cost",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=3,
        description="Capital cost of second reactor in activated sludge process",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )
    exports.add(
        obj=fs.R3.costing.capital_cost,
        name="Reactor 3 capital cost",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=3,
        description="Capital cost of third reactor in activated sludge process",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )
    exports.add(
        obj=fs.R4.costing.capital_cost,
        name="Reactor 4 capital cost",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=3,
        description="Capital cost of fourth reactor in activated sludge process",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )
    exports.add(
        obj=fs.R5.costing.capital_cost,
        name="Reactor 5 capital cost",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=3,
        description="Capital cost of fifth reactor in activated sludge process",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )
    exports.add(
        obj=fs.R6.costing.capital_cost,
        name="Reactor 6 capital cost",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=3,
        description="Capital cost of sixth reactor in activated sludge process",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )
    exports.add(
        obj=fs.R7.costing.capital_cost,
        name="Reactor 7 capital cost",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=3,
        description="Capital cost of seventh reactor in activated sludge process",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )
    exports.add(
        obj=fs.CL.costing.capital_cost,
        name="Primary clarifier capital cost",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=3,
        description="Capital cost of primary clarifier",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )
    # exports.add(
    #     obj=fs.CL2.costing.capital_cost,
    #     name="Secondary clarifier capital cost",
    #     ui_units=fs.costing.base_currency,
    #     display_units="$",
    #     rounding=3,
    #     description="Capital cost of secondary clarifier",
    #     is_input=False,
    #     is_output=True,
    #     output_category="Capital costs",
    # )
    # exports.add(
    #     obj=fs.AD.costing.capital_cost,
    #     name="Anaerobic digester capital cost",
    #     ui_units=fs.costing.base_currency,
    #     display_units="$",
    #     rounding=3,
    #     description="Capital cost of anaerobic digester",
    #     is_input=False,
    #     is_output=True,
    #     output_category="Capital costs",
    # )
    # exports.add(
    #     obj=fs.dewater.costing.capital_cost,
    #     name="Dewatering unit capital cost",
    #     ui_units=fs.costing.base_currency,
    #     display_units="$",
    #     rounding=3,
    #     description="Capital cost of dewatering",
    #     is_input=False,
    #     is_output=True,
    #     output_category="Capital costs",
    # )
    # exports.add(
    #     obj=fs.thickener.costing.capital_cost,
    #     name="Thickener capital cost",
    #     ui_units=fs.costing.base_currency,
    #     display_units="$",
    #     rounding=3,
    #     description="Capital cost of thickener",
    #     is_input=False,
    #     is_output=True,
    #     output_category="Capital costs",
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
        obj=fs.Treated.properties[0].conc_mass_comp["S_A"],
        name="S_A concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Fermentation products (acetate) concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Effluent",
    )
    exports.add(
        obj=fs.Treated.properties[0].conc_mass_comp["S_F"],
        name="S_F concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Fermentable organic substrates concentration",
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
        obj=fs.Treated.properties[0].conc_mass_comp["S_IC"],
        name="S_IC concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet inorganic carbon concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Effluent",
    )
    exports.add(
        obj=fs.Treated.properties[0].conc_mass_comp["S_K"],
        name="S_K concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet potassium concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Effluent",
    )
    exports.add(
        obj=fs.Treated.properties[0].conc_mass_comp["S_Mg"],
        name="S_Mg concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet magnesium concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Effluent",
    )
    exports.add(
        obj=fs.Treated.properties[0].conc_mass_comp["S_N2"],
        name="S_N2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet dinitrogen product (from denitrification) concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Effluent",
    )
    exports.add(
        obj=fs.Treated.properties[0].conc_mass_comp["S_NH4"],
        name="S_NH4 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet ammonium plus ammonia nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Effluent",
    )
    exports.add(
        obj=fs.Treated.properties[0].conc_mass_comp["S_NO3"],
        name="S_NO3 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet nitrate plus nitrite nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Effluent",
    )
    exports.add(
        obj=fs.Treated.properties[0].conc_mass_comp["S_O2"],
        name="S_O2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet dissolved oxygen concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Effluent",
    )
    exports.add(
        obj=fs.Treated.properties[0].conc_mass_comp["S_PO4"],
        name="S_PO4 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet inorganic soluble phosphorus concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Effluent",
    )
    exports.add(
        obj=fs.Treated.properties[0].conc_mass_comp["X_AUT"],
        name="X_AUT concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet autotrophic nitrifying biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Effluent",
    )
    exports.add(
        obj=fs.Treated.properties[0].conc_mass_comp["X_H"],
        name="X_H concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet heterotrophic biomass concentration",
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
        description="Outlet inert particulate organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Effluent",
    )
    exports.add(
        obj=fs.Treated.properties[0].conc_mass_comp["X_PAO"],
        name="X_PAO concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet phosphate-accumulating biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Effluent",
    )
    exports.add(
        obj=fs.Treated.properties[0].conc_mass_comp["X_PHA"],
        name="X_PHA concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet polyhydroxyalkanoates concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Effluent",
    )
    exports.add(
        obj=fs.Treated.properties[0].conc_mass_comp["X_PP"],
        name="X_PP concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet poly-phosphate concentration",
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

    # Primary clarifier output effluent
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
        obj=fs.CL.effluent.conc_mass_comp[0, "S_A"],
        name="S_A concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Fermentation products (acetate) concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Effluent",
    )
    exports.add(
        obj=fs.CL.effluent.conc_mass_comp[0, "S_F"],
        name="S_F concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Fermentable organic substrates concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Effluent",
    )
    exports.add(
        obj=fs.CL.effluent.conc_mass_comp[0, "S_I"],
        name="S_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet soluble inert organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Effluent",
    )
    exports.add(
        obj=fs.CL.effluent.conc_mass_comp[0, "S_IC"],
        name="S_IC concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet inorganic carbon concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Effluent",
    )
    exports.add(
        obj=fs.CL.effluent.conc_mass_comp[0, "S_K"],
        name="S_K concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet potassium concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Effluent",
    )
    exports.add(
        obj=fs.CL.effluent.conc_mass_comp[0, "S_Mg"],
        name="S_Mg concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet magnesium concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Effluent",
    )
    exports.add(
        obj=fs.CL.effluent.conc_mass_comp[0, "S_N2"],
        name="S_N2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet dinitrogen product (from denitrification) concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Effluent",
    )
    exports.add(
        obj=fs.CL.effluent.conc_mass_comp[0, "S_NH4"],
        name="S_NH4 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet ammonium plus ammonia nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Effluent",
    )
    exports.add(
        obj=fs.CL.effluent.conc_mass_comp[0, "S_NO3"],
        name="S_NO3 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet nitrate plus nitrite nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Effluent",
    )
    exports.add(
        obj=fs.CL.effluent.conc_mass_comp[0, "S_O2"],
        name="S_O2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet dissolved oxygen concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Effluent",
    )
    exports.add(
        obj=fs.CL.effluent.conc_mass_comp[0, "S_PO4"],
        name="S_PO4 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet inorganic soluble phosphorus concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Effluent",
    )
    exports.add(
        obj=fs.CL.effluent.conc_mass_comp[0, "X_AUT"],
        name="X_AUT concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet autotrophic nitrifying biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Effluent",
    )
    exports.add(
        obj=fs.CL.effluent.conc_mass_comp[0, "X_H"],
        name="X_H concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet heterotrophic biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Effluent",
    )
    exports.add(
        obj=fs.CL.effluent.conc_mass_comp[0, "X_I"],
        name="X_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet inert particulate organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Effluent",
    )
    exports.add(
        obj=fs.CL.effluent.conc_mass_comp[0, "X_PAO"],
        name="X_PAO concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet phosphate-accumulating biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Effluent",
    )
    exports.add(
        obj=fs.CL.effluent.conc_mass_comp[0, "X_PHA"],
        name="X_PHA concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet polyhydroxyalkanoates concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Effluent",
    )
    exports.add(
        obj=fs.CL.effluent.conc_mass_comp[0, "X_PP"],
        name="X_PP concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet poly-phosphate concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Effluent",
    )
    exports.add(
        obj=fs.CL.effluent.conc_mass_comp[0, "X_S"],
        name="X_S concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet slowly biodegradable substrate concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Effluent",
    )

    # Primary clarifier output underflow
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
        obj=fs.CL.underflow.conc_mass_comp[0, "S_A"],
        name="S_A concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Fermentation products (acetate) concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Underflow",
    )
    exports.add(
        obj=fs.CL.underflow.conc_mass_comp[0, "S_F"],
        name="S_F concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Fermentable organic substrates concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Underflow",
    )
    exports.add(
        obj=fs.CL.underflow.conc_mass_comp[0, "S_I"],
        name="S_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet soluble inert organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Underflow",
    )
    exports.add(
        obj=fs.CL.underflow.conc_mass_comp[0, "S_IC"],
        name="S_IC concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet inorganic carbon concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Underflow",
    )
    exports.add(
        obj=fs.CL.underflow.conc_mass_comp[0, "S_K"],
        name="S_K concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet potassium concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Underflow",
    )
    exports.add(
        obj=fs.CL.underflow.conc_mass_comp[0, "S_Mg"],
        name="S_Mg concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet magnesium concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Underflow",
    )
    exports.add(
        obj=fs.CL.underflow.conc_mass_comp[0, "S_N2"],
        name="S_N2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet dinitrogen product (from denitrification) concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Underflow",
    )
    exports.add(
        obj=fs.CL.underflow.conc_mass_comp[0, "S_NH4"],
        name="S_NH4 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet ammonium plus ammonia nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Underflow",
    )
    exports.add(
        obj=fs.CL.underflow.conc_mass_comp[0, "S_NO3"],
        name="S_NO3 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet nitrate plus nitrite nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Underflow",
    )
    exports.add(
        obj=fs.CL.underflow.conc_mass_comp[0, "S_O2"],
        name="S_O2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet dissolved oxygen concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Underflow",
    )
    exports.add(
        obj=fs.CL.underflow.conc_mass_comp[0, "S_PO4"],
        name="S_PO4 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet inorganic soluble phosphorus concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Underflow",
    )
    exports.add(
        obj=fs.CL.underflow.conc_mass_comp[0, "X_AUT"],
        name="X_AUT concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet autotrophic nitrifying biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Underflow",
    )
    exports.add(
        obj=fs.CL.underflow.conc_mass_comp[0, "X_H"],
        name="X_H concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet heterotrophic biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Underflow",
    )
    exports.add(
        obj=fs.CL.underflow.conc_mass_comp[0, "X_I"],
        name="X_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet inert particulate organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Underflow",
    )
    exports.add(
        obj=fs.CL.underflow.conc_mass_comp[0, "X_PAO"],
        name="X_PAO concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet phosphate-accumulating biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Underflow",
    )
    exports.add(
        obj=fs.CL.underflow.conc_mass_comp[0, "X_PHA"],
        name="X_PHA concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet polyhydroxyalkanoates concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Underflow",
    )
    exports.add(
        obj=fs.CL.underflow.conc_mass_comp[0, "X_PP"],
        name="X_PP concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet poly-phosphate concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Underflow",
    )
    exports.add(
        obj=fs.CL.underflow.conc_mass_comp[0, "X_S"],
        name="X_S concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Outlet slowly biodegradable substrate concentration",
        is_input=False,
        is_output=True,
        output_category="Primary Clarifier Underflow",
    )

    # Separator 2 SP2
    exports.add(
        obj=fs.SP2.recycle.flow_vol[0],
        name="ASP recycle inlet flow rate",
        ui_units=pyunits.m**3 / pyunits.day,
        display_units="m3/day",
        rounding=2,
        description="ASP recycle flow rate",
        is_input=False,
        is_output=True,
        output_category="ASP Recycle Outlet",
    )
    exports.add(
        obj=fs.SP2.recycle.conc_mass_comp[0, "S_A"],
        name="S_A concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Fermentation products (acetate) concentration",
        is_input=False,
        is_output=True,
        output_category="ASP Recycle Outlet",
    )
    exports.add(
        obj=fs.SP2.recycle.conc_mass_comp[0, "S_F"],
        name="S_F concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Fermentable organic substrates concentration",
        is_input=False,
        is_output=True,
        output_category="ASP Recycle Outlet",
    )
    exports.add(
        obj=fs.SP2.recycle.conc_mass_comp[0, "S_I"],
        name="S_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Soluble inert organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="ASP Recycle Outlet",
    )
    exports.add(
        obj=fs.SP2.recycle.conc_mass_comp[0, "S_IC"],
        name="S_IC concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inorganic carbon concentration",
        is_input=False,
        is_output=True,
        output_category="ASP Recycle Outlet",
    )
    exports.add(
        obj=fs.SP2.recycle.conc_mass_comp[0, "S_K"],
        name="S_K concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Potassium concentration",
        is_input=False,
        is_output=True,
        output_category="ASP Recycle Outlet",
    )
    exports.add(
        obj=fs.SP2.recycle.conc_mass_comp[0, "S_Mg"],
        name="S_Mg concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Magnesium concentration",
        is_input=False,
        is_output=True,
        output_category="ASP Recycle Outlet",
    )
    exports.add(
        obj=fs.SP2.recycle.conc_mass_comp[0, "S_N2"],
        name="S_N2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Dinitrogen product (from denitrification) concentration",
        is_input=False,
        is_output=True,
        output_category="ASP Recycle Outlet",
    )
    exports.add(
        obj=fs.SP2.recycle.conc_mass_comp[0, "S_NH4"],
        name="S_NH4 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Ammonium plus ammonia nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="ASP Recycle Outlet",
    )
    exports.add(
        obj=fs.SP2.recycle.conc_mass_comp[0, "S_NO3"],
        name="S_NO3 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Nitrate plus nitrite nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="ASP Recycle Outlet",
    )
    exports.add(
        obj=fs.SP2.recycle.conc_mass_comp[0, "S_O2"],
        name="S_O2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Dissolved oxygen concentration",
        is_input=False,
        is_output=True,
        output_category="ASP Recycle Outlet",
    )
    exports.add(
        obj=fs.SP2.recycle.conc_mass_comp[0, "S_PO4"],
        name="S_PO4 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inorganic soluble phosphorus concentration",
        is_input=False,
        is_output=True,
        output_category="ASP Recycle Outlet",
    )
    exports.add(
        obj=fs.SP2.recycle.conc_mass_comp[0, "X_AUT"],
        name="X_AUT concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Autotrophic nitrifying biomass concentration",
        is_input=False,
        is_output=True,
        output_category="ASP Recycle Outlet",
    )
    exports.add(
        obj=fs.SP2.recycle.conc_mass_comp[0, "X_H"],
        name="X_H concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Heterotrophic biomass concentration",
        is_input=False,
        is_output=True,
        output_category="ASP Recycle Outlet",
    )
    exports.add(
        obj=fs.SP2.recycle.conc_mass_comp[0, "X_I"],
        name="X_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inert particulate organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="ASP Recycle Outlet",
    )
    exports.add(
        obj=fs.SP2.recycle.conc_mass_comp[0, "X_PAO"],
        name="X_PAO concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Phosphate-accumulating biomass concentration",
        is_input=False,
        is_output=True,
        output_category="ASP Recycle Outlet",
    )
    exports.add(
        obj=fs.SP2.recycle.conc_mass_comp[0, "X_PHA"],
        name="X_PHA concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Polyhydroxyalkanoates concentration",
        is_input=False,
        is_output=True,
        output_category="ASP Recycle Outlet",
    )
    exports.add(
        obj=fs.SP2.recycle.conc_mass_comp[0, "X_PP"],
        name="X_PP concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Poly-phosphate concentration",
        is_input=False,
        is_output=True,
        output_category="ASP Recycle Outlet",
    )
    exports.add(
        obj=fs.SP2.recycle.conc_mass_comp[0, "X_S"],
        name="X_S concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Slowly biodegradable substrate concentration",
        is_input=False,
        is_output=True,
        output_category="ASP Recycle Outlet",
    )

    # Separator SP1
    exports.add(
        obj=fs.SP1.overflow.flow_vol[0],
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
        obj=fs.SP1.overflow.conc_mass_comp[0, "S_A"],
        name="S_A concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Fermentation products (acetate) concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Inlet",
    )
    exports.add(
        obj=fs.SP1.overflow.conc_mass_comp[0, "S_F"],
        name="S_F concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Fermentable organic substrates concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Inlet",
    )
    exports.add(
        obj=fs.SP1.overflow.conc_mass_comp[0, "S_I"],
        name="S_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Soluble inert organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Inlet",
    )
    exports.add(
        obj=fs.SP1.overflow.conc_mass_comp[0, "S_IC"],
        name="S_IC concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inorganic carbon concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Inlet",
    )
    exports.add(
        obj=fs.SP1.overflow.conc_mass_comp[0, "S_K"],
        name="S_K concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Potassium concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Inlet",
    )
    exports.add(
        obj=fs.SP1.overflow.conc_mass_comp[0, "S_Mg"],
        name="S_Mg concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Magnesium concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Inlet",
    )
    exports.add(
        obj=fs.SP1.overflow.conc_mass_comp[0, "S_N2"],
        name="S_N2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Dinitrogen product (from denitrification) concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Inlet",
    )
    exports.add(
        obj=fs.SP1.overflow.conc_mass_comp[0, "S_NH4"],
        name="S_NH4 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Ammonium plus ammonia nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Inlet",
    )
    exports.add(
        obj=fs.SP1.overflow.conc_mass_comp[0, "S_NO3"],
        name="S_NO3 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Nitrate plus nitrite nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Inlet",
    )
    exports.add(
        obj=fs.SP1.overflow.conc_mass_comp[0, "S_O2"],
        name="S_O2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Dissolved oxygen concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Inlet",
    )
    exports.add(
        obj=fs.SP1.overflow.conc_mass_comp[0, "S_PO4"],
        name="S_PO4 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inorganic soluble phosphorus concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Inlet",
    )
    exports.add(
        obj=fs.SP1.overflow.conc_mass_comp[0, "X_AUT"],
        name="X_AUT concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Autotrophic nitrifying biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Inlet",
    )
    exports.add(
        obj=fs.SP1.overflow.conc_mass_comp[0, "X_H"],
        name="X_H concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Heterotrophic biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Inlet",
    )
    exports.add(
        obj=fs.SP1.overflow.conc_mass_comp[0, "X_I"],
        name="X_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inert particulate organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Inlet",
    )
    exports.add(
        obj=fs.SP1.overflow.conc_mass_comp[0, "X_PAO"],
        name="X_PAO concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Phosphate-accumulating biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Inlet",
    )
    exports.add(
        obj=fs.SP1.overflow.conc_mass_comp[0, "X_PHA"],
        name="X_PHA concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Polyhydroxyalkanoates concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Inlet",
    )
    exports.add(
        obj=fs.SP1.overflow.conc_mass_comp[0, "X_PP"],
        name="X_PP concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Poly-phosphate concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Inlet",
    )
    exports.add(
        obj=fs.SP1.overflow.conc_mass_comp[0, "X_S"],
        name="X_S concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Slowly biodegradable substrate concentration",
        is_input=False,
        is_output=True,
        output_category="Secondary Clarifier Inlet",
    )

    # Thickener
    exports.add(
        obj=fs.thickener.inlet.flow_vol[0],
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
        obj=fs.thickener.inlet.conc_mass_comp[0, "S_A"],
        name="S_A concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Fermentation products (acetate) concentration",
        is_input=False,
        is_output=True,
        output_category="Thickener Inlet",
    )
    exports.add(
        obj=fs.thickener.inlet.conc_mass_comp[0, "S_F"],
        name="S_F concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Fermentable organic substrates concentration",
        is_input=False,
        is_output=True,
        output_category="Thickener Inlet",
    )
    exports.add(
        obj=fs.thickener.inlet.conc_mass_comp[0, "S_I"],
        name="S_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Soluble inert organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="Thickener Inlet",
    )
    exports.add(
        obj=fs.thickener.inlet.conc_mass_comp[0, "S_IC"],
        name="S_IC concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inorganic carbon concentration",
        is_input=False,
        is_output=True,
        output_category="Thickener Inlet",
    )
    exports.add(
        obj=fs.thickener.inlet.conc_mass_comp[0, "S_K"],
        name="S_K concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Potassium concentration",
        is_input=False,
        is_output=True,
        output_category="Thickener Inlet",
    )
    exports.add(
        obj=fs.thickener.inlet.conc_mass_comp[0, "S_Mg"],
        name="S_Mg concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Magnesium concentration",
        is_input=False,
        is_output=True,
        output_category="Thickener Inlet",
    )
    exports.add(
        obj=fs.thickener.inlet.conc_mass_comp[0, "S_N2"],
        name="S_N2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Dinitrogen product (from denitrification) concentration",
        is_input=False,
        is_output=True,
        output_category="Thickener Inlet",
    )
    exports.add(
        obj=fs.thickener.inlet.conc_mass_comp[0, "S_NH4"],
        name="S_NH4 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Ammonium plus ammonia nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Thickener Inlet",
    )
    exports.add(
        obj=fs.thickener.inlet.conc_mass_comp[0, "S_NO3"],
        name="S_NO3 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Nitrate plus nitrite nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Thickener Inlet",
    )
    exports.add(
        obj=fs.thickener.inlet.conc_mass_comp[0, "S_O2"],
        name="S_O2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Dissolved oxygen concentration",
        is_input=False,
        is_output=True,
        output_category="Thickener Inlet",
    )
    exports.add(
        obj=fs.thickener.inlet.conc_mass_comp[0, "S_PO4"],
        name="S_PO4 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inorganic soluble phosphorus concentration",
        is_input=False,
        is_output=True,
        output_category="Thickener Inlet",
    )
    exports.add(
        obj=fs.thickener.inlet.conc_mass_comp[0, "X_AUT"],
        name="X_AUT concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Autotrophic nitrifying biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Thickener Inlet",
    )
    exports.add(
        obj=fs.thickener.inlet.conc_mass_comp[0, "X_H"],
        name="X_H concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Heterotrophic biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Thickener Inlet",
    )
    exports.add(
        obj=fs.thickener.inlet.conc_mass_comp[0, "X_I"],
        name="X_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inert particulate organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="Thickener Inlet",
    )
    exports.add(
        obj=fs.thickener.inlet.conc_mass_comp[0, "X_PAO"],
        name="X_PAO concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Phosphate-accumulating biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Thickener Inlet",
    )
    exports.add(
        obj=fs.thickener.inlet.conc_mass_comp[0, "X_PHA"],
        name="X_PHA concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Polyhydroxyalkanoates concentration",
        is_input=False,
        is_output=True,
        output_category="Thickener Inlet",
    )
    exports.add(
        obj=fs.thickener.inlet.conc_mass_comp[0, "X_PP"],
        name="X_PP concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Poly-phosphate concentration",
        is_input=False,
        is_output=True,
        output_category="Thickener Inlet",
    )
    exports.add(
        obj=fs.thickener.inlet.conc_mass_comp[0, "X_S"],
        name="X_S concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Slowly biodegradable substrate concentration",
        is_input=False,
        is_output=True,
        output_category="Thickener Inlet",
    )

    # ASM2d/ADM1 inlet
    exports.add(
        obj=fs.translator_asm2d_adm1.inlet.flow_vol[0],
        name="ASM2D/ADM1 inlet flow rate",
        ui_units=pyunits.m**3 / pyunits.day,
        display_units="m3/day",
        rounding=2,
        description="ASM2D/ADM1 interface inlet flow rate",
        is_input=False,
        is_output=True,
        output_category="ASM2D/ADM1 Interface Inlet",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.inlet.conc_mass_comp[0, "S_A"],
        name="S_A concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Fermentation products (acetate) concentration",
        is_input=False,
        is_output=True,
        output_category="ASM2D/ADM1 Interface Inlet",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.inlet.conc_mass_comp[0, "S_F"],
        name="S_F concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Fermentable organic substrates concentration",
        is_input=False,
        is_output=True,
        output_category="ASM2D/ADM1 Interface Inlet",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.inlet.conc_mass_comp[0, "S_I"],
        name="S_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Soluble inert organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="ASM2D/ADM1 Interface Inlet",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.inlet.conc_mass_comp[0, "S_IC"],
        name="S_IC concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inorganic carbon concentration",
        is_input=False,
        is_output=True,
        output_category="ASM2D/ADM1 Interface Inlet",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.inlet.conc_mass_comp[0, "S_K"],
        name="S_K concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Potassium concentration",
        is_input=False,
        is_output=True,
        output_category="ASM2D/ADM1 Interface Inlet",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.inlet.conc_mass_comp[0, "S_Mg"],
        name="S_Mg concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Magnesium concentration",
        is_input=False,
        is_output=True,
        output_category="ASM2D/ADM1 Interface Inlet",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.inlet.conc_mass_comp[0, "S_N2"],
        name="S_N2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Dinitrogen product (from denitrification) concentration",
        is_input=False,
        is_output=True,
        output_category="ASM2D/ADM1 Interface Inlet",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.inlet.conc_mass_comp[0, "S_NH4"],
        name="S_NH4 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Ammonium plus ammonia nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="ASM2D/ADM1 Interface Inlet",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.inlet.conc_mass_comp[0, "S_NO3"],
        name="S_NO3 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Nitrate plus nitrite nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="ASM2D/ADM1 Interface Inlet",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.inlet.conc_mass_comp[0, "S_O2"],
        name="S_O2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Dissolved oxygen concentration",
        is_input=False,
        is_output=True,
        output_category="ASM2D/ADM1 Interface Inlet",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.inlet.conc_mass_comp[0, "S_PO4"],
        name="S_PO4 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inorganic soluble phosphorus concentration",
        is_input=False,
        is_output=True,
        output_category="ASM2D/ADM1 Interface Inlet",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.inlet.conc_mass_comp[0, "X_AUT"],
        name="X_AUT concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Autotrophic nitrifying biomass concentration",
        is_input=False,
        is_output=True,
        output_category="ASM2D/ADM1 Interface Inlet",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.inlet.conc_mass_comp[0, "X_H"],
        name="X_H concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Heterotrophic biomass concentration",
        is_input=False,
        is_output=True,
        output_category="ASM2D/ADM1 Interface Inlet",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.inlet.conc_mass_comp[0, "X_I"],
        name="X_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inert particulate organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="ASM2D/ADM1 Interface Inlet",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.inlet.conc_mass_comp[0, "X_PAO"],
        name="X_PAO concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Phosphate-accumulating biomass concentration",
        is_input=False,
        is_output=True,
        output_category="ASM2D/ADM1 Interface Inlet",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.inlet.conc_mass_comp[0, "X_PHA"],
        name="X_PHA concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Polyhydroxyalkanoates concentration",
        is_input=False,
        is_output=True,
        output_category="ASM2D/ADM1 Interface Inlet",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.inlet.conc_mass_comp[0, "X_PP"],
        name="X_PP concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Poly-phosphate concentration",
        is_input=False,
        is_output=True,
        output_category="ASM2D/ADM1 Interface Inlet",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.inlet.conc_mass_comp[0, "X_S"],
        name="X_S concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Slowly biodegradable substrate concentration",
        is_input=False,
        is_output=True,
        output_category="ASM2D/ADM1 Interface Inlet",
    )

    # Translator outlet
    exports.add(
        obj=fs.translator_asm2d_adm1.outlet.flow_vol[0],
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
        obj=fs.translator_asm2d_adm1.outlet.conc_mass_comp[0, "S_su"],
        name="S_su concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Monosaccharides concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.outlet.conc_mass_comp[0, "S_aa"],
        name="S_aa concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Amino acids concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.outlet.conc_mass_comp[0, "S_fa"],
        name="S_fa concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Long chain fatty acids concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.outlet.conc_mass_comp[0, "S_va"],
        name="S_va concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Total valerate concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.outlet.conc_mass_comp[0, "S_bu"],
        name="S_bu concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Total butyrate concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.outlet.conc_mass_comp[0, "S_pro"],
        name="S_pro concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Total propionate concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.outlet.conc_mass_comp[0, "S_ac"],
        name="S_ac concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Total acetate concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.outlet.conc_mass_comp[0, "S_IC"],
        name="S_IC concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inorganic carbon concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.outlet.conc_mass_comp[0, "S_IN"],
        name="S_IN concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inorganic nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.outlet.conc_mass_comp[0, "S_IP"],
        name="S_IP concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inorganic phosphorus concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.outlet.conc_mass_comp[0, "S_I"],
        name="S_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Soluble inerts concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.outlet.conc_mass_comp[0, "X_ch"],
        name="X_ch concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Carbohydrates concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.outlet.conc_mass_comp[0, "X_pr"],
        name="X_pr concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Proteins concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.outlet.conc_mass_comp[0, "X_li"],
        name="X_li concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Lipids concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.outlet.conc_mass_comp[0, "X_su"],
        name="X_su concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Sugar degraders concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.outlet.conc_mass_comp[0, "X_aa"],
        name="X_aa concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Amino acid degraders concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.outlet.conc_mass_comp[0, "X_fa"],
        name="X_fa concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Long chain fatty acid (LCFA) degraders concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.outlet.conc_mass_comp[0, "X_c4"],
        name="X_c4 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Valerate and butyrate degraders concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.outlet.conc_mass_comp[0, "X_pro"],
        name="X_pro concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Propionate degraders concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.outlet.conc_mass_comp[0, "X_ac"],
        name="X_ac concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Acetate degraders concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.outlet.conc_mass_comp[0, "X_h2"],
        name="X_h2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Hydrogen degraders concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.outlet.conc_mass_comp[0, "X_I"],
        name="X_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Particulate inerts concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.outlet.conc_mass_comp[0, "X_PHA"],
        name="X_PHA concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Polyhydroxyalkanoates concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.outlet.conc_mass_comp[0, "X_PP"],
        name="X_PP concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Polyphosphates concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.outlet.conc_mass_comp[0, "X_PAO"],
        name="X_PAO concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Phosphorus accumulating organisms concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.outlet.conc_mass_comp[0, "S_K"],
        name="S_K concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Potassium concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_asm2d_adm1.outlet.conc_mass_comp[0, "S_Mg"],
        name="S_Mg concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Magnesium concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Inlet (post-interface)",
    )

    # Translator adm1 asm2d inlet
    exports.add(
        obj=fs.translator_adm1_asm2d.inlet.flow_vol[0],
        name="ADM1/ASM2D inlet flow rate",
        ui_units=pyunits.m**3 / pyunits.day,
        display_units="m3/day",
        rounding=2,
        description="ADM1/ASM2D interface inlet flow rate",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM2D Interface Inlet",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.inlet.conc_mass_comp[0, "S_su"],
        name="S_su concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Monosaccharides concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM2D Interface Inlet",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.inlet.conc_mass_comp[0, "S_aa"],
        name="S_aa concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Amino acids concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM2D Interface Inlet",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.inlet.conc_mass_comp[0, "S_fa"],
        name="S_fa concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Long chain fatty acids concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM2D Interface Inlet",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.inlet.conc_mass_comp[0, "S_va"],
        name="S_va concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Total valerate concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM2D Interface Inlet",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.inlet.conc_mass_comp[0, "S_bu"],
        name="S_bu concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Total butyrate concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM2D Interface Inlet",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.inlet.conc_mass_comp[0, "S_pro"],
        name="S_pro concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Total propionate concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM2D Interface Inlet",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.inlet.conc_mass_comp[0, "S_ac"],
        name="S_ac concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Total acetate concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM2D Interface Inlet",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.inlet.conc_mass_comp[0, "S_IC"],
        name="S_IC concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inorganic carbon concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM2D Interface Inlet",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.inlet.conc_mass_comp[0, "S_IN"],
        name="S_IN concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inorganic nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM2D Interface Inlet",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.inlet.conc_mass_comp[0, "S_IP"],
        name="S_IP concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inorganic phosphorus concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM2D Interface Inlet",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.inlet.conc_mass_comp[0, "S_I"],
        name="S_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Soluble inerts concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM2D Interface Inlet",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.inlet.conc_mass_comp[0, "X_ch"],
        name="X_ch concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Carbohydrates concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM2D Interface Inlet",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.inlet.conc_mass_comp[0, "X_pr"],
        name="X_pr concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Proteins concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM2D Interface Inlet",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.inlet.conc_mass_comp[0, "X_li"],
        name="X_li concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Lipids concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM2D Interface Inlet",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.inlet.conc_mass_comp[0, "X_su"],
        name="X_su concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Sugar degraders concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM2D Interface Inlet",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.inlet.conc_mass_comp[0, "X_aa"],
        name="X_aa concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Amino acid degraders concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM2D Interface Inlet",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.inlet.conc_mass_comp[0, "X_fa"],
        name="X_fa concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Long chain fatty acid (LCFA) degraders concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM2D Interface Inlet",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.inlet.conc_mass_comp[0, "X_c4"],
        name="X_c4 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Valerate and butyrate degraders concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM2D Interface Inlet",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.inlet.conc_mass_comp[0, "X_pro"],
        name="X_pro concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Propionate degraders concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM2D Interface Inlet",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.inlet.conc_mass_comp[0, "X_ac"],
        name="X_ac concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Acetate degraders concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM2D Interface Inlet",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.inlet.conc_mass_comp[0, "X_h2"],
        name="X_h2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Hydrogen degraders concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM2D Interface Inlet",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.inlet.conc_mass_comp[0, "X_I"],
        name="X_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Particulate inerts concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM2D Interface Inlet",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.inlet.conc_mass_comp[0, "X_PHA"],
        name="X_PHA concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Polyhydroxyalkanoates concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM2D Interface Inlet",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.inlet.conc_mass_comp[0, "X_PP"],
        name="X_PP concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Polyphosphates concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM2D Interface Inlet",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.inlet.conc_mass_comp[0, "X_PAO"],
        name="X_PAO concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Phosphorus accumulating organisms concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM2D Interface Inlet",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.inlet.conc_mass_comp[0, "S_K"],
        name="S_K concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Potassium concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM2D Interface Inlet",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.inlet.conc_mass_comp[0, "S_Mg"],
        name="S_Mg concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Magnesium concentration",
        is_input=False,
        is_output=True,
        output_category="ADM1/ASM2D Interface Inlet",
    )

    # Translator outlet
    exports.add(
        obj=fs.translator_adm1_asm2d.outlet.flow_vol[0],
        name="Anaerobic digester outlet flow rate",
        ui_units=pyunits.m**3 / pyunits.day,
        display_units="m3/day",
        rounding=2,
        description="Anaerobic digester outlet flow rate",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Outlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.outlet.conc_mass_comp[0, "S_A"],
        name="S_A concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Fermentation products (acetate) concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Outlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.outlet.conc_mass_comp[0, "S_F"],
        name="S_F concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Fermentable organic substrates concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Outlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.outlet.conc_mass_comp[0, "S_I"],
        name="S_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Soluble inert organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Outlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.outlet.conc_mass_comp[0, "S_IC"],
        name="S_IC concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inorganic carbon concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Outlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.outlet.conc_mass_comp[0, "S_K"],
        name="S_K concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Potassium concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Outlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.outlet.conc_mass_comp[0, "S_Mg"],
        name="S_Mg concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Magnesium concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Outlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.outlet.conc_mass_comp[0, "S_N2"],
        name="S_N2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Dinitrogen product (from denitrification) concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Outlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.outlet.conc_mass_comp[0, "S_NH4"],
        name="S_NH4 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Ammonium plus ammonia nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Outlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.outlet.conc_mass_comp[0, "S_NO3"],
        name="S_NO3 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Nitrate plus nitrite nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Outlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.outlet.conc_mass_comp[0, "S_O2"],
        name="S_O2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Dissolved oxygen concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Outlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.outlet.conc_mass_comp[0, "S_PO4"],
        name="S_PO4 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inorganic soluble phosphorus concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Outlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.outlet.conc_mass_comp[0, "X_AUT"],
        name="X_AUT concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Autotrophic nitrifying biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Outlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.outlet.conc_mass_comp[0, "X_H"],
        name="X_H concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Heterotrophic biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Outlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.outlet.conc_mass_comp[0, "X_I"],
        name="X_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inert particulate organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Outlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.outlet.conc_mass_comp[0, "X_PAO"],
        name="X_PAO concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Phosphate-accumulating biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Outlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.outlet.conc_mass_comp[0, "X_PHA"],
        name="X_PHA concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Polyhydroxyalkanoatess concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Outlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.outlet.conc_mass_comp[0, "X_PP"],
        name="X_PP concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Poly-phosphate concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Outlet (post-interface)",
    )
    exports.add(
        obj=fs.translator_adm1_asm2d.outlet.conc_mass_comp[0, "X_S"],
        name="X_S concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Slowly biodegradable substrate concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Outlet (post-interface)",
    )

    # Anaerobic digestor
    exports.add(
        obj=fs.AD.vapor_outlet.flow_vol[0],
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
        obj=fs.AD.vapor_outlet.conc_mass_comp[0, "S_h2"],
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
        obj=fs.AD.vapor_outlet.conc_mass_comp[0, "S_ch4"],
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
        obj=fs.AD.vapor_outlet.conc_mass_comp[0, "S_co2"],
        name="Anaerobic digester vapor S_co2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Carbon dioxide gas concentration",
        is_input=False,
        is_output=True,
        output_category="Anaerobic Digester Vapor Outlet",
    )

    # Dewatering unit
    exports.add(
        obj=fs.dewater.underflow.flow_vol[0],
        name="Dewatered sludge outlet flow rate",
        ui_units=pyunits.m**3 / pyunits.day,
        display_units="m3/day",
        rounding=2,
        description="Dewatering unit underflow flow rate",
        is_input=False,
        is_output=True,
        output_category="Dewatered Sludge",
    )
    exports.add(
        obj=fs.dewater.underflow.conc_mass_comp[0, "S_A"],
        name="S_A concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Fermentation products (acetate) concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatered Sludge",
    )
    exports.add(
        obj=fs.dewater.underflow.conc_mass_comp[0, "S_F"],
        name="S_F concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Fermentable organic substrates concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatered Sludge",
    )
    exports.add(
        obj=fs.dewater.underflow.conc_mass_comp[0, "S_I"],
        name="S_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Soluble inert organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatered Sludge",
    )
    exports.add(
        obj=fs.dewater.underflow.conc_mass_comp[0, "S_IC"],
        name="S_IC concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inorganic carbon concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatered Sludge",
    )
    exports.add(
        obj=fs.dewater.underflow.conc_mass_comp[0, "S_K"],
        name="S_K concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Potassium concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatered Sludge",
    )
    exports.add(
        obj=fs.dewater.underflow.conc_mass_comp[0, "S_Mg"],
        name="S_Mg concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Magnesium concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatered Sludge",
    )
    exports.add(
        obj=fs.dewater.underflow.conc_mass_comp[0, "S_N2"],
        name="S_N2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Dinitrogen product (from denitrification) concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatered Sludge",
    )
    exports.add(
        obj=fs.dewater.underflow.conc_mass_comp[0, "S_NH4"],
        name="S_NH4 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Ammonium plus ammonia nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatered Sludge",
    )
    exports.add(
        obj=fs.dewater.underflow.conc_mass_comp[0, "S_NO3"],
        name="S_NO3 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Nitrate plus nitrite nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatered Sludge",
    )
    exports.add(
        obj=fs.dewater.underflow.conc_mass_comp[0, "S_O2"],
        name="S_O2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Dissolved oxygen concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatered Sludge",
    )
    exports.add(
        obj=fs.dewater.underflow.conc_mass_comp[0, "S_PO4"],
        name="S_PO4 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inorganic soluble phosphorus concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatered Sludge",
    )
    exports.add(
        obj=fs.dewater.underflow.conc_mass_comp[0, "X_AUT"],
        name="X_AUT concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Autotrophic nitrifying biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatered Sludge",
    )
    exports.add(
        obj=fs.dewater.underflow.conc_mass_comp[0, "X_H"],
        name="X_H concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Heterotrophic biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatered Sludge",
    )
    exports.add(
        obj=fs.dewater.underflow.conc_mass_comp[0, "X_I"],
        name="X_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inert particulate organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatered Sludge",
    )
    exports.add(
        obj=fs.dewater.underflow.conc_mass_comp[0, "X_PAO"],
        name="X_PAO concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Phosphate-accumulating biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatered Sludge",
    )
    exports.add(
        obj=fs.dewater.underflow.conc_mass_comp[0, "X_PHA"],
        name="X_PHA concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Polyhydroxyalkanoates concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatered Sludge",
    )
    exports.add(
        obj=fs.dewater.underflow.conc_mass_comp[0, "X_PP"],
        name="X_PP concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Poly-phosphate concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatered Sludge",
    )
    exports.add(
        obj=fs.dewater.underflow.conc_mass_comp[0, "X_S"],
        name="X_S concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Slowly biodegradable substrate concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatered Sludge",
    )

    exports.add(
        obj=fs.dewater.overflow.flow_vol[0],
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
        obj=fs.dewater.overflow.conc_mass_comp[0, "S_A"],
        name="S_A concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Fermentation products (acetate) concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatering Unit Liquid Outlet",
    )
    exports.add(
        obj=fs.dewater.overflow.conc_mass_comp[0, "S_F"],
        name="S_F concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Fermentable organic substrates concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatering Unit Liquid Outlet",
    )
    exports.add(
        obj=fs.dewater.overflow.conc_mass_comp[0, "S_I"],
        name="S_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Soluble inert organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatering Unit Liquid Outlet",
    )
    exports.add(
        obj=fs.dewater.overflow.conc_mass_comp[0, "S_IC"],
        name="S_IC concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inorganic carbon concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatering Unit Liquid Outlet",
    )
    exports.add(
        obj=fs.dewater.overflow.conc_mass_comp[0, "S_K"],
        name="S_K concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Potassium concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatering Unit Liquid Outlet",
    )
    exports.add(
        obj=fs.dewater.overflow.conc_mass_comp[0, "S_Mg"],
        name="S_Mg concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Magnesium concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatering Unit Liquid Outlet",
    )
    exports.add(
        obj=fs.dewater.overflow.conc_mass_comp[0, "S_N2"],
        name="S_N2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Dinitrogen product (from denitrification) concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatering Unit Liquid Outlet",
    )
    exports.add(
        obj=fs.dewater.overflow.conc_mass_comp[0, "S_NH4"],
        name="S_NH4 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Ammonium plus ammonia nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatering Unit Liquid Outlet",
    )
    exports.add(
        obj=fs.dewater.overflow.conc_mass_comp[0, "S_NO3"],
        name="S_NO3 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Nitrate plus nitrite nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatering Unit Liquid Outlet",
    )
    exports.add(
        obj=fs.dewater.overflow.conc_mass_comp[0, "S_O2"],
        name="S_O2 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Dissolved oxygen concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatering Unit Liquid Outlet",
    )
    exports.add(
        obj=fs.dewater.overflow.conc_mass_comp[0, "S_PO4"],
        name="S_PO4 concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inorganic soluble phosphorus concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatering Unit Liquid Outlet",
    )
    exports.add(
        obj=fs.dewater.overflow.conc_mass_comp[0, "X_AUT"],
        name="X_AUT concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Autotrophic nitrifying biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatering Unit Liquid Outlet",
    )
    exports.add(
        obj=fs.dewater.overflow.conc_mass_comp[0, "X_H"],
        name="X_H concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Heterotrophic biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatering Unit Liquid Outlet",
    )
    exports.add(
        obj=fs.dewater.overflow.conc_mass_comp[0, "X_I"],
        name="X_I concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Inert particulate organic matter concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatering Unit Liquid Outlet",
    )
    exports.add(
        obj=fs.dewater.overflow.conc_mass_comp[0, "X_PAO"],
        name="X_PAO concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Phosphate-accumulating biomass concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatering Unit Liquid Outlet",
    )
    exports.add(
        obj=fs.dewater.overflow.conc_mass_comp[0, "X_PHA"],
        name="X_PHA concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Polyhydroxyalkanoates concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatering Unit Liquid Outlet",
    )
    exports.add(
        obj=fs.dewater.overflow.conc_mass_comp[0, "X_PP"],
        name="X_PP concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Poly-phosphate concentration",
        is_input=False,
        is_output=True,
        output_category="Dewatering Unit Liquid Outlet",
    )
    exports.add(
        obj=fs.dewater.overflow.conc_mass_comp[0, "X_S"],
        name="X_S concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=5,
        description="Slowly biodegradable substrate concentration",
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

    if build_options is not None:
        if build_options["BioP"]:
            bioP = True
        else:
            bioP = False

        m = build(bio_P=bioP)

        set_operating_conditions(m, bio_P=bioP)

        for mx in m.fs.mixers:
            mx.pressure_equality_constraints[0.0, 2].deactivate()
        m.fs.MX3.pressure_equality_constraints[0.0, 2].deactivate()
        m.fs.MX3.pressure_equality_constraints[0.0, 3].deactivate()

        initialize_system(m, bio_P=bioP)

        for mx in m.fs.mixers:
            mx.pressure_equality_constraints[0.0, 2].deactivate()
        m.fs.MX3.pressure_equality_constraints[0.0, 2].deactivate()
        m.fs.MX3.pressure_equality_constraints[0.0, 3].deactivate()

        solve(m)

        # Switch to fixed KLa in R5, R6, and R7 (S_O concentration is controlled in R5)
        m.fs.R5.KLa.fix(240)
        m.fs.R6.KLa.fix(240)
        m.fs.R7.KLa.fix(84)
        m.fs.R5.outlet.conc_mass_comp[:, "S_O2"].unfix()
        m.fs.R6.outlet.conc_mass_comp[:, "S_O2"].unfix()
        m.fs.R7.outlet.conc_mass_comp[:, "S_O2"].unfix()

        # Resolve with controls in place
        solve(m)

        add_costing(m)
        m.fs.costing.initialize()

        # interval_initializer(m.fs.costing)

        assert_degrees_of_freedom(m, 0)

        solve(m)

    else:
        m = build(bio_P=False)

        set_operating_conditions(m, bio_P=False)

        for mx in m.fs.mixers:
            mx.pressure_equality_constraints[0.0, 2].deactivate()
        m.fs.MX3.pressure_equality_constraints[0.0, 2].deactivate()
        m.fs.MX3.pressure_equality_constraints[0.0, 3].deactivate()

        initialize_system(m, bio_P=False)

        for mx in m.fs.mixers:
            mx.pressure_equality_constraints[0.0, 2].deactivate()
        m.fs.MX3.pressure_equality_constraints[0.0, 2].deactivate()
        m.fs.MX3.pressure_equality_constraints[0.0, 3].deactivate()

        solve(m)

        # Switch to fixed KLa in R5, R6, and R7 (S_O concentration is controlled in R5)
        m.fs.R5.KLa.fix(240)
        m.fs.R6.KLa.fix(240)
        m.fs.R7.KLa.fix(84)
        m.fs.R5.outlet.conc_mass_comp[:, "S_O2"].unfix()
        m.fs.R6.outlet.conc_mass_comp[:, "S_O2"].unfix()
        m.fs.R7.outlet.conc_mass_comp[:, "S_O2"].unfix()

        # Resolve with controls in place
        solve(m)

        add_costing(m)
        m.fs.costing.initialize()

        interval_initializer(m.fs.costing)

        assert_degrees_of_freedom(m, 0)

        solve(m)

    return m


def solve_flowsheet(flowsheet=None):
    """
    Solves the initial flowsheet.
    """
    fs = flowsheet
    results = solve(fs)
    return results
