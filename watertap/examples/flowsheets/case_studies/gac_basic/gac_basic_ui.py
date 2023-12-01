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

import gac_basic_flowsheet as gac_fs

from pyomo.environ import units as pyunits
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

    # --- Input data ---
    # Feed conditions
    exports.add(
        obj=fs.feed.flow_vol[0],
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
