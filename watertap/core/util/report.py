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

__author__ = "Alexander V. Dudchenko"

from pyomo.common.formatting import tabular_writer
from pyomo.environ import (
    value,
    units as pyunits,
)
import sys
from idaes.core.util.units_of_measurement import report_quantity


def build_report_table(
    unit_name, data_dict, ostream=None, prefix="", use_default_units=False
):
    """Builds a report table for the unit model
    Developer should define a data dict as follows:
        data_dict = {
            "Composition": {},
            "Physical state": {},
        }
        for phase, ion in self.feed.properties[0].conc_mass_phase_comp:
            data_dict["Composition"][ion] = self.feed.properties[
                0
            ].conc_mass_phase_comp[phase, ion]
        data_dict["Physical state"]["pH"] = self.feed.pH
        data_dict["Physical state"]["Temperature"] = self.feed.properties[
            0
        ].temperature
        data_dict["Physical state"]["Pressure"] = self.feed.properties[0].pressure

    ------------------------------------------------------------------------------------
        fs.feed state

        Composition:
        Key    : Value     : Units                 : Fixed : Bounds
         Ca_2+ :   0.25800 : kilogram / meter ** 3 :  True : (0, 2000.0)
          Cl_- :   0.87000 : kilogram / meter ** 3 :  True : (0, 2000.0)
           H2O :    996.64 : kilogram / meter ** 3 : False : (0, 2000.0)
        HCO3_- :   0.38500 : kilogram / meter ** 3 :  True : (0, 2000.0)
           K_+ : 0.0090000 : kilogram / meter ** 3 :  True : (0, 2000.0)
         Mg_2+ :  0.090000 : kilogram / meter ** 3 :  True : (0, 2000.0)
          Na_+ :   0.73900 : kilogram / meter ** 3 :  True : (0, 2000.0)
        SO4_2- :    1.0110 : kilogram / meter ** 3 :  True : (0, 2000.0)

        Physical state:
        Key         : Value      : Units         : Fixed : Bounds
           Pressure : 1.0000e+05 :        pascal :  True : (100000.0, None)
        Temperature :     293.15 :        kelvin :  True : (273.15, 373.15)
                 pH :     7.0700 : dimensionless :  True : (None, None)

    ------------------------------------------------------------------------------------

    Args:
        unit_name: Name of the unit
        data_dict: Dictionary containing data to report
        ostream: Output stream for the report (default: sys.stdout)
        prefix: String prefix for formatting the report
        use_default_units: Boolean to indicate if default units should be used
    """

    def _get_fixed_state(v):
        """Get the fixed state of a variable, if none exists then return N/A"""
        try:
            return v.fixed
        except AttributeError:
            return "N/A"

    def _get_bounds(v):
        """Get the bounds of a variable, if none exists then return N/A"""
        try:
            return v.bounds
        except AttributeError:
            return "N/A"

    def get_values(k, v):
        """Get the values of a variable,  dimensions, fixed state and bounds"""
        if isinstance(v, int):
            return [
                v,
                "dimensionless",
                "N/A",
                "N/A",
            ]
        elif isinstance(v, float):
            return [
                "{:#.5g}".format(v),
                "dimensionless",
                "N/A",
                "N/A",
            ]
        elif use_default_units:
            return [
                "{:#.5g}".format(report_quantity(v).m),
                report_quantity(v).u,
                _get_fixed_state(v),
                _get_bounds(v),
            ]
        else:
            return [
                "{:#.5g}".format(value(v)),
                pyunits.get_units(v),
                _get_fixed_state(v),
                _get_bounds(v),
            ]

    if ostream is None:
        ostream = sys.stdout
    max_str_length = 84
    tab = " " * 4
    ostream.write("\n" + "-" * max_str_length + "\n")
    ostream.write(f"{prefix}{tab}{unit_name} state")
    ostream.write("\n" * 2)
    for key, sub_data in data_dict.items():

        ostream.write(f"{prefix}{tab}{key}: \n")
        tabular_writer(
            ostream,
            prefix + tab,
            ((k, v) for k, v in sub_data.items()),
            ("Value", "Units", "Fixed", "Bounds"),
            get_values,
        )
        ostream.write(f"\n")

    ostream.write("-" * max_str_length + "\n")
