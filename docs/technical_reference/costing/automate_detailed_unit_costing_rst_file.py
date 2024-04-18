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

__author__ = "Adam Atia"

import pandas as pd

df = pd.read_csv("WT_unit_costing_for_docs.csv")

unit_title_list = [i.title() for i in df["name"]]

unit_name_list = df["unit"]

title_exceptions = {
    "cstr_injection": "Completely Stirred Tank Reactor w/Injection Stream",
    "cstr": "Completely Stirred Tank Reactor",
    "dewatering": "Dewatering Unit",
    "electroNP": "ElectroN-P",
    "gac": "Granular Activated Carbon",
    "heater_chiller": "Heater/Chiller",
    "uv_aop": "UV with Advanced Oxidation Processes",
}


if __name__ == "__main__":

    # Create index file for all unit model costing docs
    with open("detailed_unit_model_costing.rst", "w") as f:
        f.write("Detailed Unit Model Costing\n")
        f.write("=" * len("Detailed Unit Model Costing"))
        f.write("\n")
        f.write(
            "Default costing methods are provided for the unit models listed below. However, users should supply their own cost relationship, if possible, instead of relying completely on the defaults.\n\n"
        )
        f.write(".. toctree::\n")
        f.write("   :maxdepth: 1\n\n")

    for i, u in enumerate(unit_title_list):

        # list = [
        #     f"This unit model is formulated as a **{model_type_list[i]}** model form.",
        #     f"See documentation for :ref:`{model_type_list[i]} Helper Methods<{model_type_ref_list[i]}>`.",
        #     f"Electricity consumption is calculated using the **{elect_func_list[i]}** helper function.",
        #     None,  # created on the fly now
        #     f"   pair: watertap.unit_models.zero_order.{unit_name_list[i]};{unit_name_list[i]}",
        #     f".. currentmodule:: watertap.unit_models.zero_order.{unit_name_list[i]}",
        #     f".. automodule:: watertap.unit_models.zero_order.{unit_name_list[i]}",
        # ]

        # append unit doc to index
        with open("detailed_unit_model_costing.rst", "a") as f:
            f.write(f"   {unit_name_list[i]}\n")

        with open(f"{unit_name_list[i]}.rst", "w", encoding="utf-8") as f:
            # write doc title based on unit name
            if unit_name_list[i] in title_exceptions:
                f.write(f"{title_exceptions[unit_name_list[i]]} Costing Method")
                f.write("\n")
                f.write("=" * len(f"{title_exceptions[unit_name_list[i]]} Costing Method"))
            else:
                f.write(f"{unit_title_list[i]} Costing Method")
                f.write("\n")
                f.write("=" * len(f"{unit_title_list[i]} Costing Method"))
            f.write("\n")

            #TODO: add parameter tables
            f.write("\nCosting Method Parameters\n")
            f.write("+" * len("Costing Method Parameters"))
            f.write(f"\n\nThe following parameters are constructed when applying the `cost_{unit_name_list[i]}` costing method in the ``watertap_costing_package``:\n\n")
            
            #TODO: add var tables
            f.write("\n\nCosting Method Variables\n")
            f.write("+" * len("Costing Method Variables"))
            f.write(f"\n\nThe following variables are constructed when applying the `cost_{unit_name_list[i]}` costing method in the ``watertap_costing_package``:\n\n")

            #TODO: add capex eqs
            f.write("\n\nCapital Cost Calculations\n")
            f.write("+" * len("Capital Cost Calculations"))

            #TODO: add opex eqs
            f.write("\n\nOperating Cost Calculations\n")
            f.write("+" * len("Operating Cost Calculations"))

            #TODO: add module directives to unit and cost method 
            f.write("\n\nCode Documentation\n")
            f.write("-" * len("Code Documentation"))

            f.write("\n\nReferences\n")
            f.write("-" * len("References"))
            f.write("\n")
