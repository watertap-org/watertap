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

__author__ = "Adam Atia"

import pandas as pd

df = pd.read_excel("WT3_unit_classification_for_doc.xlsx")

unit_name_list = [i.title() for i in df["Name"]]
model_type_list = df["model type long"]
model_type_ref_list = df["model type doc ref"]
elect_func_list = df["energy"]
cost_func_list = df["Cost function"]
zo_name_list = df["zo_unit"]

with open("index.rst", "w") as f:
    f.write("Zero-Order Unit Models\n")
    f.write("=" * len("Zero-Order Unit Models"))
    f.write("\n")
    f.write(".. toctree::\n")
    f.write("   :maxdepth: 1\n\n")

for i, u in enumerate(unit_name_list):

    list = [
        f"This unit model is formulated as a {model_type_list[i]} model form.",
        f"See documentation for :ref:`{model_type_list[i]} Helper Methods<{model_type_ref_list[i]}>`.",
        f"Electricity consumption is calculated using the {elect_func_list[i]} helper function.",
        f"Costing is calculated using the {cost_func_list[i]} method in the zero-order costing package.",
        f"   pair: watertap.unit_models.zero_order.{zo_name_list[i]};{zo_name_list[i]}",
        f".. currentmodule:: watertap.unit_models.zero_order.{zo_name_list[i]}",
        f".. automodule:: watertap.unit_models.zero_order.{zo_name_list[i]}",
    ]
    with open("index.rst", "a") as f:
        f.write(f"   {zo_name_list[i]}\n")

    with open(f"{zo_name_list[i]}.rst", "w") as f:
        f.write(f"{unit_name_list[i]} (ZO)")
        count = 0
        f.write("\n")
        f.write("=" * (5 + len(u)))
        f.write("\n")
        f.write("\nModel Type\n")
        f.write("-" * len("Model Type"))
        f.write(f"\n{list[count]}")
        count += 1
        if not (
            zo_name_list[i] == "feed_zo" or zo_name_list[i] == "gas_sparged_membrane_zo"
        ):

            f.write(f"\n{list[count]}\n")
        else:
            f.write("\n")
        count += 1
        f.write("\nElectricity Consumption\n")
        f.write("-" * len("Electricity Consumption"))
        f.write(f"\n{list[count]}")
        count += 1
        f.write(
            f"\nSee documentation for :ref:`Helper Methods for Electricity Demand<electricity_methods>`.\n"
        )
        f.write("\nCosting Method\n")
        f.write("-" * len("Costing Method"))
        f.write(f"\n{list[count]}")
        count += 1
        f.write(
            f"\nSee documentation for the :ref:`zero-order costing package<zero_order_costing>`.\n"
        )
        f.write("\n.. index::")
        f.write(f"\n{list[count]}\n")
        count += 1
        f.write(f"\n{list[count]}\n")
        count += 1
        f.write("\nClass Documentation\n")
        f.write("-" * len("Class Documentation"))
        f.write(f"\n\n{list[count]}\n")
        f.write("    :members:\n")
        f.write("    :noindex:\n")
