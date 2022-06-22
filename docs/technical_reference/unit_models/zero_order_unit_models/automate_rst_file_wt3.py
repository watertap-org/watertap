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
from watertap.unit_models.zero_order import *

df = pd.read_excel("WT3_unit_classification_for_doc.xlsx")

unit_name_list = [i.title() for i in df["Name"]]
model_type_list = df["model type long"]
model_type_ref_list = df["model type doc ref"]
elect_func_list = df["energy"]
cost_func_list = df["Cost function"]
zo_name_list = df["zo_unit"]
energy_helper_list = df["energy_helper_func"]
class_list = df["class"]

# model doc exceptions: keys= zo_unit, value= custom f string
title_exceptions = {
    "anaerobic_mbr_mec_zo": "Integrated Anaerobic Membrane Bioreactor/Microbial Electrolysis Cell",
    "CANDOP_zo": "CANDO-P",
    "co2_addition_zo": "CO2 Addition",
    "dmbr_zo": "Recirculating Dynamic Membrane Bioreactor",
    "gac_zo": "Granular Activated Carbon",
    "mabr_zo": "Membrane Aerated Biofilm Reactor",
    "mbr_zo": "Membrane Bioreactor",
    "metab_zo": "Modular Encapsulated Two-stage Anaerobic Biological Reactor",
    "municipal_wwtp_zo": "Municipal Wastewater Treatment Plant",
    "ozone_aop_zo": "Ozone with Advanced Oxidation Processes",
    "secondary_treatment_wwtp_zo": "Secondary Wastewater Treatment Plant",
    "sw_onshore_intake_zo": "Seawater Onshore Intake",
    "uv_aop_zo": "UV with Advanced Oxidation Processes",
    "uv_zo": "UV Reactor",
    "vfa_recovery_zo": "Volatile Fatty Acid (VFA) Recovery Unit",
    "waiv_zo": "Wind-Aided Intensified Evaporation Unit",
}

model_type_exceptions = {}

elec_func_exceptions = {}

costing_exceptions = {}

has_subtype = {}


# Create index file for all zero order model docs
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

    # append unit doc to index
    with open("index.rst", "a") as f:
        f.write(f"   {zo_name_list[i]}\n")

    with open(f"{zo_name_list[i]}.rst", "w") as f:

        # write doc title based on unit name
        if zo_name_list[i] in title_exceptions:
            f.write(f"{title_exceptions[zo_name_list[i]]} (ZO)")
            f.write("\n")
            f.write("=" * len(f"{title_exceptions[zo_name_list[i]]} (ZO)"))
        else:
            f.write(f"{unit_name_list[i]} (ZO)")
            f.write("\n")
            f.write("=" * (5 + len(u)))
        f.write("\n")
        count = 0

        # write Model Type section
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

        # write Electricity Consumption section
        f.write("\nElectricity Consumption\n")
        f.write("-" * len("Electricity Consumption"))
        f.write(f"\n{list[count]}")
        count += 1
        f.write(
            f"\nSee documentation for :ref:`Helper Methods for Electricity Demand<electricity_methods>`.\n"
        )

        # write Costing Method section
        f.write("\nCosting Method\n")
        f.write("-" * len("Costing Method"))
        f.write(f"\n{list[count]}")
        count += 1
        f.write(
            f"\nSee documentation for the :ref:`zero-order costing package<zero_order_costing>`.\n"
        )

        # write Additional Variables section if unit is non-basic
        if class_list[i] == "non-basic":
            f.write("\nAdditional Variables\n")
            f.write("-" * len("Additional Variables"))
            f.write("\n")

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


def grab_unit_variables(unit):
    from pyomo.environ import ConcreteModel
    from idaes.core import FlowsheetBlock
    from watertap.core.zero_order_properties import WaterParameterBlock

    m = ConcreteModel()
    m.fs = FlowsheetBlock()
    m.fs.props = WaterParameterBlock({"solute_list": ["tss", "cod", "tds"]})
