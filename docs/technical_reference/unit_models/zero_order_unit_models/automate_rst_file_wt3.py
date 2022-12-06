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
from idaes.core import declare_process_block_class
from watertap.core import ZeroOrderBaseData
from pyomo.environ import ConcreteModel, Var, Constraint
from watertap.core.wt_database import Database
from idaes.core import FlowsheetBlock
from watertap.core.zero_order_properties import WaterParameterBlock
import watertap.unit_models.zero_order as zo
from watertap.core import ZeroOrderBaseData
from watertap.core.tests.test_zero_order_base import DerivedZOBase
from watertap.core import build_pt, build_sido, build_siso, build_sido_reactive
from watertap.core import pump_electricity, constant_intensity
from pyomo.environ import Reference
import os
from glob import glob
from pathlib import Path


sidor_db_path = os.path.dirname(os.path.abspath(__file__))


def grab_unit_components(unit_class):

    m = ConcreteModel()
    m.zero_db = Database(dbpath=sidor_db_path)
    m.db = Database()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.props = WaterParameterBlock(
        solute_list=[
            "toc",
            "tss",
            "tkn",
            "cod",
            "tds",
            "nitrogen",
            "phosphates",
            "phosphorus",
            "struvite",
            "nonbiodegradable_cod",
            "hydrogen",
            "ammonium_as_nitrogen",
            "nitrate",
            "bod",
            "organic_solid",
            "organic_liquid",
            "iron",
            "peracetic_acid",
            "total_coliforms_fecal_ecoli",
        ]
    )
    unit = getattr(zo, unit_class)
    if unit_class == "MetabZO":
        p_subtype = "hydrogen"
    else:
        p_subtype = "default"
    m.fs.unit = unit(
        property_package=m.fs.props,
        database=m.db,
        process_subtype=p_subtype,
    )

    if model_type_short_list[i] == "SIDO reactive":
        zdb = m.zero_db
    else:
        zdb = m.db
    m.fs.zo_base = DerivedZOBase(
        property_package=m.fs.props,
        database=zdb,
        # process_subtype = process_subtype,
    )

    if model_type_short_list[i] == "PT":
        build_pt(m.fs.zo_base)
    elif model_type_short_list[i] == "SIDO":
        build_sido(m.fs.zo_base)
    elif model_type_short_list[i] == "SISO":
        build_siso(m.fs.zo_base)
    elif model_type_short_list[i] == "SIDO reactive":
        m.fs.zo_base._tech_type = "dummy_sidor_data"
        build_sido_reactive(m.fs.zo_base)
    else:
        pass

    if elect_func_list[i] == "pump_electricity":
        if model_type_short_list[i] != "PT":
            m.fs.zo_base._Q = Reference(m.fs.zo_base.properties_in[:].flow_vol)
        else:
            m.fs.zo_base._Q = Reference(m.fs.zo_base.properties[:].flow_vol)

        pump_electricity(m.fs.zo_base, m.fs.zo_base._Q)
    elif elect_func_list[i] == "constant_intensity":
        constant_intensity(m.fs.zo_base)
    else:
        pass

    zo_base_vars = []
    for var in m.fs.zo_base.component_data_objects(Var, descend_into=False):
        zovarname = var.name
        zo_base_vars.append(zovarname.replace("fs.zo_base.", "").split("[", 1)[0])
    zo_base_cons = []
    for con in m.fs.zo_base.component_data_objects(Constraint, descend_into=False):
        zoconame = con.name
        zo_base_cons.append(zoconame.replace("fs.zo_base.", "").split("[", 1)[0])
    added_vars = []
    added_var_docs = []
    added_var_units = []
    added_cons = []
    added_con_docs = []
    for var in m.fs.unit.component_data_objects(Var, descend_into=False):
        addedvarname = var.name
        newname = addedvarname.replace("fs.unit.", "").split("[", 1)[0]
        if newname not in zo_base_vars:
            model_var = getattr(m.fs.unit, newname)
            added_vars.append(newname)
            added_var_docs.append(model_var.doc)
            added_var_units.append(str(model_var._units).replace("'", ""))
            if added_var_units[-1] == "None":
                added_var_units[-1] = "dimensionless"
            if "**" in added_var_units[-1]:
                added_var_units[-1] = (":math:" + f"`{added_var_units[-1]}`").replace(
                    "**", "^"
                )
            else:
                added_var_units[-1] = ":math:" + f"`{added_var_units[-1]}`"

    for con in m.fs.unit.component_data_objects(Constraint, descend_into=False):
        addedconame = con.name
        connewname = addedconame.replace("fs.unit.", "").split("[", 1)[0]
        if connewname not in zo_base_cons:
            model_con = getattr(m.fs.unit, connewname)
            added_cons.append(connewname)
            added_con_docs.append(model_con.doc)
    return (
        m,
        zo_base_vars,
        added_vars,
        added_var_docs,
        added_var_units,
        added_cons,
        added_con_docs,
    )


def grab_unit_components_feed(unit_class):

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.props = WaterParameterBlock(
        solute_list=[
            "toc",
            "tkn",
            "tss",
            "cod",
            "tds",
            "nitrogen",
            "phosphates",
            "phosphorus",
            "struvite",
            "nonbiodegradable_cod",
            "hydrogen",
            "ammonium_as_nitrogen",
            "nitrate",
            "bod",
        ]
    )
    unit = getattr(zo, unit_class)
    m.fs.unit = zo.FeedZO(
        property_package=m.fs.props,
    )

    added_vars = []
    added_var_docs = []
    added_var_units = []
    added_cons = []
    added_con_docs = []
    for var in m.fs.unit.component_data_objects(Var, descend_into=False):

        addedvarname = var.name
        newname = addedvarname.replace("fs.unit.", "").split("[", 1)[0]
        model_var = getattr(m.fs.unit, newname)

        if not (newname == "properties"):
            if newname not in added_vars:
                added_vars.append(newname)

                added_var_docs.append(model_var.doc)
                added_var_units.append(str(model_var._units).replace("'", ""))
                if added_var_units[-1] == "None":
                    added_var_units[-1] = "dimensionless"
                if "**" in added_var_units[-1]:
                    added_var_units[-1] = (
                        ":math:" + f"`{added_var_units[-1]}`"
                    ).replace("**", "^")
                else:
                    added_var_units[-1] = ":math:" + f"`{added_var_units[-1]}`"

    for con in m.fs.unit.component_data_objects(Constraint, descend_into=False):
        addedconame = con.name
        connewname = addedconame.replace("fs.unit.", "").split("[", 1)[0]
        model_con = getattr(m.fs.unit, connewname)
        if connewname not in added_cons:
            added_cons.append(connewname)
            added_con_docs.append(model_con.doc)

    return (
        m,
        added_vars,
        added_var_docs,
        added_var_units,
        added_cons,
        added_con_docs,
    )


df = pd.read_excel("WT3_unit_classification_for_doc.xlsx")

unit_name_list = [i.title() for i in df["Name"]]
model_type_list = df["model type long"]
model_type_short_list = df["type"]

model_type_ref_list = df["model type doc ref"]
elect_func_list = df["energy"]
cost_func_list = df["Cost function"]
zo_name_list = df["zo_unit"]
class_name_list = df["class_name"]
energy_helper_list = df["energy_helper_func"]
class_list = df["class"]

# model doc exceptions: keys= zo_unit, value= custom f string
title_exceptions = {
    "anaerobic_mbr_mec_zo": "Integrated Anaerobic Membrane Bioreactor/Microbial Electrolysis Cell",
    "CANDOP_zo": "CANDO-P",
    "co2_addition_zo": "CO2 Addition",
    "dmbr_zo": "Recirculating Dynamic Membrane Bioreactor",
    "gac_zo": "Granular Activated Carbon",
    "hrcs_zo": "High-Rate Contact Stabilization",
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

p_subtype_exceptions = {"MetabZO": "hydrogen"}
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
        f"This unit model is formulated as a **{model_type_list[i]}** model form.",
        f"See documentation for :ref:`{model_type_list[i]} Helper Methods<{model_type_ref_list[i]}>`.",
        f"Electricity consumption is calculated using the **{elect_func_list[i]}** helper function.",
        f"Costing is calculated using the **{cost_func_list[i]}** method in the zero-order costing package.",
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

        if zo_name_list[i] == "feed_zo":
            f.write(
                "\nThe Feed (ZO) model adds volumetric flowrate and mass concentration "
                "as variables in connection with the zero-order property model's state "
                "variable, mass flow rate. "
                "This allows the user to enter feed volumetric flow rate and "
                "concentrations without calculating the mass flowrates of each "
                "component in the feed.\n"
            )

        f.write("\nModel Type\n")
        f.write("-" * len("Model Type"))
        if not (zo_name_list[i] == "feed_zo"):
            # write Model Type section
            f.write(f"\n{list[0]}")
            f.write(f"\n{list[1]}\n")
        else:
            f.write(
                "\nThe Feed (ZO) block for zero-order flowsheets contains "
                "methods for getting concentration data from the database, and it "
                "has been created to work with the zero-order property package.\n"
            )

        if not (zo_name_list[i] == "feed_zo") and not (
            zo_name_list[i] == "constructed_wetlands_zo"
        ):
            # write Electricity Consumption section
            f.write("\nElectricity Consumption\n")
            f.write("-" * len("Electricity Consumption"))
            if (elect_func_list[i] != "pump_electricity") and (
                elect_func_list[i] != "constant_intensity"
            ):
                (
                    _,
                    _,
                    _,
                    _,
                    _,
                    addedconscheck,
                    _,
                ) = grab_unit_components(class_name_list[i])

                if len(addedconscheck) > 0:
                    f.write(
                        "\nThe constraint used to calculate energy consumption is described in the Additional Constraints section below. More details can be found in the unit model class.\n"
                    )
            else:
                f.write(f"\n{list[2]}")
                f.write(
                    f"\nSee documentation for :ref:`Helper Methods for Electricity Demand<electricity_methods>`.\n"
                )

        if not (zo_name_list[i] == "feed_zo"):
            # write Costing Method section
            f.write("\nCosting Method\n")
            f.write("-" * len("Costing Method"))
            if not cost_func_list[i]:
                f.write("\nThis unit does not include costing.\n")
            else:
                f.write(f"\n{list[3]}")
                f.write(
                    f"\nSee documentation for the :ref:`zero-order costing package<zero_order_costing>`.\n"
                )

        # write Additional Variables section if unit is non-basic
        # TODO: conditional setting section to Variables if custom model type; add indices?; Add constraints section
        # if class_list[i] == "non-basic":
        if not (zo_name_list[i] == "feed_zo"):
            print(class_name_list[i])
            (
                _,
                _,
                addedvars,
                vardocs,
                varunits,
                addedcons,
                condocs,
            ) = grab_unit_components(class_name_list[i])
            if len(addedvars) > 0:
                f.write("\nAdditional Variables\n")
                f.write("-" * len("Additional Variables"))
                f.write("\n\n")
                f.write(".. csv-table::\n")
                f.write('   :header: "Description", "Variable Name", "Units"\n\n')
            for k, v in enumerate(addedvars):
                f.write(f'   "{vardocs[k]}", "{v}", "{varunits[k]}"\n')

            # write Additional Constraints section if unit is non-basic
            if len(addedcons) > 0:
                f.write("\nAdditional Constraints\n")
                f.write("-" * len("Additional Constraints"))
                f.write("\n\n")
                f.write(".. csv-table::\n")
                f.write('   :header: "Description", "Constraint Name"\n\n')
                for k, c in enumerate(addedcons):
                    f.write(f'   "{condocs[k]}", "{c}"\n')

        else:
            print(class_name_list[i])
            (
                _,
                addedvars,
                vardocs,
                varunits,
                addedcons,
                condocs,
            ) = grab_unit_components_feed(class_name_list[i])
            if len(addedvars) > 0:
                f.write("\nAdditional Variables\n")
                f.write("-" * len("Additional Variables"))
                f.write("\n\n")
                f.write(".. csv-table::\n")
                f.write('   :header: "Description", "Variable Name", "Units"\n\n')
            for k, v in enumerate(addedvars):
                f.write(f'   "{vardocs[k]}", "{v}", "{varunits[k]}"\n')

            # write Additional Constraints section if unit is non-basic
            if len(addedcons) > 0:
                f.write("\nAdditional Constraints\n")
                f.write("-" * len("Additional Constraints"))
                f.write("\n\n")
                f.write(".. csv-table::\n")
                f.write('   :header: "Description", "Constraint Name"\n\n')
                for k, c in enumerate(addedcons):
                    f.write(f'   "{condocs[k]}", "{c}"\n')

        if not (zo_name_list[i] == "feed_zo"):
            f.write("\n.. index::")
            f.write(f"\n{list[4]}\n")
            f.write(f"\n{list[5]}\n")
            f.write("\nClass Documentation\n")
            f.write("-" * len("Class Documentation"))
            f.write(f"\n\n{list[6]}\n")
            f.write("    :members:\n")
            f.write("    :noindex:\n")
        else:
            f.write("\nClass Documentation\n")
            f.write("-" * len("Class Documentation"))
            f.write(f"\n\n{list[6]}\n")
            f.write("    :members:\n")
            f.write("    :noindex:\n")
