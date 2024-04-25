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

unit_name_list = list(df["unit"])

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
    # Toggle OVERWRITE to True to overwrite unit model costing rst files
    OVERWRITE = True
    # Specify unit model costing filenames (e.g., "ion_exchange") that you want to overwrite
    # Otherwise, if list left empty, the default is to overwrite all files
    OVERWRITE_LIST = []
    try:
        with open("detailed_unit_model_costing.rst", "r") as f:
            lines = f.readlines()
    except:

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
        with open("detailed_unit_model_costing.rst", "r") as f:
            lines = f.readlines()

    for i, unit_title in enumerate(unit_title_list):
        for l, line in enumerate(lines):
            # set on_list flag to False in unit name not found in landing page list
            if unit_name_list[i] != line.strip():
                on_list = False
            # set on_list flag to True if unit name is on landing page list and break inner loop
            else:
                on_list = True
                line_index = l
                break

        # if new entry or if overwriting existing entry
        if (not on_list) or (OVERWRITE and (on_list)):

            # check that entries in OVERWRITE_LIST are valid
            if OVERWRITE and len(OVERWRITE_LIST):
                if all(
                    overwrite_item not in unit_name_list
                    for overwrite_item in OVERWRITE_LIST
                ):
                    raise ValueError(
                        "Entries provided in OVERWRITE_LIST are invalid.\nOVERWRITE_LIST should be empty or contain a list of valid unit names (e.g., 'ion_exchange')"
                    )

            # if this is a new entry, append the unit name to the landing page for all unit costing
            if not on_list:
                # append unit doc to index
                with open("detailed_unit_model_costing.rst", "a") as f:
                    f.write(f"   {unit_name_list[i]}\n")
            else:
                # either overwrite all units in the landing page OR overwrite specific units specified on OVERWRITE_LIST
                if (not len(OVERWRITE_LIST)) or (
                    any(item == unit_name_list[i] for item in OVERWRITE_LIST)
                ):
                    lines[line_index] = f"   {unit_name_list[i]}\n"
                    with open("detailed_unit_model_costing.rst", "w") as f:
                        f.writelines(lines)
                else:
                    pass

            if (
                (not on_list)
                or (not len(OVERWRITE_LIST))
                or (any(item == unit_name_list[i] for item in OVERWRITE_LIST))
            ):
                with open(f"{unit_name_list[i]}.rst", "w") as f:

                    # write doc title based on unit name
                    if unit_name_list[i] in title_exceptions:
                        f.write(f"{title_exceptions[unit_name_list[i]]} Costing Method")
                        f.write("\n")
                        f.write(
                            "="
                            * len(
                                f"{title_exceptions[unit_name_list[i]]} Costing Method\n"
                            )
                        )
                    else:
                        f.write(f"{unit_title} Costing Method")
                        f.write("\n")
                        f.write("=" * len(f"{unit_title} Costing Method\n"))

                    # TODO: add parameter tables
                    f.write("\n\nCosting Method Parameters\n")
                    f.write("+" * len("Costing Method Parameters"))
                    f.write(
                        f"\n\nThe following parameters are constructed when applying the `cost_{unit_name_list[i]}` costing method in the ``watertap_costing_package``:\n\n"
                    )
                    f.write(".. csv-table::\n")
                    f.write(
                        '   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units"\n\n'
                    )
                    f.write(
                        '   "description", ":math:`Symbol_{example}`", "parameter_name", "1", ":math:`\\text{dimensionless}`"\n'
                    )

                    # TODO: add var tables
                    f.write("\nCosting Method Variables\n")
                    f.write("+" * len("Costing Method Variables"))
                    f.write(
                        f"\n\nThe following variables are constructed when applying the `cost_{unit_name_list[i]}` costing method in the ``watertap_costing_package``:\n\n"
                    )
                    f.write(".. csv-table::\n")
                    f.write(
                        '   :header: "Description", "Symbol", "Variable Name", "Default Value", "Units"\n\n'
                    )
                    f.write(
                        '   "description", ":math:`Symbol_{example}`", "variable_name", "1", ":math:`\\text{dimensionless}`"\n'
                    )

                    # TODO: add capex eqs
                    f.write("\nCapital Cost Calculations\n")
                    f.write("+" * len("Capital Cost Calculations"))
                    f.write(
                        "\n\nDescribe capital costs..keep it concise where possible\n\n"
                    )
                    f.write("    .. math::\n\n")
                    f.write(
                        "        C_{cap,tot} = C_{cap,example1}+C_{cap,example2}+C_{cap,other}"
                    )
                    f.write("\n\n    .. math::\n\n")
                    f.write(
                        "        C_{cap,example1} = fill in equation for each component in total capex equation\n\n "
                    )

                    # TODO: add opex eqs
                    f.write("\nOperating Cost Calculations\n")
                    f.write("+" * len("Operating Cost Calculations"))
                    f.write(
                        "\n\nDescribe operating/maintenance costs..keep it concise where possible\n\n"
                    )
                    f.write("    .. math::\n\n")
                    f.write(
                        "        C_{op,tot} = C_{op,example1}+C_{op,example2}+C_{op,other}"
                    )
                    f.write("\n\n    .. math::\n\n")
                    f.write(
                        "        C_{op,example1} = fill in equation for each component in total opex equation\n\n "
                    )

                    # TODO: add module directives to unit and cost method
                    f.write("\nCode Documentation\n")
                    f.write("-" * len("Code Documentation"))
                    f.write(
                        f"\n\n* :mod:`watertap.costing.unit_models.{unit_name_list[i]}`"
                    )

                    f.write("\n\nReferences\n")
                    f.write("-" * len("References"))
                    f.write(
                        "\nAim to include at least one reference in most cases, but delete this section if no references used for cost relationships/default values"
                    )
