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
This file contains methods to convert WaterTAP naming conventions to OLI
and generate molecular weight and charge dictionaries from molecular formulae.

It calculates molecular weights using the periodic_table.csv from:
https://gist.github.com/GoodmanSciences/c2dd862cd38f21b0ad36b8f96b4bf1ee.
"""

__author__ = "Paul Vecchiarelli, Ben Knueven, Adam Atia"

from collections import namedtuple

from watertap.core.util.chemistry import (
    get_charge,
    get_charge_group,
    get_molar_mass,
)

# TODO: consider replacing some functionality with molmass: https://pypi.org/project/molmass

OLIName = namedtuple(
    "OLIName", ["oli_name", "watertap_name", "charge", "charge_group", "molar_mass"]
)


def watertap_to_oli(watertap_name: str) -> OLIName:
    """
    This method creates a named tuple
    which can be passed directly into OLI
    or into MCAS property models.

    :param watertap_name: string name of substance in WaterTAP format, i.e., B[OH]4_-

    :return OLIName: named tuple containing attributes derived from molecular formula
    """

    c = findall(r"[A-Z]", watertap_name)
    if len(c) == 0:
        raise IOError(
            f" At least 1 uppercase letter is required to specify a molecule, not '{watertap_name}'."
        )

    oli_name = get_oli_name(watertap_name)
    charge = get_charge(watertap_name)
    charge_group = get_charge_group(charge)
    molar_mass = get_molar_mass(watertap_name)
    return OLIName(oli_name, watertap_name, charge, charge_group, molar_mass)


def get_oli_name(watertap_name: str) -> str:
    """
    Converts an WaterTAP formatted name, i.e., "Na_+"
    into an OLI formatted name, i.e., "NAION".

    :param watertap_name: string name of a solute in WaterTAP format

    :return oli_name: string name of a solute in OLI format
    """

    exclude_items = ["temperature", "pressure", "volume"]
    if watertap_name.lower() in exclude_items:
        return watertap_name
    if hasattr(watertap_name, "oli_name"):
        return watertap_name
    components = watertap_name.split("_")
    if len(components) == 0:
        raise IOError(f" Unable to parse solute '{watertap_name}'.")
    if len(components) == 1:
        molecule = components[0]
    elif len(components) == 2:
        molecule = components[0] + "ION"
    oli_name = molecule.replace("[", "").replace("]", "").upper()
    return oli_name


def get_oli_names(source: dict):
    """
    Updates source dictionary with data to populate MCAS property model.

    :param source: dictionary containing WaterTAP names as keys

    :return source: dictionary with OLIName named tuples as keys
    """

    source = dict(
        map(lambda k, v: (watertap_to_oli(k), v), source.keys(), source.values())
    )
    return source


def oli_reverse_lookup(oli_name: str, names_db) -> OLIName:
    """
    Looks up WaterTAP formatted name for solute in OLI format, if listed in names_db dictionary.

    :param oli_name: string name of a solute in OLI format

    :return watertap_name: string name of a solute in WaterTAP format
    """

    if oli_name in names_db:
        return names_db[oli_name]
    else:
        raise IOError(
            f" Component {oli_name} not found in names_db."
            + " Update this dictionary to hard code additional OLI names."
        )


"""
Here follows a dictionary of OLI names and their WaterTAP counterparts.

It functions to aid reverse lookup, i.e., if a name is already in OLI format,
the necessary data can still be extracted.

TODO: method to add novel (valid) names to names_db
"""

names_db = {
    "NAION": "Na_+",
    "CLION": "Cl_-",
    "SO4ION": "SO4_2-",
    "MGION": "Mg_2+",
    "CAION": "Ca_2+",
    "KION": "K_+",
    "HCO3ION": "HCO3_-",
    "NA2CO3": "Na2CO3",
    "CO2": "CO2",
    "H2O": "H2O",
}
