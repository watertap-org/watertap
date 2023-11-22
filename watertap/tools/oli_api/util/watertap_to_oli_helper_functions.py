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
#
#################################################################################

"""
This file contains methods to convert WaterTAP naming conventions to OLI
and generate molecular weight and charge dictionaries from molecular formulae. 

It calculates molecular weights using the periodic_table.csv from:
https://gist.github.com/GoodmanSciences/c2dd862cd38f21b0ad36b8f96b4bf1ee.
"""

__author__ = "Paul Vecchiarelli, Ben Knueven"

from collections import namedtuple
from re import findall
from pathlib import Path
from os.path import join
from pandas import read_csv

# TODO: maybe replace some functionality with molmass: https://pypi.org/project/molmass

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


# TODO: merge with other helper functions in this file
def get_oli_name(watertap_name: str) -> str:
    """
    Converts an WaterTAP formatted name, i.e., "Na_+"
    into an OLI formatted name, i.e., "NAION".

    :param watertap_name: string name of a solute in WaterTAP format

    :return oli_name: string name of a solute in OLI format
    """

    components = watertap_name.split("_")
    if len(components) == 0:
        raise IOError(f" Unable to parse solute '{watertap_name}'.")
    if len(components) == 1:
        molecule = components[0]
    elif len(components) == 2:
        molecule = components[0] + "ION"
    oli_name = molecule.replace("[", "").replace("]", "").upper()
    return oli_name


def get_charge(watertap_name: str) -> int:
    """
    Gets charge from WaterTAP formatted names.

    :param watertap_name: string name of a solute in WaterTAP format

    :return charge: integer value of charge
    """

    components = watertap_name.split("_")
    if len(components) == 0:
        raise IOError(f" Unable to parse solute '{watertap_name}'.")
    if len(components) == 1:
        molecule = components[0]
        charge = 0
    elif len(components) == 2:
        molecule = components[0] + "ION"
        charge = components[1]
        charge_sign = charge[-1]
        if len(charge) > 1:
            charge_magnitude = int(charge[:-1])
        else:
            charge_magnitude = 1
        if charge_sign == "+":
            charge = charge_magnitude
        elif charge_sign == "-":
            charge = -charge_magnitude
        else:
            raise IOError(" Only + and - are valid charge indicators.")
    return charge


def get_charge_group(charge: int) -> str:
    """
    Categorizes molecule based on its charge.

    :param charge: integer value for charge

    :return group: string name for charge group
    """

    if charge == 0:
        group = "Neutrals"
    elif charge > 0:
        group = "Cations"
    elif charge < 0:
        group = "Anions"
    return group


def get_molar_mass(watertap_name: str) -> float:
    """
    Extracts atomic weight data from a periodic table file
    to generate the molar mass of a chemical substance.

    TODO: additional testing for complex solutes
    such as CH3CO2H, [UO2]2[OH]4, etc.

    :param watertap_name: string name of a solute in WaterTAP format

    :return molar_mass: float value for molar mass of solute
    """

    file_path = Path(__file__).parents[0]
    periodic_table = read_csv(join(file_path, "periodic_table.csv"))

    components = watertap_name.split("_")
    elements = findall("[A-Z][a-z]?[0-9]*", components[0])
    element_counts = {}
    for element in elements:
        if len(element) == 1:
            element_counts[element] = 1
        elif len(element) == 2 and element.isalpha():
            element_counts[element] = 1
        elif len(element) == 2 and not element.isalpha():
            element_counts[element[:-1]] = int(element[-1])
        elif len(element) == 3 and element[:-1].isalpha():
            element_counts[element[:-1]] = int(element[-1])
        elif len(element) == 3 and not element[:-1].isalpha():
            element_counts[element[:-2]] = int(element[-2:-1])
        else:
            raise IOError(f" Too many characters in {element}.")

        element_location = components[0].find(element)
        if "[" in components[0]:
            boundary = (components[0].find("["), components[0].find("]"))
            coefficient = int(components[0][boundary[1] + 1])
            if element_location > boundary[0] and element_location < boundary[1]:
                element_counts[element] *= coefficient

    molar_mass = 0
    for element in element_counts:
        atomic_mass = float(
            periodic_table["AtomicMass"][(periodic_table["Symbol"] == element)].values[
                0
            ]
        )
        molar_mass += element_counts[element] * atomic_mass
    return molar_mass


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
