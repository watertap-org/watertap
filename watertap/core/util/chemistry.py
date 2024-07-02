from re import findall
from pathlib import Path

import pandas as pd
from pyomo.environ import units as pyunits


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
        try:
            charge_sign = charge[-1]
        except IndexError:
            raise IOError(
                f"Charge sign could not be determined from the string '{watertap_name}'"
            )
        if len(charge) > 1:
            try:
                charge_magnitude = int(charge[:-1])
            except ValueError:
                raise IOError(
                    f"Charge sign could not be determined from the string '{watertap_name}'"
                )
        else:
            charge_magnitude = 1
        if charge_sign == "+":
            charge = charge_magnitude
        elif charge_sign == "-":
            charge = -charge_magnitude
        else:
            raise IOError(
                f"Only + and - are valid charge indicators and neither was provided in '{watertap_name}'."
            )
    else:
        raise IOError(
            f"Charge could not be determined from the string '{watertap_name}'"
        )
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

    periodic_table = get_periodic_table()

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
        try:
            atomic_mass = float(
                periodic_table["AtomicMass"][
                    (periodic_table["Symbol"] == element)
                ].values[0]
            )
        except IndexError:
            raise IOError(
                f"The symbol '{element}' from the component name '{components[0]}' could not be found in the periodic table."
            )

        molar_mass += element_counts[element] * atomic_mass

    if not molar_mass:
        raise IOError(f"Molecular weight data could not be found for {watertap_name}.")

    return molar_mass


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


def get_periodic_table() -> pd.DataFrame:
    parent_dir = Path(__file__).parent
    return pd.read_csv(parent_dir / "periodic_table.csv")


def get_molar_mass_quantity(watertap_name: str, units=pyunits.kg / pyunits.mol):
    """
    Extracts atomic weight data from a periodic table file
    to generate the molar mass of a chemical substance in pint units.
    Since get_molar_mass returns only the value, which has inherent units of g/mol,
    this function converts to kg/mol by default, the units used for molecular weight by convention in WaterTAP.

    :param watertap_name: string name of a solute in WaterTAP format

    :return desired_quantity: molar mass of solute in pint units. Conversion from g/mol to kg/mol by default.
    """
    molar_mass_value = get_molar_mass(watertap_name)
    inherent_quantity = molar_mass_value * pyunits.g / pyunits.mol
    desired_quantity = pyunits.convert(inherent_quantity, to_units=units)
    return desired_quantity
