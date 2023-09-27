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

"""
This file contains methods to convert WaterTAP naming conventions to OLI
and generate molecular weight and charge dictionaries from molecular formulae. 

It calculates molecular weights using the periodic_table.csv from:
https://gist.github.com/GoodmanSciences/c2dd862cd38f21b0ad36b8f96b4bf1ee.
"""

__author__ = "Paul Vecchiarelli"

from collections import namedtuple
from pandas import read_csv
from re import findall

OLIName = namedtuple("OLIName", ["oli_name", "watertap_name", "charge", "charge_group", "molar_mass"])

def watertap_to_oli(watertap_name: str) -> OLIName:
    """
    This method creates a named tuple
    which can be passed directly into OLI
    or into MCAS property models.
    
    Parameters
    ----------
    watertap_name : str
        a substance name in WaterTAP format (typical IDAES convention, i.e., B[OH]4_-)

    Raises
    ------
    RuntimeError: this generic error is raised if:
        - watertap_name is a blank string
        - an extracted charge_sign is not plus or minus
        - if more than 1 "_" is present in watertap_name
        
    Returns
    -------
    OLIName: named tuple containing:
        - oli_name: substance name converted to OLI convention
        - watertap_name: original name passed into method
        - charge: charge value
        - charge_group: categorization under one of: Cations, Anions, Neutrals
        - molar_mass: weight of 1 mole of the substance (grams per mole)

    """
    components = watertap_name.split("_")
    
    # neutral molecule
    if len(components) == 1:
        if components[0] != '':
            molecule = components[0]
            charge = 0
        else:
            raise RuntimeError("Unrecognized component " + watertap_name )
    
    # charged molecule
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
            raise RuntimeError("Unrecognized charge magnitude " + charge_magnitude )
        
    else:
        raise RuntimeError("Unrecognized name format " + watertap_name )
        
    charge_group = _charge_group_computation(charge)
        
    oli_name = molecule.replace("[", "").replace("]", "").upper()    
    molar_mass = _molar_mass_from_watertap_name(components[0])
    
    return OLIName(oli_name, watertap_name, charge, charge_group, molar_mass)

def _molar_mass_from_watertap_name(watertap_formula: str) -> float:
    """
    This function extracts atomic weight data from a periodic table file
    to generate the molar mass of a chemical substance.
    
    This function requires additional testing for complex solutes
    such as CH3CO2H, [UO2]2[OH]4, etc.
    """
    
    # change to local file location
    periodic_table = read_csv("watertap/tools/oli/periodic_table.csv")
    
    # isolate single- and double- letter elements and add to element_counts dict
    elements = findall('[A-Z][a-z]?[0-9]*', watertap_formula)
    
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
            raise RuntimeError("Too many characters in " + element)
        
        element_location = watertap_formula.find(element)
        
        # find brackets and following coefficient
        if "[" in watertap_formula:
            boundary = (watertap_formula.find("["), watertap_formula.find("]"))
            coefficient = int(watertap_formula[boundary[1]+1])
            if element_location > boundary[0] and element_location < boundary[1]:
                element_counts[element] *= coefficient
                    
    # compute molecular weight of compound
    molar_mass = 0
    for element in element_counts:
        atomic_mass = float(periodic_table["AtomicMass"][(periodic_table["Symbol"] == element)].values[0])
        molar_mass += element_counts[element]*atomic_mass
    return molar_mass

def _charge_group_computation(charge: int) -> str:
    if charge == 0:
        group = "Neutrals"
    elif charge > 0:
        group = "Cations"
    elif charge < 0:
        group = "Anions"
    return group

def oli_reverse_lookup(oli_name: str, names_db) -> OLIName:
    if oli_name in names_db:
        return watertap_to_oli(names_db[oli_name])
    else:
        raise RuntimeError(f" Component {oli_name} not found in names_db." +
                           " Please update the dictionary if you wish to hard code OLI names.")
    
"""
Here follows a dictionary of OLI names and their WaterTAP counterparts.

It functions to aid reverse lookup, i.e., if a name is already in OLI format,
the necessary data can still be extracted.

TODO: method to add novel (valid) names to names_db
"""

names_db = {"NAION": "Na_+",
            "CLION": "Cl_-",
            "SO4ION": "SO4_2-",
            "MGION": "Mg_2+",
            "CAION": "Ca_2+",
            "KION": "K_+",
            "HCO3ION": "HCO3_-",
            "NA2CO3": "Na2CO3",
            "CO2": "CO2",
            "H2O": "H2O"}