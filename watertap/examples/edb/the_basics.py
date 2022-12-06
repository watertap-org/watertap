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
    This file demonstrates the basics of working with and using the electrolyte
    database (EDB).

    (1) Before we can start, you must install MongoDB (which is installed separately)

        [See more information on the ReadTheDocs under 'Getting Started --> Installing WaterTAP']


    (2) After installing MongoDB, you will need to 'load' the database using the
        command line function 'edb load -b'. This will load the default database
        that WaterTAP is bundled with.

        [NOTE: If you need to 'reload' the database, simply use the command 'edb drop -d electrolytedb'
                in the command line. The database on MongoDB is named "electrolytedb"]

        [NOTE 2: You can invoke the command line utility with the "help" keyword to
                get more information on funtionality. Command: 'edb --help' or 'edb [arg] --help']


    (3) To use EDB in python, start by importing the interface class object 'ElectrolyteDB'


    (4) Invoke the 'ElectrolyteDB' object to connect to the database


    (5) Grab a 'base' for a configuration dictionary, and place it into a class object


    (6) Get the chemcial species/components for a simulation case. There are a number of ways
        to do this. In this example, we will grab them by finding all components that contain
        only specific elements. Then, we add those components and their associated parameters
        to the configuration dictionary being built from the 'base'.

        [NOTE: An alternative method is to provide a list of the names of components you want]


    (7) Get the set of reactions you want in your system and put into a 'base' object.
        That 'base' can be either a 'thermo' base or a 'reaction' (as in this case)
        base. IF you are adding reactions to a 'thermo' base, they should be added
        as 'inherent' reactions. IF you are adding reactions to a 'reaction' base,
        they should be added as 'equilibrium' (or other) reactions.


    (8) When using an reactor object in IDAES, you must always provide a 'reaction_config'
        to match with the 'thermo_config'. We can create a base 'reaction' config from
        the database and add reactions to that config in the same way we do for the
        'thermo_config' when adding reactions as inherent.

        [NOTE: If a reaction is added to a 'thermo_config' as 'inherent', this it should
               NOT be added to a 'reaction_config' as 'equilibrium']

"""

# ========= These imports (below) are for testing the configs from EDB ===============
# Import specific pyomo objects
from pyomo.environ import (
    ConcreteModel,
)

# Import the idaes objects for Generic Properties and Reactions
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models.properties.modular_properties.base.generic_reaction import (
    GenericReactionParameterBlock,
)

# Import the idaes object for the EquilibriumReactor unit model
from idaes.models.unit_models.equilibrium_reactor import EquilibriumReactor

# Import the core idaes objects for Flowsheets and types of balances
from idaes.core import FlowsheetBlock

# ========= These imports (above) are for testing the configs from EDB ===============


# ========================== (3) ================================
# Import ElectrolyteDB object
from watertap.edb import ElectrolyteDB

__author__ = "Austin Ladshaw"

# ========================== (4) ================================
# By default, invoking the 'ElectrolyteDB' object (with no args)
#   will attempt to connect to the local host database. You can
#   check the connection by calling the 'can_connect' function
#   and passing the 'host' and 'port' as args. If no 'host' or
#   'port' are given, then it uses the defaults.
def connect_to_edb(test_invalid_host=False):
    print("connecting to " + str(ElectrolyteDB.DEFAULT_URL))
    db = ElectrolyteDB()
    connected = db.can_connect()
    return (db, connected)


# ========================== (5) ================================
# All configuration files used in WaterTAP for electrolyte chemistry
#   require a 'base' dictionary to start. For example, we need to
#   create a 'thermo_config' dictionary to pass to the GenericProperty
#   package in IDAES. That 'thermo_config' file will always have a
#   few specific items in common with most other configuration files.
#   Thus, this operation will populate a class object that holds the
#   data assocated with that 'base' dictionary.
#
# In the EDB, there are several different 'base' structures to start
#   from. In this example, we will build from the 'default_thermo'
#   configuration base.
def grab_base_thermo_config(db):
    # Get the base and place into a result object
    base = db.get_base("default_thermo")
    return base


# ========================== (6) ================================
# Get chemical components/species for a simulation case
#       NOTE: This function here also returns a 'list' of the
#           components that it finds. This is not a built in
#           feature of the EDB, but is very useful because
#           getting reactions is dependent on the component list.
def get_components_and_add_to_idaes_config(db, base_obj, by_elements=False):
    # Going to grab all components that contain ONLY "H" and "O"
    #   Expected behavior = Will get "H2O", "H_+", and "OH_-"
    element_list = ["H", "O"]

    # Alternatively, you can pass a list of individual componets
    #   you want to grab and the EDB functions should grab explicitly
    #   those components/species you want.
    comp_list = ["H2O", "H_+", "OH_-"]

    # Just like before, this function returns a results object
    #   that contains other objects that must be iterated through
    #   in order to access the information. Then, call the 'add'
    #   function to add those components to the 'base' object
    if by_elements == True:
        res_obj_comps = db.get_components(element_names=element_list)
    else:
        res_obj_comps = db.get_components(component_names=comp_list)

    # Iterate through the results object and add the components
    #   to the base_obj
    db_comp_list = []
    for comp_obj in res_obj_comps:
        print("Adding " + str(comp_obj.name) + "")
        base_obj.add(comp_obj)
        db_comp_list.append(comp_obj.name)
    print()
    return (base_obj, db_comp_list)


# ========================== (7) ================================
# Grab the reactions associated with the list of components and add
#   them to a base object (which could be a 'thermo' base or 'reaction' base)
#
def get_reactions_return_object(db, base_obj, comp_list, is_inherent=True):
    react_obj = db.get_reactions(component_names=comp_list)
    for r in react_obj:
        print("Found reaction: " + str(r.name))
        if is_inherent == True:
            r._data["type"] = "inherent"
        base_obj.add(r)
    return base_obj


# ========================== (8) ================================
# Create a base config for reactions.
def grab_base_reaction_config(db):
    # Get the base and place into a result object
    base = db.get_base("reaction")
    return base


# This function will produce an error if the thermo config is not correct
def is_thermo_config_valid(thermo_config):
    model = ConcreteModel()
    model.fs = FlowsheetBlock(dynamic=False)
    model.fs.thermo_params = GenericParameterBlock(**thermo_config)
    return True


# This function will produce an error if the thermo config is not correct
#   or if the pairing of the thermo and reaction config are invalid
def is_thermo_reaction_pair_valid(thermo_config, reaction_config):
    model = ConcreteModel()
    model.fs = FlowsheetBlock(dynamic=False)
    model.fs.thermo_params = GenericParameterBlock(**thermo_config)
    model.fs.rxn_params = GenericReactionParameterBlock(
        property_package=model.fs.thermo_params, **reaction_config
    )
    return True


# Run script for testing
def run_the_basics_with_mockdb(db):
    base_obj = grab_base_thermo_config(db)

    (base_obj, comp_list) = get_components_and_add_to_idaes_config(db, base_obj)

    # Create a reaction config
    react_base = grab_base_reaction_config(db)

    # Add reactions to the reaction base as 'equilibrium'
    react_base = get_reactions_return_object(
        db, react_base, comp_list, is_inherent=False
    )

    # If all goes well, this function returns true
    return is_thermo_reaction_pair_valid(base_obj.idaes_config, react_base.idaes_config)


# Run script for testing
def run_the_basics_alt_with_mockdb(db):
    base_obj = grab_base_thermo_config(db)

    (base_obj, comp_list) = get_components_and_add_to_idaes_config(
        db, base_obj, by_elements=True
    )

    # Add reactions to the thermo base as 'inherent'
    base_obj = get_reactions_return_object(db, base_obj, comp_list, is_inherent=True)

    # If all goes well, this function returns true
    return is_thermo_config_valid(base_obj.idaes_config)


# Run script for testing
def run_the_basics_dummy_rxn_with_mockdb(db):
    base_obj = grab_base_thermo_config(db)

    (base_obj, comp_list) = get_components_and_add_to_idaes_config(
        db, base_obj, by_elements=True
    )

    # Add reactions to the thermo base as 'inherent'
    base_obj = get_reactions_return_object(db, base_obj, comp_list, is_inherent=True)

    # Create a reaction config
    react_base = grab_base_reaction_config(db)

    # If no reactions are in the reaction base, this will cause an error in IDAES.
    #   However, we can add a 'dummy' reaction just to satisfy the IDAES code base.
    react_obj = db.get_reactions(reaction_names=["dummy"])
    for r in react_obj:
        print("Found reaction: " + str(r.name))
        react_base.add(r)

    # IDAES will throw an exception when we try to do this if something is wrong
    thermo_config = base_obj.idaes_config
    reaction_config = react_base.idaes_config
    model = ConcreteModel()
    model.fs = FlowsheetBlock(dynamic=False)
    model.fs.thermo_params = GenericParameterBlock(**thermo_config)
    model.fs.rxn_params = GenericReactionParameterBlock(
        property_package=model.fs.thermo_params, **reaction_config
    )

    model.fs.unit = EquilibriumReactor(
        property_package=model.fs.thermo_params,
        reaction_package=model.fs.rxn_params,
        has_rate_reactions=False,
        has_equilibrium_reactions=False,
        has_heat_transfer=False,
        has_heat_of_reaction=False,
        has_pressure_change=False,
    )

    # If all goes well, this function returns true
    return is_thermo_reaction_pair_valid(base_obj.idaes_config, react_base.idaes_config)
