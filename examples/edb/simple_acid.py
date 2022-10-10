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
    This file demonstrates how to use EDB to create a simple acid problem.

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
        This time, we will grab a base that is for a Liq only problem using FpcTP
        state variables.


    (6) Get the chemcial species/components for a simulation case. There are a number of ways
        to do this. In this example, we will grab them by finding all components that contain
        only specific elements. Then, we add those components and their associated parameters
        to the configuration dictionary being built from the 'base'.

        [NOTE: An alternative method is to provide a list of the names of components you want]


    (7) Get the set of reactions you want in your system and put into a 'base' object.
        In this case, we are getting all reactions associated with a system of water
        and carbonic acid. We should get three reactions::

            H2O <--> H_+ + OH_-
            H2CO3 <--> H_+ + HCO3_-
            HCO3_- <--> H_+ + CO3_2-


    (8) When using an reactor object in IDAES, you must always provide a 'reaction_config'
        to match with the 'thermo_config'. We can create a base 'reaction' config from
        the database and add reactions to that config in the same way we do for the
        'thermo_config' when adding reactions as inherent.

        [NOTE: If a reaction is added to a 'thermo_config' as 'inherent', it should
               NOT be added to a 'reaction_config' as 'equilibrium']


    (9) [NEW Step] Build an equilibrium reactor from the 'thermo_config' and 'reaction_config'
        that were generated from the EDB.

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


# ========================== (3 & 4) ================================
# Import ElectrolyteDB object
from watertap.edb import ElectrolyteDB
from watertap.examples.edb.the_basics import (
    connect_to_edb,
    is_thermo_config_valid,
    grab_base_reaction_config,
    is_thermo_reaction_pair_valid,
)

__author__ = "Austin Ladshaw"

# ========================== (5) ================================
# Grab a new base config for our thermo, but this time we will use
#   one of the newer bases that will use the FpcTP state vars and
#   a Liq only system.
def grab_thermo_Liq_FpcTP_base(db):
    # Get the base and place into a result object
    base = db.get_base("thermo_Liq_FpcTP")
    return base


# ========================== (6) ================================
# Get chemical components/species for a simulation case
#       NOTE: This function here also returns a 'list' of the
#           components that it finds. This is not a built in
#           feature of the EDB, but is very useful because
#           getting reactions is dependent on the component list.
def get_components_and_add_to_idaes_config(db, base_obj, comp_list):
    res_obj_comps = db.get_components(component_names=comp_list)

    # Iterate through the results object and add the components
    #   to the base_obj
    for comp_obj in res_obj_comps:
        print("Adding " + str(comp_obj.name) + "")
        base_obj.add(comp_obj)
    print()
    return base_obj


# ========================== (7) ================================
# Grab the reactions associated with the list of components and add
#   them to a reaction base as equilibrium reactions
#
def add_equilibrium_reactions_to_react_base(db, react_base_obj, comp_list):
    react_obj = db.get_reactions(component_names=comp_list)
    for r in react_obj:
        print("Found reaction: " + str(r.name))
        react_base_obj.add(r)
    return react_base_obj


# ========================== (9) ================================
# Create the Pyomo model by using the thermo_config and reaction_config
#   that were generated from the EDB.
#
def build_equilibrium_model(thermo_config, reaction_config):
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
        has_equilibrium_reactions=True,
        has_heat_transfer=False,
        has_heat_of_reaction=False,
        has_pressure_change=False,
    )

    return model


# Run script for testing
def run_simple_acid_with_mockdb(db):
    base_obj = grab_thermo_Liq_FpcTP_base(db)

    # Our components for this problem are as follows:
    comp_list = ["H2O", "H_+", "OH_-", "H2CO3", "HCO3_-", "CO3_2-"]

    base_obj = get_components_and_add_to_idaes_config(db, base_obj, comp_list)

    # Create a reaction config
    react_base = grab_base_reaction_config(db)

    # Add reactions to the reaction base as 'equilibrium'
    react_base = add_equilibrium_reactions_to_react_base(db, react_base, comp_list)

    # Now, we can actually see if we created a correct model by looking
    #       for degrees of freedom, state variables, etc.
    thermo_config = base_obj.idaes_config
    reaction_config = react_base.idaes_config
    model = build_equilibrium_model(thermo_config, reaction_config)

    return model
