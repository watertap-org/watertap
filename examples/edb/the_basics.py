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

        [# WARNING: A feature of 'ElectrolyteDB' appears to be missing from WaterTAP main
        that is preventing us from following the tutorial verbatim.]


    (5) Grab a 'base' for a configuration dictionary, and place it into a class object

        [# WARNING: There are a number of ways to do this. Behavior is inconsistent with
        some of the tutorials currently.]


    (6) Get the chemcial species/components for a simulation case. There are a number of ways
        to do this. In this example, we will grab them by finding all components that contain
        only specific elements. Then, we add those components and their associated parameters
        to the configuration dictionary being built from the 'base'.

        [NOTE: An alternative method is to provide a list of the names of components you want]


    (7-a) Get the set of reactions you want in your system and put into a 'base' object.
        That 'base' can be either a 'thermo' base (as in this example) or a 'reaction'
        base. IF you are adding reactions to a 'thermo' base, they should be added
        as 'inherent' reactions. IF you are adding reactions to a 'reaction' base,
        they should be added as 'equilibrium' (or other) reactions. 

"""

# ========================== (3) ================================
# Import ElectrolyteDB object
from watertap.edb import ElectrolyteDB

# ========================== (4) ================================
# By default, invoking the 'ElectrolyteDB' object (with no args)
#   will attempt to connect to the local host database. You can
#   check the connection by calling the 'can_connect' function
#   and passing the 'host' and 'port' as args. If no 'host' or
#   'port' are given, then it uses the defaults.
def connect_to_edb():
    print("connecting to " + str(ElectrolyteDB.DEFAULT_URL))
    db = ElectrolyteDB()
    connected = db.can_connect()

    # Supposedly, we should be able to pass a check_connection arg to the
    #   initialization of the 'ElectrolyteDB' object. However, that does
    #   not appear to be a feature that exists in WaterTAP main.
    try:
        db2 = ElectrolyteDB("mongodb://some.other.host", check_connection=False)
    except:
        print("\tWARNING: If you are here, that means that the EDB code is not up-to-date")

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
#   from. In this example, we will build from the default 'thermo'
#   configuration base.
#
#   NOTE: There is a difference between the 'get_base' function and
#       the 'get_one_base' function. The 'get_base' function will
#       return a 'Result' object that contains an iterator for each
#       base requested. This means that in order to access the config
#       you just created, you would need to iterate through that
#       result object and look at the 'idaes_config' for each obj
#       in the 'result'. Alternatively, you can use the 'get_one_base'
#       function, which will directly return a single object that
#       can be directly queried for the 'idaes_config'.
def grab_base_thermo_config(db):

    # Get the base and place into a result object that
    #   needs to be iterated through to access data at
    #   each object held in the result
    res_it_obj = db.get_base("thermo")
    for obj in res_it_obj:
        obj.idaes_config  # This is where you access the data

    # Get the base and place into a singular object
    #   with no need to iterate through
    base_obj = db.get_one_base("thermo")
    base_obj.idaes_config  # This gives more direct access to the information
    return base_obj

# ========================== (6) ================================
# Get chemical components/species for a simulation case
#       NOTE: This function here also returns a 'list' of the
#           components that it finds. This is not a built in
#           feature of the EDB, but is very useful because
#           getting reactions is dependent on the component list.
def get_components_and_add_to_idaes_config(db, base_obj, by_elements=True):
    # Going to grab all components that contain ONLY "H" and "O"
    #   Expected behavior = Will get "H2O", "H_+", and "OH_-"
    element_list = ["H","O"]

    # Alternatively, you can pass a list of individual componets
    #   you want to grab and the EDB functions should grab explicitly
    #   those components/species you want.
    comp_list = ["H2O","H_+","OH_-"]

    # Just like before, this function returns a results object
    #   that contains other objects that must be iterated through
    #   in order to access the information. Then, call the 'add'
    #   function to add those components to the 'base' object
    if (by_elements==True):
        res_obj_comps = db.get_components(element_names=element_list)
    else:
        res_obj_comps = db.get_components(component_names=comp_list)

    # Iterate through the results object and add the components
    #   to the base_obj
    db_comp_list = []
    for comp_obj in res_obj_comps:
        print("Adding " + str(comp_obj.name) + "" )
        base_obj.add(comp_obj)
        db_comp_list.append(comp_obj.name)
    print()
    return (base_obj, db_comp_list)

# ========================== (7-a) ================================
# Grab the reactions associated with the list of components and add
#   them to a base object (which could be a 'thermo' base or 'reaction' base)
#
def get_reactions_return_object(db, base_obj, comp_list, is_inherent=True):
    react_obj = db.get_reactions(component_names=comp_list)
    for r in react_obj:
        print("Found reaction: " + str(r.name))
        if (is_inherent == True):
            r._data["type"] = "inherent"
        else:
            r._data["type"] = "equilibrium"
        base_obj.add(r)
    return base_obj

if __name__ == "__main__":
    (db, connected) = connect_to_edb()
    if (connected == False):
        print("\nFailed to connect!!!\n")
        exit()
    else:
        print("\n...connected\n")

    base_obj = grab_base_thermo_config(db)
    try:
        base_obj.idaes_config
    except:
        print("Error! Object does NOT contain 'idaes_config' dict!")
        exit()

    (base_obj, comp_list) = get_components_and_add_to_idaes_config(db, base_obj)

    base_obj = get_reactions_return_object(db, base_obj, comp_list)

    print(base_obj.idaes_config)

    # At this point, we should have a valid dict...
