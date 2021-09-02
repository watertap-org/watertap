###############################################################################
# ProteusLib Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/nawi-hub/proteuslib/"
#
###############################################################################

"""
This is the EDB version of test_enrtl_water_pH.py

See that module for details on what the tests accomplish
in terms of the IDAES chemistry packages.

XXX: This test currently doesn't *do* anything
"""
import pytest
from proteuslib.edb import ElectrolyteDB

# Set global database object after checking that MongoDB server is up
g_edb = None
if ElectrolyteDB.can_connect():
    g_edb = ElectrolyteDB()


def get_thermo_config(edb):
    from idaes.generic_models.properties.core.eos.enrtl import ENRTL
    from idaes.generic_models.properties.core.eos.enrtl_reference_states import (
        Unsymmetric,
    )

    base = edb.get_one_base("water_reaction")
    elements = ["H", "O"]
    components = []
    # Add the components
    for c in edb.get_components(element_names=elements):
        # Need to remove these to avoid errors when using the generated config
        c.remove("valid_phase_types")
        c.remove("enth_mol_ig_comp")
        c.remove("phase_equilibrium_form")
        c.remove("pressure_sat_comp")
        # Add specific stuff not in the DB
        if c.name == "H2O":
            c.data["relative_permittivity_liq_comp"] = "relative_permittivity_constant"
            c.set_parameter("relative_permittivity_liq_comp", 78.54)
        base.add(c)
        components.append(c.name)

    # Add the reaction
    r = list(edb.get_reactions(reaction_names=["H2O_Kw_2"]))[0]
    r.set_reaction_order("Liq", ("H2O",), ("H_+", "OH_-"))
    base.add(r)

    cfg = base.idaes_config.copy()
    cfg["phases"]["Liq"]["equation_of_state"] = ENRTL
    cfg["phases"]["Liq"]["equation_of_state_options"] = {"reference_state": Unsymmetric}
    cfg["inherent_reactions"] = cfg["equilibrium_reactions"]
    del cfg["equilibrium_reactions"]
    return cfg


def get_water_reaction_config(edb):
    from idaes.generic_models.properties.core.reactions.equilibrium_forms import (
        log_power_law_equil,
    )

    elements = ["H", "O"]
    components = [c.name for c in edb.get_components(element_names=elements)]
    base = edb.get_one_base("water_reaction")
    # Just use base units, and add a dummy equilibrium reaction config
    cfg = {
        "base_units": base.idaes_config["base_units"],
        "equilibrium_reactions": {
            "dummy": {
                "stoichiometry": {},
                "equilibrium_form": log_power_law_equil,
            }
        },
    }
    return cfg


# ====================================================================================


def get_carbonic_thermo_config(edb):
    from idaes.generic_models.properties.core.eos.enrtl import ENRTL
    from idaes.generic_models.properties.core.eos.enrtl_reference_states import (
        Unsymmetric,
    )

    base = edb.get_one_base("water_reaction")
    # The 'right' way to fetch all components
    # elements = ["H", "O", "C", "Na"]
    components = []
    # Add the components
    component_names = ("Na_+", "CO3_2-",  "HCO3_-", "H2CO3", "OH_-", "H_+", "H2O")
    for c in edb.get_components(component_names=component_names):
        # Need to remove these to avoid errors when using the generated config
        if c.name != "H2CO3":
            c.remove("valid_phase_types")
        c.remove("enth_mol_ig_comp")
        c.remove("phase_equilibrium_form")
        c.remove("pressure_sat_comp")
        # Add specific stuff not in the DB
        if c.name == "H2O":
            c.data["relative_permittivity_liq_comp"] = "relative_permittivity_constant"
            c.set_parameter("relative_permittivity_liq_comp", 78.54)
        base.add(c)
        components.append(c.name)

    # Add the reaction
    for r in edb.get_reactions(reaction_names=["H2O_Kw_2", "H2CO3_Ka1", "H2CO3_Ka2"]):
        if r.name == "H2O_Kw_2":
            r.set_reaction_order("Liq", ("H2O",), ("H_+", "OH_-"))
        elif r.name == "H2CO3_Ka1":
            r.set_reaction_order("Liq", ("H2CO3",), ("H_+", "HCO3_-"), lhs_value=-1)
            r.remove_parameter("ds_rxn_ref")
            r.set_parameter("dh_rxn_ref", 0, units="kJ/mol")
        elif r.name == "H2CO3_Ka2":
            r.set_reaction_order("Liq", ("HCO3_-",), ("H_+", "CO3_2-"), lhs_value=-1)
            r.remove_parameter("ds_rxn_ref")
            r.set_parameter("dh_rxn_ref", 0, units="kJ/mol")
        base.add(r)

    cfg = base.idaes_config.copy()
    cfg["phases"]["Liq"]["equation_of_state"] = ENRTL
    cfg["phases"]["Liq"]["equation_of_state_options"] = {"reference_state": Unsymmetric}
    cfg["inherent_reactions"] = cfg["equilibrium_reactions"]
    del cfg["equilibrium_reactions"]
    return cfg

