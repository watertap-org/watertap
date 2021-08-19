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
This test is to establish that the core chemistry packages in IDAES solve
a simple water dissociation problem and return the correct pH value.

Modified to use the Electrolyte Database -dang 08/2021
"""
import pytest
from proteuslib.edb import ElectrolyteDB
from .test_pure_water_pH import TestPureWater

g_edb = None
try:
    raise RuntimeError()
    g_edb = ElectrolyteDB()
except Exception:
    pass


def get_thermo_config(edb):
    base = edb.get_one_base("water_reaction")
    components = ["H2O", "H_+", "OH_-"]
    # Add the components
    for c in edb.get_components(components):
        # Need to remove these to avoid errors when using the generated config
        c.remove("valid_phase_types")
        c.remove("enth_mol_ig_comp")
        c.remove("phase_equilibrium_form")
        c.remove("pressure_sat_comp")
        base.add(c)
    # Add the reactions
    for r in edb.get_reactions(component_names=components):
        r.set_reaction_order('Liq', ('H2O',), ('H_+', 'OH_-'))
        r._data["type"] = "inherent"
        base.add(r)
    return base.idaes_config


def get_water_reaction_config(edb):
    components = ["H2O", "H_+", "OH_-"]
    base = edb.get_one_base("water_reaction")
    # Need to remove these to avoid errors when using the generated config
    base.remove("phases")
    base.remove("pressure_ref")
    base.remove("state_bounds")
    base.remove("state_definition")
    base.remove("temperature_ref")
    # Add the reactions
    for r in edb.get_reactions(component_names=components):
        # Set a custom reaction order
        r.set_reaction_order('Liq', ('H2O',), ('H_+', 'OH_-'))
        # Need to remove this to avoid errors when using the generated config
        r.remove_parameter("ds_rxn_ref")
        base.add(r)
    return base.idaes_config


@pytest.mark.skipif(g_edb is None, reason="Cannot connect to MongoDB")
class TestPureWaterEdb(TestPureWater):
    """Run all tests in TestPureWater, but with different configs.
    """
    if g_edb:
        thermo_config = get_thermo_config(g_edb)
        wr_config = get_water_reaction_config(g_edb)
        thermo_only_config = {}
        thermo_only_config.update(thermo_config)
        del thermo_only_config["inherent_reactions"]
