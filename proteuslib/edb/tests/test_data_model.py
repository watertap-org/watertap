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
Tests for data_model module
"""
import logging
from pprint import pprint  # for debugging
import pytest
from . import data as testdata
from typing import Dict
from pyomo.environ import units as pyunits

# Used to test IDAES config -> DataWrapper
from idaes.generic_models.properties.core.pure import Perrys
from idaes.generic_models.properties.core.pure.NIST import NIST
from idaes.generic_models.properties.core.phase_equil.forms import fugacity
from idaes.generic_models.properties.core.reactions.equilibrium_forms import log_power_law_equil
from idaes.core import Component as IComponent
from idaes.generic_models.properties.core.reactions.equilibrium_constant import (
    van_t_hoff,
)
from idaes.generic_models.properties.core.reactions.dh_rxn import constant_dh_rxn
from idaes.generic_models.properties.core.generic.generic_reaction import (
    ConcentrationForm,
)
from idaes.core.components import Solvent, Solute, Cation, Anion


from ..data_model import (
    ConfigGenerator,
    Component,
    Reaction,
    Result,
    Base,
    ThermoConfig,
)

# For validating DataWrapper contents
from ..validate import validate

# errors
from ..error import BadConfiguration, ValidationError


@pytest.fixture
def debug_logging():
    """Use fixture cleanup ability to set and unset debug logging."""
    logging.getLogger("proteuslib.edb.data_model").setLevel(logging.DEBUG)
    yield "debug"
    logging.getLogger("proteuslib.edb.data_model").setLevel(logging.INFO)


def assert_configuration_equal(a: Dict, b: Dict, fld: str):
    """Walk through and compare things in two config dicts.

    This is needed because simple comparison will also compare objects that may not have appropriate
    equality methods defined and therefore fail for the simple reason that they have different object ids.
    """
    a_section, b_section = a.get(fld, {}), b.get(fld, {})
    assert len(a_section) == len(b_section)

    for name in a_section:
        assert name in b_section
        a_data, b_data = a_section[name], b_section[name]
        for key in a_data:
            assert key in b_data
            if key == "parameter_data":
                for key2, value2 in a_data[key].items():
                    assert key2 in b_data[key]
                    if (
                        isinstance(value2, tuple) and len(value2) == 2
                    ):  # number, unit pair
                        assert b_data[key][key2][0] == pytest.approx(value2[0])


@pytest.mark.unit
def test_config_generator():
    hc = ConfigGenerator({})


@pytest.mark.unit
def test_component_ca_thermo():
    logging.getLogger("idaes.proteuslib.edb.data_model").setLevel(logging.DEBUG)
    comp = Component(testdata.Ca_thermo_data)
    generated_config = comp.idaes_config

    # for debugging
    print("Config generated from data:")
    pprint(generated_config)
    print("Expected idaes_config:")
    pprint(testdata.Ca_thermo_config)

    assert_configuration_equal(
        comp.idaes_config, testdata.Ca_thermo_config, "components"
    )
    logging.getLogger("idaes.proteuslib.edb.data_model").setLevel(logging.INFO)


@pytest.mark.unit
def test_reaction_bicarbonate():
    print(testdata.bicarbonate_reaction_data)
    react = Reaction(testdata.bicarbonate_reaction_data)
    generated_config = react.idaes_config

    # for debugging
    print("Config generated from data:")
    pprint(generated_config)
    print("Expected idaes_config:")
    pprint(testdata.bicarbonate_reaction_config)

    assert_configuration_equal(
        generated_config, testdata.bicarbonate_reaction_config, "equilibrium_reactions"
    )


@pytest.mark.unit
@pytest.mark.parametrize("data", testdata.reaction_data)
def test_reaction(data):
    reaction = Reaction(data)


@pytest.mark.unit
def test_result():
    Result(iterator=[], item_class=Component)
    Result(iterator=[], item_class=Reaction)


@pytest.mark.unit
@pytest.mark.parametrize(
    "data_wrapper_class,required",
    [
        (
            Component,
            {"name": "foo", "elements": ["H"], "type": "solvent", "parameter_data": {}},
        ),
        (Reaction, {"name": "foo", "type": "equilibrium", "parameter_data": {},
                    "components": ["H2O"]}),
    ],
)
def test_config_generator_required(data_wrapper_class, required):
    # make sure it fails with any missing required field(s)
    all_fields = list(required.keys())
    for i in range(len(all_fields)):
        except_i = all_fields[:i] + all_fields[i+1:]
        missing_a_field = {f: required[f] for f in except_i}
        print(f"Data missing field '{all_fields[i]}': {missing_a_field}")
        if data_wrapper_class is Component and all_fields[i] in ("elements", "type"):
            # 'elements' field will be filled by pre-processor for components
            _ = data_wrapper_class(missing_a_field).idaes_config
        else:
            with pytest.raises(ValidationError):
                _ = data_wrapper_class(missing_a_field).idaes_config
    # add all required fields and make sure it passes
    d = {r: required[r] for r in required}
    _ = data_wrapper_class(d).idaes_config


@pytest.mark.unit
@pytest.mark.parametrize("starting_value", [None, {}])
def test_base(starting_value):
    starting = dict.fromkeys(Component.merge_keys, starting_value)
    foo_value = 12.34
    mk0 = Component.merge_keys[0]
    starting[mk0] = {"foo": foo_value}
    b = Base(starting)
    assert b.idaes_config == starting
    # Add an empty component
    with pytest.raises(BadConfiguration):
        c = Component({})  # "name" is required
    # also required: elements, parameter_data
    c = Component({"name": "a", "elements": ["H +"], "parameter_data": {}})
    b.add(c)
    assert b.idaes_config[mk0]["foo"] == starting[mk0]["foo"]
    # Add a non-empty component
    name = "baz"
    component_data = {"name": name, "elements": [], "parameter_data": {"foo_coeff": [{"i": 0, "v": 1, "u": "g"}]}}
    c = Component(component_data)
    b.add(c)
    print(f"b.idaes_config={b.idaes_config} component_data={component_data}")
    config_foo_coeff_value = b.idaes_config[mk0][name]["parameter_data"]["foo_coeff"][0]
    data_foo_coeff_value = component_data["parameter_data"]["foo_coeff"][0]["v"]
    assert config_foo_coeff_value == data_foo_coeff_value


subst_foo, subst_bar, subst_y1, subst_y2 = "foo_obj", "bar_obj", 1, 2


class SubstituteTestGenerator(ConfigGenerator):
    pass


@pytest.mark.unit
def test_config_generator_substitute(debug_logging):
    SubstituteTestGenerator.substitute_values = {
        "a.b": {"foo": subst_foo},
        "x.*_bla": {"foo": subst_foo, "bar": subst_bar},
        "y": {"1": subst_y1, "2": subst_y2},
        "d.e.e.p": {"number_one": subst_y1},
        "d.e.e.p.e.*": {"number_two": subst_y2},
        "z.*": SubstituteTestGenerator.SUBST_UNITS,
    }
    data = {
        "a": {"b": "foo"},
        "x": {"one_bla": "bar", "two_bla": "hello", "three_bla": "foo", "ignore": "me"},
        "y": "1",
        "z": {"time": "s", "ignored_due_to_value": 0.123},
        "d": {"e": {"e": {"p": "number_one"}}},
    }
    print(f"before: {data}")
    SubstituteTestGenerator._substitute(data)
    print(f"after: {data}")
    assert data == {
        "a": {"b": "foo_obj"},
        "x": {
            "one_bla": "bar_obj",
            "two_bla": "hello",
            "three_bla": "foo_obj",
            "ignore": "me",
        },
        "y": 1,
        "z": {"time": pyunits.s, "ignored_due_to_value": 0.123},
        "d": {"e": {"e": {"p": 1}}},
    }


def test_substitute_phases(debug_logging):
    data = {"phases": {"Liq": {"type": "LiquidPhase"}}}
    SubstituteTestGenerator.substitute_values = {
        "phases.Liq.type": {"LiquidPhase": NIST}
    }
    print(f"before: {data}")
    SubstituteTestGenerator._substitute(data)
    print(f"after: {data}")
    assert data["phases"]["Liq"]["type"] == NIST


def test_component_from_idaes_config(debug_logging):
    H2O_thermo_config = {
        "components": {
            "H2O": {
                "type": IComponent,
                # Define the methods used to calculate the following properties
                "dens_mol_liq_comp": Perrys,
                "enth_mol_liq_comp": Perrys,
                "cp_mol_liq_comp": Perrys,
                "entr_mol_liq_comp": Perrys,
                "enth_mol_ig_comp": NIST,
                "pressure_sat_comp": NIST,
                "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
                # Parameter data is always associated with the methods defined above
                "parameter_data": {
                    "mw": (18.0153, pyunits.g / pyunits.mol),
                    "pressure_crit": (220.64e5, pyunits.Pa),
                    "temperature_crit": (647, pyunits.K),
                    # Comes from Perry's Handbook:  p. 2-98
                    "dens_mol_liq_comp_coeff": {
                        "1": (5.459, pyunits.kmol * pyunits.m ** -3),
                        "2": (0.30542, pyunits.dimensionless),
                        "3": (647.13, pyunits.K),
                        "4": (0.081, pyunits.dimensionless),
                    },
                    "enth_mol_form_liq_comp_ref": (-285.830, pyunits.kJ / pyunits.mol),
                    "enth_mol_form_vap_comp_ref": (0, pyunits.kJ / pyunits.mol),
                    # Comes from Perry's Handbook:  p. 2-174
                    "cp_mol_liq_comp_coeff": {
                        "1": (2.7637e5, pyunits.J / pyunits.kmol / pyunits.K),
                        "2": (-2.0901e3, pyunits.J / pyunits.kmol / pyunits.K ** 2),
                        "3": (8.125, pyunits.J / pyunits.kmol / pyunits.K ** 3),
                        "4": (-1.4116e-2, pyunits.J / pyunits.kmol / pyunits.K ** 4),
                        "5": (9.3701e-6, pyunits.J / pyunits.kmol / pyunits.K ** 5),
                    },
                    "cp_mol_ig_comp_coeff": {
                        "A": (30.09200, pyunits.J / pyunits.mol / pyunits.K),
                        "B": (
                            6.832514,
                            pyunits.J
                            * pyunits.mol ** -1
                            * pyunits.K ** -1
                            * pyunits.kiloK ** -1,
                        ),
                        "C": (
                            6.793435,
                            pyunits.J
                            * pyunits.mol ** -1
                            * pyunits.K ** -1
                            * pyunits.kiloK ** -2,
                        ),
                        "D": (
                            -2.534480,
                            pyunits.J
                            * pyunits.mol ** -1
                            * pyunits.K ** -1
                            * pyunits.kiloK ** -3,
                        ),
                        "E": (
                            0.082139,
                            pyunits.J
                            * pyunits.mol ** -1
                            * pyunits.K ** -1
                            * pyunits.kiloK ** 2,
                        ),
                        "F": (-250.8810, pyunits.kJ / pyunits.mol),
                        "G": (223.3967, pyunits.J / pyunits.mol / pyunits.K),
                        "H": (0, pyunits.kJ / pyunits.mol),
                    },
                    "entr_mol_form_liq_comp_ref": (
                        69.95,
                        pyunits.J / pyunits.K / pyunits.mol,
                    ),
                    "pressure_sat_comp_coeff": {
                        "A": (4.6543, None),  # [1], temperature range 255.9 K - 373 K
                        "B": (1435.264, pyunits.K),
                        "C": (-64.848, pyunits.K),
                    },
                },
            }
        }
    }
    logging.getLogger("proteuslib.edb.data_model").setLevel(logging.DEBUG)
    # create Component from the config
    components = Component.from_idaes_config(H2O_thermo_config)
    assert len(components) == 1
    component = components[0]
    # make sure it is a valid Component
    validate(component)
    # create a config from the Component, i.e. round-trip, and compare
    assert_configuration_equal(H2O_thermo_config, component.idaes_config, "components")
    logging.getLogger("proteuslib.edb.data_model").setLevel(logging.INFO)


def test_reaction_from_idaes_config(debug_logging):
    carbonation_reaction_config = {
        "base_units": {
            "time": pyunits.s,
            "length": pyunits.m,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "temperature": pyunits.K,
        },
        "equilibrium_reactions": {
            "CO2_to_H2CO3": {
                "stoichiometry": {
                    ("Liq", "H2O"): -1,
                    ("Liq", "CO2"): -1,
                    ("Liq", "H2CO3"): 1,
                },
                "heat_of_reaction": constant_dh_rxn,
                "equilibrium_constant": van_t_hoff,
                # "equilibrium_constant": gibbs_energy,
                "equilibrium_form": log_power_law_equil,
                "concentration_form": ConcentrationForm.molarity,
                "parameter_data": {
                    "dh_rxn_ref": (0, pyunits.kJ / pyunits.mol),
                    # NOTE: This 'ds_rxn_ref' was calculated to give 'keq' of 1.7*10**-3
                    # "ds_rxn_ref": (-53.0, pyunits.J/pyunits.mol/pyunits.K),
                    "k_eq_ref": (1.7 * 10 ** -3, None),
                    "T_eq_ref": (300, pyunits.K),
                    "reaction_order": {
                        ("Liq", "H2CO3"): 1,
                        ("Liq", "CO2"): -1,
                        ("Liq", "H2O"): 0,
                    },
                },
            }
        },
    }
    # create Component from the config
    result = Reaction.from_idaes_config(carbonation_reaction_config)
    assert len(result) == 1
    reaction = result[0]
    # make sure it is a valid Reaction
    validate(reaction)

    # in the generated config 'reaction_order' will be dropped since it is a
    # runtime addition, not something we store in the DB. So for comparison
    # purposes drop it from the original config as well
    del carbonation_reaction_config["equilibrium_reactions"]["CO2_to_H2CO3"][
        "parameter_data"
    ]["reaction_order"]

    # create a config from the Reaction, i.e. round-trip, and compare
    assert_configuration_equal(
        carbonation_reaction_config, reaction.idaes_config, "equilibrium_reactions"
    )


def test_thermoconfig_set_type():
    def _type(d):
        keys = list(d["components"].keys())
        return d["components"][keys[0]]["type"]

    # water -> solvent
    data = {
        "parameter_data": {
            "mw": [{"v": 18.0153, "u": "g/mol", "i": 0}],
        },
        "name": "H2O",
        "elements": ["O", "H"],
    }
    thermo = Component(data, validation=False)
    config = thermo.idaes_config
    assert _type(config) == Solvent

    # neutral non-water -> solute
    data = {
        "parameter_data": {
            "mw": [{"v": 96.994, "u": "g/mol", "i": 0}],
        },
        "name": "H2PO4",
        "elements": ["O", "H", "P"],
    }
    thermo = Component(data)
    config = thermo.idaes_config
    assert _type(config) == Solute

    # positive charge -> cation
    data = {
        "parameter_data": {
            "mw": [{"v": 57.002, "u": "g/mol", "i": 0}],
        },
        "name": "CaOH +",
        "elements": ["Ca", "O", "H"],
    }
    thermo = Component(data)
    config = thermo.idaes_config
    assert _type(config) == Cation

    # negative charge -> anion
    data = {
        "parameter_data": {
            "mw": [{"v": 17.008, "u": "g/mol", "i": 0}],
        },
        "name": "OH -",
        "elements": ["O", "H"],
    }
    thermo = Component(data)
    config = thermo.idaes_config
    assert _type(config) == Anion

    # already present, even if wrong, leave alone
    data = {
        "type": "solvent",
        "parameter_data": {
            "mw": [{"v": 17.008, "u": "g/mol", "i": 0}],
        },
        "name": "OH -",
        "elements": ["O", "H"],
    }
    thermo = Component(data)
    config = thermo.idaes_config
    assert _type(config) == Solvent
