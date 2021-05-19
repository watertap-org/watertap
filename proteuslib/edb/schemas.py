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
# WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
#
# This module is a work in progress. Do not use it for real work right now.
#
# WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING

"""
JSON schema embedded as variables for:
  - component
  - reaction
"""

schemas = {
    "component": {
        "$schema": "http://json-schema.org/draft-07/schema#",
        "$id": "https://nawi-hub.github.com/electrolytedb/component",
        "title": "Electrolyte database component",
        "type": "object",
        "description": "Electrolyte component schema",
        "properties": {
            "type": {"type": "string", "const": "component"},
            "name": {
                "description": "The chemical name of the component",
                "examples": ["HPO4 2-", "H3PO4", "NH4 +"],
                "type": "string",
            },
            "valid_phase_types": {"type": "string"},
            "phase_equilibrium_form": {
                "type": "object",
                "properties": {"vap": {"type": "string"}, "liq": {"type": "string"}},
            },
            "parameter_data": {
                "type": "object",
                "properties": {
                    "mw": {"$ref": "#/definitions/value_units"},
                    "pressure_crit": {"$ref": "#/definitions/value_units"},
                    "temperature_crit": {"$ref": "#/definitions/value_units"},
                },
                "patternProperties": {
                    "_coeff": {
                        "type": "object",
                        "additionalProperties": {"$ref": "#/definitions/value_units"},
                    },
                    "_ref": {"$ref": "#/definitions/value_units"},
                },
                "additionalProperties": False,
            },
        },
        "required": ["name", "type"],
        "patternProperties": {"_comp": {"type": "string"}},
        "additionalProperties": False,
        "definitions": {
            "value_units": {
                "type": "array",
                "description": "Value and units as a pair in an array",
                "items": [
                    {"description": "value", "type": "number"},
                    {"description": "units", "type": "string"},
                ],
                "additionalItems": False,
            }
        },
    },
    "reaction": {
        "$schema": "http://json-schema.org/draft-07/schema#",
        "$id": "https://nawi-hub.github.com/electrolytedb/reaction",
        "title": "Electrolyte database reaction",
        "description": "Electrolyte reaction schema",
        "type": "object",
        "properties": {
            "type": {"type": "string", "const": "reaction"},
            "reaction_type": {"type": "string", "enum": ["equilibrium"]},
            "name": {"type": "string"},
            "stoichiometry": {"type": "object"},
            "heat_of_reaction": {"type": "string"},
            "equilibrium_constant": {"type": "string"},
            "equilibrium_form": {"type": "string"},
            "concentration_form": {"type": "string"},
            "parameter_data": {
                "type": "object",
                "properties": {"reaction_order": {"type": "object"}},
                "patternProperties": {"_ref": {"$ref": "#/definitions/value_units"}},
                "additionalProperties": False,
            },
        },
        "required": ["name", "type"],
        "definitions": {
            "value_units": {
                "type": "array",
                "description": "Value and units as a pair in an array",
                "items": [
                    {"description": "value", "type": "number"},
                    {"description": "units", "type": "string"},
                ],
                "additionalItems": False,
            }
        }
    }
}

