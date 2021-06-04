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
                    "mw": {"$ref": "#/definitions/parameter"},
                    "pressure_crit": {"$ref": "#/definitions/parameter"},
                    "temperature_crit": {"$ref": "#/definitions/parameter"},
                },
                "patternProperties": {
                    "_coeff": {
                        "type": "object",
                        "additionalProperties": {"$ref": "#/definitions/parameter"},
                    },
                    "_ref": {"$ref": "#/definitions/parameter"},
                },
                "additionalProperties": False,
            },
        },
        "required": ["name", "type"],
        "patternProperties": {"_comp": {"type": "string"}},
        "additionalProperties": False,
        "definitions": {
            "parameter": {
                "type": "object",
                "description": "Value, units, etc. for a parameter",
                "properties": {
                    "v": {"description": "value", "type": "number"},
                    "u": {"description": "units", "type": "string"},
                    "i": {"description": "index", "type": "number"}
                },
                "required": ["v", "u"]
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
            "type": {"type": "string", "enum": ["equilibrium"], "description": "Type of reaction"},
            "name": {"type": "string", "description": "Name of reaction"},
            "stoichiometry": {
                "type": "object",
                "properties": {
                    "Liq": {"$ref": "#/definitions/stoichiometry"},
                    "Vap": {"$ref": "#/definitions/stoichiometry"}
                }
            },
            "heat_of_reaction": {"type": "string"},
            "equilibrium_constant": {"type": "string"},
            "equilibrium_form": {"type": "string"},
            "concentration_form": {"type": "string"},
            "parameter_data": {
                "type": "object",
                "patternProperties": {"_ref": {"$ref": "#/definitions/parameter"}},
                "additionalProperties": False,
            },
        },
        "required": ["name", "type"],
        "definitions": {
            "parameter": {
                "type": "array",
                "description": "List of parameter values",
                "items": {
                    "type": "object",
                    "description": "Value, units, etc. for a parameter",
                    "properties": {
                        "v": {"description": "value", "type": "number"},
                        "u": {"description": "units", "type": "string"},
                        "i": {"description": "index", "type": "number"}
                    },
                    "required": ["v", "u"]
                }
            },
            "stoichiometry": {
                "type": "object",
                "description": "Stoichiometry for a reaction",
                "patternProperties": {
                    "^[A-Z].*$": {
                        "type": "number",
                        "description": "Moles for the given species in the reaction. Negative for LHS, positive for RHS"
                    }
                }
            }
        }
    }
}

