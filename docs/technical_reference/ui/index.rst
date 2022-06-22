.. _ref_ui-api:

WaterTAP User Interface API
===========================
This is the Python API that will be called by the UI backend to process
changes from the UI front-end.

.. include:: <isonum.txt>

.. py:currentmodule:: watertap.ui.api

.. contents:: Contents
    :depth: 1
    :local:

Introduction
------------

This page describes an application programming interface (API) that is designed to help communicate information about a flowsheet and its variables to a user interface layer (or possibly several different kinds of user interface layers at the same time).
The API also provides the ability, intended for the UI developer, to update the variable values and run "actions" such as building and solving the flowsheet.

For how-to guides on specific tasks, see the :ref:`howto_ui-api`.

.. .. image:: /_static/search-icon.png
.. x :height: 65px
.. x  :align: left

Class and function reference
----------------------------

* :mod:`watertap.ui.api`

High-level API
^^^^^^^^^^^^^^

* :func:`export_variables`

* :func:`find_flowsheet_interfaces`

* :class:`FlowsheetInterface`

* :class:`WorkflowActions`

Lower-level API
^^^^^^^^^^^^^^^

* :func:`get_block_interface`

* :func:`set_block_interface`

* :class:`BlockInterface`

Flowsheet information format
----------------------------
The information about the flowsheet and all its subblocks and variables is encoded in JSON when it is transferred between the backend and the UI, or saved to a file.
The form of this information is shown in the JSON pseudo-schema in the header of the :mod:`watertap.ui.api` module.

Below is a more formal (JSON-Schema) schema for the information:

.. code-block:: json

    {
      "$ref": "#/definitions/Block",
      "definitions": {
        "IndexedValue": {
          "title": "IndexedValue",
          "type": "object",
          "properties": {
            "index": {
              "title": "Index",
              "type": "array",
              "items": {
                "type": "array",
                "items": {
                  "anyOf": [
                    {
                      "type": "number"
                    },
                    {
                      "type": "string"
                    }
                  ]
                }
              }
            },
            "value": {
              "title": "Value",
              "type": "array",
              "items": {
                "anyOf": [
                  {
                    "type": "number"
                  },
                  {
                    "type": "string"
                  }
                ]
              }
            }
          },
          "required": [
            "index",
            "value"
          ]
        },
        "ScalarValue": {
          "title": "ScalarValue",
          "type": "object",
          "properties": {
            "value": {
              "title": "Value",
              "anyOf": [
                {
                  "type": "number"
                },
                {
                  "type": "string"
                }
              ]
            }
          },
          "required": [
            "value"
          ]
        },
        "Variable": {
          "title": "Variable",
          "type": "object",
          "properties": {
            "value": {
              "title": "Value",
              "anyOf": [
                {
                  "$ref": "#/definitions/IndexedValue"
                },
                {
                  "$ref": "#/definitions/ScalarValue"
                }
              ]
            },
            "display_name": {
              "title": "Display Name",
              "default": "",
              "type": "string"
            },
            "description": {
              "title": "Description",
              "default": "",
              "type": "string"
            },
            "units": {
              "title": "Units",
              "default": "",
              "type": "string"
            },
            "readonly": {
              "title": "Readonly",
              "default": false,
              "type": "boolean"
            }
          }
        },
        "BlockMeta": {
          "title": "BlockMeta",
          "type": "object",
          "properties": {
            "parameters": {
              "title": "Parameters",
              "default": {},
              "type": "object"
            }
          }
        },
        "Block": {
          "title": "Block",
          "type": "object",
          "properties": {
            "variables": {
              "title": "Variables",
              "default": {},
              "type": "object",
              "additionalProperties": {
                "$ref": "#/definitions/Variable"
              }
            },
            "blocks": {
              "title": "Blocks",
              "default": {},
              "type": "object",
              "additionalProperties": {
                "$ref": "#/definitions/Block"
              }
            },
            "meta": {
              "title": "Meta",
              "default": {
                "parameters": {}
              },
              "allOf": [
                {
                  "$ref": "#/definitions/BlockMeta"
                }
              ]
            },
            "display_name": {
              "title": "Display Name",
              "default": "",
              "type": "string"
            },
            "description": {
              "title": "Description",
              "default": "",
              "type": "string"
            },
            "category": {
              "title": "Category",
              "default": "default",
              "type": "string"
            }
          }
        }
      }
    }
