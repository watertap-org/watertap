Electrolyte Database (EDB)
==========================

.. contents::
    :local:
    :depth: 1

.. _edb-api:

Overview
--------
The Electrolyte Database (EDB) stores metadata and data about chemical species, called here
`components`, and `reactions`. It is accessed through a Python API to return well-defined Python objects.

The data are stored in `MongoDB <https://mongodb.org>`_, so they can be queried in a number of ways, and the
system is extensible to new use-cases. The native storage format for MongoDB is a `JSON <https://json.org>`_ document,
and the expected structure and fields of the *component* and *reaction* data is defined by a
`JSON Schema <https://json-schema.org>`_. Validation using those schemas is built into the API (though it can be disabled).

To interface with the `IDAES Core Modeling Framework <https://idaes-pse.readthedocs.io/en/stable/user_guide/concepts.html>`_
(IDAES-CMF, which underlies WaterTAP), add components and reactions to a "base" object and fetch the result as a Python
`dict`. This result can be used to configure and build IDAES objects (`ParameterBlocks`, `ReactionBlocks`, etc.).
The API also has methods to construct component and reaction objects from IDAES configurations.

Examples for the usage of the EDB for creating these IDAES configurations can be found in the
:ref:`EDB How-to Documentation<how_to_use_edb>`

Workflows
---------
The EDB is intended to support some known workflows out of the box, with lower-level functions available when these
are not sufficient. Example workflows can be seen in the how to guides from the link above.


Python API
----------

Database API
^^^^^^^^^^^^
Connect to the database and create, read, update and delete its contents.

.. automodule:: watertap.edb.db_api
    :members: ElectrolyteDB
    :noindex:

Data object API
^^^^^^^^^^^^^^^
Data models for components and reactions, including conversion to IDAES config objects.

.. automodule:: watertap.edb.data_model
    :members: Base, Component, Reaction, Result
    :noindex:

.. _edb-cli:

.. program:: edb

edb command-line
----------------
The ``edb`` command-line program lets you interact with the database and the data schemas from a terminal.


edb base
^^^^^^^^

The work of the program is all done by subcommands.

edb base options
++++++++++++++++

.. option:: --help

    Show options and subcommands

.. option::  -v, --verbose

    Increase verbosity

.. option::  -q, --quiet

    Increase quietness

.. program:: edb-load

.. ###########################################################

edb load
^^^^^^^^

Load JSON records into the EDB.

edb load options
++++++++++++++++

.. option::  -f, --file FILENAME

    File to load  [required]

.. option::  -t, --type [component|reaction|base]

    Type of records  [required]

.. option::  -u, --url TEXT

    Database connection URL

.. option::  -d, --database TEXT

    Database name

.. option::  --validate / -n, --no-validate

    Turn on or off validation of input

.. option:: -b, --bootstrap

    Bootstrap a new database by loading in the standard base data.

.. ###########################################################

edb dump
^^^^^^^^

Dump JSON records from the EDB to a file.

edb dump options
++++++++++++++++

.. option::  -f, --file FILENAME

     File to create (will overwrite existing files!)  [required]

.. option::  -t, --type [component|reaction|base]

    Type of records (MongoDB collection name)

.. option::  -u, --url TEXT

    Database connection URL

.. option::  -d, --database TEXT

    Database name

.. ###########################################################

edb schema
^^^^^^^^^^

Show JSON schemas, in raw or readable forms, for the different record types.

edb schema options
++++++++++++++++++

.. option:: -f, --file FILENAME

    Write output to this file instead of printing to the screen

.. option::  -o, --format [json|markdown|html|html-js]

     Output format

.. option::  -t, --type [component|reaction]

    Type of records  [required]

.. option::  -u, --url TEXT

    Database connection URL

.. option::  -d, --database TEXT

    Database name

.. ###########################################################

.. edb-schemata

EDB schemas
-----------
The EDB data is encoded in `JSON <https://json.org>`_.
Naturally, the expected form of the records is specified as JSON Schema.
There are schemas for the `component` and `reaction` records. Currently, there is no schema for `base` data (this
will change soon, though).

Component schema
^^^^^^^^^^^^^^^^

.. only:: html

    .. raw:: html
        :file: schemas/component.html

.. only:: latex or epub or text

    .. include:: schemas/component.json
        :literal:

Reaction schema
^^^^^^^^^^^^^^^

.. only:: html

    .. raw:: html
        :file: schemas/reaction.html

.. only:: latex or epub or text

    .. include:: schemas/reaction.json
        :literal:


EDB base config options
-----------------------
The EDB data for the most common base configs are made available in the
standard bootstrap bundled with WaterTAP. Those bases are listed below:

+----------------------+-------------------------------------------------------------------------------------------+
|     Base             |  Description                                                                              |
+======================+===========================================================================================+
| default_thermo       | Default ThermoConfig: contains only AqueousPhase and uses FTPx state vars                 |
+----------------------+-------------------------------------------------------------------------------------------+
| thermo_Liq_FpcTP     | ThermoConfig: contains only AqueousPhase and uses FpcTP state vars                        |
+----------------------+-------------------------------------------------------------------------------------------+
| thermo_Liq_Sol_FpcTP | ThermoConfig: contains AqueousPhase + SolidPhase and uses FpcTP state vars                |
+----------------------+-------------------------------------------------------------------------------------------+
| thermo_Liq_Vap_FpcTP | ThermoConfig: contains AqueousPhase + VaporPhase and uses FpcTP state vars                |
+----------------------+-------------------------------------------------------------------------------------------+
| reaction             | ReactionConfig: Blank template for reaction configs                                       |
+----------------------+-------------------------------------------------------------------------------------------+

The naming convention for the **thermo** bases is as follows: (i) each section of the name is broken up
by the underscore character (_), (ii) the first word in the name is always **thermo** to denote the
type of base, (iii) the last word in the name always denotes the **state_vars**, and (iv) each word in
between denotes the set of **phases** for the configuration file. 
