Electrolyte Database (EDB)
==========================

.. contents::
    :local:
    :depth: 1

.. _edb-api:

Overview
--------
The Electrolyte Database (EDB) stores data about chemical species.
It is accessed through a Python API to return well-defined Python objects.

The data are stored in a database called `MongoDB <https://mongodb.org>`_.
MongoDB is open-source and free software, so you can install and use the database locally.
You may also use our open cloud-hosted version of the database.
See the :ref:`MongoDB installation instructions <install-mongodb>` for more details.

The native storage format for MongoDB is a `JSON <https://json.org>`_ object, which MongoDB calls a "document".
The expected structure of the EDB data is defined by a `JSON Schema <https://json-schema.org>`_.
Most users will not need to deal with the MongoDB documents, as there is a :ref:`Python data API <edb-data-api>` for searching the database and building IDAES "config blocks" from the *components* and *reactions*.

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

.. _edb-data-api:

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
