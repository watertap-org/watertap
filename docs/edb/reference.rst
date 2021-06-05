EDB Reference
=============

.. contents::
    :local:
    :depth: 2

.. _edb-api:

Python API
----------

Database API
^^^^^^^^^^^^
Connect to the database and create, read, update and delete its contents.

.. automodule:: proteuslib.edb.db_api
    :members: ElectrolyteDB
    :noindex:

Data object API
^^^^^^^^^^^^^^^
Data models for components and reactions, including conversion to IDAES config objects.

.. automodule:: proteuslib.edb.data_model
    :members: Base, Component, Reaction, Result
    :noindex:

.. _edb-cli:

.. program:: edb

edb
---
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
Naturally, the expected form of the records is specified as `JSON Schema <https://json-schema.org/>`_.
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
