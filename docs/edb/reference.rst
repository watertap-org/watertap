EDB Reference
=============

.. _edb-cli:

.. program:: edb

edb
---
The ``edb`` command-line program lets you interact with the database and the data schemas from a terminal.
The work of the program is all done by subcommands.

edb options
^^^^^^^^^^^

.. option:: --help

    Show options and subcommands

.. option::  -v, --verbose

    Increase verbosity

.. option::  -q, --quiet

    Increase quietness

.. program:: edb-load

.. ###########################################################

edb load
--------

Load JSON records into the EDB.

edb load options
^^^^^^^^^^^^^^^^

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

.. ###########################################################

edb dump
--------

Dump JSON records from the EDB to a file.

edb dump options
^^^^^^^^^^^^^^^^

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
----------

Show JSON schemas, in raw or readable forms, for the different record types.

edb schema options
^^^^^^^^^^^^^^^^^^

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