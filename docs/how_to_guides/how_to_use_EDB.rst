How to use the Electrolyte Database (EDB)
=========================================

For information on how to install and set up the EDB, see :ref:`install-edb`.

.. _use-cloud-edb:

Using the public cloud EDB
--------------------------
The EDB requires a running database (MongoDB) server that is loaded with the correct data.
You can install and run this server locally, but you may also use the cloud-hosted public database.
To use this database, provide the following URL for Python and command-line functions:

``mongodb+srv://<username>:<password>@nawi-edb.utpac.mongodb.net``

Where the value for <username> is ``edbnawi`` and <password> is ``edb-user``.

As with any distributed system, it is possible that the system is down or unresponsive.
To quickly test whether this URL is operating correctly you can dump the (small) 'base' collection from the command-line with::

    edb dump -u  '<INSERT-URL-HERE>' -f '-' -t base

This should produce output that starts with the following line, followed by some JSON output::

    Wrote 2 record(s) from collection 'base' to file '<stdout>'

EDB examples
------------

edb command-line program
^^^^^^^^^^^^^^^^^^^^^^^^
Below are some examples of using the ``edb`` command-line program.

|arrw| Load records from ``reaction.json`` into the `reaction` collection in the `test` database::

    edb load -f reaction.json -t reaction -d test


|arrw| Dump the contents of the `base` collection from the `edb` database into the file ``data/base.json``::

    edb dump -t base -d edb -f data/base.json

|arrw| Print the raw JSON-Schema for the `component` collection records to the screen::

    edb schema -t component -o json


|arrw| Document the schema of the `component` collection records with an interactive web page
stored at ``./docs/schemas/edb_component_schema.html``. Supporting files will also be added in the same directory::

    edb schema -t component -f ./docs/schemas/edb_component_schema.html -o html-js


.. |arrw| unicode:: U+27A2 .. nice looking arrow glyph