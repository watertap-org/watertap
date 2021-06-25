EDB Examples
============

edb command-line program
------------------------
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