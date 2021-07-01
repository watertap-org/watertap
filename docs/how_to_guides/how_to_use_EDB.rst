How to use the Electrolyte Database (EDB)
=========================================

Installation
------------

1. **Install MongoDB**. The EDB uses `MongoDB <https://www.mongodb.com/>`_ as its storage engine, so first you must install MongoDB.
   Then make sure MongoDB is running locally on your system. MongoDB comes with a decent desktop UI called Compass, but
   you may also choose to try a (free) third-party UI like `Robo3T <https://robomongo.org/>`_.
2. **Install ProteusLib**. The EDB is distributed as part of `ProteusLib <https://github.com/nawi-hub/proteuslib>`_.
   Follow the ProteusLib installation instructions.
3. **Load data**. Some electrolyte data is distributed with ProteusLib to bootstrap using the EDB.
   To load it, use the :ref:`edb command-line program <edb-cli>`::

    # Load the standard data into a database called, unimaginatively, 'edb'
    edb load --bootstrap -d edb


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