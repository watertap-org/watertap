EDB Installation
================

`In the following instructions the electrolyte database will be abbreviated "EDB".`

1. **Install MongoDB**. The EDB uses `MongoDB <https://www.mongodb.com/>`_ as its storage engine, so first you must install MongoDB.
   Then make sure MongoDB is running locally on your system. MongoDB comes with a decent desktop UI called Compass, but
   you may also choose to try a (free) third-party UI like `Robo3T <https://robomongo.org/>`_.
2. **Install ProteusLib**. The EDB is distributed as part of `ProteusLib <https://github.com/nawi-hub/proteuslib>`_.
   Follow the ProteusLib installation instructions.
3. **Load data**. Some electrolyte data is distributed with ProteusLib to bootstrap using the EDB.
   To load it, use the :ref:`edb command-line program <edb-cli>`::

    # Load the standard data into a database called, unimaginatively, 'edb'
    edb load --bootstrap -d edb