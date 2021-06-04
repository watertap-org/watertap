EDB Installation
================

`In the following instructions the electrolyte database will be abbreviated "EDB".`

1. **Install MongoDB**. The EDB uses MongoDB as its storage engine, so first you must install MongoDB. Then make
sure MongoDB is running locally on your system.

2. **Install ProteusLib**. The EDB is distributed as part of ProteusLib. Follow the ProteusLib installation instructions.

3. **Load data**. Some electrolyte data is distributed with ProteusLib to bootstrap using the EDB. To load it, use
the :ref:`edb command-line program <edb-cli>`::

    # note: --bootstrap is not yet implemented
    # Here 'edb' is being used as the database name. You can choose another name if you like.
    # This assumes that the MongoDB server is running on your system, on the standard port.
    edb load --bootstrap -d edb