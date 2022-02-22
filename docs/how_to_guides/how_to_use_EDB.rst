.. _how_to_use_edb:

How to use the Electrolyte Database (EDB)
=========================================

.. toctree::
    :maxdepth: 1

    edb_cli
    edb_api
    ../examples/edb/edb_example
    ../examples/edb/simple_acid_example
    ../examples/edb/multiple_phases_example

.. _install-mongodb:

Installation
------------
Before running any of the EDB example Jupyter Notebooks, you **must** have MongoDB :ref:`installed locally <install-mongodb-local>` and you need to load the 'bootstrap' data.
See the :ref:`EDB CLI How-to <how_to_use_edb_cli>` for information on loading the default (i.e., 'bootstrap' data).
For your own work, you can also use the public cloud database, which has the data pre-loaded, as described below.

.. _install-mongodb-local:

Local installation
^^^^^^^^^^^^^^^^^^
For local installation, download and install the MongoDB "community server", which is free.
At the time of this writing, the web page for this is `here <https://www.mongodb.com/try/download/community>`_.
Any recent version should work.
The default settings should work out of the box with the EDB client tools and API.
MongoDB now ships with a useful graphical user interface called Compass (see `Compass documentation here <https://docs.mongodb.com/compass/master/>`_), that you can use to view and even modify the database.

For more details on setting up the EDB, see :ref:`install-edb`.

.. _use-cloud-edb:

Using the Cloud-hosted Database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The EDB requires a running database (MongoDB) server that is loaded with the correct data.
You can install and run this server locally, but you may also use the cloud-hosted public database.
To use this database, provide the following URL for Python and command-line functions::

    mongodb+srv://<username>:<password>@nawi-edb.utpac.mongodb.net

For the public read-only user, the value for <username> is ``edbnawi`` and <password> is ``edb-user``.
You need to provide this URL to any command-line or API functions that connect to the database.

As with any distributed system, it is possible that the system is down or unresponsive.
To quickly test whether this URL is operating correctly you can dump the (small) 'base' collection from the command-line with::

    edb dump -u  '<INSERT-URL-HERE>' -f '-' -t base

This should produce output that starts with the following line, followed by some JSON output::

    Wrote 2 record(s) from collection 'base' to file '<stdout>'



