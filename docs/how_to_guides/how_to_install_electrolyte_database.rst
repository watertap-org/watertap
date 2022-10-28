.. _install-edb:

How to install the Electrolyte Database (EDB)
=============================================

To install the Electrolyte Database (EDB), follow the steps discussed below.

**Install MongoDB**. The EDB uses `MongoDB <https://www.mongodb.com/>`_ as its storage engine. MongoDB is a third-party application that must be installed separately. The "community edition" of MongoDB is free, and has all the functionality needed for the EDB. To download and install it, go to the `MongoDB homepage <https://www.mongodb.com/>`_, or directly to the `MongoDB Community Server download page <https://www.mongodb.com/try/download/community>`_ (see :ref:`screenshot <screenshot-mongodb-download>`). On that page, select the correct operating system and follow the instructions to install the server.

**Load the data**. Some electrolyte data is distributed with WaterTAP to "bootstrap" the EDB. To load it, use the ``edb load`` command --- part of the :ref:`EDB command-line tools <edb-cli>` --- with the bootstrap option, from a shell or command window:

.. code-block::

   # Load the standard data into the default MongoDB database, running locally
   edb load -b

**Verify the installation**. If the above command works, the MongoDB server is running and the data should be loaded. You can verify this in a couple of ways:

* `Use the command-line program` to dump out the 'base' collection (which is small) to the console. In a shell environment where the Python package has been installed, run the following command:

.. code-block::

   edb dump -f '-' -t base

The result should be a bunch of text that resembles the following:

.. code-block::

   Wrote 2 record(s) from collection 'base' to file '<stdout>'
   [{"phases": {"Liq": {"type": "AqueousPhase", "equation_of_state": "Ideal"}}, "state_definition":
   "FTPx", "state_bounds": {"flow_mol": [0, 50, 100], "temperature": [273.15, 300, 650], "pressure":
   [50000, 100000, 1000000]}, "pressure_ref": 100000, "temperature_ref": 300, "base_units": {"time": "s",
   "length": "m", "mass": "kg", "amount": "mol", "temperature": "K"}, "name": "thermo"}, {"base_units":
   {"time": "s", "length": "m", "mass": "kg", "amount": "mol", "temperature": "K"}, "phases": {"Liq":
   {"type": "AqueousPhase", "equation_of_state": "Ideal"}}, "state_definition": "FTPx", "state_bounds":
   {"flow_mol": [0, 50, 100], "temperature": [273.15, 300, 650], "pressure": [50000.0, 100000.0, 1000000.0]},
   "pressure_ref": 100000.0, "temperature_ref": 300, "name": "water_reaction"}]

* `Use MongoDB's graphical user interface`, "MongoDB Compass", to browse the data. To do this, find and start the application called "MongoDB Compass", which should have been installed when you installed the rest of the MongoDB application. Run it, and choose to connect to the server at URL ``mongodb://localhost:27017`` (this should be the default). You will get a screen like :ref:`this one <screenshot-mongodb-compass-initial>` (with the database you are going to click on next circled). Then, select the "electrolytedb" database. The result should show three collections with some records loaded in each, as in :ref:`this screen <screenshot-mongodb-compass-edb>`.

.. rubric:: Screenshots

.. _screenshot-mongodb-download:

.. figure:: ../_static/mongodb-download-page.*

   Download page for MongoDB community server (9/2021)

.. _screenshot-mongodb-compass-initial:

.. figure:: ../_static/mongodb-compass-initial.*

   MongoDB Compass Initial Screen (9/2021)

.. _screenshot-mongodb-compass-edb:

.. figure:: ../_static/mongodb-compass-electrolytedb.*

   MongoDB Compass electrolytedb Collections (9/2021)
