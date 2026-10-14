.. _how_to_test_custom_flowsheets_in_watertap_ui:

How to test custom flowsheets in the WaterTAP UI
=================================================

This guide describes how to install and run the IDAES Flowsheet Processor UI
and test flowsheets from a local WaterTAP pull request (PR). The installation
and UI launch steps follow the `IDAES Flowsheet Processor UI instructions
<https://github.com/prommis/idaes-flowsheet-processor-ui>`_.

Installation
------------

Before starting, install Conda, Node.js, npm, and Git, and clone the UI
repository:

.. code-block:: console

	git clone https://github.com/prommis/idaes-flowsheet-processor-ui.git
	cd idaes-flowsheet-processor-ui

Create the Conda environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

From the UI repository root, create and activate the
``flowsheet-processor-env`` environment:

.. code-block:: console

	conda env create --file environment.yml && conda activate flowsheet-processor-env

Install the JavaScript dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Change to the ``frontend`` directory and install the UI dependencies. This
creates ``frontend/node_modules``:

.. code-block:: console

	cd frontend
	npm ci

Install the IDAES solver dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With ``flowsheet-processor-env`` active, install the IDAES solver dependencies:

.. code-block:: console

	idaes get-extensions --verbose

Install WaterTAP with default flowsheets (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To use the default WaterTAP flowsheets, install the published package:

.. code-block:: console

	pip install watertap

To test a local WaterTAP PR instead, use the local installation described in
:ref:`howto_gui_local_watertap`.

Running the UI
--------------

Ensure that all dependencies are installed and activate the Conda environment:

.. code-block:: console

	conda activate flowsheet-processor-env

From the UI repository root, change to ``frontend`` and start the UI in your
browser. If you are already in ``frontend``, skip the ``cd`` command:

.. code-block:: console

	cd frontend
	npm run app-start

.. _howto_gui_local_watertap:

Running WaterTAP flowsheets from a local PR
-------------------------------------------

1. In your local WaterTAP repository, check out the branch containing the PR
	you want to test.
2. Activate ``flowsheet-processor-env`` and install WaterTAP from the local
	repository root (the directory containing ``pyproject.toml``):

	.. code-block:: console

		conda activate flowsheet-processor-env
		pip install .

3. Return to the UI repository's ``frontend`` directory and run:

	.. code-block:: console

		npm run app-start

Stop any running UI before installing the local WaterTAP package, then restart
it to load the updated flowsheets. Repeat the local installation and restart
after making further changes to WaterTAP.

Custom flowsheets must expose the UI interface described in
:ref:`howto_ui-api` to be available in the GUI.
