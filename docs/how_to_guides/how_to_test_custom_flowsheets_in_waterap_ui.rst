.. _how_to_test_custom_flowsheets_in_watertap_ui:

How to test custom flowsheets in the WaterTAP User Interface
============================================================

This guide describes how to install and run the IDAES Flowsheet Processor User Interface (UI)
and test flowsheets from a local WaterTAP pull request (PR). The installation
and UI launch steps follow the `IDAES Flowsheet Processor UI instructions
<https://github.com/prommis/idaes-flowsheet-processor-ui>`_.

Installation
------------

Before starting, clone the IDAES Flowsheet Processor User Interface repository and navigate to its root directory. Replace ``<idaes-flowsheet-processor-ui-path>`` below with the actual path to the UI repository:

.. code-block:: python

	git clone https://github.com/prommis/idaes-flowsheet-processor-ui.git
	cd <idaes-flowsheet-processor-ui-path>

Create the Conda environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

From the UI repository root (the directory containing ``environment.yml``), create and activate the
``flowsheet-processor-env`` environment:

.. code-block:: python

	conda env create --file environment.yml && conda activate flowsheet-processor-env

Install the JavaScript dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Change to the ``frontend`` directory within the UI repository path and install the UI dependencies. This
creates ``frontend/node_modules``:

.. code-block:: python

	cd frontend
	npm ci

Install the IDAES solver dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With ``flowsheet-processor-env`` environment active, install the IDAES solver dependencies:

.. code-block:: python

	idaes get-extensions --verbose

Install WaterTAP with default flowsheets (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To use the default WaterTAP flowsheets, install the published package:

.. code-block:: python

	pip install watertap

To test a custom WaterTAP flowsheet in your local WaterTAP instead, use the local installation described in
:ref:`how_to_gui_local_watertap`.

Running the UI
--------------

Ensure that all dependencies are installed and activate the Conda environment:

.. code-block:: python

	conda activate flowsheet-processor-env

From the UI repository root change to ``frontend`` and run ``app-start`` to start the UI. Replace ``<idaes-flowsheet-processor-ui-path>`` below with the actual path to the UI repository.

.. code-block:: python

	cd <idaes-flowsheet-processor-ui-path>/frontend
	npm run app-start

.. _how_to_gui_local_watertap:

Running custom WaterTAP flowsheets
----------------------------------

Follow the detailed steps below to test custom WaterTAP flowsheets in the GUI. 
Before starting, ensure that you have cloned the IDAES Flowsheet Processor User Interface repository 
and set up the necessary environment as described in the previous sections.	

1. Go to you local WaterTAP repository. Replace ``<local-watertap-repo-path>`` below with the actual path to your local WaterTAP repository.
2. Activate ``flowsheet-processor-env``
3. Check out the branch in your local WaterTAP repository containing the flowsheet you want to test in the UI. Replace ``<branch-containing-flowsheet>`` below with the actual branch name containing the flowsheet.
4. Install WaterTAP from the local WaterTAP repository root. This is the directory containing ``pyproject.toml``.

	.. code-block:: python

		conda activate flowsheet-processor-env
		cd <local-watertap-repo-path>
		git checkout <branch-containing-flowsheet>
		pip install .

5. Return to the UI repository's ``frontend`` directory and run the UI. Replace ``<idaes-flowsheet-processor-ui-path>`` below with the actual path to the UI repository.

	.. code-block:: python

		cd <idaes-flowsheet-processor-ui-path>/frontend
		npm run app-start

	This launches the WaterTAP GUI in your web browser.

Stop any running UI before installing the local WaterTAP package, then restart it to load the updated flowsheets. Repeat the local installation and restart
after making further changes to WaterTAP.

Custom flowsheets must expose the UI interface described in :ref:`how_to_ui-api` to be available in the GUI.

.. note::
   For detailed information on developer testing of custom flowsheets, refer to the README instructions on the `IDAES Flowsheet Processor UI repository <https://github.com/prommis/idaes-flowsheet-processor-ui>`_.
