.. _notebooks:

How to run WaterTAP with Jupyter notebooks
==========================================

This guide describes the supported methods to work with WaterTAP using Jupyter notebooks.

Local (developer) installation
------------------------------

.. important:: These instructions assume that a WaterTAP developer environment has already been configured, as described in e.g. :ref:`install-dev`.

1. Activate your WaterTAP environment. If your environment has a name different from ``watertap-dev``, replace ``watertap-dev`` with the actual name wherever applicable:

   .. code-block:: shell
    conda activate watertap-dev

1. Navigate to the directory where your local clone of the WaterTAP repository is located

1. Run the following command to register the currently active WaterTAP environment as a Jupyter kernel. This will create a dedicated WaterTAP kernel that's selectable in the Jupyter web interface, ensuring that the notebook code is run in the correct Python environment where WaterTAP is available.

   .. code-block:: shell

   python -m ipykernel --user --name "watertap-dev"

1. Then, start the Jupyter server from the current directory, navigate to the desired notebook using the file browser, and launch it
1. If prompted to select a kernel, select "watertap-dev" from the menu. Otherwise, ensure "watertap-dev" appears in the kernel box in the top-right corner of the notebook interface

Online using Binder.org
-----------------------

(Coming soon)
