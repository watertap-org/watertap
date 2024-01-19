.. _notebooks:

How to run WaterTAP with Jupyter notebooks
==========================================

This guide describes the supported methods to work with WaterTAP using Jupyter notebooks.

Local (developer) installation
------------------------------

.. important:: These instructions assume that a WaterTAP developer environment has already been configured, as described in e.g. :ref:`install-dev`.

#. Activate your WaterTAP environment. If your environment has a name different from ``watertap-dev``, replace ``watertap-dev`` with the actual name wherever applicable:

   .. code-block:: shell

      conda activate watertap-dev

#. Navigate to the directory where your local clone of the WaterTAP repository is located

#. Run the following command to register the currently active WaterTAP environment as a Jupyter kernel. This will create a dedicated WaterTAP kernel that's selectable in the Jupyter web interface, ensuring that the notebook code is run in the correct Python environment where WaterTAP is available.

   .. code-block:: shell

      python -m ipykernel install --user --name "watertap-dev"

#. Then, start the Jupyter server from the current directory, navigate to the desired notebook using the file browser, and launch it

#. If prompted to select a kernel, select ``watertap-dev`` from the menu. Otherwise, ensure ``watertap-dev`` appears in the kernel box in the top-right corner of the notebook interface

.. _notebooks-online:

Online using Binder.org
-----------------------

`Binder <https://mybinder.org>`_ is an online service providing a **short-lived temporary sandbox environment on public cloud resources** where Jupyter notebooks can be run, free of charge, without having to install any software locally.

Click on this button to launch an environment pointing to the current ``main`` branch of the WaterTAP repository:

.. image:: https://mybinder.org/badge_logo.svg
 :target: https://mybinder.org/v2/gh/watertap-org/watertap/main?labpath=tutorials%2Fintroduction.ipynb

.. important::

   A Binder environment is automatically **destroyed after a few minutes of inactivity**, which means that **any unsaved progress will be lost**. To avoid this, users should download a copy of a notebook file from the Binder environment to their local machine through the browser (in the Jupyter Lab file browser menu, right click on the notebook file, then ``Download``). For more information, see https://mybinder.readthedocs.io/en/latest/about/about.html#how-long-will-my-binder-session-last.

Customization
^^^^^^^^^^^^^

Binder uses a Git repository hosted on, e.g., GitHub to fetch the notebooks and create the runtime environment.
Users can specify a specific fork, branch, and Git ref (e.g., a particular commit hash), as well as the path to a particular directory or notebook file within the repository, that will be used
when first starting the environment.

These options can be specified interactively on the `Binder homepage <https://mybinder.org/>`_, which will create a URL that can then be shared with others to generate a (separate) instance of the environment with the same repository settings.

Alternatively, the URL can be generated manually according to the following schema::

   https://mybinder.org/v2/gh/<github-org-or-user-name>/<github-repo-name>/<git-ref>?labpath=<path-to-notebook>

Example, for branch ``mybranch`` of ``myuser``'s fork of this repository, pointing to the ``tutorials/introduction.ipynb`` notebook file::

   https://mybinder.org/v2/gh/myuser/watertap/mybranch?labpath=tutorials/introduction.ipynb
