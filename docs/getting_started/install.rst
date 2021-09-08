Installing ProteusLib
=====================

Introduction
------------

.. _about-conda:

Using Conda environments
^^^^^^^^^^^^^^^^^^^^^^^^

Conda environments are a way to create and manage multiple sets of packages and/or Python versions on the same system without incurring conflicts.
Each Conda environment is a dedicated directory, separate from other Conda environments and the operating system's own directories, containing its own collection of packages, executables, and Python installation, including the Python interpreter.
Once a Conda environment is *activated*, using the ``conda activate`` command in a terminal/console, the environment's own version of Python will be used to run commands or interactive sessions for the remainder of the session.

For these reasons, Conda environments are especially useful to install and manage multiple projects (and/or multiple *versions* of the same project) on the same computer with minimal effort,
as they provide a way to seamlessly switch between different projects without conflicts.

Using Conda environments is not mandatory to be able to install and use ProteusLib; however, it is strongly recommended.

To use Conda environments, the ``conda`` package manager is required. Refer to the `Conda installation guide <https://docs.conda.io/projects/conda/en/latest/user-guide/install/>`_ for detailed steps on how to install Conda for your operating system.

General installation
--------------------

If you are going to use ProteusLib's functionality, but *do not* plan to contribute to ProteusLib's codebase, choose this option.

#. Create a Conda environment (in this example, named ``proteuslib``) where ProteusLib and its runtime dependencies will be installed:

	.. code-block:: shell

		conda create --name proteuslib --yes python=3.8 pip=21.1

#. Activate the ``proteuslib`` environment:

	.. code-block:: shell

		conda activate proteuslib
	
	To verify that the correct environment is active, run the following command:

	.. code-block:: shell

		python -c "import sys; print(sys.executable)"
	
	If the environment was activated correctly, its name should be contained in the path displayed by the above command.

	.. important:: The ``conda activate`` command described above must be run each time a new terminal/console session is started.

#. Install ProteusLib using ``pip``:

	.. code-block:: shell

		pip install proteuslib

#. To verify that the installation was successful, open a Python interpreter and try importing some of ProteusLib's modules, e.g.:

	.. code-block:: shell

		python
		>>> from proteuslib.unit_models import *

.. _install-idaes-ext:

Installing solvers distributed through IDAES Extensions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After installing ProteusLib, the ``idaes get-extensions`` command can be used to automatically install the solvers distributed as part of the IDAES Extensions.

.. important:: Depending on your operating system, additional steps might be needed. For more information, refer to the `IDAES installation guide <https://idaes-pse.readthedocs.io/en/stable/getting_started/index.html>`_.

From the same environment where ProteusLib was installed, run:

    .. code-block:: shell

        idaes get-extensions

.. note:: Typically, the ``idaes get-extensions`` command only needs to be run once for each system, as it will install the required files into a common, system-wide location.


For ProteusLib developers
-------------------------

If you plan to contribute to ProteusLib's codebase, choose this option.

.. note:: Typically, *contributing to ProteusLib* will involve opening a Pull Request (PR) in ProteusLib's repository. For more information, refer to :ref:`developer-guide`.

#. Create a Conda environment (in this example, named ``proteuslib-dev``) where ProteusLib and all dependendencies needed for development will be installed, then activate it:

	.. code-block:: shell

		conda create --name proteuslib-dev --yes python=3.8 pip=21.1 && conda activate proteuslib-dev

	.. note:: For more information about using Conda environments, refer to the ":ref:`about-conda`" section above.

#. Clone the ProteusLib repository to your local development machine using ``git clone``, then enter the newly created ``proteuslib`` subdirectory:

	.. code-block:: shell

		git clone https://github.com/nawi-hub/proteuslib && cd proteuslib

#. Install ProteusLib and the development dependencies using ``pip`` and the ``requirements-dev.txt`` file:

	.. code-block:: shell

		pip install -r requirements-dev.txt

#. If needed, follow the steps described in the ":ref:`install-idaes-ext`" section above to install solvers distributed through IDAES Extensions.

#. To verify that the installation was successful, try running the ProteusLib test suite using ``pytest``:

	.. code-block:: shell

		pytest

#. To view/change the generated documentation, see the :ref:`documentation-mini-guide` section

Installing in existing development environments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When either the ``proteuslib`` package or one of its dependencies are installed, it should be possible to update those packages within an existing developer environment.

.. important:: In case of any issue or unexpected behavior when updating an existing environment,
    first try to see if the issues are solved if a freshly created environment is used instead.

#. Activate the environment, if not already active:

    .. code-block:: shell

        conda activate proteuslib-dev

#. Enter the directory where your local clone of the ProteusLib repository is located, and pull the latest changes using ``git pull``:

    .. code-block:: shell
        
        cd /path/to/your/clone
        git pull

#. Uninstall the version of ``proteuslib`` that's currently installed in the environment:

    .. code-block:: shell

        pip uninstall proteuslib

#. Run the ``pip install`` command targeting the ``requirements-dev.txt`` file.

    .. code-block:: shell

        pip --no-cache-dir install -r requirements-dev.txt

    .. note:: The ``--no-cache-dir`` flag is used to ensure that existing packages are not erroneously reused by pip,
        which would cause the wrong (outdated) version to be present in the environment after installation.


.. _documentation-mini-guide:

Documentation for developers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The project documentation is created and updated using the `Sphinx documentation tool <https://www.sphinx-doc.org>`_.
This tool generates nice, indexed, HTML webpages --- like this one --- from text files in the "docs" directory.
The documentation will include the docstrings you put on your modules, classes, methods, and functions as well
as additional documentation in text files in the "docs" directory. The project is set up so that Sphinx documentation
is generated automatically online for new releases. This section describes how to do this same documentation
generation locally in your development environment so you can preview what will be shown to the users.

.. _documentation-mini-guide-gen:

Generating the documentation
++++++++++++++++++++++++++++

To generate a local copy of the documentation for the first time, follow these steps:

1. Change directory to the "docs" subdirectory

2. Generate the tree of API documentation with "sphinx-apidoc". For convenience, a script has been
   provided that has all the required options.

   * On Windows, run ``.\apidoc.bat``

   * On Linux/OSX run ``./apidoc.sh``

3. Generate the HTML with Sphinx.

   * On Windows, run ``.\make html``

   * On Linux/OSX run ``make html``

After these steps are complete, you should be able to preview the HTML documentation by opening the file
located at "_build/html/index.html" in a web browser. To see the tree of API documentation that is generated
automatically from the source code, browse to the "Technical Reference" page and click on the "Modules" link at the
bottom.

.. _documentation-mini-guide-update:

Updating the documentation
++++++++++++++++++++++++++

If you make changes in your code's docstrings that you want to see reflected in the generated documentation,
you need to re-generate the API documentation using "sphinx-apidoc". To do this, simply re-run the command
given in step 2 of :ref:`documentation-mini-guide-gen`.

If you edited some documentation directly, i.e. created or modified a text file with extension `.rst`, then you
don't need to run the previous command. Regardless, you will next need to update the documentation with the
Sphinx build command given in step 3 of :ref:`documentation-mini-guide-gen`.

.. note:: The files under "docs/apidoc" are tracked in Git, otherwise they would not be available to the
ReadTheDocs builder (that doesn't know about sphinx-apidoc, strangely). Please remember to commit and push
them along with the changes in the source code.

Documenting your modules
++++++++++++++++++++++++
Full documentation for modules should be placed in the appropriate subfolder --- e.g., `property_models` or
`unit_models` --- of the `docs/technical_reference` section (and folder). See `docs/technical_reference/unit_modles/reverse_osmosis_0D.rst`
for an example.

Note that at the bottom of the file you should add the ``.. automodule::`` directive that will insert the
documentation for your module as generated from the source code (and docstrings). This generally looks like this:

    .. autodoc:: proteuslib.<package_name>.<module_name>
        :members:
        :noindex:

The ``:members:`` option says to include all the classes, functions, etc. in the module. It is important to add
the ``:noindex:`` option, otherwise Sphinx will try to generate an index entry that conflicts with the
entry that was created by the API docs (step 2 of :ref:`documentation-mini-guide-gen`), which would result
in warnings and failed builds for ReadTheDocs and the tests.