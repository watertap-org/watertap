How to use apparent and true chemical species
=============================================

In ProteusLib, most all chemical processes simulated will be considered "true"
species, i.e., species that actually exist in an aqueous solution. However, their
may be times when you want to have your system report the "apparent" species
between your inlet and/or outlet ports of a flowsheet. To do this, you need
to understand how to define a species in a given configuration file as "apparent",
which is a special component type in the **thermo-properties** configuration dictionary.

For more information on creating a **thermo-properties** configuration dictionary,
see 'How to setup simple chemistry'.


What you need to update in the thermo-properties configuration dictionary
-------------------------------------------------------------------------

1. Add the ``"Apparent"`` and ``"StateIndex"`` objects under your import statements
2. Specify specific species as ``"Apparent"`` within its component definition
3. Give a list of "true" species it will dissociate into
4. Define a "state_components" argument in the config and specify ``"StateIndex.true"``


Add the following import statements
-----------------------------------

.. testsetup::

   # quiet idaes logs
   import idaes.logger as idaeslogger
   idaeslogger.getLogger('ideas.core').setLevel('CRITICAL')
   idaeslogger.getLogger('idaes.init').setLevel('CRITICAL')

.. testcode::

    # Importing the component types
    from idaes.core.components import Apparent

    # Importing the state basis 
    from idaes.generic_models.properties.core.generic.generic_property import StateIndex
