How to use apparent and true chemical species
=============================================

:ref: `How to setup simple chemistry`

In ProteusLib, most all chemical processes simulated will be considered "true"
species, i.e., species that actually exist in an aqueous solution. However, their
may be times when you want to have your system report the "apparent" species
between your inlet and/or outlet ports of a flowsheet. To do this, you need
to understand how to define a species in a given configuration file as "apparent",
which is a special component type in the **thermo-properties** configuration dictionary.
