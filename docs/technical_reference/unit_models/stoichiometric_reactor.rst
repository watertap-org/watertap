Stoichiometric reactor
======================

.. index::
   pair: watertap.unit_models.stoichiometric_reactor;stoichiometric_reactor

The main assumptions of the implemented model are as follows:

1) Single liquid phase only (solids are not explicitly modeled)
2) Steady state only
3) Has no volume
4) Reaction phenomena is not modeled (e.g. exothermic, endothermic, reaction rate, etc.)
   User must provide amount of reagent to add, and amount of precipitant that would form

Introduction
------------
This is a basic stoichiometric reactor designed to provide a simple mass balanced method for adding reactants and removing solids from a stream, for more comprehesive reactor models please use IDAES Stoichiometric reactor model (https://idaes-pse.readthedocs.io/en/latest/reference_guides/model_libraries/generic/unit_models/stoichiometric_reactor.html). 

The stoichiometric reactor is a basic unit operation designed to aid in modeling addition of reagent to a feed stream and removal of ions from a stream through precipitation. A basic example for using this model is the lime/soda ash softening process. 
The reactor can be configured to include only dissolution of a reagent or precipitation of specific species, and both as shown in figure below.

.. figure:: ../../_static/unit_models/stoichiometric_reactor.png
    :width: 600
    :align: center
    
    Figure 1. Schematic representation of a processes considered in stoichiometric reactor

Degrees of Freedom
------------------
The degrees of freedom in a stoichiometric reactor unit model are as follows:

For inlet conditions:
    * temperature
    * pressure
    * component mass compositions

If reagents are supplied:
   * reagent mass flow rate or reagent dose

If precipitants are supplied 
   * precipitant mass flow rate 
   * mass fraction of solids in precipitant (solid) waste stream

Model Structure and usage
-------------------------
The stoichiometric reactor uses control volumes to perform the dissolution reaction and precipitation reaction, while a an IDAES separator is used to separate precipitated solids from the feed stream. The model should be used with MCAS property package.
The user needs to specify how supplied reagent, and precipitant dissolve or precipitate out of the feed stream, using ions present in the feed. 

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "time", ":math:`t`", "[0]"
   "phases", ":math:`p`", "['Liq']"
   "components", ":math:`j`", "['H2O', solutes]"
   "reagents", ":math:`\text{reagents}`",[reagent]
   "precipitants", ":math:`\text{precipitants}`",[precipitants]

Variables
----------

.. csv-table::
   :header: "Description", "Variable Name", "Index", "Units"
   
   "Reagent dose", 'reagent_dose','[reagent]','kg/:math:`\text{m}^3`'
   "Reagent flow mass", 'flow_mass_reagent','[reagent]','kg/s'
   "Flow mass of precipitant",'flow_mass_precipitate',[precipitant],'kg/s'
   "Mass concentration of precipitant",'conc_mass_precipitate',[precipitant],'kg/:math:`\text{m}^3`'
   "Fraction of solids in waste stream",  "waste_mass_frac_precipitate", None, fraction
   