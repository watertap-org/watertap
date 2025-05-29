.. _generic_separator:

Generic Separator Unit
======================

.. code-block:: python

   from watertap.unit_models.generic_separation import GenericSeparation

.. index::
   pair: watertap.unit_models.generic_separation;generic_separation

Introduction
------------
This is a generic separation unit for target solute intended to imitate a selective separation process (e.g. softening, ion exchange, etc.) 
The unit assumes that target separation can be performed allowing user to configure percent of constituent to move from 
feed stream into product stream. 

Degrees of Freedom
------------------
For inlet conditions:
    * temperature
    * pressure
    * component mass compositions

Unit overall operation:
    * percent of component to remove from feed stream 
    * amount of chemical to add for operation 
    
Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Phases", ":math:`p`", "['Liq']"
   "Components", ":math:`j`", "['H2O', solutes]"
   
Variables
----------

.. csv-table::
   :header: "Description", "Variable Name", "Index", "Units"
   
   "Removal % of target component", "component_removal_percent", "j", "%"
   "Chemical dose", "additive_dose", "None", "mg/L"
   "Chemical mass flow", "additive_mass_flow", "None", "kg/s"

Class Documentation
-------------------

* :mod:`watertap.unit_models.generic_separation`