.. _generic_desalter:

Generic Desalter Unit
=====================

.. code-block:: python

   from watertap.unit_models.generic_desalter import GenericDesalter

.. index::
   pair: watertap.unit_models.generic_desalter;generic_desalter

Introduction
------------
This is a generic desalination unit designed to separate water from other stream components. The user must specify target water recovery.
The unit assumes complete separation of ions from the liquid stream (100 % rejection). 

Degrees of Freedom
------------------
For inlet conditions:
    * temperature
    * pressure
    * component mass compositions

Unit overall operation:
    * water recovery 
    * solids concentration in brine (User must provide tracked solids list)
    * water content of the brine (User must provide tracked solids list)

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Phases", ":math:`p`", "['Liq']"
   "Components", ":math:`j`", "['H2O', solutes]"
   "Tracked solids list", ":math:`tds`", "['TDS','X', etc.]"

Variables
----------

.. csv-table::
   :header: "Description", "Variable Name", "Index", "Units"
   
   "Water recovery", "water_recovery", "None", ":math:`\text{dimensionless}`"
   "Brine water mass fraction", "brine_water_mass_percent", "None", ":math:`\text{dimensionless}`"
   "Brine solids concentration", "brine_solids_concentration", "None", "kg/:math:`\text{m}^3`"

Class Documentation
-------------------

* :mod:`watertap.unit_models.generic_desalter`