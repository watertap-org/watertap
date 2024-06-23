Generic desalination unit
=================================

.. index::
   pair: watertap.unit_models.generic_desalter;generic_desalter

Introduction
------------
This is a generic desalination unit designed to separate water from other stream components. The user must specify target water recovery 
The unit assumes full complete separation of ions from the liquid stream (100 % rejection). 

Degrees of Freedom
------------------
For inlet conditions:
    * temperature
    * pressure
    * component mass compositions

Unit overall operaiton:
    * water recovery 
    * solids concentration in brine (User must provide tracked solids list)
    * water content of the brine (User must provide tracked solids list)

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "time", ":math:`t`", "[0]"
   "phases", ":math:`p`", "['Liq']"
   "components", ":math:`j`", "['H2O', solutes]"
   "tracked_solids_list", ":math:`tds`", "['TDS','X', etc.]"

Variables
----------

.. csv-table::
   :header: "Description", "Variable Name", "Index", "Units"
   
   "Water recovery", "water_recovery", "None", "%"
   "Brine water mass fraction", "brine_water_mass_percent", "None", "%"
   "Brine solids concentration", "brine_solids_concentration", "None", "kg/:math:`\text{m}^3`"

Class Documentation
-------------------

* :mod:`watertap.unit_models.generic_desalter`