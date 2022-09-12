Feed  (ZO)
==========

Model Type
----------
This unit model is formulated as a **None** model form.
The Feed (ZO) block for zero-order flowsheets contains methods for getting concentration data from the database, and it has been created to work with the zero-order property package.

Electricity Consumption
-----------------------
This unit does not include energy consumption.

Costing Method
--------------
This unit does not include costing.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Volumetric flowrate in feed", "flow_vol", ":math:`m^3/s`"
   "Component mass concentrations", "conc_mass_comp", ":math:`kg/m^3`"

Additional Constraints
----------------------

.. csv-table::
   :header: "Description", "Constraint Name"

   "Volumetric flowrate of the feed", "flow_vol_constraint"
   "Component mass concentrations", "conc_mass_constraint"

.. index::
   pair: watertap.unit_models.zero_order.feed_zo;feed_zo

.. currentmodule:: watertap.unit_models.zero_order.feed_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.feed_zo
    :members:
    :noindex:
