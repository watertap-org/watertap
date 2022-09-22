Feed  (ZO)
==========

The Feed (ZO) model adds volumetric flowrate and mass concentration as variables in connection with the zero-order property model's state variable, mass flow rate. This allows the user to enter feed volumetric flow rate and concentrations without calculating the mass flowrates of each component in the feed.

Model Type
----------
The Feed (ZO) block for zero-order flowsheets contains methods for getting concentration data from the database, and it has been created to work with the zero-order property package.

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

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.feed_zo
    :members:
    :noindex:
