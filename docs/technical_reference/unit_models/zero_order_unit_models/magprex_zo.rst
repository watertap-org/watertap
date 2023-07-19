Magprex (ZO)
============

Model Type
----------
This unit model is formulated as a **reactive single-inlet, double-outlet** model form.
See documentation for :ref:`reactive single-inlet, double-outlet Helper Methods<sidor_methods>`.

Electricity Consumption
-----------------------
Electricity consumption is calculated using the **constant_intensity** helper function.
See documentation for :ref:`Helper Methods for Electricity Demand<electricity_methods>`.

Costing Method
--------------
Costing is calculated using the :py:meth:`~watertap.unit_models.zero_order.magprex_zo.MagprexZOData.cost_magprex` method.
For full details on costing, see documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Dosage of magnesium chloride per phosphates", "magnesium_chloride_dosage", ":math:`dimensionless`"
   "Magnesium chloride flowrate", "MgCl2_flowrate", ":math:`kg/h`"

Additional Constraints
----------------------

.. csv-table::
   :header: "Description", "Constraint Name"

   "Constraint for magnesium chloride demand based on sludge flowrate.", "MgCl2_demand"

.. index::
   pair: watertap.unit_models.zero_order.magprex_zo;magprex_zo

.. currentmodule:: watertap.unit_models.zero_order.magprex_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.magprex_zo
    :members:
    :noindex:
