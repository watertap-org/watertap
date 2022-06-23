Chlorination  (ZO)
==================

Model Type
----------
This unit model is formulated as a single-input, single-output model form.
See documentation for :ref:`single-input, single-output Helper Methods<siso_methods>`.

Electricity Consumption
-----------------------
The constraint used to calculate energy consumption is described in the Additional Constraints section below. More details can be found in the unit model class.

Costing Method
--------------
Costing is calculated using the cost_chlorination method in the zero-order costing package.
See documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Initial chlorine demand", "initial_chlorine_demand", "mg/l"
   "Chlorine contact time", "contact_time", "hr"
   "CT value for chlorination", "concentration_time", "mg*min/l"
   "Chlorine decay rate", "chlorine_decay_rate", "mg/hr/l"
   "Chlorine dose", "chlorine_dose", "mg/l"

Additional Constraints
----------------------

.. csv-table::
   :header: "Description", "Constraint Name"

   "Chlorine dose constraint", "chlorine_dose_constraint"

.. index::
   pair: watertap.unit_models.zero_order.chlorination_zo;chlorination_zo

.. currentmodule:: watertap.unit_models.zero_order.chlorination_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.chlorination_zo
    :members:
    :noindex:
