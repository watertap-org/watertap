UV Reactor (ZO)
===============

Model Type
----------
This unit model is formulated as a single-input, single-output model form.
See documentation for :ref:`single-input, single-output Helper Methods<siso_methods>`.

Electricity Consumption
-----------------------
The constraint used to calculate energy consumption is described in the Additional Constraints section below. More details can be found in the unit model class.

Costing Method
--------------
Costing is calculated using the cost_uv method in the zero-order costing package.
See documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Reduced equivalent dosage", "uv_reduced_equivalent_dose", "mJ/cm**2"
   "UV transmittance of solution at UV reactor inlet", "uv_transmittance_in", "None"

Additional Constraints
----------------------

.. csv-table::
   :header: "Description", "Constraint Name"


.. index::
   pair: watertap.unit_models.zero_order.uv_zo;uv_zo

.. currentmodule:: watertap.unit_models.zero_order.uv_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.uv_zo
    :members:
    :noindex:
