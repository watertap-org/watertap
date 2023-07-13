UV Reactor (ZO)
===============

Model Type
----------
This unit model is formulated as a **single-input, single-output** model form.
See documentation for :ref:`single-input, single-output Helper Methods<siso_methods>`.

Electricity Consumption
-----------------------
Electricity consumption is calculated using the **constant_intensity** helper function.
See documentation for :ref:`Helper Methods for Electricity Demand<electricity_methods>`.

Costing Method
--------------
Costing is calculated using the :py:meth:`~watertap.unit_models.zero_order.uv_zo.UVZOData.cost_uv` method.
For full details on costing, see documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Reduced equivalent dosage", "uv_reduced_equivalent_dose", ":math:`mJ/cm^2`"
   "UV transmittance of solution at UV reactor inlet", "uv_transmittance_in", ":math:`dimensionless`"

.. index::
   pair: watertap.unit_models.zero_order.uv_zo;uv_zo

.. currentmodule:: watertap.unit_models.zero_order.uv_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.uv_zo
    :members:
    :noindex:
