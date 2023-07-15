UV with Advanced Oxidation Processes (ZO)
=========================================

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
Costing is calculated using the :py:meth:`~watertap.unit_models.zero_order.uv_aop_zo.UVAOPZOData.cost_uv_aop` method.
For full details on costing, see documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Reduced equivalent dosage", "uv_reduced_equivalent_dose", ":math:`mJ/cm^2`"
   "UV transmittance of solution at UV reactor inlet", "uv_transmittance_in", ":math:`dimensionless`"
   "Oxidant dosage", "oxidant_dose", ":math:`mg/l`"
   "Mass flow rate of oxidant solution", "chemical_flow_mass", ":math:`kg/s`"

Additional Constraints
----------------------

.. csv-table::
   :header: "Description", "Constraint Name"

   "Chemical mass flow constraint", "chemical_flow_mass_constraint"

.. index::
   pair: watertap.unit_models.zero_order.uv_aop_zo;uv_aop_zo

.. currentmodule:: watertap.unit_models.zero_order.uv_aop_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.uv_aop_zo
    :members:
    :noindex:
