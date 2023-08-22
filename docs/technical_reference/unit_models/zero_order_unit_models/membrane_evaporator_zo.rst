Membrane Evaporator (ZO)
========================

Model Type
----------
This unit model is formulated as a **single-input, double-output** model form.
See documentation for :ref:`single-input, double-output Helper Methods<sido_methods>`.

Electricity Consumption
-----------------------
Electricity consumption is calculated using the **constant_intensity** helper function.
See documentation for :ref:`Helper Methods for Electricity Demand<electricity_methods>`.

Costing Method
--------------
Costing is calculated using the :py:meth:`~watertap.unit_models.zero_order.membrane_evaporator_zo.MembraneEvaporatorData.cost_membrane_evaporator` method.
For full details on costing, see documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Water flux through membrane", "water_flux", ":math:`m/h`"
   "Membrane area", "membrane_area", ":math:`m^2`"

Additional Constraints
----------------------

.. csv-table::
   :header: "Description", "Constraint Name"

   "Constraint for water flux.", "wat_flux"

.. index::
   pair: watertap.unit_models.zero_order.membrane_evaporator_zo;membrane_evaporator_zo

.. currentmodule:: watertap.unit_models.zero_order.membrane_evaporator_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.membrane_evaporator_zo
    :members:
    :noindex:
