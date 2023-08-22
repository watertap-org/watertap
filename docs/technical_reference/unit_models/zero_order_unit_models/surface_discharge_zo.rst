Surface Discharge (ZO)
======================

Model Type
----------
This unit model is formulated as a **pass-through** model form.
See documentation for :ref:`pass-through Helper Methods<pt_methods>`.

Electricity Consumption
-----------------------
Electricity consumption is calculated using the **pump_electricity** helper function.
See documentation for :ref:`Helper Methods for Electricity Demand<electricity_methods>`.

Costing Method
--------------
Costing is calculated using the :py:meth:`~watertap.unit_models.zero_order.surface_discharge_zo.SurfaceDischargeData.cost_surface_discharge` method.
For full details on costing, see documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Piping distance", "pipe_distance", ":math:`mi`"
   "Pipe diameter", "pipe_diameter", ":math:`in`"

.. index::
   pair: watertap.unit_models.zero_order.surface_discharge_zo;surface_discharge_zo

.. currentmodule:: watertap.unit_models.zero_order.surface_discharge_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.surface_discharge_zo
    :members:
    :noindex:
