Ion Exchange (ZO)
=================

Model Type
----------
This unit model is formulated as a **single-input, double-output** model form.
See documentation for :ref:`single-input, double-output Helper Methods<sido_methods>`.

Electricity Consumption
-----------------------
Electricity consumption is calculated using the **pump_electricity** helper function.
See documentation for :ref:`Helper Methods for Electricity Demand<electricity_methods>`.

Costing Method
--------------
Costing is calculated using the :py:meth:`~watertap.unit_models.zero_order.ion_exchange_zo.IonExchangeZOData.cost_ion_exchange` method.
For full details on costing, see documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Flowrate of NaCl addition", "NaCl_flowrate", ":math:`kg/s`"
   "Dosage of NaCl addition", "NaCl_dose", ":math:`kg/m^3`"
   "Replacement rate of ion exchange resin", "resin_demand", ":math:`kg/s`"
   "Resin replacement as a function of flow", "resin_replacement", ":math:`kg/m^3`"

Additional Constraints
----------------------

.. csv-table::
   :header: "Description", "Constraint Name"

   "None", "NaCl_constraint"
   "None", "resin_constraint"

.. index::
   pair: watertap.unit_models.zero_order.ion_exchange_zo;ion_exchange_zo

.. currentmodule:: watertap.unit_models.zero_order.ion_exchange_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.ion_exchange_zo
    :members:
    :noindex:
