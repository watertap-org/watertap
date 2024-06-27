Photothermal Membrane (ZO)
==========================

Model Type
----------
This unit model is formulated as a **single-input, double-output** model form.
See documentation for :ref:`single-input, double-output Helper Methods<sido_methods>`.

Electricity Consumption
-----------------------
The constraint used to calculate energy consumption is described in the Additional Constraints section below. More details can be found in the unit model class.

Costing Method
--------------
Costing is calculated using the :py:meth:`~watertap.unit_models.zero_order.photothermal_membrane_zo.PhotothermalMembraneData.cost_photothermal_membrane` method.
For full details on costing, see documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Water flux through membrane", "water_flux", ":math:`kg/h/m^2`"
   "Membrane area", "membrane_area", ":math:`m^2`"

Additional Constraints
----------------------

.. csv-table::
   :header: "Description", "Constraint Name"

   "Constraint for water flux.", "wat_flux"

.. index::
   pair: watertap.unit_models.zero_order.photothermal_membrane_zo;photothermal_membrane_zo

.. currentmodule:: watertap.unit_models.zero_order.photothermal_membrane_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.photothermal_membrane_zo
    :members:
    :noindex:
