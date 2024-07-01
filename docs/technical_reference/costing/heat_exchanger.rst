Heat Exchanger Costing Method
=============================

Costing Method Parameters
+++++++++++++++++++++++++

The following parameters are constructed for the unit on the FlowsheetCostingBlock (e.g., `m.fs.costing.heat_exchanger`) when applying the `cost_heat_exchanger` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units"

   "Heat exchanger unit cost", ":math:`C_{hx}`", "``unit_cost``", "300", ":math:`\text{USD}_{2020}`"
   "Material factor cost", ":math:`f_{m}`", "``material_factor_cost``", "1", ":math:`\text{dimensionless}`"

Costing Method Variables
++++++++++++++++++++++++

The following variables are used on the unit block (e.g., m.fs.unit.costing) when applying the `cost_heat_exchanger` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Heat exchanger area", ":math:`A_{hx}`", "``area``", "None", ":math:`\text{m}^2`"

Capital Cost Calculations
+++++++++++++++++++++++++

The capital cost is dependent on the heat exchanger area, :math:`A`, and the material factor, :math:`f_{m}`, as shown in the equation below.

    .. math::

        C_{cap, tot} = C_{hx} \cdot f_{m} \cdot A_{hx}

Operating Cost Calculations
+++++++++++++++++++++++++++

The operating cost includes the cost of steam when the heat exchanger is used as a steam heater. 
The steam consumption operating cost can be calculated as:

.. math::

    C_{op, tot} = C_{steam} \cdot \dot{m}_{steam}



Code Documentation
------------------

* :mod:`watertap.costing.unit_models.heat_exchanger`
