Evaporator Costing Method
=========================

Costing Method Parameters
+++++++++++++++++++++++++

The following parameters are constructed for the unit on the FlowsheetCostingBlock (e.g., `m.fs.costing.evaporator`) when applying the `cost_evaporator` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units"

   "Evaporator unit cost", ":math:`C_{evap}`", "``unit_cost``", "1000", ":math:`\text{USD}_{2020}`"
   "Material factor cost", ":math:`f_{m}`", "``material_factor_cost``", "1", ":math:`\text{dimensionless}`"

Costing Method Variables
++++++++++++++++++++++++

The following variables are used on the unit block (e.g., m.fs.unit.costing) when applying the `cost_evaporator` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Evaporator area", ":math:`A_{evap}`", "``area``", "None", ":math:`\text{m}^2`"

Capital Cost Calculations
+++++++++++++++++++++++++

The capital cost is dependent on the evaporator area, :math:`A`, and the material factor, :math:`f_{m}`, as shown in the equation below.

    .. math::

        C_{cap, tot} = C_{evap} \cdot f_{m} \cdot A_{evap}


Operating Cost Calculations
+++++++++++++++++++++++++++

There are no unique operating costs specific to the evaporator unit.

Code Documentation
------------------

* :mod:`watertap.costing.unit_models.evaporator`

