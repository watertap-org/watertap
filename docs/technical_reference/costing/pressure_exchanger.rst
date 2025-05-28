Pressure Exchanger Costing Method
==================================

Costing Method Parameters
+++++++++++++++++++++++++

The following parameters are constructed for the unit on the FlowsheetCostingBlock (e.g., ``m.fs.costing.pressure_exchanger``) when applying the ``cost_pressure_exchanger`` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units"

   "Pressure exchanger unit cost", ":math:`C_{PX}`", "``pressure_exchanger_cost``", "535", ":math:`\text{USD}_{2018}\text{/m}^3\text{/hr}`"

Costing Method Variables
++++++++++++++++++++++++

The following variables are constructed on the unit block (e.g., ``m.fs.unit.costing``) when applying the ``cost_pressure_exchanger`` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Inlet volumetric flow rate", ":math:`Q_{in}`", "``flow_in``", "[t]", ":math:`\text{m}^3\text{/hr}`"

Capital Cost Calculations
+++++++++++++++++++++++++

Capital cost is dependent upon the volumetric flow rate, :math:`Q_{in}`, as shown in the equation below.

    .. math::

        C_{cap,tot} = C_{PX} * Q_{in}

 
Operating Cost Calculations
+++++++++++++++++++++++++++

There are no operating costs unique to the pressure exchanger unit.

 
Code Documentation
------------------

* :mod:`watertap.costing.unit_models.pressure_exchanger`
