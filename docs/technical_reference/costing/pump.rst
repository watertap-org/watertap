Pump Costing Method
====================

Costing Method Parameters
+++++++++++++++++++++++++

The following parameters are constructed for the unit on the FlowsheetCostingBlock (e.g., `m.fs.costing.pump`) when applying the `cost_pump` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units"

   "**High-pressure pump**"
   "Pump unit cost", ":math:`C_{pump}`", "``cost``", "1.908", ":math:`\text{USD}_{2018}\text{/W}`"

   "**Low-pressure pump**"
   "Pump unit cost", ":math:`C_{pump}`", "``cost``", "889", ":math:`\text{USD}_{2018}\text{/L/s}`"

Costing Method Variables
++++++++++++++++++++++++

The following variables are constructed on the unit block (e.g., m.fs.unit.costing) when applying the `cost_pump` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "**High-pressure pump**"
   "Mechanical work", ":math:`W_{mech}`", "``work_mechanical``", "[t]", ":math:`\text{W}`"

   "**Low-pressure pump**"
   "Inlet volumetric flow rate", ":math:`Q_{in}`", "``flow_in``", "[t]", ":math:`\text{m}^3\text{/s}`"

Capital Cost Calculations
+++++++++++++++++++++++++

For the high pressure pump , capital cost is dependent upon the mechanical work, :math:`W_{mech}`, whereas the capital cost of
the low pressure pump is based on the volumetric flow rate :math:`Q_{in}`.

    .. math::

        C_{cap,high pressure} = C_{pump} * W_{mech}

    .. math::

        C_{cap,low pressure} = C_{pump} * Q_{in}

 
Operating Cost Calculations
+++++++++++++++++++++++++++

Electricity :math:`C_{elec}` is a variable operating cost based on the energy intensity :math:`E` of the unit process
(mechanical work for the pump), electricity price :math:`P`, electricity flow :math:`Q`, and the plant
utilization factor :math:`f_{util}`. The annual electricity costs are calculated as:

    .. math::

        C_{op, tot} = C_{elec} = E Q f_{util} P

 
Code Documentation
------------------

* :mod:`watertap.costing.unit_models.pump`
