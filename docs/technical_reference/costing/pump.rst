.. _pump_costing:

Pump Costing Method
====================

Costing Method Parameters
+++++++++++++++++++++++++

The following parameters are constructed for the unit on the FlowsheetCostingBlock (e.g., ``m.fs.costing.pump``) when applying the ``cost_pump`` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units"

   "**High-pressure pump** (cost method = ``cost_high_pressure_pump``)"
   "Pump unit cost", ":math:`C_{pump}`", "``cost``", "1.908", ":math:`\text{USD}_{2018}\text{/W}`"

   "**Low-pressure pump** (cost method = ``cost_low_pressure_pump``)"
   "Pump unit cost", ":math:`C_{pump}`", "``cost``", "889", ":math:`\text{USD}_{2018}\text{/L/s}`"

Costing Method Variables
++++++++++++++++++++++++

The following variables are constructed on the unit block (e.g., ``m.fs.unit.costing``) when applying the ``cost_pump`` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "**High-pressure pump**"
   "Mechanical work", ":math:`W_{mech}`", "``work_mechanical``", "[t]", ":math:`\text{W}`"

   "**Low-pressure pump**"
   "Inlet volumetric flow rate", ":math:`Q_{in}`", "``flow_in``", "[t]", ":math:`\text{m}^3\text{/s}`"

Capital Cost Calculations
+++++++++++++++++++++++++

For the high pressure pump, capital cost is dependent upon the mechanical work, :math:`W_{mech}`, whereas the capital cost of
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

References
----------

| T. V. Bartholomew, N. S. Siefert and M. S. Mauter (2018).
| Cost Optimization of Osmotically Assisted Reverse Osmosis.
| *Environ Sci Technol* 2018 Vol. 52 Issue 20 Pages 11813-11821.
| DOI: 10.1021/acs.est.8b02771

| A. Malek, M. N. A. Hawlader and J. C. Ho (1996).
| Design and economics of RO seawater desalination.
| *Desalination* 1996 Vol. 105 Issue 3 Pages 245-261.
| DOI: 10.1016/0011-9164(96)00081-1