Anaerobic Digester Costing Method
==================================

Costing Method Parameters
+++++++++++++++++++++++++

The following parameters are constructed for the unit on the FlowsheetCostingBlock (e.g., `m.fs.costing.anaerobic_digester`) when applying the `cost_anaerobic_digester` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units"

   "Capital cost A parameter :math:`^1`", ":math:`A`", "``capital_a_parameter``", "19.3552312e6", ":math:`\text{USD}_{2012}`"
   "Capital cost B parameter :math:`^1`", ":math:`B`", "``capital_b_parameter``", "0.6", ":math:`\text{dimensionless}`"
   "Reference flow :math:`^1`", ":math:`F_r`", "``reference_flow``", "911.054", ":math:`m^3/h`"

Costing Method Variables
++++++++++++++++++++++++

There are no costing method variables unique to the anaerobic digester.

Capital Cost Calculations
+++++++++++++++++++++++++

Capital cost is dependent upon the unit's volumetric flowrate, :math:`F`, as shown in the equation below.

    .. math::

        C_{cap,tot} = A * F/F_r^{B}

 
Operating Cost Calculations
+++++++++++++++++++++++++++

Electricity :math:`C_{elec}` is a variable operating cost based on the energy intensity :math:`E` of the unit process
(electricity consumption for the anaerobic digester), electricity price :math:`P`, electricity flow :math:`Q`, and the plant
utilization factor :math:`f_{util}`. The annual electricity costs are calculated as:

    .. math::
        C_{op, tot} = C_{elec} = E Q f_{util} P
 
Code Documentation
------------------

* :mod:`watertap.costing.unit_models.anaerobic_digester`

References
----------
[1] Eberle, Annika, Irina Tsiryapkina, Steve Peterson, Laura Vimmerstedt, Dylan Hettinger,
and Daniel Inman. 2020. An Overview of the Waste-to-Energy System Simulation
(WESyS) Model. Golden, CO: National Renewable Energy Laboratory.
NREL/TP-6A20-77166. https://www.nrel.gov/docs/fy21osti/77166.pdf.