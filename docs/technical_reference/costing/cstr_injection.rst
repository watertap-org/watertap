Completely Stirred Tank Reactor w/Injection Stream Costing Method
==================================================================

Costing Method Parameters
+++++++++++++++++++++++++

The following parameters are constructed for the unit on the FlowsheetCostingBlock (e.g., `m.fs.costing.cstr_injection`) when applying the `cost_cstr_injection` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units"

   "Capital cost A parameter :math:`^1`", ":math:`A`", "capital_a_parameter", "1.1141e3", ":math:`\text{USD}_{1997}`"
   "Capital cost B parameter :math:`^1`", ":math:`B`", "capital_b_parameter", "0.6", ":math:`\text{dimensionless}`"

Costing Method Variables
++++++++++++++++++++++++

There are no costing method variables unique to the CSTR with injection.

Capital Cost Calculations
+++++++++++++++++++++++++

Capital cost is dependent upon the unit's volume, :math:`V`, as shown in the equation below.

    .. math::

        C_{cap,tot} = A * V^{B}

 
Operating Cost Calculations
+++++++++++++++++++++++++++

Electricity :math:`C_{elec}` is a variable operating cost based on the energy intensity :math:`E` of the unit process
(electricity consumption for the CSTR with injection), electricity price :math:`P`, electricity flow :math:`Q`, and the plant
utilization factor :math:`f_{util}`. The annual electricity costs are calculated as:

    .. math::

        C_{op, tot} = C_{elec} = E Q f_{util} P

 
Code Documentation
------------------

* :mod:`watertap.costing.unit_models.cstr_injection`

References
----------
[1] Eberle, Annika, Irina Tsiryapkina, Steve Peterson, Laura Vimmerstedt, Dylan Hettinger,
and Daniel Inman. 2020. An Overview of the Waste-to-Energy System Simulation
(WESyS) Model. Golden, CO: National Renewable Energy Laboratory.
NREL/TP-6A20-77166. https://www.nrel.gov/docs/fy21osti/77166.pdf.