Clarifier Costing Method
=========================

There are three classes of clarifier costing types (circular, rectangular, and primary), each with their own parameters, variables,
and costing relationships. The default configuration is a circular clarifier, so users must manually change the clarifier type
if they wish to invoke the costing method for rectangular or primary clarifiers.

Costing Method Parameters
+++++++++++++++++++++++++

The following parameters are constructed for the unit on the FlowsheetCostingBlock (e.g., `m.fs.costing.clarifier`) when applying the `cost_clarifier` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units"

   "**Circular clarifier**"
   "Capital cost A parameter :math:`^1`", ":math:`A_{circular}`", "construction_a_parameter", "-6e-4", ":math:`\text{USD}_{2011}\text{/ft}^4`"
   "Capital cost B parameter :math:`^1`", ":math:`B_{circular}`", "construction_b_parameter", "98.952", ":math:`\text{USD}_{2011}\text{/ft}^2`"
   "Capital cost C parameter :math:`^1`", ":math:`C_{circular}`", "construction_c_parameter", "191806", ":math:`\text{USD}_{2011}`"

   "**Rectangular clarifier**"
   "Capital cost A parameter :math:`^1`", ":math:`A_{rectangular}`", "construction_a_parameter", "-2.9e-3", ":math:`\text{USD}_{2011}\text{/ft}^4`"
   "Capital cost B parameter :math:`^1`", ":math:`B_{rectangular}`", "construction_b_parameter", "169.19", ":math:`\text{USD}_{2011}\text{/ft}^2`"
   "Capital cost C parameter :math:`^1`", ":math:`C_{rectangular}`", "construction_c_parameter", "94365", ":math:`\text{USD}_{2011}`"

   "**Primary clarifier**"
   "Capital cost A parameter :math:`^2`", ":math:`A_{primary}`", "capital_a_parameter", "-2.9e-3", ":math:`\text{USD}_{2021}`"
   "Capital cost B parameter :math:`^2`", ":math:`B_{primary}`", "capital_b_parameter", "538746.398", ":math:`\text{dimensionless}`"

Costing Method Variables
++++++++++++++++++++++++

The following variables are constructed on the unit block (e.g., m.fs.unit.costing) when applying the `cost_clarifier` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "**Circular clarifier**"
   "Surface area", ":math:`A_{s}`", "surface_area", "None", ":math:`\text{ft}^2`"

   "**Rectangular clarifier**"
   "Surface area", ":math:`A_{s}`", "surface_area", "None", ":math:`\text{ft}^2`"

   "**Primary clarifier**"
   "Inlet volumetric flow rate", ":math:`Q_{in}`", "flow_in", "[t]", ":math:`\text{gal/day}`"

Capital Cost Calculations
+++++++++++++++++++++++++

For the circular and rectangular clarifiers, capital cost is dependent upon the surface area, :math:`A_{s}`, whereas the capital cost of
the primary clarifier is based on the volumetric flow rate :math:`Q_{in}`.

    .. math::

        C_{cap,circular} = A_{circular} * A_{s}^{2} + B_{circular} * A_{s} + C_{circular}

    .. math::

        C_{cap,rectangular} = A_{rectangular} * A_{s}^{2} + B_{rectangular} * A_{s} + C_{rectangular}

    .. math::

        C_{cap,primary} = A_{primary} * Q_{in}^{B_{primary}}

Operating Cost Calculations
+++++++++++++++++++++++++++

Electricity :math:`C_{elec}` is a variable operating cost based on the energy intensity :math:`E` of the unit process
(electricity consumption for the clarifier), electricity price :math:`P`, electricity flow :math:`Q`, and the plant
utilization factor :math:`f_util`. The annual electricity costs are calculated as:

    .. math::

        C_{op, tot} = C_{elec} = E Q f_{util} P

 
Code Documentation
------------------

* :mod:`watertap.costing.unit_models.clarifier`

References
----------
[1] Sharma, Jwala R., Mohammad Najafi, and Syed R. Qasim.
"Preliminary cost estimation models for construction, operation, and maintenance of water treatment plants."
Journal of Infrastructure Systems 19.4 (2013): 451-464.

[2] Byun, Jaewon, Maravelias, Christos.
Benchmark Model for Wastewater Treatment Using an Activated Sludge Process.
United States: N.p., 21 Jan, 2022. Web. doi: 10.7481/1844539.