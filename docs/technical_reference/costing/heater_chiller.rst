Heater/Chiller Costing Method
=============================

Costing Method Parameters
+++++++++++++++++++++++++

The following parameters are constructed for the unit on the FlowsheetCostingBlock (e.g., `m.fs.costing.heater_chiller`) when applying the `cost_heater_chiller` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units"

   "**Electric Heater**"
   "Heater unit cost", ":math:`C_{heater}`", "``unit_cost``", "0.066", ":math:`\text{USD}_{2018}\text{/W}`"
   "Heat generation efficiency", ":math:`HE`", "``HE``", "0.99", ":math:`\text{dimensionless}`"

   "**Chiller**"
   "Chiller unit cost", ":math:`C_{chiller}`", "``unit_cost``", "0.2", ":math:`\text{USD}_{2018}\text{/W}`"
   "Coefficient of performance", ":math:`COP`", "``COP``", "7", ":math:`\text{dimensionless}`"

Costing Method Variables
++++++++++++++++++++++++

The following variables are used on the unit block (e.g., m.fs.unit.costing) when applying the `cost_heater_chiller` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Heat duty", ":math:`E`", "``heat_duty``", "None", ":math:`\text{W}`"

Capital Cost Calculations
+++++++++++++++++++++++++

For the **Electric Heater**:

The capital cost is dependent on the heat duty, :math:`E`, and the heat generation efficiency, :math:`HE`, as shown in the equation below.

    .. math::

        C_{cap, tot} = C_{heater} \cdot \frac{E}{HE}

For the **Chiller**:

The capital cost is dependent on the effective heat duty, :math:`E`, and the coefficient of performance, :math:`COP`, as shown in the equation below.

    .. math::

        C_{cap, tot} = C_{chiller} \cdot \frac{E}{COP}

Operating Cost Calculations
+++++++++++++++++++++++++++

Electricity :math:`C_{elec}` is a variable operating cost based on the energy intensity :math:`E` of the unit process
(heat duty for the heater or the chiller), electricity price :math:`P`, electricity flow :math:`Q`, and the plant
utilization factor :math:`f_{util}`. The annual electricity costs are calculated as:

    .. math::

        C_{op, tot} = C_{elec} = E Q f_{util} P

Code Documentation
------------------

* :mod:`watertap.costing.unit_models.heater_chiller`

References
----------

- Estimated from multiple sources.
