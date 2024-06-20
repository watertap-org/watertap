Compressor Costing Method
==========================

Costing Method Parameters
+++++++++++++++++++++++++

The following parameters are constructed for the unit on the FlowsheetCostingBlock (e.g., `m.fs.costing.compressor`) when applying the `cost_compressor` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units"

   "Compressor unit cost", ":math:`C_{comp}`", "``unit_cost``", "7364", ":math:`\text{USD}_{2001}`"
   "Compressor cost exponent", ":math:`e`", "``exponent``", "0.7", ":math:`\text{dimensionless}`"

Costing Method Variables
++++++++++++++++++++++++

The following variables are used on the unit block (e.g., m.fs.unit.costing) when applying the `cost_compressor` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Vapor flow rate", ":math:`\dot{m}_{\text{Vap, H2O}}`", "``flow_mass_phase_comp``", "None", ":math:`\text{kg/s}`"
   "Pressure ratio", ":math:`PR`", "``pressure_ratio``", "None", ":math:`\text{dimensionless}`"

Capital Cost Calculations
+++++++++++++++++++++++++

The capital cost is dependent on the flow rate of water vapor, :math:`\dot{m}_{\text{Vap, H2O}}`, and the pressure ratio, :math:`PR`, as shown in the equation below.

    .. math::

        C_{cap, tot} = C_{comp} \dot{m}_{\text{Vap, H2O}} PR \left( \frac{\eta}{1 - \eta} \right)^e

where
:math:`\eta` is the compressor efficiency

Operating Cost Calculations
+++++++++++++++++++++++++++

Electricity :math:`C_{elec}` is a variable operating cost based on the energy intensity :math:`E` of the unit process
(mechanical work for the compressor), electricity price :math:`P`, electricity flow :math:`Q`, and the plant
utilization factor :math:`f_{util}`. The annual electricity costs are calculated as:

    .. math::

        C_{op, tot} = C_{elec} = E Q f_{util} P

Code Documentation
------------------

* :mod:`watertap.costing.unit_models.compressor`

References
----------

El-Sayed, Y. M. (2001). Designing desalination systems for higher productivity. Desalination, 134, 129–58.
