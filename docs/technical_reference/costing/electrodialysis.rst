Electrodialysis Costing Method
===============================

Costing Method Parameters
+++++++++++++++++++++++++

The following parameters are constructed for the unit on the FlowsheetCostingBlock (e.g., `m.fs.costing.electrodialysis`) when applying the `cost_electrodialysis` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units"

   "Membrane unit cost", ":math:`C_{mem}`", "``membrane_capital_cost``", "160", ":math:`\text{USD}_{2018}\text{/m}^2`"
   "Membrane replacement factor (fraction of membrane replaced/year)", ":math:`f_{mem,\, replace}`", "``factor_membrane_replacement``", "0.2", ":math:`\text{yr}^{-1}`"
   "Electrode unit cost", ":math:`C_{elec}`", "``stack_electrode_capital_cost``", "2100", ":math:`\text{USD}_{2018}\text{/m}^2`"
   "Electrode replacement factor (fraction of electrode replaced/year)", ":math:`f_{elec,\, replace}`", "``factor_stack_electrode_replacement``", "0.2", ":math:`\text{yr}^{-1}`"
   "Rectifier cost coefficient, a", ":math:`a_{rectifier}`", "``rectifier_cost_coeff[0]``", "508.6", ":math:`\text{dimensionless}`"
   "Rectifier cost coefficient, b", ":math:`b_{rectifier}`", "``rectifier_cost_coeff[1]``", "2810", ":math:`\text{dimensionless}`"
   "AC to DC conversion efficiency", ":math:`\eta_{AC-DC}`", "``ac_dc_conversion_efficiency``", "0.9", ":math:`\text{dimensionless}`"

Costing Method Variables
++++++++++++++++++++++++

The following variables are constructed on the unit block (e.g., m.fs.unit.costing) when applying the `cost_electrodialysis` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Cell pair number in a stack", ":math:`n_{pair}`", "``cell_pair_num``", "None", ":math:`\text{dimensionless}`"
   "Cell width", ":math:`w`", "``cell_width``", "None", ":math:`\text{m}`"
   "Cell length", ":math:`l`", "``cell_length``", "None", ":math:`\text{m}`"
   "Power", ":math:`P`", "``power``", "[t]", ":math:`\text{kW}`"

Capital Cost Calculations
+++++++++++++++++++++++++

Capital cost is dependent upon the cell pair number, cell width, and cell length, as shown in the equations below.

If the system has a rectifier (``has_rectifier=True``):

    .. math::

        C_{rectifier} = b_{rectifier} + (a_{rectifier} * P / \eta_{AC-DC})

        C_{cap,tot} = (C_{mem} * 2 * n_{pair} * w * l) + (C_{elec} * 2 * w * l) + C_{rectifier}

Otherwise:

    .. math::

        C_{cap,tot} = (C_{mem} * 2 * n_{pair} * w * l) + (C_{elec} * 2 * w * l)

 
Operating Cost Calculations
+++++++++++++++++++++++++++

Electricity :math:`C_{elec}` is a variable operating cost based on the energy intensity :math:`E` of the unit process
(power for electrodialysis), electricity price :math:`P`, electricity flow :math:`Q`, and the plant
utilization factor :math:`f_{util}`. If a rectifier is included, the energy intensity is the power required for the rectifier.
The annual electricity costs are calculated as:

    .. math::

        C_{op, tot} = C_{elec} = E Q f_{util} P

 
Code Documentation
------------------

* :mod:`watertap.costing.unit_models.electrodialysis`
