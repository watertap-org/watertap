Reverse Osmosis Costing Method
===============================

Costing Method Parameters
+++++++++++++++++++++++++

The following parameters are constructed for the unit on the FlowsheetCostingBlock (e.g., `m.fs.costing.reverse_osmosis`) when applying the `cost_reverse_osmosis` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units"

   "**Standard RO**"
   "Membrane replacement factor (fraction of membrane replaced/year)", ":math:`f_{mem,\, replace}`", "factor_membrane_replacement", "0.2", ":math:`\text{yr}^{-1}`"
   "Membrane unit cost", ":math:`C_{mem}`", "membrane_unit_cost", "30", ":math:`\text{USD}_{2012}\text{/m}^2`"

   "**High-pressure RO**"
   "Membrane replacement factor (fraction of membrane replaced/year)", ":math:`f_{mem,\, replace}`", "factor_membrane_replacement", "0.2", ":math:`\text{yr}^{-1}`"
   "Membrane unit cost", ":math:`C_{mem}`", "membrane_unit_cost", "75", ":math:`\text{USD}_{2012}\text{/m}^2`"


Costing Method Variables
++++++++++++++++++++++++

The following variables are constructed on the unit block (e.g., m.fs.unit.costing) when applying the `cost_dewatering` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Membrane area", ":math:`A_{mem}`", "area", "None", ":math:`\text{m}^2`"

Capital Cost Calculations
+++++++++++++++++++++++++

Capital cost is dependent upon the RO membrane area, :math:`A_{mem}`, as shown in the equations below.

    .. math::

        C_{cap,tot} = A_{mem} * C_{mem}

 
Operating Cost Calculations
+++++++++++++++++++++++++++

There are no operating costs unique to the RO unit.

 
Code Documentation
------------------

* :mod:`watertap.costing.unit_models.reverse_osmosis`

References
----------
Aim to include at least one reference in most cases, but delete this section if no references used for cost relationships/default values