Completely Stirred Tank Reactor Costing Method
===============================================

Costing Method Parameters
+++++++++++++++++++++++++

The following parameters are constructed for the unit on the FlowsheetCostingBlock (e.g., `m.fs.costing.cstr`) when applying the `cost_cstr` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units"

   "Capital cost A parameter :math:`^1`", ":math:`A`", "capital_a_parameter", "1246.1", ":math:`\text{USD}_{1990}`"
   "Capital cost B parameter :math:`^1`", ":math:`B`", "capital_b_parameter", "0.71", ":math:`\text{dimensionless}`"

Costing Method Variables
++++++++++++++++++++++++

There are no costing method variables unique to the CSTR.

Capital Cost Calculations
+++++++++++++++++++++++++

Capital cost is dependent upon the unit's volume, :math:`V`, as shown in the equation below.

    .. math::

        C_{cap,tot} = A * V^{B}

 
Operating Cost Calculations
+++++++++++++++++++++++++++

There are no operating costs unique to the CSTR.
 
Code Documentation
------------------

* :mod:`watertap.costing.unit_models.cstr`

References
----------
Aim to include at least one reference in most cases, but delete this section if no references used for cost relationships/default values