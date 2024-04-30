Anaerobic Digester Costing Method
==================================

Costing Method Parameters
+++++++++++++++++++++++++

The following parameters are constructed for the unit on the FlowsheetCostingBlock (e.g., `m.fs.costing.anaerobic_digester`) when applying the `cost_anaerobic_digester` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units"

   "Capital cost A parameter", ":math:`A`", "capital_a_parameter", "19.3552312e6", ":math:`\text{USD}_{2012}`"
   "Capital cost B parameter", ":math:`B`", "capital_b_parameter", "0.6", ":math:`\text{dimensionless}`"
   "Reference flow", ":math:`Q_{ref}`", "reference_flow", "911.054", ":math:`\text{m}^3\text{/hr}`"

Costing Method Variables
++++++++++++++++++++++++

There are no costing method variables unique to the anaerobic digester

Capital Cost Calculations
+++++++++++++++++++++++++

Capital cost is dependent upon the anaerobic digester's influent flow rate, :math:`Q_{in}`, as shown in the equation below.

    .. math::

        C_{cap,tot} = A * \frac{Q_{in}}{Q_{ref}}**B

 
Operating Cost Calculations
+++++++++++++++++++++++++++

There are no operating costs unique to the anaerobic digester.

 
Code Documentation
------------------

* :mod:`watertap.costing.unit_models.anaerobic_digester`

References
----------
Aim to include at least one reference in most cases, but delete this section if no references used for cost relationships/default values