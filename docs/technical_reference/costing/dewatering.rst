Dewatering Unit Costing Method
===============================

There are three classes of dewatering costing types (centrifuge, filter belt press, and filter plate press), each with their own parameters, variables,
and costing relationships. The default configuration is a centrifuge, so users must manually change the dewatering type
if they wish to invoke the costing method for filter belt press or filter plate press dewatering units.

Costing Method Parameters
+++++++++++++++++++++++++

The following parameters are constructed for the unit on the FlowsheetCostingBlock (e.g., `m.fs.costing.dewatering`) when applying the `cost_dewatering` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units"

   "**Centrifuge**"
   "Capital cost A parameter :math:`^1`", ":math:`A_{centrifuge}`", "``capital_a_parameter``", "328.03", ":math:`\text{USD}_{2007}\text{/gal/hr}`"
   "Capital cost B parameter :math:`^1`", ":math:`B_{centrifuge}`", "``capital_b_parameter``", "751295", ":math:`\text{USD}_{2007}`"

   "**Filter belt press**"
   "Capital cost A parameter :math:`^1`", ":math:`A_{fbp}`", "``capital_a_parameter``", "146.29", ":math:`\text{USD}_{2007}\text{/gal/hr}`"
   "Capital cost B parameter :math:`^1`", ":math:`B_{fbp}`", "``capital_b_parameter``", "433972", ":math:`\text{USD}_{2007}`"

   "**Filter plate press**"
   "Capital cost A parameter :math:`^1`", ":math:`A_{fpp}`", "``capital_a_parameter``", "102794", ":math:`\text{USD}_{2007}\text{/gal/hr}`"
   "Capital cost B parameter :math:`^1`", ":math:`B_{fpp}`", "``capital_b_parameter``", "0.4216", ":math:`\text{USD}_{2007}`"

Costing Method Variables
++++++++++++++++++++++++

The following variables are constructed on the unit block (e.g., m.fs.unit.costing) when applying the `cost_dewatering` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Inlet volumetric flow rate", ":math:`Q_{in}`", "``flow_in``", "[t]", ":math:`\text{gal/hr}`"

Capital Cost Calculations
+++++++++++++++++++++++++

Capital cost is dependent upon the volumetric flow rate, :math:`Q_{in}`, as shown in the equations below.

    .. math::

        C_{cap,centrifuge} = A_{centrifuge} * Q_{in} + B_{centrifuge}

    .. math::

        C_{cap,fbp} = A_{fbp} * Q_{in} + B_{fbp}

    .. math::

        C_{cap,fpp} = A_{fpp} * Q_{in} + B_{fpp}
 
Operating Cost Calculations
+++++++++++++++++++++++++++

There are no operating costs unique to the dewatering unit.

 
Code Documentation
------------------

* :mod:`watertap.costing.unit_models.dewatering`

References
----------
[1] W. McGivney, S. Kawamura, Cost estimating manual for water treatment facilities,
John Wiley & Sons, 2008. http://onlinelibrary.wiley.com/book/10.1002/9780470260036.