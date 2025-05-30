Thickener Costing Method
=========================

Costing Method Parameters
+++++++++++++++++++++++++

The following parameters are constructed for the unit on the FlowsheetCostingBlock (e.g., ``m.fs.costing.thickener``) when applying the ``cost_thickener`` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units"

   "Capital cost A parameter :math:`^1`", ":math:`A`", "``capital_a_parameter``", "4729.8", ":math:`\text{USD}_{2011}\text{/ft}`"
   "Capital cost B parameter :math:`^1`", ":math:`B`", "``capital_b_parameter``", "37068", ":math:`\text{USD}_{2007}`"

Costing Method Variables
++++++++++++++++++++++++

The following variables are constructed on the unit block (e.g., ``m.fs.unit.costing``) when applying the ``cost_thickener`` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Diameter", ":math:`d`", "diameter", "None", ":math:`\text{ft}`"

Capital Cost Calculations
+++++++++++++++++++++++++

Capital cost is dependent upon the thickener diameter, :math:`d`, as shown in the equation below.

    .. math::

        C_{cap,tot} = A * d + B

 
Operating Cost Calculations
+++++++++++++++++++++++++++

Electricity :math:`C_{elec}` is a variable operating cost based on the energy intensity :math:`E` of the unit process
(electricity consumption for the thickener), electricity price :math:`P`, electricity flow :math:`Q`, and the plant
utilization factor :math:`f_{util}`. The annual electricity costs are calculated as:

    .. math::

        C_{op, tot} = C_{elec} = E Q f_{util} P

 
Code Documentation
------------------

* :mod:`watertap.costing.unit_models.thickener`

References
----------
[1] W. McGivney, S. Kawamura, Cost estimating manual for water treatment facilities,
John Wiley & Sons, 2008. http://onlinelibrary.wiley.com/book/10.1002/9780470260036.