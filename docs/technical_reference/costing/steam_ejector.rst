Steam Ejector Costing Method
============================

Costing Method Parameters
+++++++++++++++++++++++++

The following parameters are constructed for the unit on the FlowsheetCostingBlock (e.g., ``m.fs.costing.thickener``) when applying the ``cost_thickener`` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units"

   "Base cost coefficient for steam ejector", ":math:`C_{base}`", "``base_cost``", "1949", ":math:`\text{USD}_{2020}`"
   "Cost scaling exponent", ":math:`C_{exponent}`", "``cost_exponent``", "0.3", ":math:`\text{dimensionless}`"
   "Steam cost", ":math:`C_{steam}`", "``steam_cost``", "0.008", ":math:`\text{USD}_{2018}\text{/kg}`"

Costing Method Variables
++++++++++++++++++++++++

There are no costing method variables unique to the steam ejector.

Capital Cost Calculations
+++++++++++++++++++++++++

Capital cost is dependent upon the steam and entrained vapor mass flow, :math:`Q`, as shown in the equation below.

    .. math::

        C_{cap,tot} = C_{base} * (Q_{steam} + Q_{entrained vapor}^{C_{exponent}}

 
Operating Cost Calculations
+++++++++++++++++++++++++++

The only operating cost associated with this steam ejector model is the cost of steam,
which is a function of the steam cost (on a mass basis), :math:`C_{steam}`, and the
mass flow of steam, :math:`Q_{steam}`, as shown in the equation below:

    .. math::

        C_{op, tot} = C_{steam} * Q_{steam}

 
Code Documentation
------------------

* :mod:`watertap.costing.unit_models.steam_ejector`

