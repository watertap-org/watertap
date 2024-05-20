Nanofiltration Costing Method
==============================

Costing Method Parameters
+++++++++++++++++++++++++

The following parameters are constructed for the unit on the FlowsheetCostingBlock (e.g., `m.fs.costing.nanofiltration`) when applying the `cost_nanofiltration` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units"

   "Membrane replacement factor (fraction of membrane replaced/year)", ":math:`f_{mem,\, replace}`", "factor_membrane_replacement", "0.2", ":math:`\text{yr}^{-1}`"
   "Membrane unit cost", ":math:`C_{mem}`", "membrane_unit_cost", "15", ":math:`\text{USD}_{2018}\text{/m}^2`"

Costing Method Variables
++++++++++++++++++++++++

The following variables are constructed on the unit block (e.g., m.fs.unit.costing) when applying the `cost_nanofiltration` costing method in the ``watertap_costing_package``:

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

There are no operating costs unique to the nanofiltration unit.

 
Code Documentation
------------------

* :mod:`watertap.costing.unit_models.nanofiltration`
