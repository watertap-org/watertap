Membrane Distillation Costing Method
=====================================

Costing Method Parameters
+++++++++++++++++++++++++

The following parameters are constructed for the unit on the FlowsheetCostingBlock (e.g., `m.fs.costing.membrane_distillation`) when applying the `cost_membrane_distillation` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units"

   "Membrane replacement factor", ":math:`f_{replace}`", "``factor_membrane_replacement``", "0.2", ":math:`\text{yr}^{-1}`"
   "Membrane cost", ":math:`C_{mem}`", "``membrane_unit_cost``", "56", ":math:`\text{USD}_{2018}\text{/m}^2`"

Costing Method Variables
++++++++++++++++++++++++

The following variables are constructed on the unit block (e.g., m.fs.unit.costing) when applying the `cost_membrane_distillation` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Membrane area", ":math:`A_{mem}`", "``area``", "None", ":math:`\text{m}^2`"

Capital Cost Calculations
+++++++++++++++++++++++++

The capital cost is dependent upon the membrane area, :math:`A_{mem}`, as shown in the equations below.

    .. math::

        C_{cap,tot} = A_{mem} * C_{mem}

Operating Cost Calculations
+++++++++++++++++++++++++++

The fixed operating cost is proportional to the membrane capital cost, adjusted by the membrane replacement factor.

    .. math::

        C_{op,tot} = f_{replace} * C_{mem} * A_{mem}

Code Documentation
------------------

* :mod:`watertap.costing.unit_models.membrane_distillation`

References
----------

Shamlou, E., Vidic, R., & Khanna, V. (2022). Optimization-based modeling and economic comparison of membrane distillation configurations for application in shale gas produced water treatment. Desalination, 526, 115513.