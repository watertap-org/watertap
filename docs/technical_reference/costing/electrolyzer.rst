Electrolyzer Costing Method
============================

Costing Method Parameters
+++++++++++++++++++++++++

The following parameters are constructed for the unit on the FlowsheetCostingBlock (e.g., `m.fs.costing.electrolyzer`) when applying the `cost_electrolyzer` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units"

   "Membrane replacement factor (fraction of membrane replaced/year)", ":math:`f_{mem,\, replace}`", "factor_membrane_replacement", "0.33", ":math:`\text{yr}^{-1}`"
   "Membrane unit cost", ":math:`C_{mem}`", "membrane_unit_cost", "25", ":math:`\text{USD}_{2012}\text{/m}^2`"
   "Anode unit cost", ":math:`C_{anode}`", "anode_unit_cost", "300", ":math:`\text{USD}_{2005}\text{/m}^2`"
   "Cathode unit cost", ":math:`C_{cathode}`", "cathode_unit_cost", "600", ":math:`\text{USD}_{2005}\text{/m}^2`"
   "Membrane, anode, and cathode fraction of total capital", ":math:`f_{material}`", "fraction_material_cost", "0.65", ":math:`\text{dimensionless}`"


Costing Method Variables
++++++++++++++++++++++++

The following variables are constructed on the unit block (e.g., m.fs.unit.costing) when applying the `cost_electrolyzer` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Total unit capital cost", ":math:`C_{cap}`", "capital_cost", "None", ":math:`\text{USD}_{2020}`"
   "Fixed operating costs", ":math:`C_{op}`", "fixed_operating_cost", "None", ":math:`\text{USD}_{2020}\text{/yr}`"
   "Membrane capital cost", ":math:`C_{cap,\, mem}`", "membrane_cost", "None", ":math:`\text{USD}_{2020}`"
   "Anode capital cost", ":math:`C_{cap,\, anode}`", "anode_cost", "None", ":math:`\text{USD}_{2020}`"
   "Cathode capital cost", ":math:`C_{cap,\, cathode}`", "cathode_cost", "None", ":math:`\text{USD}_{2020}`"
   "Membrane replacement cost", ":math:`C_{op,\, mem}`", "membrane_replacement_cost", "None", ":math:`\text{USD}_{2020}\text{/yr}`"

Capital Cost Calculations
+++++++++++++++++++++++++

Capital costs contribute to the majority of material costs for the anode, cathode, and membrane. Each material cost is calculated individually then summed and increased by the assumed fractional cost to estimate the total capital cost (O’Brien, 2005).

    .. math::

        & C_{cap,\, mem} = C_{mem}A_{mem} \\\\
        & C_{cap,\, anode} = C_{anode}A_{anode} \\\\
        & C_{cap,\, cathode} = C_{cathode}A_{cathode} \\\\
        & C_{cap} = \frac{C_{cap,\, mem}+C_{cap,\, anode}+C_{cap,\, cathode}}{f_{material}}

Operating Cost Calculations
+++++++++++++++++++++++++++

The electrolyzer costing model considers electricity as a variable operating costs and membrane replacement as a fixed operating costs. Electricity is costed using ``cost_flow`` function applied to the ``power`` variable on the unit model. Replacement costs for the anode and cathode are not currently considered in the costing function.

    .. math::

        C_{op} = C_{op,\, mem} = f_{mem,\, replace}C_{mem}A_{mem}
 
Code Documentation
------------------

* :mod:`watertap.costing.unit_models.electrolyzer`

References
----------
Bommaraju, T. V., & O’Brien, T. F. (2015). Brine electrolysis. Electrochemistry Encyclopedia. https://knowledge.electrochem.org/encycl/art-b01-brine.htm

Kent, J. A. (Ed.). (2007). Kent and Riegel’s Handbook of Industrial Chemistry and Biotechnology. Springer US. https://doi.org/10.1007/978-0-387-27843-8

Kumar, A., Du, F., & Lienhard, J. H. (2021). Caustic Soda Production, Energy Efficiency, and Electrolyzers. ACS Energy Letters, 6(10), 3563–3566. https://doi.org/10.1021/acsenergylett.1c01827

O’Brien, T., Bommaraju, T. V., & Hine, F. (2005). Handbook of chlor-alkali technology. Springer.

Phillips, R., Edwards, A., Rome, B., Jones, D. R., & Dunnill, C. W. (2017). Minimising the ohmic resistance of an alkaline electrolysis cell through effective cell design. International Journal of Hydrogen Energy, 42(38), 23986–23994. https://doi.org/10.1016/j.ijhydene.2017.07.184

