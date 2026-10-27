UV with Advanced Oxidation Processes Costing Method
====================================================

Costing Method Parameters
+++++++++++++++++++++++++

The following parameters are constructed for the unit on the FlowsheetCostingBlock (e.g., ``m.fs.costing.ultraviolet``) when applying the ``cost_uv_aop`` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units"

   "UV replacement factor", ":math:`f`", "``factor_lamp_replacement``", "0.33278", ":math:`\text{dimensionless}`"
   "UV reactor cost", ":math:`C_F`", "``reactor_cost``", "202.346", ":math:`\text{USD}_{2018}\text{/m}^3\text{/hr}`"
   "UV lamps, sleeves, ballasts and sensors cost", ":math:`C_l`", "``lamp_cost``", "235.5", ":math:`\text{USD}_{2018}\text{/kW}`"

Costing Method Variables
++++++++++++++++++++++++

The following variables are used on the unit block (e.g., ``m.fs.unit.costing``) when applying the ``cost_uv_aop`` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Inlet volumetric flowrate", ":math:`Q_{in}`", "``flow_vol``", "[t]", ":math:`\text{m}^3\text{/s}`"
   "Electricity demand", ":math:`E`", "``electricity_demand``", "[t]", ":math:`\text{W}`"

Capital Cost Calculations
+++++++++++++++++++++++++

Capital cost is the sum of the unit reactor cost and lamp cost:

    .. math::

        C_{cap,tot} = C_{cap,reactor}+C_{cap,lamp}

    .. math::

        C_{cap,reactor} = C_F \times Q_{in}
        C_{cap,lamp} = C_l \times E

 
Operating Cost Calculations
+++++++++++++++++++++++++++

Fixed operating cost depends on the electricity demand of lamps, :math:`E`, as shown in the equation below:

    .. math::

        C_{op,tot} = f \times C_l \times E

Code Documentation
------------------

* :mod:`watertap.costing.unit_models.uv_aop`

References
----------

| Costing relationships and parameters were developed from data in the UV System Cost Analysis Tool (UVCAT)
| https://www.nyserda.ny.gov/About/Publications/Energy-Analysis-Reports-and-Studies/Water-and-Wastewater-Technical-Reports/Optimization-of-UV-Disinfection
| Prepared by Carollo for NYSERDA; Appendix A

| Wright, H., Gaithuma, D., Greene, D., Aieta, M. (2006) 
| Integration of validation, design, and operation provides optimal implementation of UV disinfection.
| American Water Works Association, 98, 81-92
| https://doi.org/10.1002/j.1551-8833.2006.tb07780.x
