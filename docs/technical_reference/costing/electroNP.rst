ElectroN-P Costing Method
==========================

Costing Method Parameters
+++++++++++++++++++++++++

The following parameters are constructed for the unit on the FlowsheetCostingBlock (e.g., `m.fs.costing.electroNP`) when applying the `cost_electroNP` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units"

   "Hydraulic retention time", ":math:`HRT`", "HRT", "1.3333", ":math:`\text{hr}`"
   "Reactor sizing cost", ":math:`C_V`", "sizing_cost", "1000", ":math:`\text{USD}_{2020}\text{/m}^3`"
   "Magnesium chloride cost", ":math:`C_{MgCl2}`", "magnesium_chloride_cost", "0.0786", ":math:`\text{USD}_{2020}\text{/kg}`"
   "Phosphorus recovery value*", ":math:`C_{RP}`", "phosphorus_recovery_value", "-0.07", ":math:`\text{USD}_{2020}\text{/kg}`"

\* Negative value represents revenue from recovering phosphorus

Costing Method Variables
++++++++++++++++++++++++

The following variables are used on the unit block (e.g., m.fs.unit.costing) when applying the `cost_electroNP` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Inlet volumetric flow rate", ":math:`Q_{in}`", "mixed_state[0].flow_vol", "[t]", ":math:`\text{m}^3\text{/hr}`"
   "Byproduct phosphate mass flow rate", ":math:`Q_{byproduct, S_{PO4}}`", "byproduct_state[0].conc_mass_comp[0, S_PO4]", "[t]", ":math:`\text{kg/hr}`"
   "Magnesium chloride flowrate", ":math:`Q_{MgCl2}`", "magnesium_chloride_dosage", "[t]", ":math:`\text{kg/hr}`"

Capital Cost Calculations
+++++++++++++++++++++++++

Capital cost depends on the unit's inlet flow rate, :math:`Q_{in}`, as shown in the equation below:

    .. math::

        C_{cap,tot} = HRT * Q_{in} * C_V

 
Operating Cost Calculations
+++++++++++++++++++++++++++

Operating/maintenance costs consist of magnesium chloride usage cost and phosphorus recovery revenue:

    .. math::

        C_{op,tot} = C_{op,MgCl2}+C_{op,RP}

    .. math::

        C_{op,MgCl2} = Q_{MgCl2} * C_{MgCl2}
        C_{op,RP} = Q_{byproduct, S_{PO4}} * C_{RP}

 
Code Documentation
------------------

* :mod:`watertap.costing.unit_models.electroNP`