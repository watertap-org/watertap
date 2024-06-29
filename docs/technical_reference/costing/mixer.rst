Mixer Costing Method
=====================

Costing Method Parameters
+++++++++++++++++++++++++

The following parameters are constructed for the unit on the FlowsheetCostingBlock (e.g., `m.fs.costing.mixer`) when applying the `cost_mixer` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units"

   "**Standard mixer**"
   "Mixer unit cost", ":math:`C_{mix, standard}`", "``unit_cost``", "361", ":math:`\text{USD}_{2018}\text{/L/s}`"

   "**NaOCl mixer**"
   "Mixer unit cost", ":math:`C_{mix, NaOCl}`", "``cost``", "5.08", ":math:`\text{USD}_{2018}\text{/m}^{3}\text{/day}`"
   "NaOCl cost", ":math:`C_{NaOCl}`", "``unit_cost``", "0.23", ":math:`\text{USD}_{2018}\text{/kg}`"
   "NaOCl purity", ":math:`p_{NaOCl}`", "``purity``", "0.15", ":math:`\text{dimensionless}`"

   "**CaOH2 mixer**"
   "Mixer unit cost", ":math:`C_{mix, CaOH2}`", "``cost``", "873.911", ":math:`\text{USD}_{2018}\text{/kg/day}`"
   "CaOH2 cost", ":math:`C_{CaOH2}`", "``unit_cost``", "0.12", ":math:`\text{USD}_{2018}\text{/kg}`"
   "CaOH2 purity", ":math:`p_{CaOH2}`", "``purity``", "1", ":math:`\text{dimensionless}`"

Costing Method Variables
++++++++++++++++++++++++

The following variables are constructed on the unit block (e.g., m.fs.unit.costing) when applying the `cost_mixer` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Inlet volumetric flow rate", ":math:`Q_{in}`", "``flow_in``", "[t]", ":math:`\text{m}^3\text{/hr}`"
   "Lime mass flow", ":math:`M_{in}`", "``lime_kg_per_day``", "[t]", ":math:`\text{kg/day}`"
   "Dosing rate", ":math:`D`", "``dosing_rate``", "[j]", ":math:`\text{kg/s}`"

Capital Cost Calculations
+++++++++++++++++++++++++

Capital cost is dependent upon the volumetric flow rate, :math:`Q_{in}`, or
lime mass flow, :math:`M_{in}`, as shown in the equations below.

   "**Standard mixer**"

    .. math::

        C_{cap,tot} = C_{mix, standard} * Q_{in}


   "**NaOCl mixer**"

    .. math::

        C_{cap,tot} = C_{mix, NaOCl} * Q_{in}


   "**CaOH2 mixer**"

    .. math::

        C_{cap,tot} = C_{mix, CaOH2} * M_{in}

 
Operating Cost Calculations
+++++++++++++++++++++++++++

The total operating cost is dependent upon the dosing rate, :math:`D`, as shown in the equations below.

   "**Standard mixer**"

    .. math::

        C_{op,tot} = 0

   "**NaOCl mixer**"

    .. math::

        C_{op,tot} = D_{NaOCl} * C_{NaOCl} / p_{NaOCl}

   "**CaOH2 mixer**"

    .. math::

        C_{op,tot} = D_{CaOH2} * C_{CaOH2} / p_{CaOH2}

 
Code Documentation
------------------

* :mod:`watertap.costing.unit_models.mixer`
