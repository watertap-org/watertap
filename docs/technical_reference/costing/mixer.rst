Mixer Costing Method
=====================

Costing Method Parameters
+++++++++++++++++++++++++

The following parameters are constructed for the unit on the FlowsheetCostingBlock (e.g., `m.fs.costing.mixer`) when applying the `cost_mixer` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units"

   "**Standard mixer**"
   "Mixer unit cost", ":math:`C_{standard}`", "``unit_cost``", "361", ":math:`\text{USD}_{2018}\text{/L/s}`"

   "**NaOCl mixer**"
   "Mixer unit cost", ":math:`C_{NaOCl}`", "``unit_cost``", "160", ":math:`\text{USD}_{2018}\text{/L/s}`"

   "**CaOH2 mixer**"
   "Mixer unit cost", ":math:`C_{CaOH2}`", "``unit_cost``", "160", ":math:`\text{USD}_{2018}\text{/L/s}`"

Costing Method Variables
++++++++++++++++++++++++

The following variables are constructed on the unit block (e.g., m.fs.unit.costing) when applying the `cost_mixer` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "description", ":math:`Symbol_{example}`", "``variable_name``", "[t]", ":math:`\text{dimensionless}`"

Capital Cost Calculations
+++++++++++++++++++++++++

Describe capital costs..keep it concise where possible

    .. math::

        C_{cap,tot} = C_{cap,example1}+C_{cap,example2}+C_{cap,other}

    .. math::

        C_{cap,example1} = fill in equation for each component in total capex equation

 
Operating Cost Calculations
+++++++++++++++++++++++++++

Describe operating/maintenance costs..keep it concise where possible

    .. math::

        C_{op,tot} = C_{op,example1}+C_{op,example2}+C_{op,other}

    .. math::

        C_{op,example1} = fill in equation for each component in total opex equation

 
Code Documentation
------------------

* :mod:`watertap.costing.unit_models.mixer`

References
----------
Aim to include at least one reference in most cases, but delete this section if no references used for cost relationships/default values