Osmotically Assisted Reverse Osmosis Costing Method
====================================================

Costing Method Parameters
+++++++++++++++++++++++++

The following parameters are constructed when applying the `cost_osmotically_assisted_reverse_osmosis` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units"

   "description", ":math:`Symbol_{example}`", "parameter_name", "1", ":math:`\text{dimensionless}`"

Costing Method Variables
++++++++++++++++++++++++

The following variables are constructed when applying the `cost_osmotically_assisted_reverse_osmosis` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Default Value", "Units"

   "description", ":math:`Symbol_{example}`", "variable_name", "1", ":math:`\text{dimensionless}`"

Capital Cost Calculations
+++++++++++++++++++++++++

Describe capital costs..keep it concise where possible

    .. math::

        C_{cap,tot} = C_{cap,example1}+C_{cap,example2}+C_{cap,other}

    .. math::

        & C_{cap,example1} = fill in equation for each component in total capex equation

 
Operating Cost Calculations
+++++++++++++++++++++++++++

Describe operating/maintenance costs..keep it concise where possible

    .. math::

        C_{op,tot} = C_{op,example1}+C_{op,example2}+C_{op,other}

    .. math::

        & C_{op,example1} = fill in equation for each component in total opex equation

 
Code Documentation
------------------

* :mod:`watertap.unit_models.osmotically_assisted_reverse_osmosis`
* :mod:`watertap.costing.unit_models.osmotically_assisted_reverse_osmosis`

References
----------
Aim to include at least one reference in most cases, but delete this section if no references used for cost relationships/default values