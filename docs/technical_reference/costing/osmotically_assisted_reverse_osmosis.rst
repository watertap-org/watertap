Osmotically Assisted Reverse Osmosis Costing Method
====================================================

Costing Method Parameters
+++++++++++++++++++++++++

The following parameters are constructed for the unit on the FlowsheetCostingBlock (e.g., `m.fs.costing.osmotically_assisted_reverse_osmosis`) when applying the `cost_osmotically_assisted_reverse_osmosis` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units"

   "Membrane replacement factor", ":math:`f`", "factor_membrane_replacement", "0.15", ":math:`\text{dimensionless}`"
   "Membrane cost", ":math:`C_A`", "membrane_cost", "30", ":math:`\text{$/m^2}`"
   "High pressure membrane cost", ":math:`C_hA`", "high_pressure_membrane_cost", "50", ":math:`\text{$/m^2}`"

Costing Method Variables
++++++++++++++++++++++++

The following variables are used on the unit block (e.g., m.fs.unit.costing) when applying the `cost_osmotically_assisted_reverse_osmosis` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Units"

   "Membrane area", ":math:`A`", "area", ":math:`\text{m^2}`"

Capital Cost Calculations
+++++++++++++++++++++++++

There are two classes of OARO costing types (``standard`` and ``high_pressure``), the membrane cost is different with these two
methods and the default costing method is ``standard``.

For ``standard`` type, the  capital cost is defineds as:

    .. math::

        C_{cap,tot} = C_A * A

For ``high_pressure`` type, the  capital cost is defineds as:

    .. math::

        C_{cap,tot} = C_hA * A

 
Operating Cost Calculations
+++++++++++++++++++++++++++

The fixed operating cost is correlated to OARO costing types (``standard`` and ``high_pressure``) and the default is ``standard``.

For ``standard`` type, the  capital cost is defineds as:

    .. math::

        C_{op,tot} = f * C_A * A

For ``high_pressure`` type, the  capital cost is defineds as:

    .. math::

        C_{op,example1} = f * C_A * A

 
Code Documentation
------------------

* :mod:`watertap.costing.unit_models.osmotically_assisted_reverse_osmosis`

References
----------
Bartholomew, T. V., Mey, L., Arena, J. T., Siefert, N. S., & Mauter, M. S. (2017).
Osmotically assisted reverse osmosis for high salinity brine treatment. Desalination, 421, 3-11.