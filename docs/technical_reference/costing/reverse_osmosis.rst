Reverse Osmosis Costing Method
===============================

Costing Method Parameters
+++++++++++++++++++++++++

The following parameters are constructed for the unit on the FlowsheetCostingBlock (e.g., ``m.fs.costing.reverse_osmosis``) when applying the ``cost_reverse_osmosis`` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units"

   "**Standard RO**"
   "Membrane replacement factor (fraction of membrane replaced/year)", ":math:`f_{mem,\, replace}`", "``factor_membrane_replacement``", "0.2", ":math:`\text{yr}^{-1}`"
   "Membrane unit cost", ":math:`C_{mem}`", "``membrane_unit_cost``", "30", ":math:`\text{USD}_{2018}\text{/m}^2`"

   "**High-pressure RO**"
   "Membrane replacement factor (fraction of membrane replaced/year)", ":math:`f_{mem,\, replace}`", "``factor_membrane_replacement``", "0.2", ":math:`\text{yr}^{-1}`"
   "Membrane unit cost", ":math:`C_{mem}`", "``membrane_unit_cost``", "75", ":math:`\text{USD}_{2018}\text{/m}^2`"

Costing Method Variables
++++++++++++++++++++++++

The following variables are constructed on the unit block (e.g., ``m.fs.unit.costing``) when applying the ``cost_reverse_osmosis`` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Membrane area", ":math:`A_{mem}`", "``area``", "None", ":math:`\text{m}^2`"

Capital Cost Calculations
+++++++++++++++++++++++++

Capital cost is dependent upon the RO membrane area, :math:`A_{mem}`, as shown in the equation below.

    .. math::

        C_{cap,tot} = A_{mem} * C_{mem}

 
Operating Cost Calculations
+++++++++++++++++++++++++++

The membrane replacement cost :math:`C_{op,\, mem}` is calculated on an annual basis from the membrane cost and replacement frequency factor:

    .. math::

        C_{op,\, mem} = f_{mem,\, replace}C_{mem}A_{mem}

.. note:: 
    
    Energy consumption for RO units in WaterTAP is calculated via the separate :ref:`pump costing method <pump_costing>` for the :ref:`pump` unit model.
 
Code Documentation
------------------

* :mod:`watertap.costing.unit_models.reverse_osmosis`

References
----------

| T. V. Bartholomew, N. S. Siefert and M. S. Mauter (2018).
| Cost Optimization of Osmotically Assisted Reverse Osmosis.
| *Environ Sci Technol* 2018 Vol. 52 Issue 20 Pages 11813-11821.
| DOI: 10.1021/acs.est.8b02771

| Y.-Y. Lu, Y.-D. Hu, X.-L. Zhang, L.-Y. Wu and Q.-Z. Liu (2007).
| Optimum design of reverse osmosis system under different feed concentration and product specification.
| *Journal of Membrane Science* 2007 Vol. 287 Issue 2 Pages 219-229.
| 10.1016/j.memsci.2006.10.037