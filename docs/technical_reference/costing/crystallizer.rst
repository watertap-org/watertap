Crystallizer Costing Method
============================

Costing Method Parameters
+++++++++++++++++++++++++

The following parameters are constructed for the unit on the FlowsheetCostingBlock (e.g., `m.fs.costing.crystallizer`) when applying the `cost_crystallizer` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units"

   "**Mass-based Capital Costing**"
   "Reference free-on-board (FOB) capital cost :math:`^1`", ":math:`Cost_{ref}`", "fob_unit_cost", "675000", ":math:`\text{USD}_{2007}`"
   "Reference crystallizer capacity :math:`^1`", ":math:`size_{ref}`", "ref_capacity", "1", ":math:`\text{kg/s}`"
   "Crystallizer cost exponent parameter :math:`^1`", ":math:`n`", "ref_exponent", "0.53", ":math:`\text{dimensionless}`"
   "Installed equipment cost factor :math:`^2`", ":math:`IEC`", "iec_percent", "1.43", ":math:`\text{dimensionless}`"

   "**Volume-based Capital Costing**"
   "Capital cost A parameter :math:`^3`", "A", "volume_cost", "16320", ":math:`\text{USD}_{2007}\text{/ft}^3`"
   "Capital cost B parameter :math:`^3`", "B", "vol_basis_exponent", "0.47", ":math:`\text{dimensionless}`"

   "**Operating Costs**"
   "Heating steam Pressure :math:`^4`", ":math:`P_{steam}`", "steam_pressure", "3", ":math:`\text{bar}`"
   "Heating steam cost parameter :math:`^5`", ":math:`Cost_{steam}`", "steam_cost", "0.004", ":math:`\text{USD}_{2018}\text{/m}^3`"
   "Recirculation pump head height", ":math:`h_{rec}`", "pump_head_height", "1", ":math:`\text{m}`"
   "Recirculation pump efficiency", ":math:`\eta_{pump}`", "efficiency_pump", "0.7", ":math:`\text{dimensionless}`"

Costing Method Variables
++++++++++++++++++++++++

There are no costing method variables unique to the crystallizer.

Capital Cost Calculations
+++++++++++++++++++++++++

The crystallizer offers two options for computing the capital cost: mass-based costing or volume-based costing.

The mass-based capital cost is dependent upon the mass of solid crystals produced in the crystallizer, :math:`S`, as shown in the equation below.

    .. math::

        C_{cap,tot} = IEC * Cost_{ref} * (\frac{S}{size_{ref}})^{n}

The volume-based capital cost is dependent upon the unit's volume, :math:`V`, as shown in the equation below.

    .. math::

        C_{cap,tot} = A * V^{B}

 
Operating Cost Calculations
+++++++++++++++++++++++++++

The operating cost of the crystallizer is the sum of the electricity cost for the crystallizer recirculation pump, and the cost of steam for process heating. 

    .. math::

        C_{op,tot} = C_{op,electricity}+C_{op,heat}


:math:`C_{op,electricity}`  is computed with WaterTAP's standard approach for costing electricity consumption, with assumptions of :math:`h_{rec}=` 1m pump head height and :math:`\eta_{pump}` = 70% pump efficiency.


Process heat is supplied via steam at :math:`P_{steam}=` 3 bar (latent heat), and the process heating cost is computed from  the crystallizer heating requirement :math:`Q` (:math:`\text{kJ}`):


    .. math::

        C_{op,heat} = Cost_{steam} * \frac{Q}{\rho_{steam} * L_{v}}

where :math:`\rho_{steam}`  and :math:`L_v` are the density (:math:`\text{kg}\text{/m}^3`) and latent heat of condensation (:math:`\text{kJ/kg}`) of steam, respectively.

 
Code Documentation
------------------

* :mod:`watertap.costing.unit_models.crystallizer`

References
----------
[1] Woods, Donald R (2007).
Rules of Thumb in Engineering Practice.
Wiley. 2007. `DOI: 10.1002/9783527611119 <https://onlinelibrary.wiley.com/doi/book/10.1002/9783527611119>`_.


[2] Diab, Samir and Gerogiorgis, Dimitrios I (2017). 
Technoeconomic Evaluation of Multiple Mixed Suspension-Mixed Product Removal (MSMPR) Crystallizer Configurations for Continuous Cyclosporine Crystallization. 
*ACS Organic Process Research & Development*, Vol. 21, No. 10 p. 1571-1587. `DOI: 10.1021/acs.oprd.7b00225 <https://pubs.acs.org/doi/10.1021/acs.oprd.7b00225>`_.

[3] Yusuf, A et. al. (2019). 
CO2 utilization from power plant: A comparative techno-economic assessment of soda ash production and scrubbing by monoethanolamine.
*Journal of Cleaner Production*, Vol. 237, p. 117760. `DOI: 10.1016/j.jclepro.2019.117760 <https://doi.org/10.1016/j.jclepro.2019.117760>`_.

[4] Dutta, B. 
Principles of mass transfer and separation processes. PHI Learning, 2007.

[5] Panagopoulos, Argyris (2020) 
Process simulation and techno-economic assessment of a zero liquid discharge/multi-effect desalination/thermal vapor compression (ZLD/MED/TVC) system. 
*International Journal of Energy Research* , Vol. 44, No. 1, p. 473-495. `DOI: 10.1002/er.4948 <https://doi.org/10.1002/er.4948>`_.
