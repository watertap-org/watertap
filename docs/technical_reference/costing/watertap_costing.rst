.. _watertap_costing:

WaterTAP Costing Package
========================

.. currentmodule:: watertap.costing.watertap_costing

The WaterTAP costing module contains the costing package with simplified technoeconomic data. It inherits all the functionality and parameters of the :ref:`WaterTAPCostingBlockData` base class :ref:`technical_reference/costing/costing_base:Common Global Costing Parameters`.


Costing Index and Technoeconomic Factors
----------------------------------------

The following technoeconomic factors are specific to the WaterTAP Costing Package

=============================================  ====================  =======================================  ===============  ==============================================================================
                 Cost factor                     Variable                 Name                                 Default Value    Description
=============================================  ====================  =======================================  ===============  ==============================================================================
Total investment factor                           :math:`f_{toti}`    ``total_investment_factor``              1.0             Total investment factor (investment cost / equipment cost)
Maintenance-labor-chemical factor                 :math:`f_{mlc}`     ``maintenance_labor_chemical_factor``    0.03            Maintenance, labor, and chemical factor (fraction of equipment cost / year)
=============================================  ====================  =======================================  ===============  ==============================================================================

Costing Calculations
--------------------

All costing calculations are provided through the :ref:`WaterTAPCostingBlockData`: :ref:`technical_reference/costing/costing_base:Costing Calculations`.

Class Documentation
-------------------

* :class:`WaterTAPCostingData`

