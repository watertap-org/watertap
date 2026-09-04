.. _watertap_costing:

WaterTAP Costing Package
========================

.. currentmodule:: watertap.costing.watertap_costing

The WaterTAP costing package contains the costing package with simplified technoeconomic data. It inherits all the functionality and parameters of the :ref:`WaterTAPCostingBlockData` base class :ref:`technical_reference/costing/costing_base:Common Global Costing Parameters`.

.. _watertap_costing_tea_factors:

Costing Index and Technoeconomic Factors
----------------------------------------

The following technoeconomic factors are specific to the WaterTAP Costing Package:

=============================================  ====================  =======================================  ===============  ==============================================================================
                 Cost factor                     Variable                 Name                                 Default Value    Description
=============================================  ====================  =======================================  ===============  ==============================================================================
Total investment factor                           :math:`f_{toti}`    ``total_investment_factor``              1.0             Total investment factor (investment cost / equipment cost)
Maintenance-labor-chemical factor                 :math:`f_{mlc}`     ``maintenance_labor_chemical_factor``    0.03            Maintenance, labor, and chemical factor (fraction of equipment cost / year)
=============================================  ====================  =======================================  ===============  ==============================================================================

For the WaterTAP costing package, the base currency year and base period should be set using the ``base_currency_year`` and ``base_period`` configuration arguments when instantiating the costing package.


.. code-block:: python

    m.fs.costing = WaterTAPCosting(base_currency_year=2021, base_period="year")


The default costing year is 2018, but any year between 1990 and 2023 can be used. Any unit of time can be used for the base period.


.. important:: 
    Though users **could** directly set the ``base_currency`` on the flowsheet costing block (e.g., ``m.fs.costing.base_currency = pyunits.USD_2023``), this is discouraged. 
    It is recommended to use the ``base_currency_year`` configuration argument when instantiating the WaterTAP costing package to ensure consistency 
    across all costing calculations and parameters. 

Costing Calculations
--------------------

All costing calculations are provided through the :ref:`WaterTAPCostingBlockData`: :ref:`technical_reference/costing/costing_base:Costing Calculations`.

Class Documentation
-------------------

* :class:`WaterTAPCostingData`

