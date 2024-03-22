.. _watertap_costing_detailed:

Detailed WaterTAP Costing Package
=================================

.. currentmodule:: watertap.costing.watertap_costing

The detailed WaterTAP costing package contains the costing package with detailed technoeconomic data. It inherits all the functionality and parameters of the :ref:`WaterTAPCostingBlockData` base class :ref:`technical_reference/costing/costing_base:Common Global Costing Parameters`.

Costing Index and Technoeconomic Factors
----------------------------------------

The following technoeconomic factors are specific to the WaterTAP Costing Detailed Package. These parameters are chosen to give identical results as the WaterTAPCosting package which uses only simpler higher-level parameters.

======================================  ====================  =======================================  ===============  ===================================================================
            Cost factor                  Variable                 Name                                 Default Value    Description
======================================  ====================  =======================================  ===============  ===================================================================
Land cost factor                           :math:`f_{land}`    ``land_cost_percent_FCI``                0%              Unit process land cost as percentage of capital cost
Working capital cost factor                :math:`f_{work}`    ``working_capital_percent_FCI``          0%              Unit process working capital cost as percentage of capital cost
Salaries cost factor                       :math:`f_{sal}`     ``salaries_percent_FCI``                 0.20134%        Unit process salaries cost as percentage of capital cost
Maintenance cost factor                    :math:`f_{maint}`   ``maintenance_costs_percent_FCI``        1.6107%         Unit process maintenance costs as percentage of capital cost
Lab cost factor                            :math:`f_{lab}`     ``laboratory_fees_percent_FCI``          0.60403%        Unit process laboratory costs as percentage of capital cost
Insurance/taxes cost factor                :math:`f_{ins}`     ``insurance_and_taxes_percent_FCI``      0.40268%        Unit process insurance & taxes cost as percentage of capital cost
Benefits cost factor                       :math:`f_{ben}`     ``benefit_percent_of_salary``            90%             Unit process benefits cost as percentage of salary costs
Total investment factor                    :math:`f_{toti}`    ``total_investment_factor``              1.0             Total investment factor [calculated]
Maintenance-labor-chemical factor          :math:`f_{mlc}`     ``maintenance_labor_chemical_factor``    0.03            Maintenance, labor, and chemical factor [calculated]
======================================  ====================  =======================================  ===============  ===================================================================

High-level factor calculations
++++++++++++++++++++++++++++++

Total total investement factor :math:`f_{toti}` is calculated from the land cost and working capital cost factors :math:`f_{land}`, :math:`f_{work}`:

    .. math::

        f_{toti} = 1 + f_{work} + f_{land}

Total maintenance-labor-chemical factor :math:`f_{mlc}` is calculated from the salaires, maintenance, lab, insurance/taxes, and benefits cost factors:

    .. math::

        f_{mlc} = f_{sal} + f_{ben} f_{sal} + f_{maint} + f_{lab} + f_{ins}


Costing Process-Wide Costs 
--------------------------

The WaterTAPCostingDetailed class includes variables and constraints necessary to calculate process-wide costs, in addition to those provided by the WaterTAP Costing Framework :ref:`technical_reference/costing/costing_base:Costing Process-Wide Costs`.

=============================================  ====================  =====================================  ==============================================================================
                 Cost                               Variable                 Name                               Description
=============================================  ====================  =====================================  ==============================================================================
Unit capital cost                               :math:`C_{ca,u}`      ``aggregate_capital_cost``            Unit processes capital cost
Land cost                                       :math:`C_{land}`      ``land_cost``                         Cost of land for unit process
Working capital cost                            :math:`C_{work}`      ``working_capital``                   Working capital for unit process
Salary cost                                     :math:`C_{sal}`       ``salary_cost``                       Salary cost for unit process
Benefits cost                                   :math:`C_{ben}`       ``benefits_cost``                     Benefits cost for unit process
Maintenance cost                                :math:`C_{maint}`     ``maintenance_cost``                  Maintenance cost for unit process
Laboratory cost                                 :math:`C_{lab}`       ``laboratory_cost``                   Laboratory cost for unit process
Insurance & taxes cost                          :math:`C_{ins}`       ``insurance_and_taxes_cost``          Insurance & taxes for unit process
Total annualized cost                           :math:`C_{annual}`    ``total_annualized_costs``            Total cost on a annualized basis
=============================================  ====================  =====================================  ==============================================================================

Land cost is defined as:

    .. math::
   
        C_{land} = f_{land} C_{ca,u}

Working capital cost is defined as:

    .. math::
   
        C_{land} = f_{work} C_{ca,u}

Salary cost is defined as:

    .. math::
   
        C_{sal} = f_{sal} C_{ca,u}

Benefits cost is defined as:

    .. math::
   
        C_{ben} = f_{ben} f_{sal} C_{ca,u}

Maintenance cost is defined as:

    .. math::
   
        C_{maint} = f_{maint} C_{ca,u}

Land cost is defined as:

    .. math::
   
        C_{land} = f_{land} C_{ca,u}

All other costing calculations are provided through the :ref:`WaterTAPCostingBlockData`: :ref:`technical_reference/costing/costing_base:Costing Calculations`.

Class Documentation
-------------------

* :class:`WaterTAPCostingDetailedData`
