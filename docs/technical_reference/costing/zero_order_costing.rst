.. _zero_order_costing:

Zero Order Costing Package
==========================

.. currentmodule:: watertap.costing.zero_order_costing

The zero order costing module contains the costing package typically used for zero order models, though it can also be used for high-fidelity models. Technoeconomic data used for zero order models is contained in the
``.yaml`` file for that model located in the data/techno_economic folder.


Usage
-----

The ZeroOrderCosting class contains all the variables and constraints needed to cost a unit model derived from the :ref:`ZeroOrderBaseData`. It also inherits the functionality of the :ref:`WaterTAPCostingBlockData`.

The code below shows an outline of how the ZeroOrderCostingData class is intended to be used to cost zero-order type models.

.. doctest::

  from pyomo.environ import ConcreteModel

  from idaes.core import FlowsheetBlock

  from watertap.costing.zero_order_costing import ZeroOrderCosting
  from watertap.core.wt_database import Database
  from watertap.core.zero_order_properties import WaterParameterBlock
  from watertap.unit_models.zero_order import MyZOUnit


  m = ConcreteModel()
  m.db = Database()
  m.fs = FlowsheetBlock(dynamic=False)
  m.fs.params = WaterParameterBlock(solute_list=["comp_a", "comp_b", "comp_c"])
  m.fs.costing = ZeroOrderCosting()
  m.fs.unit = MyZOUnit(property_package=m.fs.params, database=m.db)

  # Add necessary statements to fix component flows prior to solve

Costing Zero Order Models
-------------------------

The ZeroOrderCostingData class includes variables and constraints necessary to calculate process-wide costs:

=============================================  ====================  =====================================  ==============================================================================
                 Cost                               Variable                 Name                               Description
=============================================  ====================  =====================================  ==============================================================================
Total capital cost                              :math:`C_{ZO,tot}`    ``total_capital_cost``                Total capital cost
Unit capital cost                               :math:`C_{ZO,u}`      ``aggregate capital_cost``            Unit processes capital cost
Total operating cost                            :math:`C_{op,tot}`    ``total_operating_cost``              Total operating cost for unit process
Total fixed operating cost                      :math:`C_{op,fix}`    ``total_fixed_operating_cost``        Total fixed operating cost for unit process
Total variable operating cost                   :math:`C_{op,var}`    ``total_variable_operating_cost``     Total variable operating cost for unit process
Land cost                                       :math:`C_{land}`      ``land_cost``                         Cost of land for unit process
Working capital cost                            :math:`C_{work}`      ``working_capital``                   Working capital for unit process
Salary cost                                     :math:`C_{sal}`       ``salary_cost``                       Salary cost for unit process
Benefits cost                                   :math:`C_{ben}`       ``benefits_cost``                     Benefits cost for unit process
Maintenance cost                                :math:`C_{maint}`     ``maintenance_cost``                  Maintenance cost for unit process
Laboratory cost                                 :math:`C_{lab}`       ``laboratory_cost``                   Laboratory cost for unit process
Insurance & taxes cost                          :math:`C_{ins}`       ``insurance_and_taxes_cost``          Insurance & taxes for unit process
Total annualized cost                           :math:`C_{annual}`    ``total_annualized_costs``            Total cost on a annualized basis
=============================================  ====================  =====================================  ==============================================================================

Calculations for each of these costs are presented below.

Costing Index and Technoeconomic Factors
----------------------------------------

Costing indices are available in ``default_case_study.yaml`` located in the data/techno_economic folder.

WaterTAP uses the CE (Chemical Engineering) Cost Index to help account for the time-value of investments and are used in the capital
and operating cost calculations. Unit process capital costs are adjusted to the year of the case study. The default year is 2018.

Other technoeconomic factors used to calculate various system metrics, capital, and operating costs are presented in the table below:

=============================================  ====================  =====================================  ===============  ==============================================================================
                 Cost factor                     Variable                 Name                               Default Value    Description
=============================================  ====================  =====================================  ===============  ==============================================================================
Plant capacity utilization factor                 :math:`f_{util}`    ``utilization_factor``                 100%               Percentage of year plant is operating
Plant lifetime                                    :math:`L`           ``plant_lifetime``                     30 yr              Expected lifetime of unit process
Electricity price                                 :math:`P`           ``electricity_cost``                   $0.0595/kWh        Electricity price in 2019 USD.
Land cost factor                                  :math:`f_{land}`    ``land_cost_percent_FCI``              0.15%              Unit process land cost as percentage of capital cost
Working capital cost factor                       :math:`f_{work}`    ``working_capital_percent_FCI``        5%                 Unit process working capital cost as percentage of capital cost
Salaries cost factor                              :math:`f_{sal}`     ``salaries_percent_FCI``               0.1%               Unit process salaries cost as percentage of capital cost
Maintenance cost factor                           :math:`f_{maint}`   ``maintenance_costs_percent_FCI``      0.8%               Unit process maintenance costs as percentage of capital cost
Lab cost factor                                   :math:`f_{lab}`     ``laboratory_fees_percent_FCI``        0.3%               Unit process laboratory costs as percentage of capital cost
Insurance/taxes cost factor                       :math:`f_{ins}`     ``insurance_and_taxes_percent_FCI``    0.2%               Unit process insurance & taxes cost as percentage of capital cost
Benefits cost factor                              :math:`f_{ben}`     ``benefit_percent_of_salary``          90%                Unit process benefits cost as percentage of salary costs
Weighted average cost of capital factor           :math:`f_{wacc}`    ``wacc``                               5%                 Weighted average cost of capital over plant lifetime
Capital recovery factor                           :math:`f_{crf}`     ``capital_recovery_factor``            6.51%              Calculated from default :math:`f_{WACC}` and :math:`L` values
=============================================  ====================  =====================================  ===============  ==============================================================================


The capital recovery factor is calculated with:

    .. math::

        f_{crf} = \frac{ f_{wacc} (1 + f_{wacc}) ^ L}{ (1 + f_{wacc}) ^ L - 1}


Capital Cost Calculations
+++++++++++++++++++++++++

In general, unit process capital costs :math:`C_{ZO,u}` for zero order unit models are a function of flow:

    .. math::

        C_{ZO,u} = A \bigg( \frac{Q_{in}}{Q_{basis}} \bigg) ^ {B}

:math:`Q_{basis}`, :math:`A`, and :math:`B` are specific to the unit model and can be found in the unit model ``.yaml`` file.
The :math:`A` value has units of USD for the costing reference year of the unit. For example, if the unit costing model
is from a reference that used 2015 USD, the units for :math:`A` are ``USD_2015``. After calculating the costs in 2015 USD, WaterTAP
adjusts the cost to the user-specified year via the Consumer Price Index.

The total capital cost of a zero order model :math:`C_{ZO,tot}` includes the land cost :math:`C_{land}` and working
capital costs :math:`C_{work}`:

    .. math::

        C_{ZO,tot} = C_{ZO,u} + C_{land} + C_{work}

Where:

    .. math::

        & C_{land} = f_{land} C_{ZO,u} \\\\
        & C_{work} = f_{work} C_{ZO,u}


Custom Capital Cost Methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are several zero order models that have costing relationships that do not follow this general form. If that is the case, a custom costing method can
be added to the unit model class to perform that calculation.

Zero order models that have custom capital costing methods include:

* Brine concentrator - ``cost_brine_concentrator()``
* CANDOP - ``cost_CANDOP()``
* Chemical addition - ``cost_chemical_addition()``
* Chlorination - ``cost_chlorination()``
* Coagulation/Flocculation - ``cost_coag_and_floc()``
* Deep well injection - ``cost_deep_well_injection()``
* DMBR - ``cost_dmbr()``
* Electrochemical nutrient removal - ``cost_electrochemical_nutrient_removal()``
* Evaporation pond - ``cost_evaporation_pond()``
* Filter press - ``cost_filter_press()``
* Fixed bed - ``cost_fixed_bed()``
* GAC - ``cost_gac()``
* Landfill - ``cost_landfill()``
* MABR - ``cost_mabr()``
* Ion exchange - ``cost_ion_exchange()``
* Iron/Manganese removal - ``cost_iron_and_manganese_removal()``
* Metab - ``cost_metab()``
* Nanofiltration - ``cost_nanofiltration()``
* Ozone - ``cost_ozonation()``
* Ozone + AOP - ``cost_ozonation_aop()``
* Photothermal membrane - ``cost_photothermal_membrane()``
* Sedimentation - ``cost_sedimentation()``
* Storage tank - ``cost_storage_tank()``
* Surface discharge - ``cost_surface_discharge()``
* UV irradiation - ``cost_uv()``
* UV + AOP - ``cost_uv_aop()``
* Well field - ``cost_well_field()``

To add a custom capital calculation method, the unit model class must register its custom costing method by setting its `default_costing_method` attribute.


Operating Cost Calculations
+++++++++++++++++++++++++++

Total operating costs for zero order models :math:`C_{op,tot}` include fixed :math:`C_{op,fix}` and variable operating costs :math:`C_{op,var}`:

    .. math::

        C_{op,tot} = C_{op,fix} + C_{op,var}

The total fixed operating costs are calculated as:

    .. math::

        C_{op,fix} = C_{sal} + C_{ben} + C_{maint} + C_{lab} + C_{ins}

Where:

    .. math::

        & C_{sal} = f_{sal} C_{ZO,u} \\\\
        & C_{ben} = f_{ben} C_{sal} \\\\
        & C_{maint} = f_{maint} C_{ZO,u} \\\\
        & C_{lab} = f_{lab} C_{ZO,u} \\\\
        & C_{ins} = f_{ins} C_{ZO,u}

Variable operating costs include any chemical additions, electricity costs, and other variable costs such as equipment
replacements.

    .. math::

        C_{op,var} = C_{chem} + C_{elec} + C_{other} + C_{op,tot}

Chemical costs are based on the chemical dosage for a given chemical addition. Default chemical costs are found in ``default_case_study.yaml``. The
annual chemical costs are calculated as:

    .. math::

        C_{chem} = \sum_{k}^{n} D_k C_k Q_{in} f_{util}`

Where :math:`D` is the dose of chemical :math:`k` and :math:`C` is the unit cost of chemical :math:`k`.

Electricity costs :math:`C_{elec}` are based on the energy intensity :math:`E` of the unit process
(see individual unit model documentation for details).
The annual electricity costs are calculated as:

    .. math::

        C_{elec} = E Q_{in} f_{util} P

LCOW Calculation
++++++++++++++++

The Levelized Cost Of Water (LCOW) [$/m3] is a metric used to assess the technoeconomics of a unit process:

    .. math::

        LCOW = \frac{ f_{crf} C_{ZO,tot} + C_{op,tot} }{Q f_{util} }

Other aggregates, like specific energy consumption, are provided through the :ref:`WaterTAPCostingBlockData`: :ref:`technical_reference/costing/costing_base:Aggregates Metrics`. 

Class Documentation
-------------------

* :class:`ZeroOrderCostingData`
