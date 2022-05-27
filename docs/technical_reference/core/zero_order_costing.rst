.. _zero_order_costing:

Zero Order Costing Package
==========================

.. index::
   pair: watertap.core.zero_order_costing;ZeroOrderCostingData

.. currentmodule:: watertap.core.zero_order_costing

The zero order costing module contains the costing package used for zero order models.

Usage
-----

The ZeroOrderCosting class contains all the variables and constraints needed to cost a unit model derived from the :ref:`ZeroOrderBaseData`.

The code below shows an outline of how the ZeroOrderCostingData class is intended to be used to cost zero-order type models.

.. testcode::

  from pyomo.environ import ConcreteModel
  
  from idaes.core import FlowsheetBlock
  from watertap.core.zero_order_costing import ZeroOrderCosting
  from watertap.core.wt_database import Database
  from watertap.core.zero_order_properties import WaterParameterBlock
  from watertap.unit_models.zero_order import MyZOUnit


  m = ConcreteModel()
  m.db = Database()
  m.fs = FlowsheetBlock(default{"dynamic": False})
  m.fs.params = WaterParameterBlock(default={"solute_list": ["comp_a", "comp_b", "comp_c"]})
  m.fs.costing = ZeroOrderCosting()
  m.fs.unit = MyZOUnit(default={"property_package": m.fs.params, "database": m.db)

  # Add necessary statements to fix component flows prior to solve

Class Documentation
-------------------

The ZeroOrderCostingData class includes variables and constraints necessary to cost the following:

* Land
* Working capital
* Total capital
* Salary
* Benefits
* Maintenance
* Laboratory
* Insurance & Taxes
* Total fixed operating cost
* Total variable operating cost
* Total operating cost

Cost calculations for many of these are calculated as a percentage of the fixed capital investment (FCI). Details are provided in the following sections

Cost Calculation Details
------------------------

Broadly speaking, costs are classified as either capital or operating costs.

Costing Indices and Factors
-----------------------------------

Costing indices are available in ``default_case_study.yaml`` located in the data/techno_economoic folder. 

WaterTAP uses the Consumer Price Index (CPI) to help account for the time-value of investments and are used in the capital
and operating cost calculations. 

There are also various assumed costing factors for each case study read in from ``case_study_basis.csv``:

* Electricity price (:math:`P`)
* Plant capacity utilization (:math:`f_{util}`)
* Land cost as percent of FCI (:math:`f_{land}`)
* Working capital as percent of FCI (:math:`f_{work}`)
* Salaries as percent of FCI (:math:`f_{sal}`)
* Maintenance costs as percent of FCI (:math:`f_{maint}`)
* Laboratory costs as percent of FCI (:math:`f_{lab}`)
* Insurance/taxes as percent of FCI (:math:`f_{ins}`)
* Benefits as percent of salary (:math:`f_{ben}`)
* Assumed plant lifetime (:math:`L`)
* Weighted Average Cost of Capital (debt interest rate) (:math:`WACC`)

Capital Cost Calculations
+++++++++++++++++++++++++++

In general, capital costs :math:`C_{ZO}` for zero order unit models are a function of flow:

    .. math::

        C_{ZO} = A \bigg( \frac{Q_{in}}{Q_{basis}} \bigg) ^ {B}
|

:math:`Q_{basis}`, :math:`A`, and :math:`B` are specific to the unit model and can be found in the unit model ``.yaml`` file.


Custom Capital Cost Methods
++++++++++++++++++++++++++++++

There are several zero order models that have costing relationships that don't follow this general form. If that is the case, a custom costing method can 
be added to this base class to perform that calculation. 

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
* Sedimentation - `1cost_sedimentation()``
* Storage tank - ``cost_storage_tank()``
* Surface discharge - ``cost_surface_discharge()``
* UV irradiation - ``cost_uv()``
* UV + AOP - ``cost_uv_aop()``
* Well field - ``cost_well_field()``

To add a custom capital calculation method, the unit module class must be part of the ``import`` statement at the top of this file
and there must be an entry in the ``unit_mapping`` dictionary that maps the unit model class to the costing method. Convention is to name the method
``cost_unit_class_name`` without the "ZO" at the end of the class. For example, if the unit model class is ``MyUnitZO``, the custom cost method would be 
``cost_my_unit`` and ``MyUnitZO: cost_my_unit`` would be an entry in the ``unit_mapping`` dictionary found in this file.




.. autoclass:: ZeroOrderCostingData
    :members:
