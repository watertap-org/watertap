.. _WaterTAPCostingBlockData:

WaterTAP Costing Base Class
===========================

.. index::
   pair: watertap.costing.costing_base;WaterTAPCostingBlockData

.. currentmodule:: watertap.costing.costing_base

The WaterTAP Costing Base class contains extensions, methods, and variables and constraints common to both WaterTAP Costing Packages, and which would be useful for creating custom costing packages for WaterTAP.


Extensions Over IDAES FlowsheetCostingBlock
-------------------------------------------

The WaterTAPCostingBlock class extends the functionality of the `IDAES Process Costing Framework <https://idaes-pse.readthedocs.io/en/stable/reference_guides/core/costing/costing_framework.html>`_ in two ways:

1. Unit models can self-register a default costing method by specifying a `default_costing_method` attribute. This allows the costing method(s) to be specified with the unit model definition.

.. testcode::

    import pyomo.environ as pyo
    import idaes.core as idc
    from watertap.costing import WaterTAPCosting

    def cost_unit_model(blk):
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

    @idc.declare_process_block_class("MyUnitModel")
    class MyUnitModelData(idc.UnitModelBlockData):

        @property
        def default_costing_method(self):
            # could point to a static method on
            # this class, could be function in 
            # a different module even
            return cost_unit_model

    m = pyo.ConcreteModel()
    m.fs = idc.FlowsheetBlock(dynamic=False)
    m.fs.costing = WaterTAPCosting()

    m.fs.my_unit = MyUnitModel()
    
    # the `default_costing_method_attribute` on the
    # unit model is checked, and the function
    # `cost_unit_model` returned then build the costing block
    m.fs.my_unit.costing = idc.UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
    )


2. The method `register_flow_type` will create a new Expression if a costing component is not already defined *and* the costing component is not constant. The default behavior in IDAES is to always create a new Var. This allows the user to specify intermediate values in `register_flow_type`. 

.. testcode::

    import pyomo.environ as pyo
    from idaes.core import FlowsheetBlock
    from watertap.costing import WaterTAPCosting

    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = WaterTAPCosting()

    m.fs.naocl_bulk_cost = pyo.Param(
        mutable=True,
        initialize=0.23,
        doc="NaOCl cost",
        units=pyo.units.USD_2018 / pyo.units.kg,
    )
    m.fs.naocl_purity = pyo.Param(
       mutable=True,
       initialize=0.15,
       doc="NaOCl purity",
       units=pyo.units.dimensionless,
    )

    # This will create an Expression m.fs.costing.naocl_cost whose expr is the second argument
    # so changes to m.fs.naocl_bulk_cost and m.fs.naocl_purity will affect the underlying
    # new Expression m.fs.costing.naocl_cost.
    m.fs.costing.register_flow_type("naocl", m.fs.naocl_bulk_cost / m.fs.naocl_purity)

    # This, however, will create a Var called m.fs.costing.caoh2_cost whose *value* is the second argument
    m.fs.costing.register_flow_type("caoh2", 0.12 * pyo.units.USD_2018 / pyo.units.kg)


Costing Process-Wide Costs
--------------------------

The WaterTAPCostingBlockData class includes variables necessary to calculate process-wide costs:

=============================================  ====================  =====================================  ==============================================================================
                 Cost                               Variable                 Name                               Description
=============================================  ====================  =====================================  ==============================================================================
Total capital cost                              :math:`C_{ca,tot}`    ``total_capital_cost``                Total capital cost -- defined by derived class
Total operating cost                            :math:`C_{op,tot}`    ``total_operating_cost``              Total operating cost for unit process -- defined by derived class
Aggregate electricity cost                      :math:`C_{el,tot}`    ``aggregate_electricity_cost``        Sum of all electricity costs
=============================================  ====================  =====================================  ==============================================================================


Common Global Costing Parameters
--------------------------------

The `_build_common_global_parameters` method builds common cost factor parameters necessary to calculate aggregates such as levelized cost of water (LCOW).
Note that the default values can be overwritten in the derived class.

=============================================  ====================  =====================================  ===============  ==============================================================================
                 Cost factor                     Variable                 Name                               Default Value    Description
=============================================  ====================  =====================================  ===============  ==============================================================================
Plant capacity utilization factor                 :math:`f_{util}`    ``utilization_factor``                 90%                Percentage of year plant is operating
Electricity price                                 :math:`P`           ``electricity_cost``                   $0.07/kWh          Electricity price in 2018 USD
Electricity carbon intensity                      :math:`f_{eci}`     ``electrical_carbon_intensity``        0.475 kg/kWh       Carbon intensity of electricity
Capital recovery factor                           :math:`f_{crf}`     ``capital_recovery_factor``            None               Calculated by derived class 
=============================================  ====================  =====================================  ===============  ==============================================================================


Aggregates Metrics
------------------

Levelized Cost of Water (LCOW)
++++++++++++++++++++++++++++++

For a given volumetric flow :math:`Q`, the LCOW, :math:`LCOW_{Q}` is calculated by the `add_LCOW` method as

    .. math::
    
        LCOW_{Q} = \frac{f_{crf}   C_{ca,tot} + C_{op,tot}}{f_{util} Q}

Specific Energy Consumption
+++++++++++++++++++++++++++

For a given volumetric flow `Q`, the specific energy consumption, :math:`E_{spec,Q}` is calculated by the `add_specific_energy_consumption` method as

    .. math::
    
        E_{spec,Q} = \frac{C_{el,tot}}{Q}

Specific Electrical Carbon Intensity
++++++++++++++++++++++++++++++++++++

For a given volumetric flow `Q`, the specific electrical carbon intensity, :math:`E^{C}_{spec,Q}` is calculated by the `add_specific_electrical_carbon_intensity` method as

    .. math::
    
        E^{C}_{spec,Q} = \frac{f_{eci} C_{el,tot}}{Q}

Annual Water Production
+++++++++++++++++++++++

For a given volumetric flow `Q`, the annual water production, :math:`W^{A}_{Q}` is calculated by the `add_annual_water_production` method as

    .. math::
    
        W^{A}_{Q} = f_{util} Q


Default Costing Methods
-----------------------

While the expectation is that unit models use the self-registration process noted above, for interoperability with IDAES unit models the WaterTAPCostingBlockData class defines two default costing methods for IDAES unit models:

* Mixer - :py:func:`watertap.costing.unit_models.mixer.cost_mixer`
* HeatExchanger - :py:func:`watertap.costing.unit_models.heat_exchanger.cost_heat_exchanger`


Class Documentation
-------------------

* :class:`WaterTAPCostingBlockData`
