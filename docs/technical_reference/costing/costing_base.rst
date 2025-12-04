.. _WaterTAPCostingBlockData:

WaterTAP Costing Framework
==========================

.. index::
   pair: watertap.costing.watertap_costing;WaterTAPCostingBlockData

.. currentmodule:: watertap.costing.watertap_costing

The WaterTAP Costing Base class and utility functions contain extensions, methods, and variables and constraints common to all WaterTAP Costing Packages, and which would be useful for creating custom costing packages for WaterTAP.


Extensions Over IDAES Costing Framework
---------------------------------------

The WaterTAP Costing Framework extends the functionality of the `IDAES Process Costing Framework <https://idaes-pse.readthedocs.io/en/stable/reference_guides/core/costing/costing_framework.html>`_ in several ways:

1. Unit models can self-register a default costing method by specifying a ``default_costing_method`` attribute. This allows the costing method(s) to be specified with the unit model definition.

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


2. The method ``register_flow_type`` will create a new Expression if a costing component is not already defined *and* the costing component is not constant. 
The default behavior in IDAES is to always create a new Var. This allows the user to specify intermediate values in ``register_flow_type``. 

.. testcode::

    import pyomo.environ as pyo
    import idaes.core as idc
    from watertap.costing import WaterTAPCosting

    m = pyo.ConcreteModel()
    m.fs = idc.FlowsheetBlock(dynamic=False)
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


3. Unit models specify one of the global indirect capital cost multipliers, 
either `TIC` or `TPEC` (defined below) when defining their capital costs. 
The costing package will then aggregate both direct and total capital costs.

.. testcode::

    import pyomo.environ as pyo
    import idaes.core as idc
    from watertap.costing import WaterTAPCosting

    def cost_my_unit_model(blk):
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )
        # Adds blk.cost_factor, an expression pointing
        # to the appropriate indirect capital cost adder
        # and blk.direct_capital_cost, which is a expression
        # defined to be blk.capital_cost / blk.cost_factor.
        # Valid strings are "TIC" and "TPEC", all others
        # will result in an indirect capital cost factor
        # of 1.
        blk.costing_package.add_cost_factor(blk, "TIC")

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost
            == blk.cost_factor * (42 * pyo.units.USD_2018)
        )

    @idc.declare_process_block_class("MyUnitModel")
    class MyUnitModelData(idc.UnitModelBlockData):
        pass

    m = pyo.ConcreteModel()
    m.fs = idc.FlowsheetBlock(dynamic=False)
    m.fs.costing = WaterTAPCosting()

    m.fs.my_unit = MyUnitModel()

    m.fs.my_unit.costing = idc.UnitModelCostingBlock(
        costing_method=cost_my_unit_model,
        flowsheet_costing_block=m.fs.costing,
    )
    m.fs.my_unit.costing.initialize()

    m.fs.my_unit.costing.cost_factor.pprint()
    m.fs.my_unit.costing.capital_cost.pprint()
    m.fs.my_unit.costing.direct_capital_cost.pprint()

.. testoutput::

    cost_factor : Size=1, Index=None
        Key  : Expression
        None : fs.costing.TIC
    capital_cost : Capital cost of unit operation
        Size=1, Index=None, Units=USD_2018
        Key  : Lower : Value : Upper : Fixed : Stale : Domain
        None :     0 :  84.0 :  None : False : False :  Reals
    direct_capital_cost : Size=1, Index=None
        Key  : Expression
        None : fs.my_unit.costing.capital_cost/fs.costing.TIC


4. A helper utility for defining global-level parameters specific to a unit model 
without changing the base costing package implementation.

.. testcode::

    import pyomo.environ as pyo
    import idaes.core as idc
    from watertap.costing import (
        WaterTAPCosting,
        register_costing_parameter_block,
        make_capital_cost_var,
    )

    def build_my_unit_model_param_block(blk):
        """
        This function builds the global parameters for MyUnitModel.

        This function should also register needed flows using the
        blk.parent_block().register_flow_type method on the costing package.
        """
        blk.fixed_capital_cost = pyo.Var(
            initialize=42,
            doc="Fixed capital cost for all of my units",
            units=pyo.units.USD_2020,
        )

    # This decorator ensures that the function
    # `build_my_unit_model_param_block` is only
    # added to the costing package once.
    # It registers it as a sub-block with the
    # name `my_unit`.
    @register_costing_parameter_block(
        build_rule=build_my_unit_model_param_block,
        parameter_block_name="my_unit",
    )
    def cost_my_unit_model(blk):
        """
        Cost an instance of MyUnitModel
        """
        # creates the `capital_cost` Var
        make_capital_cost_var(blk)
        blk.costing_package.add_cost_factor(blk, "TIC")

        # here we reference the `fixed_capital_cost` parameter
        # automatically added by the `register_costing_parameter_block`
        # decorator.
        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost
            == blk.cost_factor * blk.costing_package.my_unit.fixed_capital_cost
        )

    @idc.declare_process_block_class("MyUnitModel")
    class MyUnitModelData(idc.UnitModelBlockData):

        @property
        def default_costing_method(self):
            # could point to a static method on
            # this class, could be function in
            # a different module even
            return cost_my_unit_model

    m = pyo.ConcreteModel()
    m.fs = idc.FlowsheetBlock(dynamic=False)
    m.fs.costing = WaterTAPCosting()

    m.fs.my_unit_1 = MyUnitModel()

    # The `default_costing_method_attribute` on the
    # unit model is checked, and the function
    # `cost_my_unit_model` returned then build the costing block.
    # This method also adds the `my_unit` global parameter block,
    # so the global costing parameter m.fs.costing.my_unit.fixed_capital_cost
    # is the same for all instances of MyUnitModel.
    m.fs.my_unit_1.costing = idc.UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
    )

    m.fs.my_unit_2 = MyUnitModel()

    # Here everything as before, but the global parameter block
    # m.fs.costing.my_unit is not re-built.
    m.fs.my_unit_2.costing = idc.UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
    )


Costing Index and Technoeconomic Factors
----------------------------------------

Default costing indices are provided with the WaterTAP Costing Framework, 
but the user is free to modify these for their needs. Costs from year
A to year B are adjusted according to:

.. math::

    \text{Cost in B} = \text{Cost in A} \left( \frac{\text{Index at B}}{\text{Index at A}} \right)


WaterTAP uses the `Chemical Engineering Plant Cost Index <https://www.toweringskills.com/financial-analysis/cost-indices/>`_ (CEPCI) 
to account for the time-value of investments. Aggregated capital and operating costs are 
adjusted to the desired year for the model, accessible on the costing block as ``base_currency``. 
The default costing year is 2018, but the user can directly set the ``base_currency`` at 
the flowsheet level (e.g., ``m.fs.costing.base_currency = pyo.units.USD_2020``).


Common Global Costing Parameters
--------------------------------

The ``build_global_params`` method builds common cost factor parameters necessary to calculate aggregated metrics such as levelized cost of water (LCOW).
Note that the default values can be overwritten in the derived class.

=============================================  ====================  =====================================  ===============  ==============================================================================
                 Cost factor                     Variable                 Name                               Default Value    Description
=============================================  ====================  =====================================  ===============  ==============================================================================
Plant capacity utilization factor                 :math:`f_{util}`    ``utilization_factor``                 90%                Percentage of year plant is operating
Electricity price                                 :math:`P`           ``electricity_cost``                   $0.07/kWh          Electricity price in 2018 USD
Electricity carbon intensity                      :math:`f_{eci}`     ``electrical_carbon_intensity``        0.475 kg/kWh       Carbon intensity of electricity
Capital recovery factor                           :math:`f_{crf}`     ``capital_recovery_factor``            10%                Capital annualization (fraction of investment cost/year)
Plant lifetime                                    :math:`L`           ``plant_lifetime``                     30 years           Plant lifetime
Weighted average cost of capital                  :math:`f_{wacc}`    ``wacc``                               9.30734%           Average cost of capital over plant lifetime
Total purchased equipment cost (TPEC)             :math:`f_{TPEC}`    ``TPEC``                               4.121212           Common indirect capital cost multiplier for unit models
Total installed cost (TIC)                        :math:`f_{TIC}`     ``TIC``                                2.0                Common indirect capital cost multiplier for unit models
=============================================  ====================  =====================================  ===============  ==============================================================================

The relationship between the :math:`f_{crf}`, :math:`L`, and :math:`f_{wacc}` is as follows:

    .. math::

        f_{crf} = \frac{ f_{wacc} (1 + f_{wacc}) ^ L}{ (1 + f_{wacc}) ^ L - 1}

Therefore, at exactly two of the variables ``capital_recovery_factor``, ``plant_lifetime`` and ``wacc`` must be fixed. By default the variables ``plant_lifetime`` and ``wacc`` are fixed
and the variable ``capital_recovery_factor`` is calculated.


The process-wide costs described below rely on two other factors which must be supplyed by the derived class, the total investment factor and the maintenance-labor-chemical factor.

=============================================  ====================  =======================================  ===============  ==============================================================================
                 Cost factor                     Variable                 Name                                 Default Value    Description
=============================================  ====================  =======================================  ===============  ==============================================================================
Total investment factor                           :math:`f_{toti}`    ``total_investment_factor``              None            Total investment factor (investment cost / equipment cost)
Maintenance-labor-chemical factor                 :math:`f_{mlc}`     ``maintenance_labor_chemical_factor``    None            Maintenance, labor, and chemical factor (fraction of equipment cost / year)
=============================================  ====================  =======================================  ===============  ==============================================================================


Costing Process-Wide Costs
--------------------------

The WaterTAPCostingBlockData class includes variables necessary to calculate process-wide costs:

=============================================  ====================  =====================================  ==============================================================================
                 Cost                               Variable                 Name                               Description
=============================================  ====================  =====================================  ==============================================================================
Total capital cost                              :math:`C_{ca,tot}`    ``total_capital_cost``                Total capital cost
Unit capital cost                               :math:`C_{ca,u}`      ``aggregate_capital_cost``            Unit processes capital cost
Total operating cost                            :math:`C_{op,tot}`    ``total_operating_cost``              Total operating cost for unit process
Total fixed operating cost                      :math:`C_{op,fix}`    ``total_fixed_operating_cost``        Total fixed operating cost for unit process
Total variable operating cost                   :math:`C_{op,var}`    ``total_variable_operating_cost``     Total variable operating cost for unit process
Total annualized cost                           :math:`C_{annual}`    ``total_annualized_costs``            Total cost on a annualized basis
Aggregate electricity cost                      :math:`C_{el,tot}`    ``aggregate_electricity_cost``        Sum of all electricity costs
=============================================  ====================  =====================================  ==============================================================================


Costing Calculations
--------------------

Total annulized cost is a simple function of the annualized capital cost and the annualized operating cost:

    .. math::
 
        C_{annual} = f_{crf} C_{ca,tot} + C_{op,tot}

The total capital cost is a simple factor of the sum of the unit model capital costs:

    .. math::

        C_{ca,tot} = f_{toti} C_{ca,u}

The total operating cost is the sum of the fixed and variable operating costs:

    .. math::

        C_{op,tot} = C_{op,fix} + c_{op,var}

The total fixed operating cost :math:`C_{op,fix}` is the sum of the maintence, labor, and chemical operating costs, :math:`C_{mlc}` and the total fixed operating costs from the unit models, :math:`C_{fop,u}`:

   .. math::

        C_{op,fix} = C_{mlc} + C_{fop,u}

Where the maintenance-labor-chemical operating cost :math:`C_{mlc}` is defined as:

   .. math::

        C_{mlc} = f_{mlc} C_{ca,tot}
  
The total variable operating cost is the sum of the total variable operating cost from the unit models, :math:`C_{vop,u}` plus the sum of the flow costs, :math:`C_{flow,tot}` times the plant utilization factor :math:`f_{util}`:

   .. math::

        C_{op,var} = C_{vop,u} + f_{util} C_{flow,tot}


Aggregate Metrics
------------------

Built in methods can be used to add expressions for common aggregate metrics used in technoeconomic analyses of water systems.
The following metrics can be added to the costing block by calling the respective method:

.. csv-table::
   :header: "Method", "Default Expression Name", "Description"

   "``add_LCOW``", "``LCOW``", "Adds LCOW variable and constraint"
   "``add_specific_energy_consumption``", "``specific_energy_consumption``", "Adds specific energy consumption variable and constraint"
   "``add_specific_electrical_carbon_intensity``", "``specific_electrical_carbon_intensity``", "Adds specific electrical carbon intensity variable and constraint"
   "``add_annual_water_production``", "``annual_water_production``", "Adds annual water production variable and constraint" 

Each of these methods requires the user pass a volumetric flow rate :math:`Q` (with units of volume per time) to be used as the basis for the calculation.
Users can optionally provide custom names for the created expression via the `name` keyword argument. For example, creating an expression called ``SEC`` on ``m.fs.costing`` 
based on ``flow_rate`` would be:

.. code-block:: python

    m.fs.costing.add_specific_energy_consumption(
        flow_rate,
        name="SEC",
    )

Levelized Cost of Water (LCOW)
++++++++++++++++++++++++++++++

For a given volumetric flow :math:`Q`, the LCOW, :math:`LCOW_{Q}` is calculated by the ``add_LCOW`` method as

    .. math::
  
        LCOW_{Q} = \frac{f_{crf}   C_{ca,tot} + C_{op,tot}}{f_{util} Q}

In addition to creating the LCOW expression at the system level, the ``add_LCOW`` method will create the following indexed expressions 
to further break down the cost components contributing to the LCOW:

.. csv-table::
   :header: "Description", "Default Expression Name :sup:`1`", "Index", "Equation :sup:`2`"

    "Direct capital expenditure by flowsheet component", "``LCOW_component_direct_capex``", "Unit model flowsheet name ", ":math:`\cfrac{f_{crf} C_{dir,i}}{f_{util} Q}`"
    "Indirect capital expenditure by flowsheet component", "``LCOW_component_indirect_capex``", "Unit model flowsheet name", ":math:`\cfrac{f_{crf} C_{indir,i}}{f_{util} Q}`"
    "Fixed operating expenditure by flowsheet component", "``LCOW_component_fixed_opex``", "Unit model flowsheet name", ":math:`\cfrac{f_{crf} C_{fop,i}}{f_{util} Q}`"
    "Variable operating expenditure by flowsheet component", "``LCOW_component_variable_opex``", "Unit model flowsheet name *or* flow name", ":math:`\cfrac{f_{crf} C_{vop,i}}{f_{util} Q}`"
    "Aggregate direct capital expenditure by unit type", "``LCOW_aggregate_direct_capex``", "Unit model class name", ":math:`\cfrac{f_{crf} \sum C_{dir,u}}{f_{util} Q}`"
    "Aggregate indirect capital expenditure by unit type", "``LCOW_aggregate_indirect_capex``", "Unit model class name", ":math:`\cfrac{f_{crf} \sum C_{indir,u}}{f_{util} Q}`"
    "Aggregate fixed operating expenditure by unit type", "``LCOW_aggregate_fixed_opex``", "Unit model class name", ":math:`\cfrac{f_{crf} \sum C_{fop,u}}{f_{util} Q}`"
    "Aggregate variable operating expenditure by unit type", "``LCOW_aggregate_variable_opex``", "Unit model class name *or* flow name", ":math:`\cfrac{f_{crf} \sum C_{vop,u}}{f_{util} Q}`"

.. note::
    :sup:`1` The default expression names prepend the method argument `name` to the extended variable name; e.g., ``add_LCOW(flow_rate, name="MyLCOW")``, will result in ``MyLCOW_component_direct_capex``.

    :sup:`2` The index :math:`i` refers to individual unit model instances on the flowsheet, while :math:`u` refers to unit model classes.

Note, the difference between the "component" and "aggregate" expressions: the component expressions break down costs by individual unit model instances,
while the aggregate expressions sum costs by unit model class (e.g., Mixer, Pump, ReverseOsmosis0D, etc.). So, if there are multiple pumps on the flowsheet, the individual contributions 
to LCOW from each pump would be available in the ``LCOW_component_*`` expressions, while the total contribution from all pumps would be available as ``LCOW_aggregate_*`` expressions.
The ``LCOW_component_*`` expressions are indexed by the unit model flowsheet name. This is the name that is assigned when the unit model
is added to the flowsheet. For example, if you add a unit model as ``m.fs.unit1 = MyUnitModel()``, the name used in the LCOW component expressions will be ``fs.unit1``.
The indexes for the ``LCOW_aggregate_*`` expressions are by unit model class name, which is the string representation of the class used to define the unit model (e.g., "ReverseOsmosis0D", "Pump").

Importantly, both ``LCOW_component_variable_opex`` and ``LCOW_aggregate_variable_opex`` expressions are also indexed by flow name for registered flows.
Energy (e.g., ``electricity``) and material (e.g., ``naocl``, ``caustic``) flows registered with the costing package will have their variable operating costs
broken out in these expressions. This allows the user to see the contribution of individual flow costs to the overall LCOW.

For an example of the breakdowns presented by each of these expressions, consider a flowsheet that has two pump units (``m.fs.pump1``, ``m.fs.pump2``) and one chemical addition 
unit (``m.fs.chem_add``), and has registered ``electricity`` and ``anti_scalant`` flows. The user adds the LCOW via ``m.fs.costing.add_LCOW(flow_rate)``. The following expressions and indexes would be available on ``m.fs.costing``:

* ``LCOW_component_direct_capex``: ``fs.pump1``, ``fs.pump2``, ``fs.chem_add``
* ``LCOW_component_indirect_capex``: ``fs.pump1``, ``fs.pump2``, ``fs.chem_add``
* ``LCOW_component_fixed_opex``: ``fs.pump1``, ``fs.pump2``, ``fs.chem_add``
* ``LCOW_component_variable_opex``: ``fs.pump1``, ``fs.pump2``, ``fs.chem_add``, ``electricity``, ``anti_scalant``
* ``LCOW_aggregate_direct_capex``: ``Pump``, ``ChemicalAdditionZO``
* ``LCOW_aggregate_indirect_capex``: ``Pump``, ``ChemicalAdditionZO``
* ``LCOW_aggregate_fixed_opex``: ``Pump``, ``ChemicalAdditionZO``
* ``LCOW_aggregate_variable_opex``: ``Pump``, ``ChemicalAdditionZO``, ``electricity``, ``anti_scalant``

The contribution of flows to the LCOW is found as both an individual contribution to individual unit model breakdown *and* as separate entries. For example, ``electricity`` is counted both in
``LCOW_component_variable_opex['electricity']`` as well as part of the ``LCOW_component_variable_opex`` for each of the unit models that contribute electricity flow (e.g., ``LCOW_component_variable_opex['fs.pump1']`` includes the electricity cost for pump1). 
Similarly, ``electricity`` is counted both as ``LCOW_aggregate_variable_opex['electricity']`` as well as part of the ``LCOW_aggregate_variable_opex['Pump']``. For this reason, the system LCOW 
is the summation of all indexes in any of the component or aggregate expressions *except* those indexed by flow.

Specific Energy Consumption
+++++++++++++++++++++++++++

For a given volumetric flow `Q`, the specific energy consumption, :math:`E_{spec,Q}` is calculated by the ``add_specific_energy_consumption`` method as

    .. math::
  
        E_{spec,Q} = \frac{C_{el,tot}}{Q}

Additionally, the specific energy consumption will be broken down by unit model. An expression is created with ``_component`` prepended to the name provided by the user (or ``specific_energy_consumption`` by default).
This expression is indexed by unit model flowsheet name and is calculated as
    
    .. math::
  
        E^{component}_{spec,Q,i} = \frac{C_{el,i}}{Q}

For a flowsheet with two pump units (``m.fs.pump1``, ``m.fs.pump2``), calling ``m.fs.costing.add_specific_energy_consumption(flow_rate, name="SEC")`` would create ``m.fs.costing.SEC`` and ``m.fs.costing.SEC_component`` indexed by ``fs.pump1`` and ``fs.pump2``.

Specific Electrical Carbon Intensity
++++++++++++++++++++++++++++++++++++

For a given volumetric flow `Q`, the specific electrical carbon intensity, :math:`E^{C}_{spec,Q}` is calculated by the ``add_specific_electrical_carbon_intensity`` method as

    .. math::
  
        E^{C}_{spec,Q} = \frac{f_{eci} C_{el,tot}}{Q}

Additionally, the specific electrical carbon intensity will be broken down by unit model. An expression is created with ``_component`` prepended to the name provided by the user (or ``specific_electrical_carbon_intensity`` by default).
This expression is indexed by unit model flowsheet name and is calculated as
    
    .. math::
    
            E^{C,comp}_{spec,Q,i} = \frac{f_{eci} C_{el,i}}{Q}

For a flowsheet with two pump units (``m.fs.pump1``, ``m.fs.pump2``), calling ``m.fs.costing.add_specific_electrical_carbon_intensity(flow_rate, name="SECI")`` would create ``m.fs.costing.SECI`` and ``m.fs.costing.SECI_component`` indexed by ``fs.pump1`` and ``fs.pump2``.

Annual Water Production
+++++++++++++++++++++++

For a given volumetric flow `Q`, the annual water production, :math:`W^{A}_{Q}` is calculated by the ``add_annual_water_production`` method as

    .. math::
   
        W^{A}_{Q} = f_{util} Q


Default Costing Methods
-----------------------

While the expectation is that unit models use the self-registration process noted above, 
for interoperability with IDAES unit models the WaterTAPCostingBlockData class defines default costing methods for IDAES unit models:

* Mixer - :py:func:`watertap.costing.unit_models.mixer.cost_mixer`
* HeatExchanger - :py:func:`watertap.costing.unit_models.heat_exchanger.cost_heat_exchanger`
* CSTR - :py:func:`watertap.costing.unit_models.cstr.cost_cstr`
* Heater - :py:func:`watertap.costing.unit_models.heater_chiller.cost_heater_chiller`


Class Documentation
-------------------

* :class:`WaterTAPCostingBlockData`
