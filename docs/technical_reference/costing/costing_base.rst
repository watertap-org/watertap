.. _WaterTAPCostingBlockData:

WaterTAP Costing Framework
==========================

.. index::
   pair: watertap.costing.watertap_costing;WaterTAPCostingBlockData

.. currentmodule:: watertap.costing.watertap_costing

The WaterTAP Costing Base class and utility functions contain extensions, methods, variables, and constraints common to all WaterTAP Costing Packages, and which would be useful for creating custom costing packages for WaterTAP.
An example of using WaterTAP costing in a flowsheet is provided in the :ref:`How to use WaterTAP Costing <how_to_use_watertap_costing>` guide.

.. _extensions_over_idaes_costing_framework:

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

.. _common_global_costing_parameters:

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

        f_{crf} = \frac{ f_{wacc}\,(1 + f_{wacc}) ^ L}{ (1 + f_{wacc}) ^ L - 1}

Therefore, exactly two of the variables ``capital_recovery_factor``, ``plant_lifetime`` and ``wacc`` must be fixed. By default, ``plant_lifetime`` and ``wacc`` are fixed
and ``capital_recovery_factor`` is calculated.


The process-wide costs described below rely on two other factors that must be supplied by the derived class: the total investment factor and the maintenance-labor-chemical factor.

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
Total annualized cost                           :math:`C_{annual}`    ``total_annualized_costs``            Total cost on an annualized basis
Aggregate electricity cost                      :math:`C_{el,tot}`    ``aggregate_electricity_cost``        Sum of all electricity costs
=============================================  ====================  =====================================  ==============================================================================


Costing Calculations
--------------------

Total annualized cost is a simple function of the annualized capital cost and the annualized operating cost:

    .. math::
 
        C_{annual} = f_{crf}\,C_{ca,tot} + C_{op,tot}

The total capital cost is a simple factor of the sum of the unit model capital costs:

    .. math::

        C_{ca,tot} = f_{toti}\,C_{ca,u}

The total operating cost is the sum of the fixed and variable operating costs:

    .. math::

        C_{op,tot} = C_{op,fix} + C_{op,var}

The total fixed operating cost :math:`C_{op,fix}` is the sum of the maintenance, labor, and chemical operating costs, :math:`C_{mlc}`, and the total fixed operating costs from the unit models, :math:`C_{fop,u}`:

   .. math::

        C_{op,fix} = C_{mlc} + C_{fop,u}

Where the maintenance-labor-chemical operating cost :math:`C_{mlc}` is defined as:

   .. math::

        C_{mlc} = f_{mlc}\,C_{ca,tot}
  
The total variable operating cost is the sum of the total variable operating cost from the unit models, :math:`C_{vop,u}` plus the sum of the flow costs, :math:`C_{flow,tot}` times the plant utilization factor :math:`f_{util}`:

   .. math::

        C_{op,var} = C_{vop,u} + f_{util}\,C_{flow,tot}


Aggregate Metrics
------------------

Built-in methods can be used to add expressions for common aggregate metrics used in technoeconomic analyses of water systems.
The following methods can be used to add different metrics to the costing block:

.. csv-table::
   :header: "Method", "Default Expression Name", "Description"

    "``add_levelized_cost``", "``levelized_cost``", "Adds a levelized cost expression and its component breakdowns"
    "``add_LCOW``", "``LCOW``", "Adds LCOW expression and its component breakdowns"
    "``add_specific_energy_consumption``", "``specific_energy_consumption``", "Adds a specific energy consumption expression and its component breakdown"
    "``add_specific_electrical_carbon_intensity``", "``specific_electrical_carbon_intensity``", "Adds a specific electrical carbon intensity expression and its component breakdown"
    "``add_process_throughput``", "``annual_process_throughput``", "Adds process throughput over a specified period (default is annual)"
    "``add_annual_water_production``", "``annual_water_production``", "Adds annual water production expression"
    "``add_flow_component_breakdown``", "``*_component``", "Adds a flow component breakdown expression"

These methods accept a flow rate :math:`Q` (with units of quantity per time) as the basis for the calculation. Two arguments are available to specify the flow basis and the corresponding units for the flow basis:

- ``flow_basis`` (optional): The basis for the flow rate, either ``"volumetric"``, ``"mass"``, or ``"energy"``.
- ``flow_basis_units`` (optional): The units for the flow basis (e.g., ``flow_basis_units = pyo.units.kg``).

The flow basis is inferred from the flow rate units unless ``flow_basis`` or ``flow_basis_units`` is provided. The default units for the flow basis when ``flow_basis_units`` is not specified are:

.. csv-table::
   :header: "``flow_rate`` Inferred Basis", "Specified ``flow_basis``", "Default Units"

    "Volumetric", ``"volumetric"``, "m\ :sup:`3`"
    "Mass", ``"mass"``, "kg"
    "Energy", ``"energy"``, "kWh"

.. _aggregate_metric_LCOW:

Levelized Cost Metrics
++++++++++++++++++++++

For a given flow rate, an expression for the levelized cost is added by the ``add_levelized_cost`` method. The method has four arguments:

- ``flow_rate`` (required): the flow rate used as the basis for the levelized cost calculation
- ``name`` (optional): custom name for the created expression. If not provided, ``levelized_cost`` is used as the default name
- ``flow_basis`` (optional): Specifies the basis for the flow rate (``"volumetric"``, ``"mass"``, or ``"energy"``)
- ``flow_basis_units`` (optional): explicit units for the flow basis

The ``add_LCOW`` method is a convenience wrapper around ``add_levelized_cost`` that adds the levelized cost of water :math:`LCOW_Q` expression for a volumetric flow rate :math:`Q`, setting ``flow_basis="volumetric"`` and ``LCOW`` as the default expression name. The LCOW expression is calculated as:

    .. math::
  
        LCOW_Q = \frac{f_{crf}\,C_{ca,tot} + C_{op,tot}}{f_{util}\,Q}

In addition to creating the levelized cost expression at the system level, the ``add_levelized_cost`` method will create the following indexed expressions 
to further break down the cost components contributing to the levelized cost:

.. csv-table::
   :header: "Description", "Default Expression Name :sup:`1`", "Index", "Equation :sup:`2`"

    "Direct capital expenditure by flowsheet component", "``*_component_direct_capex``", "Unit model flowsheet name :sup:`3` ", ":math:`\cfrac{f_{crf}\,C_{dir,i}}{f_{util}\,Q}`"
    "Indirect capital expenditure by flowsheet component", "``*_component_indirect_capex``", "Unit model flowsheet name", ":math:`\cfrac{f_{crf}\,C_{indir,i}}{f_{util}\,Q}`"
    "Fixed operating expenditure by flowsheet component", "``*_component_fixed_opex``", "Unit model flowsheet name", ":math:`\cfrac{f_{crf}\,C_{fop,i}}{f_{util}\,Q}`"
    "Variable operating expenditure by flowsheet component", "``*_component_variable_opex``", "Unit model flowsheet name *or* flow name :sup:`4`", ":math:`\cfrac{f_{crf}\,C_{vop,i}}{f_{util}\,Q}`"
    "Aggregate direct capital expenditure by unit type", "``*_aggregate_direct_capex``", "Unit model class name :sup:`5`", ":math:`\cfrac{f_{crf} \sum C_{dir,u}}{f_{util}\,Q}`"
    "Aggregate indirect capital expenditure by unit type", "``*_aggregate_indirect_capex``", "Unit model class name", ":math:`\cfrac{f_{crf} \sum C_{indir,u}}{f_{util}\,Q}`"
    "Aggregate fixed operating expenditure by unit type", "``*_aggregate_fixed_opex``", "Unit model class name", ":math:`\cfrac{f_{crf} \sum C_{fop,u}}{f_{util}\,Q}`"
    "Aggregate variable operating expenditure by unit type", "``*_aggregate_variable_opex``", "Unit model class name *or* flow name", ":math:`\cfrac{f_{crf} \sum C_{vop,u}}{f_{util}\,Q}`"

.. note::
    :sup:`1` The default expression names prepend the method argument `name` to the extended variable name; e.g., ``add_levelized_cost(flow_rate, name="MyLCOW")``, will result in ``MyLCOW_component_direct_capex``.

    :sup:`2` The index :math:`i` refers to individual unit model instances on the flowsheet, while :math:`u` refers to unit model classes.

    :sup:`3` The unit model flowsheet name is the name assigned to the unit model when it is added to the flowsheet (e.g., ``m.fs.unit1 = MyUnitModel()`` would have a flowsheet name of ``"fs.unit1"``).

    :sup:`4` The flow name is the name used when registering the flow with the costing package (e.g., ``m.fs.costing.register_flow_type("foobaz", foobaz_unit_cost)`` would have a flow name of ``"foobaz"``).

    :sup:`5` The unit model class name is the string representation of the class used to define the unit model (e.g., ``"ReverseOsmosis0D"``, ``"Pump"``).

Note the difference between the "component" and "aggregate" expressions: the component expressions break down costs by individual unit model instances,
while the aggregate expressions sum costs by unit model class. So, for the levelized cost of water, if there are multiple pumps on the flowsheet, the individual contributions 
to LCOW from each pump would be available in the ``LCOW_component_*`` expressions, while the total contribution from all pumps would be available as ``LCOW_aggregate_*`` expressions.
The ``LCOW_component_*`` expressions are indexed by the string representation of the unit model flowsheet name.
The indexes for the ``LCOW_aggregate_*`` expressions are the unit model class name.

Importantly, both ``*_component_variable_opex`` and ``*_aggregate_variable_opex`` expressions are also indexed by flow name for registered flows.
Energy (e.g., ``"electricity"``) and material (e.g., ``"naocl"``, ``"caustic"``) flows registered with the costing package will have their variable operating costs
broken out in these expressions. This allows the user to see the contribution of individual flow costs to the overall levelized cost.

For an example of the breakdowns presented by each of these expressions, see the :ref:`how to use WaterTAP costing<how_to_use_watertap_costing>` guide.


.. _aggregate_metric_SEC:

Specific Energy Consumption (SEC)
+++++++++++++++++++++++++++++++++

For a given flow rate :math:`Q_B`, an expression for the specific energy consumption is added by the ``add_specific_energy_consumption`` method. The method has four arguments:

- ``flow_rate``: the flow rate of the stream for which the specific energy consumption is calculated
- ``name`` (optional): custom name for the created expression. If not provided, ``specific_energy_consumption`` is set as the default name
- ``flow_basis`` (optional): the basis for the flow rate (``"volumetric"``, ``"mass"``, or ``"energy"``)
- ``flow_basis_units`` (optional): explicit units for the flow basis

The calculation uses an hourly period, so the resulting expression has units of energy per unit of the selected flow basis (e.g., energy per cubic meter, energy per kilogram, etc.):

    .. math::
  
        \text{SEC}_{Q_B} = \frac{C_{el,tot}}{Q_B}

Here, :math:`Q_B` is the flow rate converted to the corresponding basis units :math:`B` per hour and :math:`C_{el,tot}` is the total power consumption.

Users can optionally provide custom names for the created expressions via the ``name`` keyword argument. For example, creating an expression called ``SEC`` on ``m.fs.costing`` based on ``flow_rate`` would be:

.. code-block:: python

    m.fs.costing.add_specific_energy_consumption(
        flow_rate,
        name="SEC",
    )

The method also creates a component breakdown expression with ``_component`` appended to the selected name. This expression is indexed by the unit model or registered flow associated with each electricity flow and is calculated as:
    
    .. math::
  
        \text{SEC}^{\text{component}}_{Q_B,i} = \frac{C_{el,i}}{Q_B}

Specific Electrical Carbon Intensity (SECI)
+++++++++++++++++++++++++++++++++++++++++++

For a given flow rate :math:`Q_B`, an expression for the specific electrical carbon intensity is added by the ``add_specific_electrical_carbon_intensity`` method.
For a given flow :math:`Q_B`, an expression for the specific electrical carbon intensity, :math:`\text{SECI}_{Q_B}` is added by the ``add_specific_electrical_carbon_intensity`` method as:

    .. math::
  
        \text{SECI}_{Q_B} = \frac{f_{eci}\,C_{el,tot}}{Q_B}

Additionally, the specific electrical carbon intensity will be broken down by unit model. An expression is created with ``_component`` appended to the name provided by the user (or ``specific_electrical_carbon_intensity`` by default).
This expression is indexed by unit model flowsheet name and is calculated as:
    
    .. math::
    
            \text{SECI}^{\text{component}}_{Q_B,i} = \frac{f_{eci}\,C_{el,i}}{Q_B}

Process Throughput
++++++++++++++++++

For a given flow rate :math:`Q`, an expression for process throughput over a period :math:`T` is added by the ``add_process_throughput`` method as:

    .. math::
   
        \mathrm{Throughput}^T_Q = f_{util}\,Q

The method has the following arguments:

- ``flow_rate`` (required): flow rate to be used in calculating throughput
- ``name`` (optional): name for the throughput expression (default: ``annual_process_throughput``)
- ``flow_basis`` (optional): basis for the flow rate, either ``"volumetric"``, ``"mass"``, or ``"energy"``
- ``flow_basis_units`` (optional): explicit flow units (e.g., m\ :sup:`3`, kg, kWh)
- ``period`` (optional): reporting period for throughput (e.g., year, month, day). Defaults to year

The ``add_annual_water_production`` method remains available as a convenience wrapper that creates annual volumetric throughput with the default name ``annual_water_production``. For a given volumetric flow rate :math:`Q`, the annual water production, :math:`\mathrm{W}^A_Q`, is calculated as: 

    .. math::

        \mathrm{W}^A_Q = f_{util}\,Q

Flow Breakdowns
++++++++++++++++

An additional method on the WaterTAP costing block is ``add_flow_component_breakdown``. 
This allows the user to break down the costs associated with individual registered flows for a specific flow type (e.g., electricity, chemicals) per unit (e.g., cubic meter, kilogram) of product flow :math:`Q_p`. 
For a given registered flow type :math:`x` the flow component breakdown :math:`\text{FCB}_{u}` by flow source :math:`u` is calculated as

    .. math::

        \text{FCB}_{u} = \frac{F_{x,u}\,M_f}{Q_p}

Where :math:`F_{x,u}` is the flow of :math:`x` from source :math:`u`, :math:`M_f` is an optional multiplier, and :math:`Q_p` is a specified flow rate in the chosen flow basis units or inferred units from the flow rate.
:math:`M_f` must have units that, when multiplied with the units for :math:`F_{x,u}`, result in a rate (i.e., units per time). For example, if the flow rate was electricity (units of kW),
the multiplier could be an electrical carbon intensity (units of kg/kWh) and the resulting units would be kg/hr.

The method has three required arguments and five optional arguments:

- ``flow_name`` (required): string for a registered flow type
- ``name`` (required): base name appended with ``_component`` for expression name
- ``flow_rate`` (required): flow rate to be used for normalization
- ``flow_basis`` (optional): flow basis, either ``"volumetric"``, ``"mass"``, or ``"energy"``
- ``flow_basis_units`` (optional): explicit units for the flow rate
- ``period`` (optional): time period for normalization (default is ``base_period``)
- ``utilization_factor`` (optional): utilization factor for the flow (default is the costing block's ``utilization_factor``)
- ``multiplier`` (optional): multiplier for the flow (default is 1.0)

To create a breakdown of costs for ``bazchem`` used per hour per cubic meter of water, the following will create 
an expression ``m.fs.costing.bazchem_flow_component`` indexed to every unit that is associated with a flow of ``bazchem``.

.. code-block:: python

    m.fs.costing.add_flow_component_breakdown(
        "bazchem", "bazchem_flow", flow_rate, period=pyunits.hour
    )

.. note::

    The ``add_flow_component_breakdown`` method will try to find the unit associated with each registered flow automatically. If it can't, the logger will print a warning and the index for the unidentified flow will be the name of the flow expression.

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
