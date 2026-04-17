.. _how_to_use_watertap_costing:

How to use WaterTAP Costing
===========================

Overview
--------

This guide shows you how to add the WaterTAP costing package, add unit model costing, create a custom costing model, and access costing results. 

How To
------

.. testsetup::

    # quiet idaes logs
    import idaes.logger as idaeslogger
    idaeslogger.getLogger('ideas.core').setLevel('CRITICAL')
    idaeslogger.getLogger('ideas.core.util.scaling').setLevel('CRITICAL')
    idaeslogger.getLogger('idaes.init').setLevel('CRITICAL')

    from pyomo.environ import (
        ConcreteModel,
        Constraint,
        Var,
        Expression,
        Objective,
        TransformationFactory,
        assert_optimal_termination,
        value,
    )
    from pyomo.environ import units as pyunits
    from pyomo.network import Arc

    from idaes.core import FlowsheetBlock, UnitModelCostingBlock
    import idaes.core.util.scaling as iscale
    from idaes.core.util.model_statistics import degrees_of_freedom
    from idaes.core.util.initialization import propagate_state
    from idaes.models.unit_models import Feed, Product, StateJunction

    from watertap.costing import WaterTAPCosting
    from watertap.unit_models.reverse_osmosis_1D import (
        ReverseOsmosis1D,
        ConcentrationPolarizationType,
        MassTransferCoefficient,
        PressureChangeType,
    )
    from watertap.costing.util import (
        register_costing_parameter_block,
        make_capital_cost_var,
        make_fixed_operating_cost_var,
    )
    from watertap.unit_models.pressure_changer import Pump, EnergyRecoveryDevice
    from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock

    from watertap.core.solvers import get_solver

    solver = get_solver()

The :ref:`WaterTAP costing package<watertap_costing>` can be added to any WaterTAP flowsheet but is not required to run a WaterTAP model.
For this How-To, we will build a flowsheet with:

- Feed block
- Pump 1
- Chemical addition unit
- Pump 2 
- Reverse Osmosis unit
- ERD unit
- Product block
- Brine block


How-To Add WaterTAP Costing to a Flowsheet
===========================================

Below is the code to build this flowsheet. 


.. testcode:: 


    def build():

        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = SeawaterParameterBlock()
        m.fs.costing = WaterTAPCosting()

        m.fs.feed = Feed(property_package=m.fs.properties)
        m.fs.pump1 = Pump(property_package=m.fs.properties)

        m.fs.chem_addition = StateJunction(property_package=m.fs.properties)
        m.fs.chem_addition.chem_dose = Var(
            initialize=10,
            bounds=(0, None),
            units=pyunits.mg / pyunits.liter,
            doc="Chemical addition dose",
        )
        m.fs.chem_addition.chem_mass_flow = Expression(
            expr=pyunits.convert(
                m.fs.chem_addition.chem_dose * m.fs.feed.properties[0].flow_vol_phase["Liq"],
                to_units=pyunits.kg / pyunits.s,
            ),
            doc="Mass flow of chemical addition",
        )

        m.fs.RO = ReverseOsmosis1D(
            property_package=m.fs.properties,
            has_pressure_change=True,
            pressure_change_type=PressureChangeType.calculated,
            mass_transfer_coefficient=MassTransferCoefficient.calculated,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
            transformation_scheme="BACKWARD",
            transformation_method="dae.finite_difference",
            finite_elements=10,
            has_full_reporting=True,
        )

        m.fs.pump2 = Pump(property_package=m.fs.properties)
        m.fs.ERD = EnergyRecoveryDevice(property_package=m.fs.properties)
        m.fs.product = Product(property_package=m.fs.properties)
        m.fs.brine = Product(property_package=m.fs.properties)

        return m


Adding the WaterTAP costing package to the flowsheet is done by simply by creating an instance of the ``WaterTAPCosting`` class and assigning it to a flowsheet attribute. Convention is to call it ``m.fs.costing``.
This is referred to as the "flowsheet costing block" (contrasted with a "unit model costing block" discussed later). At this point, the flowsheet costing block is not doing anything. To get costing results, we need to add costing blocks to each unit model.

How-To Create Custom Costing Method
===================================

Note that the chemical addition unit is not an existing WaterTAP model but a `state junction <https://idaes-pse.readthedocs.io/en/stable/reference_guides/model_libraries/generic/unit_models/statejunction.html>`_ (passthrough) 
with an added dose variable and expression to calculate the mass flow of the chemical. 

.. This unit model is being used here for illustrative purposes to show how to add costing to a unit model that does not have a built in costing method.

Unlike other unit models on this flowsheet (like RO, which has a built-in costing method), this chemical addition unit model does not have a built-in costing method. If we want our 
system results to reflect the cost of chemical addition, we need to create a custom costing method.

.. Above you can see we assigned the custom costing method ``chem_addition_costing`` to the chemical addition unit model property ``default_costing_method``. 

The code below shows how to build the custom costing method. This is the general structure of costing models for existing WaterTAP unit models that are in the *watertap/costing/unit_models* directory.

.. testcode:: 

    def build_chem_addition_cost_param_block(blk):

        blk.chemical_capex_base = Var(
            initialize=5e4,
            units=pyunits.USD_2020 / (pyunits.Mgallons / pyunits.day),
            doc="Base capital cost for chemical addition",
        )

        blk.chemical_capex_exponent = Var(
            initialize=0.7,
            units=pyunits.dimensionless,
            doc="Exponent for chemical addition capital cost scaling",
        )

        blk.factor_equip_replacement = Var(
            initialize=0.1,
            units=pyunits.year**-1,
            doc="Fraction of chemical addition equipment replaced per year",
        )

        blk.chemical_unit_cost = Var(
            initialize=0.1,
            units=pyunits.USD_2023 / pyunits.kg,
            doc="Unit cost of chemical addition",
        )
        costing_pkg = blk.parent_block()
        costing_pkg.register_flow_type("chemical", blk.chemical_unit_cost)


    @register_costing_parameter_block(
        build_rule=build_chem_addition_cost_param_block,
        parameter_block_name="chem_addition",
    )
    def chem_addition_costing(blk):

        make_capital_cost_var(blk)
        blk.costing_package.add_cost_factor(blk, "TIC")
        make_fixed_operating_cost_var(blk)

        capex_base = pyunits.convert(
            blk.costing_package.chem_addition.chemical_capex_base
            * blk.unit_model.properties[0].flow_vol_phase["Liq"]
            * blk.costing_package.base_currency**-1,
            to_units=pyunits.dimensionless,
        )

        blk.capital_cost_constraint = Constraint(
            expr=blk.capital_cost
            == blk.cost_factor
            * capex_base**blk.costing_package.chem_addition.chemical_capex_exponent
        )
        blk.fixed_operating_cost_constraint = Constraint(
            expr=blk.fixed_operating_cost
            == blk.costing_package.chem_addition.factor_equip_replacement * blk.capital_cost
        )

        blk.costing_package.cost_flow(blk.unit_model.chem_mass_flow, "chemical")


Custom costing methods generally consist of two functions:

1. A function to build a costing parameter block (``build_chem_addition_cost_param_block``). This function defines the parameters needed for the costing method and registers any flow types needed for variable cost calculations. In this example, the this function will:

    - Create variables for chemical addition capital cost calculation (``chemical_capex_base`` and ``chemical_capex_exponent``)
    - Create a variable for calculating fixed operating cost (``factor_equip_replacement``)
    - Create a variable for the unit cost of the chemical (``chemical_unit_cost``)
    - Register a flow type of "chemical" with the costing package with assigned cost ``chemical_unit_cost``. This will be used to calculate operating costs based on the mass flow of chemical addition.

2. A function to build the costing model (``chem_addition_costing``), which is decorated with the ``@register_costing_parameter_block`` decorator. This function defines the costing variables and constraints needed to calculate capital and operating costs, and also defines the variable cost calculations using the ``cost_flow`` method of the costing package.

    The first argument to the ``@register_costing_parameter_block`` decorator is the function that builds the costing parameter block, and the second argument is the name of the parameter block on the flowsheet costing block.
    This is the name that will be used to access the parameters for this costing method from the flowsheet costing block.
    In this example, the parameter block is accessed as ``m.fs.costing.chem_addition``.

    Within the costing method function, we first create the necessary costing variables (capital cost and fixed operating cost). Then we define the constraints that calculate capital and operating cost based on the parameters defined in the parameter block. 
    
    .. Finally, we use the ``cost_flow`` method to define variable costs based on the chemical mass flow.

    The ``cost_flow`` method is used to aggregate flows of the same type across multiple units (most commonly this is done with chemical and electricty flows).

    Costing methods that calculate capital costs must provide a capital cost factor ("TIC", "TPEC", or ``None``) to be used to calculate direct and indirect capital costs.

.. important::

    For proper aggregation of capital and operating costs, the flowsheet costing block requires the following naming conventions:
        - The capital cost variable must be named ``capital_cost`` and constraint ``capital_cost_constraint``.
        - The fixed operating cost variable must be named ``fixed_operating_cost`` and constraint ``fixed_operating_cost_constraint``.

How-To Add Unit Model Costing
=============================

At this point in our flowsheet build, we have our flowsheet costing block ``m.fs.costing`` and our custom costing method for chemical addition defined, but we have not yet added costing to any of our unit models.
To add costing to a unit model, we create a ``UnitModelCostingBlock`` to the unit model's ``costing`` attribute and specify the flowsheet costing block as an argument.
For the WaterTAP unit models on this flowsheet (RO, pumps, and ERD), that is all that is needed because their custom costing methods are already assigned to their ``default_costing_method`` attribute. 
However, if we were to do the same for the chemical addition unit, an error would be raised because the flowsheet costing block cannot find the custom costing method for chemical addition. 
In this case, we must pass the costing method via the ``costing_method`` argument when creating the ``UnitModelCostingBlock``.

The following snippet shows how to add costing to the unit models on our flowsheet, including the custom costing method for chemical addition.
Costing proceeds in the following steps:

1. Add a ``UnitModelCostingBlock`` to the ``costing`` attribute of each unit model, passing the flowsheet costing block as an argument. For the chemical addition unit, also pass the custom costing method via the ``costing_method`` argument.
2. After all unit model costing blocks have been added, call the ``cost_process`` method on the flowsheet costing block to build the costing model. This method constructs the process-level costing components based on the registered unit operations and flows.
3. Add any desired system-level costing metrics (LCOW and SEC in this example).

.. testcode::

    def add_costing(m):

        # Add unit model costing blocks 
        m.fs.pump1.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.pump2.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.RO.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.ERD.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        # Pass the custom costing method to the chemical addition costing block via the costing_method argument
        m.fs.chem_addition.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=chem_addition_costing
        )

        # Build the process costing model
        m.fs.costing.cost_process()
        # Add system-level costing metrics
        m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol_phase["Liq"])
        m.fs.costing.add_specific_energy_consumption(
            m.fs.product.properties[0].flow_vol_phase["Liq"], name="SEC"
        )

At this point in the build, you can change values of flowsheet-level costing variables (``electricity_price`` in this example), or change the ``base_currency``.
In the example below, we also create an ``Objective`` to minimize LCOW.

.. testcode::

    # Add objective
    m.fs.obj = Objective(expr=m.fs.costing.LCOW)

    # Change values of flowsheet-level costing variables if desired
    m.fs.costing.electricity_cost.fix(0.1)
    m.fs.costing.base_currency = pyunits.USD_2023

    # Change values of unit model costing parameters if desired
    m.fs.costing.chem_addition.chemical_unit_cost.fix(0.1)
    m.fs.costing.reverse_osmosis.membrane_unit_cost.fix(45)

How-To Access Costing Results
==============================

After building the costing model and solving the optimization problem, system-level costing results are on the flowsheet costing block (``m.fs.costing``) and unit-level costing results are on the ``costing`` attribute of each unit model (e.g. ``m.fs.chem_addition.costing``).

Examples of system-level costing results:

- Total capital cost: ``m.fs.costing.total_capital_cost``
- Total operating cost: ``m.fs.costing.total_operating_cost``
- Total electricity required: ``m.fs.costing.aggregate_flow_electricity``
- Total cost of electricity: ``m.fs.costing.aggregate_flow_costs["electricity"]``
- Total cost of chemicals: ``m.fs.costing.aggregate_flow_costs["chemical"]``
- LCOW: ``m.fs.costing.LCOW``
- SEC: ``m.fs.costing.SEC``

Accessing the values of these variables/expressions is done using the ``value`` function from Pyomo.

.. code::

    # Example of accessing system-level costing results
    total_capital_cost = value(m.fs.costing.total_capital_cost)
    total_operating_cost = value(m.fs.costing.total_operating_cost)
    LCOW = value(m.fs.costing.LCOW)
    SEC = value(m.fs.costing.SEC)

LCOW and SEC are further broken down into the contributions from each unit model, unit model class, and flow type in the following expressions.

- Direct capital expenditures by flowsheet component: ``m.fs.costing.LCOW_component_direct_capex``
- Indirect capital expenditures by flowsheet component: ``m.fs.costing.LCOW_component_indirect_capex``
- Fixed operating expenditures by flowsheet component: ``m.fs.costing.LCOW_component_fixed_opex``
- Variable operating expenditures by flowsheet component: ``m.fs.costing.LCOW_component_variable_opex``
- Aggregate direct capital expenditures by unit model type: ``m.fs.costing.LCOW_aggregate_direct_capex``
- Aggregate indirect capital expenditures by unit model type: ``m.fs.costing.LCOW_aggregate_indirect_capex``
- Aggregate fixed operating expenditures by unit model type: ``m.fs.costing.LCOW_aggregate_fixed_opex``
- Aggregate variable operating expenditures by unit model type: ``m.fs.costing.LCOW_aggregate_variable_opex``
- Specific energy consumption by unit model: ``m.fs.costing.SEC_component``


Values of these expressions can similarly be accessed using the ``value`` function from Pyomo.

.. code::

    # Example of accessing breakdown of LCOW results
    pump1_direct_capex = value(m.fs.costing.LCOW_component_direct_capex["fs.pump1"])
    pump2_direct_capex = value(m.fs.costing.LCOW_component_direct_capex["fs.pump2"])
    electricity_variable_opex = value(m.fs.costing.LCOW_aggregate_variable_opex["electricity"])
    chemical_variable_opex = value(m.fs.costing.LCOW_aggregate_variable_opex["chemical"])

    # Example of accessing breakdown of SEC results
    pump1_SEC = value(m.fs.costing.SEC_component["fs.pump1"])
    pump2_SEC = value(m.fs.costing.SEC_component["fs.pump2"])

Further descriptions of these breakdowns, how each expression is indexed, and the equations are presented in the :ref:`LCOW <aggregate_metric_LCOW>` and :ref:`SEC <aggregate_metric_SEC>` sections of the :ref:`costing package<watertap_costing>` documentation.
