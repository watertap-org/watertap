.. _how_to_use_watertap_costing:

How to use WaterTAP Costing
===========================

Overview
--------

This guide shows you how to add the WaterTAP costing package, add unit model costing, create a custom costing model, and access costing results.
Additional details on the WaterTAP costing package, including equations and default parameter values, can be found in the :ref:`official documentation<WaterTAPCostingBlockData>`.

.. Note that this guide will use the primary costing package :ref:`the primary costing package<watertap_costing>` but identical steps could be taken with :ref:`the detailed costing package<watertap_costing_detailed>`.

How To
------

.. testsetup:: python

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

    def real_build():

        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = SeawaterParameterBlock()

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
                m.fs.chem_addition.chem_dose
                * m.fs.feed.properties[0].flow_vol_phase["Liq"],
                to_units=pyunits.kg / pyunits.s,
            ),
            doc="Mass flow of chemical addition",
        )
        # 

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

        # Add system recovery constraint
        m.fs.recovery = Var(initialize=0.5, bounds=(0, 1), doc="Recovery fraction")
        m.fs.recovery_constraint = Constraint(
            expr=m.fs.recovery * m.fs.feed.properties[0].flow_vol_phase["Liq"]
            == m.fs.product.properties[0].flow_vol_phase["Liq"]
        )

        # Make connections
        m.fs.feed_to_pump1 = Arc(source=m.fs.feed.outlet, destination=m.fs.pump1.inlet)
        m.fs.pump1_to_chem_addition = Arc(
            source=m.fs.pump1.outlet, destination=m.fs.chem_addition.inlet
        )
        m.fs.chem_addition_to_pump2 = Arc(
            source=m.fs.chem_addition.outlet, destination=m.fs.pump2.inlet
        )
        m.fs.pump2_to_RO = Arc(source=m.fs.pump2.outlet, destination=m.fs.RO.inlet)
        m.fs.RO_to_ERD = Arc(source=m.fs.RO.retentate, destination=m.fs.ERD.inlet)
        m.fs.RO_to_product = Arc(source=m.fs.RO.permeate, destination=m.fs.product.inlet)
        m.fs.ERD_to_brine = Arc(source=m.fs.ERD.outlet, destination=m.fs.brine.inlet)

        TransformationFactory("network.expand_arcs").apply_to(m)

        # Set scaling
        flow_vol = 1 * pyunits.liter / pyunits.second
        conc = 35 * pyunits.gram / pyunits.liter
        rho = 1000 * pyunits.kg / pyunits.m**3

        mass_flow_water = pyunits.convert(flow_vol * rho, to_units=pyunits.kg / pyunits.s)
        mass_flow_salt = pyunits.convert(flow_vol * conc, to_units=pyunits.kg / pyunits.s)
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1 / value(mass_flow_water), index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp",
            1 / value(mass_flow_salt),
            index=("Liq", "TDS"),
        )

        for pump in [m.fs.pump1, m.fs.pump2]:
            iscale.set_scaling_factor(pump.control_volume.work, 1e-3)
            iscale.set_scaling_factor(
                pump.control_volume.properties_out[0].flow_vol_phase["Liq"], 1
            )
            iscale.set_scaling_factor(pump.work_fluid[0], 1)

        iscale.set_scaling_factor(m.fs.RO.area, 1e-2)
        iscale.set_scaling_factor(m.fs.RO.width, 1)
        iscale.set_scaling_factor(m.fs.RO.length, 1)

        iscale.calculate_scaling_factors(m)

        # Set operating conditions
        m.fs.feed.properties.calculate_state(
            var_args={
                ("flow_vol_phase", "Liq"): flow_vol,
                ("conc_mass_phase_comp", ("Liq", "TDS")): conc,
                ("pressure", None): 101325,
                ("temperature", None): 298,
            },
            hold_state=True,
        )

        m.fs.pump1.control_volume.properties_out[0].pressure.fix(1 * pyunits.bar)
        m.fs.pump1.efficiency_pump.fix(0.8)

        m.fs.chem_addition.chem_dose.fix(10)

        m.fs.pump2.control_volume.deltaP.fix(60 * pyunits.bar)
        m.fs.pump2.efficiency_pump.fix(0.8)

        m.fs.RO.A_comp.fix(4.2e-12)
        m.fs.RO.B_comp.fix(3.5e-8)
        m.fs.RO.feed_side.channel_height.fix(1e-3)
        m.fs.RO.feed_side.spacer_porosity.fix(0.9)
        m.fs.RO.permeate.pressure[0].fix(101325)
        m.fs.RO.width.fix(5)
        m.fs.RO.area.fix(30)

        m.fs.ERD.efficiency_pump.fix(0.95)
        m.fs.ERD.control_volume.properties_out[0].pressure.fix(101325)

        m.fs.feed.initialize()
        propagate_state(m.fs.feed_to_pump1)

        m.fs.pump1.initialize()
        propagate_state(m.fs.pump1_to_chem_addition)
        m.fs.chem_addition.initialize()
        propagate_state(m.fs.chem_addition_to_pump2)
        m.fs.pump2.initialize()
        propagate_state(m.fs.pump2_to_RO)

        m.fs.RO.initialize()
        propagate_state(m.fs.RO_to_ERD)
        propagate_state(m.fs.RO_to_product)

        m.fs.ERD.initialize()
        propagate_state(m.fs.ERD_to_brine)

        m.fs.product.initialize()
        m.fs.brine.initialize()

        results = solver.solve(m, tee=True)
        assert_optimal_termination(results)

        m.fs.pump1.control_volume.properties_out[0].pressure.unfix()
        m.fs.pump2.control_volume.properties_out[0].pressure.unfix()
        m.fs.RO.area.unfix()
        m.fs.RO.width.unfix()
        # m.fs.RO.length.fix(1)

        m.fs.recovery.fix(0.5)
        results = solver.solve(m, tee=True)
        assert_optimal_termination(results)
    
    m = real_build()


The :ref:`WaterTAP costing package<watertap_costing>` can be added to any WaterTAP flowsheet but is not required to run a WaterTAP model.
For this guide, we will consider a flowsheet with these units in the following order:

- Feed block
- Pump 1
- Chemical addition unit
- Pump 2 
- Reverse Osmosis unit
- ERD unit
- Product block
- Brine block


How-To Add WaterTAP Costing to a Flowsheet
*******************************************

Below is a code snippet to create the flowsheet.


.. code:: python


    def build():

        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = SeawaterParameterBlock()
        # add the WaterTAP costing package to the flowsheet
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


Adding the WaterTAP costing package to the flowsheet is done by simply by creating an instance of the ``WaterTAPCosting`` class and assigning it to a flowsheet attribute. Convention is to name this attribute ``m.fs.costing``.
This is referred to as the "flowsheet costing block" (contrasted with a "unit model costing block" discussed later). At this point, the flowsheet costing block is not doing anything. 
It contains instructions to aggregate the costs from the individual unit model costing blocks into overall flowsheet-level costs.
To get costing results, we need to add costing blocks to each unit model.

How-To Create Custom Costing Method 
********************************************

Note that the chemical addition unit in this flowsheet is not an existing WaterTAP model but a :doc:`state junction <idaes:reference_guides/model_libraries/generic/unit_models/statejunction>` (passthrough model) 
with an added ``dose`` variable and expression ``chem_mass_flow`` to calculate the mass flow of the chemical. 
Unlike other unit models on this flowsheet (like RO, which has a built-in costing method defined in the unit model file), this chemical addition unit model does not have a built-in costing method. If we want our 
system results to reflect the cost of chemical addition, we need to create a custom costing method:sup:`1`.

The code below shows an example of how to build a custom costing method. This is the general structure of :ref:`costing methods for existing WaterTAP unit models<detailed_unit_model_costing>` that are in the *watertap/costing/unit_models* directory.

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
    - Create a variable for calculating fixed operating cost (``factor_equip_replacement``) as fraction of the capital cost per year
    - Create a variable for the unit cost of the chemical (``chemical_unit_cost``)
    - Register a "chemical" flow type with the costing package and assign cost as ``chemical_unit_cost``. This will be used to calculate operating costs based on the mass flow of the chemical addition.

2. A function to build the costing model (``chem_addition_costing``), which is decorated with the ``@register_costing_parameter_block`` decorator. This function creates the costing variables and constraints needed to calculate capital and operating costs, and also defines the variable cost calculations using the ``cost_flow`` method of the costing package:sup:`2`.

    - The first argument to the ``@register_costing_parameter_block`` decorator is the function that builds the costing parameter block, and the second argument is the desired name for the parameter block on the flowsheet costing block.
    This is the name that will be used to access the parameters for this costing method from the flowsheet costing block.
    In this example, the parameter block is named "chem_addition" and is accessed with ``m.fs.costing.chem_addition``.
    - Within the costing method function, we first create the necessary costing variables (capital cost and fixed operating cost:sup:`3`). Then we define the constraints that calculate capital and operating cost based on the parameters defined in the parameter block. 
    - The ``cost_flow`` method is used to aggregate flows of the same type across multiple units (most commonly this is done with chemical and electricty flows).
    - Costing methods that calculate capital costs must provide a capital cost factor ("TIC", "TPEC", or ``None``) to be used to calculate direct and indirect capital costs.

.. important::

    :sup:`1` Users can create a custom costing method for *any* unit model in the same way. The flowsheet costing block will use the costing method passed via the ``costing_method`` argument when creating the ``UnitModelCostingBlock`` before using the costing method defined by the unit model's ``default_costing_method`` attribute.

    :sup:`2` The flow cost set via the ``register_flow_type`` method multiplied by the flow passed to the ``cost_flow`` method *must* be convertable to cost / time units because it is considered a variable operating cost for that flow type.

    :sup:`3` For proper aggregation of capital and operating costs, the flowsheet costing block requires the following naming conventions:

        - The capital cost variable must be named ``capital_cost`` and constraint ``capital_cost_constraint``.
        - The fixed operating cost variable must be named ``fixed_operating_cost`` and constraint ``fixed_operating_cost_constraint``.

        For this reason, the imported utility functions ``make_capital_cost_var`` and ``make_fixed_operating_cost_var`` should be used.

How-To Add Unit Model Costing
******************************

At this point in our flowsheet build, we have our flowsheet costing block ``m.fs.costing`` and our custom costing method for chemical addition defined, but we have not yet added costing to any of our unit models.
To add costing to a unit model, we assign a ``UnitModelCostingBlock`` to the unit model's ``costing`` attribute and specify the flowsheet costing block as an argument.
For the WaterTAP unit models on this flowsheet (RO, pumps, and ERD), that is all that is needed because their custom costing methods are already assigned to their ``default_costing_method`` attribute. 
However, if we were to do the same for the chemical addition unit, an error would be raised because the flowsheet costing block cannot find the custom costing method for chemical addition. 
In this case, we must pass the costing method via the ``costing_method`` argument when creating the ``UnitModelCostingBlock``.

The following snippet shows how to add costing to the unit models on our flowsheet, including the custom costing method for chemical addition.
Costing proceeds in the following steps:

1. Assign a ``UnitModelCostingBlock`` to the ``costing`` attribute of each unit model, passing the flowsheet costing block as an argument. For the chemical addition unit, also pass the custom costing method via the ``costing_method`` argument.
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


.. testsetup::

    add_costing(m)

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

.. testsetup::

    results = solver.solve(m)
    assert_optimal_termination(results)

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

.. testcode:: python

    from pyomo.environ import value

    # Example of accessing system-level costing results
    total_capital_cost = value(m.fs.costing.total_capital_cost)
    total_operating_cost = value(m.fs.costing.total_operating_cost)
    LCOW = value(m.fs.costing.LCOW)
    SEC = value(m.fs.costing.SEC)


In this example, we have registered flows of ``"electricity"`` and ``"chemical"`` for costing purposes. All registered flows are accessed through the ``aggregate_flow_*`` variable(s) on the costing block, 
where the ``*`` is replaced by the name of the registered flow (e.g., ``electricity`` or ``chemical``).

.. testcode:: python

    # Example of accessing aggregate flow results
    electricity_flow = value(m.fs.costing.aggregate_flow_electricity) # kW
    chemical_flow = value(m.fs.costing.aggregate_flow_chemical) # kg/s

Costs for these flows are accessed via the ``aggregate_flow_costs`` indexed variable on the costing block. Similarly, the index name is the name of the flow.

.. testcode:: python

    # Example of accessing aggregate flow costs
    electricity_cost = value(m.fs.costing.aggregate_flow_costs["electricity"])
    chemical_cost = value(m.fs.costing.aggregate_flow_costs["chemical"])


LCOW and SEC are further broken down into the contributions from each unit model, unit model class, and flow type in the following expressions.
Further descriptions of these breakdowns, how each expression is indexed, and the equations are presented in the :ref:`LCOW <aggregate_metric_LCOW>` and :ref:`SEC <aggregate_metric_SEC>` sections of the :ref:`costing package<watertap_costing>` documentation.

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

.. testcode:: python

    # Example of accessing breakdown of LCOW results
    pump1_direct_capex = value(m.fs.costing.LCOW_component_direct_capex["fs.pump1"])
    pump2_direct_capex = value(m.fs.costing.LCOW_component_direct_capex["fs.pump2"])
    electricity_variable_opex = value(m.fs.costing.LCOW_aggregate_variable_opex["electricity"])
    chemical_variable_opex = value(m.fs.costing.LCOW_aggregate_variable_opex["chemical"])

    # Example of accessing breakdown of SEC results
    pump1_SEC = value(m.fs.costing.SEC_component["fs.pump1"])
    pump2_SEC = value(m.fs.costing.SEC_component["fs.pump2"])


Alternatively, users can use the ``.display()`` method:sup:`1` on these components to view their values directly.

.. testcode:: python

    m.fs.costing.LCOW.display()
    m.fs.costing.SEC.display()
    m.fs.costing.aggregate_flow_costs.display()
    m.fs.costing.LCOW_component_direct_capex.display()
    m.fs.costing.LCOW_component_indirect_capex.display()
    m.fs.costing.LCOW_component_fixed_opex.display()
    m.fs.costing.LCOW_component_variable_opex.display()
    m.fs.costing.SEC_component.display()

The output below is only for demonstration purposes and may differ depending on the specific model and parameters used.

.. code-block:: none

    LCOW : Size=1
        Key  : Value
        None : 0.7431456904640147
    SEC : Size=1
        Key  : Value
        None : 2.6082564743928325
    aggregate_flow_costs : Size=2, Index=fs.costing.used_flows, Units=USD_2023/a
        Key         : Lower : Value              : Upper : Fixed : Stale : Domain
        chemical :  None : 31.557600005338784 :  None : False : False :  Reals
        electricity :  None :  3811.372904540384 :  None : False : False :  Reals
    LCOW_component_direct_capex : Size=5
        Key              : Value
                fs.pump1 : 2.5358855332579296e-09
                fs.pump2 :     0.1333160045572335
                fs.RO :   0.026985738454967153
                fs.ERD :   0.008982409079407897
        fs.chem_addition :  0.0011924289089755807
    LCOW_component_indirect_capex : Size=5
        Key              : Value
                fs.pump1 : 2.5358855332579296e-09
                fs.pump2 :     0.1333160045572335
                fs.RO :   0.026985738454967153
                fs.ERD :   0.008982409079407897
        fs.chem_addition :  0.0011924289089755807
    LCOW_component_fixed_opex : Size=5
        Key              : Value
                fs.pump1 : 1.5215313199547577e-09
                fs.pump2 :     0.0799896027343401
                fs.RO :    0.07016291998291461
                fs.ERD :   0.005389445447644737
        fs.chem_addition :    0.00310031516333651
    LCOW_component_variable_opex : Size=6
        Key                                         : Value
                                        fs.pump1 :  7.33995321187609e-09
                                        fs.pump2 :    0.3858743713038191
                                            fs.RO :                   0.0
                                            fs.ERD :  -0.14432414010246428
                                fs.chem_addition :                   0.0
        fs.feed.properties[0.0].flow_vol_phase[Liq] : 0.0019999999999999996
    SEC_component : Size=3
        Key      : Value
        fs.pump1 : 7.925672357944015e-08
        fs.pump2 :     4.166666666666666
        fs.ERD :   -1.5584102715305574

.. note::

    :sup:`1` The contents of many Pyomo objects (e.g., ``Var``, ``Constraint``, ``Param``, ``Expression``) can be accessed with the ``.display()`` method. If the object is indexed, the ``.display()`` method will show all indices and their corresponding values, but this method cannot be used on individual indexes (i.e., ``m.fs.costing.aggregate_flow_costs["electricity"].display()`` will raise an error.)
