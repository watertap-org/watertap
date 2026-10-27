.. _how_to_use_watertap_costing:

How to use the WaterTAP costing package
=======================================

This guide will demonstrate how to use the WaterTAP costing package for a full flowsheet. It pulls together the smaller steps from other how-to guides, including

- :ref:`How to add costing packages to a flowsheet<how_to_add_watertap_costing_to_flowsheet>`,
- :ref:`How to create a custom costing method<how_to_create_custom_costing_method>`, 
- :ref:`How to cost a chemical flow<how_to_cost_a_flow>`, and
- :ref:`How to access costing results<how_to_access_costing_results>`

Additional details on the WaterTAP costing package, including equations and default parameter values, can be found in the :ref:`official WaterTAP costing documentation<WaterTAPCostingBlockData>`.


.. testsetup::

    from pyomo.environ import (
        ConcreteModel,
        Constraint,
        Param,
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

    from watertap.costing import WaterTAPCosting, WaterTAPCostingDetailed
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
    from watertap.unit_models.zero_order import ChemicalAdditionZO
    from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock

    from watertap.core.solvers import get_solver

    solver = get_solver()

    # quiet idaes logs
    import idaes.logger as idaeslogger
    idaeslogger.getLogger('idaes.core').setLevel('CRITICAL')
    idaeslogger.getLogger('idaes.core.util.scaling').setLevel('CRITICAL')
    idaeslogger.getLogger('idaes.init').setLevel('CRITICAL')


    def build_chem_addition_cost_param_block(blk):

        blk.chemical_capex_slope = Var(
            initialize=1.23e4,
            units=pyunits.USD_2020 / (pyunits.Mgallons / pyunits.day),
            doc="Base capital cost for chemical addition",
        )

        blk.chemical_capex_intercept = Var(
            initialize=5e3,
            units=pyunits.USD_2020,
            doc="Exponent for chemical addition capital cost scaling",
        )

        blk.factor_equip_replacement = Var(
            initialize=0.067,
            units=pyunits.year**-1,
            doc="Fraction of chemical addition equipment replaced per year",
        )


    def build_bazchem_cost_param_block(blk):

        blk.unit_cost = Var(
            initialize=0.089,
            units=pyunits.USD_2023 / pyunits.kg,
            doc="Unit cost of bazchem",
        )

        costing_pkg = blk.parent_block()
        costing_pkg.register_flow_type("bazchem", blk.unit_cost)


    @register_costing_parameter_block(
        build_rule=build_chem_addition_cost_param_block,
        parameter_block_name="chem_addition",
    )
    @register_costing_parameter_block(
        build_rule=build_bazchem_cost_param_block,
        parameter_block_name="bazchem",
    )
    def chem_addition_costing(blk):

        make_capital_cost_var(blk)
        blk.costing_package.add_cost_factor(blk, "TIC")
        make_fixed_operating_cost_var(blk)

        blk.capital_cost_constraint = Constraint(
            expr=blk.capital_cost
            == blk.cost_factor
            * pyunits.convert(
                blk.costing_package.chem_addition.chemical_capex_intercept,
                to_units=blk.costing_package.base_currency,
            )
            + pyunits.convert(
                blk.costing_package.chem_addition.chemical_capex_slope
                * blk.unit_model.properties[0].flow_vol_phase["Liq"],
                to_units=blk.costing_package.base_currency,
            )
        )
        blk.fixed_operating_cost_constraint = Constraint(
            expr=blk.fixed_operating_cost
            == blk.costing_package.chem_addition.factor_equip_replacement * blk.capital_cost
        )

        blk.costing_package.cost_flow(blk.unit_model.chem_mass_flow, "bazchem")


    def real_build():

        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = SeawaterParameterBlock()
        m.fs.costing = WaterTAPCosting()

        m.fs.feed = Feed(property_package=m.fs.properties)
        m.fs.pump1 = Pump(property_package=m.fs.properties)

        m.fs.chem_addition = ChemicalAdditionZO(
            property_package=m.fs.properties, process_subtype="default"
        )
        m.fs.chem_addition.chem_mass_flow = Expression(
            expr=pyunits.convert(
                m.fs.chem_addition.chemical_dosage
                * m.fs.chem_addition.properties[0].flow_vol_phase["Liq"],
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
        m.fs.RO.chem_dose = Param(
            initialize=690,
            units=pyunits.mg / pyunits.L,
            mutable=True,
            doc="Chemical dose for RO",
        )
        m.fs.RO.chem_mass_flow = Expression(
            expr=pyunits.convert(
                m.fs.RO.chem_dose
                * m.fs.RO.feed_side.properties[0, 0].flow_vol_phase["Liq"],
                to_units=pyunits.kg / pyunits.s,
            ),
            doc="Mass flow of chemical for RO",
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

        m.fs.pump1.control_volume.properties_out[0].pressure.setlb(
            10 * pyunits.bar
        )  # so the numbers in breakdown are reasonable
        m.fs.pump1.control_volume.properties_out[0].pressure.fix(10 * pyunits.bar)
        m.fs.pump1.efficiency_pump.fix(0.8)

        m.fs.chem_addition.chemical_dosage.fix(10)
        m.fs.chem_addition.ratio_in_solution.fix(1)
        m.fs.chem_addition.solution_density.fix(1000)

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

        results = solver.solve(m, tee=False)
        assert_optimal_termination(results)

        m.fs.pump1.control_volume.properties_out[0].pressure.unfix()
        m.fs.pump2.control_volume.properties_out[0].pressure.unfix()
        m.fs.RO.area.unfix()
        m.fs.RO.width.unfix()

        m.fs.recovery.fix(0.5)
        results = solver.solve(m, tee=False)
        assert_optimal_termination(results)

        return m


    def real_add_costing(m):

        m.fs.pump1.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.pump2.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.RO.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.ERD.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

        m.fs.chem_addition.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing, costing_method=chem_addition_costing
        )

        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol_phase["Liq"])
        m.fs.costing.add_specific_energy_consumption(
            m.fs.product.properties[0].flow_vol_phase["Liq"], name="SEC"
        )
        m.fs.costing.add_flow_component_breakdown(
            "bazchem",
            m.fs.product.properties[0].flow_vol_phase["Liq"],
            period=pyunits.hr,
        )
        m.fs.costing.add_flow_component_breakdown(
            "electricity",
            m.fs.product.properties[0].flow_vol_phase["Liq"],
            period=pyunits.hr,
        )
        # Add objective
        m.fs.obj = Objective(expr=m.fs.costing.LCOW)

        # Change values of flowsheet-level costing variables if desired
        m.fs.costing.electricity_cost.fix(0.12)
        m.fs.costing.base_currency = pyunits.USD_2023

        # Change values of unit model costing parameters if desired
        m.fs.costing.chem_addition.chemical_capex_slope.fix(2.34e4)
        m.fs.costing.bazchem.unit_cost.fix(0.42)
        m.fs.costing.reverse_osmosis.membrane_cost.fix(24)
        m.fs.costing.reverse_osmosis.factor_membrane_replacement.fix(0.02)


    m = real_build()
    real_add_costing(m)
    results = solver.solve(m)
    assert_optimal_termination(results)
    m.fs.del_component(m.fs.obj)


The :ref:`WaterTAP costing package<watertap_costing>` can be added to any WaterTAP flowsheet but is not required to run a WaterTAP model.
For this guide, we will consider a flowsheet with these units in the following order:

- Feed block
- Pump 1
- Chemical addition unit
- Pump 2 
- Reverse osmosis unit with custom chemical addition
- ERD unit
- Product block
- Brine block


Adding costing to a WaterTAP model
*******************************************

This section applies the steps outlined in the :ref:`how to add WaterTAP costing to a flowsheet guide<how_to_add_watertap_costing_to_flowsheet>` for all unit models on this flowsheet.

Adding the WaterTAP costing package to a flowsheet can be done at any point in the flowsheet build prior to adding unit model costing blocks. 
Below is a build function to create the flowsheet. The costing package is added after added the property package.


.. code-block:: python


    def build():

        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = SeawaterParameterBlock()

        # Add the WaterTAP costing package to the flowsheet
        m.fs.costing = WaterTAPCosting()

        # Add unit models to flowsheet
        m.fs.feed = Feed(property_package=m.fs.properties)
        m.fs.pump1 = Pump(property_package=m.fs.properties)

        m.fs.chem_addition = ChemicalAdditionZO(
            property_package=m.fs.properties, process_subtype="default"
        )

        # Create an expression for the mass flow of the chemical on the unit model block.
        # This will be used in our custom costing method.

        m.fs.chem_addition.chem_mass_flow = Expression(
            expr=pyunits.convert(
                m.fs.chem_addition.chemical_dosage[0]
                * m.fs.chem_addition.properties[0].flow_vol_phase["Liq"],
                to_units=pyunits.kg / pyunits.s,
            ),
            doc="Mass flow of chemical addition",
        )

        # Create a unit model that also uses a chemical that we want to track
        # for costing purposes.
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
        m.fs.RO.chem_dose = Param(
            initialize=690,
            units=pyunits.mg / pyunits.L,
            mutable=True,
            doc="Chemical dose for RO",
        )
        m.fs.RO.chem_mass_flow = Expression(
            expr=pyunits.convert(
                m.fs.RO.chem_dose
                * m.fs.RO.feed_side.properties[0, 0].flow_vol_phase["Liq"],
                to_units=pyunits.kg / pyunits.s,
            ),
            doc="Mass flow of chemical for RO",
        )

        m.fs.pump2 = Pump(property_package=m.fs.properties)
        m.fs.ERD = EnergyRecoveryDevice(property_package=m.fs.properties)
        m.fs.product = Product(property_package=m.fs.properties)
        m.fs.brine = Product(property_package=m.fs.properties)

        return m


Adding the WaterTAP costing package to the flowsheet is done by simply by creating an instance of the ``WaterTAPCosting`` class and assigning it to a flowsheet attribute. 
Convention is to name this attribute ``m.fs.costing``. 
This is referred to as the "flowsheet costing block" (contrasted with a "unit model costing block" discussed later). 
At this point, the flowsheet costing block only contains instructions for aggregating costs from individual unit model costing blocks into overall flowsheet-level costs and metrics.
To get costing results, we need to add costing blocks to each unit model.

Create custom costing method
************************************

All WaterTAP unit models in the `library <https://github.com/watertap-org/watertap/tree/main/watertap/costing/unit_models>`_ have a default costing method assigned. 
For the purposes of this guide, we assume use of a custom costing method for the chemical addition unit model.
See :ref:`how to create a custom costing method<how_to_create_custom_costing_method>` for instructions on how to create a custom costing method.


Add unit model costing
******************************

The example presented here uses the steps outlined in :ref:`how to add WaterTAP costing to a flowsheet <how_to_add_watertap_costing_to_flowsheet>` to add costing to all unit models on the flowsheet, including using the custom costing method for chemical addition.

At this point in our flowsheet build, we have our flowsheet costing block ``m.fs.costing``, our custom costing method for chemical addition defined, and default costing methods for the other unit models, but we have not yet added costing to any of our unit models.
Adding unit model costing is required because the flowsheet costing block will use these individual costing models to aggregate capital and operating costs at the system level.

To add costing to a unit model, we assign a ``UnitModelCostingBlock`` to the unit model's ``costing`` attribute and specify the flowsheet costing block as an argument.
For all the WaterTAP unit models on this flowsheet (chemical addition, RO, pumps, and ERD), that is all that is required because their custom costing methods are already assigned to their ``default_costing_method`` attribute. 
However, if we want to use our own custom costing method for chemical addition (i.e., bypass the default costing method), we can specify via the ``costing_method`` argument when adding the unit model costing block to ``m.fs.chem_addition``.

Costing proceeds in the following steps:

1. Assign a ``UnitModelCostingBlock`` to the ``costing`` attribute of each unit model, passing the flowsheet costing block as an argument. For the chemical addition unit, also pass the custom costing method via the ``costing_method`` argument.
2. After all unit model costing blocks have been added, call the ``cost_process()`` method on the flowsheet costing block to build the costing model. This method constructs the process-level costing components based on the registered unit operations and flows.
3. Add any desired system-level costing metrics (LCOW and SEC in this example) or flow breakdowns (the breakdown for bazchem across chemical addition and RO units).

The following snippet shows how to add costing to the unit models on our flowsheet, including using the custom costing method for chemical addition.

.. code-block:: python

    def add_costing(m):

        # Add unit model costing blocks 
        # The system costing block (m.fs.costing) is passed as the flowsheet costing block to each unit model costing block
        m.fs.pump1.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.RO.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.pump2.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.ERD.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

        # Pass the custom costing method to the chemical addition costing block via the costing_method argument
        m.fs.chem_addition.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=chem_addition_costing
        )
        # Cost the flow of bazchem for RO
        # Note this can only be done after registering the bazchem flow type
        m.fs.costing.cost_flow(m.fs.RO.chem_mass_flow, "bazchem")

        # Build the process costing model
        m.fs.costing.cost_process()

        # Designate flow rate for costing metrics
        flow_rate = m.fs.product.properties[0].flow_vol_phase["Liq"]

        # Add system-level costing metrics
        m.fs.costing.add_LCOW(flow_rate)
        m.fs.costing.add_specific_energy_consumption(flow_rate, name="SEC")

        # Add flow component breakdown for electricity and bazchem
        m.fs.costing.add_flow_component_breakdown(
            "electricity",
            flow_rate,
            period=pyunits.hr,
        )
        m.fs.costing.add_flow_component_breakdown(
            "bazchem",
            flow_rate,
            period=pyunits.hr,
        )

At this point in the build, you can change values of flowsheet-level costing variables (``electricity_price`` in this example) and change the ``base_currency`` 
or other costing variables/parameters. In the example below, we also add an ``Objective`` to minimize LCOW.

.. testcode::

    # Add objective
    m.fs.obj = Objective(expr=m.fs.costing.LCOW)

    # Change values of flowsheet-level costing variables if desired
    m.fs.costing.electricity_cost.fix(0.12)
    m.fs.costing.base_currency = pyunits.USD_2023

    # Change values of unit model costing parameters if desired
    m.fs.costing.chem_addition.chemical_capex_slope.fix(2.34e4)
    m.fs.costing.bazchem.unit_cost.fix(0.42)
    m.fs.costing.reverse_osmosis.membrane_cost.fix(24)
    m.fs.costing.reverse_osmosis.factor_membrane_replacement.fix(0.02)


Access costing results
******************************

The approaches for accessing results are also presented in the :ref:`how to access costing results guide<how_to_access_costing_results>`.

After building and solving the model, system-level costing results are on the flowsheet costing block (``m.fs.costing``) and unit-level costing results are on the ``costing`` attribute of each unit model (e.g. ``m.fs.chem_addition.costing``).
Accessing the values for any variable, expression, parameter, etc. can be done using the ``value`` function from Pyomo *or* by "calling" the component directly (i.e., placing ``()`` after the component name). 
Examples of both of these approaches are presented below.


Unit model costing results
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Costing results for individual unit models are on the unit model costing block. These are the individual ``UnitModelCostingBlock`` that we assigned to the ``costing`` attribute of each 
unit model in the ``add_costing`` function presented in this guide. So, to access the capital cost for each unit model:

.. testcode::

    from pyomo.environ import value

    # Example of accessing results with value()
    pump1_capital_cost = value(m.fs.pump1.costing.capital_cost)
    chem_addition_capital_cost = value(m.fs.chem_addition.costing.capital_cost)
    pump2_capital_cost = value(m.fs.pump2.costing.capital_cost)
    
    # Example of accessing results by directly calling the component
    ro_capital_cost = m.fs.RO.costing.capital_cost()
    erd_capital_cost = m.fs.ERD.costing.capital_cost()
    chem_addition_capital_cost = m.fs.chem_addition.costing.capital_cost()

Likewise, if the unit model has any fixed operating costs, they can be accessed in a similar manner:

.. testcode::

    ro_fixed_operating_cost = value(m.fs.RO.costing.fixed_operating_cost)
    chem_addition_fixed_operating_cost = value(m.fs.chem_addition.costing.fixed_operating_cost)

Importantly, the currency units for these values are only converted to the ``base_currency`` and ``base_period`` defined on the flowsheet costing block if the model developer 
did so in the costing method. For this reason, it is recommended to make this conversion when creating costing methods to ensure consistency in the reported values.
If you are unsure, the units for costing variables (or any variable) can be accessed using ``pyunits.get_units(var)`` from the units package and printing the output.

.. testcode::

    ro_capex_units = pyunits.get_units(m.fs.RO.costing.capital_cost)
    chem_addition_capex_units = pyunits.get_units(m.fs.chem_addition.costing.capital_cost)
    pump_capex_units = pyunits.get_units(m.fs.pump1.costing.capital_cost)
    erd_capex_units = pyunits.get_units(m.fs.ERD.costing.capital_cost)

    print(f"Base currency: {m.fs.costing.base_currency}")
    print(f"RO capital cost units: {ro_capex_units}")
    print(f"Chemical addition capital cost units: {chem_addition_capex_units}")
    print(f"Pump capital cost units: {pump_capex_units}")
    print(f"ERD capital cost units: {erd_capex_units}")

.. testoutput::

    Base currency: USD_2023
    RO capital cost units: USD_2018
    Chemical addition capital cost units: USD_2018
    Pump capital cost units: USD_2018
    ERD capital cost units: USD_2018


System costing results
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Costing results for all units with a costing block are aggregated to the system level and converted to be in the same currency year (as defined by ``base_currency``). Operating costs are likewise converted to 
units of ``base_currency / base_period``. These results are found on the flowsheet costing block (``m.fs.costing``) and for this example includes the following:

- Total capital cost: ``m.fs.costing.total_capital_cost``
- Total operating cost: ``m.fs.costing.total_operating_cost``
- Total electricity required: ``m.fs.costing.aggregate_flow_electricity``
- Total cost of electricity: ``m.fs.costing.aggregate_flow_costs["electricity"]``
- Total cost of bazchem: ``m.fs.costing.aggregate_flow_costs["bazchem"]``
- LCOW: ``m.fs.costing.LCOW``
- SEC: ``m.fs.costing.SEC``


.. testcode::

    total_capital_cost = value(m.fs.costing.total_capital_cost)
    total_operating_cost = value(m.fs.costing.total_operating_cost)
    LCOW = m.fs.costing.LCOW()
    SEC = m.fs.costing.SEC()


In this example, we have registered flows of ``"electricity"`` and ``"bazchem"`` for costing purposes. Results for all registered flows are accessed through the ``aggregate_flow_*`` variable(s) on the costing block, 
where the ``*`` is replaced by the name of the registered flow (e.g., ``electricity`` or ``bazchem``).

.. testcode::

    # Example of accessing aggregate flow results
    electricity_flow = value(m.fs.costing.aggregate_flow_electricity) # kW
    bazchem_flow = value(m.fs.costing.aggregate_flow_bazchem) # kg/s

Costs for these flows are accessed via the ``aggregate_flow_costs`` indexed variable on the costing block. The index is the registered name of the flow.

.. testcode::

    # Example of accessing aggregate flow costs
    electricity_cost = value(m.fs.costing.aggregate_flow_costs["electricity"]) # $/year
    bazchem_cost = m.fs.costing.aggregate_flow_costs["bazchem"]() # $/year


:ref:`LCOW<aggregate_metric_LCOW>` and :ref:`SEC<aggregate_metric_SEC>` are further broken down into the contributions from each unit model, unit model class, and flow type in the following expressions.

- Direct capital expenditures by flowsheet component: ``m.fs.costing.LCOW_component_direct_capex``
- Indirect capital expenditures by flowsheet component: ``m.fs.costing.LCOW_component_indirect_capex``
- Fixed operating expenditures by flowsheet component: ``m.fs.costing.LCOW_component_fixed_opex``
- Variable operating expenditures by flowsheet component: ``m.fs.costing.LCOW_component_variable_opex``
- Aggregate direct capital expenditures by unit model type: ``m.fs.costing.LCOW_aggregate_direct_capex``
- Aggregate indirect capital expenditures by unit model type: ``m.fs.costing.LCOW_aggregate_indirect_capex``
- Aggregate fixed operating expenditures by unit model type: ``m.fs.costing.LCOW_aggregate_fixed_opex``
- Aggregate variable operating expenditures by unit model type: ``m.fs.costing.LCOW_aggregate_variable_opex``
- Specific energy consumption by unit model: ``m.fs.costing.SEC_component``


.. testcode::

    # Example of accessing breakdown of LCOW results
    pump1_direct_capex = value(m.fs.costing.LCOW_component_direct_capex["fs.pump1"])
    pump2_direct_capex = value(m.fs.costing.LCOW_component_direct_capex["fs.pump2"])
    electricity_variable_opex = m.fs.costing.LCOW_aggregate_variable_opex["electricity"]()
    chemical_variable_opex = m.fs.costing.LCOW_aggregate_variable_opex["bazchem"]()

    # Example of accessing breakdown of SEC results
    pump1_SEC = value(m.fs.costing.SEC_component["fs.pump1"])
    pump2_SEC = value(m.fs.costing.SEC_component["fs.pump2"])

Below is an example of how to display the LCOW and SEC breakdowns using the ``display()`` method on various expressions.

.. code-block:: python

    print("Example LCOW breakdowns")
    m.fs.costing.LCOW.display()
    m.fs.costing.aggregate_flow_costs.display()
    m.fs.costing.LCOW_component_direct_capex.display()
    m.fs.costing.LCOW_component_variable_opex.display()
    m.fs.costing.LCOW_aggregate_variable_opex.display()

    print("\nExample SEC breakdowns")
    m.fs.costing.SEC.display()
    m.fs.costing.SEC_component.display()

    print("\nExample flow breakdowns")
    m.fs.costing.electricity_component.display()
    m.fs.costing.bazchem_component.display()

The example output below will differ depending on the specific model and parameters used.

The component-based LCOW breakdowns are indexed by each unit model on the flowsheet. However,
the aggregation-based LCOW breakdowns are indexed by each unit model type and flow type. 

Each registered flow appears in multiple places within the LCOW breakdowns. For example, for the ``LCOW_component_variable_opex`` breakdown, the electricity flow 
is included in the results for the `fs.pump1`, `fs.pump2`, and `fs.ERD` indexes; the sum of these values (i.e., the total aggregation of electricity) equals the ``LCOW_aggregate_variable_opex["electricity"]`` entry.
Additionally, ``LCOW_aggregate_variable_opex["Pump"]`` corresponds to the sum of the variable costs for all pump unit models. 
For this reason, the system LCOW is the summation of all indexes in any of the component or aggregate expressions *except* those indexed by a material/energy flow (see :ref:`how to cost a flow<how_to_cost_a_flow>`).

.. code-block:: none

    Example LCOW breakdowns
    LCOW : Size=1
        Key  : Value
        None : 1.14016307973748
    aggregate_flow_costs : Size=2, Index=fs.costing.used_flows, Units=USD_2023/a
        Key         : Lower : Value             : Upper : Fixed : Stale : Domain
            bazchem :  None : 132.5419200224229 :  None : False : False :  Reals
        electricity :  None : 7472.211129411062 :  None : False : False :  Reals
    LCOW_component_direct_capex : Size=5
        Key              : Value
                fs.pump1 : 0.019967960566987127
                fs.pump2 :  0.13331600455723355
                fs.RO :  0.01454986992197548
                fs.ERD : 0.008982485001933636
        fs.chem_addition :  0.04963718498491835
    LCOW_component_variable_opex : Size=5
        Key              : Value
                fs.pump1 :  0.09907875980955234
                fs.pump2 :   0.6614989222351185
                fs.RO :                  0.0
                fs.ERD : -0.28701751725730146
        fs.chem_addition : 0.008399999999999998
    LCOW_aggregate_variable_opex : Size=6
        Key                  : Value
                        Pump :   0.7605776820446709
            ReverseOsmosis1D :                  0.0
        EnergyRecoveryDevice : -0.28701751725730146
        ChemicalAdditionZO : 0.008399999999999998
                electricity :   0.4735601647873693
                    bazchem : 0.008399999999999998
    
    Example SEC breakdowns
    SEC : Size=1
        Key  : Value
        None : 2.982873118845955
    SEC_component : Size=3
        Key      : Value
        fs.pump1 :  0.6240798767717449
        fs.pump2 :   4.166666666666667
        fs.ERD : -1.8078734245924564

    Example flow breakdowns
    electricity_component : Size=3
        Key      : Value
        fs.pump1 :  0.6240798818300959
        fs.pump2 :   4.166666666666666
        fs.ERD : -1.8078734308552329
    bazchem_component : Size=2
        Key              : Value
        fs.chem_addition :  5.555555555555555e-06
                fs.RO : 0.00038333333333333324