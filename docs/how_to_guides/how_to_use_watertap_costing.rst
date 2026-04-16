.. _how_to_use_watertap_costing:

How to use WaterTAP Costing
===========================

Overview
--------

This guide shows you how to add the WaterTAP costing package, add unit model costing, access costing results, and create a 
custom costing model. 

How To
------

.. testsetup::

    # quiet idaes logs
    import idaes.logger as idaeslogger
    idaeslogger.getLogger('ideas.core').setLevel('CRITICAL')
    idaeslogger.getLogger('ideas.core.util.scaling').setLevel('CRITICAL')
    idaeslogger.getLogger('idaes.init').setLevel('CRITICAL')


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


Below is the code to build this flowsheet. 


.. code:: 


    def build():

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
                m.fs.chem_addition.chem_dose * m.fs.feed.properties[0].flow_vol_phase["Liq"],
                to_units=pyunits.kg / pyunits.s,
            ),
            doc="Mass flow of chemical addition",
        )
        m.fs.chem_addition.default_costing_method = chem_addition_costing

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


Note that the chemical addition unit is a state junction with a dose variable added and expression to calculate the mass flow of the chemical
Unlike other unit models on this flowsheet, this chemical addition unit model does not have a built in costing method, so we will create and assign it a custom costing method.
Above you can see we assigned the custom costing method ``chem_addition_costing`` to the chemical addition unit model property ``default_costing_method``. 
The code below shows how to build the custom costing method.

.. code::

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


.. code::


    def add_costing(m):
        m.fs.costing = WaterTAPCosting()
        m.fs.costing.base_currency = pyunits.USD_2023

        m.fs.pump1.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.pump2.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.RO.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.ERD.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

        m.fs.chem_addition.default_costing_method = chem_addition_costing
        m.fs.chem_addition.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing
        )

        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol_phase["Liq"])
        m.fs.costing.add_specific_energy_consumption(
            m.fs.product.properties[0].flow_vol_phase["Liq"], name="SEC"
        )
        m.fs.obj = Objective(expr=m.fs.costing.LCOW)