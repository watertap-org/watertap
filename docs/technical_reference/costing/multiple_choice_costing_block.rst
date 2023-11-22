.. _multiple_choice_costing_block:

Multiple Choice Unit Model Costing Block
========================================

.. currentmodule:: watertap.costing.multiple_choice_costing_block

The MultiUnitModelCostingBlock class is a wrapper around the IDAES
`UnitModelCostingBlock` which allows the modeller to easily explore the implications of
different unit model costing choices on overall flowsheet costs without rebuilding
a new flowsheet or replacing Pyomo components.

MultiUnitModelCostingBlock objects instead construct every single costing relationship
specified for the unit model, and controls which one is active through the indexed
mutable parameter `costing_block_selector`, which takes the value `1` when the
associated costing block is active and the value `0` otherwise. The helper method
`select_costing_block` handles activating a single costing block whilst deactivating
all other costing blocks.

Usage
-----

The MultiUnitModelCostingBlock allows the user to pre-specify all costing methods
on-the-fly and change between them without rebuilding the flowsheet or even the base
costing relationships on the costing package.

The code below demonstrates its use on a reverse osmosis unit model.

.. testcode::

    import pyomo.environ as pyo
    from idaes.core import FlowsheetBlock

    import watertap.property_models.NaCl_prop_pack as props
    from watertap.unit_models.reverse_osmosis_0D import (
        ReverseOsmosis0D,
        ConcentrationPolarizationType,
        MassTransferCoefficient,
        PressureChangeType,
    )
    from watertap.costing import WaterTAPCosting, MultiUnitModelCostingBlock
    from watertap.costing.unit_models.reverse_osmosis import (
        cost_reverse_osmosis,
    )

    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.NaClParameterBlock()
    m.fs.costing = WaterTAPCosting()

    m.fs.RO = ReverseOsmosis0D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
    )

    def my_reverse_osmosis_costing(blk):
        blk.variable_operating_cost = pyo.Var(
            initialize=42,
            units=blk.costing_package.base_currency / blk.costing_package.base_period,
            doc="Unit variable operating cost",
        )
        blk.variable_operating_cost_constraint = pyo.Constraint(
            expr=blk.variable_operating_cost
            == 42 * blk.costing_package.base_currency / blk.costing_package.base_period
        )

    m.fs.RO.costing = MultiUnitModelCostingBlock(
        # always needed, just like for UnitModelCostingBlock
        flowsheet_costing_block=m.fs.costing,

        # The keys to the costing-block are a user-defined name.
        # The values are either the costing method itself or another
        # dictionary with the keys "costing_method", specifying the
        # costing method, and optionally "costing_method_arguments",
        # which defines the keyword arguments into the costing method
        costing_blocks={
            "normal_pressure": cost_reverse_osmosis,
            "high_pressure": {
                "costing_method": cost_reverse_osmosis,
                "costing_method_arguments": {"ro_type": "high_pressure"},
            },
            "custom_costing": my_reverse_osmosis_costing,
        },

        # This argument is optional, if it is not specified the block
        # utilizes the first method, in this case "normal_pressure"
        initial_costing_block="normal_pressure",
    )

    # set an area
    m.fs.RO.area.set_value(100)

    # create the system-level aggregates
    m.fs.costing.cost_process()

    # initialize the unit level costing blocks and the flowsheet costing block
    m.fs.costing.initialize()

    # Since the `initial_costing_block` was active, its aggregates are used:
    assert ( pyo.value(m.fs.RO.costing.capital_cost) == 
        pyo.value(m.fs.RO.costing.costing_blocks["normal_pressure"].capital_cost)
    )
    assert ( pyo.value(m.fs.costing.aggregate_capital_cost) == 
        pyo.value(m.fs.RO.costing.costing_blocks["normal_pressure"].capital_cost)
    )
    assert pyo.value(m.fs.RO.costing.variable_operating_cost) == 0
    assert pyo.value(m.fs.costing.aggregate_variable_operating_cost) == 0

    # We can activate the "high_pressure" costing block:
    m.fs.RO.costing.select_costing_block("high_pressure")

    # Need re-initialize to have the new values in the aggregates
    m.fs.costing.initialize()

    assert ( pyo.value(m.fs.RO.costing.capital_cost) == 
        pyo.value(m.fs.RO.costing.costing_blocks["high_pressure"].capital_cost)
    )
    assert ( pyo.value(m.fs.costing.aggregate_capital_cost) == 
        pyo.value(m.fs.RO.costing.costing_blocks["high_pressure"].capital_cost)
    )
    assert pyo.value(m.fs.RO.costing.variable_operating_cost) == 0
    assert pyo.value(m.fs.costing.aggregate_variable_operating_cost) == 0

    # We can activate the "custom_costing" costing block:
    m.fs.RO.costing.select_costing_block("custom_costing")

    # Need re-initialize to have the new values in the aggregates
    m.fs.costing.initialize()

    # No capital cost for block "custom_costing", but it does have
    # a "variable" operating cost 
    assert pyo.value(m.fs.RO.costing.capital_cost) == 0
    assert pyo.value(m.fs.costing.aggregate_capital_cost) == 0
    assert pyo.value(m.fs.RO.costing.variable_operating_cost) == 42
    assert pyo.value(m.fs.costing.aggregate_variable_operating_cost) == 42


Class Documentation
-------------------

* :class:`MultiUnitModelCostingBlock`
