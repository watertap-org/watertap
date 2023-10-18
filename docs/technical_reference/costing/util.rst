.. _costing_utils:

Utility Functions for Costing
=============================

.. module:: watertap.costing.util
   :noindex:

The module :mod:`watertap.costing.util` contains several utility functions for
building unit model costing methods.

Utility for Defining Global-Level Parameters Specific to a Unit Model
---------------------------------------------------------------------

The ``register_costing_parameter_block()`` is a decorator for unit model costing
functions. It works by automatically inserting a specified parameter block onto
the costing package. It includes additional safeguards for ensuring a costing parameter
block is not overwritten during model creation. Below is an example of the
functionality:

.. testcode::

    import pyomo.environ as pyo
    import idaes.core as idc
    from watertap.costing import (
        WaterTAPCosting,
        register_costing_parameter_block,
        make_capital_cost_var,
    )

    @idc.declare_process_block_class("MyUnitModel")
    class MyUnitModelData(idc.UnitModelBlockData):

        def build(self):
            super().build()
            self.foo_flow = pyo.Var(
                initialize=10,
                units=pyo.units.kg / pyo.units.s,
                doc="foo",
            )

        @property
        def default_costing_method(self):
            # could point to a static method on
            # this class, could be function in 
            # a different module even
            return cost_my_unit_model

    def build_my_unit_model_param_block(blk):
        """
        This function builds the global parameters and custom flows for MyUnitModel.
        """
        blk.fixed_capital_cost = pyo.Var(
            initialize=42,
            doc="Fixed capital cost for all of my units",
            units=pyo.units.USD_2020,
        )

        blk.foo_flow_cost = pyo.Var(
            initialize=3,
            doc="Foo flow cost",
            units=pyo.units.USD_2016 / pyo.units.kg,
        )
        blk.parent_block().register_flow_type("foo", blk.foo_flow_cost)

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

        # here we reference the `fixed_capital_cost` parameter
        # automatically added by the `register_costing_parameter_block`
        # decorator.
        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.costing_package.my_unit.fixed_capital_cost,
            name="fixed capital cost constraint",
        )

        blk.costing_package.cost_flow(blk.unit_model.foo_flow, "foo")

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

    # Here everythin as before, but the global parameter block
    # m.fs.costing.my_unit is not re-built. 
    m.fs.my_unit_2.costing = idc.UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
    )

.. autofunction:: register_costing_parameter_block 
   :noindex:


Utilities for Common Variable Creation
--------------------------------------

The IDAES costing framework utilizes specific names on a unit model's costing block to
calculate aggregates like `aggregate_capital_cost` and `aggregate_fixed_operating_cost`.
The ``make_capital_cost_var()`` and ``make_fixed_operating_cost_var()`` are utilities
to create a Pyomo ``Var`` with these standard names.

.. autofunction:: make_capital_cost_var
   :noindex:

.. autofunction:: make_fixed_operating_cost_var
   :noindex:


Utilities for Common Cost Calculations
--------------------------------------

Some costing functions/methods share similar structure. For example, costing for
nanofiltration and reverse osmosis is largely a function of the membrane total area.
Therefore, WaterTAP has a few utility functions to capture these similar features.

.. autofunction:: cost_by_flow_volume
   :noindex:

.. autofunction:: cost_membrane
   :noindex:

.. autofunction:: cost_rectifier
   :noindex:
