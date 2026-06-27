.. _how_to_cost_a_flow:

How to cost a flow in WaterTAP costing
=======================================

.. testsetup:: python

    from pyomo.environ import (
        ConcreteModel,
        Param,
        Var,
        Expression,
        assert_optimal_termination,
    )
    from pyomo.environ import units as pyunits

    from idaes.core import FlowsheetBlock

    from watertap.costing import WaterTAPCosting
    from watertap.core.solvers import get_solver

    solver = get_solver()

    # quiet idaes logs
    import idaes.logger as idaeslogger
    idaeslogger.getLogger('idaes.core').setLevel('CRITICAL')
    idaeslogger.getLogger('idaes.core.util.scaling').setLevel('CRITICAL')
    idaeslogger.getLogger('idaes.init').setLevel('CRITICAL')

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = WaterTAPCosting()

    bazchem_cost = 0.42 * pyunits.USD_2020 / pyunits.kg
    m.fs.costing.register_flow_type("bazchem", bazchem_cost)
    assert "bazchem" in m.fs.costing.flow_types

    # Test different components for costing flows
    m.fs.bazchem_flow_mass1 = Var(
        initialize=1,
        units=pyunits.kg / pyunits.s,
        bounds=(0, None),
        doc="Mass flow of bazchem",
    )
    m.fs.bazchem_flow_mass1.fix()
    m.fs.bazchem_flow_mass2 = Expression(expr=2 * pyunits.kg / pyunits.s)
    m.fs.bazchem_flow_mass3 = Param(initialize=3, units=pyunits.kg / pyunits.s)
    m.fs.bazchem_flow_mass4 = 4 * pyunits.kg / pyunits.s

    m.fs.costing.cost_flow(m.fs.bazchem_flow_mass1, "bazchem")
    m.fs.costing.cost_flow(m.fs.bazchem_flow_mass2, "bazchem")
    m.fs.costing.cost_flow(m.fs.bazchem_flow_mass3, "bazchem")
    m.fs.costing.cost_flow(m.fs.bazchem_flow_mass4, "bazchem")


    m.fs.costing.cost_process()

    m.fs.costing.initialize()
    results = solver.solve(m)
    assert_optimal_termination(results)



In the WaterTAP costing package, variable operational costs for the system are calculated by collecting all the material and energy flows from each unit model on the flowsheet. 
These commonly include flows of power (electricity) and various chemicals.
If you are using a unit model with a default costing method, the necessary flows are already costed in that unit model costing method.
However, if you are using a new costing model that has a material or energy flow, you will need to add that to the flowsheet.

In this guide, we demonstrate how to register and cost a flow at the flowsheet level. 
Frequently these steps are done within a unit model costing method; 
see :ref:`how to create a custom costing method<how_to_create_custom_costing_method>` for an example of how to cost this flow within a unit model costing method.

Costing a flow in WaterTAP has two steps:

1. Register the flow type with the flowsheet costing block.
2. Cost the flow expression in a costing method.

In this example, we will register and cost a flow of a new chemical "bazchem" and assume there is a flowsheet costing block called ``m.fs.costing``.

The first step is to register the flow type.
Many flows, like ``"electricity"``, are already registered by default.
To register the flow type, use the ``register_flow_type`` method of the flowsheet costing block.
This method takes two arguments: the name of the flow type and the cost of the flow.

.. code-block:: python

    bazchem_cost = 0.42 * pyunits.USD_2020 / pyunits.kg

    m.fs.costing.register_flow_type("bazchem", bazchem_cost)

Here, we have registered a flow type called "bazchem" with a cost of 0.42 $/kg.
It is important to include the units in the cost expression to ensure that the costing calculations are done with the correct units.
This method will create a variable on the flowsheet costing block ``m.fs.costing.bazchem_cost`` that will be used to calculate "bazchem" flow costs, 
and adds ``"bazchem"`` to the ``flow_types`` set on the flowsheet costing block.

Next, we need to cost the flow expression. 
Assume that we have a mass-based flow of "bazchem" in our unit model that we want a cost for called ``m.fs.unit.bazchem_flow_mass`` with units of kg/min.
To cost this flow, we use the ``cost_flow`` method of the flowsheet costing block.
This method takes two arguments: the modeling component that represents the flow and the name of the flow type (which we just registered as "bazchem" in step 1).

.. code-block:: python

    m.fs.costing.cost_flow(m.fs.unit.bazchem_flow_mass, "bazchem")

The first argument can be anything that represents the flow (variable, parameter, expression, etc.) as long as the product of the flow expression and flow cost 
are in units of currency per time (e.g., $/year).
