.. _how_to_access_costing_results:

How to access costing results
===============================

This guide provides instructions on how to access costing results from a flowsheet with WaterTAP costing.
See :ref:`how to add the costing packages<how_to_add_watertap_costing_to_flowsheet>` for instructions on adding costing to a model.

After building and solving a flowsheet, there are results associated with the flowsheet costing block and each individual unit model costing block.
In this example, we assume the costing models are structured as follows:

.. code-block:: python

    m.fs.costing = WaterTAPCosting()
    m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

We call ``m.fs.costing`` the "flowsheet costing block" and ``m.fs.unit.costing`` the "unit model costing block".

Results for only the unit model are found on the unit model costing block and could include:

- Capital cost: ``m.fs.unit.costing.capital_cost``
- Operating cost: ``m.fs.unit.costing.fixed_operating_cost``

The specific results available will depend on if they are included in the costing method for the unit model.

Results from all the unit model costing blocks are aggregated to the flowsheet costing block and could include:

- Total capital cost: ``m.fs.costing.total_capital_cost``
- Total operating cost: ``m.fs.costing.total_operating_cost``
- Total electricity required: ``m.fs.costing.aggregate_flow_electricity``
- Total cost of electricity: ``m.fs.costing.aggregate_flow_costs["electricity"]``
- LCOW: ``m.fs.costing.LCOW``
- SEC: ``m.fs.costing.SEC``


Accessing the values for any variable, expression, parameter, etc. can be done using the ``value`` function from Pyomo *or* by "calling" the component directly (i.e., placing ``()`` after the component name).
In either case, if the component is indexed, the index(es) must be specified to access the value for a specific index.

.. code-block:: python

    from pyomo.environ import value

    # Accessing costing results using the value function
    unit_capital_cost = value(m.fs.unit.costing.capital_cost)
    total_capital_cost = value(m.fs.costing.total_capital_cost)
    total_operating_cost = value(m.fs.costing.total_operating_cost)
    total_electricity = value(m.fs.costing.aggregate_flow_electricity)

    # Accessing costing results by "calling" the component directly
    unit_operating_cost = m.fs.unit.costing.fixed_operating_cost()
    # Access costs of an indexed component by specifying the index
    total_electricity_cost = m.fs.costing.aggregate_flow_costs["electricity"]()
    LCOW = m.fs.costing.LCOW()
    SEC = m.fs.costing.SEC()

Alternatively, users can use the ``display`` method on these components to view their values directly. In addition to the value of the component, 
the ``display`` method will also show additional information such as units, bounds, and other metadata associated with the modeling component.

Below is an example of how to display costing results using the ``display()`` method.

.. code-block:: python

    m.fs.costing.LCOW.display()
    m.fs.costing.total_capital_cost.display()
    m.fs.costing.aggregate_flow_costs.display()