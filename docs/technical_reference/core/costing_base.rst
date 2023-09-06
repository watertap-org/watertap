.. _costing_base:

WaterTAP Costing Block
======================

.. currentmodule:: watertap.core.costing_base


Usage
-----

The WaterTAPCostingBlock class extends the functionality of the `IDAES Process Costing Framework <https://idaes-pse.readthedocs.io/en/stable/reference_guides/core/costing/costing_framework.html>`_ in two ways:

#. Unit models can self-register a default costing method by specifying a `default_costing_method` attribute. This allows the costing method(s) to be specified with the unit model definition.


#. The method `register_flow_type` will create a new Expression if a costing component is not already defined *and* the costing component is not constant. The default behavior in IDAES is to always create a new Var. This allows the user to specify intermediate values in `register_flow_type`. 

.. testcode::

    import pyomo.environ as pyo
    from idaes.core import FlowsheetBlock
    from watertap.costing import WaterTAPCosting

    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock()
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


