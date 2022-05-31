.. _zero_order_costing:

Zero Order Costing Package
==========================

.. index::
   pair: watertap.core.zero_order_costing;ZeroOrderCostingData

.. currentmodule:: watertap.core.zero_order_costing

The zero order costing module contains the costing package used for zero order models. Technoeconomic data used for zero order models is contained in the 
``.yaml`` file for that model located in the data/techno_economic folder.


Usage
-----

The ZeroOrderCosting class contains all the variables and constraints needed to cost a unit model derived from the :ref:`ZeroOrderBaseData`.

The code below shows an outline of how the ZeroOrderCostingData class is intended to be used to cost zero-order type models.

.. testcode::

  from pyomo.environ import ConcreteModel
  
  from idaes.core import FlowsheetBlock

  from watertap.core.zero_order_costing import ZeroOrderCosting
  from watertap.core.wt_database import Database
  from watertap.core.zero_order_properties import WaterParameterBlock
  from watertap.unit_models.zero_order import MyZOUnit


  m = ConcreteModel()
  m.db = Database()
  m.fs = FlowsheetBlock(default{"dynamic": False})
  m.fs.params = WaterParameterBlock(default={"solute_list": ["comp_a", "comp_b", "comp_c"]})
  m.fs.costing = ZeroOrderCosting()
  m.fs.unit = MyZOUnit(default={"property_package": m.fs.params, "database": m.db)

  # Add necessary statements to fix component flows prior to solve

.. autoclass:: ZeroOrderCostingData
    :members:
