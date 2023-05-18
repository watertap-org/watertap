.. _pt_methods:

Pass Through Unit Helper Methods
================================

.. index::
   pair: watertap.core.zero_order_pt;build_pt

.. currentmodule:: watertap.core.zero_order_pt

The :func:`build_pt` method is intended to be used to rapidly construct a standard set of material balance equations for zero-order type models in which there is no change in flow rate or concentration of the solutes. This is most commonly seen in pressure and temperature change equipment.

Usage
-----

.. doctest::

  from idaes.core import declare_process_block_class
  from watertap.core import build_pt, ZeroOrderBaseData

  @declare_process_block_class("PumpZO")
  class PumpZOData(ZeroOrderBaseData):

      CONFIG = ZeroOrderBaseData.CONFIG()

      def build(self):
          super().build()

          self._tech_type = "pump"

          build_pt(self)

Model Structure
---------------

The build_pt method constructs a simple representation of unit operation with a single inlet (named `inlet`) and outlet (named `outlet`) which have the same material state. A single `StateBlock` is constructed which is used for both the inlet and outlet Ports.

Variables
---------

The build_pt method creates no additional variables beyond those created by the `StateBlocks`.

Constraints
-----------

The build_pt method writes no additional constraints beyond those created by the `StateBlocks`.

Module Documentation
--------------------

* :mod:`watertap.core.zero_order_pt`
