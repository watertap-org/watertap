.. _electricity_methods:

Helper Methods for Electricity Demand
=====================================

.. currentmodule:: watertap.core.zero_order_electricity

A number of methods are available for implementing common forms for calculating electricity intensity in zero-order type models. These helper functions are intended to be called after all necessary material and chemical flows have been added to the unit model.

Usage
-----

.. doctest::

  from idaes.core import declare_process_block_class
  from watertap.core import build_pt, constant_intensity, ZeroOrderBaseData

  @declare_process_block_class("PumpZO")
  class PumpZOData(ZeroOrderBaseData):

      CONFIG = ZeroOrderBaseData.CONFIG()

      def build(self):
          super().build()

          self._tech_type = "pump"

          build_pt(self)

          # Add variables and constraints for constant electricity intensity
          constant_intensity(self)

Module Documentation
--------------------

* :mod:`watertap.core.zero_order_electricity`
