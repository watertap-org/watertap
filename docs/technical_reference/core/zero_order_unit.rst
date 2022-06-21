.. _ZeroOrderBaseData:

Zero Order Unit Model Base Class
================================

.. index::
   pair: watertap.core.zero_order_base;ZeroOrderBaseData

.. currentmodule:: watertap.core.zero_order_base

The zero-order unit model base class is intended to form the basis for a library of generic zero-order type (fixed performance) unit models. It defines a standard API for these types of unit operations and a set of common helper methods for building and working with these.

Usage
-----

The ZeroOrderBaseData class contains no variables or constraints of its own - rather it is intended to be used as the starting point for construction of zero-order type models for specific unit operations.

The code below shows an outline of how the ZeroOrderBaseData class is intended to be used to develop custom zero-order type models.

.. doctest::

  from idaes.core import declare_process_block_class
  from watertap.core import ZeroOrderBaseData

  @declare_process_block_class("MyZOUnit")
  class MyZOUnitData(ZeroOrderBaseData):

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "my_technology"

        # Add necessary StateBlocks, Ports, Variables and Constraints
        # A library of methods for common model form is available,
        # or users can custom code their own components

Property Package Requirements
-----------------------------

The ZeroOrderBase class makes a number of assumptions about the structure of the associated property package, and contains a number of checks to ensure the property package meets these requirements. The requirements for the property package are:

1. A single phase named `Liq`.
2. A single `Solvent` species named `H2O` (package must define a `solvent_set` with length one that contains "H2O").
3. All other species are defined as `Solutes` (package must define a `solute_set` and the `component_list` must be the union of `solvent_set` and `solute_set`).
4. The property package must define `flow_vol`, `conc_mass_comp`, `pressure` and `temperature` as properties.

Class Documentation
-------------------

* :class:`ZeroOrderBaseData`
