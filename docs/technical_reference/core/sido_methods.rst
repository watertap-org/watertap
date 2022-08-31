.. _sido_methods:

Single Inlet - Double Outlet Helper Methods
===========================================

.. index::
   pair: watertap.core.zero_order_sido;build_sido

.. currentmodule:: watertap.core.zero_order_sido

The `build_sido` method is intended to be used to rapidly construct a standard set of material balance equations for zero-order type models with a single inlet and two outlets.

Usage
-----

.. doctest::

  from idaes.core import declare_process_block_class
  from watertap.core import build_sido, ZeroOrderBaseData

  @declare_process_block_class("NanofiltrationZO")
  class NanofiltrationZOData(ZeroOrderBaseData):

      CONFIG = ZeroOrderBaseData.CONFIG()

      def build(self):
          super().build()

          self._tech_type = "nanofiltration"

          build_sido(self)

Model Structure
---------------

The build_sido method constructs a simple representation of unit operation with a single inlet (named `inlet`) and two outlets (named `treated` and `byproduct`). A `StateBlock` is constructed for the inlet and each outlet with a `Port` associated with each of these.

Variables
---------

The build_sido method creates the following variables in addition to those created by the `StateBlocks`.

=============================== ========================= =============== ======================================================================
Variable                        Name                      Indices         Notes
=============================== ========================= =============== ======================================================================
:math:`r_{t}`                   recovery_frac_mass_H2O    time            Fraction of mass flow of water in inlet that goes to treated stream.
:math:`f_{t,j}`                 removal_frac_mass_comp    time, component Fraction of mass flow of each component that goes to byproduct stream.
=============================== ========================= =============== ======================================================================

Constraints
-----------

The build_sido method writes the following constraints which relate the inlet state to those in the treated and byproduct streams.
First, a mass recovery equation for water is written to relate the flowrate at the treated outlet to that at the inlet:

`water_recovery_equation(t)`:

.. math:: r_t \times M_{inlet,t,H2O} = M_{treated,t,H2O}

where :math:`M_{t,H2O}` is mass flowrate of water at time :math:`t`.

Next, a mass balance for water is written to relate the flowrate in the byproduct stream to that in the inlet and treated streams.

`flow_balance(t)`:

.. math:: M_{inlet,t,H2O} = M_{treated,t,H2O} + M_{byproduct,t,H2O}

Removal constraints are then written for each solute to relate the solute concentration in the byproduct stream to that in the inlet stream.

`solute_removal_equation(t, j)`:

.. math:: f_{t, j} \times M_{inlet,t,j} = M_{byproduct,t,j}

A mass balance constraint is then written for each solute.

`solute_treated_equation(t, j)`:

.. math:: M_{inlet,t,j} = M_{treated,t,j} + M_{byproduct,t,j}

Module Documentation
--------------------

* :mod:`watertap.core.zero_order_sido`
