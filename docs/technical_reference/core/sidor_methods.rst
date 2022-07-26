.. _sidor_methods:

Reactive Single Inlet - Double Outlet Helper Methods
====================================================

.. index::
   pair: watertap.core.zero_order_sido_reactive;build_sido_reactive

.. currentmodule:: watertap.core.zero_order_sido_reactive

The `build_sido_reactive` method is intended to be used to rapidly construct a standard set of material balance equations for zero-order type models with a single inlet and two outlets that involve chemical reactions (i.e. units that involve both reaction and separation).

Usage
-----

.. doctest::

  from idaes.core import declare_process_block_class
  from watertap.core import build_sido_reactive, ZeroOrderBaseData

  @declare_process_block_class("ReactiveSeparationZO")
  class ReactiveSeparationZOData(ZeroOrderBaseData):

      CONFIG = ZeroOrderBaseData.CONFIG()

      def build(self):
          super().build()

          self._tech_type = "reactive_separation"

          build_sido_reactive(self)

Model Structure
---------------

The build_sido_reactive method constructs a simple representation of unit operation with a single inlet (named `inlet`) and two outlets (named `treated` and `byproduct`). A `StateBlock` is constructed for the inlet and each outlet with a `Port` associated with each of these.

Variables
---------

The build_sido_reactive method creates the following variables in addition to those created by the `StateBlocks`.

=============================== ========================= =============== ======================================================================
Variable                        Name                      Indices         Notes
=============================== ========================= =============== ======================================================================
:math:`r_{t}`                   recovery_frac_mass_H2O    time            Fraction of mass flow of water in inlet that goes to treated stream.
:math:`f_{t,j}`                 removal_frac_mass_comp    time, component Fraction of mass flow of each component that goes to byproduct stream.
:math:`X_{t, r}`                reaction_conversion       time, reaction  Fractional conversion of key reactant via reaction r.
:math:`\xi_{t, r}`              extent_of_reaction        time, reaction  Extent of reaction r in kg/s.
=============================== ========================= =============== ======================================================================

Parameters
----------

The build_sido_reactive method introduces a single indexed parameter.

===================  ================= ==================== =========================================================================
Parameter            Name              Indices              Notes
===================  ================= ==================== =========================================================================
:math:`\nu_{r, j}`   generation_ratio  reaction, component  Mass conversion ratio for component j in reaction r w.r.t. key reactant.
===================  ================= ==================== =========================================================================

Note that :math:`\nu_{r, j}` should be negative for reactants and positive for products. Generation ratio can either be specified directly as a parameter, or by providing stoichiometric coefficients and molecular weights for each component.

Expressions
-----------

The build_sido_reactive method creates an Expression for the generation of each species in each reaction. 

===================  ==================== ========================== ===============================================
Parameter            Name                 Indices                    Notes
===================  ==================== ========================== ===============================================
:math:`G_{t, r, j}`  generation_rxn_comp  time, reaction, component  Mass generation for component j in reaction r.
===================  ==================== ========================== ===============================================

.. math:: G_{t, r, j} = \nu_{r,j} \times \xi_{t,r}

Constraints
-----------

The build_sido method writes the following constraints which relate the inlet state to those in the treated and byproduct streams.
First, a constraint for the extent of reaction is written:

.. math:: \xi_{t, r} = X_{t, r} \times M_{in,t,key}

where :math:`M_{t,key}` is mass flowrate of the key reactant for the reaction at time :math:`t`.

Next, a mass recovery equation for water is written to relate the flowrate at the treated outlet to that at the inlet. Note that recovery is applied after any reactions have occurred.

`water_recovery_equation(t)`:

.. math:: r_t \times (M_{inlet,t,H2O} + \sum{G_{t, r, H2O}}_r) = M_{treated,t,H2O}

Next, a mass balance for water is written to relate the flowrate in the byproduct stream to that in the inlet and treated streams.

`flow_balance(t)`:

.. math:: M_{inlet,t,H2O} + \sum{G_{t, r, H2O}}_r = M_{treated,t,H2O} + M_{byproduct,t,H2O}

Removal constraints are then written for each solute to relate the solute concentration in the byproduct stream to that in the inlet stream. Again, removal is applied after any reactions have occured.

`solute_removal_equation(t, j)`:

.. math:: f_{t, j} \times (M_{inlet,t,j} + \sum{G_{t, r, j}}_r) = M_{byproduct,t,j}

A mass balance constraint is then written for each solute.

`solute_treated_equation(t, j)`:

.. math:: M_{inlet,t,j} + \sum{G_{t, r, j}}_r = M_{treated,t,j} + M_{byproduct,t,j}

Module Documentation
--------------------

* :mod:`watertap.core.zero_order_sido_reactive`
