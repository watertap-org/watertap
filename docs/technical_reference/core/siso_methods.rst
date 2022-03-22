Single Inlet - Single Outlet Helper Methods
===========================================

.. index::
   pair: watertap.core.zero_order_siso;build_siso

.. currentmodule:: watertap.core.zero_order_siso

The `build_siso` method is intended to be used to rapidly construct a standard set of material balance equations for zero-order type models with a single inlet and single outlet.

Usage
-----

.. testcode::

  from idaes.core import declare_process_block_class
  from watertap.core import build_siso, ZeroOrderBaseData

  @declare_process_block_class("UVZO")
  class UVZOData(ZeroOrderBaseData):

      CONFIG = ZeroOrderBaseData.CONFIG()

      def build(self):
          super().build()

          self._tech_type = "uv"

          build_siso(self)

Model Structure
---------------

The build_siso method constructs a simple representation of unit operation with one inlet (named `inlet`) and one outlet (named `treated`). A `StateBlock` is constructed for each inlet and outlet with an associated `Port`.

Variables
---------

The build_siso method creates the following variables in addition to those created by the `StateBlocks`.

=============================== ========================= ============== =====================================================================
Variable                        Name                      Indices        Notes
=============================== ========================= ============== =====================================================================
:math:`r_{t}`                   recovery_frac_mass_H2O    time           Fraction of mass flow of water in inlet that goes to treated stream.
:math:`f_{t,j}`                 removal_frac_mass_solute  time, solutes  Fraction of mass flow of each solute that is removed from the inlet stream.
=============================== ========================= ============== =====================================================================

recovery_frac_mass_H2O is intended to be fixed to 1 (e.g., UV reactor which yields product stream without water losses), but the user can optionally set this to some fraction.

Constraints
-----------

The build_siso method writes the following constraints which relate the inlet state to the treated stream.
First, a water recovery equation is written for water to relate the flowrate at the treated outlet to that at the inlet:

`water_recovery_equation(t)`:

.. math:: r_t \times M_{inlet,t,H2O} = M_{treated,t,H2O}

where :math:`M_{t,H2O}` is mass flowrate of water at time :math:`t`.

Note, a mass balance for water is ignored since build_siso is intended to only account for constituent removal/conversion at the treated outlet.
Thus, a mass balance constraint is only written for each solute.

`solute_treated_equation(t, j)`:

.. math:: (1 - f_{t, j}) \times M_{inlet,t,j} = M_{treated,t,j}

Module Documentation
--------------------

.. automodule:: watertap.core.zero_order_siso.build_siso
    :members:
    :noindex:
