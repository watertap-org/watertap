Single Inlet - Double Outlet Helper Methods
===========================================

.. index::
   pair: watertap.core.zero_order_sido;build_sido

.. currentmodule:: watertap.core.zero_order_sido

The `build_sido` method is intended to be used to rapidly construct a standard set of material balance equations for zero-order type models with a signle inlet and two outlets.

Usage
-----

.. testcode::

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

=============================== ==================== ============== ===================================================================
Variable                        Name                 Indices        Notes
=============================== ==================== ============== ===================================================================
:math:`r_{t}`                   recovery_vol         time           Fraction of volumetric flow in inlet that goes to treated
:math:`f_{t,j}`                 removal_mass_solute  time, solutes  Fraction of mass flow of each solute that goes to byproduct stream
:math:`\Delta P_{treated,t}`    deltaP_treated       time           Pressure difference between inlet and treated stream. Optional.
:math:`\Delta P_{byproduct,t}`  deltaP_byproduct     time           Pressure difference between inlet and byproduct stream. Optional.
=============================== ==================== ============== ===================================================================

The :math:`\Delta P_{treated,t}` and :math:`\Delta P_{byproduct,t}` terms are optional and are only created if the `self._has_deltaP_treated` and `self._has_deltaP_byproduct` attributes are `True` respectively.

Constraints
-----------

The build_sido method writes the following constraints which relate the inlet state to those in the treated and byproduct streams. As mentioned previously, these constraints are formulated such that they will be linear if `flow_vol`, `conc_mass_comp`, `pressure` and `temperature` are the state variables used by the property package and the recovery and removal fractions are fixed.

First, a volumetric flow recovery equation is written to relate the flowrate at the treated to that at the inlet:

`water_recovery_equation(t)`:

.. math:: r_t \times Q_{inlet,t} = Q_{treated,t}

where :math:`Q_t` is volumetric flowrate at time :math:`t`.

Next, a volumetric flow balance is written to relate the flowrate in the byproduct stream to that in the inlet and treated.

`flow_balance(t)`:

.. math:: Q_{inlet,t} = Q_{treated,t} + Q_{byproduct,t}

Removal constraints are then written for each solute to relate the solute concentration in the byproduct stream to that in the inlet stream.

`solute_removal_equation(t, j)`:

.. math:: f_{t, j} \times C_{inlet,t, j} = (1-r_t) \times C_{byproduct,t,j}

where :math:`C_{t,j}` is mass concentration of solute :math:`j` at time :math:`t`.

A similar constraint is then written for each solute to relate the solute concentration in the treated stream to that in the inlet stream.

`solute_treated_equation(t, j)`:

.. math:: (1-f_{t, j}) \times C_{inlet,t, j} = r_t \times C_{treated,t,j}

Pressure balances are then written for each steam, including the :math:`\Delta P` terms if present.

`treated_pressure_constraint(t)`:

.. math:: P_{inlet,t} = P_{treated,t} [+ \Delta P_{treated,t}]

`byproduct_pressure_constraint(t)`:

.. math:: P_{inlet,t} = P_{byproduct,t} [+ \Delta P_{byproduct,t}]

where :math:`P_t` is the pressure at time :math:`t`.

Finally, temperature equality constraints are written to equate the temperature in the treated and byproduct streams to that in the inlet.

`treated_temperature_equality(t)`:

.. math:: T_{inlet,t} = T_{treated,t}

`byproduct_temperature_equality(t)`:

.. math:: T_{inlet,t} = T_{byproduct,t}

where :math:`T_t` is the temperature at time :math:`t`.

Module Documentation
--------------------

.. automodule:: watertap.core.zero_order_sido.build_sido
