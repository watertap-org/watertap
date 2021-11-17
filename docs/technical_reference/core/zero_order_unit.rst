Zero Order Unit Model Base Class
================================

.. index::
   pair: watertap.core.zero_order_unit;SITOBaseData

.. currentmodule:: watertap.core.zero_order_unit

The zero-order unit model base class is intended to form the basis for a library of generic single inlet-two outlet (SITO) units models using recovery and removal fractions to partition the incoming flow between the two outlets. This base class is not intended to be used by itself, and requires the user to define a derived class before it can be used in a flowsheet.

Usage
-----

.. testcode::

  # Import IDAES Core components
  from idaes.core import declare_process_block_class

  # Import base class
  from watertap.core.zero_order_unit import SITOBaseData

  # Create derived class with decorator
  @declare_process_block_class("DerivedModel")
  class DerivedModelData(SITOBaseData):

    def build(self):

      # Create attributes to indicate whether pressure change terms should be created
      # This must happen before calling super().build()
      self._has_deltaP_outlet = True
      self._has_deltaP_waste = True

      # Call super().build()
      super().build()

Model Structure
---------------

The SITO base class constructs a simple representation of unit operation with a single inlet (named `inlet`) and two outlets (named `outlet` and `waste`). A `StateBlock` is constructed for the inlet and each outlet with a `Port` associated with each of these.

Property Package Requirements
-----------------------------

The SITO base class makes a number of assumptions about the structure of the associated property package, and contains a number of checks to ensure the property package meets these requirements. The requirements for the property package are:

1. A single phase named `Liq`.
2. A single `Solvent` species named `H2O` (package must define a `solvent_set` with length that contains "H2O").
3. All other species are defined as `Solutes` (package must define a `solute_set` and the `component_list` must be the union of `solvent_set` and `solute_set`).
4. The property package must define `flow_vol`, `conc_mass_comp`, `pressure` and `temperature` as properties.

The SITO base class is formulated in such a way that if `flow_vol`, `conc_mass_comp`, `pressure` and `temperature` are the state variables used by the property package, then the resulting unit model will be linear if the removal and recovery fractions are fixed (or bilinear if they are not). The :ref:`ideal water <technical_reference/core/water_props:Ideal Water Properties>` property package has been designed specifically to be used with models deriving from the SITO base class and meets all of the above requirements.

Variables
---------

The SITO base class create the following variables in addition to those created by the `StateBlocks`.

============================ ==================== ============== ================================================================
Variable                     Name                 Indices        Notes
============================ ==================== ============== ================================================================
:math:`r_{t}`                recovery_vol         time           Fraction of volumetric flow in inlet that goes to outlet
:math:`f_{t,j}`              removal_mass_solute  time, solutes  Fraction of mass flow of each solute that goes to waste stream
:math:`\Delta P_{outlet,t}`  deltaP_outlet        time           Pressure difference between inlet and outlet stream. Optional.
:math:`\Delta P_{waste,t}`   deltaP_waste         time           Pressure difference between inlet and waste stream. Optional.
============================ ==================== ============== ================================================================

The :math:`\Delta P_{outlet,t}` and :math:`\Delta P_{waste,t}` terms are optional and are only created if the `self._has_deltaP_outlet` and `self._has_deltaP_waste` attributes are `True` respectively.

Constraints
-----------

The SITO base class writes the following constraints which relate the inlet state to those in the outlet and waste streams. As mentioned previously, these constraint are formulated such that they will be linear if `flow_vol`, `conc_mass_comp`, `pressure` and `temperature` are the state variables used by the property package and the recovery and removal fractions are fixed.

First, a volumetric flow recovery equation is written to relate the flowrate at the outlet to that at the inlet:

`water_recovery_equation(t)`:

.. math:: r_t \times Q_{inlet,t} = Q_{outlet,t}

where :math:`Q_t` is volumetric flowrate at time :math:`t`.

Next, a volumetric flow balance is written to relate the flowrate in the waste stream to that in the inlet and outlet.

`flow_balance(t)`:

.. math:: Q_{inlet,t} = Q_{outlet,t} + Q_{waste,t}

Removal constraints are then written for each solute to relate the solute concentration in the waste stream to that in the inlet stream.

`solute_removal_equation(t, j)`:

.. math:: f_{t, j} \times C_{inlet,t, j} = (1-r_t) \times C_{waste,t,j}

where :math:`C_{t,j}` is mass concentration of solute :math:`j` at time :math:`t`.

A similar constraint is then written for each solute to relate the solute concentration in the outlet stream to that in the inlet stream.

`solute_outlet_equation(t, j)`:

.. math:: (1-f_{t, j}) \times C_{inlet,t, j} = r_t \times C_{outlet,t,j}

Pressure balances are then written for each steam, including the :math:`\Delta P` terms if present.

`outlet_pressure_constraint(t)`:

.. math:: P_{inlet,t} = P_{outlet,t} [+ \Delta P_{outlet,t}]

`waste_pressure_constraint(t)`:

.. math:: P_{inlet,t} = P_{waste,t} [+ \Delta P_{waste,t}]

where :math:`P_t` is the pressure at time :math:`t`.

Finally, temperature equality constraints are written to equate the temperature in the outlet and waste streams to that in the inlet.

`outlet_temperature_equality(t)`:

.. math:: T_{inlet,t} = T_{outlet,t}

`waste_temperature_equality(t)`:

.. math:: T_{inlet,t} = T_{waste,t}

where :math:`T_t` is the temperature at time :math:`t`.

Class Documentation
-------------------

.. autoclass:: SITOBaseData
   :members:
