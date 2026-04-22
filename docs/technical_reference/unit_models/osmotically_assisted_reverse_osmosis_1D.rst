.. _OARO_1D:

Osmotically Assisted Reverse Osmosis (1D)
=========================================

.. code-block:: python

   from watertap.unit_models.osomtically_assisted_reverse_osmosis_1D import OsmoticallyAssistedReverseOsmosis1D

This osmotically assisted reverse osmosis (OARO) unit model
   * is 1-dimensional
   * supports a single liquid phase only
   * supports steady-state only
   * supports flat-sheet module designs
   * is based on the solution-diffusion model and film theory
   * assumes isothermal conditions
   * assumes the feed-side flows in the forward direction
   * assumes the permeate-side flows in the backwards direction

.. index::
   pair: watertap.unit_models.osmotically_assisted_reverse_osmosis_1D;osomtically_assisted_reverse_osmosis_1D

.. currentmodule:: watertap.unit_models.osmotically_assisted_reverse_osmosis_1D

Degrees of Freedom
------------------
Aside from the feed-side and permeate-side inlet state variables (i.e. temperature, pressure, component flowrates), the OARO model has
at least 4 degrees of freedom that should be fixed for the unit to be fully specified. Unlike RO, which only
accounts for concentration polarization on the feed side, the OARO model includes a structural parameter
variable, which is used to calculate the membrane-interface concentration on the permeate side.

Typically, the following variables are fixed for the OARO model, in addition to state variables at the inlet:
    * membrane water permeability, A
    * membrane salt permeability, B
    * membrane area
    * structural parameter

On the other hand, configuring the OARO unit to calculate concentration polarization effects, mass transfer
coefficient, and pressure drop would result in 6 additional degrees of freedom. In this case, in addition to the
previously fixed variables, we typically fix the following variables to fully specify the unit:

    * feed-side spacer porosity
    * feed-side channel height
    * permeate-side space porosity
    * permeate-side channel height
    * feed-side Reynolds number *or* water mass recovery

Model Structure
------------------
This OARO model consists of a separate MembraneChannel1DBlock for the feed-side and the permeate-side of the OARO unit.

* The feed-side includes StateBlocks indexed by time and space (feed_side.properties[t, x]) which are used for mass, energy, and momentum balances, and additional StateBlocks for the conditions at the membrane interface (feed_side.properties_interface[t, x]).
* The permeate-side includes StateBlocks indexed by time and space (permeate-side.properties[t, x]) which are used for mass, energy, and momentum balances, and additional StateBlocks for the conditions at the membrane interface (permeate-side.properties_interface[t, x]).

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Inlet/outlet", ":math:`x`", "['in', 'out']"
   "Phases", ":math:`p`", "['Liq']"
   "Components", ":math:`j`", "['H2O', 'NaCl']*"

\*Solute depends on the imported property model; example shown here is for the NaCl property model.

.. _1doaro_variables:

Variables
----------

Refer to the :any:`0doaro_variables` section in the  0D OARO model.

.. _1doaro_equations:

Equations
-----------

Refer to the :any:`0doaro_equations` section in the  0D OARO model.

Class Documentation
-------------------

* :mod:`watertap.unit_models.osmotically_assisted_reverse_osmosis_1D`
