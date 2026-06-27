.. _steam_heater_0D:

Steam Heater 0D
===============

.. code-block:: python

   from watertap.unit_models.steam_heater_0D import SteamHeater0D

.. index::
   pair: watertap.unit_models.steam_heater_0D;steam_heater_0D

The Steam Heater 0D unit model is a 0D heat exchanger model based on the
`IDAES heat exchanger <https://idaes-pse.readthedocs.io/en/latest/reference_guides/model_libraries/generic/unit_models/heat_exchanger.html>`_.
The model is intended for cases where the hot side represents steam that
condenses while transferring heat to a cold-side process stream.

The model extends the IDAES heat exchanger with WaterTAP-specific constraints
for the condensing hot side:

* vapor flow in the hot-side outlet is fixed to its lower bound
* hot-side outlet pressure is constrained to the saturation pressure of the
  hot-side outlet state plus a non-negative pressure margin

The pressure margin is represented by ``pressure_deltaP``. It is fixed to zero
by default, and users can call ``set_subcooling_margin`` or
``release_subcooling_margin`` to fix or unfix the margin.

Degrees of Freedom
------------------

The degrees of freedom depend on the heat exchanger configuration and on whether
the model is being used to estimate steam flow, heat transfer area, or cold-side
flow. Typical specifications include:

* hot-side inlet temperature and pressure
* hot-side inlet liquid flow
* cold-side inlet temperature, pressure, and component mass flows
* overall heat transfer coefficient
* either heat transfer area or another design/operating variable that closes the
  heat exchanger degrees of freedom

Model Structure
---------------

The model inherits the IDAES ``HeatExchangerData`` structure. Users provide hot-
and cold-side property packages through the heat exchanger configuration. Common
WaterTAP usage names the two sides ``hot`` and ``cold``:

.. code-block:: python

   m.fs.unit = SteamHeater0D(
       hot_side_name="hot",
       cold_side_name="cold",
       hot={"property_package": m.fs.steam_properties},
       cold={"property_package": m.fs.cold_side_properties},
   )

Sets
----

.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Phases", ":math:`p`", "['Liq', 'Vap']*"
   "Components", ":math:`j`", "\*"

\*Phases and components depend on the imported property packages.

Variables
---------

The Steam Heater 0D model adds the following variable:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Pressure margin above hot-side saturation pressure", ":math:`\Delta P_{sat}`", "``pressure_deltaP``", "[t]", ":math:`\text{Pa}`"

Each hot-side and cold-side property block also contains the state variables
defined by the selected property packages. Common variables of interest include:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Phase-component mass flow", ":math:`M_{p,j}`", "``flow_mass_phase_comp``", "[p, j]", "\*"
   "Temperature", ":math:`T`", "``temperature``", "[t]", "\*"
   "Pressure", ":math:`P`", "``pressure``", "[t]", "\*"
   "Saturation pressure", ":math:`P_{sat}`", "``pressure_sat``", "[t]", "\*"

\*Units depend on the imported property packages.

Equations and Relationships
---------------------------

In addition to the standard IDAES heat exchanger equations, Steam Heater 0D adds
the following relationships on the hot side:

.. csv-table::
   :header: "Description", "Equation", "Model object"

   "Total condensation at the hot-side outlet", ":math:`M^{hot,out}_{Vap,j,t} = M^{hot,out,lb}_{Vap,j,t}`", "``hot_side.properties_out[t].flow_mass_phase_comp['Vap', j]``"
   "Hot-side outlet pressure set by saturation pressure and margin", ":math:`P^{hot,out}_{t} = P^{hot,out}_{sat,t} + \Delta P_{sat,t}`", "``outlet_pressure_sat``"

The model uses the default heat exchanger costing method from
``watertap.costing.unit_models.heat_exchanger``.

Class Documentation
-------------------

.. currentmodule:: watertap.unit_models.steam_heater_0D

.. autoclass:: SteamHeater0D
    :members:
    :noindex:
