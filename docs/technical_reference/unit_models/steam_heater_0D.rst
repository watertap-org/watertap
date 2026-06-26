.. _steam_heater_0D:

Steam Heater 0D
===============

.. code-block:: python

   from watertap.unit_models.steam_heater_0D import SteamHeater0D, Mode

.. index::
   pair: watertap.unit_models.steam_heater_0D;steam_heater_0D

The Steam Heater 0D unit model is a 0D heat exchanger model based on the
`IDAES heat exchanger <https://idaes-pse.readthedocs.io/en/latest/reference_guides/model_libraries/generic/unit_models/heat_exchanger.html>`_.
The model is intended for cases where the hot side represents steam that
condenses while transferring heat to a cold-side process stream.

The model extends the IDAES heat exchanger with WaterTAP-specific constraints
for the condensing hot side:

* vapor flow in the hot-side outlet is fixed to its lower bound
* hot-side inlet vapor and liquid flow are balanced to the hot-side outlet liquid flow
* hot-side outlet pressure is constrained to be greater than or equal to the
  saturation pressure of the hot-side outlet state

Configuration
-------------

The Steam Heater 0D model supports the following additional configuration options:

.. csv-table::
   :header: "Configuration option", "Description"

   "``mode``", "Initialization mode. Options are ``Mode.HEATER`` and ``Mode.CONDENSER``."
   "``estimate_cooling_water``", "When ``mode`` is ``Mode.CONDENSER``, optionally estimate the cold-side inlet flow rate from a specified cold-side outlet temperature."

The configuration options affect the initialization routine. The unit model
equations are built from the base heat exchanger equations plus the condensing
hot-side constraints listed above.

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

When using ``Mode.HEATER``, the initialization routine can unfix the hot-side
inlet vapor flow so the model can estimate the required steam flow. When using
``Mode.CONDENSER`` with ``estimate_cooling_water=True``, the initialization
routine can estimate the cold-side inlet flow from a fixed cold-side outlet
temperature.

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
       mode=Mode.HEATER,
   )

Equations and Relationships
---------------------------

In addition to the standard IDAES heat exchanger equations, Steam Heater 0D adds
the following relationships on the hot side:

.. csv-table::
   :header: "Description", "Constraint name"

   "Condensed hot-side material balance", "``outlet_liquid_mass_balance``"
   "Hot-side outlet pressure lower bounded by saturation pressure", "``outlet_pressure_sat``"

The model uses the default heat exchanger costing method from
``watertap.costing.unit_models.heat_exchanger``.

Class Documentation
-------------------

.. currentmodule:: watertap.unit_models.steam_heater_0D

.. autoclass:: SteamHeater0D
    :members:
    :noindex:

.. autoclass:: Mode
    :members:
    :noindex:
