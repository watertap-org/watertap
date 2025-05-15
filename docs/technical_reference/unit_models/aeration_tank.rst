.. _aeration_tank:

Aeration Tank
=============

.. code-block:: python

   from watertap.unit_models.aeration_tank import AerationTank

.. index::
   pair: watertap.unit_models.aeration_tank;aeration_tank

This aeration tank unit model inherits from the :ref:`CSTR with injection <CSTR_injection>`.
The model makes the following assumptions:

   * oxygen is injected into the tank
   * is 0-dimensional
   * supports a single liquid phase only
   * supports steady-state only

Degrees of Freedom
------------------
Aside from the inlet feed state variables (i.e. temperature, pressure, component flowrates), the aeration tank model has
four degree of freedoms, and additional degrees of freedom may need to be specified depending on the configuration options.

    * volume OR hydraulic retention time
    * injection rates for all components
    * lumped mass transfer coefficient for oxygen (KLa)*
    * dissolved oxygen concentration at equilibrium (:math:`S_{O, eq}`)*

If heat transfer is included:
    * heat duty

If pressure change is included:
    * change in pressure (ΔP)

\*These degrees of freedom should only be included if the oxygen injection rate is not specified and needs to be calculated

Model Structure
---------------
The aeration tank unit model consists of a single ControlVolume0D (named control_volume) with one inlet port (named inlet) and one outlet port (named outlet).

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Inlet/outlet", ":math:`x`", "['in', 'out']"
   "Phases", ":math:`p`", "['Liq']"
   "Components", ":math:`j`", "['H2O', 'NaCl']*"

\*Solute depends on the imported property model; example shown here is for the NaCl property model.

Variables
---------
.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Hydraulic retention time", ":math:`HRT`", "hydraulic_retention_time", "[t]", ":math:`\text{s}`"
   "Dissolved oxygen concentration at equilibrium", ":math:`S_{O, eq}`", "S_O_eq", "None", ":math:`\text{kg/}\text{m}^3`"
   "Component injection", ":math:`I`", "injection", "[t, p, j]", ":math:`\text{kg/hr}`"
   "Electricity intensity", ":math:`E_{I}`", "energy_electric_flow_vol_inlet", "None", ":math:`\text{kWh/}\text{m}^3`"

Equations and Relationships
---------------------------

.. csv-table::
   :header: "Description", "Equation"

   "Aeration tank retention time", ":math:`HRT = V / Q_{in}`"
   "Oxygen mass transfer", ":math:`I_{O} = KLa * V * (S_{O, eq} - C_{O})`"
   "Electricity consumption (without aeration)", ":math:`E = E_{I} * Q_{in}`"
   "Electricity consumption (with aeration)", ":math:`E = \frac{S_{O, eq}}{1.8} * V * KLa`"

Class Documentation
-------------------
.. currentmodule:: watertap.unit_models.aeration_tank

.. autoclass:: AerationTank
    :members:
    :noindex:
