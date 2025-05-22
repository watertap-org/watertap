.. _ERD:

Energy Recovery Device
======================

.. code-block:: python

   from watertap.unit_models.pressure_changer import EnergyRecoveryDevice

.. index::
   pair: watertap.unit_models.energy_recovery_device;energy_recovery_device

This energy recovery device unit model represents a reversible pump or Pelton turbine.
It is a class of pressure changer which is based on the `IDAES pressure changer <https://idaes-pse.readthedocs.io/en/latest/reference_guides/model_libraries/generic/unit_models/pressure_changer.html>`_
and makes the following assumptions:

   * supports a single liquid phase only
   * supports steady-state only
   * supports isothermal pump only

Degrees of Freedom
------------------

Energy recovery device units generally have two degrees of freedom:

* outlet pressure, :math:`P_{out}` or pressure change across the inlet and outlet,:math:`\Delta P`
* energy recovery device efficiency

Model Structure
---------------

The energy recovery device unit model consists of a single ``ControlVolume0D`` (named ``control_volume``)
with one Inlet Port (named ``inlet``) and one Outlet Port (named ``outlet``).
This model inherits `IDAES pressure changer <https://idaes-pse.readthedocs.io/en/latest/reference_guides/model_libraries/generic/unit_models/pressure_changer.html>`_
with the configuration ``compressor`` fixed to ``False``, representing that the unit should be considered as an expander with pressure drop.

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Phases", ":math:`p`", "['Liq']"
   "Components", ":math:`j`", "['H2O', 'NaCl']*"

\*Solute depends on the imported property model; example shown here is for the NaCl property model.

Variables
---------

The variables are the same as those listed in the "Variable" tab of
`IDAES pressure changer <https://idaes-pse.readthedocs.io/en/latest/reference_guides/model_libraries/generic/unit_models/pressure_changer.html>`_.

Equations and Relationships
---------------------------
The constraints are the same as those listed in
`IDAES pressure changer <https://idaes-pse.readthedocs.io/en/latest/reference_guides/model_libraries/generic/unit_models/pressure_changer.html>`_.

Class Documentation
-------------------
.. currentmodule:: watertap.unit_models.pressure_changer

.. autoclass:: EnergyRecoveryDevice
    :members:
    :noindex:
