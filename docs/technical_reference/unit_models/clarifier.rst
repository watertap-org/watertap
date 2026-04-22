.. _clarifier:

Clarifier
=========

.. code-block:: python

   from watertap.unit_models.clarifier import Clarifier

.. index::
   pair: watertap.unit_models.clarifier;clarifier

This clarifier unit model is based on the `IDAES separator <https://idaes-pse.readthedocs.io/en/latest/reference_guides/model_libraries/generic/unit_models/separator.html>`_
and makes the following assumptions:

   * supports a single liquid phase only
   * supports steady-state only

Degrees of Freedom
------------------
Clarifier units have a number of degrees of freedom based on the separation type chosen, where the default configuration is `componentFlow`.

* If `split_basis` = 'componentFlow', degrees of freedom are generally :math:`(no. outlets-1) \times no. components`
* If `split_basis` = 'phaseFlow', degrees of freedom are generally :math:`(no. outlets-1) \times no. phases`
* If `split_basis` = 'phaseComponentFlow', degrees of freedom are generally :math:`(no. outlets-1) \times no. phases \times no. components`
* If `split_basis` = 'totalFlow', degrees of freedom are generally :math:`(no. outlets-1) \times no. phases \times no. components`


Model Structure
---------------
The clarifier unit model does not use ControlVolumes, and instead writes a set of material, energy and momentum balances to split the inlet stream into a number of outlet streams.
There is a single inlet port (named inlet) and a user-defined number of outlet ports.

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
.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Surface area", ":math:`A_{s}`", "surface_area", "None", ":math:`\text{m}^2`"
   "Electricity intensity", ":math:`E_{I}`", "energy_electric_flow_vol_inlet", "None", ":math:`\text{kWh/}\text{m}^3`"

Equations and Relationships
---------------------------

.. csv-table::
   :header: "Description", "Equation"

   "Electricity consumption", ":math:`E = E_{I} * Q_{in}`"

Class Documentation
-------------------
.. currentmodule:: watertap.unit_models.clarifier

.. autoclass:: Clarifier
    :members:
    :noindex:
