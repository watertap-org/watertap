.. _CSTR:

Continuously Stirred Tank Reactor
=================================

.. code-block:: python

   from watertap.unit_models.cstr import CSTR

.. index::
   pair: watertap.unit_models.cstr;cstr

This CSTR unit model is based on the `IDAES CSTR <https://idaes-pse.readthedocs.io/en/latest/reference_guides/model_libraries/generic/unit_models/cstr.html>`_
and makes the following assumptions:

   * is 0-dimensional
   * supports a single liquid phase only
   * supports steady-state only

Degrees of Freedom
------------------
Aside from the inlet feed state variables (i.e. temperature, pressure, component flowrates), the CSTR model has
one degree of freedom that should be fixed for the unit to be fully specified.

    * volume OR hydraulic retention time

Model Structure
---------------
The CSTR unit model consists of a single ControlVolume0D (named control_volume) with one inlet port (named inlet) and one outlet port (named outlet).

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
   :header: "Description", "Symbol", "Variable Name", "Index", "Value", "Units"

   "Hydraulic retention time", ":math:`HRT`", "hydraulic_retention_time", "[t]", "4", ":math:`\text{s}`"

Equations and Relationships
---------------------------

.. csv-table::
   :header: "Description", "Equation"

   "CSTR retention time", ":math:`HRT = V / Q_{in}`"

Class Documentation
-------------------
.. currentmodule:: watertap.unit_models.cstr

.. autoclass:: CSTR
    :members:
    :noindex:
