Nanofiltration (0D)
====================
This nanofiltration (NF) unit model is suitable for non-predictive case studies where the user wishes to specify constant recovery fractions
whilst preserving electroneutrality in the permeate stream (optional). A single recovery fraction is assumed for all solvents, whilst
individual recovery fractions can be set for neutral and monovalent solutes, with the exception of one ion specified for maintaining
electroneutrality. By default, complete rejection of all multi-valent solutes is assumed (however a single multi-valent ion recovery fraction
is provided).

Retentate pressure is assumed to be related to feed pressure with an optional pressure drop whilst permeate pressure is assumed to be
a degree of freedom, and temperature equality is assumed. Temperature and pressure constraints can be removed with configuration arguments.

This model assumes supports a single liquid phase only and assumes steady-state.

.. index::
   pair: watertap.unit_models.nanofiltration_0D;nanofiltration_0D

.. currentmodule:: watertap.unit_models.nanofiltration_0D

Degrees of Freedom
------------------
The ``Nanofiltration0D`` model has the following degrees of freedom

   * inlet state
   * permeate pressure
   * solvent recovery fraction (``solvent_recovery``)
   * solute recovery fractions (``solute recovery``) for all solutes EXCEPT the one identified as as the electroneutrality species.

The following additional degrees of freedom  may exist depending on configuration options

   * retentate pressure drop (``deltaP``)
   * multivalent ion recovery (``multivalent_recovery``). Always present but fixed to 1e-10 by default.

Model Structure
---------------
The Nanofiltration0D model consists of three state blocks (``properties_in``, ``properties_retentate`` and ``properties_permeate``).

.. _NF0D_variables:

Variables
---------

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Solvent recovery", ":math:`R_{solvent}`", "solvent_recovery", None, ":math:`\text{dimensionless}`"
   "Solute recovery", ":math:`R_j`", "solute_recovery", [j], ":math:`\text{dimensionless}`"
   "Multivalent ion recovery", ":math:`R_{multi}`", "multivalent_recovery", None, ":math:`\text{dimensionless}`"
   "Retentate pressure drop", ":math:`\deltaP`", "deltaP", [t], ":math:`\text{pressure}`"

.. _NF0D_equations:

Equations
---------

Here :math:`F` represents component flowrate and :math:`Z_j` is the charge on species j.

.. csv-table::
   :header: "Description", "Equation"

   "Component material balances", ":math:`F_{in, j} = F_{retentate, j} + F_{permeate, j}`"
   "Solvent recovery",":math:`R_{solvent} \times F_{in, j} = F_{permeate_j}`"
   "Solute recovery",":math:`R_j \times F_{in, j} = F_{permeate_j}`"
   "Multi-valent recovery",":math:`R_{multi} \times F_{in, j} = F_{permeate_j}`"
   "Permeate electroneutrality",":math:`0 = F_{permeate, j} \times Z_{j}`"
   "Retentate pressure balance",":math:`P_{in} = P_{retentate} - \deltaP`"
   "Retentate temperature equality", ":math:`T_{in} = T_{retentate}`"
   "Permeate temperature equality", ":math:`T_{in} = T_{permeate}`"

Class Documentation
-------------------

* :mod:`watertap.unit_models.nanofiltration_0D`
