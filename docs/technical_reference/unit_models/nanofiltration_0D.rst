Nanofiltration (0D)
====================
.. code-block:: python

   from watertap.unit_models import Nanofiltration0D

This nanofiltration (NF) unit model is suitable for non-predictive case studies where the user wishes to specify constant rejection fractions
whilst preserving electroneutrality in the permeate stream (optional). A single recovery fraction is assumed for all solvents, whilst
individual rejection fractions can be set for all solutes, with the exception of one ion specified for maintaining
electroneutrality. Solutes are assigned default rejection value by grouping them into either two categories, passing and excluded, with
default rejection values for each category. Solutes are categorised using either a user-provided list of species which freely pass the
membrane, or by assuming all neutral and monovalent species pass the membrane.

Retentate pressure is assumed to be related to feed pressure with an optional pressure drop whilst permeate pressure is assumed to be
a degree of freedom, and temperature equality is assumed. Temperature and pressure constraints can be removed with configuration arguments.

This model supports a single liquid phase only and assumes steady-state.

.. index::
   pair: watertap.unit_models.nanofiltration_0D;nanofiltration_0D

.. currentmodule:: watertap.unit_models.nanofiltration_0D

Degrees of Freedom
------------------
The ``Nanofiltration0D`` model has the following degrees of freedom

   * inlet state
   * permeate pressure
   * solvent recovery fraction (``recovery_solvent``)
   * solute rejection fractions (``rejection_comp``) for all solutes except the one identified as the electroneutrality species.
   * membrane area (``area``)

The following additional degrees of freedom may exist depending on configuration options:

   * retentate pressure drop (``deltaP``)

Model Structure
---------------
The Nanofiltration0D model consists of three state blocks (``properties_in``, ``properties_retentate`` and ``properties_permeate``).

.. _NF0D_variables:

Variables
---------

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Solvent recovery", ":math:`Q_{solvent}`", "``recovery_solvent``", None, ":math:`\text{dimensionless}`"
   "Solute rejection", ":math:`r_j`", "``rejection_comp``", [j], ":math:`\text{dimensionless}`"
   "Retentate pressure drop", ":math:`\Delta P`", "``deltaP``", [t], ":math:`\text{Pa}`"
   "Membrane area", ":math:`A_m`", "``area``", "None", ":math:`\text{m}^2`"

.. _NF0D_equations:

Equations
---------

Here :math:`F` represents component flowrate and :math:`Z_j` is the charge on species j.

.. csv-table::
   :header: "Description", "Equation"

   "Component material balances", ":math:`F_{in, j} = F_{retentate, j} + F_{perm, j}`"
   "Solvent recovery",":math:`Q_{solvent} \times F_{in, j} = F_{p, j}`"
   "Solute rejection",":math:`r_j = 1 - \frac{C_{perm, j}}{C_{feed_j}}`"
   "Permeate electroneutrality",":math:`0 = F_{perm, j} \times Z_{j}`"
   "Retentate pressure balance",":math:`P_{in} = P_{retentate} - \Delta P`"
   "Retentate temperature equality", ":math:`T_{in} = T_{retentate}`"
   "Permeate temperature equality", ":math:`T_{in} = T_{perm}`"
   "Solvent mass transfer", ":math:`M_{p, solv} = A_m J_{solv} \rho_{solvent}`"
   "Component recovery rate",":math:`R_j = \frac{M_{p,j}}{M_{f,in,j}}`"
   "Volumetric recovery rate",":math:`R_{vol} = \frac{Q_{p}}{Q_{f,in}}`"



References
----------

| O. Labban, C. Liu, T. H. Chong and J. H. Lienhard V (2017)
| Fundamentals of low-pressure nanofiltration: Membrane characterization, modeling, and understanding the multi-ionic interactions in water softening
| *Journal of Membrane Science* 2017 Vol. 521 Pages 18-32
| doi: 10.1016/j.memsci.2016.08.062


Class Documentation
-------------------

* :mod:`watertap.unit_models.nanofiltration_0D`
