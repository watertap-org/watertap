.. _electroNP:

ElectroNP (ZO)
==============

.. code-block:: python

   from watertap.unit_models.electroNP_ZO import ElectroNPZO

This electrochemical nutrient removal (ElectroNP) unit model
   * is a zero-order model that enables the user to specify performance in terms of solute rejection
   * supports a single liquid phase only
   * supports steady-state only
   * assumes isothermal conditions

.. index::
   pair: watertap.unit_models.electroNP_ZO;electroNP_ZO

.. currentmodule:: watertap.unit_models.electroNP_ZO

Degrees of Freedom
------------------
The zero-order electroNP model has at least 6 degrees of freedom that should be fixed for the unit to be fully specified.
Typically, the following variables are fixed:

   * inlet volumetric flow rate
   * inlet temperature
   * inlet pressure
   * electricity intensity
   * dosage of magnesium chloride per phosphorus removal

There is 1 degree of freedom for each solute in a given property model:

   * inlet solute concentration

Model Structure
------------------
This electroNP model is built based on the `IDAES Separator module <https://idaes-pse.readthedocs.io/en/stable/reference_guides/model_libraries/generic/unit_models/separator.html>`_. It consists of one inlet port and two outlet ports as
treated and byproducts. The ``SplittingType`` is set to ``componentFlow``.

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Phases", ":math:`p`", "['Liq']"
   "Components", ":math:`j`", "['H2O', 'Solute']*"

\*Solute depends on the imported property model, and more than one solute can be provided.

Variables
----------

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Mass recovery fraction of water", ":math:`R_{H2O}`", "recovery_frac_mass_H2O", "[t]", ":math:`\text{dimensionless}`"
   "Removal fraction of component", ":math:`r_j`", "removal_frac_mass_comp", "[t, 'byproduct', j]", ":math:`\text{dimensionless}`"
   "Electricity consumption", ":math:`E`", "electricity", "[t]", ":math:`\text{kW}`"
   "Electricity intensity with respect to phosphorus removal", ":math:`EI`", "energy_electric_flow_mass", "None", ":math:`\text{kWh/kg}`"
   "Dosage of magnesium chloride per phosphorus removal", ":math:`D_{MgCl2}`", "magnesium_chloride_dosage", "None", ":math:`\text{dimensionless}`"
   "Magnesium chloride flowrate", ":math:`Q_{MgCl2}`", "magnesium_chloride_dosage", "[t]", ":math:`\text{kg/hr}`"

The electroNP model also has two mutable parameters:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Phosphorus removal fraction", ":math:`r_P`", "P_removal", "None", ":math:`\text{dimensionless}`"
   "Nitrogen removal fraction", ":math:`r_N`", "N_removal", "None", ":math:`\text{dimensionless}`"

**NOTE: Removal fractions for other components are fixed to 0**

Equations
-----------

.. csv-table::
   :header: "Description", "Equation"

   "Water removal", ":math:`r_{H2O} = 1 - R_{H2O}`"
   "Phosphorus removal",":math:`r_{S_{PO4}} = r_P`"
   "Nitrogen removal",":math:`r_{S_{NH4}} = r_N`"
   "Energy consumption",":math:`E = EI Q_{byproduct, S_{PO4}}`"
   "MgCl2 demand",":math:`Q_{MgCl2} = D_{MgCl2} Q_{byproduct, S_{PO4}}`"

Class Documentation
-------------------

* :mod:`watertap.unit_models.electroNP_ZO`