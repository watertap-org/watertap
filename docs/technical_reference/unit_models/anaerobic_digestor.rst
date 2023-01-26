Anaerobic Digestor
==================

Introduction
------------

Anaerobic biological processes have been used in food and beverage production for centuries.
However, modern advances pertaining to anaerobic conversions have introduced various forms of
high-rate treatment processes that have proven to be particularly useful in wastewater treatment plants.
High organic loading rates and low sludge production gives anaerobic processes a significant advantage over other
biological unit operations, and, more importantly, the positive net energy production from the produced biogas
can help replace fossil fuel sources, lowering greenhouse gas emissions. While a wide variety of anaerobic
digestion models have been developed over the years, their widespread use has been made impossible by the
models' very specific natures. This generic anaerobic digestion model seeks to overcome this restriction by
limiting itself to only the major biological processes and excluding many of the more niche relationships in this
first version of the model. Likewise, while this implementation may not be as accurate as models tailored to their
specific applications, its simplicity makes it applicable for a wide variety of anaerobic processes, providing a
common basis by which future model development and validation studies can be compared.

In this implementation of the model, the user MUST provide two property packages - one for the liquid phase and
another for the vapor phase. This ADM1 model is based on the standard IDAES CSTR with the addition of mass transfer
terms and an extra port for the vapor outlet. The model relies on the following key assumptions:

   * supports steady-state only
   * liquid phase property package has a single phase named Liq
   * vapor phase property package has a single phase named Vap
   * liquid and vapor phase properties need not have the same component lists

.. figure:: ../../_static/unit_models/anaerobic_digestor.png
    :width: 400
    :align: center

    Figure 1. Schematic representation of an anaerobic digestor

.. index::
   pair: watertap.unit_models.anaerobic_digestor;anaerobic_digestor

.. currentmodule:: watertap.unit_models.anaerobic_digestor

Degrees of Freedom
------------------
Aside from the inlet feed state variables (i.e. temperature, pressure, component flowrates), the ADM1 model has
at least 5 degrees of freedom that should be fixed for the unit to be fully specified.

Typically, the following variables are fixed, in addition to state variables at the inlet:
    * cation concentration
    * anion concentration
    * liquid volume
    * vapor volume
    * liquid outlet temperature

Control Volumes
---------------

This model has two separate 0D control volumes for the liquid and vapor phases.

* Liquid
* Vapor

Ports
-----

This model provides three ports:

* inlet
* liquid_outlet
* vapor_outlet

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Inlet/outlet", ":math:`x`", "['in', 'out']"
   "Phases", ":math:`p`", "['Liq', 'Vap']"
   "Components", ":math:`j`", "['H2O', 'S_su', 'S_aa', 'X_fa', 'X_va', 'X_bu', 'X_pro', 'X_ac', 'S_h2', 'S_ch4', 'S_IC', 'S_IN', 'S_I', 'X_c', 'X_ch', 'X_pr', 'X_li', 'X_su', 'X_aa', 'X_fa', 'X_c4', 'X_pro', 'X_ac', 'X_h2', 'X_I', 'S_cat', 'S_an', 'S_co2']"
   "Ion", ":math:`j`", "['S_cat', 'S_an'] \  :sup:`*`"

**Notes**
 :sup:`*` Ion" is a subset of "Component" and uses the same symbol j.

.. _ADM1_variables:

Variables
----------

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Liquid phase mass transfer term", ":math:`J`", "liquid_phase.mass_transfer_term_j", "[t]", ":math:`\text{kg/s}`"
   "Liquid volume", ":math:`V_ad, liq`", "volume_AD", "[t]", ":math:`\text{m}^3`"
   "Vapor volume", ":math:`V_ad, vap`", "volume_vapor", "[t]", ":math:`\text{m}^3`"

.. _ADM1_equations:

Equations
-----------


Nomenclature
------------


Class Documentation
-------------------

* :mod:`watertap.unit_models.anaerobic_digester`
