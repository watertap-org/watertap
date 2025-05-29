.. _thickener:

Thickener
=========

.. code-block:: python

   from watertap.unit_models.thickener import Thickener

.. index::
   pair: watertap.unit_models.thickener;thickener

The main assumptions of the implemented model are as follows:

1) Single liquid phase only
2) Steady state only
3) Has no volume

Introduction
------------
Sludge thickening is a process commonly used in wastewater treatment plants to increase the sludge concentration of a stream
as well as reduce its volume by removing free water. In doing so, the thickener reduces the load of downstream processes,
such as digestion and dewatering. This implementation of the thickener unit is based on the `IDAES separator unit <https://idaes-pse.readthedocs.io/en/stable/reference_guides/model_libraries/generic/unit_models/separator.html>`_.

Degrees of Freedom
------------------
The degrees of freedom in a thickener unit model are the inlet feed state variables:

    * temperature
    * pressure
    * component mass compositions

Model Structure
---------------
The thickener unit model does not use ControlVolumes, and instead writes a set of material, energy and momentum
balances to split the inlet stream into two outlet streams. Thickener models have a single inlet Port
(named inlet) and two outlet Ports (named overflow and underflow).

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Phases", ":math:`p`", "['Liq']"
   "Particulate Components", ":math:`j`", "['X_I', 'X_S', 'X_P', 'X_BH', 'X_BA', 'X_ND']"
   "Non-particulate Components", ":math:`j`", "['H2O', 'S_I', 'S_S', 'S_O', 'S_NO', 'S_NH', 'S_ND', 'S_ALK']"

NOTE: These components are defined in the :ref:`ASM1 Property Package <ASM1>` documentation.

Parameters
----------
.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Value", "Units"

   "Percentage of suspended solids in the underflow", ":math:`p_{thick}`", "p_thick", "None", "0.07", ":math:`\text{dimensionless}`"
   "Percentage of suspended solids removed", ":math:`TSS_{rem}`", "TSS_rem", "None", "0.98", ":math:`\text{dimensionless}`"


Equations and Relationships
---------------------------

.. csv-table::
   :header: "Description", "Equation"

   "Suspended solid concentration", ":math:`C_{TSS} = 0.75 (C_{X_{I}} + C_{X_{IP}} + C_{X_{BH}} + C_{X_{BA}} + C_{X_{S}})`"
   "Thickening factor", ":math:`f_{thick} = p_{thick} (\frac{10}{C_{TSS}})`"
   "Remove factor", ":math:`f_{q_{du}} = \frac{TSS_{rem}}{100 * f_{thick}}`"
   "Overflow particulate fraction", ":math:`split_{particulate} = 1 - TSS_{rem}`"
   "Overflow soluble fraction", ":math:`split_{soluble} = 1 - f_{q_{du}}`"

Class Documentation
-------------------
.. currentmodule:: watertap.unit_models.thickener

.. autoclass:: Thickener
    :members:
    :noindex:

References
----------
J. Alex, L. Benedetti, J.B. Copp, K.V. Gernaey, U. Jeppsson, I. Nopens, M.N. Pons, C. Rosen, J.P. Steyer & P. A. Vanrolleghem
Benchmark Simulation Model no. 2 (BSM2)