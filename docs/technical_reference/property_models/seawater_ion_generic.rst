Seawater Ion Generic Property Package
=====================================

This package implements property relationships for artificial seawater as provided in `Islam et al. (2021) <https://www.sciencedirect.com/science/article/abs/pii/S1383586621009412>`_.

This artificial seawater property package:
   * supports only 'H2O', 'Na_+', 'Ca_2+', 'Mg_2+', 'Cl\_-', 'SO4_2-'
   * supports only liquid phase
   * is formulated on a molar basis
   * does not support dynamics

This package uses the eNRTL equation of state and has been used to predict mineral scaling in full treatment
where reverse osmosis is the primary desalination technology. Additionally, this property package is
formulated as an extension of `IDAES's modular property framework <https://idaes-pse.readthedocs.io/en/stable/explanations/components/property_package/general/index.html#generic-property-package-framework>`_ (GenericParameterBlock).
In other words, this property package is a configuration file that gets passed into IDAES's generic/modular property framework.
Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Components", ":math:`j`", "['H2O', 'Na_+', 'Ca_2+', 'Mg_2+', 'Cl\_-', 'SO4_2-']"
   "Phases", ":math:`p`", "['Liq']"

State variables
---------------
.. csv-table::
   :header: "Description", "Symbol", "Variable", "Index", "Units"

   "Component molar flowrate", ":math:`N_j`", "flow_mol_phase_comp", "[p, j]", ":math:`\text{mole/s}`"
   "Temperature", ":math:`T`", "temperature", "None", ":math:`\text{K}`"
   "Pressure", ":math:`P`", "pressure", "None", ":math:`\text{Pa}`"

Properties
----------
.. csv-table::
   :header: "Description", "Symbol", "Variable", "Index", "Units"

   "Component mole fraction", ":math:`y_j`", "mole_frac_phase_comp", "[p, j]", ":math:`\text{dimensionless}`"
   "Component mole flowrate", ":math:`N_j`", "flow_mol_phase_comp", "[p, j]", ":math:`\text{mole/s}`"
   "Interaction parameter", ":math:`τ`", "Liq_tau", "[j_1, j_2, j_3]", ":math:`\text{dimensionless}`"
   "Component molar volume", ":math:`V_j`", "vol_mol_liq_comp", "[j]", ":math:`\text{m}^3\text{/mole}`"

Relationships
-------------
.. csv-table::
   :header: "Description", "Equation"

   "Component mole flowrate", ":math:`N_j = \frac{M_j}{MW_j}`"
   "Component mole fraction", ":math:`y_j = \frac{N_j}{\sum_{j} N_j}`"



Scaling
-------
This artificial seawater property package includes support for scaling, such as default or user-supplied scaling factors for all variables.

The default scaling factors are as follows:

   * 1e1 for 'Na_+' mole flowrate
   * 1e3 for 'Ca_2+' mole flowrate
   * 1e2 for 'Mg_2+' mole flowrate
   * 1e2 for 'SO4_2-' mole flowrate
   * 1e1 for 'Cl\_-' mole flowrate
   * 1e-1 for 'H2O' mole flowrate
   * 1e2 for 'Na_+' mole fraction
   * 1e4 for 'Ca_2+' mole fraction
   * 1e3 for 'Mg_2+' mole fraction
   * 1e3 for 'SO4_2-' mole fraction
   * 1e2 for 'Cl\_-' mole fraction
   * 1 for 'H2O' mole fraction
   * 1e1 for 'NaCl' apparent mole flowrate
   * 1e2 for 'Na2SO4' apparent mole flowrate
   * 1e2 for 'CaCl2' apparent mole flowrate
   * 1e3 for 'CaSO4' apparent mole flowrate
   * 1e2 for 'MgCl2' apparent mole flowrate
   * 1e3 for 'MgSO4' apparent mole flowrate
   * 1e-1 for 'H2O' apparent mole flowrate
   * 1e3 for 'NaCl' apparent mole fraction
   * 1e4 for 'Na2SO4' apparent mole fraction
   * 1e4 for 'CaCl2' apparent mole fraction
   * 1e5 for 'CaSO4' apparent mole fraction
   * 1e4 for 'MgCl2' apparent mole fraction
   * 1e5 for 'MgSO4' apparent mole fraction
   * 1 for 'H2O' apparent mole fraction

Scaling factors for other variables can be calculated based on their relationships with the user-supplied or default scaling factors.
   
Reference
---------

M.R.Islam, I.Hsieh, B.Lin, A.K.Thakur, C.Chen, M.Malmali, Molecular thermodynamics for scaling prediction: Case of membrane distillation, Separation and Purification Technology, 2021,Vol. 276. https://www.sciencedirect.com/science/article/abs/pii/S1383586621009412

