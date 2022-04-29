ASM1 Property Package
=====================

This package implements property relationships for the treatment of wastewater using an activated sludge biological reactor as provided in `Henze, M. et al. (1987) <https://belinra.inrae.fr/doc_num.php?explnum_id=4467>`_.

This Activated Sludge Model no.1 (ASM1) property package:
   * supports only 'H2O' and 'Sludge' as Components
   * supports only liquid phase
   * is formulated on a mass basis

Sets
----
.. csv-table::
  :header: "Description", "Symbol", "Indices"

  "Components", ":math:`j`", "['H2O', 'Sludge']"
  "Phases", ":math:`p`", "['Liq']"

State variables
---------------
.. csv-table::
   :header: "Description", "Symbol", "Variable", "Index", "Units"

   "Total volumetric flowrate", ":math:`Q`", "flow_vol", "None", ":math:`\text{m}^3\text{/s}`"
   "Temperature", ":math:`T`", "temperature", "None", ":math:`\text{K}`"
   "Pressure", ":math:`P`", "pressure", "None", ":math:`\text{Pa}`"
   "Component mass concentrations", ":math:`C_j`", "conc_mass_comp", "[p]", ":math:`\text{kg/}\text{m}^3`"
   "Alkalinity in molar concentration", ":math:`A`", "alkalinity", "[p]", ":math:`\text{kmol/m}^{3}`"

Parameters
----------
.. csv-table::
 :header: "Description", "Symbol", "Parameter", "Index", "Units"

   "Y_A (0.24)", ":math:`Y_A`", "Y_A", "[p]", ":math:`\text{g COD formed/}\text{g N oxidized}`"
   "Y_H (0.67)", ":math:`Y_H`", "Y_H", "[p]", ":math:`\text{g COD formed/}\text{g COD oxidized}`"
   "f_P (0.08)", ":math:`f_P`", "f_P", "[p]", ":math:`\text{dimensionless}`"
   "i_XB (0.08)", ":math:`i_{XB}`", "i_XB", "[p]", ":math:`\text{g N/}\text{g COD in biomass}`"
   "i_XP (0.06)", ":math:`i_{XP}`", "i_XP", "[p]", ":math:`\text{g N/}\text{g COD in particulate products}`"
   "µ_H (4.0)", ":math:`µ_H`", "µ_H", "[p]", ":math:`\text{d}^{-1}`"
   "K_S (10.0)", ":math:`K_S`", "K_S", "[p]", ":math:`\text{g COD/}\text{m}^{3}`"
   "K_OH (0.2)", ":math:`K_{O,H}`", "K_OH", "[p]", ":math:`\text{g -COD/}\text{m}^{3}`"
   "K_NO (0.5)", ":math:`K_{NO}`", "K_NO", "[p]", ":math:`\text{g NO}_{3}\text{-N/}\text{m}^{3}`"
   "b_H (0.3)", ":math:`b_H`", "b_H", "[p]", ":math:`\text{d}^{-1}`"
   "η_g (0.8)", ":math:`η_g`", "η_g", "[p]", ":math:`\text{dimensionless}`"
   "η_h (0.8)", ":math:`η_h`", "η_h", "[p]", ":math:`\text{dimensionless}`"
   "k_h (3.0)", ":math:`k_h`", "k_h", "[p]", ":math:`\text{g slowly biodegradable COD/}\text{g COD . d}`"
   "K_X (0.1)", ":math:`K_X`", "K_X", "[p]", ":math:`\text{g slowly biodegradable COD/}\text{g COD}`"
   "µ_A (0.5)", ":math:`µ_A`", "µ_A", "[p]", ":math:`\text{d}^{-1}`"
   "K_NH (1.0)", ":math:`K_{NH}`", "K_NH", "[p]", ":math:`\text{g NH}_{3}\text{-N/}\text{m}^{3}`"
   "b_A (0.05)", ":math:`b_A`", "b_A", "[p]", ":math:`\text{d}^{-1}`"
   "K_OA (0.4)", ":math:`K_{O,A}`", "K_OA", "[p]", ":math:`\text{g -COD/}\text{m}^{3}`"
   "k_a (0.05)", ":math:`k_a`", "k_a", "[p]", ":math:`\text{m}^{3}\text{/}\text{g COD . d}`"


Properties
----------
.. csv-table::
  :header: "Description", "Symbol", "Variable", "Index", "Units"

  "Fluid specific heat capacity", ":math:`c_p`", "cp", "None", ":math:`\text{J/kg/K}`"
  "Mass density", ":math:`\rho`", "dens_mass", "[p]", ":math:`\text{kg/}\text{m}^3`"
  "Soluble inert organic matter, S_I", ":math:`S_I`", "S_I", "[p]", ":math:`\text{kg/m}^{3}`"
  "Readily biodegradable substrate S_S", ":math:`S_S`", "S_S", "[p]", ":math:`\text{kg/m}^{3}`"
  "Particulate inert organic matter, X_I", ":math:`X_I`", "X_I", "[p]", ":math:`\text{kg/m}^{3}`"
  "Slowly biodegradable substrate X_S", ":math:`X_S`", "X_S", "[p]", ":math:`\text{kg/m}^{3}`"
  "Active heterotrophic biomass X_B,H", ":math:`X_{B,H}`", "X_BH", "[p]", ":math:`\text{kg/m}^{3}`"
  "Active autotrophic biomass X_B,A", ":math:`X_{B,A}`", "X_BA", "[p]", ":math:`\text{kg/m}^{3}`"
  "Particulate products arising from biomass decay, X_P", ":math:`X_P`", "X_P", "[p]", ":math:`\text{kg/m}^{3}`"
  "Oxygen, S_O", ":math:`S_O`", "S_O", "[p]", ":math:`\text{kg/m}^{3}`"
  "Nitrate and nitrite nitrogen, S_NO", ":math:`S_{NO}`", "S_NO", "[p]", ":math:`\text{kg/m}^{3}`"
  "NH4 :math:`^{+}` + NH :math:`_{3}` Nitrogen, S_NH", ":math:`S_{NH}`", "S_NH", "[p]", ":math:`\text{kg/m}^{3}`"
  "Soluble biodegradable organic nitrogen, S_ND", ":math:`S_{ND}`", "S_ND", "[p]", ":math:`\text{kg/m}^{3}`"
  "Particulate biodegradable organic nitrogen, X_ND", ":math:`X_{ND}`", "X_ND", "[p]", ":math:`\text{kg/m}^{3}`"
  "Alkalinity, S_ALK", ":math:`S_{ALK}`", "S_ALK", "[p]", ":math:`\text{kg/m}^{3}`"

Relationships
-------------
.. csv-table::
   :header: "Description", "Equation"

   "Aerobic growth of heterotrophs", ":math:`ρ_1 = µ_{H}(\frac{S_{S}}{K_{S}+S_{S}})(\frac{S_{O}}{K_{O,H}+S_{O}})X_{B,H}`"
   "Anoxic growth of heterotrophs", ":math:`ρ_2 = µ_{H}(\frac{S_{S}}{K_{S}+S_{S}})(\frac{K_{O,H}}{K_{O,H}+S_{O}})(\frac{S_{NO}}{K_{NO}+S_{NO}})η_{g}X_{B,H}`"
   "Aerobic growth of autotrophs", ":math:`ρ_3 = µ_{A}(\frac{S_{NH}}{K_{NH}+S_{NH}})(\frac{S_{O}}{K_{O,A}+S_{O}})X_{B,A}`"
   "Decay of heterotrophs", ":math:`ρ_4 = b_{H}X_{B,H}`"
   "Decay of autotrophs", ":math:`ρ_5 = b_{H}X_{B,H}`"
   "Ammonification of soluble organic nitrogen", ":math:`ρ_6 = k_{a}S_{ND}X_{B,H}`"
   "Hydrolysis of entrapped organics", ":math:`ρ_7 = k_{H}(\frac{X_{S}/X_{B,H}}{K_{X}+(X_{S}/X_{B,H})})[(\frac{S_{O}}{K_{O,H}+S_{O}})+η_{h}(\frac{K_{O,H}}{K_{O,H}+S_{O}})(\frac{S_{NO}}{K_{NO}+S_{NO}})]X_{B,H}`"
   "Hydrolysis of entrapped organic nitrogen", ":math:`ρ_7 = k_{H}(\frac{X_{S}/X_{B,H}}{K_{X}+(X_{S}/X_{B,H})})[(\frac{S_{O}}{K_{O,H}+S_{O}})+η_{h}(\frac{K_{O,H}}{K_{O,H}+S_{O}})(\frac{S_{NO}}{K_{NO}+S_{NO}})]X_{B,H}(X_{ND}/X_{S})`"


Scaling
-------
This ASM1 property package includes support for scaling, such as providing
default or calculating scaling factors for almost all variables. The only variables
that do not have scaling factors are the component mass flowrate and the user will
receive a warning if these are not set.

The user can specify the scaling factors for component mass flowrates with the following:

.. testsetup::

  from pyomo.environ import ConcreteModel
  from idaes.core import FlowsheetBlock

.. testcode::

  # relevant imports
  import watertap.property_models.coagulation_prop_pack as props    # Needs to be replaced with ASM prop pack
  from idaes.core.util.scaling import calculate_scaling_factors

  # relevant assignments
  m = ConcreteModel()
  m.fs = FlowsheetBlock(default={"dynamic": False})
  m.fs.properties = props.CoagulationParameterBlock()               # Needs to be replaced with ASM

  # set scaling for component mass flowrate
  m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq','H2O'))
  m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq','TDS'))
  m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq','TSS'))
  m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e3, index=('Liq','Sludge'))

  # calculate scaling factors
  calculate_scaling_factors(m.fs)

The default scaling factors are as follows:

  * 1e-2 for temperature
  * 1e-6 for pressure
  * 1e-3 for mass density

Scaling factors for other variables can be calculated based on their relationships
with the user-supplied or default scaling factors.


References
----------
[1] Henze, M., Grady, C.P.L., Gujer, W., Marais, G.v.R., Matsuo, T.,
"Activated Sludge Model No. 1", 1987, IAWPRC Task Group on Mathematical Modeling
for Design and Operation of Biological Wastewater Treatment.
https://belinra.inrae.fr/doc_num.php?explnum_id=4467

[2] Alex, J. et al. Benchmark Simulation Model no.1 (BSM1). Lund University, 2008, 5-6.
https://www.iea.lth.se/publications/Reports/LTH-IEA-7229.pdf