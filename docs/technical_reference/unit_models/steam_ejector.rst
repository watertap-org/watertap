Steam Jet Ejector
================================
This Steam Jet Ejector unit model:
   * Simulates the performance of a steam jet ejector for thermal vapor compression.
   * Uses semi-empirical correlations for entrainment ratio, pressure correction factor (PCF), and temperature correction factor (TCF), based on El-Dessouky (1997).
   * Operates in steady-state only.
   * Assumes the discharge mixture pressure equals its saturation pressure.
   
.. index::
   pair: watertap.unit_models.steam_ejector;steam_ejector

.. currentmodule:: watertap.unit_models.steam_ejector

Degrees of Freedom
-------------------
In addition to the inlet state variables (i.e., temperature, pressure, and component flowrates for motive steam and entrained vapor), the Steam Ejector model has at least 1 degree of freedom that must be fixed for the unit to be fully specified. Typically, the following variables are fixed:

    * Entrainment ratio
    * Compression ratio

Model Structure
------------------
This Steam Ejector model consists of state blocks for the properties of the motive steam inlet, entrained vapor inlet, and discharge mixture. It incorporates semi-empirical equations to model key performance parameters.



Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Inlet/Outlet", ":math:`x`", "['in', 'out']"
   "Phases", ":math:`p`", "['Liq', 'Vap']"
   "Components", ":math:`j`", "['H2O']"

Performance Metrics
--------------------
.. csv-table::
   :header: "Metric", "Equation"

   "Entrainment Ratio", ":math:`Ra = \frac{\dot{m}_{motive}}{\dot{m}_{entrained}}`"
   "Compression Ratio", ":math:`CR = \frac{P_s}{P_{ev}}`"

Variables
----------
.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Units", "Bounds"

   "Entrainment Ratio", ":math:`Ra`", "entrainment_ratio", "Dimensionless", "<4"
   "Compression Ratio", ":math:`CR`", "compression_ratio", "Dimensionless", ">1.89"
   "Pressure Correction Factor", ":math:`PCF`", "PCF", "Dimensionless", "N/A"
   "Temperature Correction Factor", ":math:`TCF`", "TCF", "Dimensionless", "N/A"
   "Motive Steam Pressure", ":math:`P_m`", "properties_motive_steam[0].pressure", "kPa", "[100, 3500]"
   "Entrained Vapor Pressure", ":math:`P_{ev}`", "properties_entrained_vapor[0].pressure", "kPa", "N/A"
   "Discharge Mixture Pressure", ":math:`P_s`", "properties_discharge_mix[0].pressure", "kPa", "N/A"

Equations
---------
.. csv-table::
   :header: "Description", "Equation"

   "Pressure Correction Factor", ":math:`PCF = 3 \times 10^{-7} P_m^2 - 0.0009 P_m + 1.6101`"
   "Temperature Correction Factor", ":math:`TCF = 2 \times 10^{-8} T_{ev}^2 - 0.0006 T_{ev} + 1.0047`"
   "Entrainment Ratio Model", ":math:`Ra \times TCF = 0.296 \frac{P_s^{1.19}}{P_{ev}^{1.04}} \left(\frac{P_m}{P_{ev}}\right)^{0.015} PCF`"
   "Entrainment Ratio Definition", ":math:`Ra = \frac{\dot{m}_{motive}}{\dot{m}_{entrained}}`"
   "Compression Ratio", ":math:`CR = \frac{P_s}{P_{ev}}`"
   

Class Documentation
-------------------
* :mod:`watertap.unit_models.steam_ejector`

References
----------
El-Dessouky, H., Modeling and simulation of thermal vapor compression desalination plant. Symposium on Desalination of Seawater with Nuclear Energy, Taejon, Republic of Korea, 26-30 May, 1997.