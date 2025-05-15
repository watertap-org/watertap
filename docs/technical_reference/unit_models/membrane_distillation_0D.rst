.. _MD_0D:

Membrane Distillation (0D)
=========================================

.. code-block:: python

   from watertap.unit_models.MD.membrane_distillation_0D import MembraneDistillation0D

This Membrane Distillation (MD) unit model:
   * supports the following configurations: 

     - **DCMD** (Direct Contact Membrane Distillation)
     - **VMD** (Vacuum Membrane Distillation)
     - **GMD** (Permeate Gap/Conductive Gap Membrane Distillation)

   * is 0-dimensional
   * supports steady-state only
   * assumes heat loss in equipment is negligible
   * assumes permeate exits the membrane pores with zero salinity
   * assumes no concentration polarization for the cold channel
   * assumes complete vapor condensation for the cold channel (in DCMD and GMD)
   * accounts for vapor expansion in VMD
   * assumes linear temperature change across gap channel (in GMD)
   * assumes no pressure change and temperature polarization in VMD vacuum channel


Degrees of Freedom
------------------
In addition to the hot channel and cold channel inlet state variables (i.e, temperature, pressure, and component flowrates) for the **DCMD** and **GMD** configurations, the MD model has at least **4 degrees of freedom** for all configurations that should be fixed for the unit to be fully specified. Typically, the following variables are fixed:

- Membrane permeability coefficient
- Membrane thickness
- Membrane thermal conductivity
- Recovery *or* membrane area

**Additional degress of freedom**:

- **VMD** introduces vacuum pressure at the cold side.
- **GMD** introduces gap thermal conductivity and gap thickness.

Configuring the MD unit to calculate temperature polarization, concentration polarization, mass transfer
coefficient, and pressure drop would result in five additional degrees of freedom. In this case, in addition to the
previously fixed variables, we typically fix the following variables to fully specify the unit:

    * Hot channel spacer porosity
    * Hot channel height
    * Cold channel spacer porosity (in DCMD and GMD)
    * Cold channel height (in DCMD and GMD)
    * Membrane length *or* membrane width

Model Structure
---------------
The MD model consists of a separate `MDchannel0Dblock` for each channel depending on the configuration:

- **DCMD**: Includes **hot channel** and **cold channel**.
- **VMD**: Includes **hot channel** and **vacuum (cold) channel**.
- **GMD**: Includes **hot channel**, **gap channel**, and **cold channel**.

- **hot and cold channels in all configurations** includes bulk properties at the inlet and outlet (`properties_in` and `properties_out`) which are used for mass, energy, and momentum balances
- **hot channel in all configurations, cold channel in DCMD and GMD, and gap channel in GMD** includes 2 StateBlocks for the conditions at the membrane interface and gap interface
- **hot channel in all configurations, cold channel in DCMD, and gap channel in GMD** includes Vapor properties at the membrane interface (for **DCMD** and **VMD** configurations).

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Inlet/outlet", ":math:`x`", "['in', 'out']"
   "Phases", ":math:`p`", "['Liq', 'Vap']"
   "Components", ":math:`j`", "['H2O', solute]*"

\*Solute depends on the imported property model.

.. _0MD_variables:


Variables
---------
.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Membrane permeability coefficient", ":math:`B_0`", "permeability_coef", "[t]", ":math:`\text{kg/m/Pa/s}`"
   "Membrane thickness", ":math:`\sigma`", "membrane_thickness", "None", ":math:`\text{m}`"
   "Membrane thermal conductivity", ":math:`k_m`", "membrane_tc", "None", ":math:`\text{W/K/m}`"
   "Mass flux across membrane", ":math:`J`", "flux_mass", "[t, x]", ":math:`\text{kg/s}\text{/m}^2`"
   "Conduction heat flux across membrane", ":math:`q_{cond}`", "flux_conduction_heat", "[t, x]", ":math:`\text{W}\text{/m}^2`"
   "Evaporation heat flux from hot channel", ":math:`q_{evap}`", "flux_enth_hot", "[t, x]", ":math:`\text{W}\text{/m}^2`"
   "Condensation heat flux to cold channel", ":math:`q_{conden}`", "flux_enth_cold", "[t, x]", ":math:`\text{W}\text{/m}^2`"
   "Membrane area", ":math:`A_m`", "area", "None", ":math:`\text{m}^2`"
   "Recovery rate", ":math:`R`", "recovery_mass", "[t]", ":math:`\text{dimensionless}`"

**Additional Variables for VMD**:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Vapor expansion heat flux", ":math:`q_{exp}`", "flux_expansion_heat", "[t, x]", ":math:`\text{W}\text{/m}^2`"

**Additional Variables for GMD**:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Gap thermal conductivity", ":math:`k_{gap}`", "gap_thermal_conductivity", "None", ":math:`\text{W/K/m}`"
   "Gap thickness", ":math:`\sigma_{gap}`", "gap_thickness", "None", ":math:`\text{m}`"
   "gap conduction heat flux", ":math:`q_{gap}`", "flux_conduction_heat_gap", "[t, x]", ":math:`\text{W}\text{/m}^2`"

The following variables are only built when specific configuration key-value pairs are selected.

if ``has_pressure_change`` is set to ``True``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Pressure drop", ":math:`ΔP`", "deltaP", "[t]", ":math:`\text{Pa}`"

if ``temperature_polarization_type`` is set to ``TemperaturePolarizationType.fixed``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Hot channel Convective heat transfer coefficient", ":math:`h_{conv,h}`", "hot_ch.h_conv", "[t]", ":math:`\text{W/K}\text{/m}^2`"
   "Cold channel Convective heat transfer coefficient", ":math:`h_{conv,c}`", "hot_ch.h_conv", "[t]", ":math:`\text{W/K}\text{/m}^2`"

if ``temperature_polarization_type`` is set to ``TemperaturePolarizationType.calculated``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Hot channel Prandtl number", ":math:`Pr_h`", "hot_ch.N_Pr", "[t]", ":math:`\text{dimensionless}`"
   "Cold channel Prandtl number", ":math:`Pr_c`", "cold_ch.N_Pr", "[t]", ":math:`\text{dimensionless}`"
   "Hot channel Nusselt number", ":math:`Nu_h`", "hot_ch.N_Nu", "[t]", ":math:`\text{dimensionless}`"
   "Cold channel Nusselt number", ":math:`Nu_c`", "cold_ch.N_Nu", "[t]", ":math:`\text{dimensionless}`"


if ``concentration_polarization_type`` is set to ``ConcentrationPolarizationType.fixed``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Concentration polarization modulus in hot channel", ":math:`CP_{mod,h}`", "hot_ch.cp_modulus", "[t, j]", ":math:`\text{dimensionless}`"


if ``concentration_polarization_type`` is set to ``ConcentrationPolarizationType.calculated``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Mass transfer coefficient in hot channel", ":math:`k_h`", "hot_ch.K", "[t, x, j]", ":math:`\text{m/s}`"

if ``temperature_polarization_type`` is set to ``TemperaturePolarizationType.calculated``:
or ``mass_transfer_coefficient`` is set to ``MassTransferCoefficient.calculated``
or ``pressure_change_type`` is set to ``PressureChangeType.calculated``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Hot channel height", ":math:`h_{ch,h}`", "hot_ch.channel_height", "None", ":math:`\text{m}`"
   "Hot channel Hydraulic diameter", ":math:`d_{h,h}`", "cold_ch.dh", "None", ":math:`\text{m}`"
   "Hot channel Spacer porosity", ":math:`\epsilon_{sp,h}`", "hot_ch.spacer_porosity", "None", ":math:`\text{dimensionless}`"
   "Hot channel Reynolds number", ":math:`Re_{h}`", "hot_ch.N_Re", "[t, x]", ":math:`\text{dimensionless}`"
   "Cold channel height", ":math:`h_{ch,c}`", "cold_ch.channel_height", "None", ":math:`\text{m}`"
   "Cold channel Hydraulic diameter", ":math:`d_{h,c}`", "cold_ch.dh", "None", ":math:`\text{m}`"
   "Cold channel Spacer porosity", ":math:`\epsilon_{sp,c}`", "cold_ch.spacer_porosity", "None", ":math:`\text{dimensionless}`"
   "Cold channel Reynolds number", ":math:`Re_{c}`", "cold_ch.N_Re", "[t, x]", ":math:`\text{dimensionless}`"


if ``mass_transfer_coefficient`` is set to ``MassTransferCoefficient.calculated``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Schmidt number", ":math:`Sc_h`", "hot_ch.N_Sc", "[t, x]", ":math:`\text{dimensionless}`"
   "Sherwood number", ":math:`Sh_h`", "hot_ch.N_Sh", "[t, x]", ":math:`\text{dimensionless}`"
   "Schmidt number", ":math:`Sc_c`", "cold_ch.N_Sc", "[t, x]", ":math:`\text{dimensionless}`"
   "Sherwood number", ":math:`Sh_c`", "cold_ch.N_Sh", "[t, x]", ":math:`\text{dimensionless}`"

if ``temperature_polarization_type`` is set to ``TemperaturePolarizationType.calculated``:
or ``mass_transfer_coefficient`` is set to ``MassTransferCoefficient.calculated``
or ``pressure_change_type`` is **NOT** set to ``PressureChangeType.fixed_per_stage``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Membrane length", ":math:`L`", "length", "None", ":math:`\text{m}`"
   "Membrane width", ":math:`W`", "width", "None", ":math:`\text{m}`"

if ``pressure_change_type`` is set to ``PressureChangeType.fixed_per_unit_length``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Average pressure drop per unit length of hot channel", ":math:`(\frac{ΔP}{Δx})_{avg,h}`", "hot_ch.dP_dx", "[t]", ":math:`\text{Pa/m}`"
   "Average pressure drop per unit length of cold channel", ":math:`(\frac{ΔP}{Δx})_{avg,c}`", "cold_ch.dP_dx", "[t]", ":math:`\text{Pa/m}`"

if ``pressure_change_type`` is set to ``PressureChangeType.calculated``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Hot channel velocity", ":math:`v_h`", "hot_ch.velocity", "[t, x]", ":math:`\text{m/s}`"
   "Hot channel Friction factor", ":math:`f_h`", "hot_ch.friction_factor_darcy", "[t, x]", ":math:`\text{dimensionless}`"
   "Pressure drop per unit length of hot channel at inlet/outlet", ":math:`(ΔP/Δx)_h`", "hot_ch.dP_dx", "[t, x]", ":math:`\text{Pa/m}`"
   "Cold channel velocity", ":math:`v_c`", "cold_ch.velocity", "[t, x]", ":math:`\text{m/s}`"
   "Pressure drop per unit length of cold channel at inlet/outlet", ":math:`(ΔP/Δx)_c`", "cold_ch.dP_dx", "[t, x]", ":math:`\text{Pa/m}`"

.. _0MD_equations:

Equations
---------

if ``MD_configuration_type`` is set to ``MDconfigurationType.DCMD``:

.. csv-table:: 
   :header: "Description", "Equation"

   "Vapor flux across membrane", ":math:`J(t, x) = \frac{B_0(t)}{\sigma} \times \left( P_{\text{sat, hot}}(t, x) - P_{\text{sat, cold}}(t, x) \right)`"
   "Condensation heat flux to cold channel", ":math:`q_{\text{conden}}(t, x) = J(t, x) \times \widehat{H}_{\text{c}}(t, x, Vap)`"
   "Average condensation flux to cold channel", ":math:`\overline{q}_{\text{conden}}(t) = \frac{1}{2} \sum_{x} q_{\text{conden}}(t, x)`"
   "Cold channel convective heat transfer", ":math:`h_{\text{conv}, c}(t, x) \left( T_{\text{interface}, c}(t, x) - T_{\text{bulk}, c}(t, x) \right) = q_{\text{cond}}(t, x) + q_{\text{conden}}(t, x) - J(t, x) \cdot \widehat{H}_{\text{bulk}, c}(t, x, Liq)`"
   "Hot channel convective heat transfer", ":math:`h_{\text{conv}, h}(t, x) \left( T_{\text{bulk}, h}(t, x) - T_{\text{interface}, h}(t, x) \right) = q_{\text{cond}}(t, x) + q_{\text{evap}}(t, x) - J(t, x) \cdot \widehat{H}_{\text{bulk, h}}(t, x, Liq)`"
   "Conduction heat flux across membrane", ":math:`q_{\text{cond}}(t, x) = \frac{k_{\text{m}}}{\sigma} \left( T_{\text{interface}, h}(t, x) - T_{\text{interface}, c}(t, x) \right)`"
   "Mass transfer balance between hot and cold channel", ":math:`\dot{m}_{\text{cold}}(t, x, p, j) = -\dot{m}_{\text{hot}}(t, x, p, j)`"
   "Conductive heat transfer to cold channel", ":math:`q_{\text{cond, hot}}(t, x) = -q_{\text{cond, cold}}(t, x)`"

if ``MD_configuration_type`` is set to ``MDconfigurationType.VMD``:

.. csv-table:: 
   :header: "Description", "Equation"

   "Vapor flux across membrane", ":math:`J(t, x) = \frac{B_0(t)}{\sigma} \times \left( P_{\text{sat, hot}}(t, x) - P_{\text{vaccuum, cold}}(t, x) \right)`"
   "Hot channel convective heat transfer", ":math:`h_{\text{conv}, h}(t, x) \left( T_{\text{bulk}, h}(t, x) - T_{\text{interface}, h}(t, x) \right) = q_{\text{exp}}(t, x) + q_{\text{evap}}(t, x) - J(t, x) \cdot \widehat{H}_{\text{bulk, h}}(t, x, Liq)`"
   "Vapor expansion heat flux", ":math:`q_{\text{exp}}(t, x) = \frac{R \cdot T}{M} \ln\left( \frac{P_f}{P_p} \right) \cdot J(t, x)`"
   "Mass transfer from vapor phase to vacuum channel", ":math:`\dot{m}_{\text{cold}}(t, x, Vap, j) = -\dot{m}_{\text{hot}}(t, x, Liq, j)`"
   "Conductive heat transfer to cold channel", ":math:`q_{\text{cond, hot}}(t, x) = -q_{\text{exp}}(t, x)`"
   "Cold channel inlet temperature", ":math:`T_{\text{cold, in}}(t) = T_{\text{hot, in}}(t)`"

if ``MD_configuration_type`` is set to ``MDconfigurationType.GMD``:

.. csv-table:: 
   :header: "Description", "Equation"

   "Vapor flux across membrane", ":math:`J(t, x) = \frac{B_0(t)}{\sigma} \times \left( P_{\text{sat, hot}}(t, x) - P_{\text{sat, gap}}(t, x) \right)`"
   "Cold channel convective heat transfer", ":math:`h_{\text{conv}, c}(t, x) \left( T_{\text{interface}, c}(t, x) - T_{\text{bulk}, c}(t, x) \right) = q_{\text{gap}}(t, x)`"
   "Hot channel convective heat transfer", ":math:`h_{\text{conv}, h}(t, x) \left( T_{\text{bulk}, h}(t, x) - T_{\text{interface}, h}(t, x) \right) = q_{\text{cond}}(t, x) + q_{\text{evap}}(t, x) - J(t, x) \cdot \widehat{H}_{\text{bulk, h}}(t, x, Liq)`"
   "Conduction heat flux across membrane", ":math:`q_{\text{cond}}(t, x) = \frac{k_{\text{m}}}{\sigma} \left( T_{\text{interface}, h}(t, x) - T_{\text{interface}, gap}(t, x) \right)`"
   "Conduction heat flux across gap", ":math:`q_{\text{gap}}(t, x) = \frac{k_{\text{m}}}{\sigma} \left( T_{\text{interface}, gap}(t, x) - T_{\text{interface}, c}(t, x) \right)`"
   "Mass transfer balance between hot and gap channel", ":math:`\dot{m}_{\text{gap}}(t, x, Liq, H2O) = -\dot{m}_{\text{hot}}(t, x, Liq, H2O)`"
   "Conductive heat transfer between channels", ":math:`q_{\text{cold}}(t, x) = -q_{\text{hot}}(t, x) - ΔH_{\text{hot}}(t, x) - ΔH_{\text{gap}}(t, x)`"

Common in all configurations:

.. csv-table:: 
   :header: "Description", "Equation"

   "Average flux across membrane", ":math:`J_{avg, j} = \frac{1}{2}\sum_{x} J_{x, j}`"
   "Hot channel membrane-interface solute concentration", ":math:`C_{\text{interface, j, h}}(t, x) = C_{\text{bulk, j, h}}(t, x) \times \exp\left( \frac{J(t, x)}{\rho_{\text{solvent}} \times k_h(t, x, j)} \right)`"
   "Evaporation heat flux from hot channel", ":math:`q_{\text{evap}}(t, x) = J(t, x) \times \widehat{H}_{\text{h}}(t, x, Vap)`"
   "Average evaporation flux from hot channel", ":math:`\overline{q}_{\text{evap}}(t) = \frac{1}{2} \sum_{x} q_{\text{evap}}(t, x)`"
   "Convective heat transfer coefficient", ":math:`h_{\text{conv},(t, x)} = \frac{\kappa_{(t, x)} \cdot \text{Nu}_{(t, x)}}{d_h}`"
   "Effectiveness", ":math:`\epsilon(t) = \frac{T_{\text{cold, first}}(t) - T_{\text{c, last}}(t)}{T_{\text{h, first}}(t) - T_{\text{c, last}}(t)}`"
   "Thermal efficiency", ":math:`\eta(t) = \frac{q_{\text{evap,total}}(t)}{q_{\text{evap,total}}(t) + q_{\text{cond,total}}(t)}`"
   "Concentration polarization modulus",":math:`CP_{mod} = C_{interface}/C_{bulk}`"
   "Mass transfer coefficient",":math:`k_h = \frac{D Sh}{d_h}`"
   "Sherwood number",":math:`Sh[t, x] == 0.2 * (Re[t, x] ** 0.57) * (Pr[t, x] ** 0.4)`"
   "Schmidt number",":math:`Sc = \frac{\mu}{\rho D}`"
   "Reynolds number",":math:`Re = \frac{\rho v_f d_h}{\mu}`"
   "Hydraulic diameter",":math:`d_h = \frac{4\epsilon_{sp}}{2/h_{ch} + (1-\epsilon_{sp})8/h_{ch}}`"
   "Cross-sectional area",":math:`A_c = h_{ch}W\epsilon_{sp}`"
   "Membrane area",":math:`A_m = LW`"
   "Pressure drop",":math:`ΔP = (\frac{ΔP}{Δx})_{avg}L`"
   "Hot channel velocity",":math:`v_h = Q_h/A_c`"
   "Friction factor",":math:`f = 0.42+\frac{189.3}{Re}`"
   "Pressure drop per unit length",":math:`\frac{ΔP}{Δx} = \frac{1}{2d_h}f\rho v_h^{2}`"
   "Recovery rate",":math:`R = \frac{M_{p}}{M_{h,in}}`"
  


Class Documentation
-------------------

* :mod:`watertap.unit_models.MD.membrane_distillation_0D`
* :mod:`watertap.unit_models.MD.membrane_distillation_base`
* :mod:`watertap.unit_models.MD.MD_channel_0D`
* :mod:`watertap.unit_models.MD.MD_channel_base`