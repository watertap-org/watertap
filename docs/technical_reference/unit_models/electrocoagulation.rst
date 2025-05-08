.. _EC_0D:

Electrocoagulation (0D)
=======================

.. index::
   pair: watertap.unit_models.electrocoagulation;electrocoagulation

The main assumptions of the implemented model are (partially adopted from Dubrawski, et al. 2014):

1) Model dimensionality is limited to a 0D control volume
2) Single liquid phase only
3) Steady state only
4) Single solvent (water) only
5) Plug flow only
6) The system is insulated and adiabatic
7) No passivation on electrode surfaces
8) Negligible internal circuit resistance
9) Stoichiometric electrochemical reactions occur at cathode and anode
10) Parallel plate electrodes
11) Each electrode is the same material and same size

Introduction
------------

Electrocoagulation (EC) is a water treatment process that uses electrical current to 
destabilize and aggregate suspended particles in water. The process involves the generation 
of coagulant species via electrochemical reactions at the electrodes. As the cathode is oxidized,
metal ions are released in to the water matrix and form hydroxide species, which then interact 
with and enmesh suspended particles. After formation and agglomeration, the flocculated material 
is either settled (or floated) out of the treated water stream.

EC is an electrochemically complex process. The performance and technoeconomics is influenced by 
many factors including the composition of the water matrix, the applied current density, the electrode material, 
and other aspects of the reactor design. Because energy consumption can be a significant component
of the overall cost of the process, this model presents three different approaches to estimate the
overpotentials associated with the electrochemical reactions.

Model Configurations
---------------------

The EC model includes different configuration options for the electrode material, reactor material, and the overpotential calculation:

- Electrode material: aluminum (default) and iron. 
- Overpotential calculation: fixed (default), regression approximation, and a detailed calculation.
- Reactor material: carbon steel, stainless steel (default) and PVC.

Selecting either electrode material will properly set the WaterTAP model parameters for electrode density 
(``density_electrode_material``), molecular weight (``mw_electrode_material``), charge transfer number
(``charge_transfer_number``), and the stoichiometric coefficient of the electrochemical reaction (``stoich_coeff``).
Any of these parameters can be changed by the user after the model build. 

If the user uses the defaut overpotential calculation, the overpotential variable (``overpotential``) is a degree of freedom and must be fixed.

If the user selects the regression approximation for the overpotential calculation, the model will set default values
for the overpotential regression coefficients (``overpotential_k1`` and ``overpotential_k2``).

If the user also selects the detailed calculation for the overpotential calculation, the model will also set default values 
for the anode cell potential (``anode_cell_potential_std``), the anode entropy change (``anode_entropy_change_std``), 
the anodic exchange current density (``anodic_exchange_current_density``), and the cathodic exchange current density 
(``cathodic_exchange_current_density``) relevant to each electrode material.

Overpotential Calculation
^^^^^^^^^^^^^^^^^^^^^^^^^

The overpotential is the additional voltage required over the ohmic potential to drive the electrochemical reaction.
The electrocoagulation model determines the total cell voltage required to drive the electrochemical reactions according to:

.. math::
    E_{cell} = E_{ohmic} + E_{over}

Where :math:`E_{cell}` is the total cell voltage, :math:`E_{ohmic}` is the ohmic potential required, and 
:math:`E_{over}` is the overpotential. 

The WaterTAP electrocoagulation model provides three options for calculating the overpotential, outlined as follows.

Fixed
++++++

If this overpotential calculation is selected, the user must provide a fixed value for the overpotential variable 
(``overpotential``) in volts. This value is used directly in the calculation of the total cell voltage.

Regression Approximation
++++++++++++++++++++++++

This overpotential calculation uses a regression adapted from Eq. 18 in Gu et al. (2009) to determine the overpotential:

.. math::
    E_{over} = k_1 \text{ln}\left(i \right) + k_2

Where :math:`E_{over}` is the overpotential, :math:`i` is the current density (mA/cm\ :superscript:`2`),
:math:`k_1` is a regression coefficient (mV), and :math:`k_2` is a regression coefficient (mV). The values for the regression coefficients 
have default values, but users should adjust them based on experimental data for their specific system.

Detailed Calculation
++++++++++++++++++++

If the detailed calculation is selected, the model will use the Nernst equation to calculate the overpotential based on the standard cell potential,
the entropy change, and the exchange current densities for the anodic and cathodic reactions, and a Tafel slope parameter to estimate the
activation overpotential. 
In general, the overpotential is calculated as follows:

.. math::
    E_{over} = |E_c - E_a| + \varphi_a + |\varphi_c| + \psi_a + |\psi_c|

Where :math:`E_c` is the non-equilibrium electrode potential at the cathode, :math:`E_a` is the non-equilibrium electrode potential at the anode,
:math:`\varphi_a` is the anodic activation overpotential, :math:`\varphi_c` is the cathodic activation overpotential,
:math:`\psi_a` is the anodic concentration overpotential, and :math:`\psi_c` is the cathodic concentration overpotential.
The electrocoagulation model assumes the concentration overpotential is negligible 
(i.e., that the electrochemical reactions are not mass transfer limited) and :math:`\psi_c = \psi_a = 0`.

The non-equilibrium electrode potentials at the cathode and anode are calculated via the Nernst equation:

.. math::
    E_a = E_{a}^0 + \frac{\Delta S_a (T - T_0)}{z_a F} + \frac{RT}{z_a F} \text{ln}\left( C_{i} \right)


.. math::
    E_c = E_{c}^0 + \frac{\Delta S_c (T - T_0)}{z_c F} + \frac{RT}{z_cF} \text{ln}\left( p_{H_2} \left( C_{OH}\right)^2 \right)


Where :math:`E_{i}^0` is the standard cell potential, :math:`R` is the universal gas constant (8.314 J/(mol K)),
:math:`T` is the temperature (K), :math:`z_i` is the number of electrons transferred in the electrochemical reaction, 
:math:`F` is the Faraday constant (96,485 C/mol), :math:`\Delta S_i` is the entropy change for the reaction (J/(mol K)),
:math:`C_{i}` is the concentration of the reactant species (mol/L), :math:`C_{OH}` is the hydroxide concentration (mol/L), 
and :math:`p_{H_2}` is the partial pressure of hydrogen gas (atm).

The anodic and cathodic activation overpotentials are calculated using the Tafel equation:

.. math::
    \varphi_a = b_a \text{ln}\left( \frac{i}{i_{a0}} \right)

.. math::
    \varphi_c = b_c \text{ln}\left( \frac{i}{i_{c0}} \right)

Where :math:`i_{a0}` and :math:`i_{c0}` are the anodic and cathodic exchange current densities (A/m\ :superscript:`2`),
:math:`b_a` and :math:`b_c` are the anodic and cathodic Tafel slope parameters (V), and :math:`i` is the current density (A/m\ :superscript:`2`). 


Ports
-----

The model provides three ports (Pyomo notation in parenthesis):

* Inlet port (``inlet``)
* Outlet port (``outlet``)
* Byproduct port (``byproduct``)


Sets
----

The table below outlines example Sets that could be used with the electrocoagulation model.
Any component can be included as long as it is properly configured into the property package.

.. csv-table::
   :header: "Description", "Symbol", "Example Indices"

   "Time", ":math:`t`", "``[0]``"
   "Phases", ":math:`p`", "``['Liq']``"
   "Components", ":math:`j`", "``['H2O', 'Cation_+', 'Anion_-', 'Inert']``"


.. _EC_variables:

Model Components
-----------------

The electrocoagulation model includes variables, parameters, and expressions that are common to 
all configurations. These are provided in the table below.

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"
   
   **Variables**
   "Inlet temperature", ":math:`T`", "``temperature``", "``[t]``", ":math:`\text{K}`"
   "Inlet pressure", ":math:`p`", "``pressure``", "``[t]``", ":math:`\text{Pa}`"
   "Component mass flow rate", ":math:`M_j`", "``flow_mass_phase_comp``", "``[t, p, j]``", ":math:`\text{kg s}^{-1}`"
   "Phase volumetric flow rate", ":math:`q_j`", "``flow_vol_phase``", "``[t, p]``", ":math:`\text{m}^{3} \text{ s}^{-1}`"
   "Coagulant dose", ":math:`D_c`", "``coagulant_dose``", None, ":math:`\text{g L}^{-1}`"
   "Electrode thickness", ":math:`d_{electrode}`", "``electrode_thickness``", None, ":math:`\text{m}`"
   "Electrode mass", ":math:`m_{electrode}`", "``electrode_mass``", None, ":math:`\text{kg}`"
   "Electrode volume", ":math:`V_{electrode}`", "``electrode_volume``", None, ":math:`\text{m}^3`"
   "Electrode gap", ":math:`d_{gap}`", "``electrode_gap``", None, ":math:`\text{m}`"
   "Electrolysis time", ":math:`t_{elec}`", "``electrolysis_time``", None, ":math:`\text{min}`"
   "Current density", ":math:`i`", "``current_density``", None, ":math:`\text{A m}^{-2}`"
   "Applied current", ":math:`I`", "``applied_current``", None, ":math:`\text{A}`"
   "Ohmic resistance", ":math:`R_{ohmic}`", "``ohmic_resistance``", None, ":math:`\Omega \text{ m}^{2}`"
   "Charge loading rate", ":math:`CLR`", "``charge_loading_rate``", None, ":math:`\text{C L}^{-1}`"
   "Current efficiency", ":math:`\eta`", "``current_efficiency``", None, ":math:`\text{dimensionless}`"
   "Overpotential", ":math:`E_{over}`", "``overpotential``", None, ":math:`\text{V}`"
   "Cell voltage", ":math:`E_{cell}`", "``cell_voltage``", None, ":math:`\text{V}`"
   "Anode area", ":math:`A_{anode}`", "``anode_are``", None, ":math:`\text{m}^2`"
   "Cathode area", ":math:`A_{cathode}`", "``cathode_area``", None, ":math:`\text{m}^2`"
   "Volume of electrocoagulation reactor", ":math:`V_{r}`", "``cell_volume``", None, ":math:`\text{m}^3`"
   "Total floc basin volume (flotation + sedimentation)", ":math:`V_{floc}`", "``floc_basin_vol``", None, ":math:`\text{m}^3`"
   "Floc basin retention time", ":math:`t_{floc}`", "``floc_retention_time``", None, ":math:`\text{min}`"

   **Parameters**
   "Component removal efficiency on mass basis", ":math:`\eta_{j}`", "``removal_frac_mass_comp``", ``[j]``, ":math:`\text{dimensionless}`"
   "Water recovery on mass basis", ":math:`\eta_{w}`", "``recovery_frac_mass_water``", None, ":math:`\text{dimensionless}`"
   "Conversion factor for mg/L TDS to S/m", ":math:`x`", "``tds_to_cond_conversion``", None, ":math:`\text{mg m }\text{L}^{-1}\text{ S}^{-1}`"
   "Standard temperature", ":math:`T_0`", "``standard_temperature``", None, ":math:`\text{K}`"
   "Electrode molecular weight", ":math:`MW`", "``mw_electrode_material``", None, ":math:`\text{kg mol}^{-1}`"
   "Stoichiometric coefficient for electrode material", ":math:`\nu`", "``stoich_coeff``", None, ":math:`\text{dimensionless}`"
   "Charge transfer number", ":math:`z`", "``charge_transfer_number``", None, ":math:`\text{dimensionless}`"
   "Electrode density", ":math:`\rho_{electrode}`", "``density_electrode_material``", None, ":math:`\text{kg m}^{-3}`"
   "Fractional increase in water temperature from inlet to outlet", ":math:`x_T`", "``frac_increase_temperature``", None, ":math:`\text{dimensionless}`"

   **Expressions**
   "Conductivity", ":math:`\kappa`", "``conductivity``", None, ":math:`\text{S m}^{-1}`"
   "Electrode area total", ":math:`A_{electrode}`", "``electrode_area_total``", None, ":math:`\text{m}^2`"
   "Total power required", ":math:`P_{tot}`", "``power_required``", None, ":math:`\text{W}`"
   "Power density Faradaic", ":math:`p_{F}`", "``power_density_faradaic``", None, ":math:`\mu\text{W m}^{-2}`"
   "Power density total", ":math:`p_{total}`", "``power_density_total``", None, ":math:`\mu\text{W m}^{-2}`"


If ``overpotential_calculation`` is set to ``regression``, the following variables are also created:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   **Variables**
   "Overpotential regression coefficient 1", ":math:`k_1`", "``overpotential_k1``", None, ":math:`\text{mV}`"
   "Overpotential regression coefficient 2", ":math:`k_2`", "``overpotential_k2``", None, ":math:`\text{mV}`"

If ``overpotential_calculation`` is set to ``detailed``, the following variables, parameters, and expressions are also created:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units", "Default Value"

   **Variables**
   "Anodic Tafel slope", ":math:`b_a`", "``tafel_slope_anode``", None, ":math:`\text{V}`"
   "Cathodic Tafel slope", ":math:`b_c`", "``tafel_slope_cathode``", None, ":math:`\text{V}`"

   **Parameters**
   "Anodic non-equilibrium cell potential, standard @ 25C", ":math:`E_{a}^0`", "``anode_cell_potential_std``", None, ":math:`\text{V}`", -0.5
   "Anodic entropy change", ":math:`\frac{\Delta S_a}{z_aF}`", "``anode_entropy_change_std``", None, ":math:`\text{V K}^{-1}`", 1e-4
   "Anodic exchange current density", ":math:`i_{a0}`", "``anodic_exchange_current_density``", None, ":math:`\text{A m}^{-2}`", 2e-5
   "Cathodic non-equilibrium cell potential, standard @ 25C", ":math:`E_{c}^0`", "``cathode_cell_potential_std``", None, ":math:`\text{V}`", -0.83
   "Cathode entropy change", ":math:`\frac{\Delta S_c}{z_cF}`", "``cathode_entropy_change_std``", None, ":math:`\text{V K}^{-1}`", -0.000836
   "Cathode surface pH", ":math:`pH`", "``cathode_surface_pH``", None, ":math:`\text{dimensionless}`", 11

   **Expressions**
   "Anode cell potential via Nernst equation", ":math:`E_a`", "``anode_cell_potential``", None, ":math:`\text{V}`"
   "Cathodic cell potential via Nernst equation", ":math:`E_c`", "``cathode_cell_potential``", None, ":math:`\text{V}`"
   "Anodic activation overpotential", ":math:`\varphi_a`", "``anode_overpotential``", None, ":math:`\text{V}`"
   "Cathodic activation overpotential", ":math:`\varphi_c`", "``cathode_overpotential``", None, ":math:`\text{V}`"

Degrees of Freedom
--------------------

Aside from the inlet feed state variables (temperature, pressure, component molar flowrate),
the user must specify 8-9 degrees of freedom to fully specify the model, depending on the configuration.

The following degrees of freedom should be specified regardless of the configuration:

- ``electrode_thickness``
- ``electrode_gap``
- ``electrolysis_time``
- ``floc_retention_time``

The following degrees of freedom are fixed dependent on the configuration:

- ``overpotential`` (if ``overpotential_calculation`` is set to ``fixed``)
- ``overpotential_k1`` and ``overpotential_k2`` (if ``overpotential_calculation`` is set to ``regression``)
- ``tafel_slope_anode`` and ``tafel_slope_cathode`` (if ``overpotential_calculation`` is set to ``detailed``)

Then, the user can select combinations of three of the following variables to have a fully specified model.
The specific combination would be dependent on what the user knows about the system and their modeling objectives.

- ``current_density``
- ``applied_current``
- ``current_efficiency``
- ``cell_voltage``
- ``coagulant_dose``
- ``charge_loading_rate``
- ``anode_area`` or ``cathode_area`` 

Solution Component Information
------------------------------
The electrocoagulation model is designed to work with WaterTAP's
multi-component aqueoous solution (MCAS) property package.
The inlet solute list must contain ``TDS`` because the model 
uses the TDS concentration to calculate the conductivity of the solution.
Because the removal efficiency is defined on a mass basis, MCAS must 
be configured to use mass as the material flow basis.

An example configuration is provided below:

.. code-block::

    ec_feed = {
        "solute_list": ["TDS", "Ca_2+", "Mg_2+"],
        "mw_data": {
            "TDS": 58.44e-3,
            "Ca_2+": 40.08e-3,
            "Mg_2+": 24.31e-3,
        },
        "material_flow_basis": MaterialFlowBasis.mass,
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(**ec_feed)
    m.fs.unit = Electrocoagulation(
        property_package=m.fs.properties,
        electrode_material="iron",
        overpotential_calculation="detailed",
    )

Equations and Relationships
---------------------------

.. csv-table::
    :header: "Description", "Equation"
    
    **Common** 
    "Conductivity", ":math:`\kappa = C_{TDS} / x`"
    "Total electrode area", ":math:`A_{electrode} = A_{anode} + A_{cathode}`"
    "Power required", ":math:`P_{tot} = E_{cell} I`"
    "Power density Faradaic", ":math:`p_{F} = \frac{E_{over}I}{A_{anode}}`"
    "Power density total", ":math:`p_{tot} = \frac{P_{tot}}{A_{anode}}`"
    "Effluent temperature", ":math:`T_{out} = x_T T_{in}`"
    "Water recovery", ":math:`M_{H_2O, out} = M_{H_2O, in} \eta_w`"
    "Water mass balance", ":math:`M_{H_2O, out} = M_{H_2O, in} - M_{H_2O, byprod}`"
    "Component mass balance", ":math:`M_{j, out} = M_{j, in} - M_{j, byprod}`"
    "Component removal efficiency", ":math:`M_{j, byprod} = \eta_j M_{j, in}`"
    "Charge loading rate", ":math:`CLR = \frac{I}{q_{liq}}`"
    "Floc reactor volume", ":math:`V_{floc} = q_{liq} t_{floc}`"
    "Faraday's Law", ":math:`D_c = \frac{I \eta MW}{q_{liq} z F}`"
    "Anode area required", ":math:`A_{anode} = \frac{I}{i}`"
    "Cathode area required", ":math:`A_{cathode} = A_{anode}`"
    "Cell voltage required", ":math:`E_{cell} = E_{over} + \frac{I R_{ohmic}}{A_{anode}}`"
    "Electrode volume", ":math:`V_{electrode} = \left( A_{anode} + A_{cathode} \right)  d_{electrode}`"
    "Electrode mass", ":math:`m_{electrode} = V_{electrode} \rho_{electrode}`"
    "Reactor volume", ":math:`V_{cell} = q_{liq} t_{elec}`"
    "Ohmic resistance", ":math:`R_{ohmic} = \frac{d_{gap}}{\kappa}`"

    **Regression**
    "Overpotential regression", ":math:`E_{over} = k_1 \text{ln}(i) + k_2`"

    **Detailed**
    "Anodic cell potential", ":math:`E_a = E_{a}^0 + \frac{\Delta S_a (T - T_0)}{z_a F} + \frac{RT}{z_a F} \text{ln}(C_{i})`"
    "Cathodic cell potential", ":math:`E_c = E_{c}^0 + \frac{\Delta S_c (T - T_0)}{z_c F} + \frac{RT}{z_cF} \text{ln}(p_{H_2} (C_{OH})^2)`"
    "Anodic activation overpotential", ":math:`\varphi_a = b_a \text{ln}(i / i_{a0})`"
    "Cathodic activation overpotential", ":math:`\varphi_c = b_c \text{ln}(i / i_{c0})`"
    "Overpotential", ":math:`E_{over} = |E_c - E_a| + \varphi_a + |\varphi_c|`"


References
----------

| K. L. Dubrawski, C. Du and M. Mohseni (2014)
| General Potential-Current Model and Validation for Electrocoagulation
| Electrochimica Acta 2014 Vol. 129 Pages 187-195
| DOI: 10.1016/j.electacta.2014.02.089

| Z. Gu, Z. Liao, M. Schulz, J. R. Davis, J. C. Baygents and J. Farrell (2009)
| Estimating Dosing Rates and Energy Consumption for Electrocoagulation Using Iron and Aluminum Electrodes
| Industrial & Engineering Chemistry Research 2009 Vol. 48 Issue 6 Pages 3112-3117
| DOI: 10.1021/ie801086c

| Bratsch, S. G. (1989). 
| Standard Electrode Potentials and Temperature Coefficients in Water at 298.15 K. 
| Journal of Physical and Chemical Reference Data, 18(1), 1-21. 
| DOI: 10.1063/1.555839 

| Zhang, F., Yang, C., Zhu, H., Li, Y., & Gui, W. (2020). 
| An integrated prediction model of heavy metal ion concentration for iron electrocoagulation process. 
| Chemical Engineering Journal, 391, 123628. 
| DOI: 10.1016/j.cej.2019.123628 