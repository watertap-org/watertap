Ion Exchange (0D)
=================

.. index::
   pair: watertap.unit_models.ion_exchange_0D;ion_exchange_0D

The main assumptions of the implemented model are as follows:

1) Model dimensionality is limited to a 0D control volume
2) Single liquid phase only
3) Steady state only
4) Single solute and single solvent (water) only
5) Plug flow conditions
6) Isothermal conditions
7) Favorable Langmuir *or* Freundlich isotherm

Introduction
------------

Ion exchange is the reversible transfer of one or more solutes between a fluid phase and a sorbent.
This process is becoming increasingly popular in drinking water treatment applications where it is
used for water softening and demineralization. This implementation of the fixed-bed ion exchange model
accounts for process equilibrium, kinetics, and hydrodynamics to predict performance, bed and column geometry, and capital/operating costs.
The ion exchange process operates as a cycle with four steps:

(1) Service
(2) Backwashing
(3) Regeneration
(4) Rinsing

Critical to predicting performance of an ion exchange process is having an estimate for the breakthrough time,
or the duration of treatment before the solute begins exiting the column at a concentration unacceptable to the operator.
At this time, the mass transfer zone is approaching the end of the ion exchange bed, the resin is nearing exhaustion,
and the regeneration cycle can begin. Fundamental to this model is the assumption that the isotherm between the solute
and the resin is favorable, and thus the mass transfer zone is shallow.

.. note:: 
    If using ``single-use`` configuration for ``regenerant``, the backwashing, regeneration, and rinsing steps are not modeled and all associated costs for these steps are zero.

Isotherm Configurations
^^^^^^^^^^^^^^^^^^^^^^^

The model requires the user has either Langmuir or Freundlich isotherm equilibrium parameters for their specific system.
Variables in the following equations and their corollary in the WaterTAP model are defined in a following section.

Langmuir
++++++++

For the Langmuir isotherm, the Langmuir parameter :math:`La` (``langmuir`` in the WaterTAP model) is:

.. math::
    La = \frac{1}{1 + K C_0}

Where :math:`K` is an equilibrium constant derived from experimental data, and :math:`C_0` is the influent concentration of the target ion. 
(Note: This equation is not included in the model). :math:`La` is used in the dimensionless form of the Langmuir isotherm:

.. math::
    La = \frac{X (1 - Y)}{Y (1 - X)}

:math:`Y` is the ratio of equilibrium to total resin capacity (``resin_eq_capacity`` and ``resin_max_capacity``, respectively in the WaterTAP model).
For a favorable isotherm (a core assumption of the model), :math:`La` is less than one.

Freundlich
++++++++++

For the Freundlich isotherm, the model assumes the user has fit breakthrough data to the Clark model.  
The general solution evaluated at 50% breakthrough is:

.. math::
    \frac{C_b}{C_0} = \frac{1}{\bigg(1 + (2^{n - 1} - 1)\text{exp}\bigg[\frac{k_T Z (n - 1)}{BV_{50} u_{bed}} (BV_{50} - BV)\bigg]\bigg)^{\frac{1}{n-1}}}

The form often fit to breakthrough data is:

.. math::
    \frac{C_b}{C_0} = \frac{1}{A \text{exp}\big[\frac{-r Z}{u_{bed}} BV\big]^{\frac{1}{n-1}}}

The full derivation for both equations is provided in Croll et al. (2023).

Ports
-----

The model provides three ports (Pyomo notation in parenthesis):

* Inlet port (inlet)
* Outlet port (outlet)
* Regeneration port (regen)

Sets
----

The table below outlines example Sets that could be used with the ion exchange model.
"Components" is a subset of "Ions" and uses the same symbol ``j``. 
They can include any ion as long as the ion is configured into the property package.
``target_ion_set`` includes the component to be removed via the ion exchange process. 
The current model implementation is only for a single component, but ``target_ion_set`` is included for future development of a multi-component model.

.. csv-table::
   :header: "Description", "Symbol", "Example Indices"

   "Time", ":math:`t`", "``[0]``"
   "Phases", ":math:`p`", "``['Liq']``"
   "Components", ":math:`j`", "``['H2O', 'Cation_+', 'Anion_-', 'Inert']``"
   "Ions", ":math:`j`", "``['Cation_+', 'Anion_-']``"
   "Target Ion", ":math:`j`", "``['Cation_+']``"

In this example, the influent stream contains ``H2O`` (always included), ``Cation_+``, ``Anion_-``, and an uncharged component ``Inert``. 
The user would specify the concentration of each as part of the property package in the model build.
The charged components are included in "Ions", a subset of "Compoenents". The model is configured as a cation exchange process since ``target_ion_set`` contains a positively
charged component, ``Cation_+``.


.. _IX_variables:

Model Components
----------------

The ion exchange model includes many variables, parameters, and expressions that are common to both the
``langmuir`` and ``freundlich`` isotherm configurations. These are provided in the table below.

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"
   
   **Variables**
   "Inlet temperature", ":math:`T`", "``temperature``", "``[t]``", ":math:`\text{K}`"
   "Inlet pressure", ":math:`p`", "``pressure``", "``[t]``", ":math:`\text{Pa}`"
   "Component molar flow rate", ":math:`N_j`", "``flow_mol_phase_comp``", "``[t, 'Liq', 'H2O']``", ":math:`\text{mol/s}`"
   "Control volume mass transfer term", ":math:`\dot{m}_j`", "``process_flow.mass_transfer_term``", "``[t, 'Liq', j]``", ":math:`\text{mol/s}`"
   "Service flow rate through resin bed in bed volumes per hour", ":math:`SFR`", "``service_flow_rate``", "None", ":math:`\text{hr}^{-1}`"
   "Linear velocity through bed", ":math:`u_{bed}`", "``vel_bed``", "None", ":math:`\text{m/s}`"
   "Interstitial velocity through bed", ":math:`u_{inter}`", "``vel_inter``", "None", ":math:`\text{m/s}`"
   "Number of operational columns", ":math:`n_{op}`", "``number_columns``", "None", ":math:`\text{dimensionless}`"
   "Number of redundant columns", ":math:`n_{red}`", "``number_columns_redund``", "None", ":math:`\text{dimensionless}`"
   "Bed depth", ":math:`Z`", "``bed_depth``", "None", ":math:`\text{m}`"
   "Column height", ":math:`H_{col}`", "``col_height``", "None", ":math:`\text{m}`"
   "Column diameter", ":math:`D_{col}`", "``col_diam``", "None", ":math:`\text{m}`"
   "Column height to diameter ratio", ":math:`R_{HD}`", "``col_height_to_diam_ratio``", "None", ":math:`\text{dimensionless}`"
   "Total bed volume", ":math:`V_{res, tot}`", "``bed_vol_tot``", "None", ":math:`\text{m}^3`"
   "Resin bead diameter", ":math:`d`", "``resin_diam``", "None", ":math:`\text{m}`"
   "Resin bulk density", ":math:`\rho_{b}`", "``resin_bulk_dens``", "None", ":math:`\text{kg/L}`"
   "Resin surface area per volume", ":math:`a_{s}`", "``resin_surf_per_vol``", "None", ":math:`\text{m}^{-1}`"
   "Bed porosity", ":math:`\varepsilon`", "``bed_porosity``", "None", ":math:`\text{dimensionless}`"
   "Number of cycles before regenerant disposal", ":math:`N_{regen}`", "``regen_recycle``", "None", ":math:`\text{dimensionless}`"
   "Relative breakthrough concentration at breakthrough time ", ":math:`X`", "``c_norm``", "``target_ion_set``", ":math:`\text{dimensionless}`"
   "Breakthrough time", ":math:`t_{break}`", "``t_breakthru``", "None", ":math:`\text{s}`"
   "Empty Bed Contact Time (EBCT)", ":math:`EBCT`", "``ebct``", "None", ":math:`\text{s}`"
   "Reynolds number", ":math:`Re`", "``N_Re``", "None", ":math:`\text{dimensionless}`"
   "Schmidt number", ":math:`Sc`", "``N_Sc``", "``target_ion_set``", ":math:`\text{dimensionless}`"
   "Sherwood number", ":math:`Sh`", "``N_Sh``", "``target_ion_set``", ":math:`\text{dimensionless}`"
   "Peclet particle number", ":math:`Pe_{p}`", "``N_Pe_particle``", "None", ":math:`\text{dimensionless}`"
   "Peclet bed number", ":math:`Pe_{bed}`", "``N_Pe_bed``", "None", ":math:`\text{dimensionless}`"
   
   **Parameters**
   "Regeneration time", ":math:`t_{regen}`", "``t_regen``", "None", ":math:`\text{s}`"
   "Backwash time", ":math:`t_{bw}`", "``t_bw``", "None", ":math:`\text{s}`" 
   "Backwash loading rate", ":math:`u_{bw}`", "``bw_rate``", "None", ":math:`\text{m/hr}`" 
   "Number of bed volumes for rinse step", ":math:`N_{rinse}`", "``rinse_bv``", "None", ":math:`\text{dimensionless}`" 
   "Pump efficiency", ":math:`\eta`", "``pump_efficiency``", "None", ":math:`\text{dimensionless}`" 
   "Service-to-regeneration flow ratio", ":math:`R`", "``service_to_regen_flow_ratio``", "None", ":math:`\text{dimensionless}`" 
   "Pressure drop equation intercept", ":math:`p_{drop,A}`", "``p_drop_A``", "None", ":math:`\text{dimensionless}`" 
   "Pressure drop equation B", ":math:`p_{drop,B}`", "``p_drop_B``", "None", ":math:`\text{dimensionless}`" 
   "Pressure drop equation C", ":math:`p_{drop,C}`", "``p_drop_C``", "None", ":math:`\text{dimensionless}`" 
   "Bed expansion fraction equation intercept", ":math:`H_{expan,A}`", "``bed_expansion_frac_A``", "None", ":math:`\text{dimensionless}`" 
   "Bed expansion fraction equation B parameter", ":math:`H_{expan,B}`", "``bed_expansion_frac_B``", "None", ":math:`\text{dimensionless}`" 
   "Bed expansion fraction equation C parameter", ":math:`H_{expan,C}`", "``bed_expansion_frac_C``", "None", ":math:`\text{dimensionless}`" 

    **Expressions**
   "Fraction of bed depth increase during backwashing", ":math:`X_{expan}`", "``bed_expansion_frac``", "None", ":math:`\text{dimensionless}`" 
   "Additional column sidewall height required for bed expansion", ":math:`H_{expan}`", "``bed_expansion_h``", "None", ":math:`\text{dimensionless}`" 
   "Backwashing volumetric flow rate", ":math:`Q_{bw}`", "``bw_flow``", "None", ":math:`\text{m}^{3}\text{/s}`" 
   "Rinse time", ":math:`t_{rinse}`", "``t_rinse``", "None", ":math:`\text{s}`" 
   "Rinse volumetric flow rate", ":math:`Q_{rinse}`", "``rinse_flow``", "None", ":math:`\text{m}^{3}\text{/s}`" 
   "Regen + Rinse + Backwash time", ":math:`t_{waste}`", "``t_waste``", "None", ":math:`\text{s}`" 
   "Cycle time", ":math:`t_{cycle}`", "``t_cycle``", "None", ":math:`\text{s}`" 
   "Bed volume of one unit", ":math:`V_{res}`", "``bed_vol``", "None", ":math:`\text{m}^{3}`"
   "Column volume of one unit", ":math:`V_{col}`", "``col_vol_per``", "None", ":math:`\text{m}^{3}`" 
   "Total column volume", ":math:`V_{col, tot}`", "``col_vol_tot``", "None", ":math:`\text{m}^{3}`" 
   "Bed volumes of throughput at breakthrough", ":math:`BV`", "``bv_calc``", "None", ":math:`\text{dimensionless}`" 
   "Regeneration solution tank volume", ":math:`V_{regen}`", "``regen_tank_vol``", "None", ":math:`\text{m}^{3}`" 
   "Pressure drop through resin bed", ":math:`p_{drop}`", "``pressure_drop``", "None", ":math:`\text{psi}`" 
   "Power of main booster pump", ":math:`P_{main}`", "``main_pump_power``", "None", ":math:`\text{kW}`" 
   "Regen pump power", ":math:`P_{regen}`", "``regen_pump_power``", "None", ":math:`\text{kW}`" 
   "Backwash pump power", ":math:`P_{bw}`", "``bw_pump_power``", "None", ":math:`\text{kW}`" 
   "Rinse pump power", ":math:`P_{rinse}`", "``rinse_pump_power``", "None", ":math:`\text{kW}`" 


If ``isotherm`` is set to ``langmuir``, the model includes the following components:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   **Variables**
   "Langmuir equilibrium parameter for resin/ion system", ":math:`La`", "``langmuir``", "``target_ion_set``", ":math:`\text{dimensionless}`"
   "Maximum resin capacity", ":math:`q_{max}`", "``resin_max_capacity``", "None", ":math:`\text{mol/kg}`"
   "Equilibrium resin capacity", ":math:`q_{eq}`", "``resin_eq_capacity``", "None", ":math:`\text{mol/kg}`"
   "Unused resin capacity", ":math:`q_{un}`", "``resin_unused_capacity``", "None", ":math:`\text{mol/kg}`"
   "Sorbed mass of ion", ":math:`M_{out}`", "``mass_removed``", "``target_ion_set``", ":math:`\text{mol}`"
   "Number of transfer units", ":math:`N`", "``num_transfer_units``", "None", ":math:`\text{dimensionless}`"
   "Dimensionless time", ":math:`\tau`", "``dimensionless_time``", None, ":math:`\text{dimensionless}`"
   "Partition ratio", ":math:`\Lambda`", "``partition_ratio``", "None", ":math:`\text{dimensionless}`"
   "Fluid mass transfer coefficient", ":math:`k_{f}`", "``fluid_mass_transfer_coeff``", "``target_ion_set``", ":math:`\text{m/s}`"
   "Mass removed during service", ":math:`M_{rem,j}`", "``mass_removed``", "``target_ion_set``", ":math:`\text{mol}`"
   


If ``isotherm`` is set to ``freundlich``, the model includes the following components:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   **Variables**
   "Freundlich isotherm exponent for resin/ion system", ":math:`n`", "``freundlich_n``", "None", ":math:`\text{dimensionless}`"
   "Bed volumes at breakthrough", ":math:`BV`", "``bv``", "None", ":math:`\text{dimensionless}`"
   "Bed volumes at 50% influent conc.", ":math:`BV_{50}`", "``bv_50``", "None", ":math:`\text{dimensionless}`"
   "Mass transfer coefficient", ":math:`k_T`", "``mass_transfer_coeff``", "None", ":math:`\text{s}^{-1}`"
   "Average relative breakthrough concentration at breakthrough time", ":math:`X_{avg}`", "``c_norm_avg``", "None", ":math:`\text{dimensionless}`"
   "Relative breakthrough conc. for trapezoids", ":math:`X_{trap,k}`", "``c_traps``", "``k``", ":math:`\text{dimensionless}`"
   "Breakthrough times for trapezoids", ":math:`t_{trap,k}`", "``tb_traps``", "``k``", ":math:`\text{s}`"
   "Area of trapezoids", ":math:`A_{trap,k}`", "``traps``", "``k``", ":math:`\text{dimensionless}`"


Degrees of Freedom
------------------

Aside from the inlet feed state variables (temperature, pressure, component molar flowrate), the user must specify an additional 9 degrees of freedom
for both the ``langmuir`` and ``freundlich`` isotherm model configurations to achieve a fully specified model (i.e., zero degrees of freedom).
Depending on the data available to the user and the objectives of the modeling exercise, different combinations of variables can be fixed to achieve 
zero degrees of freedom.

For either model configuration, the user can fix the following variables:

* ``resin_diam``
* ``resin_bulk_dens``
* ``bed_porosity``
* ``service_flow_rate`` (alternatively, ``vel_bed``)
* ``bed_depth``
* ``number_columns``


Langmuir DOF 
^^^^^^^^^^^^

If ``isotherm`` is set to ``langmuir``, the additional variables to fix are:

* ``langmuir`` 
* ``resin_max_capacity``
* ``dimensionless_time`` (can be fixed to default value of 1)


Freundlich DOF
^^^^^^^^^^^^^^

If ``isotherm`` is set to ``freundlich``, the additional variables to fix are:

* ``freundlich_n``
* ``bv`` 
* ``c_norm``
* one of ``bv_50`` or ``mass_transfer_coeff`` as determined from Clark model equations



Solution Component Information
------------------------------
The IonExchange0D model is designed to work with WaterTAP's 
Multi-component aqueous solution (MCAS) property package. 
In addition to providing a list of solute ions, users must 
provide parameter information for each ion including molecular weight,
diffusivity data, and charge data. An example of how this 
data is used to build a model is provided below.

.. code-block::

    target_ion = "Ca_2+"
    ion_props = {
        "solute_list": [target_ion],
        "diffusivity_data": {("Liq", target_ion): 9.2e-10},
        "mw_data": {"H2O": 0.018, target_ion: 0.04},
        "charge": {target_ion: 2},
    }
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(**ion_props)
    ix_config = {
        "property_package": m.fs.properties,
        "target_ion": target_ion,
    }
    m.fs.ix = IonExchange0D(**ix_config)


.. .. code-block::

Equations and Relationships
---------------------------

.. csv-table::
   :header: "Description", "Equation"

    **Common**
   "Service flow rate", ":math:`SFR = \frac{Q_{p, in}}{V_{res, tot}}`"
   "Total bed volume", ":math:`V_{res, tot} = V_{bed}n_{op}`"
   "Flow through bed constraint", ":math:`\frac{Z}{u_{bed}} = \frac{V_{res, tot}}{Q_{p, in}}`"
   "Total resin volume required", ":math:`V_{res, tot} = Z \pi \frac{D_{col}^2}{4} n_{op}`"
   "Volume of single column", ":math:`V_{col} = H_{col} \frac{V_{bed}}{Z}`"
   "Total column volume required", ":math:`V_{col, tot} = n_{op}V_{col}`"
   "Column height to diameter ratio", ":math:`R_{HD} = \frac{H_{col}}{D_{col}}`"
   "Column height", ":math:`H_{col} = Z + H_{distributor} + H_{underdrain} + H_{expan}`"
   "Interstitial velocity", ":math:`u_{inter} = \frac{u_{bed}}{\varepsilon}`"
   "Contact time", ":math:`t_{contact} = EBCT \varepsilon`"
   "Empty bed contact time", ":math:`EBCT = \frac{Z}{u_{bed}}`"
   "Regeneration tank volume", ":math:`V_{regen} = t_{regen} (Q_{p, in} / R)`"
   "Bed expansion fraction from backwashing (T = 20C)", ":math:`X_{expan} = H_{expan,A} + H_{expan,B}u_{bw} + H_{expan,C}u_{bw}^{2}`"
   "Bed expansion from backwashing", ":math:`H_{expan} = X_{expan}Z`"
   "Regen volumetric flow rate", ":math:`Q_{regen} = \frac{Q_{p, in}N_{regen}}{R}`"
   "Backwashing flow rate", ":math:`Q_{bw} = u_{bw} \frac{V_{bed}}{Z}n_{op}`"
   "Rinse flow rate", ":math:`Q_{rinse} = u_{bed} \frac{V_{bed}}{Z}n_{op}`"
   "Main pump power", ":math:`P_{main} = \frac{g \rho_{in} p_{drop}Q_{p, in}}{\eta} \Big( \frac{t_{break}}{t_{cycle}} \Big)`"
   "Regen pump power", ":math:`P_{regen} = \frac{g \rho_{in} p_{drop}Q_{regen}}{\eta} \Big( \frac{t_{regen}}{t_{cycle}} \Big)`"
   "Rinse pump power", ":math:`P_{rinse} = \frac{g \rho_{in} p_{drop}Q_{rinse}}{\eta} \Big( \frac{t_{rinse}}{t_{cycle}} \Big)`"
   "Backwash pump power", ":math:`P_{bw} = \frac{g \rho_{in} p_{drop}Q_{bw}}{\eta} \Big( \frac{t_{bw}}{t_{cycle}} \Big)`"
   "Pressure drop (T = 20C)", ":math:`p_{drop} = Z(p_{drop,A} + p_{drop,B}u_{bed} + p_{drop,C}u_{bed}^{2})`"
   "Rinse time", ":math:`t_{rinse} = EBCT N_{rinse}`"
   "Waste time", ":math:`t_{waste} = t_{regen} + t_{bw} + t_{rinse}`"
   "Cycle time", ":math:`t_{cycle} = t_{break} + t_{waste}`"
   "Reynolds number", ":math:`Re = \frac{u_{bed}d}{\mu}`"
   "Schmidt number", ":math:`Sc = \frac{\mu}{D}`"
   "Sherwood number", ":math:`Sh = 2.4 \varepsilon^{0.66} Re^{0.34} Sc^{0.33}`"
   "Bed Peclet number", ":math:`Pe_{bed} = Pe_{p} \frac{Z}{d}`"
   "Particle Peclet number", ":math:`Pe_{p} = 0.05 Re^{0.48}`"
   "Resin surface area per vol", ":math:`a_{s} = 6 \frac{1-\varepsilon}{d}`"

    **Langmuir**
   "Langmuir isotherm", ":math:`\frac{C_{b}}{C_{0}} (1-\frac{q_{eq}}{q_{max}}) = La (1-\frac{C_{b}}{C_{0}})\frac{q_{eq}}{q_{max}}`"
   "Constant pattern solution for Langmuir isotherm", ":math:`N(\tau - 1) = 1 + \frac{\log{(C_{b}/C_{0})} - La \log{(1 - C_{b}/C_{0})}}{1 - La}`"
   "Resin capacity mass balance", ":math:`q_{max} = q_{avail} + q_{eq}`"
   "Partition ratio", ":math:`\Lambda = \frac{q_{eq} \rho_{b}}{C_{0}}`"
   "Fluid mass transfer coeff", ":math:`k_{f} = \frac{D Sh}{d}`"
   "Number of mass-transfer units", ":math:`N = \frac{k_{f}a_{s}Z}{u_{bed}}`"
   "Dimensionless time", ":math:`\tau = (\frac{u_{inter}t_{break} \varepsilon}{Z} - \varepsilon) / \Lambda`"
   "Height of transfer unit", ":math:`HTU = \frac{u_{bed}}{\rho_{b}k}`"
   "Rate coefficient", ":math:`k = 6 \frac{(1-\varepsilon)k_{f}}{\rho_{b}d}`"
   "Mass removed", ":math:`M_{rem,j} = V_{res,tot}q_{eq} \rho_{b}`"
   "Mass transfer term", ":math:`\dot{m}_j = -M_{rem,j} / t_{break}`"

    **Freundlich**
   "Breakthrough concentration", ":math:`X = \frac{C_b}{C_0}`"
   "Bed volumes at breakthrough concentration", ":math:`BV = \frac{t_{break} u_{bed}}{Z}`"
   "Clark equation with fundamental constants", ":math:`X = \frac{1}{\bigg(1 + (2^{n - 1} - 1)\text{exp}\bigg[\frac{k_T Z (n - 1)}{BV_{50} u_{bed}} (BV_{50} - BV)\bigg]\bigg)^{\frac{1}{n-1}}}`"
   "Evenly spaced c_norm for trapezoids", ":math:`X_{trap,k} = X_{trap,min} + (k - 1) \frac{X - X_{trap,min}}{n_{trap} - 1}`"
   "Breakthru time calculation for trapezoids", ":math:`t_{trap,k} = - \log{\frac{X_{trap,k}^{n-1}-1}{A}} / k_T`"
   "Area of trapezoids", ":math:`A_{trap,k} = \frac{t_{trap,k} - t_{trap,k - 1}}{t_{trap,n_{trap}}} \frac{X_{trap,k} + X_{trap,k - 1}}{2}`"
   "Average relative effluent concentration", ":math:`X_{avg} = \sum{A_{trap,k}}`"
   "Mass transfer term", ":math:`\dot{m}_j = -(1 - X_{avg}) N_j`"


Costing Method
--------------

The following is a list of variables and/or parameters that are created when applying the ion exchange costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Default Value", "Units", "Notes"

   "Anion exchange resin cost", ":math:`c_{res}`", "``anion_exchange_resin_cost``", "205", ":math:`\text{\$/}\text{ft}^{3}`", "Assumes strong base polystyrenic gel-type Type II. From EPA-WBS cost model."
   "Cation exchange resin cost", ":math:`c_{res}`", "``cation_exchange_resin_cost``", "153", ":math:`\text{\$/}\text{ft}^{3}`", "Assumes strong acid polystyrenic gel-type. From EPA-WBS cost model."
   "Regenerant dose per volume of resin", ":math:`D_{regen}`", "``regen_dose``", "300", ":math:`\text{kg/}\text{m}^3`", "Mass of regenerant chemical per cubic meter of resin volume"
   "Ion exchange column cost equation A coeff", ":math:`C_{col,A}`", "``vessel_A_coeff``", "1596.499", ":math:`\text{\$}`", "Carbon steel w/ stainless steel internals. From EPA-WBS cost model."
   "Ion exchange column cost equation B coeff", ":math:`C_{col,b}`", "``vessel_b_coeff``", "0.459496", ":math:`\text{dimensionless}`", "Carbon steel w/ stainless steel internals. From EPA-WBS cost model."
   "Backwash/rinse tank cost equation A coeff", ":math:`C_{bw,A}`", "``backwash_tank_A_coeff``", "308.9371", ":math:`\text{\$}`", "Steel tank. From EPA-WBS cost model."
   "Backwash/rinse tank cost equation B coeff", ":math:`C_{bw,b}`", "``backwash_tank_b_coeff``", "0.501467", ":math:`\text{dimensionless}`", "Steel tank. From EPA-WBS cost model."
   "Regeneration solution tank cost equation A coeff", ":math:`C_{regen,A}`", "``regen_tank_A_coeff``", "57.02158", ":math:`\text{\$}`", "Stainless steel tank. From EPA-WBS cost model."
   "Regeneration solution tank cost equation B coeff", ":math:`C_{regen,b}`", "``regen_tank_b_coeff``", "0.729325", ":math:`\text{dimensionless}`", "Stainless steel tank. From EPA-WBS cost model."
   "Fraction of resin replaced per year", ":math:`f_{res}`", "``annual_resin_replacement_factor``", "0.05", ":math:`\text{yr}^{-1}`", "Estimated 4-5% per year. From EPA-WBS cost model."
   "Minimum hazardous waste disposal cost", ":math:`f_{haz,min}`", "``hazardous_min_cost``", "3240", ":math:`\text{\$/}\text{yr}`", "Minimum cost per hazardous waste shipment. From EPA-WBS cost model."
   "Unit cost for hazardous waste resin disposal", ":math:`f_{haz,res}`", "``hazardous_resin_disposal``", "347.10", ":math:`\text{\$/}\text{ton}`", "From EPA-WBS cost model."
   "Unit cost for hazardous waste regeneration solution disposal", ":math:`f_{haz,regen}`", "``hazardous_regen_disposal``", "3.64", ":math:`\text{\$/}\text{gal}`", "From EPA-WBS cost model."
   "Number of cycles the regenerant can be reused before disposal", ":math:`f_{recycle}`", "``regen_recycle``", "1", ":math:`\text{dimensionless}`", "Can optionally be set by the user to investigate more efficient regen regimes."
   "Costing factor to account for total installed cost installation of equipment", ":math:`f_{TIC}`", "``total_installed_cost_factor``", "1.65", ":math:`\text{dimensionless}`", "Costing factor to account for total installed cost of equipment"
   "Unit cost of NaCl", ":math:`c_{regen}`", "``costing.nacl``", "0.09", ":math:`\text{\$/}\text{kg}`", "Assumes solid NaCl. From CatCost v 1.0.4"
   "Unit cost of HCl", ":math:`c_{regen}`", "``costing.hcl``", "0.17", ":math:`\text{\$/}\text{kg}`", "Assumes 37% solution HCl. From CatCost v 1.0.4"
   "Unit cost of NaOH", ":math:`c_{regen}`", "``costing.naoh``", "0.59", ":math:`\text{\$/}\text{kg}`", "Assumes 30% solution NaOH. From iDST"
   "Unit cost of Methanol (MeOH)", ":math:`c_{regen}`", "``costing.meoh``", "3.395", ":math:`\text{\$/}\text{kg}`", "Assumes 100% pure MeOH. From ICIS"

Capital Cost Calculations
^^^^^^^^^^^^^^^^^^^^^^^^^

Capital costs for ion exchange in the ``watertap_costing_package`` are the summation of the 
total cost of the resin, columns, backwashing tank, and regeneration solution tank:

Resin is costed based on the total volume of resin required for the system, where :math:`c_{res}` is the cost per volume of resin (either cation or anion exchange resin):

.. math::
    C_{resin} = V_{res,tot} c_{res}

Vessel cost as a function of volume was fit to a power function to determine capital cost of each column:

.. math::
    C_{col} = C_{col,A} V_{col}^{C_{col,b}}
   

The backwashing tank is assumed to include backwash and rinsing volumes. The total volume of this tank is:

.. math::
    V_{bw} = Q_{bw} t_{bw} + Q_{rinse} t_{rinse}

Backwashing tank cost as a function of volume was fit to a power function to determine capital cost of the backwashing tank:

.. math::
    C_{bw} = C_{bw,A} V_{bw}^{C_{bw,b}}
   
Regeneration tank cost as a function of volume was fit to a power function to determine capital cost of the regeneration tank:

.. math::
    C_{regen} = C_{regen,A} V_{regen}^{C_{regen,b}}

And the total capital cost for the ion exchange system is the summation of these:

.. math::
    C_{tot} = ((C_{resin} + C_{col}) (n_{op} + n_{red}) + C_{bw} + C_{regen}) f_{TIC}

A total installed cost (:math:`f_{TIC}`) factor of 1.65 is applied to account for installation costs. 

.. note::
    If using ``single_use`` option for ``regenerant`` configuration keyword, the capital for the regeneration tank is zero.

Operating Cost Calculations
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The operating costs for ion exchange includes the annual resin replacement cost, regeneration solution flow, energy consumption for booster pumps, 
and any hazardous waste handling costs.

Generally, the largest operating cost is the cost of the regeneration solution. The type of regeneration solution used is set via the 
optional model configuration keyword ``regenerant``. Costing data is available for the following regenerant chemicals:

* NaCl
* HCl
* NaOH
* MeOH

If the user does not provide a value for this option, the model defaults to a NaCl regeneration solution. The dose of regenerant needed
is set by the parameter ``regen_dose`` in kg regenerant per cubic meter of resin volume. The mass flow of regenerant solution [kg/yr] is:

.. math::
    \dot{m}_{regen} = \frac{D_{regen} V_{res} (n_{op} + n_{red})}{t_{cycle} f_{recycle}}

Annual resin replacement cost is:

.. math::
    C_{op,res} = V_{res} (n_{op} + n_{red}) f_{res} c_{res}

If the spent resin and regenerant contains hazardous material, the user designates this by the model configuration keyword ``hazardous_waste``. If set to ``True``, hazardous
disposal costs are calculated as a function of the annual mass of resin replaced and regenerant consumed:

.. math::
    C_{op,haz} = f_{haz,min} + \bigg( M_{res} (n_{op} + n_{red}) f_{res} \bigg)  f_{haz,res} + \dot{v}_{regen} f_{haz,regen}

Where :math:`M_{res}` is the resin mass for a single bed and :math:`\dot{v}_{regen}` is the volumetric flow of regenerant solution. If ``hazardous_waste`` is set to ``False``,
:math:`C_{op,haz} = 0`

The total energy consumed by the unit is the summation of the power required for each of the booster pump, backwashing pump, regeneration pump, and rinsing pump. Each is scaled 
by the total time required for each step:

.. math::
    P_{tot} = \cfrac{P_{main} t_{break} + P_{bw} t_{bw} + P_{regen} t_{regen} + P_{rinse} t_{rinse}}{t_{cycle}} 

If the user chooses ``single_use`` for the ``regenerant`` configuration keyword, there is no cost for regeneration solution:

.. math::
    \dot{m}_{regen} = \dot{v}_{regen} = 0

Instead, the model assumes the entire volume of resin for the operational columns is replaced at the end of each service cycle by calculating the 
volumetric "flow" of resin:

.. math::
    \dot{v}_{resin} = \frac{V_{res, tot}}{t_{break}} 

And then operational cost of replacing the entire bed is:

.. math::
    C_{op,res} = \dot{v}_{resin} c_{res}

If ``hazardous_waste`` is set to ``True``, the hazardous waste disposal costs are: 

.. math::
    C_{op,haz} = f_{haz,min} + ( \dot{v}_{resin} \rho_{b} n_{op})  f_{haz,res}

Otherwise, :math:`C_{op,haz} = 0` as before. 

Lastly, the total energy consumed by the unit for ``single_use`` configuration includes the booster pump, backwashing pump, and rinsing pump:

.. math::
    P_{tot} = \cfrac{P_{main} t_{break} + P_{bw} t_{bw} + P_{rinse} t_{rinse}}{t_{cycle}} 

References
----------

| LeVan, M. D., Carta, G., & Yon, C. M. (2019).
| Section 16: Adsorption and Ion Exchange.
| Perry's Chemical Engineers' Handbook, 9th Edition.

| Crittenden, J. C., Trussell, R. R., Hand, D. W., Howe, K. J., & Tchobanoglous, G. (2012).
| Chapter 16: Ion Exchange.
| MWH's Water Treatment (pp. 1263-1334): John Wiley & Sons, Inc.

| DOWEX Ion Exchange Resins Water Conditioning Manual
| https://www.lenntech.com/Data-sheets/Dowex-Ion-Exchange-Resins-Water-Conditioning-Manual-L.pdf

| Inamuddin, & Luqman, M. (2012).
| Ion Exchange Technology I: Theory and Materials.

| Vassilis J. Inglezakis and Stavros G. Poulopoulos
| Adsorption, Ion Exchange and Catalysis: Design of Operations and Environmental Applications (2006).
| doi.org/10.1016/B978-0-444-52783-7.X5000-9

| Michaud, C.F. (2013)
| Hydrodynamic Design, Part 8: Flow Through Ion Exchange Beds
| Water Conditioning & Purification Magazine (WC&P)
| https://wcponline.com/2013/08/06/hydrodynamic-design-part-8-flow-ion-exchange-beds/

| Clark, R. M. (1987). 
| Evaluating the cost and performance of field-scale granular activated carbon systems. 
| Environ Sci Technol, 21(6), 573-580. 
| doi:10.1021/es00160a008

| Croll, H. C., Adelman, M. J., Chow, S. J., Schwab, K. J., Capelle, R., Oppenheimer, J., & Jacangelo, J. G. (2023). 
| Fundamental kinetic constants for breakthrough of per- and polyfluoroalkyl substances at varying empty bed contact times: 
| Theoretical analysis and pilot scale demonstration. 
| Chemical Engineering Journal, 464. 
| doi:10.1016/j.cej.2023.142587

| United States Environmental Protection Agency. (2021). Work Breakdown Structure-Based Cost Models
| https://www.epa.gov/sdwa/drinking-water-treatment-technology-unit-cost-models
