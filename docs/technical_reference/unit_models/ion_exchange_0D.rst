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
7) Constant separation-factor treatment with favorable Langmuir isotherm

Introduction
------------

Ion exchange is the reversible transfer of one or more solutes between a fluid phase and a sorbent.
This process is becoming increasingly popular in drinking water treatment applications where it is
used for water softening and demineralization. This implementation of the fixed-bed ion exchange model
accounts for process equilibrium, kinetics, and hydrodynamics to predict performance, bed and column geometry, and capital/operating costs.
The ion exchange process operates as a cycle with four steps: (1) service, i.e., treatment, (2) backwashing, (3) regeneration, (4) rinsing.
Critical to predicting performance of an ion exchange process is having an estimate for the breakthrough time,
or the duration of treatment before the solute begins exiting the column at a concentration unacceptable to the operator.
At this time, the mass transfer zone is approaching the end of the ion exchange bed, the resin is nearing exhaustion,
and the regeneration cycle can begin. Fundamental to this model is the assumption that the isotherm between the solute
and the resin is favorable, and thus the mass transfer zone is shallow.


Ports
-----

The model provides three ports (Pyomo notation in parenthesis):

* Inlet port (inlet)
* Outlet port (outlet)
* Regen port (regen)

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Phases", ":math:`p`", "['Liq']"
   "Components", ":math:`j`", "['H2O', 'Na_+', 'Ca_2+', '\Cl_-', 'Mg2+', 'SO4_2-', '\PFAS_-', 'Hardness_2+']"
   "Ions", ":math:`j`", "['Na_+', 'Ca_2+', '\Cl_-', 'Mg2+', 'SO4_2-', '\PFAS_-', 'Hardness_2+']*"

\*"Ion" is a subset of "Component" and uses the same symbol j.
**NOTE:The "Components" and "Ions" lists can include any ion as long as the ion is configured into the property package.**


Degrees of Freedom
------------------
Aside from the inlet feed state variables (i.e., temperature, pressure, component molar flowrate),

the Ion Exchange (0D) model has at least an additional 27 degrees of freedom that
the user must specify. The table below gives an outline of these.

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Inlet temperature", ":math:`T`", "temperature", "[t]", ":math:`\text{K}`"
   "Inlet pressure", ":math:`P`", "pressure", "[t]", ":math:`\text{Pa}`"
   "Component molar flow rate", ":math:`N_j`", "flow_mol_phase_comp", "[t, 'Liq', 'H2O']", ":math:`\text{mol/s}`"
   "Langmuir equilibrium parameter for resin/ion system", ":math:`La`", "langmuir", "``target_ion_set``", ":math:`\text{dimensionless}`"
   "Maximum resin capacity", ":math:`q_{max}`", "resin_max_capacity", "None", ":math:`\text{mol/kg}`"
   "Service flow rate through resin bed in bed volumes per hour", ":math:`SFR`", "service_flow_rate", "None", ":math:`\text{hr}^{-1}`"
   "Number of operational columns", ":math:`n_{op}`", "number_columns", "None", ":math:`\text{dimensionless}`"
   "Number of redundant columns", ":math:`n_{red}`", "number_columns_redund", "None", ":math:`\text{dimensionless}`"
   "Bed depth", ":math:`Z`", "bed_depth", "None", ":math:`\text{m}`"
   "Resin bead diameter", ":math:`d`", "resin_diam", "None", ":math:`\text{m}`"
   "Resin bulk density", ":math:`\rho_{b}`", "resin_bulk_dens", "None", ":math:`\text{kg/L}`"
   "Bed porosity", ":math:`\epsilon`", "bed_porosity", "None", ":math:`\text{dimensionless}`"
   "Dimensionless time", ":math:`\tau`", "dimensionless_time", None, ":math:`\text{dimensionless}`"
   "Regenerant dose per volume of resin", ":math:`C_{regen}`", "regen_dose", "None", ":math:`\text{kg/}\text{m}^3`"
   "Number of cycles before regenerant disposal", ":math:`N_{regen}`", "regen_recycle", "None", ":math:`\text{dimensionless}`"
   "Regeneration time", ":math:`t_{regen}`", "t_regen", "None", ":math:`\text{s}`"
   "Backwash time", ":math:`t_{bw}`", "t_bw", "None", ":math:`\text{s}`"
   "Backwash loading rate", ":math:`u_{bw}`", "bw_rate", "None", ":math:`\text{m/hr}`"
   "Number of bed volumes for rinse step", ":math:`N_{rinse}`", "rinse_bv", "None", ":math:`\text{dimensionless}`"
   "Pump efficiency", ":math:`\eta`", "pump_efficiency", "None", ":math:`\text{dimensionless}`"
   "Pressure drop equation intercept", ":math:`P_{drop,A}`", "p_drop_A", "None", ":math:`\text{dimensionless}`"
   "Pressure drop equation B", ":math:`P_{drop,B}`", "p_drop_B", "None", ":math:`\text{dimensionless}`"
   "Pressure drop equation C", ":math:`P_{drop,C}`", "p_drop_C", "None", ":math:`\text{dimensionless}`"
   "Bed expansion fraction eq intercept", ":math:`H_{expan,A}`", "bed_expansion_frac_A", "None", ":math:`\text{dimensionless}`"
   "Bed expansion fraction equation B parameter", ":math:`H_{expan,B}`", "bed_expansion_frac_B", "None", ":math:`\text{dimensionless}`"
   "Bed expansion fraction equation C parameter", ":math:`H_{expan,C}`", "bed_expansion_frac_C", "None", ":math:`\text{dimensionless}`"
   "Service-to-regeneration flow ratio", ":math:`R`", "service_to_regen_flow_ratio", "None", ":math:`\text{dimensionless}`"


**Users must provide values for and 'fix' the following variables to solve the model with DOF=0: 'pressure', 'temperature', 'flow_mol_phase_comp', 'langmuir', 'resin_max_capacity', 'service_flow_rate', 'number_columns', and 'bed_depth'. The other variables can simply be fixed to their default values ('.fix()').**

**NOTE: Variables for 'temperature', 'pressure', and 'flow_mol_phase_comp' come from the associated property package as state variables and are accessed via {port_name}.{state_var_name}**

.. _IX_variables:

Variables
---------

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "**Resin Variables**"
   "Maximum resin capacity", ":math:`q_{max}`", "resin_max_capacity", "None", ":math:`\text{mol/kg}`"
   "Usable resin capacity at equilibrium", ":math:`q_{eq}`", "resin_eq_capacity", "None", ":math:`\text{mol/kg}`"
   "Available resin capacity at equilibrium", ":math:`q_{avail}`", "resin_unused_capacity", "None", ":math:`\text{dimensionless}`"
   "Resin bead diameter", ":math:`d`", "resin_diam", "None", ":math:`\text{m}`"
   "Resin bulk density", ":math:`\rho_{b}`", "resin_bulk_dens", "None", ":math:`\text{kg/L}`"
   "Resin particle density", ":math:`\rho_{p}`", "resin_particle_dens", "None", ":math:`\text{dimensionless}`"
   "Separation factor", ":math:`\alpha`", "separation_factor", "``target_ion_set``", ":math:`\text{dimensionless}`"
   "Resin surface area per volume", ":math:`a_{s}`", "resin_surf_per_vol", "None", ":math:`\text{m}^{-1}`"
   "Langmuir equilibrium parameter for resin/ion system", ":math:`La`", "langmuir", "``target_ion_set``", ":math:`\text{dimensionless}`"

   "**Bed/Column Variables**"
   "Ratio of bed depth to column diameter", ":math:`X`", "bed_depth_to_diam_ratio", "None", ":math:`\text{dimensionless}`"
   "Bed volume of one unit", ":math:`V_{bed}`", "bed_vol", "None", ":math:`\text{m}^{3}`"
   "Total bed volume", ":math:`V_{tot}`", "bed_vol_tot", "None", ":math:`\text{m}^{3}`"
   "Bed depth", ":math:`Z`", "bed_depth", "None", ":math:`\text{m}`"
   "Bed porosity", ":math:`\epsilon`", "bed_porosity", "None", ":math:`\text{dimensionless}`"
   "Column height", ":math:`H`", "col_height", "None", ":math:`\text{m}`"
   "Column diameter", ":math:`D_{col}`", "col_diam", "None", ":math:`\text{m}`"
   "Column volume of one unit", ":math:`V_{col}`", "col_vol_per", "None", ":math:`\text{m}^{3}`"
   "Total column volume", ":math:`V_{col, tot}`", "col_vol_tot", "None", ":math:`\text{m}^{3}`"
   "Number of operational columns", ":math:`n_{op}`", "number_columns", "None", ":math:`\text{dimensionless}`"
   "Number of redundant columns", ":math:`n_{red}`", "number_columns_redund", "None", ":math:`\text{dimensionless}`"
   "Underdrain height", ":math:`H_{underdrain}`", "underdrain_h", "None", ":math:`\text{m}`"
   "Distributor height", ":math:`H_{distributor}`", "distributor_h", "None", ":math:`\text{m}`"

   "**Kinetic Variables**"
   "Partition ratio", ":math:`\Lambda`", "partition_ratio", "None", ":math:`\text{dimensionless}`"
   "Fluid mass transfer coefficient", ":math:`k_{f}`", "fluid_mass_transfer_coeff", "``target_ion_set``", ":math:`\text{m/s}`"
   "Rate coefficient based on fluid-phase concentration driving force", ":math:`k`", "rate_coeff", "``target_ion_set``", ":math:`\text{m}^{3}\text{kg*s}`"
   "Number of transfer units", ":math:`N`", "num_transfer_units", "None", ":math:`\text{dimensionless}`"
   "Height of a transfer unit", ":math:`HTU`", "HTU", "``target_ion_set``", ":math:`\text{m}`"
   "Position of breakthrough on constant-pattern wave", ":math:`lh`", "lh", "None", ":math:`\text{dimensionless}`"
   "Influent mass of ion", ":math:`M_{in}`", "mass_in", "``target_ion_set``", ":math:`\text{mol}`"
   "Sorbed mass of ion", ":math:`M_{out}`", "mass_removed", "``target_ion_set``", ":math:`\text{mol}`"
   "Effluent mass of ion", ":math:`M_{rem}`", "mass_out", "``target_ion_set``", ":math:`\text{mol}`"

   "**Hydrodynamic Variables**"
   "Service flow rate through resin bed in bed volumes per hour", ":math:`SFR`", "service_flow_rate", "None", ":math:`\text{hr}^{-1}`"
   "Velocity through resin bed", ":math:`u_{bed}`", "vel_bed", "None", ":math:`\text{m/s}`"
   "Interstitial velocity", ":math:`u_{inter}`", "vel_inter", "None", ":math:`\text{m/s}`"
   "Holdup percent", ":math:`holdup`", "holdup", "None", ":math:`\text{dimensionless}`"
   "Pressure drop through resin bed", ":math:`P_{drop}`", "pressure_drop", "None", ":math:`\text{psi}`"
   "Pressure drop equation intercept", ":math:`P_{drop,A}`", "p_drop_A", "None", ":math:`\text{dimensionless}`"
   "Pressure drop equation B", ":math:`P_{drop,B}`", "p_drop_B", "None", ":math:`\text{dimensionless}`"
   "Pressure drop equation C", ":math:`P_{drop,C}`", "p_drop_C", "None", ":math:`\text{dimensionless}`"

   "**Time Variables**"
   "Rinse time", ":math:`t_{rinse}`", "t_rinse", "None", ":math:`\text{s}`"
   "Dimensionless time", ":math:`\tau`", "dimensionless_time", None, ":math:`\text{dimensionless}`"
   "Breakthrough time", ":math:`t_{breakthru}`", "t_breakthru", "None", ":math:`\text{s}`"
   "Cycle time", ":math:`t_{cycle}`", "t_cycle", "None", ":math:`\text{s}`"
   "Contact time", ":math:`t_{contact}`", "t_contact", "None", ":math:`\text{s}`"
   "Regen + Rinse + Backwash time", ":math:`t_{waste}`", "t_waste", "None", ":math:`\text{s}`"
   "Regeneration time", ":math:`t_{regen}`", "t_regen", "None", ":math:`\text{s}`"
   "Backwash time", ":math:`t_{bw}`", "t_bw", "None", ":math:`\text{s}`"

   "**Dimensionless Variables**"
   "Reynolds number", ":math:`Re`", "Re", "None", ":math:`\text{dimensionless}`"
   "Schmidt number", ":math:`Sc`", "Sc", "``target_ion_set``", ":math:`\text{dimensionless}`"
   "Sherwood number", ":math:`Sh`", "Sh", "``target_ion_set``", ":math:`\text{dimensionless}`"
   "Peclet particle number", ":math:`Pe_{p}`", "Pe_p", "None", ":math:`\text{dimensionless}`"
   "Peclet bed number", ":math:`Pe_{bed}`", "Pe_bed", "None", ":math:`\text{dimensionless}`"
   "Ratio of breakthrough concentration to influent concentration", ":math:`C_{b}/C_{0}`", "c_norm", "``target_ion_set``", ":math:`\text{dimensionless}`"

   "**Regeneration Variables**"
   "Service-to-regeneration flow ratio", ":math:`R`", "service_to_regen_flow_ratio", "None", ":math:`\text{dimensionless}`"
   "Number of cycles before regenerant disposal", ":math:`N_{regen}`", "regen_recycle", "None", ":math:`\text{dimensionless}`"
   "Regenerant dose per volume of resin", ":math:`C_{regen}`", "regen_dose", "None", ":math:`\text{kg/}\text{m}^3`"

   "**Backwashing Variables**"
   "Backwashing volumetric flow rate", ":math:`Q_{bw}`", "bw_flow", "None", ":math:`\text{m}^{3}\text{/s}`"
   "Backwash loading rate", ":math:`u_{bw}`", "bw_rate", "None", ":math:`\text{m/hr}`"
   "Fraction of bed depth increase during backwashing", ":math:`X_{expan}`", "bed_expansion_frac", "None", ":math:`\text{dimensionless}`"
   "Additional column sidewall height required for bed expansion", ":math:`H_{expan}`", "bed_expansion_h", "None", ":math:`\text{dimensionless}`"
   "Bed expansion fraction eq intercept", ":math:`H_{expan,A}`", "bed_expansion_frac_A", "None", ":math:`\text{dimensionless}`"
   "Bed expansion fraction equation B parameter", ":math:`H_{expan,B}`", "bed_expansion_frac_B", "None", ":math:`\text{dimensionless}`"
   "Bed expansion fraction equation C parameter", ":math:`H_{expan,C}`", "bed_expansion_frac_C", "None", ":math:`\text{dimensionless}`"

   "**Rinsing Variables**"
   "Rinse volumetric flow rate", ":math:`Q_{rinse}`", "rinse_flow", "None", ":math:`\text{m}^{3}\text{/s}`"
   "Number of bed volumes for rinse step", ":math:`N_{rinse}`", "rinse_bv", "None", ":math:`\text{dimensionless}`"
   "Power of main booster pump", ":math:`P_{main}`", "main_pump_power", "None", ":math:`\text{kW}`"
   "Regen pump power", ":math:`P_{regen}`", "regen_pump_power", "None", ":math:`\text{kW}`"
   "Backwash pump power", ":math:`P_{bw}`", "bw_pump_power", "None", ":math:`\text{kW}`"
   "Rinse pump power", ":math:`P_{rinse}`", "rinse_pump_power", "None", ":math:`\text{kW}`"
   "Assumed efficiency for all pumps", ":math:`\eta`", "pump_efficiency", "None", ":math:`\text{dimensionless}`"


Solution Component Information
------------------------------
In addition to providing a list of solute ions, the users will
need to provide parameter information for each ion including molecular weight,
diffusivity data, and charge data.

To provide this information to the unit model, users must add
dictionaries to the initialization of the unit model. These dictionaries must have the
following format.

.. code-block::

   def get_ix_in(ions):
    diff_data = {
        "Na_+": 1.33e-9,
        "Ca_2+": 9.2e-10,
        "Cl_-": 2.03e-9,
        "Mg_2+": 0.706e-9,
        "SO4_2-": 1.06e-9,
        "PFAS_-": 0.49e-9,
        "Hardness_2+": 0.706e-9,
    }
    mw_data = {
        "Na_+": 23e-3,
        "Ca_2+": 40e-3,
        "Cl_-": 35e-3,
        "Mg_2+": 24e-3,
        "SO4_2-": 96e-3,
        "PFAS_-": 414.1e-3,
        "Hardness_2+": 100.0869e-3,
    }
    charge_data = {
        "Na_+": 1,
        "Ca_2+": 2,
        "Cl_-": -1,
        "Mg_2+": 2,
        "SO4_2-": -2,
        "PFAS_-": -1,
        "Hardness_2+": 2,
    }
    ix_in = {
        "solute_list": [],
        "diffusivity_data": {},
        "mw_data": {"H2O": 18e-3},
        "charge": {},
    }
    for ion in ions:
        ix_in["solute_list"].append(ion)
        ix_in["diffusivity_data"][("Liq", ion)] = diff_data[ion]
        ix_in["mw_data"][ion] = mw_data[ion]
        ix_in["charge"][ion] = charge_data[ion]
    return ix_in

**NOTE: 'ions' is an ion_set, which is a configuration argument of the property package as shown below**


.. code-block::

        ions = m.fs.unit.config.property_package.ion_set

**NOTE: The above example assumes you have already constructed a pyomo model named 'm' and attached an IDAES flowsheet named 'fs' to it.**

Equations and Relationships
---------------------------

.. csv-table::
   :header: "Description", "Equation"

   "Separation factor", ":math:`\alpha = \frac{1}{La}`"
   "Langmuir isotherm", ":math:`\alpha \frac{C_{b}}{C_{0}} (1-\frac{q_{eq}}{q_{max}}) = (1-\frac{C_{b}}{C_{0}})\frac{q_{eq}}{q_{max}}`"
   "Reynolds number", ":math:`Re = \frac{u_{bed}d}{\mu}`"
   "Schmidt number", ":math:`Sc = \frac{\mu}{D}`"
   "Sherwood number", ":math:`Sh = \frac{1.09}{\epsilon}Re^{0.33}Sc^{0.33}`"
   "Bed Peclet number", ":math:`Pe_{bed} = Pe_{p} \frac{Z}{d}`"
   "Particle Peclet number", ":math:`Pe_{p} = 0.05 Re^{0.48}`"
   "Resin capacity mass balance", ":math:`q_{max} = q_{avail} + q_{eq}`"
   "Interstitial velocity", ":math:`u_{inter} = \frac{u_{bed}}{\epsilon}`"
   "Resin surface area per vol", ":math:`a_{s} = 6 \frac{1-\epsilon}{d}`"
   "Contact time", ":math:`t_{contact} = \frac{Z}{u_{inter}}`"
   "Service flow rate", ":math:`SFR = \frac{Q_{p, in}}{V_{tot}}`"
   "Flow through bed constraint", ":math:`\frac{Z \epsilon}{u_{bed}} = \frac{V_{bed} \epsilon}{Q_{p, in} / n_{op}}`"
   "Total bed volume", ":math:`V_{tot} = V_{bed}n_{op}`"
   "Column height", ":math:`H = Z + H_{distributor} + H_{underdrain} + H_{expan}`"
   "Column volume calculated from bed volume", ":math:`V_{col} = H \frac{V_{bed}}{Z}`"
   "Column volume calculated from column diameter", ":math:`V_{col} = \pi (\frac{D_{col}}{2})^{2} H`"
   "Column diameter calculation", ":math:`(\frac{D_{col}}{2})^{2} = (\frac{H}{2X})^{2}`"
   "Fluid mass transfer coeff", ":math:`k_{f} = \frac{D Sh}{d}`"
   "Rate coefficient", ":math:`k = 6 \frac{(1-\epsilon)k_{f}}{\rho_{b}d}`"
   "Height of transfer unit", ":math:`HTU = \frac{u_{bed}}{\rho_{b}k}`"
   "Partition ratio", ":math:`\Lambda = \frac{q_{eq} \rho_{b}}{ñ_{in}}`"
   "Left hand side of constant pattern solution", ":math:`lh = N(\tau - 1)`"
   "Right hand side of constant pattern solution", ":math:`lh = 1 + \frac{\log{(C_{b}/C_{0})} - La \log{(1 - C_{b}/C_{0})}}{1 - La}`"
   "Dimensionless time", ":math:`\tau = (\frac{u_{inter}t_{breakthru} \epsilon}{Z} - \epsilon) / \Lambda`"
   "Number of mass-transfer units", ":math:`N = \frac{k_{f}a_{s}Z}{u_{bed}}`"
   "Flow conservation", ":math:`Q_{p, in} - \frac{Q_{bw}t_{bw} + Q_{rinse}t_{rinse}}{t_{cycle}} = Q_{p, out} - \frac{Q_{regen}t_{regen}}{t_{cycle}}`"
   "Influent total mass of ion", ":math:`M_{in} = Q_{p, in}t_{breakthru}ñ_{in}`"
   "Removed total mass of ion", ":math:`M_{rem} = V_{bed}q_{eq}n_{op} \rho_{b}`"
   "Mass of ion in effluent", ":math:`M_{out} = M_{in} - M_{rem}`"
   "Steady-state effluent concentration (for target ion)", ":math:`ñ_{out} = \frac{M_{out}}{Q_{p, in}t_{breakthru}}`"
   "Steady-state effluent concentration", ":math:`ñ_{out} = ñ_{in}`"
   "Steady-state regen concentration (for target ion)", ":math:`ñ_{regen} = \frac{M_{rem}N_{regen}}{Q_{p, regen}t_{regen}}`"
   "Steady-state regen concentration", ":math:`ñ_{regen} = 0`"
   "Cycle time", ":math:`t_{cycle} = t_{breakthru} + t_{waste}`"
   "Waste time", ":math:`t_{waste} = t_{regen} + t_{bw} + t_{rinse}`"
   "Regen volumetric flow rate", ":math:`Q_{p, regen} = \frac{Q_{p, in}N_{regen}}{R}`"
   "Regen pump power", ":math:`P_{regen} = \frac{9.81 \rho_{in} 0.70325P_{drop}Q_{p, regen}}{\eta}`"
   "Bed expansion fraction from backwashing (T = 20C)", ":math:`X_{expan} = H_{expan,A} + H_{expan,B}u_{bw} + H_{expan,C}u_{bw}^{2}`"
   "Bed expansion from backwashing", ":math:`H_{expan} = X_{expan}Z`"
   "Backwashing flow rate", ":math:`Q_{bw} = u_{bw} \frac{V_{bed}}{Z}n_{op}`"
   "Backwash pump power", ":math:`P_{bw} = \frac{9.81 \rho_{in} 0.70325P_{drop}Q_{bw}}{\eta}`"
   "Rinse time", ":math:`t_{rinse} t_{contact} + N_{rinse}`"
   "Rinse flow rate", ":math:`Q_{rinse} = u_{bed} \frac{V_{bed}}{Z}n_{op}`"
   "Rinse pump power", ":math:`P_{rinse} = \frac{9.81 \rho_{in} 0.70325P_{drop}Q_{rinse}}{\eta}`"
   "Main pump power", ":math:`P_{main} = \frac{9.81 \rho_{in} 0.70325P_{drop}Q_{p, in}}{\eta}`"
   "Pressure drop (T = 20C)", ":math:`P_{drop} = Z(P_{drop,A} + P_{drop,B}u_{bed} + P_{drop,C}u_{bed}^{2})`"
   "Total column volume required", ":math:`V_{col, tot} = n_{op}V_{col}`"


References
----------
Hand, D. W., Crittenden, J. C., & Thacker, W. E. (1984). Simplified models for design of fixed-bed adsorption systems.
Journal of Environmental Engineering, 110(2), 440-456.

Crittenden, J., Rhodes, R., Hand, D., Howe, K., & Tchobanoglous, G. (2012). MWHs Water Treatment. Principles and Design.
EditorialJohn Wiley & Sons.

LeVan, M. D., Carta, G., & Yon, C. M. (2019). Section 16: Adsorption and Ion Exchange. Perry's Chemical Engineers' Handbook, 9th Edition.

Inamuddin, & Luqman, M. (2012). Ion Exchange Technology I: Theory and Materials.
