Ion Exchange 0D
===============

.. note::

    Documentation for the Ion Exchange (0D) model is undergoing refinement.

Introduction
------------
[Add intro and image]

The main assumptions of the implemented model are as follows:

1) Model dimensionality is limited to a 0D control volume
2) Single liquid phase only
3) Steady state only


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
   "Ion", ":math:`j`", "['Na_+', 'Ca_2+', '\Cl_-', 'Mg2+', 'SO4_2-', '\PFAS_-', 'Hardness_2+'] \  :sup:`1`"

 :sup:`1` "Ion" is a subset of "Component" and uses the same symbol j.

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
   "Langmuir isotherm coefficient", ":math:`?`", "langmuir", "[j]", ":math:`\text{dimensionless}`"
   "Resin max capacity", ":math:`n^{ref}_{i}`", "resin_max_capacity", "None", ":math:`\text{mol/kg}`"
   "Service flow rate (in bed volumes)", ":math:`?`", "service_flow_rate", "None", ":math:`\text{hr}^{-1}`"
   "Number of operational columns", ":math:`?`", "number_columns", "None", ":math:`\text{dimensionless}`"
   "Number of redundant columns", ":math:`?`", "number_columns_redund", "None", ":math:`\text{dimensionless}`"
   "Bed depth", ":math:`?`", "bed_depth", "None", ":math:`\text{m}`"
   "Resin bead diameter", ":math:`?`", "resin_diam", "None", ":math:`\text{m}`"
   "Resin bulk density", ":math:`?`", "resin_bulk_dens", "None", ":math:`\text{kg/L}`"
   "Bed porosity", ":math:`?`", "bed_porosity", "None", ":math:`\text{dimensionless}`"
   "Dimensionless time", ":math:`?`", "dimensionless_time", None, ":math:`\text{dimensionless}`"
   "Regenerant dose per volume of resin", ":math:`?`", "regen_dose", "None", ":math:`\text{kg/}\text{m}^3`"
   "Number of cycles before regenerant disposal", ":math:`?`", "regen_recycle", "None", ":math:`\text{dimensionless}`"
   "Regeneration time", ":math:`?`", "t_regen", "None", ":math:`\text{s}`"
   "Backwash time", ":math:`?`", "t_bw", "None", ":math:`\text{s}`"
   "Backwash loading rate", ":math:`?`", "bw_rate", "None", ":math:`\text{m/hr}`"
   "Number of bed volumes for rinse step", ":math:`?`", "rinse_bv", "None", ":math:`\text{dimensionless}`"
   "Pump efficiency", ":math:`?`", "pump_efficiency", "None", ":math:`\text{dimensionless}`"
   "Pressure drop equation intercept", ":math:`?`", "p_drop_A", "None", ":math:`\text{dimensionless}`"
   "Pressure drop equation B", ":math:`?`", "p_drop_B", "None", ":math:`\text{dimensionless}`"
   "Pressure drop equation C", ":math:`?`", "p_drop_C", "None", ":math:`\text{dimensionless}`"
   "Bed expansion fraction eq intercept", ":math:`?`", "bed_expansion_frac_A", "None", ":math:`\text{dimensionless}`"
   "Bed expansion fraction equation B parameter", ":math:`?`", "bed_expansion_frac_B", "None", ":math:`\text{dimensionless}`"
   "Bed expansion fraction equation C parameter", ":math:`?`", "bed_expansion_frac_C", "None", ":math:`\text{dimensionless}`"
   "Service-to-regeneration flow ratio", ":math:`?`", "service_to_regen_flow_ratio", "None", ":math:`\text{dimensionless}`"

**Users must provide values for and 'fix' these variables to solve the model with DOF=0. However, users may also leave variables unfixed for optimization purposes.**

**NOTE: Variables for 'temperature', 'pressure', and 'flow_mol_phase_comp' come from the associated property package as state variables and are accessed via {port_name}.{state_var_name}**


Additional Variables and Parameters
-----------------------------------

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "**Parameters**"
   "Diffusivity of ion through resin bead", ":math:`?`", "diff_ion_resin", "None", ":math:`\text{m}^{2}\text{/s}`"
   "Underdrain height", ":math:`?`", "underdrain_h", "None", ":math:`\text{m}`"
   "Distributor height", ":math:`?`", "distributor_h", "None", ":math:`\text{m}`"
   "Conversion for pressure drop in psi to m", ":math:`?`", "p_drop_psi_to_m", "None", ":math:`\text{m/psi}`"
   "Holdup equation A parameter", ":math:`?`", "holdup_A", "None", ":math:`\text{dimensionless}`"
   "Holdup equation B parameter", ":math:`?`", "holdup_B", "None", ":math:`\text{dimensionless}`"
   "Holdup equation exponent", ":math:`?`", "holdup_exp", "None", ":math:`\text{dimensionless}`"
   "Peclet particle equation A parameter", ":math:`?`", "Pe_p_A", "None", ":math:`\text{dimensionless}`"
   "Peclet particle equation exponent", ":math:`?`", "Pe_p_exp", "None", ":math:`\text{dimensionless}`"
   "Sherwood equation A parameter", ":math:`?`", "Sh_A", "None", ":math:`\text{dimensionless}`"
   "Sherwood equation exp", ":math:`?`", "Sh_exp", "None", ":math:`\text{dimensionless}`"


   "**Bed/Column Variables**"
   "Minimum ratio of bed depth to diameter", ":math:`?`", "bed_depth_to_diam_ratio", "None", ":math:`\text{dimensionless}`"
   "Bed volume of one unit", ":math:`?`", "bed_vol", "None", ":math:`\text{m}^{3}`"
   "Total bed volume", ":math:`?`", "bed_vol_tot", "None", ":math:`\text{m}^{3}`"
   "Column height", ":math:`?`", "col_height", "None", ":math:`\text{m}`"
   "Column diameter", ":math:`?`", "col_diam", "None", ":math:`\text{m}`"
   "Column volume", ":math:`?`", "col_vol_per", "None", ":math:`\text{m}^{3}`"
   "**Resin Variables**"
   "Resin equilibrium capacity", ":math:`n_{i}`", "resin_eq_capacity", "None", ":math:`\text{mol/kg}`"
   "Resin available capacity", ":math:`?`", "resin_unused_capacity", "None", ":math:`\text{dimensionless}`"
   "Resin particle density", ":math:`?`", "resin_particle_dens", "None", ":math:`\text{dimensionless}`"
   "Separation factor", ":math:`r`", "separation_factor", "[j]", ":math:`\text{dimensionless}`"
   "Resin surface area per volume", ":math:`?`", "resin_surf_per_vol", "None", ":math:`\text{m}^{-1}`"
   "**Kinetic Variables**"
   "Partition ratio", ":math:`?`", "partition_ratio", "None", ":math:`\text{dimensionless}`"
   "Fluid mass transfer coefficient", ":math:`?`", "fluid_mass_transfer_coeff", "[j]", ":math:`\text{m/s}`"
   "Rate coefficient", ":math:`?`", "rate_coeff", "[j]", ":math:`\text{m}^{3}\text{kg*s}`"
   "Breakthrough time", ":math:`?`", "t_breakthru", "None", ":math:`\text{s}`"
   "Cycle time", ":math:`?`", "t_cycle", "None", ":math:`\text{s}`"
   "Contact time", ":math:`?`", "t_contact", "None", ":math:`\text{s}`"
   "Regen + Rinse + Backwash time", ":math:`?`", "t_waste", "None", ":math:`\text{s}`"
   "Number of transfer units", ":math:`?`", "num_transfer_units", "None", ":math:`\text{dimensionless}`"
   "Height of a transfer unit", ":math:`?`", "HTU", "[j]", ":math:`\text{m}`"
   "Position of breakthrough on constant-pattern wave", ":math:`?`", "lh", "None", ":math:`\text{dimensionless}`"
   "Influent mass of ion", ":math:`?`", "mass_in", "[j]", ":math:`\text{mol}`"
   "Sorbed mass of ion", ":math:`?`", "mass_removed", "[j]", ":math:`\text{mol}`"
   "Effluent mass of ion", ":math:`?`", "mass_out", "[j]", ":math:`\text{mol}`"
   "**Hydrodynamic Variables**"
   "Velocity through resin bed", ":math:`?`", "vel_bed", "None", ":math:`\text{m/s}`"
   "Interstitial velocity", ":math:`?`", "vel_inter", "None", ":math:`\text{m/s}`"
   "Holdup percent", ":math:`?`", "holdup", "None", ":math:`\text{dimensionless}`"
   "Pressure drop across column", ":math:`?`", "pressure_drop", "None", ":math:`\text{psi}`"
   "**Dimensionless Variables**"
   "Reynolds number", ":math:`?`", "Re", "None", ":math:`\text{dimensionless}`"
   "Schmidt number", ":math:`?`", "Sc", "[j]", ":math:`\text{dimensionless}`"
   "Sherwood number", ":math:`?`", "Sh", "[j]", ":math:`\text{dimensionless}`"
   "Peclet particle number", ":math:`?`", "Pe_p", "None", ":math:`\text{dimensionless}`"
   "Peclet bed number", ":math:`?`", "Pe_bed", "None", ":math:`\text{dimensionless}`"
   "Dimensionless concentration", ":math:`c^{*}_{i}`", "c_norm", "[j]", ":math:`\text{dimensionless}`"
   "**Backwashing**"
   "Backwashing volumetric flow rate", ":math:`?`", "bw_flow", "None", ":math:`\text{m}^{3}\text{/s}`"
   "Fraction of bed depth increase during backwashing", ":math:`?`", "bed_expansion_frac", "None", ":math:`\text{dimensionless}`"
   "Additional column sidewall height required for bed expansion", ":math:`?`", "bed_expansion_h", "None", ":math:`\text{dimensionless}`"
   "**Rinse**"
   "Rinse volumetric flow rate", ":math:`?`", "rinse_flow", "None", ":math:`\text{m}^{3}\text{/s}`"
   "Rinse time", ":math:`?`", "t_rinse", "None", ":math:`\text{s}`"
   "Main pump power", ":math:`?`", "main_pump_power", "None", ":math:`\text{kW}`"
   "Regen pump power", ":math:`?`", "regen_pump_power", "None", ":math:`\text{kW}`"
   "Backwash pump power", ":math:`?`", "bw_pump_power", "None", ":math:`\text{kW}`"
   "Rinse pump power", ":math:`?`", "rinse_pump_power", "None", ":math:`\text{kW}`"


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

   "Langmuir isotherm", ":math:`r = \frac{c^{*}_{i}(1-n^{*}_i})}{n^{*}_{i}(1-c^{*}_i}`"
   "Reynolds number", ":math:``"
   "Schmidt number", ":math:``"
   "Sherwood number", ":math:``"
   "Bed Peclet number", ":math:``"
   "Particle Peclet number", ":math:``"
   "Resin capacity mass balance", ":math:``"
   "Interstitial velocity", ":math:``"
   "Column holdup", ":math:``"
   "Resin surface area per vol", ":math:``"
   "Contact time", ":math:``"
   "Service flow rate", ":math:``"
   "Flow through bed constraint", ":math:``"
   "Total bed volume", ":math:``"
   "Column height", ":math:``"
   "Column volume calculated from bed volume", ":math:``"
   "Column volume calculated from column diameter", ":math:``"
   "Column diameter calculation", ":math:``"
   "Fluid mass transfer coeff", ":math:``"
   "Rate coefficient", ":math:``"
   "Height of transfer unit - HTU", ":math:``"
   "Partition ratio", ":math:``"
   "Left hand side of constant pattern solution", ":math:``"
   "Right hand side of constant pattern solution", ":math:``"
   "Dimensionless time", ":math:``"
   "Number of mass-transfer units", ":math:``"
   "Flow conservation", ":math:``"
   "Influent total mass of ion", ":math:``"
   "Removed total mass of ion", ":math:``"
   "Mass of ion in effluent", ":math:``"
   "Steady-state effluent concentration", ":math:``"
   "Steady-state regen concentration", ":math:``"
   "Cycle time", ":math:``"
   "Waste time", ":math:``"
   "Regen volumetric flow rate", ":math:``"
   "Regen pump power", ":math:``"
   "Bed expansion fraction from backwashing", ":math:``"
   "Bed expansion from backwashing", ":math:``"
   "Backwashing flow rate", ":math:``"
   "Backwash pump power", ":math:``"
   "Rinse time", ":math:``"
   "Rinse flow rate", ":math:``"
   "Rinse pump power", ":math:``"
   "Main pump power", ":math:``"
   "Pressure drop", ":math:``"
   "Total column volume required", ":math:``"


References
----------
LeVan, M. D., Carta, G., & Yon, C. M. (2019).
Section 16: Adsorption and Ion Exchange.
Perry's Chemical Engineers' Handbook, 9th Edition.

Crittenden, J. C., Trussell, R. R., Hand, D. W., Howe, K. J., & Tchobanoglous, G. (2012).
Chapter 16: Ion Exchange.
MWH's Water Treatment (pp. 1263-1334): John Wiley & Sons, Inc.

DOWEX Ion Exchange Resins Water Conditioning Manual
https://www.lenntech.com/Data-sheets/Dowex-Ion-Exchange-Resins-Water-Conditioning-Manual-L.pdf

Inamuddin, & Luqman, M. (2012).
Ion Exchange Technology I: Theory and Materials.

Michaud, C.F. (2013)
Hydrodynamic Design, Part 8: Flow Through Ion Exchange Beds
Water Conditioning & Purification Magazine (WC&P)
https://wcponline.com/2013/08/06/hydrodynamic-design-part-8-flow-ion-exchange-beds/
