Reverse Osmosis Unit (0D)
=========================
This reverse osmosis (RO) unit model:
   * is a 0-dimensional model
   * supports a single liquid phase only
   * supports steady-state only

.. index::
   pair: proteuslib.unit_models.reverse_osmosis_0D;ReverseOsmosis0D

.. currentmodule:: proteuslib.unit_models.reverse_osmosis_0D

Example
-------
The example below shows how to setup a simple RO unit model.

.. code-block:: python

    # Import NaCl property model
    import proteuslib.property_models.NaCl_prop_pack as props
    # Import utility tool for calculating scaling factors
    import idaes.core.util.scaling as calculate_scaling_factors

    # Create a concrete model, flowsheet, and NaCl property parameter block.
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = props.NaClParameterBlock()

    # Add an RO unit to the flowsheet.
    m.fs.unit = ReverseOsmosis0D(default={"property_package": m.fs.properties})

    # Specify system variables.
    m.fs.unit.inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl'].fix(0.035) # mass flow rate of NaCl
    m.fs.unit.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O'].fix(0.965)  # mass flow rate of water
    m.fs.unit.inlet.pressure[0].fix(50e5)  # feed pressure
    m.fs.unit.inlet.temperature[0].fix(298.15)  # feed temperature
    m.fs.unit.deltaP.fix(-3e5)  # membrane pressure drop
    m.fs.unit.area.fix(50)  # membrane area
    m.fs.unit.A_comp.fix(4.2e-12)  # membrane water permeability
    m.fs.unit.B_comp.fix(3.5e-8)  # membrane salt permeability
    m.fs.unit.permeate.pressure[0].fix(101325)  # permeate pressure

    # Set scaling factors for component mass flowrates.
    m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq','H2O'))
    m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq','TDS'))

    # Calculate scaling factors.
    calculate_scaling_factors(m.fs)

    # Initialize the model.
    m.fs.unit.initialize(optarg={'nlp_scaling_method': 'user-scaling'})


Configuration Options
---------------------
In addition to the core configuration options that are normally included for an IDAES process block, the RO unit model
contains options that account for concentration polarization, the mass transfer coefficient, and pressure drop.

.. csv-table::
   :header: "Keyword", "Value"

   "concentration_polarization_type", "ConcentrationPolarizationType.none"
   "concentration_polarization_type", "ConcentrationPolarizationType.fixed"
   "concentration_polarization_type", "ConcentrationPolarizationType.calculated"


Degrees of Freedom
------------------
This reverse osmosis (RO) unit model:
   * is a 0-dimensional model
   * supports a single liquid phase only
   * supports steady-state only

Model Structure
------------------
Feed-side: 0-D control volume
Permeate-side: state block

Variables
----------
.. csv-table::
   :header: "Description", "Symbol", "Variable", "Index", "Units"

   "Component mass fraction", ":math:`x_j`", "mass_frac_phase_comp", "[p, j]", ":math:`\text{dimensionless}`"
   "Mass density of seawater", ":math:`\rho`", "dens_mass_phase", "[p]", ":math:`\text{kg/}\text{m}^3`"
   "Mass density of pure water", ":math:`\rho_w`", "dens_mass_w_phase", "[p]", ":math:`\text{kg/}\text{m}^3`"
   "Phase volumetric flowrate", ":math:`Q_p`", "flow_vol_phase", "[p]", ":math:`\text{m}^3\text{/s}`"
   "Volumetric flowrate", ":math:`Q`", "flow_vol", "None", ":math:`\text{m}^3\text{/s}`"
   "Mass concentration", ":math:`C_j`", "conc_mass_phase_comp", "[p, j]", ":math:`\text{kg/}\text{m}^3`"
   "Dynamic viscosity", ":math:`\mu`", "visc_d_phase", "[p]", ":math:`\text{Pa}\cdotp\text{s}`"
   "Osmotic coefficient", ":math:`\phi`", "osm_coeff", "None", ":math:`\text{dimensionless}`"
   "Specific enthalpy", ":math:`\widehat{H}`", "enth_mass_phase", "[p]", ":math:`\text{J/kg}`"
   "Enthalpy flow", ":math:`H`", "enth_flow", "None", ":math:`\text{J/s}`"
   "Saturation pressure", ":math:`P_v`", "pressure_sat", "None", ":math:`\text{Pa}`"
   "Specific heat capacity", ":math:`c_p`", "cp_phase", "[p]", ":math:`\text{J/kg/K}`"
   "Thermal conductivity", ":math:`\kappa`", "therm_cond_phase", "[p]", ":math:`\text{W/m/K}`"
   "Latent heat of vaporization", ":math:`h_{vap}`", "dh_vap", "None", ":math:`\text{J/kg}`"


Constraints
-----------
.. csv-table::
   :header: "Description", "Symbol", "Variable", "Index", "Units"

   "Component mole flowrate", ":math:`N_j`", "flow_mol_phase_comp", "[p, j]", ":math:`\text{mole/s}`"
   "Component mole fraction", ":math:`y_j`", "mole_frac_phase_comp", "[p, j]", ":math:`\text{dimensionless}`" 
   "Molality", ":math:`Cm`", "molality_comp", "['TDS']", ":math:`\text{mole/kg}`"
   "Osmotic pressure", ":math:`\pi`", "pressure_osm", "None", ":math:`\text{Pa}`"

Relationships
-------------
.. csv-table::
   :header: "Description", "Equation"

   "Component mass fraction", ":math:`x_j = \frac{M_j}{\sum_{j} M_j}`"
   "Mass density", "Equation 8 in Sharqawy et al. (2010)"
   "Volumetric flowrate", ":math:`Q = \frac{\sum_{j} M_j}{\rho}`"
   "Mass concentration", ":math:`C_j = x_j \cdotp \rho`"
   "Dynamic viscosity", "Equations 22 and 23 in Sharqawy et al. (2010)"
   "Osmotic coefficient", "Equation 49 in Sharqawy et al. (2010)"
   "Specific enthalpy", "Equations 43 and 55 in Sharqawy et al. (2010)"
   "Enthalpy flow", ":math:`H = \sum_{j} M_j \cdotp \widehat{H}`"
   "Component mole flowrate", ":math:`N_j = \frac{M_j}{MW_j}`"
   "Component mole fraction", ":math:`y_j = \frac{N_j}{\sum_{j} N_j}`"
   "Molality", ":math:`Cm = \frac{x_{TDS}}{(1-x_{TDS}) \cdotp MW_{TDS}}`"
   "Osmotic pressure", ":math:`\pi = \phi \cdotp Cm \cdotp \rho_w \cdotp R \cdotp T` [See note below]"
   "Saturation pressure", "Equations 5 and 6 in Nayar et al. (2016)"
   "Specific heat capacity", "Equation 9 in Sharqawy et al. (2010)"
   "Thermal conductivity", "Equation 13 in Sharqawy et al. (2010)"
   "Latent heat of vaporization", "Equations 37 and 55 in Sharqawy et al. (2010)"


Class Documentation
-------------------

.. autoclass:: ReverseOsmosis0D
   :members:

.. autoclass:: ReverseOsmosisData
   :members:


