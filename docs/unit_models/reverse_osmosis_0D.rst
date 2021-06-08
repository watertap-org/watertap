Reverse Osmosis Unit (0D)
=========================
This reverse osmosis (RO) unit model
   * is 0-dimensional
   * supports a single liquid phase only
   * supports steady-state only
   * is based on the solution-diffusion model and film theory

.. index::
   pair: proteuslib.unit_models.reverse_osmosis_0D;reverse_osmosis_0D

.. currentmodule:: proteuslib.unit_models.reverse_osmosis_0D

Example A: Setting up the RO model
----------------------------------
The example below shows how to setup a simple RO unit model.

.. code-block:: python

    # Import concrete model from Pyomo
    from pyomo.environ import ConcreteModel
    # Import flowsheet block from IDAES core
    from idaes.core import FlowsheetBlock
    # Import NaCl property model
    import proteuslib.property_models.NaCl_prop_pack as props
    # Import utility tool for calculating scaling factors
    from idaes.core.util.scaling import calculate_scaling_factors
    # Import RO model
    from proteuslib.unit_models.reverse_osmosis_0D import ReverseOsmosis0D


    # Create a concrete model, flowsheet, and NaCl property parameter block.
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = props.NaClParameterBlock()

    # Add an RO unit to the flowsheet.
    m.fs.unit = ReverseOsmosis0D(default={"property_package": m.fs.properties})

    # Specify system variables.
    m.fs.unit.inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl'].fix(0.035)  # mass flow rate of NaCl
    m.fs.unit.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O'].fix(0.965)  # mass flow rate of water
    m.fs.unit.inlet.pressure[0].fix(50e5)  # feed pressure
    m.fs.unit.inlet.temperature[0].fix(298.15)  # feed temperature
    m.fs.unit.area.fix(50)  # membrane area
    m.fs.unit.A_comp.fix(4.2e-12)  # membrane water permeability
    m.fs.unit.B_comp.fix(3.5e-8)  # membrane salt permeability
    m.fs.unit.permeate.pressure[0].fix(101325)  # permeate pressure

    # Set scaling factors for component mass flowrates.
    m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
    m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'NaCl'))

    # Calculate scaling factors.
    calculate_scaling_factors(m)

    # Initialize the model.
    m.fs.unit.initialize(optarg={'nlp_scaling_method': 'user-scaling'})

Configuration Options
---------------------
In addition to the core configuration options that are normally included for an IDAES control volume, the RO unit model
contains options that account for concentration polarization, the mass transfer coefficient, and pressure drop.

**Default options denoted by** **

Configuration Keyword: ``concentration_polarization_type``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. csv-table::
    :header: "Configuration Options", "Description"

    "``ConcentrationPolarizationType.none`` **", "Simplifying assumption to ignore concentration polarization"
    "``ConcentrationPolarizationType.fixed``", "Specify an estimated value for the concentration polarization modulus"
    "``ConcentrationPolarizationType.calculated``", "Allow model to perform calculation of membrane-interface concentration"

Configuration Keyword: ``mass_transfer_coefficient``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. csv-table::
    :header: "Configuration Options", "Description"

    "``MassTransferCoefficient.none`` **", "Mass transfer coefficient not used in calculations"
    "``MassTransferCoefficient.fixed``", "Specify an estimated value for the mass transfer coefficient in the feed channel"
    "``MassTransferCoefficient.calculated``", "Allow model to perform calculation of mass transfer coefficient"


Configuration Keyword: ``pressure_change_type``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Notes:**
    * **to use any of the** ``pressure_change_type`` **options to account for pressure drop, the configuration keyword** ``has_pressure_change`` **must also be set to** ``True`` **.**
    * **if a value is specified for pressure drop, it should be negative.**

.. csv-table::
    :header: "Configuration Options", "Description"

    "``PressureChangeType.fixed_per_stage`` **", "Specify an estimated value for pressure drop across the membrane feed channel"
    "``PressureChangeType.fixed_per_unit_length``", "Specify an estimated value for pressure drop per unit length across the membrane feed channel"
    "``PressureChangeType.calculated``", "Allow model to perform calculation of pressure drop across the membrane feed channel"

Example B: Configure the RO model to account for concentration polarization and pressure drop
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    # Import concrete model from Pyomo
    from pyomo.environ import ConcreteModel
    # Import flowsheet block from IDAES core
    from idaes.core import FlowsheetBlock
    # Import NaCl property model
    import proteuslib.property_models.NaCl_prop_pack as props
    # Import utility tool for calculating scaling factors
    import idaes.core.util.scaling as calculate_scaling_factors
    #Import RO model and configuration classes
    from proteuslib.unit_models.reverse_osmosis_0D import (ReverseOsmosis0D,
                                                           ConcentrationPolarizationType,
                                                           MassTransferCoefficient,
                                                           PressureChangeType)

    # Create a concrete model, flowsheet, and NaCl property parameter block.
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = props.NaClParameterBlock()

    # Add an RO unit to the flowsheet and specify configuration options to calculate effects of
    # concentration polarization and pressure drop.
    m.fs.unit = ReverseOsmosis0D(default={"property_package": m.fs.properties,
                                          "has_pressure_change": True,
                                          "concentration_polarization_type": ConcentrationPolarizationType.calculated,
                                          "mass_transfer_coefficient": MassTransferCoefficient.calculated,
                                          "pressure_change_type": PressureChangeType.calculated})

Degrees of Freedom
------------------
Aside from the feed temperature, feed pressure, and component mass flow rates at the inlet, the RO model typically has
at least 4 degrees of freedom that should be fixed for the unit to be fully specified.

In Example A, the following variables were fixed, in addition to state variables at the inlet:
    * membrane water permeability, A
    * membrane salt permeability, B
    * permeate pressure
    * membrane area

The degrees of freedom will depend on which RO configuration options are selected. For example, setting
``has_pressure_change= True`` adds 1 degree of freedom. In this case, the pressure drop ``deltaP`` would need to be fixed
to eliminate that degree of freedom.

On the other hand, in Example B, configuring the RO unit to calculate concentration polarization effects, mass transfer
coefficient, and pressure drop would result in 3 more degrees of freedom than Example A. In this case, in addition to the
previously fixed variables, we typically fix the following variables to fully specify the unit:
    * feed-spacer porosity
    * feed-channel height
    * membrane length *or* membrane width

Model Structure
------------------
#TODO: Feed-side: 0-D control volume
Permeate-side: state block

Variables
----------
#TODO:

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
#TODO:

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


