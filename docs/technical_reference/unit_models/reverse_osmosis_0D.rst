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
    * **if a value is specified for pressure change, it should be negative.**

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
This RO model consists of a ControlVolume0DBlock for the feed-channel of the membrane and a StateBlock for the permeate channel.

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Inlet/outlet", ":math:`x`", "['in', 'out']"
   "Phases", ":math:`p`", "['Liq']"
   "Components", ":math:`j`", "['H2O', 'NaCl']*"

* Solute depends on the imported property model; example shown here is for the NaCl property model.


Variables
----------

.. csv-table::
   :header: "Description", "Symbol", "Variable", "Index", "Units"

   "Solvent permeability coefficient", ":math:`A`", "A_comp", "[t, j]", ":math:`\text{m/Pa/s}`"
   "Solvent permeability coefficient", ":math:`B`", "B_comp", "[t, j]", ":math:`\text{m/s}`"
   "Mass density of pure water", ":math:`\rho_w`", "dens_solvent", "[p]", ":math:`\text{kg/}\text{m}^3`"
   "Mass flux across membrane", ":math:`J`", "flux_mass_io_phase_comp", "[t, x, p, j]", ":math:`\text{kg/s}\text{/m}^2`"
   "Membrane area", ":math:`A_m`", "area", "None", ":math:`\text{m}^2`"
   "Component recovery rate", ":math:`R_j`", "recovery_mass_phase_comp", "[t, p, j]", ":math:`\text{dimensionless}`"
   "Volumetric recovery rate", ":math:`R`", "recovery_vol_phase", "[t, p]", ":math:`\text{dimensionless}`"


Constraints
-----------
#TODO:

.. .. csv-table::
..   :header: "Description", "Equation"




Class Documentation
-------------------


.. autoclass:: ReverseOsmosis0D
   :members:
   :noindex:

.. autoclass:: ReverseOsmosisData
   :members:
   :noindex:


