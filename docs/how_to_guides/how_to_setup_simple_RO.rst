How to setup a simple RO model
------------------------------------------------

The example below shows how to setup and initialize a simple RO unit model.

.. testsetup::

   # quiet idaes logs
   import idaes.logger as idaeslogger
   idaeslogger.getLogger('ideas.core').setLevel('CRITICAL')
   idaeslogger.getLogger('idaes.init').setLevel('CRITICAL')

.. code-block::

    # Import concrete model from Pyomo
    from pyomo.environ import ConcreteModel
    # Import flowsheet block from IDAES core
    from idaes.core import FlowsheetBlock
    # Import NaCl property model
    import watertap.property_models.NaCl_prop_pack as props
    # Import utility tool for calculating scaling factors
    from idaes.core.util.scaling import calculate_scaling_factors
    # Import RO model
    from watertap.unit_models.reverse_osmosis_0D import ReverseOsmosis0D


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
