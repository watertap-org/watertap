How to setup a simple RO model
------------------------------------------------

The example below shows how to setup and initialize a simple RO unit model.

.. testsetup::

   # quiet idaes logs
   import idaes.logger as idaeslogger
   idaeslogger.getLogger('ideas.core').setLevel('CRITICAL')
   idaeslogger.getLogger('idaes.init').setLevel('CRITICAL')

.. testcode::
    
    from pyomo.environ import ConcreteModel
    from idaes.core import FlowsheetBlock
    from idaes.core.util.scaling import calculate_scaling_factors
    from watertap.property_models.NaCl_prop_pack import NaClParameterBlock
    from watertap.unit_models.reverse_osmosis_0D import ReverseOsmosis0D
    from watertap.unit_models.reverse_osmosis_0D import ConcentrationPolarizationType
    from watertap.unit_models.reverse_osmosis_0D import MassTransferCoefficient


    # Create a concrete model, flowsheet, and NaCl property parameter block.
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = NaClParameterBlock()

    # Add an RO unit to the flowsheet.
    m.fs.unit = ReverseOsmosis0D(
        property_package=m.fs.properties,
        concentration_polarization_type=ConcentrationPolarizationType.none,
        mass_transfer_coefficient=MassTransferCoefficient.none,
        has_pressure_change=False,
    )

    # Specify system variables.
    # mass flow rate of NaCl (kg/s)
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(0.035)
    # mass flow rate of water (kg/s)
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(0.965)
    m.fs.unit.inlet.pressure[0].fix(50e5)  # feed pressure (Pa)
    m.fs.unit.inlet.temperature[0].fix(298.15)  # feed temperature (K)
    m.fs.unit.area.fix(50)  # membrane area (m^2)
    m.fs.unit.A_comp.fix(4.2e-12)  # membrane water permeability (m/Pa/s)
    m.fs.unit.B_comp.fix(3.5e-8)  # membrane salt permeability (m/s)
    m.fs.unit.permeate.pressure[0].fix(101325)  # permeate pressure (Pa)

    # Set scaling factors for component mass flowrates.
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1e2, index=("Liq", "NaCl"))

    # Calculate scaling factors.
    calculate_scaling_factors(m)

    # Initialize the model.
    m.fs.unit.initialize()
