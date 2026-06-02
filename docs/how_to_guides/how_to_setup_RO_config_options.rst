How to use configuration options in the RO model
------------------------------------------------
In addition to the core configuration options that are normally included for an IDAES control volume, the RO unit model
contains options that account for concentration polarization, the mass transfer coefficient, and pressure drop.

Example: Configure the RO model to account for concentration polarization and pressure drop
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. testcode::

    # Import RO model and configuration classes
    from watertap.unit_models.reverse_osmosis_0D import ReverseOsmosis0D
    from watertap.unit_models.reverse_osmosis_0D import ConcentrationPolarizationType
    from watertap.unit_models.reverse_osmosis_0D import MassTransferCoefficient
    from watertap.unit_models.reverse_osmosis_0D import PressureChangeType

    # Import concrete model from Pyomo
    from pyomo.environ import ConcreteModel

    # Import flowsheet block from IDAES core
    from idaes.core import FlowsheetBlock

    # Import NaCl property model
    import watertap.property_models.NaCl_prop_pack as props

    # Create a concrete model, flowsheet, and NaCl property parameter block.
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.NaClParameterBlock()

    # Add an RO unit to the flowsheet and specify configuration options to calculate effects of
    # concentration polarization and pressure drop.
    m.fs.unit = ReverseOsmosis0D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        pressure_change_type=PressureChangeType.calculated,
    )
