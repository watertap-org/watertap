.. _unit_model_utilities:

Unit Model Utilities
======================

This section contains documentation for utility functions that can be used to help build and/or initialize unit models.


Calculate Operating Pressure
-----------------------------

.. code-block:: python

   from watertap.core.util import calculate_operating_pressure

This function will estimate the operating pressure of a unit model based on the osmotic pressure of the inlet stream. It can accept the following arguments:

    * ``state_block``: The state block of the RO feed that has the non-pressure state variables set to desired values. Can accept state blocks from NaCl or seawater property models.
    * ``over_pressure_factor``: The amount of operating pressure above the brine osmotic pressure represented as a fraction (default=1.15)
    * ``water_recovery_mass``: The mass-based fraction of inlet H2O that becomes permeate (default=0.5)
    * ``salt_passage``: The mass-based fraction of inlet salt that becomes permeate (default=0)
    * ``solver``: Solver object to be used (default=None)

An separate model of the provided ``state_block`` is created and the operating pressure is calculated based on the provided ``water_recovery_mass`` and ``salt_passage`` values. 
The calculated operating pressure is returned in Pascals.

An example usage is provided below.

.. code-block:: python

   from pyomo.environ import ConcreteModel
   from idaes.core import FlowsheetBlock
   from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
   from watertap.unit_models import ReverseOsmosis0D
   from watertap.core.solvers import get_solver
   from watertap.core.util import calculate_operating_pressure

    # Create a state block with the desired inlet conditions
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SeawaterParameterBlock()

    m.fs.RO = ReverseOsmosis0D(
        property_package=m.fs.properties,
        concentration_polarization_type="none",
        mass_transfer_coefficient="none",
        has_pressure_change=False,
    )

    m.fs.RO.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(0.965)
    m.fs.RO.inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(0.035)
    m.fs.RO.inlet.temperature.fix(273 + 25)

    solver = get_solver()
    operating_pressure = calculate_operating_pressure(
        state_block=m.fs.RO.inlet,
        over_pressure_factor=1.15,
        water_recovery_mass=0.5,
        salt_passage=0,
        solver=solver,
    )

    m.fs.RO.inlet.pressure.fix(operating_pressure)
