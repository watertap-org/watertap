.. _how_to_run_models_in_a_py_script:

How to run models in a Python script
====================================

In most circumstances, code involving WaterTAP models should run without issues in any way Python code can normally be run (assuming that the :ref:`installation procedure<install>` has been completed successfully).

However, when writing a script, it is recommended (and in some cases, required; see below for details) to structure the script in the following way:

#. Start from an empty Python source file, e.g. ``run_my_model.py``
#. At the beginning of the file, in the **global scope**, define the imports for external modules (and, if applicable, other externally defined objects).
   In this context, defining something at the global scope means that it is written at the leftmost level of the Python file, i.e. without any indentation.

   .. code-block::

      import pyomo.environ as pyo
      import idaes
      from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP

#. Next, also in the global scope, add any constants, functions, and/or classes:

   .. code-block::
   
       def my_function():
           return 42
   
   
       MY_CONSTANT = 1

#. Next, define a function (typically called ``main()``), where all other code (creating a model, variable, model initialization, solve, ...) will be defined.

   .. code-block::
   
      def main():
   
          m = pyo.ConcreteModel()
          m.x = pyo.Var()
          m.c = pyo.Constraint(expr=(0, m.x, 1))
          m.o = pyo.Objective(expr=m.x)
   
          nlp = PyomoNLP(m)

#. At the end of the file, add the following code snippet. This will cause the ``main()`` function to be called only when the Python file is invoked directly:

   .. code-block::

      if __name__ == "__main__":
          main()

#. Finally, save and run the Python file:

   .. code-block:: bash

      python run_my_model.py

Example: Python file with recommended structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. testcode::

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
   # import the solver
   from watertap.core.solvers import get_solver


   # Put all the model constructors, initialization, and solver in a separate function
   def main():
       # Create a concrete model, flowsheet, and NaCl property parameter block.
       m = ConcreteModel()
       m.fs = FlowsheetBlock(dynamic=False)
       m.fs.properties = props.NaClParameterBlock()
       # Add an RO unit to the flowsheet.
       m.fs.unit = ReverseOsmosis0D(property_package=m.fs.properties)

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

       # Get default watertap solver
       solver = get_solver()

       # Initialize the model passing default solver options
       m.fs.unit.initialize(optarg=solver.options)

       # Solve the model (using the tee=True option to display solver info)
       solver.solve(m, tee=True)


   # Call that function in the "__main__" for the script
   if __name__ == "__main__":
       main()

Example: the same code without recommended structure (may cause errors on Windows)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
   # import the solver
   from watertap.core.solvers import get_solver

   # Create a concrete model, flowsheet, and NaCl property parameter block.
   m = ConcreteModel()
   m.fs = FlowsheetBlock(dynamic=False)
   m.fs.properties = props.NaClParameterBlock()
   # Add an RO unit to the flowsheet.
   m.fs.unit = ReverseOsmosis0D(property_package=m.fs.properties)

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

   # Get default watertap solver
   solver = get_solver()

   # Initialize the model passing default solver options
   m.fs.unit.initialize(optarg=solver.options)

   # Solve the model (using the tee=True option to display solver info)
   solver.solve(m, tee=True)

If code other than imports and constant/function/class definitions is run in the global scope (i.e. not defined inside a function), it is likely to cause errors when run on Windows.
See `issue #387 <https://github.com/watertap-org/watertap/issues/387>`_ for more details.
