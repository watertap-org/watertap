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

.. code-block::

   import pyomo.environ as pyo
   import idaes
   from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP


   def my_function():
       return 42


   MY_CONSTANT = 1


   def main():

       m = pyo.ConcreteModel()
       m.x = pyo.Var()
       m.c = pyo.Constraint(expr=(0, m.x, 1))
       m.o = pyo.Objective(expr=m.x)

       nlp = PyomoNLP(m)


   if __name__ == "__main__":
       main()


Example: the same code without recommended structure (may cause errors on Windows)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block::

   import pyomo.environ as pyo
   import idaes
   from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP


   def my_function():
       return 42


   MY_CONSTANT = 1


   m = pyo.ConcreteModel()
   m.x = pyo.Var()
   m.c = pyo.Constraint(expr=(0, m.x, 1))
   m.o = pyo.Objective(expr=m.x)

   nlp = PyomoNLP(m)


If code other than imports and constant/function/class definitions is run in the global scope (i.e. not defined inside a function), it is likely to cause errors when run on Windows.
See `issue #387 <https://github.com/watertap-org/watertap/issues/387>`_ for more details.
