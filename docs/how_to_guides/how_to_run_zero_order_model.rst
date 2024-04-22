How to run a zero-order model
-----------------------------

The example script shown below is for the dual media filtration zero-order model. This example uses the default parameters defined in the watertap YAML database.

.. testcode::

   from pyomo.environ import ConcreteModel

   from idaes.core import FlowsheetBlock
   from watertap.core.solvers import get_solver

   from watertap.core.wt_database import Database
   from watertap.core.zero_order_properties import WaterParameterBlock
   from watertap.unit_models.zero_order import DualMediaFiltrationZO


   def main():

       # Create a Pyomo model and initialize the YAML database
       model = ConcreteModel()
       model.db = Database()

       # Create an IDAES flowsheet and define the solutes
       model.fs = FlowsheetBlock(dynamic=False)
       model.fs.params = WaterParameterBlock(solute_list=["nonvolatile_toc", "toc", "tss"])

       # Setup the zero-order model and define inlet flows
       model.fs.unit = DualMediaFiltrationZO(property_package=model.fs.params, database=model.db)
       model.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(10)
       model.fs.unit.inlet.flow_mass_comp[0, "nonvolatile_toc"].fix(1)
       model.fs.unit.inlet.flow_mass_comp[0, "toc"].fix(1)
       model.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(1)

       # Load default parameters from the YAML database
       model.fs.unit.load_parameters_from_database()

       # Access the solver and solve the model
       solver = get_solver()
       solver.solve(model)

       # Display a report of the results
       model.fs.unit.report()


   if __name__ == "__main__":
       main()

The output from running the example is shown below.

.. code-block:: text

   ====================================================================================
   Unit : fs.unit                                                             Time: 0.0
   ------------------------------------------------------------------------------------
       Unit Performance

       Variables:

       Key                              : Value      : Units         : Fixed : Bounds
                     Electricity Demand :     2373.6 :          watt : False : (0, None)
                  Electricity Intensity : 1.8259e+05 :        pascal :  True : (None, None)
       Solute Removal [nonvolatile_toc] :    0.20000 : dimensionless :  True : (0, None)
                   Solute Removal [toc] :    0.20000 : dimensionless :  True : (0, None)
                   Solute Removal [tss] :    0.97000 : dimensionless :  True : (0, None)
                         Water Recovery :    0.99000 : dimensionless :  True : (0.0, 1.0000001)

   ------------------------------------------------------------------------------------
       Stream Table
                                                  Units            Inlet   Treated  Byproduct
       Volumetric Flowrate                   meter ** 3 / second 0.013000 0.011530 0.0014700
       Mass Concentration H2O              kilogram / meter ** 3   769.23   858.63    68.027
       Mass Concentration nonvolatile_toc  kilogram / meter ** 3   76.923   69.384    136.05
       Mass Concentration toc              kilogram / meter ** 3   76.923   69.384    136.05
       Mass Concentration tss              kilogram / meter ** 3   76.923   2.6019    659.86
   ====================================================================================

The zero-order models rely on default model parameter values specified in YAML files. Users should supply their own values, if possible, instead of relying on the default parameter values. To use a local YAML database file, define the folder path to the file using the :code:`dbpath` parameter for the :code:`Database()` class. The line below defines the current working directory as the path for the YAML file.

.. code-block::

   model.db = Database(dbpath=".")

When using a local YAML file, the file name must be the same as the model name without the ZO letters. It must also use snake case style; for example, a YAML file for the :code:`DualMediaFiltrationZO()` model must be named :code:`dual_media_filtration.yaml`. See below for the contents of this YAML file.

.. code-block:: yaml

   # Contents of YAML file named dual_media_filtration.yaml
   # This defines parameter values for the DualMediaFiltrationZO() zero-order model

   default:
     energy_electric_flow_vol_inlet:
       value: 0.050718512
       units: kWh/m^3
     capital_cost:
       basis: flow_vol
       cost_factor: None
       reference_state:
         value: 4732.0
         units: m^3/hr
       capital_a_parameter:
         value: 12.17829669e6
         units: USD_2014
       capital_b_parameter:
         value: 0.5862
         units: dimensionless
     recovery_frac_mass_H2O:
       value: 0.99
       units: dimensionless
       reference:
     default_removal_frac_mass_comp:
       value: 0
       units: dimensionless
     removal_frac_mass_comp:
       nonvolatile_toc:
         value: 0.2
         units: dimensionless
         constituent_longform: Nonvolatile TOC
       toc:
         value: 0.2
         units: dimensionless
         constituent_longform: Total Organic Carbon (TOC)
       tss:
         value: 0.97
         units: dimensionless
         constituent_longform: Total Suspended Solids (TSS)
