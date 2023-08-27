How to use loopTool to explore flowsheets
=========================================

.. index::
   pair: watertap.tools.analysis_tools.loop_tool;loop_tool

.. currentmodule:: watertap.tools.analysis_tools.loop_tool

The loopTool is a wrapper for the parameter sweep (PS) tool set, and is designed to simplify setting up parametric sweeps, enabling sweeping over discrete design choices, and providing structured data management.  

The loopTool uses the full features of PS tool set, including standard PS tool and differential PS tool, but brings in the ability to run iteratively over different build options, initialization options, and solve options that might be required for full flowsheet analysis. 

A common example is solving processes with multiple stages such as LSRRO flow sheet, or options where different design choices need to be evaluated such as different pressure exchanger types in RO. Other examples could be exploring different initialization guess, or setting up different solve constraints. All of these options can be run in nested configuration (e.g. for every stage simulate N number of build configurations).

The loopTool uses .h5 format for structured data management, and a companion tool dataImporter provides a simple interface to explore these files (coming soon). The h5 format enables saving data in a file using directories, such that each unique simulation run with parameter sweep is stored in its own directories, enabling one file to store a large number of different simulations. The loopTool does not support storing data in CSV format. 

Setting up loopTool
-----------------------------
The loopTool setup involves creating the following:

1) A flowsheet with build, initialize (optional), and optimize function
2) A .yaml file with simulation configuration 
3) A .py file that connects loopTool to the flowsheet, directories for where to save data, and yaml file, as well as other solve or PS options  

Setting up a flowsheet for use with loopTool
----------------------------------------------------
The loopTool requires a similar functional setup as PS tool kit. Here we will setup RO_with_energy_recovery.py example flowsheet for use with PS and loopTool. The RO_with_energy_recovery flowsheet has an option 
that selects the type of ERD used, by passing erd_type option at build (either a pressure_exchanger, pump_as_turbine, or None). The user likely would want to run a parameter sweep across these options, and we will setup the loopTool to do so. 

*The user will need to set up the following three functions:*

1. The *build_function* that builds the flowsheet (ro_build)- This function should build the flowsheet, and accept any kwargs for its configuration. Here we pass explicitly erd_type, but also include *kwargs* in case we want to pass in other options in the future. When loopTool runs, it will pass selected *kwargs* into this function before running the sweep. (This function ideally should not initialize the model but should build all relevant vars, parameters, and constraints)
2. The *initialization_function* (ro_init)- This function should initialize the model and set it up for optimization or simulation run.
   
   If the user enables update_sweep_params_before_init option (default: False), the Parameter Sweep tool (and loopTool) will update the parameters on the model tree that are being swept across before calling initialize as well as before solving the model. In the current example, we will use the default setting. 
   
   *Your initialization function should also be used to prepare the model for solving or optimization run, and thus should call on any additional functions to do so. In our example, a function “optimize_set_up” needs to be called before we run actual optimization. This function unfixes variables we want to optimize and fixes those that should be fixed. For PS/LoopTool, it does not matter if you are planning to do optimization (>0 DOFS) or simulation (0 DOFs) and as such you can setup the model to do either here in.*

3. The *optimize_function* this function will run your optimization or simulation, and should return the result. In general, you are free to include tests to see if result is optimal, but PS/loopTool will do that check before saving data from the solution. A common error is to have an optimization function that does not return a result, resulting in all of the solves failing during run, even if the model is successfully solved. 
   
   *The optimize function can include complex solution steps beyond just running the solve call. For example, a user can specify in the optimize function to first solve a model using a continuous approach, followed by a resolve with discrete choices. The user just needs to remember that PS tool will update the sweep parameters only before the call, but not verify that they were applied.*

An example of the functions setup for use with loopTool for RO_with_energy_recovery flowsheet example.  

.. code-block::

   import watertap.examples.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery as ro_erd
   def ro_build(erd_type=None, **kwargs):
      m = ro_erd.build(erd_type)
      return m

   def ro_init(m, solver=None, **kwargs):
      ro_erd.set_operating_conditions(
         m, water_recovery=0.5, over_pressure=0.3, solver=solver
      )
      ro_erd.initialize_system(m, solver=solver)
      ro_erd.optimize_set_up(m)

   def ro_solve(m, solver=None, **kwargs):
      result = solver.solve(m)
      return result

Setting up the .yaml configuration file 
--------------------------------------------------
The loopTool uses a yaml file to generate all sweep configurations, define variables to sweep over, and other options. 
This enables user to have multiple setup files that can be used to run simulations with out having to change any underlying code, except changing which yaml file the loopTool uses. 

**The loopTool .yaml configuration file accepts two types of options:**

   * *default options* - these configure simulation defaults that either define PS tool behavior, loopTool behavior, or default keys that are passed into build, initialize, or optimize functions unless they are overridden by loop options 
   * *loop options* - these define options that build iterative loops, and will override existing default values or be included with them. 

**default options:**

   * *initialize_before_sweep* –(default: False) option to force run initialization_function before every optimize_function call 
   * *update_sweep_params_before_init* –(default: False)  option if set to true, will result in parameters being swept over to be updated before the initlalize call. 
   * *number_of_subprocesses* –(default: 1)  Number of logical cores to use if running in parallel manner (not used wit MPI), default: 1
   * *parallel_back_end* – (default: MultiProcessing)  Parallel back end to use if number_of_subprocesses > 1. (refer to Parallel manager doc for further info)
   * *build_defaults* – default arguments that are passed for every sweep. For example in ro_erd example, shown here, this could be erd_type, but in other flowsheets, this could define any default options, from file locations, to number of stages etc. The defaults will be passed along with loop values, unless loop value overrides them. 
   * *init_defaults* - defaults for initialization_function, same behavior as build-defaults but used for initialization_function
   * *optimize_defaults* - defaults for optimize_function, same behavior as build or init defaults, but for optimize function. 


**loop options:**

   * *build_loop* - defines list of keywords arguments to loop over for the *build_function*
   * *init_loop* - defines list of keywords arguments to loop over for the *initialization_function*
   * *optimize_loop* - defines list of keywords arguments to loop over for the *optimize_function*
   * *cases* - defines specific cases to run for given build_loop, init_loop, or optimize_loop.  
   * *sweep_param_loop* - defines sweep parameters to loop over 
   * *diff_param_loop* - defines differential parameter sweep, the first set of sweep_params define differential sweep parameters, and second key in has to be *sweep_reference_params* - which defines the parameters to use for generation of references simulations, from which differential simulations are performed. 

**Defining sweep_param_loop and diff_param_loop:**

The sweep_param_loop and diff_param_loop define the parameters to perform sweeps over using PS tool kit, and support standard configurations in the following structure.

These configuration arguments are designed to be defined and set up in a yaml file. The general structure for yaml file structure should be as follows:

.. code-block::

   # for parametric sweeps
   # for LinearSamples, UniformSample, GeomSample, ReverseGeomSample, LatinHypercubeSample
   sweep_param_loop:
      sweep_param_name:  
         type: Any type supported by PS toolkit (LinearSample, UniformSample, etc)
         param: key on the flowsheet that can be found useing m.find_component (e.g. fs.costing.reverse_osmosis.membrane_cost)
         lower_limit: lower value for sampling
         upper_limit: upper value for sampling
         num_samples: number of samples to run
		min_num_samples: number of expected solved samples (this used in backup management. 
   
   # for NormalSample
   sweep_param_loop:
      sweep_param_name:  
         type: Any type supported by PS toolkit (LinearSample, UniformSample, etc)
         param: key on the flowsheet that can be found useing m.find_component (e.g. fs.costing.reverse_osmosis.membrane_cost)
         mean: mean value of normal sample
         sd: standard deviation
         num_samples: number of samples to run 

The sweep_param_loop parameters will be iterated over one by one, to sweep over multiple parameter at once you can define a group as follows:

.. code-block::

   # for parametric sweeps
   # for LinearSamples, UniformSample, GeomSample, ReverseGeomSample, LatinHypercubeSample
   sweep_param_loop:
      sweep_group_name:           
         sweep_param_name_1:
            type: Any type supported by PS toolkit (LinearSample, UniformSample, etc.)
            param: key on the flowsheet that can be found using m.find_component (e.g. fs.costing.reverse_osmosis.membrane_cost)
            lower_limit: lower value for sampling
            upper_limit: upper value for sampling
            num_samples: number of samples to run 
         sweep_param_name_2:
            type: Any type supported by PS toolkit (LinearSample, UniformSample, etc.)
            param: key on the flowsheet that can be found using m.find_component (e.g. fs.costing.reverse_osmosis.membrane_cost)
            lower_limit: lower value for sampling
            upper_limit: upper value for sampling
            num_samples: number of samples to run 
         sweep_param_name_N:
            type: Any type supported by PS toolkit (LinearSample, UniformSample, etc.)
            param: key on the flowsheet that can be found using m.find_component (e.g. fs.costing.reverse_osmosis.membrane_cost)
            lower_limit: lower value for sampling
            upper_limit: upper value for sampling
            num_samples: number of samples to run 


**Defining diff_param_loop**

The diff_param_loop defines a run for parameter_sweep_differential tool and requires the same parametrization. Namely the differential spec, and parameter sweep values. In the yaml file, when setting up diff_param_loop, the first set of parameters will define the differential spec, while values after *sweep_reference_params* will define the parameters for reference sweep.

.. code-block::

   # for differenatial sweeps
   diff_param_loop:
      # percentile diff_type
      diff_param_name:  
         diff_mode: perentile
         diff_sample_type: UniformSample
         nominal_lb: lower nomial value
         nominal_ub: upper nomial value
         num_samples: number of samples, if you have no differnece between relative_lb and relative_ub, set to 1
         param: key on the flowsheet that can be found useing m.find_component (e.g. fs.costing.reverse_osmosis.membrane_cost)
         relative_lb: lower percnetile bound to sample 
         relative_ub: upper percnetile bound to sample, can  be same as relatibe lb, to sample single value      
      # sum or product diff_type 
      diff_param_name:  
         diff_mode: sum or product
         diff_sample_type: UniformSample
         num_samples: number of samples, if you have no differnece between relative_lb and relative_ub, set to 1
         param: key on the flowsheet that can be found useing m.find_component (e.g. fs.costing.reverse_osmosis.membrane_cost)
         relative_lb: lower realtive value to sample 
         relative_ub: upper realtive value to sample
      sweep_reference_params: # normal sweep-param input, to generate refernece simulations      
         sweep_param_name:  
            type: Any type supported by PS toolkit (LinearSample, UniformSample, etc)
            param: key on the flowsheet that can be found useing m.find_component (e.g. fs.costing.reverse_osmosis.membrane_cost)
            mean: mean value of normal sample
            sd: standard deviation
            num_samples: number of samples to run 

**General .yaml configuration file structure:**

The general structure starts with analysis name, which will be used in the file name when saving, followed by default options and configurations, and finally by loops options. The general structure is shown below: 

.. code-block::

   analysis_name:
      reinitialize_before_sweep: (optional) False or True
      update_sweep_params_before_init: (optional) False or True
      number_of_subprocesses: (optional) int > 1
      parallel_back_end: (optional) MultiProcessing, RayIo, etc.
      build_defaults: (optional)
         default_a: value_a
      init_defaults: (optional)
         init_a: value_a
      optimize_defaults: (optional)
         opt_a: value_a 
      build_loop: 
         buld_arg:
            - value_a
            - vlaue_b
            - etc
            init_loop:
               init_arg: 
                  - value_a
                  - vlaue_b
                  - etc
               etc_loop:
                  etc_arg:
                     - etc_val
                  param_sweep_loop:
                     sweep_param_name_a:
                        sweep_param_configs
                     sweep_param_name_b:
                        sweep_param_configs
                     sweep_parm_etc:
                        etc

**Examples for RO with ERD**

Here we setup a simple run on RO erd flowsheet, requesting loopTool, to run PS tool over 2 RO erd_type configurations, for each erd_configuration, we run a linear sweep over membrane cost, a linear sweep over factor_membrane_replacment, and map sweep over NaCl loading and RO recovery, the map sweep will generate a mesh grid using NaCl loading and RO recovery, producing a 9 samples total. The example of map sweep can include as many or as few parameters as user desires.  

.. code-block::

   ro_erd_type_analysis:
      reinitialize_before_sweep: False # We don't need to reinit before each solve
      build_loop: # We are only gonna loop over each build function
         erd_type:
            - pressure_exchanger
            - pump_as_turbine
         sweep_param_loop
            membrane_cost:  # Runs over differnt membrnae costs
               type: LinearSample
               param: fs.costing.reverse_osmosis.membrane_cost
               lower_limit: 10
               upper_limit: 30
               num_samples: 3
            factor_membrane_replacement:  # Runs over membrane_replacment costs, generating 10 steps
               type: LinearSample
               param: fs.costing.reverse_osmosis.factor_membrane_replacement
               lower_limit: 0.1
               upper_limit: 0.2
               num_samples: 3 
            map_sweep: # Will run meshgrid sweep over feed_mass_nacl and ro_recovery, generating 100 samples
               feed_mass_nacl:  # Runs over salt mass flow rate only, generating 10 steps
                  type: LinearSample
                  param: fs.feed.properties[0].flow_mass_phase_comp[Liq,NaCl]
                  lower_limit: 0.03
                  upper_limit: 0.04
                  num_samples: 3 
               ro_recovery:  # Runs over ro recovery, generating 10 steps
                  type: LinearSample
                  param: fs.RO.recovery_mass_phase_comp[0,Liq,H2O]
                  lower_limit: 0.3
                  upper_limit: 0.5
                  num_samples: 3 

Example using cases, here we assume ro_build function takes in 2 options, a water type, and erd type.
The loop tool will ran a parameter sweep over each case.

.. code-block::

   ro_erd_analysis_simple:
      build_loop:
         cases:
            BGW_with_erd:
               water_type: BGW
               erd_type: pump_as_turbine
            BGW_without_ERD: 
               water_type: BGW
               erd_type: null # none in yaml is null
            SW_with_ERD:
               water_type: SW
               erd_type: pump_as_turbine 
         sweep_param_loop:
               membrane_cost:  # Runs over different membrane costs
                  type: LinearSample
                  param: fs.costing.reverse_osmosis.membrane_cost
                  lower_limit: 10
                  upper_limit: 30
                  num_samples: 3

Example when we don't do any build loops, init loops etc, this will simply run membrane_cost sweep. This is the same as using PS tool directly, but it gives simple data management structure if user wants to sweep over many variables in the model.

.. code-block::

   ro_erd_analysis_simple:
       sweep_param_loop:
            membrane_cost:  # Runs over differnt membrnae costs
               type: LinearSample
               param: fs.costing.reverse_osmosis.membrane_cost
               lower_limit: 10
               upper_limit: 30
               num_samples: 3

Example of setting up differential sweep. Briefly, this type of analysis can provide insight into how reducing membrane cost can reduce RO costs, even when exact membrane costs and replacement factors are unknown. 

.. code-block::

   ro_diff_analysis:
      reinitialize_before_sweep: False # We don't need to reinit before each solve
      build_loop: # We are only gonna loop over each build function
         erd_type:
            - pressure_exchanger
            - pump_as_turbine
         diff_param_loop: #the params below will be run iteratively, differentiating from the simulation created with sweep_reference_params
            membrane_cost:  # Runs over differnt membrnae costs
               diff_mode: percentile
               diff_sample_type: UniformSample
               param: fs.costing.reverse_osmosis.membrane_cost
               relative_lb: -0.01
               relative_ub: -0.01
               nominal_lb: 10
               nominal_ub: 30
               num_samples: 1
            factor_membrane_replacement:  # Runs over differnt membrnae costs
               diff_mode: percentile
               diff_sample_type: UniformSample
               param: fs.costing.reverse_osmosis.membrane_cost
               relative_lb: -0.01
               relative_ub: -0.01
               nominal_lb: 0.1
               nominal_ub: 0.2
               num_samples: 1  
            sweep_reference_params: # The parameters below will be used to generate reference design and will be run simultaneously (e.g. all of them will be changed to generate a new hypothetical design here, we only provide one parameter, membrane cost but could provide many others)
               membrane_cost: 
                  type: UniformSample
                  param: fs.costing.reverse_osmosis.membrane_cost
                  lower_limit: 10
                  upper_limit: 30                  
               num_samples: 10

Setting up the loopTool
-------------------------------------

To loopTool can be executed by passing in the flowsheet functions, yaml file, and save locations into the loopTool, as well as additional optional arguments:

   * loop_file (required): .yaml config file that contains iterative loops to run
   * solver (optional): solver to use in model, default uses watertap solver
   * build_function (required): function to build unit model
   * initialize_function (optional): function for initialization of the unit model
   * optimize_function (required): function for solving model
   * probe_function (optional): Function to probe if a solution should be attempted or not
   * save_name (required): name to use when saving the file with
   * save_dir (required): directory to save the file in
   * number_of_subprocesses (optional): user defined number of subprocesses to use for parallel run should be 1 or greater 
   * parallel_back_end (optional): backend to use for a parallel run if MPI is not used and number_of_subprocesses > 1
   * custom_do_param_sweep (optional): custom param function (refer to parameter sweep tool)
   * custom_do_param_sweep_kwargs (optional): custom parm kwargs (refer to parameter sweep tool)
   * execute_simulations (optional): sets if looptool should execute simulations upon setup, otherwise user can call build_run_dict, and run_simulations call manually
   * h5_backup (optional): Set location for backup file, if set to False, no backup will be created, otherwise the backup will be auto-created

Example of code for setting up our example of RO with ERD

.. code-block::

   # this imports the function created in example for RO_with_energy_recovery above
   import ro_erd as ro_setup
   # import the loopTool and utility function for getting working directory
   from watertap.tools.analysis_tools.loop_tool.loop_tool import loopTool, get_working_dir

   if __name__=='__main__': we will excute script here, required for safe execution of parallel scripts
      # We assume our yaml file is in same directory as the this script
      
      lp = loopTool(
        "ro_erd_sweep.yaml", # on of the yaml files we created in example above
        build_function=ro_setup.ro_build,  # our build_function
        initialize_function=ro_setup.ro_init,  # our initialize function
        optimize_function=ro_setup.ro_solve, # our solve function
        saving_dir=get_working_dir(), # this gets working directory for script so we can save files in same dirctory
        save_name="ro_with_erd", # this will be the name of this loop run
    )

The above code example will run the RO_erd flowsheet through loops specified with our yaml file. 

The loopTool will create an output folder in the directory as found by get_working_dir(). 

Upon succesfull run, there will be an output folder, with an h5 File, that will have a name of save_name_analysisType_analysis_name, save_name is specified in loopTool, while analysis_name is specified in yaml file (this means you can have multiple analysis types in same yaml file, and each analysis will be saved in its own file).

Output data structure and data protection
-----------------------------------------

The loopTool will store data in h5 file, with a structure similar to that of the yaml file being run. For example for the ro_erd_type_analysis example, the h5 file structure would be as follows:

.. code-block::

   ro_erd_type_analysis
   |-erd_type
      |-pressure_exchanger
      |  |-membrane_cost # contains standard PS h5 output
      |  |  |-outputs 
      |  |  |-solve_successful 
      |  |  |-sweep_params 
      |  |
      |  |-factor_membrane_replacement # contains standard PS h5 output
      |  |  |-outputs 
      |  |  |-solve_successful 
      |  |  |-sweep_params 
      |  |
      |  |-map_sweep # contains standard PS h5 output
      |     |-outputs 
      |     |-solve_successful 
      |     |-sweep_params 
      |
      |-pump_as_turbine
         |-membrane_cost # contains standard PS h5 output
         |  |-outputs 
         |  |-solve_successful 
         |  |-sweep_params 
         |
         |-factor_membrane_replacement # contains standard PS h5 output
         |  |-outputs 
         |  |-solve_successful 
         |  |-sweep_params 
         |
         |-map_sweep # contains standard PS h5 output
            |-outputs 
            |-solve_successful 
            |-sweep_params 

We can readely access the data for processing by useing h5py. 
To visually explore the h5 file, use HDFviewer ( https://www.hdfgroup.org/downloads/hdfview/ )
All the parameters will use their parm names or object keys as reference (e.g. if you set up yaml file to sweep over 'RO_recovery' for which param is 'fs.ro.ro_recovery', the file will store data for RO_recovery under key 'fs.ro.ro_recovery'.

.. code-block::

   import h5py
    h5file = h5py.File(
        "ro_with_erd_analysisType_ro_erd_type_analysis.h5", "r"
    )

    # get data for cost from  pressure_exchanger erd device, not the dir path strucure
    data = h5file[
        "ro_erd_type_analysis/erd_type/pressure_exchanger/membrane_cost/outputs/fs.costing.LCOW/value"
    ][()] # wil create a list

**Backup management**

The loopTool includes a naïve data management schema to prevent overwriting existing files and minimize simulation runs through the creation of backups from existing simulation file. When the loopTool starts it will check if a file with the same name already exists, if it does, it will rename that file to include a data and time (file_name_M_D-H_M-S.h5.bak).  

The backup file will be used to check if existing completed simulations exists before running a simulation. 
It will check if the backup file contains a complete simulation for the current run, (this only checks the number of successfully solved samples but does not check if the sweep_parameters match). If all of the simulations were successfully solved in backup file, it would copy over the data from the backup file into the new file. Otherwise, it will re-run the simulations. 
The user can specify the expected number of samples if it differs from num_samples by passing additional *expected_num_samples*, or if there is minimum number of samples expected using *min_num_samples*. 
User can also specify to force a re-run or to not re-run by passing *force_rerun*.
Example of use as follows:

.. code-block::


   # single paramter
   ro_erd_analysis_simple:
       sweep_param_loop:
            membrane_cost:  # Runs over differnt membrnae costs
               type: LinearSample
               param: fs.costing.reverse_osmosis.membrane_cost
               lower_limit: 10
               upper_limit: 30
               num_samples: 3               
               min_num_samples: int (if set, loopTool will only re-run if there is less samples in the backup file than the set value)
               expected_num_samples: int (if set, loopTool will only re-run if the expected_num_samples differs from the number of samples saved in the backup file)
               force_rerun: False or True (if set to True, will force to always rerun given parameter sweep, if False, will not re-run the parameter sweep) 
	         map_sweep: # Will run meshgrid sweep over feed_mass_nacl and ro_recovery, generating 100 samples
               feed_mass_nacl:  # Runs over salt mass flow rate only, generating 10 steps
                  type: LinearSample
                  param: fs.feed.properties[0].flow_mass_phase_comp[Liq,NaCl]
                  lower_limit: 0.03
                  upper_limit: 0.04
                  num_samples: 3 
               ro_recovery:  # Runs over ro recovery, generating 10 steps
                  type: LinearSample
                  param: fs.RO.recovery_mass_phase_comp[0,Liq,H2O]
                  lower_limit: 0.3
                  upper_limit: 0.5
                  num_samples: 3 
               min_num_samples: int (if set, loopTool will only re-run if there is less samples in the backup file than the set value)
               expected_num_samples: int (if set, loopTool will only re-run if the expected_num_samples differs from the number of samples saved in the backup file)
               force_rerun: False or True (if set to True, will force to always rerun given parameter sweep, if False, will not re-run the parameter sweep) 

This feature is designed to help with getting complete simulation sets, without needing to re-run all the looped options. For example, you might run a certain build loop, where only in one build option, some solutions failed. After fixing the reason, you can rerun the loopTool, and it will only rerun those build options that failed to solve. 

The user can disable the backup generation by setting h5_backup argument in loopTool to False, or, ideally, by deleting the existing h5 file if the simulation results in it are not needed. The loopTool will avoid overwriting existing files, and if a file with the same name exists, it will error out. 
 


