This is an example workflow for analysis using loop tool and psPlotKit. 

The workflow is split into 4 steps:

    1) Setup flowsheet
    2) Configure loop tool and analysis .yaml files
    3) Run analysis and generate output .h5 files
    4) Process data and generate results using psPlotKit

In this workflow, the loop tool and .yaml files create a tractable way of running complex analysis while leveraging advanced capabilities of parameter-sweep tool and .h5 storage format. 

Loop tool allows you to configure analysis that run distinct analysis under consistent assumptions and configurations. Examples could be:

    1) Running model across different designs (number of stages in LSRRO, or different pre-treatment options of flowsheet designs)

    2) Running specific cases (Feed water compositions, or process designs options)

The loop tool stores all data in a single .h5 file per simulation configuration. For example in a multi-stage analysis all stage results would be stored in the h5 file, in addition, every Variable, Parameter, and Expression value created during the build step will be stored in the file, so you can access all model data and results after simulation/optimization is ran. 

The loop tool allows you to run parallel simulations, accelerating analysis. This can be accomplished via two pathways:

    1) Set number_of_subprocesses > 1 to run parameter sweep simulations in parallel, this is useful when running large maps for a single case. Each processor will create a copy of the model and solve it for each parameter sweep options.

    2) Set num_loop_workers >1 to run each case in its own thread. This is useful if running on local machine and solving for many flowsheet designs/options (such as number_of_stages), in this case each scenario will be run in its own thread, and sequentially execute the requested parameter sweep. 

Finally, you can use psPlotKit to access and process the generated .h5 files and generate figures of interests such as maps and cost breakdowns. 

To follow the example workflow, go through following steps:

    1) Inspect analysis_setup.py - shows how to setup lsrro flowsheet for use with looptool and executes the analysis 

    2) Inspect lssro_stage_sweep.yaml - shows how to configure map and case sweeps across different number of stages

    3) Execute analysis_setup.py - this will run the simulations and generate an output folder that will contain two files:
    3.a) lsrro_stage_sweep_analysisType_map_sweep.h5 - contains results for map sweep
    3.b) lsrro_stage_sweep_analysisType_case_sweep.h5 - contains results for case sweep
    4) Inspect and execute map_plot.py - uses psPlotKit to plot cost for cost optimal LSRRO design and optimal number of stages
    5) Inspect and execute case_cost_breakdown_plot.py - uses psPlotKit to plot cost breakdown for two feed cases across water recoveries.

The above 3 python files contain additional details on the workflows and methods. 

