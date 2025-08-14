.. _howto_ui-api:

Add a flowsheet to the UI
==========================
.. py:currentmodule:: idaes_flowsheet_processor.api

This API is intended for model developers who would like to connect their flowsheets to the UI.
Developers can select which variables to "export" to the UI for each component of the model, 
and provide extra metadata (display name, description) for them. For flowsheets, they should also 
specify how to build and solve the flowsheets.

For reference, see :class:`FlowsheetInterface` and :class:`FlowsheetExport` in the `idaes_flowsheet_processor.api` module.

----

In some Python module, define the function ``export_to_ui``, which will look
similar to this::

    from idaes_flowsheet_processor.api import FlowsheetInterface, FlowsheetCategory
    def export_to_ui():
        return FlowsheetInterface(
            name="NF-DSPM-DE",
            do_export=export_variables,
            do_build=build_flowsheet,
            do_solve=solve_flowsheet,
            get_diagram=get_diagram,
            requires_idaes_solver=True,
            category=FlowsheetCategory.wastewater,
            build_options={
                "Bypass": {
                    "name": "bypass option",
                    "display_name": "With Bypass",
                    "values_allowed": ['false', 'true'],
                    "value": "false"
                }
            }
        )

There are 3 required functions:

1. ``do_export`` - This function defines the variables that will be displayed on the UI.

There are two ways to export the variables (which can be combined, if you really want to), using
the `exports` variable, which is an instance of :class:`FlowsheetExport`.
The first way is to use the Python API to call ``exports.add()`` (:meth:`FlowsheetExport.add`) for each variable. For example::

    def export_variables(flowsheet=None, exports=None, build_options=None, **kwargs):
        fs = flowsheet
        exports.add(
            obj=fs.feed.properties[0].flow_vol_phase["Liq"],
            name="Volumetric flow rate",
            ui_units=pyunits.L / pyunits.hr,
            display_units="L/h",
            rounding=2,
            description="Inlet volumetric flow rate",
            is_input=True,
            input_category="Feed",
            is_output=False,
            output_category="Feed",
        )
        exports.add(
            obj=fs.NF.pump.outlet.pressure[0],
            name="NF pump pressure",
            ui_units=pyunits.bar,
            display_units="bar",
            rounding=2,
            description="NF pump pressure",
            is_input=True,
            input_category="NF design",
            is_output=True,
            output_category="NF design",
        )
        exports.add(
            obj=fs.NF.product.properties[0].flow_vol_phase["Liq"],
            name="NF product volume flow",
            ui_units=pyunits.L / pyunits.hr,
            display_units="L/h",
            rounding=2,
            description="NF design",
            is_input=False,
            input_category="NF design",
            is_output=True,
            output_category="NF design",
        )

The second way to export variables is to create a comma-separated values (CSV) file with the same information, and
read that in with ``exports.from_csv()`` (:meth:`FlowsheetExport.from_csv`). For example::

    def export_variables(flowsheet=None, exports=None, build_options=None, **kwargs):
        exports.from_csv(file="nf_exports.csv", flowsheet=flowsheet)

By default, the file is located in the same directory as the Python module.
The format of the file is documented in the :meth:`FlowsheetExport.from_csv` method, but it basically puts the
API keywords as columns in a table. For example, the CSV table for the API calls above would look like:

.. csv-table:: nf_exports.csv
    :header: "obj", "name", "descriptions", "ui_units", "display_units", "rounding", "is_input", "input_category", "is_output", "output_category"

    "fs.feed.properties[0].flow_vol_phase['Liq']","Volumetric flow rate","Volumetric flow rate","units.L / units.hr","L/h",2,true,"Feed",false,""
    "fs.NF.pump.outlet.pressure[0]","NF pump pressure","Nanofiltration pump outlet pressure","units.bar","bar",true,"NF design",true,"NF design"
    "fs.NF.product.properties[0].flow_vol_phase['Liq']","NF product volume flow rate","Nanofiltration product volume flow rate","units.L / units.hr","L/h",2,false,"",true,"NF design"

The raw text version is::

    "obj", "name", "descriptions", "ui_units", "display_units", "rounding", "is_input", "input_category", "is_output", "output_category"
    "fs.feed.properties[0].flow_vol_phase['Liq']","Volumetric flow rate","Volumetric flow rate","units.L / units.hr","L/h",2,true,"Feed",false,""
    "fs.NF.pump.outlet.pressure[0]","NF pump pressure","Nanofiltration pump outlet pressure","units.bar","bar",true,"NF design",true,"NF design"
    "fs.NF.product.properties[0].flow_vol_phase['Liq']","NF product volume flow rate","Nanofiltration product volume flow rate","units.L / units.hr","L/h",2,false,"",true,"NF design"

2. ``do_build`` - This function defines the build function for a flowsheet. See example below::

    from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.metab.metab import (
        build,
        set_operating_conditions,
        initialize_system,
        solve,
        add_costing,
        adjust_default_parameters,
    )
    def build_flowsheet():
        # build and solve initial flowsheet
        m = build()

        set_operating_conditions(m)
        assert_degrees_of_freedom(m, 0)
        assert_units_consistent(m)

        initialize_system(m)

        results = solve(m)
        assert_optimal_termination(results)

        add_costing(m)
        assert_degrees_of_freedom(m, 0)
        m.fs.costing.initialize()

        adjust_default_parameters(m)

        results = solve(m)
        assert_optimal_termination(results)
        return m


3. ``do_solve`` - This function defines the solve function for a flowsheet. See example below::

    from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.metab.metab import solve
    def solve_flowsheet(flowsheet=None):
        fs = flowsheet
        results = solve(fs)
        return results

Additionally, there are optional parameters to assign a category, provide build options,
and provide a diagram function among others. See additional examples below.

Build function using build options::

    def build_flowsheet(build_options=None, **kwargs):
        # build and solve initial flowsheet
        if build_options is not None:
            if build_options["Bypass"].value == "true": #build with bypass
                solver = get_solver()
                m = nf_with_bypass.build()
                nf_with_bypass.initialize(m, solver)
                nf_with_bypass.unfix_opt_vars(m)
            else: # build without bypass
                solver = get_solver()
                m = nf.build()
                nf.initialize(m, solver)
                nf.add_objective(m)
                nf.unfix_opt_vars(m)
        else: # build without bypass
            solver = get_solver()
            m = nf.build()
            nf.initialize(m, solver)
            nf.add_objective(m)
            nf.unfix_opt_vars(m)
        return m

Custom diagram function::

    def get_diagram(build_options):
        if build_options["Bypass"].value == "true":
            return "nf_with_bypass_ui.png"
        else:
            return "nf_ui.png"

Enable UI to discover flowsheet - In order for the UI to discover a flowsheet, an
entrypoint must be defined in setup.py with the path to the export file. For examples, see below::

    entry_points={
        "watertap.flowsheets": [
            "nf = watertap.examples.flowsheets.nf_dspmde.nf_ui",
            "metab = watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.metab.metab_ui",
        ]


For a complete overview of all arguments, see :class:`FlowsheetInterface`.

Testing flowsheet interfaces
-------------------------------

.. note::
    The following requires the WaterTAP testing dependencies to be installed. This is done by default in a developer environment, or can be installed manually by running ``pip install "watertap[testing]"``.

To verify that the flowsheet interfaces developed and distributed with WaterTAP function correctly, it is possible to use pytest to run a series of standardized tests against each interface. The ``idaes-flowsheets`` pytest plugin that enables this functionality is installed together with the rest of the WaterTAP testing dependencies. However, it must be activated by providing additional command-line flags when invoking pytest.

.. code-block:: shell

   # run tests for all flowsheet interfaces registered under the `watertap.flowsheets` entry point group
   pytest --idaes-flowsheets --entry-point-group watertap.flowsheets

   # run tests for one or more importable Python modules (useful when developing a new flowsheet interface as it also supports modules not yet registered as an entry point)
   pytest --idaes-flowsheets --modules watertap.flowsheets.gac.gac_ui watertap.flowsheets.mvc.mvc_single_stage_ui

.. note::
   By default, the ``idaes-flowsheets`` pytest plugin will collect and run its tests *in addition* to the normally discovered pytest tests (e.g. test functions defined in ``test_*.py`` files throughout the current working directory). To disable normal (Python) pytest collection and only run ``idaes-flowsheets`` tests, use the ``-p no:python`` flag::

        pytest --idaes-flowsheets --modules watertap.flowsheets.gac.gac_ui watertap.flowsheets.mvc.mvc_single_stage_ui -p no:python