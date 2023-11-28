.. _howto_ui-api:

How-to guide for the UI API
===========================
.. py:currentmodule:: watertap.ui.fsapi

.. contents:: Contents
    :depth: 2
    :local:

Overview
--------

.. image:: /_static/terminal-icon.png
    :height: 30px
    :align: left

This API is intended for model developers who would like to connect their flowsheets to the UI.
Developers can select which variables to "export" to the UI for each component of the model, 
and provide extra metadata (display name, description) for them. For flowsheets, they should also 
specify how to build and solve the flowsheets.

See also: :ref:`the UI flowsheet API reference section <ref_ui-fsapi>`.


Model Developers
----------------

.. image:: /_static/terminal-icon.png
    :height: 30px
    :align: left

Add an interface to your flowsheet
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In some Python module, define the function ``export_to_ui``, which will look
similar to this::

    from watertap.ui.fsapi import FlowsheetInterface, FlowsheetCategory
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

1. ``do_export`` - This function defines the variables that will be displayed on the UI. See example below::

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
