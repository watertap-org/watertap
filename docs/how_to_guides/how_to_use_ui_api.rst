.. _howto_ui-api:

Add a flowsheet to the UI
==========================
.. py:currentmodule:: watertap.ui.fsapi

This API is intended for model developers who would like to connect their flowsheets to the UI.
Developers can select which variables to "export" to the UI for each component of the model, 
and provide extra metadata (display name, description) for them. For flowsheets, they should also 
specify how to build and solve the flowsheets.

For reference, see :class:`FlowsheetInterface` and :class:`FlowsheetExport` in the `watertap.ui.fsapi` module.

----

In some Python module, define the function ``export_to_ui``, which will look
similar to this::

    from watertap.ui.fsapi import FlowsheetInterface, FlowsheetCategory
    def export_to_ui():
        return FlowsheetInterface(
            name="OARO",
            do_export=export_variables,
            do_build=build_flowsheet,
            do_solve=solve_flowsheet,
            build_options={
                "NumberOfStages": {
                    "name": "NumberOfStages",
                    "display_name": "Number of stages",
                    "values_allowed": "int",
                    "value": 3,  # default value
                    "max_val": 8,  # optional
                    "min_val": 1,  # optional
                },
                "SystemRecovery": {
                    "name": "SystemRecovery",
                    "display_name": "System Mass Recovery",
                    "values_allowed": "float",
                    "value": 0.5,  # default value
                    "max_val": 1,  # optional
                    "min_val": 0,  # optional
                },
            },
        )
There are 3 required functions: 

1. ``do_export`` - This function defines the variables that will be displayed on the UI.

There are two ways to export the variables (which can be combined, if you really want to), using
the `exports` variable, which is an instance of :class:`FlowsheetExport`.
The first way is to use the Python API to call ``exports.add()`` (:meth:`FlowsheetExport.add`) for each variable. For example::

    def export_variables(flowsheet=None, exports=None, build_options=None, **kwargs):
        fs = flowsheet
        exports.add(
            obj=fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"],
            name="Water mass flowrate",
            ui_units=pyunits.kg / pyunits.s,
            display_units="kg/s",
            rounding=3,
            description="Inlet water mass flowrate",
            is_input=True,
            input_category="Feed",
            is_output=False,
        )
        exports.add(
            obj=fs.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"],
            name="NaCl mass flowrate",
            ui_units=pyunits.kg / pyunits.s,
            display_units="kg/s",
            rounding=3,
            description="Inlet NaCl mass flowrate",
            is_input=True,
            input_category="Feed",
            is_output=False,
        )

The second way to export variables is to create a comma-separated values (CSV) file with the same information, and
read that in with ``exports.from_csv()`` (:meth:`FlowsheetExport.from_csv`). For example::

    def export_variables(flowsheet=None, exports=None, build_options=None, **kwargs):
        exports.from_csv(file="oaro_exports.csv", flowsheet=flowsheet)

By default, the file is located in the same directory as the Python module.
The format of the file is documented in the :meth:`FlowsheetExport.from_csv` method, but it basically puts the
API keywords as columns in a table. For example, the CSV table for the API calls above would look like:

.. csv-table:: nf_exports.csv
    :header: "obj", "name", "descriptions", "ui_units", "display_units", "rounding", "is_input", "input_category", "is_output"

    "fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]","Water mass flowrate","Inlet water mass flowrate","units.kg / units.w","kg/s",3,true,"Feed",false
    "fs.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"]","NaCl mass flowrate","Inlet NaCl mass flowrate","units.kg / units.s","kg/s",3,true,"Feed",false

The raw text version is::

    "obj", "name", "descriptions", "ui_units", "display_units", "rounding", "is_input", "input_category", "is_output"
    "fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]","Water mass flowrate","Inlet water mass flowrate","units.kg / units.s","kg/s",3,true,"Feed",false
    "fs.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"]","NaCl mass flowrate","Inlet NaCl mass flowrate","units.kg / units.s","kg/s",3,true,"Feed",false

2. ``do_build`` - This function defines the build function for a flowsheet. See example below::

    from watertap.flowsheets.oaro.oaro_multi import (
        build,
        set_operating_conditions,
        initialize_system,
        optimize_set_up,
        solve,
    )
    def build_flowsheet():
        m = build(number_of_stages=number_of_stages, erd_type=erd_type)
        set_operating_conditions(m)
        initialize_system(
            m,
            number_of_stages,
            solvent_multiplier=0.5,
            solute_multiplier=0.7,
            solver=solver,
        )

        optimize_set_up(
            m, number_of_stages=number_of_stages, water_recovery=system_recovery
        )

        results = solve(m, solver=solver)
        assert_optimal_termination(results)

        return m


3. ``do_solve`` - This function defines the solve function for a flowsheet. See example below::

    from watertap.flowsheets.oaro.oaro_multi import solve
    def solve_flowsheet(flowsheet=None):
        fs = flowsheet
        results = solve(fs)
        return results

Additionally, there are optional parameters to assign a category, provide build options, 
and provide a diagram function among others. See additional examples below.

Build function using build options::

    def build_flowsheet(build_options=None, **kwargs):
        if build_options is not None:
            # get solver
            solver = get_solver()

            # build, set, and initialize
            m = build(
                number_of_stages=build_options["NumberOfStages"].value, erd_type=erd_type
            )
            set_operating_conditions(m)
            initialize_system(
                m,
                number_of_stages=build_options["NumberOfStages"].value,
                solvent_multiplier=0.5,
                solute_multiplier=0.7,
                solver=solver,
            )

            optimize_set_up(
                m,
                number_of_stages=build_options["NumberOfStages"].value,
                water_recovery=build_options["SystemRecovery"].value,
            )

            # display
            solve(m, solver=solver)
        else:
            # get solver
            solver = get_solver()

            # build, set, and initialize
            m = build(number_of_stages=3, erd_type=erd_type)
            set_operating_conditions(m)
            initialize_system(
                m,
                number_of_stages=3,
                solvent_multiplier=0.5,
                solute_multiplier=0.7,
                solver=solver,
            )

            optimize_set_up(
                m,
                number_of_stages=3,
                water_recovery=0.5,
            )

            # display
            solve(m, solver=solver)
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
            "RO = watertap.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery_ui",
            "OARO = watertap.flowsheets.oaro.oaro_multi_ui",
        ]


For a complete overview of all arguments, see :class:`FlowsheetInterface`.
