#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
from idaes_flowsheet_processor.api import FlowsheetInterface, FlowsheetCategory
from watertap.flowsheets.generic_desalination_train import generic_train
from pyomo.environ import units as pyunits
from idaes.core.solvers import get_solver


__author__ = "Alexander V. Dudchenko"


def export_to_ui():
    return FlowsheetInterface(
        name="Generic treatment train",
        do_export=export_variables,
        do_build=build_flowsheet,
        do_solve=solve_flowsheet,
        get_diagram=get_diagram,
        requires_idaes_solver=True,
        category=FlowsheetCategory.wastewater,
    )


def export_variables(flowsheet=None, exports=None, build_options=None, **kwargs):

    # print("model", model)
    fs = flowsheet  # model.fs
    # --- Input data ---
    # Feed conditions
    desalination_recovery = "Desalination recovery"
    water_management_costs = f"Water management costs"
    optional_inputs = f"Optional inputs"
    optional_costs = f"Optional cost"
    optional_inputs = f"Optional inputs"
    valorizer_inputs = "Valorization"
    desal_unit_costs = "Desalination unit cost"
    system_outputs = "System outputs"
    exports.add(
        obj=fs.costing.LCOW,
        name="LCOW Total",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m^3",
        rounding=3,
        description="LCOW",
        is_input=False,
        input_category="LCOW",
        is_output=True,
        output_category="LCOW",
    )

    for unit in fs.costing.LCOW_unit:
        exports.add(
            obj=fs.costing.LCOW_unit[unit],
            name=f'LCOW {unit.replace("_", " ")}',
            ui_units=fs.costing.base_currency / pyunits.m**3,
            display_units="$/m^3",
            rounding=3,
            description="LCOW",
            is_input=False,
            input_category="LCOW",
            is_output=True,
            output_category="LCOW",
            chart_type="stacked_bar_with_net",
            chart_group="process_cost",
        )
    for proc in fs.process_order:
        block = proc["process_block"]
        process_name = proc["process_name"]
        process_name_nice = proc["process_name"].replace("_", " ")

        if proc["process_type"] == "desalter":
            exports.add(
                obj=block.desalter.base_cost,
                name=f"{process_name_nice} base cost",
                ui_units=fs.costing.base_currency / pyunits.m**3,
                display_units="$/m^3",
                rounding=3,
                description=water_management_costs,
                is_input=True,
                input_category=water_management_costs,
                is_output=True,
                output_category=water_management_costs,
            )
            exports.add(
                obj=block.desalter.water_recovery,
                name=f"{process_name_nice} water recovery",
                ui_units=pyunits.dimensionless,
                display_units="%",
                rounding=3,
                description=desalination_recovery,
                is_input=True,
                input_category=desalination_recovery,
                is_output=True,
                output_category=desalination_recovery,
            )
            exports.add(
                obj=block.desalter.LCOW,
                name=f"{process_name_nice} unit cost",
                ui_units=fs.costing.base_currency / pyunits.m**3,
                display_units="$/m^3",
                rounding=3,
                description=desal_unit_costs,
                is_input=False,
                input_category=desal_unit_costs,
                is_output=True,
                output_category=desal_unit_costs,
            )

        if proc["process_type"] == "separator" or proc["process_type"] == "valorizer":
            if proc["process_type"] == "separator":
                exports.add(
                    obj=block.separator.base_cost,
                    name=f"{process_name_nice} base cost",
                    ui_units=fs.costing.base_currency / pyunits.m**3,
                    display_units="$/m^3",
                    rounding=3,
                    description=water_management_costs,
                    is_input=True,
                    input_category=water_management_costs,
                    is_output=True,
                    output_category=water_management_costs,
                )
            if proc["process_type"] == "valorizer":
                for (phase, ion), obj in fs.feed.properties[
                    0
                ].conc_mass_phase_comp.items():

                    if ion == "X":
                        exports.add(
                            obj=obj,
                            name="Feed {}".format(ion),
                            ui_units=pyunits.mg / pyunits.L,
                            display_units="mg/L",
                            rounding=0,
                            description="{} concentration".format(ion),
                            is_input=True,
                            input_category=valorizer_inputs,
                            is_output=True,
                            output_category=valorizer_inputs,
                        )
                exports.add(
                    obj=block.separator.mass_base_cost,
                    name=f"{process_name_nice} base cost",
                    ui_units=fs.costing.base_currency / pyunits.kg,
                    display_units="$/kg",
                    rounding=3,
                    description=valorizer_inputs,
                    is_input=True,
                    input_category=valorizer_inputs,
                    is_output=True,
                    output_category=valorizer_inputs,
                )

                for ion in block.separator.product_value.keys():
                    if ion != "TDS" and ion != "H2O":
                        exports.add(
                            obj=block.separator.component_removal_percent[ion],
                            name=f"Valorizer {ion} recovery %",
                            ui_units=pyunits.dimensionless,
                            display_units="%",
                            rounding=3,
                            description=valorizer_inputs,
                            is_input=True,
                            input_category=valorizer_inputs,
                            is_output=True,
                            output_category=valorizer_inputs,
                        )
                        exports.add(
                            obj=block.separator.product_value[ion],
                            name=f"Valorizer {ion} value",
                            ui_units=fs.costing.base_currency / pyunits.kg,
                            display_units="$/kg",
                            rounding=3,
                            description=valorizer_inputs,
                            is_input=True,
                            input_category=valorizer_inputs,
                            is_output=True,
                            output_category=valorizer_inputs,
                        )
                exports.add(
                    obj=block.separator.treatment_LCOW,
                    name=f"Valorizer unit cost",
                    ui_units=fs.costing.base_currency / pyunits.m**3,
                    display_units="$/m^3",
                    rounding=3,
                    description=valorizer_inputs,
                    is_input=False,
                    input_category=valorizer_inputs,
                    is_output=True,
                    output_category=valorizer_inputs,
                )
                exports.add(
                    obj=block.separator.ion_removal_LCOW,
                    name=f"Valorizer levelized revenue (LROW)",
                    ui_units=fs.costing.base_currency / pyunits.m**3,
                    display_units="$/m^3",
                    rounding=3,
                    description=valorizer_inputs,
                    is_input=False,
                    input_category=valorizer_inputs,
                    is_output=True,
                    output_category=valorizer_inputs,
                )
                for (phase, ion), obj in block.product_mass_flow_feed_basis.items():
                    if ion != "H2O" and ion != "TDS":
                        exports.add(
                            obj=obj,
                            name="Valorizer {} production".format(ion),
                            ui_units=pyunits.kg / pyunits.m**3,
                            display_units="kg/m^3  of feed",
                            rounding=3,
                            description="{} mass flow".format(ion),
                            is_input=False,
                            input_category=valorizer_inputs,
                            is_output=True,
                            output_category=valorizer_inputs,
                        )

    exports.add(
        obj=fs.feed.base_cost,
        name="Feed source cost",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m^3",
        rounding=3,
        description=optional_costs,
        is_input=True,
        input_category=optional_costs,
        is_output=True,
        output_category=optional_costs,
    )
    exports.add(
        obj=fs.product.base_cost,
        name="Product distribution cost",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m^3",
        rounding=3,
        description=optional_costs,
        is_input=True,
        input_category=optional_costs,
        is_output=True,
        output_category=optional_costs,
    )
    for proc in fs.process_order:
        block = proc["process_block"]
        process_name = proc["process_name"]
        process_name_nice = proc["process_name"].replace("_", " ")
        if proc["process_type"] == "desalter":
            exports.add(
                obj=block.desalter.recovery_cost,
                name=f"{process_name_nice} rate cost LCOW/WR",
                ui_units=fs.costing.base_currency / pyunits.m**3,
                display_units="$/m^3/%",
                rounding=3,
                description=optional_costs,
                is_input=True,
                input_category=optional_costs,
                is_output=True,
                output_category=optional_costs,
            )
            exports.add(
                obj=block.desalter.recovery_cost_offset,
                name=f"{process_name_nice} rate cost offset ",
                ui_units=pyunits.dimensionless,
                display_units="%",
                rounding=2,
                description=optional_costs,
                is_input=True,
                input_category=optional_costs,
                is_output=False,
                output_category=optional_costs,
            )
    exports.add(
        obj=fs.disposal.base_cost,
        name="Disposal cost",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m^3",
        rounding=3,
        description=water_management_costs,
        is_input=True,
        input_category=water_management_costs,
        is_output=True,
        output_category=water_management_costs,
    )

    for (phase, ion), obj in fs.feed.properties[0].conc_mass_phase_comp.items():
        if ion != "H2O" and ion != "X":
            exports.add(
                obj=obj,
                name="Feed {}".format(ion),
                ui_units=pyunits.mg / pyunits.L,
                display_units="mg/L",
                rounding=0,
                description="{} concentration".format(ion),
                is_input=True,
                input_category=optional_inputs,
                is_output=True,
                output_category=optional_inputs,
            )

    exports.add(
        obj=fs.water_recovery,
        name="System water recovery",
        ui_units=pyunits.dimensionless,
        display_units="%",
        rounding=3,
        description=optional_inputs,
        is_input=True,
        input_category=optional_inputs,
        is_output=True,
        output_category=optional_inputs,
    )

    for proc in fs.process_order:
        block = proc["process_block"]
        process_name = proc["process_name"]
        process_name_nice = proc["process_name"].replace("_", " ")
        if proc["process_type"] == "desalter":
            if process_name != "Desal_3":
                exports.add(
                    obj=block.desalter.brine_solids_concentration,
                    name=f"{process_name_nice} outlet TDS",
                    ui_units=pyunits.g / pyunits.L,
                    display_units="g/L",
                    rounding=3,
                    description=optional_inputs,
                    is_input=True,
                    input_category=optional_inputs,
                    is_output=True,
                    output_category=optional_inputs,
                )
            else:
                exports.add(
                    obj=block.desalter.brine_solids_concentration,
                    name=f"{process_name_nice} sludge TDS",
                    ui_units=pyunits.g / pyunits.L,
                    display_units="g/L",
                    rounding=3,
                    description=optional_inputs,
                    is_input=False,
                    input_category=optional_inputs,
                    is_output=True,
                    output_category=optional_inputs,
                )
                exports.add(
                    obj=block.desalter.brine_water_mass_percent,
                    name=f"{process_name_nice} sludge water content",
                    ui_units=pyunits.dimensionless,
                    display_units="%",
                    rounding=2,
                    description=optional_inputs,
                    is_input=True,
                    input_category=optional_inputs,
                    is_output=True,
                    output_category=optional_inputs,
                )

    for (phase, ion), obj in fs.disposal.properties[0].conc_mass_phase_comp.items():
        if ion != "H2O" and ion != "X":
            exports.add(
                obj=obj,
                name="Waste {} concentration".format(ion),
                ui_units=pyunits.g / pyunits.L,
                display_units="g/L",
                rounding=3,
                description=system_outputs,
                is_input=False,
                input_category=system_outputs,
                is_output=True,
                output_category=system_outputs,
            )
            if ion == "TDS":
                exports.add(
                    obj=fs.disposal_mass_flow_feed_basis[phase, ion],
                    name="Waste solids production",  # .format(ion),
                    ui_units=pyunits.kg / pyunits.m**3,
                    display_units="kg/m^3 of feed",
                    rounding=3,
                    description=system_outputs,
                    is_input=False,
                    input_category=system_outputs,
                    is_output=True,
                    output_category=system_outputs,
                )


def build_flowsheet(build_options=None, **kwargs):
    # build and solve initial flowsheet
    solver = get_solver()
    m = generic_train.build()
    generic_train.initialize(m, solver)
    return m


def get_diagram(build_options):
    return "fig_with_costs.png"


def solve_flowsheet(flowsheet):
    solver = get_solver()
    results = generic_train.solve(flowsheet, solver)
    return results
