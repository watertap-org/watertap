#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

from pyomo.environ import units as pyunits
from idaes_flowsheet_processor.api import FlowsheetInterface
import watertap.flowsheets.multistage_RO.multistage_RO as multistage_RO
import watertap.flowsheets.multistage_RO.utils as utils
from watertap.core.solvers import get_solver

build = multistage_RO.build_n_stage_system
solve = utils.solve


def export_to_ui():
    return FlowsheetInterface(
        name="Multi-Stage Reverse Osmosis",
        do_export=export_variables,
        do_build=build_flowsheet,
        do_solve=solve_flowsheet,
        build_options={
            "recovery": {
                "name": "recovery",
                "display_name": "System Recovery (%)",
                "values_allowed": "float",
                "value": 50,  # default value
                "max_val": 95,  # optional
                "min_val": 20,  # optional
            },
            "flow_vol": {
                "name": "flow_vol",
                "display_name": "Feed Flow Rate (L/s)",
                "values_allowed": "float",
                "value": 1,  # default value
                "max_val": 10,  # optional
                "min_val": 1,  # optional
            },
            "salinity": {
                "name": "salinity",
                "display_name": "Feed Concentration (g/L)",
                "values_allowed": "float",
                "value": 35,  # default value
                "max_val": 100,  # optional
                "min_val": 5,  # optional
            },
            "n_stages": {
                "name": "n_stages",
                "display_name": "Number of stages",
                "values_allowed": "int",
                "value": 2,  # default value
                "max_val": 5,  # optional
                "min_val": 1,  # optional
            },
            "temperature": {
                "name": "temperature",
                "display_name": "Feed Temperature (°C)",
                "values_allowed": "float",
                "value": 25,  # default value
                "max_val": 100,  # optional
                "min_val": 1,  # optional
            },
            "max_perm_conc": {
                "name": "max_perm_conc",
                "display_name": "Maximum Product Concentration (g/L)",
                "values_allowed": "float",
                "value": 0.5,  # default value
                "max_val": 10,  # optional
                "min_val": 0.1,  # optional
            },
            "prop_pack": {
                "name": "prop_pack",
                "display_name": "Property Model",
                "values_allowed": ["NaCl", "Seawater"],
                "value": "Seawater",  # default value
            },
            "max_pressure": {
                "name": "max_pressure",
                "display_name": "Maximum Pumping Pressure (bar)",
                "values_allowed": "float",
                "value": 200,  # default value
                "max_val": 400,  # optional
                "min_val": 50,  # optional
            },
            "add_erd": {
                "name": "add_erd",
                "display_name": "Include ERD?",
                "values_allowed": ["False", "True"],
                "value": "True",  # default value
            },
            "A_comp": {
                "name": "A_comp",
                "display_name": "Water Permeability (L/m²/h/bar)",
                "values_allowed": "float",
                "value": 1.5,  # default value
                "max_val": 10,  # optional
                "min_val": 0,  # optional
            },
            "B_comp": {
                "name": "B_comp",
                "display_name": "Salt Permeability (L/m²/h)",
                "values_allowed": "float",
                "value": 0.1,  # default value
                "max_val": 1,  # optional
                "min_val": 0,  # optional
            },
            "stage_2_booster": {
                "name": "stage_2_booster",
                "display_name": "Include Stage 2 Booster?",
                "values_allowed": ["False", "True"],
                "value": "False",  # default value
            },
            "stage_3_booster": {
                "name": "stage_3_booster",
                "display_name": "Include Stage 3 Booster?",
                "values_allowed": ["False", "True"],
                "value": "False",  # default value
            },
        },
    )


def export_variables(flowsheet=None, exports=None, build_options=None, **kwargs):
    fs = flowsheet
    comp = fs.properties.solute_set.at(1)

    exports.add(
        obj=fs.costing.SEC,
        name="SEC",
        ui_units=pyunits.kWh / pyunits.m**3,
        display_units="kWh/m3 of product water",
        rounding=3,
        description="Specific energy consumption (SEC)",
        is_input=False,
        is_output=True,
        output_category="System metrics",
    )
    exports.add(
        obj=fs.costing.LCOW,
        name="LCOW",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="USD/m3 of product water",
        rounding=3,
        description="Levelized cost of water (LCOW)",
        is_input=False,
        is_output=True,
        output_category="System metrics",
    )
    exports.add(
        obj=fs.costing.total_capital_cost,
        name="CAPEX",
        ui_units=fs.costing.base_currency,
        display_units="USD",
        rounding=2,
        description="Total capital expenditure (CAPEX)",
        is_input=False,
        is_output=True,
        output_category="System metrics",
    )
    exports.add(
        obj=fs.costing.total_operating_cost,
        name="OPEX",
        ui_units=fs.costing.base_currency / pyunits.year,
        display_units="USD/year",
        rounding=2,
        description="Total operating expenditure (OPEX)",
        is_input=False,
        is_output=True,
        output_category="System metrics",
    )
    exports.add(
        obj=fs.total_power,
        name="Total Power Required",
        ui_units=pyunits.kW,
        display_units="kW",
        rounding=2,
        description="Total power required",
        is_input=False,
        is_output=True,
        output_category="System metrics",
    )

    # Feed
    exports.add(
        obj=fs.system_recovery,
        name="System Recovery",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=3,
        description="System volumetric recovery",
        is_input=True,
        input_category="System",
        is_output=True,
        output_category="System metrics",
    )
    exports.add(
        obj=fs.feed.properties[0].flow_vol_phase["Liq"],
        name="Feed volumetric flowrate",
        ui_units=pyunits.liter / pyunits.s,
        display_units="L/s",
        rounding=3,
        description="Inlet volumetric flowrate",
        is_input=False,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.properties[0].conc_mass_phase_comp["Liq", comp],
        name=f"Feed {comp} Concentration",
        ui_units=pyunits.g / pyunits.liter,
        display_units="g/L",
        rounding=3,
        description=f"Feed {comp} Concentration",
        is_input=False,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    # Product
    exports.add(
        obj=fs.product.properties[0].flow_vol,
        name="Product Flow Rate",
        ui_units=pyunits.liter / pyunits.s,
        display_units="L/s",
        rounding=2,
        description="Product water volumetric flow rate",
        is_input=False,
        is_output=True,
        output_category="Product",
    )
    exports.add(
        obj=fs.product.properties[0].conc_mass_phase_comp["Liq", comp],
        name=f"Product {comp} Concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=3,
        description=f"Product water {comp} concentration",
        is_input=False,
        is_output=True,
        output_category="Product",
    )
    # Disposal
    exports.add(
        obj=fs.disposal.properties[0].flow_vol,
        name="Brine Flow Rate",
        ui_units=pyunits.liter / pyunits.s,
        display_units="L/s",
        rounding=2,
        description="Brine volumetric flow rate",
        is_input=False,
        is_output=True,
        output_category="Disposal",
    )
    exports.add(
        obj=fs.disposal.properties[0].conc_mass_phase_comp["Liq", comp],
        name=f"Brine {comp} Concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=3,
        description=f"Brine {comp} concentration",
        is_input=False,
        is_output=True,
        output_category="Disposal",
    )

    # Unit model data, feed pump
    for n, stage in fs.stage.items():
        if stage.has_pump:
            exports.add(
                obj=stage.pump.control_volume.properties_out[0].pressure,
                name=f"Stage {n} Pump Operating Pressure",
                ui_units=pyunits.bar,
                display_units="bar",
                rounding=2,
                description=f"Pump operating pressure",
                is_input=True,
                input_category=f"Stage {n}",
                is_output=True,
                output_category=f"Stage {n}",
            )
            exports.add(
                obj=stage.pump.work_mechanical[0],
                name=f"Stage {n} Pump Power",
                ui_units=pyunits.kW,
                display_units="kW",
                rounding=2,
                description=f"Pump power",
                is_input=True,
                input_category=f"Stage {n}",
                is_output=True,
                output_category=f"Stage {n}",
            )
            exports.add(
                obj=stage.pump.efficiency_pump[0],
                name=f"Stage {n} Pump Efficiency",
                ui_units=pyunits.dimensionless,
                display_units="fraction",
                rounding=2,
                description=f"Pump efficiency",
                is_input=True,
                input_category=f"Stage {n}",
                is_output=True,
                output_category=f"Stage {n}",
            )
            # System metrics
        exports.add(
            obj=stage.RO.area,
            name=f"Stage {n} Area",
            ui_units=pyunits.m**2,
            display_units="m2",
            rounding=2,
            description=f"Membrane area",
            is_input=False,
            is_output=True,
            output_category=f"Stage {n}",
        )
        exports.add(
            obj=stage.RO.width,
            name=f"Stage {n} Width",
            ui_units=pyunits.m,
            display_units="m",
            rounding=2,
            description=f"Stage width",
            is_input=True,
            input_category=f"Stage {n}",
            is_output=True,
            output_category=f"Stage {n}",
        )
        exports.add(
            obj=stage.RO.length,
            name=f"Stage {n} Length",
            ui_units=pyunits.m,
            display_units="m",
            rounding=2,
            description=f"Stage length",
            is_input=True,
            input_category=f"Stage {n}",
            is_output=True,
            output_category=f"Stage {n}",
        )

        # Unit model data, RO
        exports.add(
            obj=stage.RO.A_comp[0, "H2O"],
            name=f"Stage {n} Water Permeability Coefficient",
            ui_units=pyunits.L / pyunits.hr / pyunits.m**2 / pyunits.bar,
            display_units="LMH/bar",
            rounding=2,
            description=f"Water permeability coefficient",
            is_input=True,
            input_category=f"Stage {n}",
            is_output=True,
            output_category=f"Stage {n}",
        )

        exports.add(
            obj=stage.RO.B_comp[0, comp],
            name=f"Stage {n} Salt Permeability Coefficient",
            ui_units=pyunits.L / pyunits.hr / pyunits.m**2,
            display_units="LMH",
            rounding=2,
            description=f"Salt permeability coefficient",
            is_input=True,
            input_category=f"Stage {n}",
            is_output=True,
            output_category=f"Stage {n}",
        )
        exports.add(
            obj=stage.RO.feed_side.channel_height,
            name=f"Stage {n} Channel Height",
            ui_units=pyunits.mm,
            display_units="mm",
            rounding=2,
            description=f"Feed-side channel height",
            is_input=True,
            input_category=f"Stage {n}",
            is_output=False,
        )
        exports.add(
            obj=stage.RO.feed_side.spacer_porosity,
            name=f"Stage {n} Spacer Porosity",
            ui_units=pyunits.dimensionless,
            display_units="fraction",
            rounding=2,
            description=f"Feed-side spacer porosity",
            is_input=True,
            input_category=f"Stage {n}",
            is_output=False,
        )
        exports.add(
            obj=stage.RO.permeate.pressure[0],
            name=f"Stage {n} Permeate Pressure",
            ui_units=pyunits.bar,
            display_units="bar",
            rounding=2,
            description=f"Permeate-side pressure",
            is_input=True,
            input_category=f"Stage {n}",
            is_output=False,
        )
        exports.add(
            obj=stage.average_flux_LMH,
            name=f"Stage {n} Average Flux",
            ui_units=pyunits.liter / pyunits.hr / pyunits.m**2,
            display_units="LMH",
            rounding=2,
            description=f"Average flux of stage {n} RO unit",
            is_input=False,
            input_category=f"Stage {n}",
            is_output=True,
            output_category=f"Stage {n}",
        )
        exports.add(
            obj=stage.recovery,
            name=f"Stage {n} Volumetric Recovery",
            ui_units=pyunits.dimensionless,
            display_units="fraction",
            rounding=2,
            description=f"Volumetric recovery",
            is_input=False,
            input_category=f"Stage {n}",
            is_output=True,
            output_category=f"Stage {n}",
        )
        exports.add(
            obj=stage.RO.recovery_mass_phase_comp[0, "Liq", "H2O"],
            name=f"Stage {n} Gravimetric Recovery",
            ui_units=pyunits.dimensionless,
            display_units="fraction",
            rounding=2,
            description=f"Gravimetric recovery",
            is_input=False,
            input_category=f"Stage {n}",
            is_output=True,
            output_category=f"Stage {n}",
        )

    # Unit model data, ERD
    if fs.add_erd:
        exports.add(
            obj=fs.ERD.efficiency_pump[0],
            name="ERD Efficiency",
            ui_units=pyunits.dimensionless,
            display_units="fraction",
            rounding=2,
            description="ERD efficiency",
            is_input=True,
            input_category="ERD",
            is_output=True,
            output_category="ERD",
        )
        exports.add(
            obj=fs.ERD.control_volume.properties_out[0].pressure,
            name="ERD Operating Pressure",
            ui_units=pyunits.bar,
            display_units="bar",
            rounding=2,
            description="ERD operating pressure",
            is_input=True,
            input_category="ERD",
            is_output=False,
        )

    # System costing
    exports.add(
        obj=fs.costing.utilization_factor,
        name="Utilization factor",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Utilization factor - [annual use hours/total hours in year]",
        is_input=True,
        input_category="Costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.TIC,
        name="Total Installed Cost Factor",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=1,
        description="Total Installed Cost (TIC) Factor",
        is_input=True,
        input_category="Costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.TPEC,
        name="Total Purchased Equipment Cost Factor",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=1,
        description="Total Purchased Equipment Cost (TPEC)",
        is_input=True,
        input_category="Costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.total_investment_factor,
        name="Total Investment Factor",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Total investment factor [investment cost/equipment cost]",
        is_input=True,
        input_category="Costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.maintenance_labor_chemical_factor,
        name="Maintenance-labor-chemical Factor",
        ui_units=1 / pyunits.year,
        display_units="fraction/year",
        rounding=2,
        description="Maintenance-labor-chemical factor [fraction of investment cost/year]",
        is_input=True,
        input_category="Costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.capital_recovery_factor,
        name="Capital Annualization Factor",
        ui_units=1 / pyunits.year,
        display_units="fraction/year",
        rounding=2,
        description="Capital annualization factor [fraction of investment cost/year]",
        is_input=True,
        input_category="Costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.electricity_cost,
        name="Electricity Cost",
        ui_units=fs.costing.base_currency / pyunits.kWh,
        display_units="$/kWh",
        rounding=3,
        description="Electricity cost",
        is_input=True,
        input_category="Costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.reverse_osmosis.membrane_cost,
        name="Membrane Unit Cost",
        ui_units=fs.costing.base_currency / pyunits.m**2,
        display_units="$/m2",
        rounding=3,
        description="Membrane cost",
        is_input=True,
        input_category="Costing",
        is_output=False,
    )


def build_flowsheet(build_options=None, **kwargs):

    solver = get_solver()

    if build_options is None:

        m = multistage_RO.run_n_stage_system()
        m = multistage_RO.set_system_recovery(m, 0.5)

    else:
        pump_dict = {1: True}
        for n in [2, 3]:
            if build_options[f"stage_{n}_booster"].value == "True":
                pump_dict[n] = True
            else:
                pump_dict[n] = False

        ro_op_dict = {
            "A_comp": build_options["A_comp"].value
            * pyunits.liter
            / pyunits.m**2
            / pyunits.hour
            / pyunits.bar,
            "B_comp": build_options["B_comp"].value
            * pyunits.liter
            / pyunits.m**2
            / pyunits.hour,
        }

        m = multistage_RO.run_n_stage_system(
            prop_pack=build_options["prop_pack"].value,
            flow_vol=build_options["flow_vol"].value,
            salinity=build_options["salinity"].value,
            temperature=build_options["temperature"].value,
            n_stages=build_options["n_stages"].value,
            add_erd=build_options["add_erd"].value,
            max_pressure=build_options["max_pressure"].value * 1e5,  # convert bar to Pa
            max_perm_conc=build_options["max_perm_conc"].value,
            pump_dict=pump_dict,
            ro_op_dict=ro_op_dict,
        )

        recovery = build_options["recovery"].value / 100  # convert % to fraction
        m = multistage_RO.set_system_recovery(m, recovery)

    _ = utils.solve(m, solver=solver)
    return m


def solve_flowsheet(flowsheet=None):
    fs = flowsheet
    results = utils.solve(model=fs.model())
    return results
