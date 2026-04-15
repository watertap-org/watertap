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


def export_to_ui():
    return FlowsheetInterface(
        name="Multi-Stage Reverse Osmosis",
        do_export=export_variables,
        do_build=build_flowsheet,
        do_solve=solve_flowsheet,
        build_options={
            "n_stages": {
                "name": "n_stages",
                "display_name": "Number of stages",
                "values_allowed": "int",
                "value": 3,  # default value
                "max_val": 5,  # optional
                "min_val": 1,  # optional
            },
        #     "FeedFlowRate": {
        #         "name": "FeedFlowRate",
        #         "display_name": "Feed Flow Rate (m3/s)",
        #         "values_allowed": "float",
        #         "value": 1e-3,  # default value
        #         "max_val": 1,  # optional
        #         "min_val": 1e-5,  # optional
        #     },
        #     "FeedNaClConcentration": {
        #         "name": "FeedNaClConcentration",
        #         "display_name": "Feed NaCl Concentration (kg/m3)",
        #         "values_allowed": "float",
        #         "value": 70,  # default value
        #         "max_val": 100,  # optional
        #         "min_val": 1,  # optional
        #     },
        #     "ROFiniteElements": {
        #         "name": "ROFiniteElements",
        #         "display_name": "RO Finite Elements",
        #         "values_allowed": "float",
        #         "value": 10,  # default value
        #         "max_val": 20,  # optional
        #         "min_val": 1,  # optional
        #     },
        #     "NaClSolubilityLimit": {
        #         "name": "NaClSolubilityLimit",
        #         "display_name": "NaCl Solubility Limit",
        #         "values_allowed": ["False", "True"],
        #         "value": "True",  # default value
        #     },
        #     "ConcentrationPolarization": {
        #         "name": "ConcentrationPolarization",
        #         "display_name": "Calculate Concentration Polzarization",
        #         "values_allowed": ["False", "True"],
        #         "value": "True",  # default value
        #     },
        #     "ROPressureDrop": {
        #         "name": "ROPressureDrop",
        #         "display_name": "Calculate Pressure Drop",
        #         "values_allowed": ["False", "True"],
        #         "value": "True",  # default value
        #     },
        #     "ACase": {
        #         "name": "ACase",
        #         "display_name": "Water Permeability",
        #         "values_allowed": [
        #             x.name for x in [ACase.optimize, ACase.fixed, ACase.single_optimum]
        #         ],
        #         "value": ACase.fixed.name,  # default value
        #     },
        #     "BCase": {
        #         "name": "BCase",
        #         "display_name": "Salt Permeability",
        #         "values_allowed": [
        #             x.name for x in [BCase.optimize, BCase.single_optimum]
        #         ],
        #         "value": BCase.optimize.name,  # default value
        #     },
        #     "ABTradeoff": {
        #         "name": "ABTradeoff",
        #         "display_name": "Water & Salt Permeability Equality Constraints",
        #         "values_allowed": [
        #             x.name
        #             for x in [
        #                 ABTradeoff.inequality_constraint,
        #                 ABTradeoff.equality_constraint,
        #                 ABTradeoff.none,
        #             ]
        #         ],
        #         "value": ABTradeoff.none.name,  # default value
        #     },
        #     "BMax": {
        #         "name": "BMax",
        #         "display_name": "Maximum NaCl Permeability",
        #         "values_allowed": "float",
        #         "value": 3.5e-6,  # default value
        #         "max_val": 1e-5,  # optional
        #         "min_val": 1e-7,  # optional
        #     },
        },
    )


def export_variables(flowsheet=None, exports=None, build_options=None, **kwargs):
    fs = flowsheet
    comp = fs.properties.solute_set.at(1)
    # --- Input data ---
    # Feed conditions
    exports.add(
        obj=fs.feed.properties[0].flow_vol_phase["Liq"],
        name="Feed water volumetric flowrate",
        ui_units=pyunits.m**3 / pyunits.s,
        display_units="m^3/s",
        rounding=3,
        description="Inlet water volumetric flowrate",
        is_input=True,
        input_category="Feed",
        is_output=False,
    )
    exports.add(
        obj=fs.feed.properties[0].conc_mass_phase_comp["Liq", comp],
        name=f"Feed {comp} concentration",
        ui_units=pyunits.kg / pyunits.m**3,
        display_units="g/L",
        rounding=3,
        description=f"Inlet {comp} concentration",
        is_input=True,
        input_category="Feed",
        is_output=False,
    )
    # exports.add(
    #     obj=fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"],
    #     name="Feed water mass flowrate",
    #     ui_units=pyunits.kg / pyunits.s,
    #     display_units="kg/s",
    #     rounding=3,
    #     description="Inlet water mass flowrate",
    #     is_input=True,
    #     input_category="Feed",
    #     is_output=False,
    # )
    # exports.add(
    #     obj=fs.feed.properties[0].flow_mass_phase_comp["Liq", comp],
    #     name=f"Feed {comp} mass flowrate",
    #     ui_units=pyunits.kg / pyunits.s,
    #     display_units="kg/s",
    #     rounding=3,
    #     description=f"Inlet {comp} mass flowrate",
    #     is_input=True,
    #     input_category="Feed",
    #     is_output=False,
    # )

    # Unit model data, feed pump
    for n, stage in fs.stage.items():
        if stage.has_pump:
            exports.add(
                obj=stage.pump.efficiency_pump[0],
                name=f"Stage {n} pump efficiency",
                ui_units=pyunits.dimensionless,
                display_units="fraction",
                rounding=2,
                description=f"Stage {n} pump efficiency",
                is_input=True,
                input_category=f"Stage {n} Pump",
                is_output=True,
            )
            exports.add(
                obj=stage.pump.control_volume.properties_out[0].pressure,
                name=f"Stage {n} pump operating pressure",
                ui_units=pyunits.bar,
                display_units="bar",
                rounding=2,
                description=f"Stage {n} pump operating pressure",
                is_input=True,
                input_category=f"Stage {n} Pump",
                is_output=True,
            )
            # System metrics
        exports.add(
            obj=stage.RO.area,
            name=f"Stage {n} RO membrane area",
            ui_units=pyunits.m**2,
            display_units="m^2",
            rounding=2,
            description=f"Stage {n} RO membrane area",
            is_input=False,
            is_output=True,
            output_category="System metrics",
        )

        # Unit model data, RO
        exports.add(
            obj=stage.RO.A_comp[0, "H2O"],
            name=f"Stage {n} RO water permeability coefficient",
            ui_units=pyunits.L / pyunits.hr / pyunits.m**2 / pyunits.bar,
            display_units="LMH/bar",
            rounding=2,
            description=f"Stage {n} RO water permeability coefficient",
            is_input=True,
            input_category=f"Stage {n} Reverse Osmosis",
            is_output=False,
        )

        exports.add(
            obj=stage.RO.B_comp[0, comp],
            name=f"Stage {n} RO salt permeability coefficient",
            ui_units=pyunits.L / pyunits.hr / pyunits.m**2,
            display_units="LMH",
            rounding=2,
            description=f"Stage {n} RO salt permeability coefficient",
            is_input=True,
            input_category=f"Stage {n} Reverse Osmosis",
            is_output=False,
        )
        exports.add(
            obj=stage.RO.feed_side.channel_height,
            name=f"Stage {n} RO feed-side channel height",
            ui_units=pyunits.mm,
            display_units="mm",
            rounding=2,
            description=f"Stage {n} RO feed-side channel height",
            is_input=True,
            input_category=f"Stage {n} Reverse Osmosis",
            is_output=False,
        )
        exports.add(
            obj=stage.RO.feed_side.spacer_porosity,
            name=f"Stage {n} RO feed-side spacer porosity",
            ui_units=pyunits.dimensionless,
            display_units="fraction",
            rounding=2,
            description=f"Stage {n} RO feed-side spacer porosity",
            is_input=True,
            input_category=f"Stage {n} Reverse Osmosis",
            is_output=False,
        )
        exports.add(
            obj=stage.RO.permeate.pressure[0],
            name=f"Stage {n} RO permeate-side pressure",
            ui_units=pyunits.bar,
            display_units="bar",
            rounding=2,
            description=f"Stage {n} RO permeate-side pressure",
            is_input=True,
            input_category=f"Stage {n} Reverse Osmosis",
            is_output=False,
        )
        exports.add(
            obj=stage.RO.width,
            name=f"Stage {n} RO stage width",
            ui_units=pyunits.m,
            display_units="m",
            rounding=2,
            description=f"Stage {n} RO stage width",
            is_input=True,
            input_category=f"Stage {n} Reverse Osmosis",
            is_output=False,
        )
        exports.add(
            obj=stage.RO.flux,
            name=f"Stage {n} RO flux",
            ui_units=pyunits.liter / pyunits.hr / pyunits.m**2,
            display_units="LMH",
            rounding=2,
            description=f"Flux of Stage {n} RO unit",
            is_input=False,
            input_category=f"Stage {n} Reverse Osmosis",
            is_output=True,
            # output_category="System metrics",
        )
        exports.add(
            obj=stage.RO.recovery_mass_phase_comp[0, "Liq", "H2O"],
            name=f"Stage {n} RO water mass recovery",
            ui_units=pyunits.dimensionless,
            display_units="fraction",
            rounding=2,
            description=f"Water mass recovery of Stage {n} RO unit",
            is_input=False,
            input_category=f"Stage {n} Reverse Osmosis",
            is_output=True,
            output_category="System metrics",
        )

    # Unit model data, ERD
    if fs.add_erd:
        exports.add(
            obj=fs.ERD.efficiency_pump[0],
            name="ERD Pump efficiency",
            ui_units=pyunits.dimensionless,
            display_units="fraction",
            rounding=2,
            description="Efficiency of energy recovery device",
            is_input=True,
            input_category="Energy Recovery Device",
            is_output=False,
        )
        exports.add(
            obj=fs.ERD.control_volume.properties_out[0].pressure,
            name="ERD operating pressure",
            ui_units=pyunits.bar,
            display_units="bar",
            rounding=2,
            description="Operating pressure of energy recovery device",
            is_input=True,
            input_category="Energy Recovery Device",
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
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.TIC,
        name="Total Installed Cost",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=1,
        description="Total Installed Cost (TIC)",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.TPEC,
        name="Total Purchased Equipment Cost",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=1,
        description="Total Purchased Equipment Cost (TPEC)",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.total_investment_factor,
        name="Total investment factor",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Total investment factor [investment cost/equipment cost]",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.maintenance_labor_chemical_factor,
        name="Maintenance-labor-chemical factor",
        ui_units=1 / pyunits.year,
        display_units="fraction/year",
        rounding=2,
        description="Maintenance-labor-chemical factor [fraction of investment cost/year]",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.capital_recovery_factor,
        name="Capital annualization factor",
        ui_units=1 / pyunits.year,
        display_units="fraction/year",
        rounding=2,
        description="Capital annualization factor [fraction of investment cost/year]",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.electricity_cost,
        name="Electricity cost",
        ui_units=fs.costing.base_currency / pyunits.kWh,
        display_units="$/kWh",
        rounding=3,
        description="Electricity cost",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )

    # Feed
    # exports.add(
    #     obj=fs.feed.properties[0].flow_vol_phase["Liq"],
    #     name="Feed volumetric flow rate",
    #     ui_units=pyunits.m**3 / pyunits.hour,
    #     display_units="m3/hr",
    #     rounding=2,
    #     description="Inlet volumetric flow rate",
    #     is_input=False,
    #     is_output=True,
    #     output_category="Feed",
    # )
    # exports.add(
    #     obj=fs.feed.properties[0].conc_mass_phase_comp["Liq", comp],
    #     name="Feed NaCl concentration",
    #     ui_units=pyunits.g / pyunits.L,
    #     display_units="g/L",
    #     rounding=2,
    #     description="Inlet NaCl concentration",
    #     is_input=False,
    #     is_output=True,
    #     output_category="Feed",
    # )

    # Product
    exports.add(
        obj=fs.product.properties[0].flow_vol,
        name="Product volumetric flow rate",
        ui_units=pyunits.m**3 / pyunits.hr,
        display_units="m3/h",
        rounding=2,
        description="Product water volumetric flow rate",
        is_input=False,
        is_output=True,
        output_category="Product",
    )
    exports.add(
        obj=fs.product.properties[0].conc_mass_phase_comp["Liq", comp],
        name=f"Product {comp} concentration",
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
        name="Waste brine volumetric flow rate",
        ui_units=pyunits.m**3 / pyunits.hr,
        display_units="m3/h",
        rounding=2,
        description="Waste brine volumetric flow rate",
        is_input=False,
        is_output=True,
        output_category="Disposal",
    )
    exports.add(
        obj=fs.disposal.properties[0].conc_mass_phase_comp["Liq", comp],
        name=f"Waste brine {comp} concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=3,
        description=f"Waste brine {comp} concentration",
        is_input=False,
        is_output=True,
        output_category="Disposal",
    )

    exports.add(
        obj=fs.costing.SEC,
        name="Specific energy consumption",
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
        name="Levelized cost of water",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of product water",
        rounding=3,
        description="Levelized cost of water (LCOW)",
        is_input=False,
        is_output=True,
        output_category="System metrics",
    )


def build_flowsheet(build_options=None, **kwargs):
    # build and solve initial flowsheet
    # m = build()
    m = multistage_RO.build_n_stage_system(n_stages=build_options["n_stages"].value, add_erd=True)

    # the UI sets `capital_recovery_factor`, so unfix `wacc`
    m.fs.costing.wacc.unfix()
    m.fs.costing.capital_recovery_factor.fix()

    solver = get_solver()

    m = multistage_RO.run_n_stage_system(n_stages=build_options["n_stages"].value, add_erd=True)

    # build, set, and initialize
    # m = build(erd_type=erd_type)
    # set_operating_conditions(m)
    # initialize_system(m, solver=solver)

    # NOTE: GUI aims to solve simulation-based flowsheet
    # optimize_set_up(m)

    # display
    results = utils.solve(m, solver=solver)
    return m


def solve_flowsheet(flowsheet=None):
    fs = flowsheet
    results = utils.solve(model=fs.model())
    return results
