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
from watertap.core.solvers import get_solver
from watertap.ui.fsapi import FlowsheetInterface
from watertap.examples.flowsheets.lsrro.lsrro import (
    build,
    set_operating_conditions,
    initialize,
    optimize_set_up,
    solve,
    ACase,
    BCase,
    ABTradeoff,
)
from pyomo.environ import units as pyunits


def export_to_ui():
    return FlowsheetInterface(
        name="LSRRO",
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
            "WaterRecovery": {
                "name": "WaterRecovery",
                "display_name": "Volumetric Water Recovery",
                "values_allowed": "float",
                "value": 0.5,  # default value
                "max_val": 1,  # optional
                "min_val": 0,  # optional
            },
            "FeedFlowRate": {
                "name": "FeedFlowRate",
                "display_name": "Feed Flow Rate (m3/s)",
                "values_allowed": "float",
                "value": 1e-3,  # default value
                "max_val": 1,  # optional
                "min_val": 1e-5,  # optional
            },
            "FeedNaClConcentration": {
                "name": "FeedNaClConcentration",
                "display_name": "Feed NaCl Concentration (kg/m3)",
                "values_allowed": "float",
                "value": 70,  # default value
                "max_val": 100,  # optional
                "min_val": 1,  # optional
            },
            "ROFiniteElements": {
                "name": "ROFiniteElements",
                "display_name": "RO Finite Elements",
                "values_allowed": "float",
                "value": 10,  # default value
                "max_val": 20,  # optional
                "min_val": 1,  # optional
            },
            "NaClSolubilityLimit": {
                "name": "NaClSolubilityLimit",
                "display_name": "NaCl Solubility Limit",
                "values_allowed": ["False", "True"],
                "value": "True",  # default value
            },
            "ConcentrationPolarization": {
                "name": "ConcentrationPolarization",
                "display_name": "Calculate Concentration Polzarization",
                "values_allowed": ["False", "True"],
                "value": "True",  # default value
            },
            "ROPressureDrop": {
                "name": "ROPressureDrop",
                "display_name": "Calculate Pressure Drop",
                "values_allowed": ["False", "True"],
                "value": "True",  # default value
            },
            "BMax": {
                "name": "BMax",
                "display_name": "Maximum NaCl Permeability",
                "values_allowed": "float",
                "value": 3.5e-6,  # default value
                "max_val": 1e-5,  # optional
                "min_val": 1e-7,  # optional
            },
        },
    )


def export_variables(flowsheet=None, exports=None, build_options=None, **kwargs):
    fs = flowsheet
    # --- Input data ---
    # Feed conditions
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
    exports.add(
        obj=fs.feed.pressure[0],
        name="Pressure",
        ui_units=pyunits.Pa,
        display_units="Pa",
        rounding=3,
        description="Inlet pressure",
        is_input=True,
        input_category="Feed",
        is_output=False,
    )
    exports.add(
        obj=fs.feed.temperature[0],
        name="Temperature",
        ui_units=pyunits.Pa,
        display_units="Pa",
        rounding=3,
        description="Inlet temperature",
        is_input=True,
        input_category="Feed",
        is_output=False,
    )

    # Unit model data, OARO
    for idx, pump in fs.PrimaryPumps.items():
        exports.add(
            obj=pump.control_volume.properties_out[0].pressure,
            name=f"Pump {idx} outlet pressure",
            ui_units=pyunits.Pa,
            display_units="Pa",
            rounding=2,
            description=f"Outlet pressure of pump {idx}",
            is_input=True,
            input_category="Primary Pumps",
            is_output=False,
        )
        exports.add(
            obj=pump.efficiency_pump[0],
            name=f"Pump {idx} efficiency",
            ui_units=pyunits.dimensionless,
            display_units="fraction",
            rounding=2,
            description=f"Efficiency of pump {idx}",
            is_input=True,
            input_category="Primary Pumps",
            is_output=False,
        )
    for idx, pump in fs.BoosterPumps.items():
        exports.add(
            obj=pump.efficiency_pump[0],
            name=f"Pump {idx} efficiency",
            ui_units=pyunits.dimensionless,
            display_units="fraction",
            rounding=2,
            description=f"Efficiency of pump {idx}",
            is_input=True,
            input_category="Booster Pumps",
            is_output=False,
        )
    for idx, stage in fs.ROUnits.items():
        exports.add(
            obj=stage.A_comp,
            name=f"RO {idx} water permeability coefficient",
            ui_units=pyunits.m / pyunits.Pa / pyunits.s,
            display_units="m/Pa/s",
            rounding=2,
            description=f"Membrane water permeability",
            is_input=True,
            input_category="Reverse Osmosis",
            is_output=False,
        )
        exports.add(
            obj=stage.B_comp,
            name=f"RO {idx} salt permeability coefficient",
            ui_units=pyunits.m / pyunits.s,
            display_units="m/s",
            rounding=2,
            description=f"Membrane salt permeability",
            is_input=True,
            input_category="Reverse Osmosis",
            is_output=False,
        )
        exports.add(
            obj=stage.area,
            name=f"RO {idx} area",
            ui_units=pyunits.m**2,
            display_units="m2",
            rounding=2,
            description=f"Membrane area",
            is_input=True,
            input_category="Reverse Osmosis",
            is_output=False,
        )
        exports.add(
            obj=stage.width,
            name=f"RO {idx} width",
            ui_units=pyunits.m,
            display_units="m",
            rounding=2,
            description=f"Membrane width",
            is_input=True,
            input_category="Reverse Osmosis",
            is_output=False,
        )
        exports.add(
            obj=stage.mixed_permeate[0].pressure,
            name=f"RO {idx} outlet pressure",
            ui_units=pyunits.Pa,
            display_units="Pa",
            rounding=2,
            description=f"Membrane outlet pressure",
            is_input=True,
            input_category="Reverse Osmosis",
            is_output=False,
        )
        if (
            build_options["ROPressureDrop"]
            or build_options["ConcentrationPolarization"]
        ):
            exports.add(
                obj=stage.feed_side.channel_height,
                name=f"RO {idx} channel height",
                ui_units=pyunits.m,
                display_units="m",
                rounding=2,
                description=f"Channel height in membrane stage",
                is_input=True,
                input_category="Reverse Osmosis",
                is_output=False,
            )
            exports.add(
                obj=stage.feed_side.spacer_porosity,
                name=f"RO {idx} space porosity",
                ui_units=pyunits.dimensionless,
                display_units="fraction",
                rounding=2,
                description=f"Spacer porosity in membrane stage",
                is_input=True,
                input_category="Reverse Osmosis",
                is_output=False,
            )
    for idx, erd in fs.EnergyRecoveryDevices.values():
        exports.add(
            obj=erd.efficiency_pump,
            name=f"ERD {idx} efficiency",
            ui_units=pyunits.dimensionless,
            display_units="fraction",
            rounding=2,
            description=f"Energy recovery device {idx} efficiency",
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
    exports.add(
        obj=fs.costing.reverse_osmosis.factor_membrane_replacement,
        name="Membrane replacement factor",
        ui_units=pyunits.yr**-1,
        display_units="1/yr",
        rounding=3,
        description="Membrane replacement factor [fraction of membrane replaced/year]",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.reverse_osmosis.membrane_cost,
        name="Membrane cost",
        ui_units=pyunits.USD_2018 / (pyunits.m**2),
        display_units="$/m2",
        rounding=3,
        description="Membrane cost",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.reverse_osmosis.high_pressure_membrane_cost,
        name="High pressure membrane cost",
        ui_units=pyunits.USD_2018 / (pyunits.m**2),
        display_units="$/m2",
        rounding=3,
        description="High pressure membrane cost",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.high_pressure_pump.cost,
        name="High pressure pump cost",
        ui_units=pyunits.USD_2018 / pyunits.watt,
        display_units="$/W",
        rounding=3,
        description="High pressure pump cost",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.energy_recovery_device.pressure_exchanger_cost,
        name="ERD pressure exchanger cost",
        ui_units=pyunits.USD_2018 / (pyunits.meter**3 / pyunits.hours),
        display_units="$/m3/hr",
        rounding=3,
        description="Pressure exchanger cost",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )

    # --- Output data ---
    # Feed
    exports.add(
        obj=fs.feed.properties[0].flow_vol_phase["Liq"],
        name="Volumetric flow rate",
        ui_units=pyunits.m**3 / pyunits.hour,
        display_units="m3/hr",
        rounding=2,
        description="Inlet volumetric flow rate",
        is_input=False,
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"],
        name="NaCl concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=2,
        description="Inlet NaCl concentration",
        is_input=False,
        is_output=True,
        output_category="Feed",
    )

    # Product
    exports.add(
        obj=fs.product.properties[0].flow_vol,
        name="Volumetric flow rate",
        ui_units=pyunits.m**3 / pyunits.hr,
        display_units="m3/h",
        rounding=2,
        description="Outlet product water volumetric flow rate",
        is_input=False,
        is_output=True,
        output_category="Product",
    )
    exports.add(
        obj=fs.product.properties[0].conc_mass_phase_comp["Liq", "NaCl"],
        name="NaCl concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=3,
        description="Outlet product water NaCl concentration",
        is_input=False,
        is_output=True,
        output_category="Product",
    )

    # Disposal
    exports.add(
        obj=fs.disposal.properties[0].flow_vol,
        name="Volumetric flow rate",
        ui_units=pyunits.m**3 / pyunits.hr,
        display_units="m3/h",
        rounding=2,
        description="Outlet product water volumetric flow rate",
        is_input=False,
        is_output=True,
        output_category="Disposal",
    )
    exports.add(
        obj=fs.disposal.properties[0].conc_mass_phase_comp["Liq", "NaCl"],
        name="NaCl concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=3,
        description="Outlet product water NaCl concentration",
        is_input=False,
        is_output=True,
        output_category="Disposal",
    )

    # System metrics
    exports.add(
        obj=fs.NumberOfStages,
        name="Number of stages",
        ui_units=pyunits.dimensionless,
        display_units=" ",
        rounding=0,
        description="Number of stages",
        is_input=False,
        is_output=True,
        output_category="System metrics",
    )
    exports.add(
        obj=fs.water_recovery,
        name="Water volumetric recovery",
        ui_units=pyunits.dimensionless,
        display_units=" ",
        rounding=2,
        description="System water volumetric recovery",
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
    exports.add(
        obj=fs.costing.specific_energy_consumption,
        name="Specific energy consumption",
        ui_units=pyunits.kWh / pyunits.m**3,
        display_units="kWh/m3 of product water",
        rounding=3,
        description="Specific energy consumption (SEC)",
        is_input=False,
        is_output=True,
        output_category="System metrics",
    )
    electricity_cost = (
        fs.costing.aggregate_flow_costs["electricity"]
        * fs.costing.utilization_factor
        / fs.costing.annual_water_production
    )
    exports.add(
        obj=electricity_cost,
        name="Electricity cost",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of product water",
        rounding=3,
        description="Electricity cost",
        is_input=False,
        is_output=True,
        output_category="System metrics",
    )
    total_area = sum(fs.ROUnits[i].area for i in fs.NonFinalStages) + fs.RO.area
    exports.add(
        obj=total_area,
        name="Total membrane area",
        ui_units=pyunits.m**2,
        display_units="m^2",
        rounding=2,
        description="Total membrane area",
        is_input=False,
        is_output=True,
        output_category="System metrics",
    )


def build_flowsheet(build_options=None, **kwargs):
    if build_options is not None:
        # get solver
        solver = get_solver()

        # build, set, and initialize
        m = build(
            number_of_stages=build_options["NumberOfStages"].value,
            has_NaCl_solubility_limit=build_options["NaClSolubilityLimit"].value,
            has_calculated_concentration_polarization=build_options[
                "ConcentrationPolarization"
            ].value,
            has_calculated_ro_pressure_drop=build_options["ROPressureDrop"].value,
            number_of_RO_finite_elements=build_options["ROFiniteElements"].value,
            B_max=build_options["BMax"].value,
        )
        set_operating_conditions(
            m,
            Cin=build_options["FeedNaClConcentration"].value,
            Qin=build_options["FeedFlowRate"].value,
        )

        initialize(m)

        solve(m, solver=solver)

        # Configuration options taken from user settings and run_lsrro_case
        optimize_set_up(
            m,
            set_default_bounds_on_module_dimensions=True,
            water_recovery=build_options["WaterRecovery"].value,
            Cbrine=None,
            A_case=ACase.optimize,
            B_case=BCase.optimize,
            AB_tradeoff=ABTradeoff.equality_constraint,
            A_value=None,
            permeate_quality_limit=500e-6,
            AB_gamma_factor=1,
            B_max=build_options["BMax"].value,
        )

        solve(m, raise_on_failure=False, tee=False, solver=solver)

    else:
        # get solver
        solver = get_solver()

        # build, set, and initialize
        m = build(
            number_of_stages=3,
            has_NaCl_solubility_limit=True,
            has_calculated_concentration_polarization=True,
            has_calculated_ro_pressure_drop=True,
            number_of_RO_finite_elements=10,
            B_max=None,
        )
        set_operating_conditions(m, Cin=None, Qin=None)

        initialize(m)
        solve(m, solver=solver)

        # Configuration options taken from def optimize_set_up
        optimize_set_up(
            m,
            set_default_bounds_on_module_dimensions=True,
            water_recovery=None,
            Cbrine=None,
            A_case=ACase.fixed,
            B_case=BCase.optimize,
            AB_tradeoff=ABTradeoff.none,
            A_value=None,
            permeate_quality_limit=None,
            AB_gamma_factor=None,
            B_max=None,
        )

        # display
        solve(m, solver=solver)

    return m


def solve_flowsheet(flowsheet=None):
    fs = flowsheet
    results = solve(fs)
    return results
