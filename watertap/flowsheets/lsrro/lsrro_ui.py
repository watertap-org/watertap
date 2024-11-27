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
from idaes_flowsheet_processor.api import FlowsheetInterface
from watertap.flowsheets.lsrro.lsrro import (
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
            "ACase": {
                "name": "ACase",
                "display_name": "Water Permeability",
                "values_allowed": [
                    x.name for x in [ACase.optimize, ACase.fixed, ACase.single_optimum]
                ],
                "value": ACase.fixed.name,  # default value
            },
            "BCase": {
                "name": "BCase",
                "display_name": "Salt Permeability",
                "values_allowed": [
                    x.name for x in [BCase.optimize, BCase.single_optimum]
                ],
                "value": BCase.optimize.name,  # default value
            },
            "ABTradeoff": {
                "name": "ABTradeoff",
                "display_name": "Water & Salt Permeability Equality Constraints",
                "values_allowed": [
                    x.name
                    for x in [
                        ABTradeoff.inequality_constraint,
                        ABTradeoff.equality_constraint,
                        ABTradeoff.none,
                    ]
                ],
                "value": ABTradeoff.none.name,  # default value
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
        description="Inlet salinity",
        is_input=False,
        is_output=True,
        output_category="Feed",
    )

    # Unit model data, OARO
    for idx, pump in fs.PrimaryPumps.items():
        exports.add(
            obj=pump.control_volume.properties_out[0].pressure,
            name=f"Pump {idx} outlet pressure",
            ui_units=pyunits.bar,
            display_units="bar",
            rounding=1,
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
        name = f"Stage {idx}"
        exports.add(
            obj=stage.A_comp[0, "H2O"],
            name=f"Stage {idx} water permeability coefficient",
            ui_units=pyunits.L / pyunits.hr / pyunits.m**2 / pyunits.bar,
            display_units="LMH/bar",
            rounding=2,
            description=f"Membrane water permeability",
            is_input=True,
            input_category=f"{name}",
            is_output=False,
        )
        exports.add(
            obj=stage.B_comp[0, "NaCl"],
            name=f"Stage {idx} salt permeability coefficient",
            ui_units=pyunits.L / pyunits.hr / pyunits.m**2,
            display_units="LMH",
            rounding=2,
            description=f"Membrane salt permeability",
            is_input=True,
            input_category=f"{name}",
            is_output=False,
        )
        exports.add(
            obj=stage.area,
            name=f"Stage {idx} membrane area",
            ui_units=pyunits.m**2,
            display_units="m2",
            rounding=2,
            description=f"Membrane area",
            is_input=True,
            input_category=f"{name}",
            is_output=True,
            output_category="Membrane area",
        )
        exports.add(
            obj=stage.width,
            name=f"Stage {idx} width",
            ui_units=pyunits.m,
            display_units="m",
            rounding=2,
            description=f"Membrane width",
            is_input=True,
            input_category=f"{name}",
            is_output=False,
        )
        exports.add(
            obj=stage.mixed_permeate[0].pressure,
            name=f"Stage {idx} permeate pressure",
            ui_units=pyunits.bar,
            display_units="bar",
            rounding=1,
            description=f"Membrane permeate pressure",
            is_input=True,
            input_category=f"{name}",
            is_output=False,
        )
        if (
            build_options["ROPressureDrop"]
            or build_options["ConcentrationPolarization"]
        ):
            exports.add(
                obj=stage.feed_side.channel_height,
                name=f"Stage {idx} channel height",
                ui_units=pyunits.mm,
                display_units="mm",
                rounding=1,
                description=f"Channel height in membrane stage",
                is_input=True,
                input_category=f"{name}",
                is_output=False,
            )
            exports.add(
                obj=stage.feed_side.spacer_porosity,
                name=f"Stage {idx} space porosity",
                ui_units=pyunits.dimensionless,
                display_units="fraction",
                rounding=2,
                description=f"Spacer porosity in membrane stage",
                is_input=True,
                input_category=f"{name}",
                is_output=False,
            )
    for idx, erd in fs.EnergyRecoveryDevices.items():
        exports.add(
            obj=erd.efficiency_pump[0],
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

    # System constraints
    exports.add(
        obj=fs.water_recovery,
        name="Water recovery",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Water recovery",
        is_input=True,
        input_category="System constraints",
        is_output=True,
        output_category="System constraints",
    )
    exports.add(
        obj=fs.ro_max_pressure,
        name="Maximum allowable RO pressure",
        ui_units=pyunits.bar,
        display_units="bar",
        rounding=2,
        description="Maximum allowable pressure for RO",
        is_input=True,
        input_category="System constraints",
        is_output=True,
        output_category="System constraints",
    )
    exports.add(
        obj=fs.lsrro_max_pressure,
        name="Maximum allowable LSRRO pressure",
        ui_units=pyunits.bar,
        display_units="bar",
        rounding=2,
        description="Maximum allowable pressure for LSRRO",
        is_input=True,
        input_category="System constraints",
        is_output=True,
        output_category="System constraints",
    )

    # --- Output data ---
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
    # Primary Pumps
    exports.add(
        obj=fs.total_pump_work,
        name="Pump work",
        ui_units=pyunits.kW,
        display_units="kW",
        rounding=3,
        description="Total work by primary pumps",
        is_input=False,
        is_output=True,
        output_category="Primary Pumps",
    )
    exports.add(
        obj=fs.recovered_pump_work,
        name="Pump work recovered",
        ui_units=pyunits.kW,
        display_units="kW",
        rounding=3,
        description="Total work recovered by primary pumps",
        is_input=False,
        is_output=True,
        output_category="Primary Pumps",
    )
    exports.add(
        obj=fs.net_pump_work,
        name="Net pump work",
        ui_units=pyunits.kW,
        display_units="kW",
        rounding=3,
        description="Net pump work by primary pumps",
        is_input=False,
        is_output=True,
        output_category="Primary Pumps",
    )
    energy_recovery = -fs.recovered_pump_work / fs.total_pump_work * 100
    exports.add(
        obj=energy_recovery,
        name="Energy recovery",
        ui_units=pyunits.dimensionless,
        display_units="%",
        rounding=3,
        description="Energy recovery by primary pumps",
        is_input=False,
        is_output=True,
        output_category="Primary Pumps",
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
        obj=fs.water_recovery * 100,
        name="Water volumetric recovery",
        ui_units=pyunits.dimensionless,
        display_units="%",
        rounding=2,
        description="System water volumetric recovery",
        is_input=False,
        is_output=True,
        output_category="System metrics",
    )
    exports.add(
        obj=fs.mass_water_recovery * 100,
        name="Mass water recovery rate",
        ui_units=pyunits.dimensionless,
        display_units="%",
        rounding=2,
        description="System water mass recovery",
        is_input=False,
        is_output=True,
        output_category="System metrics",
    )
    exports.add(
        obj=fs.system_salt_rejection * 100,
        name="Salt rejection",
        ui_units=pyunits.dimensionless,
        display_units="%",
        rounding=2,
        description="System salt rejection",
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
    total_area = sum(fs.ROUnits[i].area for i in fs.NonFinalStages)
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
    exports.add(
        obj=fs.annual_feed,
        name="Annual feed flow",
        ui_units=pyunits.m**3 / pyunits.yr,
        display_units="m3/yr",
        rounding=3,
        description="Annual feed water flow rate",
        is_input=False,
        is_output=True,
        output_category="System metrics",
    )
    exports.add(
        obj=fs.costing.annual_water_production,
        name="Annual water production",
        ui_units=pyunits.m**3 / pyunits.yr,
        display_units="m3/yr",
        rounding=3,
        description="Annual feed water production",
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
        output_category="LCOW breakdown",
    )
    total_capex = (
        fs.costing.total_capital_cost
        * fs.costing.capital_recovery_factor
        / fs.costing.annual_water_production
    )
    exports.add(
        obj=total_capex,
        name="Total CAPEX LCOW",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of product water",
        rounding=3,
        description="Total capital expenses levelized cost of water",
        is_input=False,
        is_output=True,
        output_category="LCOW breakdown",
    )
    exports.add(
        obj=fs.costing.primary_pump_capex_lcow,
        name="Primary pump CAPEX LCOW",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of product water",
        rounding=3,
        description="Primary pump capital expenses levelized cost of water",
        is_input=False,
        is_output=True,
        output_category="LCOW breakdown",
    )
    exports.add(
        obj=fs.costing.booster_pump_capex_lcow,
        name="Booster pump CAPEX LCOW",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of product water",
        rounding=3,
        description="Booster pump capital expenses levelized cost of water",
        is_input=False,
        is_output=True,
        output_category="LCOW breakdown",
    )
    exports.add(
        obj=fs.costing.erd_capex_lcow,
        name="Energy recovery device CAPEX LCOW",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of product water",
        rounding=3,
        description="ERD capital expenses levelized cost of water",
        is_input=False,
        is_output=True,
        output_category="LCOW breakdown",
    )
    exports.add(
        obj=fs.costing.membrane_capex_lcow,
        name="Membrane CAPEX LCOW",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of product water",
        rounding=3,
        description="Membrane capital expenses levelized cost of water",
        is_input=False,
        is_output=True,
        output_category="LCOW breakdown",
    )
    exports.add(
        obj=fs.costing.indirect_capex_lcow,
        name="Indirect CAPEX LCOW",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of product water",
        rounding=3,
        description="Indirect capital expenses levelized cost of water",
        is_input=False,
        is_output=True,
        output_category="LCOW breakdown",
    )
    exports.add(
        obj=fs.costing.electricity_lcow,
        name="Electricity LCOW",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of product water",
        rounding=3,
        description="Electricity levelized cost of water",
        is_input=False,
        is_output=True,
        output_category="LCOW breakdown",
    )
    exports.add(
        obj=fs.costing.membrane_replacement_lcow,
        name="Membrane replacement LCOW",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of product water",
        rounding=3,
        description="Membrane replacement levelized cost of water",
        is_input=False,
        is_output=True,
        output_category="LCOW breakdown",
    )
    exports.add(
        obj=fs.costing.chemical_labor_maintenance_lcow,
        name="Chemical-labor-maintenance CAPEX LCOW",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of product water",
        rounding=3,
        description="CLM levelized cost of water",
        is_input=False,
        is_output=True,
        output_category="LCOW breakdown",
    )
    exports.add(
        obj=fs.costing.pumping_energy_aggregate_lcow,
        name="Pumping energy aggregate LCOW",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of product water",
        rounding=3,
        description="System pumping energy levelized cost of water",
        is_input=False,
        is_output=True,
        output_category="LCOW breakdown",
    )
    total_opex = fs.costing.total_operating_cost / fs.costing.annual_water_production
    exports.add(
        obj=total_opex,
        name="Total OPEX LCOW",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of product water",
        rounding=3,
        description="Total operating expenses levelized cost of water",
        is_input=False,
        is_output=True,
        output_category="LCOW breakdown",
    )

    # Stage metrics
    for idx, pump in fs.PrimaryPumps.items():
        exports.add(
            obj=pyunits.convert(
                pump.control_volume.properties_out[0].pressure, to_units=pyunits.bar
            ),
            name=f"Feed pressure - stage {idx}",
            ui_units=pyunits.bar,
            display_units="bar",
            rounding=1,
            description=f"Feed pressure of stage {idx}",
            is_input=False,
            is_output=True,
            output_category="Primary pump feed pressure",
        )

    for idx, stage in fs.ROUnits.items():
        exports.add(
            obj=stage.rejection_phase_comp[0, "Liq", "NaCl"] * 100,
            name=f"Observed NaCl rejection - stage {idx}",
            ui_units=pyunits.dimensionless,
            display_units="%",
            rounding=2,
            description=f"Observed salt rejection of stage {idx}",
            is_input=False,
            is_output=True,
            output_category="RO observed rejection",
        )

    for idx, stage in fs.ROUnits.items():
        A_comp = pyunits.convert(
            stage.A_comp[0, "H2O"],
            to_units=pyunits.L / pyunits.hr / pyunits.m**2 / pyunits.bar,
        )
        exports.add(
            obj=A_comp,
            name=f"Water permeability coefficient - stage {idx}",
            ui_units=pyunits.L / pyunits.hr / pyunits.m**2 / pyunits.bar,
            display_units="LMH/bar",
            rounding=2,
            description=f"A value of stage {idx}",
            is_input=False,
            is_output=True,
            output_category="RO water permeability coefficient (A)",
        )

    for idx, stage in fs.ROUnits.items():
        B_comp = pyunits.convert(
            stage.B_comp[0, "NaCl"], to_units=pyunits.L / pyunits.hr / pyunits.m**2
        )
        exports.add(
            obj=B_comp,
            name=f"Salt permeability coefficient - stage {idx}",
            ui_units=pyunits.L / pyunits.hr / pyunits.m**2,
            display_units="LMH",
            rounding=2,
            description=f"B value of stage {idx}",
            is_input=False,
            is_output=True,
            output_category="RO salt permeability coefficient (A)",
        )

    for idx, stage in fs.ROUnits.items():
        water_flux = pyunits.convert(
            stage.flux_mass_phase_comp_avg[0, "Liq", "H2O"],
            to_units=pyunits.kg / pyunits.hr / pyunits.m**2,
        )
        exports.add(
            obj=water_flux,
            name=f"Average water flux - stage {idx}",
            ui_units=pyunits.kg / pyunits.hr / pyunits.m**2,
            display_units="kg/hr-m2",
            rounding=2,
            description=f"Average water flux of stage {idx}",
            is_input=False,
            is_output=True,
            output_category="RO water flux",
        )

    for idx, stage in fs.ROUnits.items():
        salt_flux = pyunits.convert(
            stage.flux_mass_phase_comp_avg[0, "Liq", "NaCl"],
            to_units=pyunits.g / pyunits.m**2 / pyunits.h,
        )
        exports.add(
            obj=salt_flux,
            name=f"Average NaCl flux - stage {idx}",
            ui_units=pyunits.g / pyunits.m**2 / pyunits.h,
            display_units="gMH",
            rounding=1,
            description=f"Average salt flux of stage {idx}",
            is_input=False,
            is_output=True,
            output_category="RO NaCl flux",
        )

    for idx, stage in fs.ROUnits.items():
        exports.add(
            obj=stage.feed_side.N_Re[0, 0],
            name=f"Inlet Reynolds number - stage {idx}",
            ui_units=pyunits.dimensionless,
            display_units=" ",
            rounding=2,
            description=f"Reynolds number of stage {idx}",
            is_input=False,
            is_output=True,
            output_category="RO Reynolds number",
        )

    for idx, stage in fs.ROUnits.items():
        exports.add(
            obj=stage.feed_side.velocity[0, 0],
            name=f"Inlet crossflow velocity - stage {idx}",
            ui_units=pyunits.m / pyunits.s,
            display_units="m/s",
            rounding=2,
            description=f"Inlet velocity of stage {idx}",
            is_input=False,
            is_output=True,
            output_category="RO Crossflow velocity",
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
        if build_options["ACase"].value == ACase.fixed:
            optimize_set_up(
                m,
                set_default_bounds_on_module_dimensions=True,
                Cbrine=None,
                A_case=build_options["ACase"].value,
                B_case=build_options["BCase"].value,
                AB_tradeoff=build_options["ABTradeoff"].value,
                A_value=4.2e-12,
                permeate_quality_limit=500e-6,
                AB_gamma_factor=1,
                B_max=build_options["BMax"].value,
            )
        else:
            optimize_set_up(
                m,
                set_default_bounds_on_module_dimensions=True,
                Cbrine=None,
                A_case=build_options["ACase"].value,
                B_case=build_options["BCase"].value,
                AB_tradeoff=build_options["ABTradeoff"].value,
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
            A_value=4.2e-12,
            B_case=BCase.optimize,
            AB_tradeoff=ABTradeoff.none,
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
