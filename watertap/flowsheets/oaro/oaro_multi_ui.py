#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
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
from watertap.flowsheets.oaro.oaro_multi import (
    build,
    set_operating_conditions,
    initialize_system,
    optimize_set_up,
    solve,
    ERDtype,
)
from pyomo.environ import units as pyunits


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


def export_variables(flowsheet=None, exports=None, build_options=None, **kwargs):
    fs = flowsheet
    # --- Input data ---
    # Feed conditions
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

    # Unit model data, OARO
    for idx, stage in fs.OAROUnits.items():
        exports.add(
            obj=stage.A_comp[0, "H2O"],
            name=f"OARO {idx} water permeability coefficient",
            ui_units=pyunits.L / pyunits.hr / pyunits.m**2 / pyunits.bar,
            display_units="LMH/bar",
            rounding=2,
            description=f"Water permeability coefficient of OARO {idx}",
            is_input=True,
            input_category="OARO",
            is_output=False,
        )
        exports.add(
            obj=stage.B_comp[0, "NaCl"],
            name=f"OARO {idx} salt permeability coefficient",
            ui_units=pyunits.L / pyunits.hr / pyunits.m**2,
            display_units="LMH",
            rounding=2,
            description=f"Salt permeability coefficient of OARO {idx}",
            is_input=True,
            input_category="OARO",
            is_output=False,
        )
        exports.add(
            obj=stage.feed_side.channel_height,
            name=f"OARO {idx} feed-side channel height",
            ui_units=pyunits.mm,
            display_units="mm",
            rounding=2,
            description=f"Feed-side channel height of OARO {idx}",
            is_input=True,
            input_category="OARO",
            is_output=False,
        )
        exports.add(
            obj=stage.feed_side.spacer_porosity,
            name=f"OARO {idx} feed-side spacer porosity",
            ui_units=pyunits.dimensionless,
            display_units="fraction",
            rounding=2,
            description=f"Feed-side spacer porosity of OARO {idx}",
            is_input=True,
            input_category="OARO",
            is_output=False,
        )
        exports.add(
            obj=stage.permeate_side.channel_height,
            name=f"OARO {idx} permeate-side channel height",
            ui_units=pyunits.mm,
            display_units="mm",
            rounding=2,
            description=f"Permeate-side channel height of OARO {idx}",
            is_input=True,
            input_category="OARO",
            is_output=False,
        )
        exports.add(
            obj=stage.permeate_side.spacer_porosity,
            name=f"OARO {idx} permeate-side spacer porosity",
            ui_units=pyunits.dimensionless,
            display_units="fraction",
            rounding=2,
            description=f"Permeate-side spacer porosity of OARO {idx}",
            is_input=True,
            input_category="OARO",
            is_output=False,
        )
        exports.add(
            obj=stage.feed_inlet.pressure[0],
            name=f"OARO stage{idx} feed inlet operating pressure",
            ui_units=pyunits.bar,
            display_units="bar",
            rounding=1,
            description=f"OARO stage{idx} feed inlet operating pressure",
            is_input=True,
            input_category="OARO",
            is_output=True,
            output_category="OARO metrics",
        )
        exports.add(
            obj=stage.area,
            name=f"OARO {idx} membrane area",
            ui_units=pyunits.m**2,
            display_units="m2",
            rounding=2,
            description=f"Membrane area of OARO {idx}",
            is_input=True,
            input_category="OARO",
            is_output=True,
            output_category="OARO metrics",
        )

    # Unit model, RO
    exports.add(
        obj=fs.RO.A_comp[0, "H2O"],
        name="Water permeability coefficient",
        ui_units=pyunits.L / pyunits.hr / pyunits.m**2 / pyunits.bar,
        display_units="LMH/bar",
        rounding=2,
        description="Water permeability coefficient",
        is_input=True,
        input_category="Reverse Osmosis",
        is_output=False,
    )

    exports.add(
        obj=fs.RO.B_comp[0, "NaCl"],
        name="Salt permeability coefficient",
        ui_units=pyunits.L / pyunits.hr / pyunits.m**2,
        display_units="LMH",
        rounding=2,
        description="Salt permeability coefficient",
        is_input=True,
        input_category="Reverse Osmosis",
        is_output=False,
    )
    exports.add(
        obj=fs.RO.feed_side.channel_height,
        name="Feed-side channel height",
        ui_units=pyunits.mm,
        display_units="mm",
        rounding=2,
        description="Feed-side channel height",
        is_input=True,
        input_category="Reverse Osmosis",
        is_output=False,
    )
    exports.add(
        obj=fs.RO.feed_side.spacer_porosity,
        name="Feed-side spacer porosity",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Feed-side spacer porosity",
        is_input=True,
        input_category="Reverse Osmosis",
        is_output=False,
    )
    exports.add(
        obj=fs.RO.permeate.pressure[0],
        name="Permeate-side pressure",
        ui_units=pyunits.bar,
        display_units="bar",
        rounding=2,
        description="Permeate-side pressure",
        is_input=True,
        input_category="Reverse Osmosis",
        is_output=False,
    )
    exports.add(
        obj=fs.RO.inlet.pressure[0],
        name="RO feed inlet operating pressure",
        ui_units=pyunits.bar,
        display_units="bar",
        rounding=1,
        description="RO feed inlet operating pressure",
        is_input=True,
        input_category="Reverse Osmosis",
        is_output=True,
        output_category="RO metrics",
    )
    exports.add(
        obj=fs.RO.area,
        name=f"RO membrane area",
        ui_units=pyunits.m**2,
        display_units="m2",
        rounding=2,
        description=f"Membrane area of RO",
        is_input=True,
        input_category="Reverse Osmosis",
        is_output=True,
        output_category="RO metrics",
    )

    # Unit model data, primary pumps
    for idx, pump in fs.PrimaryPumps.items():
        exports.add(
            obj=pump.efficiency_pump[0],
            name=f"Primary pump {idx} efficiency",
            ui_units=pyunits.dimensionless,
            display_units="fraction",
            rounding=2,
            description=f"Efficiency of primary pump {idx}",
            is_input=True,
            input_category="Primary Pumps",
            is_output=False,
        )

    # Unit model data, recycle pumps
    for idx, pump in fs.RecyclePumps.items():
        exports.add(
            obj=pump.efficiency_pump[0],
            name=f"Recycle pump {idx} efficiency",
            ui_units=pyunits.dimensionless,
            display_units="fraction",
            rounding=2,
            description=f"Efficiency of recycle pump {idx}",
            is_input=True,
            input_category="Recycle Pumps",
            is_output=False,
        )

    # Unit model data, ERD
    for idx, erd in fs.EnergyRecoveryDevices.items():
        exports.add(
            obj=erd.efficiency_pump[0],
            name=f"ERD {idx} pump efficiency",
            ui_units=pyunits.dimensionless,
            display_units="fraction",
            rounding=2,
            description=f"Efficiency of energy recovery device {idx}",
            is_input=True,
            input_category="Energy Recovery Device",
            is_output=False,
        )
        exports.add(
            obj=erd.control_volume.properties_out[0].pressure,
            name=f"ERD {idx} operating pressure",
            ui_units=pyunits.bar,
            display_units="bar",
            rounding=2,
            description=f"Operating pressure of energy recovery device {idx}",
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
        obj=fs.mass_water_recovery,
        name="Water mass recovery",
        ui_units=pyunits.dimensionless,
        display_units=" ",
        rounding=2,
        description="System water mass recovery",
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
    total_area = sum(fs.OAROUnits[i].area for i in fs.NonFinalStages) + fs.RO.area
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

    # Design variables
    for idx, stage in fs.OAROUnits.items():
        exports.add(
            obj=stage.recovery_mass_phase_comp[0, "Liq", "H2O"],
            name=f"OARO stage{idx} water recovery",
            ui_units=pyunits.dimensionless,
            display_units=" ",
            rounding=2,
            description=f"OARO stage{idx} water mass recovery",
            is_input=False,
            is_output=True,
            output_category="OARO metrics",
        )
        exports.add(
            obj=stage.permeate_inlet.pressure[0],
            name=f"OARO stage{idx} permeate inlet operating pressure",
            ui_units=pyunits.bar,
            display_units="bar",
            rounding=1,
            description=f"OARO stage{idx} permeate inlet operating pressure",
            is_input=False,
            is_output=True,
            output_category="OARO metrics",
        )
        exports.add(
            obj=stage.flux_mass_phase_comp_avg[0, "Liq", "H2O"] / stage.dens_solvent,
            name=f"OARO stage{idx} average water flux",
            ui_units=pyunits.L / pyunits.m**2 / pyunits.h,
            display_units="L/m2/h",
            rounding=2,
            description=f"OARO stage{idx} average water flux",
            is_input=False,
            is_output=True,
            output_category="OARO metrics",
        )
        exports.add(
            obj=stage.flux_mass_phase_comp_avg[0, "Liq", "NaCl"],
            name=f"OARO stage{idx} average salt flux",
            ui_units=pyunits.g / pyunits.m**2 / pyunits.h,
            display_units="g/m2/h",
            rounding=2,
            description=f"OARO stage{idx} average salt flux",
            is_input=False,
            is_output=True,
            output_category="OARO metrics",
        )

    exports.add(
        obj=fs.RO.recovery_mass_phase_comp[0, "Liq", "H2O"],
        name="RO water recovery",
        ui_units=pyunits.dimensionless,
        display_units=" ",
        rounding=2,
        description="RO water mass recovery",
        is_input=False,
        is_output=True,
        output_category="RO metrics",
    )
    exports.add(
        obj=fs.RO.flux_mass_phase_comp_avg[0, "Liq", "H2O"] / fs.RO.dens_solvent,
        name="RO average water flux",
        ui_units=pyunits.L / pyunits.m**2 / pyunits.h,
        display_units="L/m2/h",
        rounding=2,
        description="RO average water flux",
        is_input=False,
        is_output=True,
        output_category="RO metrics",
    )
    exports.add(
        obj=fs.RO.flux_mass_phase_comp_avg[0, "Liq", "NaCl"],
        name="RO average salt flux",
        ui_units=pyunits.g / pyunits.m**2 / pyunits.h,
        display_units="g/m2/h",
        rounding=2,
        description="RO average salt flux",
        is_input=False,
        is_output=True,
        output_category="RO metrics",
    )


def build_flowsheet(erd_type=ERDtype.pump_as_turbine, build_options=None, **kwargs):
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


def solve_flowsheet(flowsheet=None):
    fs = flowsheet
    results = solve(fs)
    return results
