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
from watertap.core.solvers import get_solver
from idaes_flowsheet_processor.api import FlowsheetInterface
from watertap.flowsheets.seawater_RO_desalination.seawater_RO_desalination import (
    build,
    set_operating_conditions,
    initialize_system,
    solve,
    add_costing,
    initialize_costing,
)
from pyomo.environ import units as pyunits


def export_to_ui():
    return FlowsheetInterface(
        name="SWRO",
        do_export=export_variables,
        do_build=build_flowsheet,
        do_solve=solve_flowsheet,
        build_options={
            "ERD_type": {
                "name": "ERD_type",
                "display_name": "Energy Recovery Device Type",
                "values_allowed": ["pressure_exchanger", "pump_as_turbine"],
                "value": "pressure_exchanger",  # default value
            },
        },
    )


def export_variables(flowsheet=None, exports=None, build_options=None, **kwargs):
    fs = flowsheet
    # --- Input data ---
    # Feed conditions
    exports.add(
        obj=fs.feed.flow_vol[0],
        name="Feed volumetric flow rate",
        ui_units=pyunits.m**3 / pyunits.day,
        display_units="m3/day",
        rounding=2,
        description="Inlet volumetric flow rate",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.conc_mass_comp[0, "tds"],
        name="Feed TDS concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=2,
        description="Total dissolved solids concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.conc_mass_comp[0, "tss"],
        name="Feed TSS concentration",
        ui_units=pyunits.g / pyunits.m**3,
        display_units="g/m3",
        rounding=2,
        description="Total suspended solids concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    # --- Pretreatment ---
    # TODO: Consider whether ZO removal fractions should be exposed - currently they are not
    # Unit model data, ferric chloride addition
    exports.add(
        obj=fs.pretreatment.ferric_chloride_addition.chemical_dosage,
        name="Ferric chloride dosage",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=1,
        description="Chemical dosage of ferric chloride",
        is_input=True,
        input_category="Pretreatment",
        is_output=False,
    )
    exports.add(
        obj=fs.pretreatment.ferric_chloride_addition.solution_density,
        name="Ferric chloride solution density",
        ui_units=pyunits.kg / pyunits.m**3,
        display_units="kg/m3",
        rounding=1,
        description="Density of ferric chloride",
        is_input=True,
        input_category="Pretreatment",
        is_output=False,
    )
    exports.add(
        obj=fs.pretreatment.ferric_chloride_addition.ratio_in_solution,
        name="Ferric chloride ratio in solution",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=1,
        description="Ratio of ferric chloride in solution",
        is_input=True,
        input_category="Pretreatment",
        is_output=False,
    )
    # Unit model data, chlorination
    exports.add(
        obj=fs.pretreatment.chlorination.contact_time,
        name="Chlorination contact time",
        ui_units=pyunits.hr,
        display_units="hr",
        rounding=1,
        description="Contact time of chlorination",
        is_input=True,
        input_category="Pretreatment",
        is_output=False,
    )
    exports.add(
        obj=fs.pretreatment.chlorination.concentration_time,
        name="Chlorination concentration time",
        ui_units=pyunits.mg * pyunits.min / pyunits.L,
        display_units="mg*min/L",
        rounding=1,
        description="Concentration time of chlorination",
        is_input=True,
        input_category="Pretreatment",
        is_output=False,
    )
    exports.add(
        obj=fs.pretreatment.chlorination.chlorine_decay_rate,
        name="Chlorine decay rate",
        ui_units=pyunits.hr,
        display_units="hr",
        rounding=1,
        description="Decay rate of chlorine",
        is_input=True,
        input_category="Pretreatment",
        is_output=False,
    )
    # Unit model data, storage tank
    exports.add(
        obj=fs.pretreatment.storage_tank_1.storage_time,
        name="Tank 1 storage time",
        ui_units=pyunits.hr,
        display_units="hr",
        rounding=1,
        description="Storage time of tank 1",
        is_input=True,
        input_category="Pretreatment",
        is_output=False,
    )
    exports.add(
        obj=fs.pretreatment.storage_tank_1.surge_capacity,
        name="Tank 1 surge capacity",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=1,
        description="Surge capacity of tank 1",
        is_input=True,
        input_category="Pretreatment",
        is_output=False,
    )
    # Unit model data, antiscalant addition
    exports.add(
        obj=fs.pretreatment.anti_scalant_addition.chemical_dosage,
        name="Antiscalant dosage",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=1,
        description="Chemical dosage of antiscalant",
        is_input=True,
        input_category="Pretreatment",
        is_output=False,
    )
    exports.add(
        obj=fs.pretreatment.anti_scalant_addition.solution_density,
        name="Antiscalant density",
        ui_units=pyunits.kg / pyunits.m**3,
        display_units="kg/m3",
        rounding=1,
        description="Density of antiscalant",
        is_input=True,
        input_category="Pretreatment",
        is_output=False,
    )
    exports.add(
        obj=fs.pretreatment.anti_scalant_addition.ratio_in_solution,
        name="Antiscalant ratio in solution",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=1,
        description="Ratio of Antiscalant in solution",
        is_input=True,
        input_category="Pretreatment",
        is_output=False,
    )
    # --- Desalination ---
    # Unit model data, pump
    exports.add(
        obj=fs.desalination.P1.efficiency_pump,
        name="Pump 1 efficiency",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=1,
        description="Efficiency of pump 1 in the desalination circuit",
        is_input=True,
        input_category="Desalination",
        is_output=False,
    )
    exports.add(
        obj=fs.desalination.P1.control_volume.properties_out[0].pressure,
        name="Pump 1 pressure",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=1,
        description="Pressure of pump 1 in the desalination circuit",
        is_input=True,
        input_category="Desalination",
        is_output=False,
    )
    exports.add(
        obj=fs.desalination.P2.efficiency_pump,
        name="Pump 2 efficiency",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=1,
        description="Efficiency of pump 2 in the desalination circuit",
        is_input=True,
        input_category="Desalination",
        is_output=False,
    )
    # Unit model data, RO
    exports.add(
        obj=fs.desalination.RO.A_comp[0, "H2O"],
        name="RO water permeability coefficient",
        ui_units=pyunits.L / pyunits.hr / pyunits.m**2 / pyunits.bar,
        display_units="LMH/bar",
        rounding=2,
        description="Water permeability coefficient",
        is_input=True,
        input_category="Desalination",
        is_output=False,
    )

    exports.add(
        obj=fs.desalination.RO.B_comp[0, "TDS"],
        name="RO salt permeability coefficient",
        ui_units=pyunits.L / pyunits.hr / pyunits.m**2,
        display_units="LMH",
        rounding=2,
        description="Salt permeability coefficient",
        is_input=True,
        input_category="Desalination",
        is_output=False,
    )
    exports.add(
        obj=fs.desalination.RO.feed_side.channel_height,
        name="RO feed-side channel height",
        ui_units=pyunits.mm,
        display_units="mm",
        rounding=2,
        description="Feed-side channel height",
        is_input=True,
        input_category="Desalination",
        is_output=False,
    )
    exports.add(
        obj=fs.desalination.RO.feed_side.spacer_porosity,
        name="RO feed-side spacer porosity",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Feed-side spacer porosity",
        is_input=True,
        input_category="Desalination",
        is_output=False,
    )
    exports.add(
        obj=fs.desalination.RO.permeate.pressure[0],
        name="RO permeate-side pressure",
        ui_units=pyunits.bar,
        display_units="bar",
        rounding=2,
        description="Permeate-side pressure",
        is_input=True,
        input_category="Desalination",
        is_output=False,
    )
    exports.add(
        obj=fs.desalination.RO.width,
        name="RO stage width",
        ui_units=pyunits.m,
        display_units="m",
        rounding=2,
        description="Stage width",
        is_input=True,
        input_category="Desalination",
        is_output=False,
    )
    exports.add(
        obj=fs.desalination.RO.recovery_mass_phase_comp[0, "Liq", "H2O"],
        name="RO water mass recovery",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Water mass recovery of RO unit",
        is_input=True,
        input_category="Desalination",
        is_output=True,
        output_category="System metrics",
    )

    # Unit model data, pressure exchanger
    if build_options["ERD_type"].value == "pressure_exchanger":
        exports.add(
            obj=fs.desalination.PXR.efficiency_pressure_exchanger,
            name="Pressure exchanger efficiency",
            ui_units=pyunits.dimensionless,
            display_units="fraction",
            rounding=2,
            description="Efficiency of pressure exchanger",
            is_input=True,
            input_category="Desalination",
            is_output=False,
        )
        exports.add(
            obj=fs.desalination.P2.efficiency_pump,
            name="Pump 2 efficiency",
            ui_units=pyunits.dimensionless,
            display_units="fraction",
            rounding=1,
            description="Efficiency of pump 2 in the desalination circuit",
            is_input=True,
            input_category="Desalination",
            is_output=False,
        )

    # Unit model data, energy recovery device
    if build_options["ERD_type"].value == "pump_as_turbine":
        exports.add(
            obj=fs.desalination.ERD.efficiency_pump,
            name="Energy recovery device efficiency",
            ui_units=pyunits.dimensionless,
            display_units="fraction",
            rounding=2,
            description="Efficiency of energy recovery device",
            is_input=True,
            input_category="Desalination",
            is_output=False,
        )

    # --- Posttreatment ---
    # Unit model data, storage tanks
    exports.add(
        obj=fs.posttreatment.storage_tank_2.storage_time,
        name="Tank 2 storage time",
        ui_units=pyunits.hr,
        display_units="hr",
        rounding=1,
        description="Storage time of tank 2",
        is_input=True,
        input_category="Posttreatment",
        is_output=False,
    )
    exports.add(
        obj=fs.posttreatment.storage_tank_2.surge_capacity,
        name="Tank 2 surge capacity",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=1,
        description="Surge capacity of tank 2",
        is_input=True,
        input_category="Posttreatment",
        is_output=False,
    )
    exports.add(
        obj=fs.posttreatment.storage_tank_3.storage_time,
        name="Tank 3 storage time",
        ui_units=pyunits.hr,
        display_units="hr",
        rounding=1,
        description="Storage time of tank 3",
        is_input=True,
        input_category="Posttreatment",
        is_output=False,
    )
    exports.add(
        obj=fs.posttreatment.storage_tank_3.surge_capacity,
        name="Tank 3 surge capacity",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=1,
        description="Surge capacity of tank 3",
        is_input=True,
        input_category="Posttreatment",
        is_output=False,
    )
    # Unit model data, UV AOP
    exports.add(
        obj=fs.posttreatment.uv_aop.uv_reduced_equivalent_dose,
        name="UV reduced equivalent dose",
        ui_units=pyunits.mJ / pyunits.cm**2,
        display_units="mJ/cm^2",
        rounding=1,
        description="Ultraviolet reduced equivalent dose",
        is_input=True,
        input_category="Posttreatment",
        is_output=False,
    )
    exports.add(
        obj=fs.posttreatment.uv_aop.uv_transmittance_in,
        name="UV transmittance",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=1,
        description="Ultraviolet transmittance",
        is_input=True,
        input_category="Posttreatment",
        is_output=False,
    )
    exports.add(
        obj=fs.posttreatment.uv_aop.oxidant_dose,
        name="UV AOP oxidant dose",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=1,
        description="Oxidant dose for ultraviolet advanced oxidation process",
        is_input=True,
        input_category="Posttreatment",
        is_output=False,
    )
    # Unit model data, lime addition
    exports.add(
        obj=fs.posttreatment.lime_addition.chemical_dosage,
        name="Lime chemical dosage",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=1,
        description="Chemical dosage of lime",
        is_input=True,
        input_category="Posttreatment",
        is_output=False,
    )
    exports.add(
        obj=fs.posttreatment.lime_addition.solution_density,
        name="Lime solution density",
        ui_units=pyunits.kg / pyunits.m**3,
        display_units="mg/L",
        rounding=1,
        description="Solution density of lime",
        is_input=True,
        input_category="Posttreatment",
        is_output=False,
    )
    exports.add(
        obj=fs.posttreatment.lime_addition.ratio_in_solution,
        name="Lime ratio in solution",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=1,
        description="Ratio of lime in solution",
        is_input=True,
        input_category="Posttreatment",
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
    # Municipal
    exports.add(
        obj=fs.municipal.properties[0].flow_mass_comp["H2O"],
        name="Municipal water mass flow",
        ui_units=pyunits.kg / pyunits.s,
        display_units="kg/s",
        rounding=2,
        description="Municipal water mass flow",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.municipal.properties[0].flow_mass_comp["tds"],
        name="Municipal tds mass flow",
        ui_units=pyunits.kg / pyunits.s,
        display_units="kg/s",
        rounding=2,
        description="Municipal total dissolved solids mass flow",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    # Landfill
    exports.add(
        obj=fs.landfill.properties[0].flow_mass_comp["H2O"],
        name="Landfill water mass flow",
        ui_units=pyunits.kg / pyunits.s,
        display_units="kg/s",
        rounding=2,
        description="Landfill water mass flow",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.landfill.properties[0].flow_mass_comp["tds"],
        name="Landfill tds mass flow",
        ui_units=pyunits.kg / pyunits.s,
        display_units="kg/s",
        rounding=2,
        description="Landfill total dissolved solids mass flow",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.landfill.properties[0].flow_mass_comp["tss"],
        name="Landfill tss mass flow",
        ui_units=pyunits.kg / pyunits.s,
        display_units="kg/s",
        rounding=2,
        description="Landfill total suspended solids mass flow",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    # Disposal
    exports.add(
        obj=fs.disposal.properties[0].flow_mass_comp["H2O"],
        name="Disposal water mass flow",
        ui_units=pyunits.kg / pyunits.s,
        display_units="kg/s",
        rounding=2,
        description="Disposal water mass flow",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.disposal.properties[0].flow_mass_comp["TDS"],
        name="Disposal tds mass flow",
        ui_units=pyunits.kg / pyunits.s,
        display_units="kg/s",
        rounding=2,
        description="Disposal total dissolved solids mass flow",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )

    # System metrics
    exports.add(
        obj=fs.LCOW,
        name="Levelized cost of water",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of product water",
        rounding=3,
        description="Levelized cost of water (LCOW)",
        is_input=False,
        is_output=True,
        output_category="System metrics",
    )
    total_capital_cost = (
        pyunits.convert(fs.zo_costing.total_capital_cost, to_units=pyunits.USD_2018)
        + fs.ro_costing.total_capital_cost
    )
    exports.add(
        obj=total_capital_cost,
        name="Total capital cost",
        ui_units=pyunits.USD_2018,
        display_units="USD_2018",
        rounding=3,
        description="Total capital cost in USD_2018",
        is_input=False,
        is_output=True,
        output_category="System metrics",
    )
    total_operating_cost = (
        pyunits.convert(
            fs.zo_costing.total_fixed_operating_cost,
            to_units=pyunits.USD_2018 / pyunits.year,
        )
        + pyunits.convert(
            fs.zo_costing.total_variable_operating_cost,
            to_units=pyunits.USD_2018 / pyunits.year,
        )
        + fs.ro_costing.total_operating_cost
    )
    exports.add(
        obj=total_operating_cost,
        name="Total operating cost",
        ui_units=pyunits.USD_2018 / pyunits.yr,
        display_units="USD_2018/yr",
        rounding=3,
        description="Total capital cost in USD_2018 per year",
        is_input=False,
        is_output=True,
        output_category="System metrics",
    )


def build_flowsheet(build_options=None, **kwargs):
    if build_options is not None:
        if build_options["ERD_type"].value == "pump_as_turbine":
            m = build(erd_type="pump_as_turbine")
        else:
            m = build(erd_type="pressure_exchanger")

    set_operating_conditions(m)
    initialize_system(m)
    # Solve flowsheet after initializing system
    solve(m)

    add_costing(m)
    initialize_costing(m)
    # Solve flowsheet with costing
    solve(m)

    return m


def solve_flowsheet(flowsheet=None):
    fs = flowsheet
    results = solve(fs)
    return results
