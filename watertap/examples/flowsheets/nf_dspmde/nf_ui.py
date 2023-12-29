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
from watertap.ui.fsapi import FlowsheetInterface, FlowsheetCategory
from watertap.examples.flowsheets.nf_dspmde import nf
from watertap.examples.flowsheets.nf_dspmde import nf_with_bypass
from watertap.unit_models.nanofiltration_DSPMDE_0D import ConcentrationPolarizationType
from pyomo.environ import units as pyunits
from idaes.core.solvers import get_solver


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
                "values_allowed": ["false", "true"],
                "value": "false",
            },
            "ConcentrationPolarization": {
                "name": "ConcentrationPolarization",
                "display_name": "Concentration Polarization Type",
                "values_allowed": ["calculated", "none"],
                "value": "calculated",
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
        ui_units=pyunits.L / pyunits.hr,
        display_units="L/h",
        rounding=2,
        description="Inlet volumetric flow rate",
        is_input=True,
        input_category="Feed",
        is_output=False,
        output_category="Feed",
    )
    for (phase, ion), obj in fs.feed.properties[0].conc_mass_phase_comp.items():
        if ion != "H2O":
            exports.add(
                obj=obj,
                name="{} concentration".format(ion),
                ui_units=pyunits.mg / pyunits.L,
                display_units="mg/L",
                rounding=2,
                description="{} concentration".format(ion),
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
        obj=fs.NF.nfUnit.area,
        name="NF area",
        ui_units=pyunits.m**2,
        display_units="m^2",
        rounding=2,
        description="NF pump pressure",
        is_input=True,
        input_category="NF design",
        is_output=True,
        output_category="NF design",
    )
    exports.add(
        obj=fs.NF.nfUnit.recovery_vol_phase[0.0, "Liq"],
        name="NF water recovery".format(ion),
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=3,
        description="NF design",
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
    exports.add(
        obj=fs.NF.nf_flux,
        name="NF water flux",
        ui_units=pyunits.dimensionless,
        display_units="LMH",
        rounding=2,
        description="NF design",
        is_input=False,
        input_category="NF design",
        is_output=True,
        output_category="NF design",
    )
    exports.add(
        obj=fs.NF.nfUnit.radius_pore,
        name="Pore size",
        ui_units=pyunits.nm,
        display_units="nm",
        rounding=2,
        description="NF membrane props.",
        is_input=True,
        input_category="NF membrane props.",
        is_output=True,
        output_category="NF membrane props.",
    )
    exports.add(
        obj=fs.NF.nfUnit.membrane_thickness_effective,
        name="Effective membrane thickness",
        ui_units=pyunits.nm,
        display_units="nm",
        rounding=2,
        description="NF membrane props.",
        is_input=True,
        input_category="NF membrane props.",
        is_output=True,
        output_category="NF membrane props.",
    )
    exports.add(
        obj=fs.NF.nfUnit.membrane_charge_density[0],
        name="Charge",
        ui_units=pyunits.mol / pyunits.m**3,
        display_units="mol/m^3",
        rounding=2,
        description="NF membrane props.",
        is_input=True,
        input_category="NF membrane props.",
        is_output=True,
        output_category="NF membrane props.",
    )
    exports.add(
        obj=fs.NF.nfUnit.dielectric_constant_pore[0],
        name="Dielectric constant for pore",
        ui_units=pyunits.dimensionless,
        display_units="",
        rounding=2,
        description="NF membrane props.",
        is_input=True,
        input_category="NF membrane props.",
        is_output=True,
        output_category="NF membrane props.",
    )

    exports.add(
        obj=fs.costing.nanofiltration.membrane_cost,
        name="Membrane cost",
        ui_units=fs.costing.base_currency / pyunits.m**2,
        display_units="$/m^2",
        rounding=2,
        description="NF CAPEX",
        is_input=True,
        input_category="NF CAPEX",
        is_output=True,
        output_category="NF CAPEX",
    )
    exports.add(
        obj=fs.costing.nanofiltration.factor_membrane_replacement,
        name="Membrane replacment rate",
        ui_units=pyunits.year**-1,
        display_units="fraction/year",
        rounding=3,
        description="NF CAPEX",
        is_input=True,
        input_category="NF CAPEX",
        is_output=True,
        output_category="NF CAPEX",
    )
    exports.add(
        obj=fs.costing.electricity_cost,
        name="Electricity cost",
        ui_units=fs.costing.base_currency / pyunits.kWh,
        display_units="$/kWh",
        rounding=2,
        description="NF OPEX",
        is_input=True,
        input_category="NF OPEX",
        is_output=True,
        output_category="NF OPEX",
        is_readonly=False,
    )
    exports.add(
        obj=fs.costing.utilization_factor,
        name="Plant capacity utilization",
        ui_units=pyunits.dimensionless,
        display_units="fraction of uptime",
        rounding=2,
        description="NF OPEX",
        is_input=True,
        input_category="NF OPEX",
        is_output=True,
        output_category="NF OPEX",
    )
    exports.add(
        obj=fs.costing.factor_maintenance_labor_chemical,
        name="Maintenance-labor-chemical factor",
        ui_units=pyunits.year**-1,
        display_units="fraction of equipment cost/year",
        rounding=4,
        description="NF OPEX",
        is_input=True,
        input_category="NF OPEX",
        is_output=True,
        output_category="NF OPEX",
    )
    exports.add(
        obj=fs.costing.disposal_cost,
        name="Disposal cost",
        ui_units=pyunits.USD_2020 / pyunits.m**3,
        display_units="$/m^3",
        rounding=4,
        description="System constraints",
        is_input=True,
        input_category="System constraints",
        is_output=True,
        output_category="System constraints",
    )

    exports.add(
        obj=fs.product.max_hardness,
        name="Product hardness as CaCO3",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=2,
        description="System constraints",
        is_input=True,
        input_category="System constraints",
        is_output=False,
        output_category="System constraints",
    )

    exports.add(
        obj=fs.product.properties[0].total_hardness,
        name="Product hardness",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=2,
        description="System streams",
        is_input=False,
        input_category="System streams",
        is_output=True,
        output_category="System streams",
    )
    exports.add(
        obj=fs.product.properties[0].flow_vol_phase["Liq"],
        name="Product volume flow",
        ui_units=pyunits.L / pyunits.hr,
        display_units="L/h",
        rounding=2,
        description="System streams",
        is_input=False,
        input_category="System streams",
        is_output=True,
        output_category="System streams",
    )
    exports.add(
        obj=fs.feed.properties[0].total_hardness,
        name="Feed hardness",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=2,
        description="System streams",
        is_input=False,
        input_category="System streams",
        is_output=True,
        output_category="System streams",
    )

    exports.add(
        obj=fs.disposal.properties[0].total_hardness,
        name="Disposal hardness",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=2,
        description="System streams",
        is_input=False,
        input_category="System streams",
        is_output=True,
        output_category="System streams",
    )
    exports.add(
        obj=fs.disposal.properties[0].flow_vol_phase["Liq"],
        name="Disposal volume flow",
        ui_units=pyunits.L / pyunits.hr,
        display_units="L/h",
        rounding=2,
        description="System streams",
        is_input=False,
        input_category="System streams",
        is_output=True,
        output_category="System streams",
    )

    exports.add(
        obj=fs.costing.LCOW,
        name="System cost",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m^3",
        rounding=2,
        description="Process cost and operating metrics",
        is_input=False,
        input_category="Process cost and operating metrics",
        is_output=True,
        output_category="Process cost and operating metrics",
    )
    exports.add(
        obj=fs.costing.specific_energy_consumption,
        name="System-level specific energy consumption",
        ui_units=pyunits.hr * pyunits.kW / pyunits.m**3,
        display_units="kWh/m^3",
        rounding=4,
        description="Process cost and operating metrics",
        is_input=False,
        input_category="Process cost and operating metrics",
        is_output=True,
        output_category="Process cost and operating metrics",
    )
    try:
        if build_options["Bypass"].value == "true":
            exports.add(
                obj=fs.by_pass_splitter.split_fraction[0, "bypass"],
                name="NF bypass",
                ui_units=pyunits.dimensionless,
                display_units="fraction",
                rounding=4,
                description="Bypass design",
                is_input=True,
                input_category="Bypass design",
                is_output=True,
                output_category="Bypass design",
            )
    except Exception as e:
        print(f"error adding bypass: {e}")

    for (t, phase, ion), obj in fs.NF.nfUnit.rejection_intrinsic_phase_comp.items():
        exports.add(
            obj=obj,
            name="{} intrinsic rejection".format(ion),
            ui_units=pyunits.dimensionless,
            display_units="fraction",
            rounding=5,
            description="NF int. rejection",
            is_input=False,
            input_category="NF intrinsic rejection",
            is_output=True,
            output_category="NF intrinsic rejection",
        )
    for (t, phase, ion), obj in fs.NF.nfUnit.rejection_observed_phase_comp.items():
        exports.add(
            obj=obj,
            name="{} obs. rejection".format(ion),
            ui_units=pyunits.dimensionless,
            display_units="fraction",
            rounding=5,
            description="NF observed rejection",
            is_input=False,
            input_category="NF observed rejection",
            is_output=True,
            output_category="NF observed rejection",
        )


def build_flowsheet(build_options=None, **kwargs):
    # build and solve initial flowsheet
    if build_options is not None:
        if build_options["Bypass"].value == "true":  # build with bypass
            solver = get_solver()
            m = nf_with_bypass.build()
            concentrationType = build_options["ConcentrationPolarization"].value
            # print(f'setting concentration polarization type to {concentrationType}')
            m.fs.NF.nfUnit.config.concentration_polarization_type = (
                ConcentrationPolarizationType[concentrationType]
            )
            nf_with_bypass.initialize(m, solver)
            nf_with_bypass.unfix_opt_vars(m)
            nf.add_objective(m)
        else:  # build without bypass
            solver = get_solver()
            m = nf.build()
            concentrationType = build_options["ConcentrationPolarization"].value
            # print(f'setting concentration polarization type to {concentrationType}')
            m.fs.NF.nfUnit.config.concentration_polarization_type = (
                ConcentrationPolarizationType[concentrationType]
            )
            nf.initialize(m, solver)
            nf.add_objective(m)
            nf.unfix_opt_vars(m)
    else:  # build without bypass
        solver = get_solver()
        m = nf.build()
        nf.initialize(m, solver)
        nf.add_objective(m)
        nf.unfix_opt_vars(m)

    return m


def get_diagram(build_options):
    if build_options["Bypass"].value == "true":
        return "nf_with_bypass_ui.png"
    else:
        return "nf_ui.png"


def solve_flowsheet(flowsheet=None):
    fs = flowsheet
    solver = get_solver()
    results = nf.optimize(fs, solver)
    return results
