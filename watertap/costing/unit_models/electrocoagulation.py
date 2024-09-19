import pyomo.environ as pyo
from ..util import (
    register_costing_parameter_block,
    make_capital_cost_var,
    make_fixed_operating_cost_var,
)

# Costing equations from:


def build_electrocoagulation_cost_param_block(blk):

    costing = blk.parent_block()

    blk.ec_reactor_cap_base = pyo.Var(
        initialize=1.15e4,
        units=pyo.units.USD_2020,
        doc="Reactor capital cost base cost",
    )

    blk.ec_reactor_cap_exp = pyo.Var(
        initialize=0.45,
        units=pyo.units.dimensionless,
        doc="Reactor capital cost exponent",
    )

    blk.ec_reactor_cap_material_coeff = pyo.Var(
        initialize=3.4,
        units=pyo.units.dimensionless,
        doc="Reactor capital cost material coeff (3.4 for stainless steel; 0.062 for PVC)",
    )

    blk.ec_reactor_cap_safety_factor = pyo.Var(
        initialize=2.5,
        units=pyo.units.dimensionless,
        doc="Reactor capital cost safety factor",
    )

    blk.ec_admin_lab_cap_base = pyo.Var(
        initialize=69195,
        units=pyo.units.USD_2010,
        doc="Lab + administration + building base cost",
    )

    blk.ec_admin_lab_cap_exp = pyo.Var(
        initialize=0.5523,
        units=pyo.units.dimensionless,
        doc="Lab + administration + building exponent",
    )

    blk.ec_power_supply_base = pyo.Var(
        initialize=260000,
        units=pyo.units.USD_2021,
        doc="DC power supply + transforer + electrical connection base cost",
    )

    blk.ec_admin_lab_op_base = pyo.Var(
        initialize=88589,
        units=pyo.units.USD_2010 / pyo.units.year,
        doc="Admin + lab + building operational base cost",
    )

    blk.ec_admin_lab_op_exp = pyo.Var(
        initialize=0.4589,
        units=pyo.units.dimensionless,
        doc="Admin + lab + building operational cost exponent",
    )

    blk.sludge_handling_cost = pyo.Var(
        initialize=0,  # 0.025 $/kg from NMSU
        units=costing.base_currency / pyo.units.kg,
        doc="Cost per kg for handling generated solids",
    )

    blk.ec_labor_maint_factor = pyo.Var(
        initialize=0.063,
        units=costing.base_currency / pyo.units.m**3,
        doc="Labor and maintenance cost factor",
    )

    blk.current_per_reactor = pyo.Var(
        initialize=3000,
        units=pyo.units.ampere,
        doc="Current required per reactor/power supply",
    )

    blk.number_redundant_reactors = pyo.Var(
        initialize=2,
        units=pyo.units.dimensionless,
        doc="Number of redundant EC reactors",
    )

    blk.electrode_material_cost = pyo.Var(
        initialize=2,
        units=costing.base_currency / pyo.units.kg,
        doc="Cost of electrode per kg",
    )

    blk.fix_all_vars()


def build_aluminum_plate_cost_param_block(blk):
    costing = blk.parent_block()
    blk.cost = pyo.Param(
        initialize=2.23,
        units=costing.base_currency / pyo.units.kg,
        mutable=True,
        doc="Cost of aluminum plate per kg",
    )
    costing.register_flow_type("aluminum", blk.cost)


def build_iron_plate_cost_param_block(blk):
    costing = blk.parent_block()
    blk.cost = pyo.Param(
        initialize=3.41,
        units=costing.base_currency / pyo.units.kg,
        mutable=True,
        doc="Cost of iron plate per kg",
    )
    costing.register_flow_type("iron", blk.cost)


@register_costing_parameter_block(
    build_rule=build_aluminum_plate_cost_param_block,
    parameter_block_name="aluminum",
)
@register_costing_parameter_block(
    build_rule=build_iron_plate_cost_param_block,
    parameter_block_name="iron",
)
@register_costing_parameter_block(
    build_rule=build_electrocoagulation_cost_param_block,
    parameter_block_name="electrocoagulation",
)
def cost_electrocoagulation(blk):

    ec_params = blk.costing_package.electrocoagulation
    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, "TPEC")
    make_fixed_operating_cost_var(blk)

    ec = blk.unit_model

    flow_mgd = pyo.units.convert(
        ec.properties_in[0].flow_vol, to_units=pyo.units.Mgallons / pyo.units.day
    )
    flow_m3_yr = pyo.units.convert(
        ec.properties_in[0].flow_vol, to_units=pyo.units.m**3 / pyo.units.year
    )

    blk.annual_sludge_flow = pyo.units.convert(
        sum(
            ec.properties_waste[0].flow_mass_phase_comp["Liq", j] if j != "H2O" else 0
            for j in ec.properties_waste[0].params.component_list
        ),
        to_units=pyo.units.kg / pyo.units.year,
    )

    base_currency = blk.config.flowsheet_costing_block.base_currency

    blk.number_chambers_system = pyo.Param(
        initialize=3,
        mutable=True,
        units=pyo.units.dimensionless,
        doc="Number total chambers for system - EC chamber > flotation chamber > sedimentation chamber. All made of same material.",
    )

    blk.number_EC_reactors = pyo.Var(
        initialize=3,
        units=pyo.units.dimensionless,
        doc="Number EC cells and power supplies",
    )

    blk.capital_cost_reactor = pyo.Var(
        initialize=1e4,
        units=pyo.units.dimensionless,
        bounds=(0, None),
        doc="Cost of EC reactor",
    )

    blk.capital_cost_electrodes = pyo.Var(
        initialize=1e4,
        units=pyo.units.dimensionless,
        bounds=(0, None),
        doc="Cost of EC electrodes",
    )

    blk.capital_cost_power_supply = pyo.Var(
        initialize=1e6,
        units=pyo.units.dimensionless,
        bounds=(0, None),
        doc="Cost of EC power supply",
    )

    blk.capital_cost_admin_lab = pyo.Var(
        initialize=1e4,
        units=pyo.units.dimensionless,
        bounds=(0, None),
        doc="Cost of administration + lab + building, etc.",
    )

    blk.annual_labor_maintenance = pyo.Var(
        initialize=1e4,
        units=base_currency / pyo.units.year,
        bounds=(0, None),
        doc="Annual labor + maintenance cost",
    )

    blk.annual_sludge_management = pyo.Var(
        initialize=1e4,
        units=base_currency / pyo.units.year,
        bounds=(0, None),
        doc="Annual sludge management cost",
    )

    blk.annual_admin_lab = pyo.Var(
        initialize=1e4,
        units=base_currency / pyo.units.year,
        bounds=(0, None),
        doc="Annual administration + lab cost",
    )

    if ec.config.reactor_material == "pvc":
        ec_params.ec_reactor_cap_material_coeff.fix(0.062)

    elif ec.config.reactor_material == "stainless_steel":
        ec_params.ec_reactor_cap_material_coeff.fix(3.4)

    if ec.config.electrode_material == "aluminum":
        ec_params.electrode_material_cost.fix(2.23)

    elif ec.config.electrode_material == "iron":
        ec_params.electrode_material_cost.fix(3.41)

    blk.number_EC_reactors_constr = pyo.Constraint(
        expr=blk.number_EC_reactors
        == ec.applied_current / ec_params.current_per_reactor
        + ec_params.number_redundant_reactors
    )

    blk.capital_cost_reactor_constraint = pyo.Constraint(
        expr=blk.capital_cost_reactor
        == (
            (
                ec_params.ec_reactor_cap_base
                * (
                    ec.reactor_volume
                    * blk.number_EC_reactors
                    * blk.number_chambers_system
                )
                ** ec_params.ec_reactor_cap_exp
            )
            * ec_params.ec_reactor_cap_material_coeff
        )
        * ec_params.ec_reactor_cap_safety_factor
    )

    blk.capital_cost_electrodes_constraint = pyo.Constraint(
        expr=blk.capital_cost_electrodes
        == (ec.electrode_mass * ec.number_electrode_pairs * 2 * blk.number_EC_reactors)
        * ec_params.electrode_material_cost
    )

    blk.capital_cost_power_supply_constraint = pyo.Constraint(
        expr=blk.capital_cost_power_supply
        == ec_params.ec_power_supply_base * blk.number_EC_reactors
    )

    blk.capital_cost_other_constraint = pyo.Constraint(
        expr=blk.capital_cost_admin_lab
        == ec_params.ec_admin_lab_cap_base * flow_mgd**ec_params.ec_admin_lab_cap_exp
    )

    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == blk.capital_cost_reactor
        + blk.capital_cost_electrodes
        + blk.capital_cost_power_supply
        + blk.capital_cost_admin_lab
    )

    blk.annual_labor_maintenance_constraint = pyo.Constraint(
        expr=blk.annual_labor_maintenance
        == flow_m3_yr * ec_params.ec_labor_maint_factor
    )

    blk.annual_sludge_management_constraint = pyo.Constraint(
        expr=blk.annual_sludge_management
        == blk.annual_sludge_flow * ec_params.sludge_handling_cost
    )

    blk.annual_admin_lab_constraint = pyo.Constraint(
        expr=blk.annual_admin_lab
        == ec_params.ec_admin_lab_op_base * flow_mgd**ec_params.ec_admin_lab_op_exp
    )

    blk.fixed_operating_cost_constraint = pyo.Constraint(
        expr=blk.fixed_operating_cost
        == blk.annual_labor_maintenance
        + blk.annual_sludge_management
        + blk.annual_admin_lab
    )

    blk.annual_electrode_replacement_mass_flow = pyo.Expression(
        expr=pyo.units.convert(
            ec.metal_loading * flow_m3_yr, to_units=pyo.units.kg / pyo.units.year
        )
    )

    blk.costing_package.cost_flow(
        blk.annual_electrode_replacement_mass_flow, ec.config.electrode_material
    )

    blk.costing_package.cost_flow(ec.power_required, "electricity")
