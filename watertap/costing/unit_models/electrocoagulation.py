#################################################################################
# WaterTAP Copyright (c) 2020-2025, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import pyomo.environ as pyo
from ..util import (
    register_costing_parameter_block,
    make_capital_cost_var,
    make_fixed_operating_cost_var,
)

# Costing data & equations adapted from:

# W. McGivney and S. Kawamura (2008)
# Cost Estimating Manual for Water Treatment Facilities
# DOI: 10.1002/9780470260036

# Power supply cost estimation from magna-power.com
# Linear equation fit to SL and MT series cost data

# R. Smith (2005)
# Chemical Process: Design and Integration, 1st ed.
# ISBN: 978-0471486817

# For material cost ratios (carbon steel = 1):
# https://www.engineeringtoolbox.com/piping-materials-cost-ratios-d_864.html
# https://web.mit.edu/course/3/3.11/www/modules/props.pdf

# A. R. Anuf, K. Ramaraj, V. S. Sivasankarapillai, R. Dhanusuraman, J. P. Maran, G. Rajeshkumar (2022)
# Optimization of electrocoagulation process for treatment of
#   rice mill effluent using response surface methodology
# DOI: 10.1016/j.jwpe.2022.103074


def build_aluminum_cost_param_block(blk):
    costing = blk.parent_block()
    blk.cost = pyo.Param(
        initialize=2.23,
        units=pyo.units.USD_2021 / pyo.units.kg,  # Anuf et al. (2022)
        mutable=True,
        doc="Cost of aluminum plate per kg",
    )
    costing.register_flow_type("aluminum", blk.cost)


def build_iron_cost_param_block(blk):
    costing = blk.parent_block()
    blk.cost = pyo.Param(
        initialize=3.41,
        units=pyo.units.USD_2021 / pyo.units.kg,  # Anuf et al. (2022)
        mutable=True,
        doc="Cost of iron plate per kg",
    )
    costing.register_flow_type("iron", blk.cost)


def build_electrocoagulation_cost_param_block(blk):

    costing = blk.parent_block()

    # Parameters from Table 2.1 for "Agitated Reactor" in Smith (2005)
    blk.reactor_capital_cost_base = pyo.Var(
        initialize=1.15e4,
        units=pyo.units.USD_2000,
        doc="Reactor capital cost base parameter",
    )

    blk.reactor_capital_cost_exponent = pyo.Var(
        initialize=0.45,
        units=pyo.units.dimensionless,
        doc="Reactor capital cost exponent",
    )

    # Material cost ratios can be found at:
    # https://www.engineeringtoolbox.com/piping-materials-cost-ratios-d_864.html
    # https://web.mit.edu/course/3/3.11/www/modules/props.pdf
    # Table 2.2 and Table 2.3 in Smith (2005)
    blk.reactor_material_coeff = pyo.Var(
        initialize=1,  # carbon steel = 1 is default
        units=pyo.units.dimensionless,
        doc="Reactor capital cost material coefficient",
    )

    blk.reactor_capital_safety_factor = pyo.Var(
        initialize=2.5,
        units=pyo.units.dimensionless,
        doc="Reactor capital cost safety factor",
    )

    # see: https://magna-power.com/products/magnadc
    blk.power_supply_capital_slope = pyo.Var(
        initialize=0.51972,
        units=pyo.units.USD_2020 / pyo.units.watt,
        doc="DC power supply + transformer + electrical connection base cost",
    )

    # Flocculator costs rom Figure 5.5.22 in McGivney & Kawamura (2008)
    # linear equation refit to power equation (C = A*x**b)
    blk.floc_capital_cost_base = pyo.Var(
        initialize=1.0757e06,
        units=pyo.units.USD_2007,
        doc="Flocculator capital cost equation base parameter",
    )

    blk.floc_capital_cost_exponent = pyo.Var(
        initialize=9.5139e-01,
        units=pyo.units.dimensionless,
        doc="Flocculator capital cost equation exponent",
    )

    blk.sludge_handling_cost = pyo.Var(
        initialize=0,  # 0.025 $/kg from NMSU
        units=costing.base_currency / pyo.units.kg,
        doc="Cost per kg for handling generated solids",
    )

    blk.electrode_material_cost = pyo.Var(
        initialize=2,
        units=costing.base_currency / pyo.units.kg,
        doc="Cost of electrode per kg",
    )

    blk.electrode_material_cost_safety_factor = pyo.Var(
        initialize=2,
        units=pyo.units.dimensionless,
        doc="Electrode material cost safety factor",
    )

    blk.fix_all_vars()


@register_costing_parameter_block(
    build_rule=build_aluminum_cost_param_block,
    parameter_block_name="aluminum",
)
@register_costing_parameter_block(
    build_rule=build_iron_cost_param_block,
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

    blk.annual_sludge_flow = pyo.units.convert(
        sum(
            (
                ec.properties_byproduct[0].flow_mass_phase_comp["Liq", j]
                if j != "H2O"
                else 0
            )
            for j in ec.properties_byproduct[0].params.component_list
        ),
        to_units=pyo.units.kg / pyo.units.year,
    )

    base_currency = blk.costing_package.base_currency
    base_period = blk.costing_package.base_period

    blk.capital_cost_reactor = pyo.Var(
        initialize=1e4,
        units=base_currency,
        bounds=(0, None),
        doc="Cost of EC reactor",
    )

    blk.capital_cost_electrodes = pyo.Var(
        initialize=1e4,
        units=base_currency,
        bounds=(0, None),
        doc="Cost of EC electrodes",
    )

    blk.capital_cost_power_supply = pyo.Var(
        initialize=1e6,
        units=base_currency,
        bounds=(0, None),
        doc="Cost of EC power supply",
    )

    blk.capital_cost_floc_reactor = pyo.Var(
        initialize=1e4,
        units=base_currency,
        bounds=(0, None),
        doc="Cost of floc. basin",
    )

    blk.annual_sludge_management = pyo.Var(
        initialize=1e4,
        units=base_currency / pyo.units.year,
        bounds=(0, None),
        doc="Annual sludge management cost",
    )

    if ec.config.reactor_material == "carbon_steel":
        ec_params.reactor_material_coeff.fix(1)

    elif ec.config.reactor_material == "pvc":
        ec_params.reactor_material_coeff.fix(0.56)

    elif ec.config.reactor_material == "stainless_steel":
        ec_params.reactor_material_coeff.fix(3.4)

    if ec.config.electrode_material == "aluminum":
        ec_params.electrode_material_cost.fix(2.23)

    elif ec.config.electrode_material == "iron":
        ec_params.electrode_material_cost.fix(3.41)

    capital_cost_expr = 0

    ec_reactor_dim = pyo.units.convert(
        ec.reactor_volume * pyo.units.m**-3, to_units=pyo.units.dimensionless
    )
    blk.capital_cost_reactor_constraint = pyo.Constraint(
        expr=blk.capital_cost_reactor
        == pyo.units.convert(
            (
                (
                    ec_params.reactor_capital_cost_base
                    * ec_reactor_dim**ec_params.reactor_capital_cost_exponent
                )
                * ec_params.reactor_material_coeff
            )
            * ec_params.reactor_capital_safety_factor,
            to_units=base_currency,
        )
    )

    capital_cost_expr += blk.capital_cost_reactor

    blk.capital_cost_electrodes_constraint = pyo.Constraint(
        expr=blk.capital_cost_electrodes
        == pyo.units.convert(
            ec.electrode_mass * ec_params.electrode_material_cost,
            to_units=base_currency,
        )
        * ec_params.electrode_material_cost_safety_factor
    )

    capital_cost_expr += blk.capital_cost_electrodes

    blk.capital_cost_power_supply_constraint = pyo.Constraint(
        expr=blk.capital_cost_power_supply
        == pyo.units.convert(
            ec_params.power_supply_capital_slope * ec.power_required,
            to_units=base_currency,
        )
    )

    capital_cost_expr += blk.capital_cost_power_supply

    blk.capital_cost_floc_constraint = pyo.Constraint(
        expr=blk.capital_cost_floc_reactor
        == pyo.units.convert(
            ec_params.floc_capital_cost_base
            * pyo.units.convert(
                ec.floc_basin_vol * pyo.units.Mgallons**-1,
                to_units=pyo.units.dimensionless,
            )
            ** ec_params.floc_capital_cost_exponent,
            to_units=base_currency,
        )
    )
    capital_cost_expr += blk.capital_cost_floc_reactor

    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == pyo.units.convert(capital_cost_expr, to_units=base_currency)
    )

    blk.annual_sludge_management_constraint = pyo.Constraint(
        expr=blk.annual_sludge_management
        == pyo.units.convert(
            blk.annual_sludge_flow * ec_params.sludge_handling_cost,
            to_units=base_currency / base_period,
        )
    )

    blk.fixed_operating_cost_constraint = pyo.Constraint(
        expr=blk.fixed_operating_cost == blk.annual_sludge_management
    )

    blk.annual_electrode_replacement_mass_flow = pyo.Expression(
        expr=pyo.units.convert(
            ec.coagulant_dose * ec.properties_in[0].flow_vol_phase["Liq"],
            to_units=pyo.units.kg / pyo.units.year,
        )
    )

    blk.costing_package.cost_flow(
        blk.annual_electrode_replacement_mass_flow, ec.config.electrode_material
    )

    blk.costing_package.cost_flow(ec.power_required, "electricity")
