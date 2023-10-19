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

import pyomo.environ as pyo
from idaes.core.util.constants import Constants
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.misc import StrEnum
from ..util import register_costing_parameter_block, make_capital_cost_var


class CrystallizerCostType(StrEnum):
    default = "default"
    mass_basis = "mass_basis"
    volume_basis = "volume_basis"


def build_crystallizer_cost_param_block(blk):

    blk.steam_pressure = pyo.Var(
        initialize=3,
        units=pyo.units.bar,
        doc="Steam pressure (gauge) for crystallizer heating: 3 bar default based on Dutta example",
    )

    blk.efficiency_pump = pyo.Var(
        initialize=0.7,
        units=pyo.units.dimensionless,
        doc="Crystallizer pump efficiency - assumed",
    )

    blk.pump_head_height = pyo.Var(
        initialize=1,
        units=pyo.units.m,
        doc="Crystallizer pump head height -  assumed, unvalidated",
    )

    # Crystallizer operating cost information from literature
    blk.fob_unit_cost = pyo.Var(
        initialize=675000,
        doc="Forced circulation crystallizer reference free-on-board cost (Woods, 2007)",
        units=pyo.units.USD_2007,
    )

    blk.ref_capacity = pyo.Var(
        initialize=1,
        doc="Forced circulation crystallizer reference crystal capacity (Woods, 2007)",
        units=pyo.units.kg / pyo.units.s,
    )

    blk.ref_exponent = pyo.Var(
        initialize=0.53,
        doc="Forced circulation crystallizer cost exponent factor (Woods, 2007)",
        units=pyo.units.dimensionless,
    )

    blk.iec_percent = pyo.Var(
        initialize=1.43,
        doc="Forced circulation crystallizer installed equipment cost (Diab and Gerogiorgis, 2017)",
        units=pyo.units.dimensionless,
    )

    blk.volume_cost = pyo.Var(
        initialize=16320,
        doc="Forced circulation crystallizer cost per volume (Yusuf et al., 2019)",
        units=pyo.units.USD_2007,  ## TODO: Needs confirmation, but data is from Perry apparently
    )

    blk.vol_basis_exponent = pyo.Var(
        initialize=0.47,
        doc="Forced circulation crystallizer volume-based cost exponent (Yusuf et al., 2019)",
        units=pyo.units.dimensionless,
    )

    blk.steam_cost = pyo.Var(
        initialize=0.004,
        units=pyo.units.USD_2018 / (pyo.units.meter**3),
        doc="Steam cost, Panagopoulos (2019)",
    )

    blk.NaCl_recovery_value = pyo.Var(
        initialize=0,
        units=pyo.units.USD_2018 / pyo.units.kg,
        doc="Unit recovery value of NaCl",
    )

    costing = blk.parent_block()
    costing.register_flow_type("steam", blk.steam_cost)
    costing.register_flow_type("NaCl", blk.NaCl_recovery_value)


def cost_crystallizer(blk, cost_type=CrystallizerCostType.default):
    """
    Function for costing the FC crystallizer by the mass flow of produced crystals.
    The operating cost model assumes that heat is supplied via condensation of saturated steam (see Dutta et al.)

    Args:
        cost_type: Option for crystallizer cost function type - volume or mass basis
    """
    if (
        cost_type == CrystallizerCostType.default
        or cost_type == CrystallizerCostType.mass_basis
    ):
        cost_crystallizer_by_crystal_mass(blk)
    elif cost_type == CrystallizerCostType.volume_basis:
        cost_crystallizer_by_volume(blk)
    else:
        raise ConfigurationError(
            f"{blk.unit_model.name} received invalid argument for cost_type:"
            f" {cost_type}. Argument must be a member of the CrystallizerCostType Enum."
        )


def _cost_crystallizer_flows(blk):
    blk.costing_package.cost_flow(
        pyo.units.convert(
            (
                blk.unit_model.magma_circulation_flow_vol
                * blk.unit_model.dens_mass_slurry
                * Constants.acceleration_gravity
                * blk.costing_package.crystallizer.pump_head_height
                / blk.costing_package.crystallizer.efficiency_pump
            ),
            to_units=pyo.units.kW,
        ),
        "electricity",
    )

    blk.costing_package.cost_flow(
        pyo.units.convert(
            (blk.unit_model.work_mechanical[0] / _compute_steam_properties(blk)),
            to_units=pyo.units.m**3 / pyo.units.s,
        ),
        "steam",
    )

    blk.costing_package.cost_flow(
        blk.unit_model.solids.flow_mass_phase_comp[0, "Sol", "NaCl"],
        "NaCl",
    )


@register_costing_parameter_block(
    build_rule=build_crystallizer_cost_param_block,
    parameter_block_name="crystallizer",
)
def cost_crystallizer_by_crystal_mass(blk):
    """
    Mass-based capital cost for FC crystallizer
    """
    make_capital_cost_var(blk)
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == pyo.units.convert(
            (
                blk.costing_package.crystallizer.iec_percent
                * blk.costing_package.crystallizer.fob_unit_cost
                * (
                    sum(
                        blk.unit_model.solids.flow_mass_phase_comp[0, "Sol", j]
                        for j in blk.unit_model.config.property_package.solute_set
                    )
                    / blk.costing_package.crystallizer.ref_capacity
                )
                ** blk.costing_package.crystallizer.ref_exponent
            ),
            to_units=blk.costing_package.base_currency,
        )
    )
    _cost_crystallizer_flows(blk)


@register_costing_parameter_block(
    build_rule=build_crystallizer_cost_param_block,
    parameter_block_name="crystallizer",
)
def cost_crystallizer_by_volume(blk):
    """
    Volume-based capital cost for FC crystallizer
    """
    make_capital_cost_var(blk)
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == pyo.units.convert(
            (
                blk.costing_package.crystallizer.volume_cost
                * (
                    (
                        pyo.units.convert(
                            blk.unit_model.volume_suspension
                            * (
                                blk.unit_model.height_crystallizer
                                / blk.unit_model.height_slurry
                            ),
                            to_units=(pyo.units.ft) ** 3,
                        )
                    )
                    / pyo.units.ft**3
                )
                ** blk.costing_package.crystallizer.vol_basis_exponent
            ),
            to_units=blk.costing_package.base_currency,
        )
    )
    _cost_crystallizer_flows(blk)


def _compute_steam_properties(blk):
    """
    Function for computing saturated steam properties for thermal heating estimation.

    Args:
        pressure_sat:   Steam gauge pressure in bar

    Out:
        Steam thermal capacity (latent heat of condensation * density) in kJ/m3
    """
    pressure_sat = blk.costing_package.crystallizer.steam_pressure
    # 1. Compute saturation temperature of steam: computed from El-Dessouky expression
    tsat_constants = [
        42.6776 * pyo.units.K,
        -3892.7 * pyo.units.K,
        1000 * pyo.units.kPa,
        -9.48654 * pyo.units.dimensionless,
    ]
    psat = (
        pyo.units.convert(pressure_sat, to_units=pyo.units.kPa)
        + 101.325 * pyo.units.kPa
    )
    temperature_sat = tsat_constants[0] + tsat_constants[1] / (
        pyo.log(psat / tsat_constants[2]) + tsat_constants[3]
    )

    # 2. Compute latent heat of condensation/vaporization: computed from Sharqawy expression
    t = temperature_sat - 273.15 * pyo.units.K
    enth_mass_units = pyo.units.J / pyo.units.kg
    t_inv_units = pyo.units.K**-1
    dh_constants = [
        2.501e6 * enth_mass_units,
        -2.369e3 * enth_mass_units * t_inv_units**1,
        2.678e-1 * enth_mass_units * t_inv_units**2,
        -8.103e-3 * enth_mass_units * t_inv_units**3,
        -2.079e-5 * enth_mass_units * t_inv_units**4,
    ]
    dh_vap = (
        dh_constants[0]
        + dh_constants[1] * t
        + dh_constants[2] * t**2
        + dh_constants[3] * t**3
        + dh_constants[4] * t**4
    )
    dh_vap = pyo.units.convert(dh_vap, to_units=pyo.units.kJ / pyo.units.kg)

    # 3. Compute specific volume: computed from Affandi expression (Eq 5)
    t_critical = 647.096 * pyo.units.K
    t_red = temperature_sat / t_critical  # Reduced temperature
    sp_vol_constants = [
        -7.75883 * pyo.units.dimensionless,
        3.23753 * pyo.units.dimensionless,
        2.05755 * pyo.units.dimensionless,
        -0.06052 * pyo.units.dimensionless,
        0.00529 * pyo.units.dimensionless,
    ]
    log_sp_vol = (
        sp_vol_constants[0]
        + sp_vol_constants[1] * (pyo.log(1 / t_red)) ** 0.4
        + sp_vol_constants[2] / (t_red**2)
        + sp_vol_constants[3] / (t_red**4)
        + sp_vol_constants[4] / (t_red**5)
    )
    sp_vol = pyo.exp(log_sp_vol) * pyo.units.m**3 / pyo.units.kg

    # 4. Return specific energy: density * latent heat
    return dh_vap / sp_vol
