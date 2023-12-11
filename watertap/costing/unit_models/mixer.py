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
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.misc import StrEnum
from ..util import cost_by_flow_volume, register_costing_parameter_block


class MixerType(StrEnum):
    default = "default"
    NaOCl = "NaOCl"
    CaOH2 = "CaOH2"


def build_mixer_cost_param_block(blk):

    blk.unit_cost = pyo.Var(
        initialize=361,
        doc="Mixer cost",
        units=pyo.units.USD_2018 / (pyo.units.liters / pyo.units.second),
    )


def cost_mixer(blk, mixer_type=MixerType.default, **kwargs):
    """
    Mixer costing method

    Args:
        mixer_type: MixerType Enum indicating mixer type,
            default = MixerType.default

        `**kwargs`: Additional keywords for the MixerType, e.g., NaOCl
            and CaOH2 mixers expect the `dosing_rate` keyword
            argument.
    """
    if mixer_type == MixerType.default:
        cost_default_mixer(blk, **kwargs)
    elif mixer_type == MixerType.NaOCl:
        cost_naocl_mixer(blk, **kwargs)
    elif mixer_type == MixerType.CaOH2:
        cost_caoh2_mixer(blk, **kwargs)
    else:
        raise ConfigurationError(
            f"{blk.unit_model.name} received invalid argument for mixer_type:"
            f" {mixer_type}. Argument must be a member of the MixerType Enum."
        )


@register_costing_parameter_block(
    build_rule=build_mixer_cost_param_block,
    parameter_block_name="mixer",
)
def cost_default_mixer(blk):
    """
    Default mixer costing method

    TODO: describe equations
    """
    cost_by_flow_volume(
        blk,
        blk.costing_package.mixer.unit_cost,
        pyo.units.convert(
            blk.unit_model.mixed_state[0].flow_vol,
            pyo.units.liter / pyo.units.second,
        ),
    )


def build_naocl_mixer_cost_param_block(blk):

    blk.unit_cost = pyo.Var(
        initialize=5.08,
        doc="NaOCl mixer cost",
        units=pyo.units.USD_2018 / (pyo.units.m**3 / pyo.units.day),
    )


def build_naocl_cost_param_block(blk):

    blk.cost = pyo.Param(
        initialize=0.23,
        doc="NaOCl cost",
        units=pyo.units.USD_2018 / pyo.units.kg,
    )
    blk.purity = pyo.Param(
        mutable=True,
        initialize=0.15,
        doc="NaOCl purity",
        units=pyo.units.dimensionless,
    )

    costing = blk.parent_block()
    costing.register_flow_type("NaOCl", blk.cost / blk.purity)


@register_costing_parameter_block(
    build_rule=build_naocl_cost_param_block, parameter_block_name="naocl"
)
@register_costing_parameter_block(
    build_rule=build_naocl_mixer_cost_param_block,
    parameter_block_name="naocl_mixer",
)
def cost_naocl_mixer(blk, dosing_rate):
    """
    NaOCl mixer costing method

    TODO: describe equations

    Args:
        dosing_rate: An expression in [mass/time] for NaOCl dosage
    """
    cost_by_flow_volume(
        blk,
        blk.costing_package.naocl_mixer.unit_cost,
        pyo.units.convert(
            blk.unit_model.inlet_stream_state[0].flow_vol,
            pyo.units.m**3 / pyo.units.day,
        ),
    )
    blk.costing_package.cost_flow(
        pyo.units.convert(dosing_rate, pyo.units.kg / pyo.units.s), "NaOCl"
    )


def build_caoh2_cost_param_block(blk):
    blk.cost = pyo.Param(
        mutable=True,
        initialize=0.12,
        doc="CaOH2 cost",
        units=pyo.units.USD_2018 / pyo.units.kg,
    )
    blk.purity = pyo.Param(
        mutable=True,
        initialize=1,
        doc="CaOH2 purity",
        units=pyo.units.dimensionless,
    )

    costing = blk.parent_block()
    costing.register_flow_type("CaOH2", blk.cost / blk.purity)


def build_caoh2_mixer_cost_param_block(blk):

    blk.unit_cost = pyo.Var(
        initialize=792.8 * 2.20462 / 2.0,
        doc="Ca(OH)2 mixer cost",
        units=pyo.units.USD_2018 / (pyo.units.kg / pyo.units.day),
    )


@register_costing_parameter_block(
    build_rule=build_caoh2_cost_param_block, parameter_block_name="caoh2"
)
@register_costing_parameter_block(
    build_rule=build_caoh2_mixer_cost_param_block,
    parameter_block_name="caoh2_mixer",
)
def cost_caoh2_mixer(blk, dosing_rate):
    """
    CaOH2 mixer costing method

    TODO: describe equations

    Args:
        dosing_rate: An expression in [mass/time] for CaOH2 dosage
    """
    stream = blk.unit_model.lime_stream
    blk.lime_kg_per_day = pyo.Expression(
        expr=pyo.units.convert(
            dosing_rate,
            pyo.units.kg / pyo.units.day,
        )
    )
    cost_by_flow_volume(
        blk,
        blk.costing_package.caoh2_mixer.unit_cost,
        blk.lime_kg_per_day,
    )
    blk.costing_package.cost_flow(
        pyo.units.convert(dosing_rate, pyo.units.kg / pyo.units.s), "CaOH2"
    )
