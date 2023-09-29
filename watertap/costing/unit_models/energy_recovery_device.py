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


class EnergyRecoveryDeviceType(StrEnum):
    pressure_exchanger = "pressure_exchanger"


def cost_energy_recovery_device(
    blk,
    energy_recovery_device_type=EnergyRecoveryDeviceType.pressure_exchanger,
    cost_electricity_flow=True,
):
    """
    Energy recovery device costing method

    TODO: describe equations

    Args:
        energy_recovery_device_type: EnergyRecoveryDeviceType Enum indicating ERD type.
            Defaults to EnergyRecoveryDeviceType.pressure_exchanger.
        cost_electricity_flow: bool, if True, the ERD's work_mechanical will
            be converted to kW and costed as an electricity. Defaults to True.
    """
    if energy_recovery_device_type == EnergyRecoveryDeviceType.pressure_exchanger:
        cost_pressure_exchanger_erd(blk)
    else:
        raise ConfigurationError(
            f"{blk.unit_model.name} received invalid argument for energy_recovery_device_type:"
            f" {energy_recovery_device_type}. Argument must be a member of the EnergyRecoveryDeviceType Enum."
        )


def build_energy_recovery_device_cost_param_block(blk):

    blk.pressure_exchanger_cost = pyo.Var(
        initialize=535,
        doc="Pressure exchanger cost",
        units=pyo.units.USD_2018 / (pyo.units.meter**3 / pyo.units.hours),
    )


@register_costing_parameter_block(
    build_rule=build_energy_recovery_device_cost_param_block,
    parameter_block_name="energy_recovery_device",
)
def cost_pressure_exchanger_erd(blk, cost_electricity_flow=True):
    """
    ERD pressure exchanger costing method

    TODO: describe equations

    Args:
        cost_electricity_flow (:obj:`bool`): if True, the ERD's work_mechanical will
            be converted to kW and costed as an electricity. Defaults to True.
    """
    t0 = blk.flowsheet().time.first()
    cost_by_flow_volume(
        blk,
        blk.costing_package.energy_recovery_device.pressure_exchanger_cost,
        pyo.units.convert(
            blk.unit_model.control_volume.properties_in[t0].flow_vol,
            (pyo.units.meter**3 / pyo.units.hours),
        ),
    )
    if cost_electricity_flow:
        blk.costing_package.cost_flow(
            pyo.units.convert(
                blk.unit_model.work_mechanical[t0], to_units=pyo.units.kW
            ),
            "electricity",
        )
