###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################

import pytest

from pyomo.environ import Block

from idaes.core.util.exceptions import ConfigurationError
from watertap.costing.watertap_costing_package import WaterTAPCostingData


def _get_config_testing_block():
    blk = Block()
    blk.unit_model = Block()
    return blk


@pytest.mark.component
def test_cost_reverse_osmosis_configuration_error():
    blk = _get_config_testing_block()
    with pytest.raises(
        ConfigurationError,
        match="unit_model received invalid argument for ro_type:"
        " foo. Argument must be a member of the ROType Enum.",
    ):
        WaterTAPCostingData.cost_reverse_osmosis(blk, ro_type="foo")


@pytest.mark.component
def test_cost_pump_configuration_error():
    blk = _get_config_testing_block()
    with pytest.raises(
        ConfigurationError,
        match="unit_model received invalid argument for pump_type:"
        " foo. Argument must be a member of the PumpType Enum.",
    ):
        WaterTAPCostingData.cost_pump(blk, pump_type="foo")


@pytest.mark.component
def test_cost_energy_recovery_device_configuration_error():
    blk = _get_config_testing_block()
    with pytest.raises(
        ConfigurationError,
        match="unit_model received invalid argument for energy_recovery_device_type:"
        " foo. Argument must be a member of the EnergyRecoveryDeviceType Enum.",
    ):
        WaterTAPCostingData.cost_energy_recovery_device(
            blk, energy_recovery_device_type="foo"
        )


@pytest.mark.component
def test_cost_mixer_configuration_error():
    blk = _get_config_testing_block()
    with pytest.raises(
        ConfigurationError,
        match="unit_model received invalid argument for mixer_type:"
        " foo. Argument must be a member of the MixerType Enum.",
    ):
        WaterTAPCostingData.cost_mixer(blk, mixer_type="foo")


@pytest.mark.component
def test_cost_crystallizer_configuration_error():
    blk = _get_config_testing_block()
    with pytest.raises(
        ConfigurationError,
        match="unit_model received invalid argument for cost_type:"
        " foo. Argument must be a member of the CrystallizerCostType Enum.",
    ):
        WaterTAPCostingData.cost_crystallizer(blk, cost_type="foo")
