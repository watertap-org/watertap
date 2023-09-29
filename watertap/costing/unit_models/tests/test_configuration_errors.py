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

import pytest

import pyomo.environ as pyo

from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.misc import add_object_reference
import idaes.core as idc

from watertap.costing.watertap_costing_package import WaterTAPCosting

from watertap.costing.unit_models.reverse_osmosis import cost_reverse_osmosis
from watertap.costing.unit_models.pump import cost_pump
from watertap.costing.unit_models.energy_recovery_device import (
    cost_energy_recovery_device,
)
from watertap.costing.unit_models.mixer import cost_mixer
from watertap.costing.unit_models.crystallizer import cost_crystallizer


def _get_config_testing_block():
    m = pyo.ConcreteModel()
    m.fs = idc.FlowsheetBlock(dynamic=False)
    m.fs.costing = WaterTAPCosting()

    m.fs.unit = pyo.Block()
    m.fs.unit.costing = pyo.Block()
    add_object_reference(m.fs.unit.costing, "costing_package", m.fs.costing)
    add_object_reference(m.fs.unit.costing, "unit_model", m.fs.unit)

    return m.fs.unit.costing


@pytest.mark.component
def test_cost_reverse_osmosis_configuration_error():
    blk = _get_config_testing_block()
    with pytest.raises(
        ConfigurationError,
        match="unit received invalid argument for ro_type:"
        " foo. Argument must be a member of the ROType Enum.",
    ):
        cost_reverse_osmosis(blk, ro_type="foo")


@pytest.mark.component
def test_cost_pump_configuration_error():
    blk = _get_config_testing_block()
    with pytest.raises(
        ConfigurationError,
        match="unit received invalid argument for pump_type:"
        " foo. Argument must be a member of the PumpType Enum.",
    ):
        cost_pump(blk, pump_type="foo")


@pytest.mark.component
def test_cost_energy_recovery_device_configuration_error():
    blk = _get_config_testing_block()
    with pytest.raises(
        ConfigurationError,
        match="unit received invalid argument for energy_recovery_device_type:"
        " foo. Argument must be a member of the EnergyRecoveryDeviceType Enum.",
    ):
        cost_energy_recovery_device(blk, energy_recovery_device_type="foo")


@pytest.mark.component
def test_cost_mixer_configuration_error():
    blk = _get_config_testing_block()
    with pytest.raises(
        ConfigurationError,
        match="unit received invalid argument for mixer_type:"
        " foo. Argument must be a member of the MixerType Enum.",
    ):
        cost_mixer(blk, mixer_type="foo")


@pytest.mark.component
def test_cost_crystallizer_configuration_error():
    blk = _get_config_testing_block()
    with pytest.raises(
        ConfigurationError,
        match="unit received invalid argument for cost_type:"
        " foo. Argument must be a member of the CrystallizerCostType Enum.",
    ):
        cost_crystallizer(blk, cost_type="foo")
