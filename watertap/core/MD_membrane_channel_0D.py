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

from pyomo.environ import (
    NegativeReals,
    Set,
    Var,
)
from idaes.core import (
    declare_process_block_class,
    FlowDirection,
)
from idaes.core.util import scaling as iscale
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.misc import add_object_reference
from idaes.core.base.control_volume0d import ControlVolume0DBlockData
import idaes.logger as idaeslog
from membrane_channel0d import MembraneChannel0DBlock
from temperature_polarization_mixn import (
    TemperaturePolarizationMixin,
    CONFIG_Template as Base_CONFIG_Template,
)

CONFIG_Template = Base_CONFIG_Template()


@declare_process_block_class("MDMembraneChannel0DBlock")
class MDMembraneChannel0DBlockData(
    TemperaturePolarizationMixin, MembraneChannel0DBlock
):
    def add_state_blocks(self, has_phase_equilibrium=False):

        super().add_state_blocks(has_phase_equilibrium)
        self._add_vapor_stateblock(has_phase_equilibrium)
