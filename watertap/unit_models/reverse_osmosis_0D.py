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

# Import Pyomo libraries
from pyomo.environ import (
    Var,
    Set,
    NonNegativeReals,
    NegativeReals,
    Reference,
    units as pyunits,
    exp,
    value,
    check_optimal_termination,
)

from idaes.core import declare_process_block_class
from idaes.core.util.misc import add_object_reference
from watertap.core import MembraneChannel0D, ConcentrationPolarizationType, MassTransferCoefficient, PressureChangeType
from watertap.unit_models.reverse_osmosis_base import ReverseOsmosisBaseData
import idaes.logger as idaeslog


__author__ = "Tim Bartholomew, Adam Atia"


@declare_process_block_class("ReverseOsmosis0D")
class ReverseOsmosisData(ReverseOsmosisBaseData):
    """
    Standard RO Unit Model Class:
    - zero dimensional model
    - steady state only
    - single liquid phase only
    """

    CONFIG = ReverseOsmosisBaseData.CONFIG()

    def _add_membrane_channel(self):
        # Build membrane channel control volume
        self.membrane_channel = MembraneChannel0D(
            default={
                "dynamic": False,
                "has_holdup": False,
                "property_package": self.config.property_package,
                "property_package_args": self.config.property_package_args,
            }
        )
