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
    NonNegativeReals,
    NegativeReals,
    value,
    check_optimal_termination,
    Set,
)
from pyomo.common.config import ConfigValue, In

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    useDefault,
)
from idaes.core.util.misc import add_object_reference
import idaes.logger as idaeslog

from watertap.core import MembraneChannel1D
from watertap.core.membrane_channel1d import CONFIG_Template
from watertap.unit_models.reverse_osmosis_base import ReverseOsmosisBaseData

__author__ = "Adam Atia"

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("ReverseOsmosis1D")
class ReverseOsmosis1DData(ReverseOsmosisBaseData):
    """Standard 1D Reverse Osmosis Unit Model Class."""

    CONFIG = CONFIG_Template()

    def _process_config(self):
        if self.config.transformation_method is useDefault:
            _log.warning(
                "Discretization method was "
                "not specified for the "
                "reverse osmosis module. "
                "Defaulting to finite "
                "difference method."
            )
            self.config.transformation_method = "dae.finite_difference"

        if self.config.transformation_scheme is useDefault:
            _log.warning(
                "Discretization scheme was "
                "not specified for the "
                "reverse osmosis module."
                "Defaulting to backward finite "
                "difference."
            )
            self.config.transformation_scheme = "BACKWARD"

    def _add_membrane_channel(self):
        # Build 1D Membrane Channel
        self.membrane_channel = MembraneChannel1D(
            default={
                "dynamic": self.config.dynamic,
                "has_holdup": self.config.has_holdup,
                "area_definition": self.config.area_definition,
                "property_package": self.config.property_package,
                "property_package_args": self.config.property_package_args,
                "transformation_method": self.config.transformation_method,
                "transformation_scheme": self.config.transformation_scheme,
                "finite_elements": self.config.finite_elements,
                "collocation_points": self.config.collocation_points,
            }
        )

    def build(self):
        """
        Build 1D RO model (pre-DAE transformation).

        Args:
            None

        Returns:
            None
        """
        # Check configuration errors
        self._process_config()

        # Call UnitModel.build to setup dynamics
        super().build()

        # ==========================================================================
        """ Add references to control volume geometry."""
