#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

# Import Pyomo libraries
from pyomo.environ import (
    Var,
    NonNegativeReals,
    NegativeReals,
    value,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In
# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    useDefault,
)
from idaes.core.util import scaling as iscale
from idaes.core.util.misc import add_object_reference
import idaes.logger as idaeslog

from watertap.core import (  # noqa # pylint: disable=unused-import
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    MembraneChannel1DBlock,
    PressureChangeType,
)
from watertap.core.membrane_channel1d import CONFIG_Template
from watertap.unit_models.pressure_changer import Pump
from watertap.unit_models.reverse_osmosis_1D import ReverseOsmosis1D, ReverseOsmosis1DData
from watertap.unit_models.reverse_osmosis_base import _add_has_full_reporting
from watertap.costing.unit_models.ccro import cost_ccro

__author__ = "Adam Atia, Bernard Knueven"

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("CCRO1D")
class CCRO1DData(ReverseOsmosis1DData):
    """
    CCRO 1D Unit Model Class.
    """
    # this model is meant to be identical to RO1D but with different costing

    CONFIG = CONFIG_Template()


    # CONFIG.declare(
    #     "multi",
    #     ConfigValue(
    #         # domain=Pump,
    #         default=None,
    #         description="",
    #         doc="""""",
    #     ),
    # )

    _add_has_full_reporting(CONFIG)
    
    def build(self):
        super().build()

        # if self.config.feed_pump is not None:
        #     self.feed_pump = add_object_reference(self, "feed_pump", self.config.feed_pump)
        
    # def calculate_scaling_factors(self):
    #     if iscale.get_scaling_factor(self.dens_solvent) is None:
    #         sf = iscale.get_scaling_factor(
    #             self.feed_side.properties[0, 0].dens_mass_phase["Liq"]
    #         )
    #         iscale.set_scaling_factor(self.dens_solvent, sf)

    #     super().calculate_scaling_factors()

    #     # these variables should have user input, if not there will be a warning
    #     if iscale.get_scaling_factor(self.width) is None:
    #         sf = iscale.get_scaling_factor(self.width, default=1, warning=True)
    #         iscale.set_scaling_factor(self.width, sf)

    #     if iscale.get_scaling_factor(self.length) is None:
    #         sf = iscale.get_scaling_factor(self.length, default=10, warning=True)
    #         iscale.set_scaling_factor(self.length, sf)

    #     for (t, x, p, j), v in self.mass_transfer_phase_comp.items():
    #         sf = (
    #             iscale.get_scaling_factor(
    #                 self.feed_side.properties[t, x].get_material_flow_terms(p, j)
    #             )
    #             / iscale.get_scaling_factor(self.feed_side.length)
    #         ) * value(self.nfe)
    #         if iscale.get_scaling_factor(v) is None:
    #             iscale.set_scaling_factor(v, sf)
    #         v = self.feed_side.mass_transfer_term[t, x, p, j]
    #         if iscale.get_scaling_factor(v) is None:
    #             iscale.set_scaling_factor(v, sf)

    #     if hasattr(self, "deltaP"):
    #         for v in self.deltaP.values():
    #             if iscale.get_scaling_factor(v) is None:
    #                 iscale.set_scaling_factor(v, 1e-4)

    @property
    def default_costing_method(self):
        return cost_ccro