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
"""
Feed block for zero-order flowsheets with methods for getting concentration
data from database
"""

from pyomo.environ import units as pyunits

from idaes.generic_models.unit_models.feed import FeedData
from idaes.core import declare_process_block_class
import idaes.logger as idaeslog

# Some more inforation about this module
__author__ = "Andrew Lee"

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("FeedZO")
class FeedZOData(FeedData):
    """
    Zero-Order feed block.
    """

    CONFIG = FeedData.CONFIG()

    def build(self):
        super().build()

    def load_feed_data_from_database(self, overwrite=False):
        # Get database and water source from property package
        db = self.config.property_package.config.database
        water_source = self.config.property_package.config.water_source

        # Get feed data from database
        data = db.get_source_data(water_source)

        for t in self.flowsheet().time:
            if overwrite or not self.outlet.flow_vol[t].fixed:
                try:
                    val = data["default_flow"]["value"]
                    units = getattr(pyunits, data["default_flow"]["units"])
                    self.outlet.flow_vol[t].fix(val*units)
                except KeyError:
                    _log.info(
                        f"{self.name} no default flowrate was definined "
                        f"in database water source. Value was not fixed.")

        for (t, j), v in self.outlet.conc_mass_comp.items():
            if overwrite or not v.fixed:
                try:
                    val = data["solutes"][j]["value"]
                    units = getattr(pyunits, data["solutes"][j]["units"])
                    v.fix(val*units)
                except KeyError:
                    _log.info(f"{self.name} component {j} is not defined in "
                              f"database water source. Value was not fixed.")
