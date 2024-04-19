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

from idaes.core.util.exceptions import InitializationError


class InitializationMixin:
    """
    Class for catching `InitializationError` in UnitModel.initialize
    """

    def initialize(self, *args, **kwargs):
        try:
            return super().initialize(*args, **kwargs)
        except InitializationError:
            for blk in self._initialization_order:
                blk.activate()
            raise
