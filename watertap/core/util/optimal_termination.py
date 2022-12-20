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

from pyomo.environ import check_optimal_termination
import idaes.logger as idaeslog

# Set up logger
_log = idaeslog.getLogger(__name__)


def optimal_termination(self):
    if check_optimal_termination(self):
        _log.warning(
            "The solver failed to converge to an optimal solution."
            "This suggests that the user provided infeasible inputs or that the model "
            "is poorly scaled, poorly initialized, or degenerate."
        )
