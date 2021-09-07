###############################################################################
# ProteusLib Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/nawi-hub/proteuslib/"
#
###############################################################################

from pyomo.environ import (
    Block, Constraint, Expression, Var, Param, Reals, NonNegativeReals, units as pyunits)
from idaes.core.util.exceptions import ConfigurationError
import proteuslib.flowsheets.RO_with_energy_recovery.financials as financials
import proteuslib.flowsheets.full_treatment_train.example_flowsheets.financials as financials

def build_costing_flowsheet(m, module=financials, **kwargs):
    '''
    m : model
    module : financials module
    '''

    # Add the costing parameter block
    module.add_costing_param_block(m.fs)

    # call get_costing for each unit model


    # call get_system_costing for whole flowsheet


if __name__ == "__main__":
    build_costing_flowsheet