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
import pytest
from pyomo.environ import value
from proteuslib.flowsheets.full_treatment_train.analysis import (flowsheet_single_stage,
                                                                 flowsheet_two_stage,
                                                                 )

@pytest.mark.component
def test_flowsheet_single_stage():
    desal_kwargs = flowsheet_single_stage.desal_kwargs
    m = flowsheet_single_stage.optimize_flowsheet(system_recovery=0.5, **desal_kwargs)
    assert value(m.fs.costing.LCOW) == pytest.approx(0.5036, rel=1e-3)


@pytest.mark.component
def test_flowsheet_two_stage():
    desal_kwargs = flowsheet_two_stage.desal_kwargs
    m = flowsheet_two_stage.optimize_flowsheet(system_recovery=0.65, **desal_kwargs)
    assert False
    assert value(m.fs.costing.LCOW) == pytest.approx(0.5036, rel=1e-3)
