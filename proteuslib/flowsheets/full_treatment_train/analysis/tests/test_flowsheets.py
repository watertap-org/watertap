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
from proteuslib.flowsheets.full_treatment_train.analysis import (flowsheet_NF,
                                                                 flowsheet_NF_no_bypass,
                                                                 flowsheet_single_stage,
                                                                 flowsheet_two_stage,
                                                                 flowsheet_NF_two_stage,
                                                                 flowsheet_softening_two_stage)

@pytest.mark.component
def test_flowsheet_NF():
    m = flowsheet_NF.solve_flowsheet(True)
    assert (value(m.fs.tb_pretrt_to_desal.properties_in[0].flow_mass_phase_comp['Liq', 'H2O'])
            == pytest.approx(0.7745, rel=1e-3))
    assert (value(m.fs.tb_pretrt_to_desal.properties_in[0].flow_mass_phase_comp['Liq', 'Ca'])
            == pytest.approx(2.151e-4, rel=1e-3))

# @pytest.mark.component
# def test_flowsheet_NF_no_bypass():
#     m = flowsheet_NF_no_bypass.solve_flowsheet()
#     assert (value(m.fs.tb_pretrt_to_desal.properties_in[0].flow_mass_phase_comp['Liq', 'H2O'])
#             == pytest.approx(0.2923, rel=1e-3))
#     assert (value(m.fs.tb_pretrt_to_desal.properties_in[0].flow_mass_phase_comp['Liq', 'Ca'])
#             == pytest.approx(2.414e-5, rel=1e-3))

@pytest.mark.component
def test_flowsheet_single_stage():
    desal_kwargs = flowsheet_single_stage.desal_kwargs
    m = flowsheet_single_stage.optimize_flowsheet(system_recovery=0.5, **desal_kwargs)
    assert value(m.fs.costing.LCOW) == pytest.approx(0.4996, rel=1e-3)


@pytest.mark.component
def test_flowsheet_two_stage():
    desal_kwargs = flowsheet_two_stage.desal_kwargs
    m = flowsheet_two_stage.optimize_flowsheet(system_recovery=0.65, **desal_kwargs)
    assert value(m.fs.costing.LCOW) == pytest.approx(0.5435, rel=1e-3)


@pytest.mark.component
def test_flowsheet_NF_two_stage():
    desal_kwargs = flowsheet_two_stage.desal_kwargs
    m = flowsheet_NF_two_stage.optimize_flowsheet(system_recovery=0.70, **desal_kwargs)
    # No longer true with NF.area.lb = 0.1 and pseudo-equality constraint
    # assert value(m.fs.costing.LCOW) == pytest.approx(0.5518, rel=1e-3)
    assert value(m.fs.costing.LCOW) == pytest.approx(0.5502, rel=1e-3)


@pytest.mark.component
def test_flowsheet_softening_two_stage():
    desal_kwargs = flowsheet_two_stage.desal_kwargs
    m = flowsheet_softening_two_stage.optimize_flowsheet(system_recovery=0.80, **desal_kwargs)
    assert value(m.fs.costing.LCOW) == pytest.approx(0.9277, rel=1e-3)
