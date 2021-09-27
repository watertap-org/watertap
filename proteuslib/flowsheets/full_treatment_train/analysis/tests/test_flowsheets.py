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
                                                                 flowsheet_softening,
                                                                 flowsheet_softening_two_stage)

@pytest.mark.component
def test_flowsheet_NF():
    m = flowsheet_NF.optimize_flowsheet(system_recovery=0.5)
    assert (value(m.fs.tb_pretrt_to_desal.properties_in[0].flow_mass_phase_comp['Liq', 'H2O'])
            == pytest.approx(0.7629, rel=1e-3))
    assert (value(m.fs.tb_pretrt_to_desal.properties_in[0].flow_mass_phase_comp['Liq', 'Ca'])
            == pytest.approx(1.139e-4, rel=1e-3))


@pytest.mark.component
def test_flowsheet_NF_no_bypass():
    m = flowsheet_NF_no_bypass.solve_flowsheet()
    assert (value(m.fs.tb_pretrt_to_desal.properties_in[0].flow_mass_phase_comp['Liq', 'H2O'])
            == pytest.approx(0.29225, rel=1e-3))
    assert (value(m.fs.tb_pretrt_to_desal.properties_in[0].flow_mass_phase_comp['Liq', 'Ca'])
            == pytest.approx(2.413e-5, rel=1e-3))
    m.fs.NF.recovery_mass_phase_comp[0, 'Liq', 'H2O'].fix(0.5)
    flowsheet_NF_no_bypass.simulate(m, **{'unfix_nf_area': True})
    assert (value(m.fs.tb_pretrt_to_desal.properties_in[0].flow_mass_phase_comp['Liq', 'H2O'])
            == pytest.approx(0.4823, rel=1e-3))


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
    assert value(m.fs.costing.LCOW) == pytest.approx(0.5502, rel=1e-3)


@pytest.mark.component
def test_flowsheet_softening():
    m = flowsheet_softening.solve_flowsheet()
    assert value(m.fs.costing.LCOW) == pytest.approx(0.3837, rel=1e-3)


@pytest.mark.component
def test_flowsheet_softening_two_stage():
    desal_kwargs = flowsheet_two_stage.desal_kwargs
    m = flowsheet_softening_two_stage.optimize_flowsheet(system_recovery=0.80, **desal_kwargs)
    assert value(m.fs.costing.LCOW) == pytest.approx(0.9277, rel=1e-3)


@pytest.mark.component
def test_flowsheet_1DRO():
    desal_kwargs = flowsheet_two_stage.desal_kwargs
    desal_kwargs['RO_type'] = '1D'
    m = flowsheet_two_stage.optimize_flowsheet(system_recovery=0.70, **desal_kwargs)
    assert value(m.fs.costing.LCOW) == pytest.approx(0.5963, rel=1e-3)
