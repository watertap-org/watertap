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
from proteuslib.flowsheets.full_treatment_train.model_components import property_models


@pytest.mark.component
def test_feed_TDS():
    m = property_models.solve_specify_feed(base='TDS')
    assert value(m.fs.stream[0].mass_frac_phase_comp['Liq', 'TDS']) == pytest.approx(0.035, rel=1e-3)


@pytest.mark.component
def test_feed_ion():
    m = property_models.solve_specify_feed(base='ion')
    assert value(m.fs.stream[0].mass_frac_phase_comp['Liq', 'Ca']) == pytest.approx(3.82e-4, rel=1e-3)
    assert value(m.fs.stream[0].mass_frac_phase_comp['Liq', 'Cl']) == pytest.approx(2.032e-2, rel=1e-3)
    assert value(m.fs.stream[0].mass_frac_phase_comp['Liq', 'H2O']) == pytest.approx(0.9646, rel=1e-3)
    assert value(m.fs.stream[0].mass_frac_phase_comp['Liq', 'Mg']) == pytest.approx(1.394e-3, rel=1e-3)
    assert value(m.fs.stream[0].mass_frac_phase_comp['Liq', 'Na']) == pytest.approx(1.112e-2, rel=1e-3)
    assert value(m.fs.stream[0].mass_frac_phase_comp['Liq', 'SO4']) == pytest.approx(2.136e-3, rel=1e-3)


@pytest.mark.component
def test_feed_salt():
    m = property_models.solve_specify_feed(base='salt')
    assert value(m.fs.stream[0].mass_frac_phase_comp['Liq', 'CaSO4']) == pytest.approx(1.298e-3, rel=1e-3)
    assert value(m.fs.stream[0].mass_frac_phase_comp['Liq', 'H2O']) == pytest.approx(0.9646, rel=1e-3)
    assert value(m.fs.stream[0].mass_frac_phase_comp['Liq', 'MgCl2']) == pytest.approx(4.251e-3, rel=1e-3)
    assert value(m.fs.stream[0].mass_frac_phase_comp['Liq', 'MgSO4']) == pytest.approx(1.529e-3, rel=1e-3)
    assert value(m.fs.stream[0].mass_frac_phase_comp['Liq', 'NaCl']) == pytest.approx(2.827e-2, rel=1e-3)
