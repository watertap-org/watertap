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
import proteuslib.flowsheets.full_treatment_train.example_models.unit_separator as unit_separator
import proteuslib.flowsheets.full_treatment_train.example_models.unit_0DRO as unit_0DRO


@pytest.mark.component
def test_unit_separator_RO_example():
    m = unit_separator.run_RO_example()
    assert value(m.fs.RO.permeate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.4825, rel=1e-3)
    assert value(m.fs.RO.permeate.flow_mass_phase_comp[0, 'Liq', 'TDS']) == pytest.approx(3.5e-4, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.4825, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'TDS']) == pytest.approx(3.465e-2, rel=1e-3)


@pytest.mark.component
def test_unit_0DRO_simple():
    m = unit_0DRO.run_simple_example()
    assert value(m.fs.RO.permeate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.3389, rel=1e-3)
    assert value(m.fs.RO.permeate.flow_mass_phase_comp[0, 'Liq', 'TDS']) == pytest.approx(7.875e-5, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.6261, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'TDS']) == pytest.approx(3.492e-2, rel=1e-3)


@pytest.mark.component
def test_unit_0DRO_detailed():
    m = unit_0DRO.run_detailed_example()
    m.fs.RO.report()
    assert value(m.fs.RO.permeate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.2653, rel=1e-3)
    assert value(m.fs.RO.permeate.flow_mass_phase_comp[0, 'Liq', 'TDS']) == pytest.approx(8.476e-5, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.6997, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'TDS']) == pytest.approx(3.492e-2, rel=1e-3)
