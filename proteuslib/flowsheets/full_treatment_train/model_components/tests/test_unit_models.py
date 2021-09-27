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
from proteuslib.flowsheets.full_treatment_train.model_components import unit_separator, unit_0DRO, unit_1DRO, unit_ZONF


@pytest.mark.component
def test_unit_separator_SepRO():
    m = unit_separator.solve_SepRO(base='TDS')
    assert value(m.fs.RO.permeate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.4825, rel=1e-3)
    assert value(m.fs.RO.permeate.flow_mass_phase_comp[0, 'Liq', 'TDS']) == pytest.approx(3.5e-4, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.4825, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'TDS']) == pytest.approx(3.465e-2, rel=1e-3)


@pytest.mark.component
def test_unit_separator_SepNF_ion():
    m = unit_separator.solve_SepNF(base='ion')
    assert value(m.fs.NF.mixed_state[0].mass_frac_phase_comp['Liq', 'Na']) == pytest.approx(11122e-6, rel=1e-3)
    assert value(m.fs.NF.mixed_state[0].mass_frac_phase_comp['Liq', 'Cl']) == pytest.approx(20317e-6, rel=1e-3)
    assert value(m.fs.NF.split_fraction[0, 'permeate', 'Cl']) == pytest.approx(0.7753, rel=1e-3)
    assert value(m.fs.NF.permeate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.8682, rel=1e-3)
    assert value(m.fs.NF.permeate.flow_mass_phase_comp[0, 'Liq', 'Na']) == pytest.approx(1.000e-2, rel=1e-3)
    assert value(m.fs.NF.permeate.flow_mass_phase_comp[0, 'Liq', 'Cl']) == pytest.approx(1.575e-2, rel=1e-3)
    assert value(m.fs.NF.retentate.flow_mass_phase_comp[0, 'Liq', 'Ca']) == pytest.approx(3.438e-4, rel=1e-3)
    assert value(m.fs.NF.retentate.flow_mass_phase_comp[0, 'Liq', 'Mg']) == pytest.approx(1.255e-3, rel=1e-3)
    assert value(m.fs.NF.retentate.flow_mass_phase_comp[0, 'Liq', 'SO4']) == pytest.approx(1.922e-3, rel=1e-3)


@pytest.mark.component
def test_unit_separator_SepNF_salt():
    m = unit_separator.solve_SepNF(base='salt')
    assert value(m.fs.NF.permeate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.8682, rel=1e-3)
    assert value(m.fs.NF.permeate.flow_mass_phase_comp[0, 'Liq', 'NaCl']) == pytest.approx(2.544e-2, rel=1e-3)
    assert value(m.fs.NF.retentate.flow_mass_phase_comp[0, 'Liq', 'CaSO4']) == pytest.approx(1.168e-3, rel=1e-3)
    assert value(m.fs.NF.retentate.flow_mass_phase_comp[0, 'Liq', 'MgSO4']) == pytest.approx(1.376e-3, rel=1e-3)
    assert value(m.fs.NF.retentate.flow_mass_phase_comp[0, 'Liq', 'MgCl2']) == pytest.approx(3.401e-3, rel=1e-3)
    assert value(m.fs.NF.retentate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(9.645e-2, rel=1e-3)


@pytest.mark.component
def test_unit_0DRO_simple():
    m = unit_0DRO.solve_RO(base='TDS', level='simple')
    assert value(m.fs.RO.permeate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.3389, rel=1e-3)
    assert value(m.fs.RO.permeate.flow_mass_phase_comp[0, 'Liq', 'TDS']) == pytest.approx(7.875e-5, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.6261, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'TDS']) == pytest.approx(3.492e-2, rel=1e-3)


@pytest.mark.component
def test_unit_0DRO_detailed():
    m = unit_0DRO.solve_RO(base='TDS', level='detailed')
    assert value(m.fs.RO.permeate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.2653, rel=1e-3)
    assert value(m.fs.RO.permeate.flow_mass_phase_comp[0, 'Liq', 'TDS']) == pytest.approx(8.476e-5, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.6997, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'TDS']) == pytest.approx(3.492e-2, rel=1e-3)


@pytest.mark.component
def test_unit_1DRO_detailed():
    m = unit_1DRO.solve_RO(base='TDS', level='detailed')
    assert value(m.fs.RO.mixed_permeate[0].flow_mass_phase_comp['Liq', 'H2O']) == pytest.approx(0.2541, rel=1e-3)
    assert value(m.fs.RO.mixed_permeate[0].flow_mass_phase_comp['Liq', 'TDS']) == pytest.approx(8.577e-5, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.7109, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'TDS']) == pytest.approx(3.491e-2, rel=1e-3)


@pytest.mark.component
def test_unit_ZONF():
    m = unit_ZONF.solve_ZONF(base='ion')
    assert value(m.fs.NF.permeate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.8350, rel=1e-3)
    assert value(m.fs.NF.permeate.flow_mass_phase_comp[0, 'Liq', 'Ca']) == pytest.approx(6.897e-5, rel=1e-3)
    assert value(m.fs.NF.permeate.flow_mass_phase_comp[0, 'Liq', 'Cl']) == pytest.approx(1.475e-2, rel=1e-3)
