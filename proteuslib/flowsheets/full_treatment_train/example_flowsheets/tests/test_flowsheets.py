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
from proteuslib.flowsheets.full_treatment_train.example_flowsheets import flowsheet_limited


@pytest.mark.component
def test_flowsheet_limited_NF_no_bypass_1():
    m = flowsheet_limited.solve_build_flowsheet_limited_NF(
        has_bypass=False, has_desal_feed=False, is_twostage=False,
        NF_type='Sep', NF_base='ion',
        RO_type='Sep', RO_base='TDS', RO_level='simple')
    assert value(m.fs.NF.retentate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(9.646e-2, rel=1e-3)
    assert value(m.fs.NF.retentate.flow_mass_phase_comp[0, 'Liq', 'Ca']) == pytest.approx(3.438e-4, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.4341, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'TDS']) == pytest.approx(2.589e-2, rel=1e-3)


@pytest.mark.component
def test_flowsheet_limited_NF_no_bypass_2():
    m = flowsheet_limited.solve_build_flowsheet_limited_NF(
        has_bypass=False, has_desal_feed=False, is_twostage=False,
        NF_type='Sep', NF_base='salt',
        RO_type='Sep', RO_base='TDS', RO_level='simple')
    assert value(m.fs.NF.retentate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(9.646e-2, rel=1e-3)
    assert value(m.fs.NF.retentate.flow_mass_phase_comp[0, 'Liq', 'CaSO4']) == pytest.approx(1.168e-3, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.4341, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'TDS']) == pytest.approx(2.631e-2, rel=1e-3)


def test_flowsheet_limited_NF_no_bypass_3():
    m = flowsheet_limited.solve_build_flowsheet_limited_NF(
        has_bypass=False, has_desal_feed=False, is_twostage=False,
        NF_type='ZO', NF_base='ion',
        RO_type='Sep', RO_base='TDS', RO_level='simple')
    assert value(m.fs.NF.retentate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.1296, rel=1e-3)
    assert value(m.fs.NF.retentate.flow_mass_phase_comp[0, 'Liq', 'Ca']) == pytest.approx(3.130e-4, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.4175, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'TDS']) == pytest.approx(2.435e-2, rel=1e-3)


def test_flowsheet_limited_NF_bypass_1():
    m = flowsheet_limited.solve_build_flowsheet_limited_NF(
        has_bypass=True, has_desal_feed=False, is_twostage=False,
        NF_type='Sep', NF_base='ion',
        RO_type='Sep', RO_base='TDS', RO_level='simple')
    assert value(m.fs.NF.retentate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(8.682e-2, rel=1e-3)
    assert value(m.fs.NF.retentate.flow_mass_phase_comp[0, 'Liq', 'Ca']) == pytest.approx(3.094e-4, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.4389, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'TDS']) == pytest.approx(2.680e-2, rel=1e-3)


def test_flowsheet_limited_NF_bypass_2():
    m = flowsheet_limited.solve_build_flowsheet_limited_NF(
        has_bypass=True, has_desal_feed=False, is_twostage=False,
        NF_type='Sep', NF_base='ion',
        RO_type='0D', RO_base='TDS', RO_level='simple')
    assert value(m.fs.NF.retentate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(8.682e-2, rel=1e-3)
    assert value(m.fs.NF.retentate.flow_mass_phase_comp[0, 'Liq', 'Ca']) == pytest.approx(3.094e-4, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.4916, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'TDS']) == pytest.approx(2.700e-2, rel=1e-3)


def test_flowsheet_limited_NF_bypass_3():
    m = flowsheet_limited.solve_build_flowsheet_limited_NF(
        has_bypass=True, has_desal_feed=False, is_twostage=False,
        NF_type='ZO', NF_base='ion',
        RO_type='0D', RO_base='TDS', RO_level='detailed')
    assert value(m.fs.NF.retentate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(3.318e-2, rel=1e-3)
    assert value(m.fs.NF.retentate.flow_mass_phase_comp[0, 'Liq', 'Ca']) == pytest.approx(2.748e-4, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.6045, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'TDS']) == pytest.approx(2.805e-2, rel=1e-3)


def test_flowsheet_limited_NF_bypass_twostage_1():
    m = flowsheet_limited.solve_build_flowsheet_limited_NF(
        has_bypass=True, has_desal_feed=False, is_twostage=True,
        NF_type='ZO', NF_base='ion',
        RO_type='0D', RO_base='TDS', RO_level='detailed')
    assert value(m.fs.NF.retentate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(3.318e-2, rel=1e-3)
    assert value(m.fs.NF.retentate.flow_mass_phase_comp[0, 'Liq', 'Ca']) == pytest.approx(2.748e-4, rel=1e-3)
    assert value(m.fs.RO2.retentate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.4229, rel=1e-3)
    assert value(m.fs.RO2.retentate.flow_mass_phase_comp[0, 'Liq', 'TDS']) == pytest.approx(2.793e-2, rel=1e-3)
