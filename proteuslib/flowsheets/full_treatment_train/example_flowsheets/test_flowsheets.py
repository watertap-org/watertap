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
from proteuslib.flowsheets.full_treatment_train.example_flowsheets import flowsheet_limited, SepNF_0DRO


@pytest.mark.component
def test_flowsheet_limited_NF_no_bypass_1():
    m = flowsheet_limited.solve_build_flowsheet_limited_NF(has_bypass=False, NF_type='ZO', NF_base='ion')
    assert value(m.fs.RO.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.8350, rel=1e-3)
    assert value(m.fs.RO.inlet.flow_mass_phase_comp[0, 'Liq', 'TDS']) == pytest.approx(2.460e-2, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.4175, rel=1e-3)


@pytest.mark.component
def test_flowsheet_limited_NF_no_bypass_2():
    m = flowsheet_limited.solve_build_flowsheet_limited_NF(has_bypass=False, NF_type='Sep', NF_base='salt')
    assert value(m.fs.RO.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.8682, rel=1e-3)
    assert value(m.fs.RO.inlet.flow_mass_phase_comp[0, 'Liq', 'TDS']) == pytest.approx(2.658e-2, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.4341, rel=1e-3)


@pytest.mark.component
def test_flowsheet_limited_NF_no_bypass_3():
    m = flowsheet_limited.solve_build_flowsheet_limited_NF(has_bypass=False, NF_type='Sep', NF_base='ion')
    assert value(m.fs.RO.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.8682, rel=1e-3)
    assert value(m.fs.RO.inlet.flow_mass_phase_comp[0, 'Liq', 'TDS']) == pytest.approx(2.615e-2, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.4341, rel=1e-3)


@pytest.mark.component
def test_flowsheet_limited_NF_bypass_1():
    m = flowsheet_limited.solve_build_flowsheet_limited_NF(has_bypass=True, NF_type='Sep', NF_base='ion')
    assert value(m.fs.mixer.pretreatment.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.6512, rel=1e-3)
    assert value(m.fs.mixer.pretreatment.flow_mass_phase_comp[0, 'Liq', 'Na']) == pytest.approx(7.507e-3, rel=1e-3)
    assert value(m.fs.mixer.bypass.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.2412, rel=1e-3)
    assert value(m.fs.mixer.bypass.flow_mass_phase_comp[0, 'Liq', 'Cl']) == pytest.approx(5.079e-3, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.4462, rel=1e-3)


@pytest.mark.component
def test_flowsheet_limited_NF_bypass_2():
    m = flowsheet_limited.solve_build_flowsheet_limited_NF(has_bypass=True, NF_type='Sep', NF_base='salt')
    assert value(m.fs.mixer.pretreatment.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.6512, rel=1e-3)
    assert value(m.fs.mixer.pretreatment.flow_mass_phase_comp[0, 'Liq', 'NaCl']) == pytest.approx(1.908e-2, rel=1e-3)
    assert value(m.fs.mixer.bypass.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.2412, rel=1e-3)
    assert value(m.fs.mixer.bypass.flow_mass_phase_comp[0, 'Liq', 'NaCl']) == pytest.approx(7.068e-3, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.4462, rel=1e-3)


@pytest.mark.component
def test_SepNF_0DRO_simple_example():
    m = SepNF_0DRO.run_flowsheet_simple_example()
    assert value(m.fs.mixer.pretreatment.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.6512, rel=1e-3)
    assert value(m.fs.mixer.pretreatment.flow_mass_phase_comp[0, 'Liq', 'NaCl']) == pytest.approx(1.908e-2, rel=1e-3)
    assert value(m.fs.mixer.bypass.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.2412, rel=1e-3)
    assert value(m.fs.mixer.bypass.flow_mass_phase_comp[0, 'Liq', 'NaCl']) == pytest.approx(7.068e-3, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.5192, rel=1e-3)
