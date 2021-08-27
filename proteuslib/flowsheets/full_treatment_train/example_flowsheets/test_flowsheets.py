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
import proteuslib.flowsheets.full_treatment_train.example_flowsheets.SepNF_SepRO_NoBypass as SepNF_SepRO_NoBypass
import proteuslib.flowsheets.full_treatment_train.example_flowsheets.SepNF_SepRO as SepNF_SepRO


@pytest.mark.component
def test_SepNF_SepRO_NoBypass_NF_salt_basis_example():
    m = SepNF_SepRO_NoBypass.run_flowsheet_example('salt')
    assert value(m.fs.RO.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.8682, rel=1e-3)
    assert value(m.fs.RO.inlet.flow_mass_phase_comp[0, 'Liq', 'TDS']) == pytest.approx(2.658e-2, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.4341, rel=1e-3)


@pytest.mark.component
def test_SepNF_SepRO_NoBypass_NF_ion_basis_example():
    m = SepNF_SepRO_NoBypass.run_flowsheet_example('ion')
    assert value(m.fs.RO.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.8682, rel=1e-3)
    assert value(m.fs.RO.inlet.flow_mass_phase_comp[0, 'Liq', 'TDS']) == pytest.approx(2.615e-2, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.4341, rel=1e-3)


@pytest.mark.component
def test_SepNF_SepRO_NF_salt_basis_example():
    m = SepNF_SepRO.run_flowsheet_example('salt')
    assert value(m.fs.mixer.pretreatment.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.6512, rel=1e-3)
    assert value(m.fs.mixer.pretreatment.flow_mass_phase_comp[0, 'Liq', 'NaCl']) == pytest.approx(1.908e-2, rel=1e-3)
    assert value(m.fs.mixer.bypass.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.2412, rel=1e-3)
    assert value(m.fs.mixer.bypass.flow_mass_phase_comp[0, 'Liq', 'NaCl']) == pytest.approx(7.068e-3, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.4462, rel=1e-3)


@pytest.mark.component
def test_SepNF_SepRO_NF_ion_basis_example():
    m = SepNF_SepRO.run_flowsheet_example('ion')
    assert value(m.fs.mixer.pretreatment.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.6512, rel=1e-3)
    assert value(m.fs.mixer.pretreatment.flow_mass_phase_comp[0, 'Liq', 'Na']) == pytest.approx(7.507e-3, rel=1e-3)
    assert value(m.fs.mixer.bypass.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.2412, rel=1e-3)
    assert value(m.fs.mixer.bypass.flow_mass_phase_comp[0, 'Liq', 'Cl']) == pytest.approx(5.079e-3, rel=1e-3)
    assert value(m.fs.RO.retentate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(0.4462, rel=1e-3)
