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

"""
    Simple unit tests for example flowsheet of SepRO with Chlorination.

    NOTE: That flowsheet is not meant to be viewed as a final product, but
    a sample of how to incorporate more complex chemistry into a simple
    flowsheet.
"""
import pytest

from proteuslib.flowsheets.full_treatment_train.flowsheet_components.chemistry.PreTreatment_Simple_Softening import (
    run_softening_example)

from proteuslib.flowsheets.full_treatment_train.flowsheet_components.chemistry.PostTreatment_SimpleNaOCl_Chlorination import (
    run_chlorination_example, run_chlorination_constrained_outlet_example)

from proteuslib.flowsheets.full_treatment_train.flowsheet_components.chemistry.SepRO_plus_Chlorination import (
    run_SepRO_Chlorination_flowsheet_example, run_SepRO_Chlorination_flowsheet_with_outlet_constraint_example)

from proteuslib.flowsheets.full_treatment_train.flowsheet_components.chemistry.ZeroDRO_plus_Chlorination import (
    run_0DRO_Chlorination_flowsheet_example, run_0DRO_Chlorination_flowsheet_optimization_example)

__author__ = "Austin Ladshaw"

@pytest.mark.component
def test_Simple_Softening():
    model = run_softening_example()
    assert model.fs.softening_unit.outlet.mole_frac_comp[0,'CaCO3'].value == \
            pytest.approx(5.99873e-05, rel=1e-3)
    assert model.fs.softening_unit.outlet.flow_mol[0].value == \
            pytest.approx(10.0, rel=1e-2)

@pytest.mark.component
def test_SimpleNaOCl_Chlorination():
    model = run_chlorination_example()
    assert model.fs.simple_naocl_unit.outlet.mole_frac_comp[0,'Cl_-'].value == \
            pytest.approx(0.0, rel=1e-3)
    assert model.fs.simple_naocl_unit.outlet.mole_frac_comp[0,'H2O'].value == \
            pytest.approx(0.9999989745408085, rel=1e-3)
    assert model.fs.simple_naocl_unit.outlet.mole_frac_comp[0,'HOCl'].value == \
            pytest.approx(5.8124920339055656e-08, rel=1e-3)
    assert model.fs.simple_naocl_unit.outlet.mole_frac_comp[0,'H_+'].value == \
            pytest.approx(5.640767226351033e-11, rel=1e-3)
    assert model.fs.simple_naocl_unit.outlet.mole_frac_comp[0,'Na_+'].value == \
            pytest.approx(4.836107279417401e-07, rel=1e-3)
    assert model.fs.simple_naocl_unit.outlet.mole_frac_comp[0,'OCl_-'].value == \
            pytest.approx(4.254858076026844e-07, rel=1e-3)
    assert model.fs.simple_naocl_unit.outlet.mole_frac_comp[0,'OH_-'].value == \
            pytest.approx(5.8181328011319155e-08, rel=1e-3)

@pytest.mark.component
def test_SimpleNaOCl_Chlorination_with_outlet_constraint():
    model = run_chlorination_constrained_outlet_example()
    assert model.fs.simple_naocl_unit.free_chlorine.value == \
            pytest.approx(2.0, rel=1e-3)
    assert model.fs.simple_naocl_unit.dosing_rate.value == \
            pytest.approx(0.38017796967633016, rel=1e-3)

@pytest.mark.component
def test_SepRO_Chlorination_flowsheet_example():
    model = run_SepRO_Chlorination_flowsheet_example()
    assert model.fs.RO.permeate.flow_mass_phase_comp[0, 'Liq', 'H2O'].value == \
            pytest.approx(0.4825, rel=1e-3)
    assert model.fs.RO.permeate.flow_mass_phase_comp[0, 'Liq', 'TDS'].value == \
            pytest.approx(0.00035000000000000005, rel=1e-3)

    assert model.fs.RO_to_Chlor.outlet.flow_mol[0].value == \
            pytest.approx(26.81154870624049, rel=1e-3)
    assert model.fs.RO_to_Chlor.outlet.pressure[0].value == \
            pytest.approx(101324.99999999999, rel=1e-3)
    assert model.fs.RO_to_Chlor.outlet.temperature[0].value == \
            pytest.approx(298.15, rel=1e-3)

    assert model.fs.RO_to_Chlor.outlet.mole_frac_comp[0 , 'Cl_-'].value == \
            pytest.approx(0.00022352870214977842, rel=1e-3)
    assert model.fs.RO_to_Chlor.outlet.mole_frac_comp[0 , 'H2O'].value == \
            pytest.approx(0.9995519753742446, rel=1e-3)
    assert model.fs.RO_to_Chlor.outlet.mole_frac_comp[0 , 'Na_+'].value == \
            pytest.approx(0.00022401231287774744, rel=1e-3)
    assert model.fs.RO_to_Chlor.outlet.mole_frac_comp[0 , 'OCl_-'].value == \
            pytest.approx(4.836107279690193e-07, rel=1e-3)

    assert model.fs.simple_naocl_unit.outlet.mole_frac_comp[0,'Cl_-'].value == \
            pytest.approx(0.00022352870213716955, rel=1e-3)
    assert model.fs.simple_naocl_unit.outlet.mole_frac_comp[0,'H2O'].value == \
            pytest.approx(0.9995519171365305, rel=1e-3)
    assert model.fs.simple_naocl_unit.outlet.mole_frac_comp[0,'HOCl'].value == \
            pytest.approx(5.8124923290629284e-08, rel=1e-3)
    assert model.fs.simple_naocl_unit.outlet.mole_frac_comp[0,'H_+'].value == \
            pytest.approx(5.640806579383571e-11, rel=1e-3)
    assert model.fs.simple_naocl_unit.outlet.mole_frac_comp[0,'Na_+'].value == \
            pytest.approx(0.00022401231286511128, rel=1e-3)
    assert model.fs.simple_naocl_unit.outlet.mole_frac_comp[0,'OCl_-'].value == \
            pytest.approx(4.2548580465111037e-07, rel=1e-3)
    assert model.fs.simple_naocl_unit.outlet.mole_frac_comp[0,'OH_-'].value == \
            pytest.approx(5.8181331356423126e-08, rel=1e-3)

@pytest.mark.component
def test_SepRO_Chlorination_flowsheet_with_outlet_constraint():
    model = run_SepRO_Chlorination_flowsheet_with_outlet_constraint_example()
    assert model.fs.simple_naocl_unit.free_chlorine.value == \
            pytest.approx(2.0, rel=1e-3)
    assert model.fs.simple_naocl_unit.dosing_rate.value == \
            pytest.approx(1.0196136841188141, rel=1e-3)


@pytest.mark.component
def test_0DRO_Chlorination_flowsheet_with_outlet_constraint_with_seq_decomp():
    model = run_0DRO_Chlorination_flowsheet_example(True)
    assert model.fs.simple_naocl_unit.free_chlorine.value == \
            pytest.approx(2.0, rel=1e-3)
    assert model.fs.simple_naocl_unit.dosing_rate.value == \
            pytest.approx(0.5604052608884579, rel=1e-3)

@pytest.mark.component
def test_0DRO_Chlorination_flowsheet_with_outlet_constraint_without_seq_decomp():
    model = run_0DRO_Chlorination_flowsheet_example(False)
    assert model.fs.simple_naocl_unit.free_chlorine.value == \
            pytest.approx(2.0, rel=1e-3)
    assert model.fs.simple_naocl_unit.dosing_rate.value == \
            pytest.approx(0.5604052608884579, rel=1e-3)


@pytest.mark.component
def test_0DRO_Chlorination_flowsheet_with_simple_optimization():
    model = run_0DRO_Chlorination_flowsheet_optimization_example()
    assert model.fs.simple_naocl_unit.free_chlorine.value == \
            pytest.approx(2.000000454910428, rel=1e-3)
    assert model.fs.simple_naocl_unit.dosing_rate.value == \
            pytest.approx(0.4211328239720586, rel=1e-3)
