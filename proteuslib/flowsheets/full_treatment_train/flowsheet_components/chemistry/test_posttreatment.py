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

from proteuslib.flowsheets.full_treatment_train.flowsheet_components.chemistry.posttreatment_ideal_naocl_chlorination_block import (
    run_ideal_naocl_mixer_example, run_ideal_naocl_chlorination_example, run_chlorination_block_example,
    build_translator_from_RO_to_chlorination_block)
from proteuslib.flowsheets.full_treatment_train.model_components import property_models

__author__ = "Austin Ladshaw"

@pytest.mark.component
def test_ideal_naocl_mixer():
    model = run_ideal_naocl_mixer_example()
    assert model.fs.ideal_naocl_mixer_unit.dosing_rate.value == \
            pytest.approx(0.4, rel=1e-3)
    assert model.fs.ideal_naocl_mixer_unit.outlet.flow_mol[0].value == \
            pytest.approx(10.00001074690966, rel=1e-3)
    assert model.fs.ideal_naocl_mixer_unit.outlet.mole_frac_comp[0, "OCl_-"].value == \
            pytest.approx(5.373448801531194e-07, rel=1e-3)

@pytest.mark.component
def test_ideal_naocl_chlorination():
    model = run_ideal_naocl_chlorination_example()
    assert model.fs.ideal_naocl_chlorination_unit.free_chlorine.value == \
            pytest.approx(1.893849296278168, rel=1e-3)
    assert model.fs.ideal_naocl_chlorination_unit.outlet.mole_frac_comp[0,'OCl_-'].value == \
            pytest.approx(4.254858076026844e-07, rel=1e-3)
    assert model.fs.ideal_naocl_chlorination_unit.outlet.mole_frac_comp[0,'H_+'].value == \
            pytest.approx(5.640767226351033e-11, rel=1e-3)

@pytest.mark.component
def test_ideal_naocl_chlorination_full_block():
    model = run_chlorination_block_example(fix_free_chlorine=True)
    assert model.fs.ideal_naocl_mixer_unit.dosing_rate.value == \
            pytest.approx(0.9504457542440085, rel=1e-3)
    assert model.fs.ideal_naocl_mixer_unit.outlet.flow_mol[0].value == \
            pytest.approx(25.000025535888078, rel=1e-3)
    assert model.fs.ideal_naocl_mixer_unit.outlet.mole_frac_comp[0, "OCl_-"].value == \
            pytest.approx(5.107172398859054e-07, rel=1e-3)

    assert model.fs.ideal_naocl_chlorination_unit.free_chlorine.value == \
            pytest.approx(2, rel=1e-3)
    assert model.fs.ideal_naocl_chlorination_unit.outlet.mole_frac_comp[0,'OCl_-'].value == \
            pytest.approx(4.508813505652909e-07, rel=1e-3)
    assert model.fs.ideal_naocl_chlorination_unit.outlet.mole_frac_comp[0,'H_+'].value == \
            pytest.approx(5.479684673084312e-11, rel=1e-3)

@pytest.mark.component
def test_addition_of_translator():
    model = run_ideal_naocl_mixer_example()
    property_models.build_prop(model, base='TDS')
    build_translator_from_RO_to_chlorination_block(model)
    assert hasattr(model.fs, 'RO_to_Chlor')
