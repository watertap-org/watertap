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

from proteuslib.flowsheets.full_treatment_train.chemistry_flowsheets.pretreatment_stoich_softening_block import (
    run_stoich_softening_mixer_example)
from proteuslib.flowsheets.full_treatment_train.example_models import property_models

__author__ = "Austin Ladshaw"

@pytest.mark.component
def test_stoich_softening_mixer():
    model = run_stoich_softening_mixer_example()
    assert model.fs.stoich_softening_mixer_unit.dosing_rate.value == \
            pytest.approx(0.9953385276732939, rel=1e-3)
    assert model.fs.stoich_softening_mixer_unit.outlet.flow_mol[0].value == \
            pytest.approx(10.000013433637829, rel=1e-3)
    assert model.fs.stoich_softening_mixer_unit.outlet.mole_frac_comp[0, "NaCl"].value == \
            pytest.approx(0.01072288805932624, rel=1e-3)
    assert model.fs.stoich_softening_mixer_unit.outlet.mole_frac_comp[0, "Ca(OH)2"].value == \
            pytest.approx(1.3433619782885838e-06, rel=1e-3)
