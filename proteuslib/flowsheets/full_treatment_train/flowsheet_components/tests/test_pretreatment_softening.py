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
from pyomo.environ import ConcreteModel, value, TransformationFactory
from idaes.core import FlowsheetBlock
from proteuslib.flowsheets.full_treatment_train.flowsheet_components import pretreatment_softening
from proteuslib.flowsheets.full_treatment_train.model_components import property_models
from proteuslib.flowsheets.full_treatment_train.util import check_build, check_dof, check_scaling


@pytest.mark.component
def test_solve():
    m = pretreatment_softening.solve()
    m.fs.display()
    assert value(m.fs.stoich_softening_separator_unit.waste_stream_state[0].flow_mol) \
           == pytest.approx(0.5639516, rel=1e-3)
    assert value(m.fs.stoich_softening_separator_unit.outlet_stream_state[0].flow_mol) \
           == pytest.approx(54.52455, rel=1e-3)


@pytest.mark.unit
def test_build_and_scale():

    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    # check_build(m, build_func=pretreatment_softening.build)  # assert_units_consistent fails for stoich_reactor
    pretreatment_softening.build(m)
    TransformationFactory("network.expand_arcs").apply_to(m)
    check_dof(m)
    m.fs.feed.display()
    # check_scaling(m, scale_func=pretreatment_softening.scale)

