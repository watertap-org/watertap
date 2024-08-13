#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################


import pytest
import itertools

import pyomo.environ as pyo
from watertap.core.membrane_channel0d import _0DPropertyHelper


class Test0DPropertyHelper:

    @pytest.fixture(scope="class")
    def model(self):
        m = pyo.ConcreteModel()
        m.time_periods = pyo.RangeSet(0, 10)
        m.properties_in = pyo.Var(m.time_periods)
        m.properties_out = pyo.Var(m.time_periods)

        m.properties_forward = _0DPropertyHelper(m, reverse=False)
        m.properties_backward = _0DPropertyHelper(m, reverse=True)
        return m

    @pytest.mark.unit
    def test_len(self, model):
        assert len(model.properties_forward) == 2 * len(model.time_periods)
        assert len(model.properties_backward) == 2 * len(model.time_periods)

    @pytest.mark.unit
    def test_iter(self, model):
        found_keys = []
        for k in model.properties_forward:
            found_keys.append(k)
        assert found_keys == [(t, x) for x in (0, 1) for t in model.time_periods]

        found_keys = []
        for k in model.properties_backward:
            found_keys.append(k)
        assert found_keys == [(t, x) for x in (0, 1) for t in model.time_periods]

    @pytest.mark.unit
    def test_contains(self, model):
        assert (11, 0) not in model.properties_forward
        assert (11, 1) not in model.properties_forward
        assert (11, 0.5) not in model.properties_forward
        assert (5, 0.5) not in model.properties_forward

        assert (5, 0) in model.properties_forward
        assert (5, 1) in model.properties_forward

    @pytest.mark.unit
    def test_getitem(self, model):
        assert model.properties_forward[5, 0] is model.properties_in[5]
        assert model.properties_forward[5, 1] is model.properties_out[5]

        assert model.properties_backward[5, 0] is model.properties_out[5]
        assert model.properties_backward[5, 1] is model.properties_in[5]

        for v1, v2 in zip(model.properties_forward[:, 0], model.properties_in[:]):
            assert v1 is v2
        for v1, v2 in zip(model.properties_backward[:, 0], model.properties_out[:]):
            assert v1 is v2
        for v1, v2 in zip(
            model.properties_forward[5, :],
            (model.properties_in[5], model.properties_out[5]),
        ):
            assert v1 is v2
        for v1, v2 in zip(
            model.properties_backward[5, :],
            (model.properties_out[5], model.properties_in[5]),
        ):
            assert v1 is v2

        for v1, v2 in zip(
            model.properties_backward[:, :],
            itertools.chain(model.properties_out[:], model.properties_in[:]),
        ):
            assert v1 is v2

        for v1, v2 in zip(
            model.properties_forward[...],
            itertools.chain(model.properties_in[:], model.properties_out[:]),
        ):
            assert v1 is v2

        with pytest.raises(KeyError):
            model.properties_forward[11, 0]
        with pytest.raises(KeyError):
            model.properties_forward[11, 1]
        with pytest.raises(KeyError):
            model.properties_forward[11, 0.5]
        with pytest.raises(KeyError):
            model.properties_forward[5, 0.5]
        with pytest.raises(KeyError):
            model.properties_forward[5, 0, 2]
        with pytest.raises(KeyError):
            model.properties_forward[5]
        with pytest.raises(IndexError):
            model.properties_forward[slice(0, 10, 2), 0]

    @pytest.mark.unit
    def test_keys(self, model):
        for k1, k2 in zip(
            model.properties_forward.keys(),
            itertools.chain(
                zip(model.properties_in.keys(), [0] * len(model.time_periods)),
                zip(model.properties_out.keys(), [1] * len(model.time_periods)),
            ),
        ):
            assert k1 == k2
        for k1, k2 in zip(
            model.properties_backward.keys(),
            itertools.chain(
                zip(model.properties_out.keys(), [0] * len(model.time_periods)),
                zip(model.properties_in.keys(), [1] * len(model.time_periods)),
            ),
        ):
            assert k1 == k2

    @pytest.mark.unit
    def test_items(self, model):
        for i1, i2 in zip(
            model.properties_forward.items(),
            itertools.chain(
                zip(model.properties_in.items(), [0] * len(model.time_periods)),
                zip(model.properties_out.items(), [1] * len(model.time_periods)),
            ),
        ):
            i2k = (i2[0][0], i2[1])
            i2v = i2[0][1]

            i1k = (i1[0], i1[1])
            i1v = i1[2]

            assert i1k == i2k
            assert i1v is i2v

        for i1, i2 in zip(
            model.properties_backward.items(),
            itertools.chain(
                zip(model.properties_out.items(), [0] * len(model.time_periods)),
                zip(model.properties_in.items(), [1] * len(model.time_periods)),
            ),
        ):
            i2k = (i2[0][0], i2[1])
            i2v = i2[0][1]

            i1k = (i1[0], i1[1])
            i1v = i1[2]

            assert i1k == i2k
            assert i1v is i2v

    @pytest.mark.unit
    def test_values(self, model):
        for v1, v2 in zip(
            model.properties_forward.values(),
            itertools.chain(
                model.properties_in.values(), model.properties_out.values()
            ),
        ):
            assert v1 is v2
        for v1, v2 in zip(
            model.properties_backward.values(),
            itertools.chain(
                model.properties_out.values(), model.properties_in.values()
            ),
        ):
            assert v1 is v2
