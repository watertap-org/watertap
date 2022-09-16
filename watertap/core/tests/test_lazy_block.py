###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################
"""
Tests for LazyBlock
"""

import pyomo.environ as pyo
import idaes.core as idc

from watertap.core.lazy_block import LazyBlock, LazyBlockMixin


@idc.declare_process_block_class("TestProcessBlock")
class TestProcessBlockData(LazyBlockMixin, idc.ProcessBlockData):
    def build(self):
        def build_scalar_lazy_block(b):
            b.x = pyo.Var()
            b.y = pyo.Var()
            b.c = pyo.Constraint(rule=b.x + b.y <= 1)

        self.scalar_lazy_block = LazyBlock(rule=build_scalar_lazy_block)

        def build_indexed_lazy_block(b, idx):
            b.x = pyo.Var()
            b.y = pyo.Var(list(range(idx)))
            b.c = pyo.Constraint(rule=b.x + sum(b.y.values()) <= 1)

        self.indexed_lazy_block = LazyBlock([1, 2], rule=build_indexed_lazy_block)


def test_lazy_block():
    m = pyo.ConcreteModel()
    m.fs = idc.FlowsheetBlock()
    m.fs.p = TestProcessBlock()

    assert len(m.fs.p._lazy_blocks) == 2

    assert len(list(m.component_data_objects((pyo.Var, pyo.Constraint)))) == 0

    # accessing adds the block, and happens
    # seamlessly
    m.fs.p.scalar_lazy_block.x.value = 2
    m.fs.p.scalar_lazy_block.y.value = 1

    assert len(m.fs.p._lazy_blocks) == 1

    assert len(list(m.component_data_objects((pyo.Var, pyo.Constraint)))) == 3

    # indexing works fine, and
    # constructs the whole set
    m.fs.p.indexed_lazy_block[2].y[1].value = 5
    m.fs.p.indexed_lazy_block[2].y[0].value = 3

    assert len(m.fs.p._lazy_blocks) == 0

    assert len(list(m.component_data_objects((pyo.Var, pyo.Constraint)))) == 3 + 3 + 4
