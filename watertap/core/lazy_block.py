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

from pyomo.core.base.block import Block, _BlockData
from pyomo.core.base.global_set import UnindexedComponent_index


__all__ = ["LazyBlockMixin", "LazyBlock"]


class _LazyBlockData(_BlockData):
    pass


class LazyBlock(Block):
    """
    A LazyBlock component, to be used with
    a class which inherits from LazyBlockMixin

    The construction of the LazyBlock will be
    delayed until the first time the block
    is accessed.

    If indexed, all indices are constructed at once.
    """

    _ComponentDataClass = _LazyBlockData

    def __new__(cls, *args, **kwds):
        if cls != LazyBlock:
            return super(LazyBlock, cls).__new__(cls)
        if args == ():
            return ScalarLazyBlock.__new__(ScalarLazyBlock)
        else:
            return IndexedLazyBlock.__new__(IndexedLazyBlock)


class ScalarLazyBlock(_LazyBlockData, LazyBlock):
    def __init__(self, *args, **kwds):
        _LazyBlockData.__init__(self, self)
        LazyBlock.__init__(self, *args, **kwds)
        self._data[None] = self
        self._index = UnindexedComponent_index


class IndexedLazyBlock(LazyBlock):
    pass


class LazyBlockMixin:
    """
    Enables the use of LazyBlocks on a
    parent block via inheriting from this
    class

    example:

        .. code-block::

            class MyClassData(LazyBlockMixin, ProcessBlockData):

                def build(self):

                    def my_lazy_block_rule(blk):
                        blk.x = Var()
                    self.lazy_block = LazyBlock(rule=my_lazy_block_rule)

    """

    def __init__(self, *args, **kwargs):
        self._lazy_blocks = {}
        super().__init__(*args, **kwargs)

    def add_component(self, name, val):
        if isinstance(val, LazyBlock):
            self._lazy_blocks[name] = val
        else:
            super().add_component(name, val)

    def __getattr__(self, name):
        val = self._lazy_blocks.pop(name, None)
        if val is None:
            return super().__getattr__(name)
        else:
            super().add_component(name, val)
            return val
