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
Model for the data that is created and consumed by the user interface API.
"""
import logging
from typing import List, Union, Optional, Dict, Tuple, Type
from pydantic import BaseModel, Extra

_log = logging.getLogger(__name__)


class IndexedValue(BaseModel):
    index: List[List[Union[float, str]]]  # do not change order!
    value: List[Union[float, str]]  # do not change order!
    bounds: List[Tuple[Optional[float], Optional[float]]]


class ScalarValue(BaseModel):
    value: Union[float, str]  # do not change order!
    bounds: Tuple[Optional[float], Optional[float]]


class Variable(BaseModel):
    display_name = ""
    description = ""
    units = ""
    readonly = False
    value: Optional[Union[IndexedValue, ScalarValue]]


class BlockMeta(BaseModel):
    parameters: Dict = {}

    class Config:
        arbitrary_types_allowed = True
        extra = Extra.allow


class Block(BaseModel):
    display_name = ""
    description = ""
    category = "default"
    variables: Dict[str, Variable] = {}
    blocks: Dict[str, "Block"] = {}
    meta: BlockMeta = BlockMeta()

    def get_sole_subblock(self):
        keys = list(self.blocks.keys())
        if len(keys) != 1:
            adj = "many" if len(keys) > 1 else "few"
            raise ValueError(f"Too {adj} sub-blocks. expected=1 got={len(keys)}")
        return keys[0]
