"""
Model for the data that is created and consumed by the user interface API.
"""
from typing import List, Union, Optional, Dict, Tuple, Type
from pydantic import BaseModel, Extra


class IndexedValue(BaseModel):
    index: List[List[Union[float, str]]]  # do not change order!
    value: List[Union[float, str]]  # do not change order!


class ScalarValue(BaseModel):
    value: Union[float, str]  # do not change order!


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
