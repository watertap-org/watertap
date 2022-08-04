"""
Simple flowsheet interface API
"""
from typing import Callable
from pydantic import BaseModel


class Action(BaseModel):
    """An action to take."""

    name: str
    func: Callable


class ModelObject(BaseModel):
    """A variable, expression, or parameter."""

    name: str
    id: str
    is_input: bool
    is_output: bool
    # TODO: more fields
