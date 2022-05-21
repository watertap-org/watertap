"""
Utility functions for the ``api`` module.
"""
# standard library
import inspect
import json
import functools
from string import Template
from typing import Dict, Union, Optional, IO

#: Set this logger from the api module
util_logger = None


def log_meth(meth):
    @functools.wraps(meth)
    def wrapper(*args, **kwargs):
        name = _get_method_classname(meth)
        util_logger.debug(f"Begin {name}")
        try:
            result = meth(*args, **kwargs)
        except Exception as e:
            error_msg = f"Error in {name}"
            util_logger.exception(error_msg)
            raise e
        util_logger.debug(f"End {name}")
        return result

    return wrapper


def _get_method_classname(m):
    """Get class name for method, assuming method is bound and class has __dict__."""
    for k, v in inspect.getmembers(m):
        if k == "__qualname__":
            return v
    return "<unknown>"


# End logging


def open_file_or_stream(fos, attr, **kwargs) -> IO:
    """Open a file or use the existing stream. Avoids adding this logic to every
      function that wants to provide multiple ways of specifying a file.

    Args:
        fos: File or stream
        attr: Attribute to check on the ``fos`` object to see if it is a stream,
          e.g. "write" or "read"
        kwargs: Additional keywords passed to the ``open`` call. Ignored if the input
          is a stream.

    Returns:
        Opened stream object
    """
    if hasattr(fos, attr):
        output = fos
    else:
        output = open(fos, **kwargs)
    return output

