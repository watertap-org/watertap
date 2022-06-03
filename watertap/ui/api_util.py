"""
Utility functions for the ``api`` module.
"""
from pathlib import Path
from typing import IO


def open_file_or_stream(fos, attr="tell", **kwargs) -> IO:
    """Open a file or use the existing stream. Avoids adding this logic to every
      function that wants to provide multiple ways of specifying a file.

    Args:
        fos: File or stream
        attr: Attribute to check on the ``fos`` object to see if it is a stream
        kwargs: Additional keywords passed to the ``open`` call. Ignored if the input
          is a stream.

    Returns:
        Opened stream object
    """
    if isinstance(fos, Path):
        output = open(fos, **kwargs)
    elif hasattr(fos, attr):
        output = fos
    else:
        output = open(fos, **kwargs)
    return output
