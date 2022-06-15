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
Utility functions for the ``api`` module.
"""
from operator import itemgetter
from pathlib import Path
from typing import IO, Dict, List, Tuple, Union


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


def flatten_tree(
    tree: Dict, tuple_keys: bool = False, copy_value: bool = True, sort: bool = True
) -> List[Tuple[Union[List[str], str], Dict]]:
    """Flatten a tree of blocks.

    Args:
        tree: The tree of blocks. See :mod:`watertap.ui.api` for details on the
          format. Should start like: ``{ "blocks": { "Flowsheet": { ... } } }``
        tuple_keys: Controls whether the flattened keys should be a tuple of
         strings or a single dotted string.
        copy_value: If True, make a copy of the value and remove the "blocks"
          from it. Otherwise, the values are references to the input.
        sort: If True, sort the result by keys before returning it.

    Returns:
        List of tuples of (key, value), where key is the full path to the block and
        the value is the value of the block.
    """
    flattened = []

    def flatten_subtree(path, subtree):
        blocks = subtree["blocks"]
        for name, val in blocks.items():
            if tuple_keys:
                full_name = tuple(list(path) + [name])
            else:
                full_name = f"{path}.{name}"
            if copy_value:
                item = val.copy()
                del item["blocks"]
            else:
                item = val
            flattened.append((full_name, item))
            flatten_subtree(full_name, val)

    # start with first (only) block at first level
    root_blocks = tree["blocks"]
    root_key = list(root_blocks.keys())[0]
    root_block = root_blocks[root_key]
    if tuple_keys:
        flatten_subtree((root_key,), root_block)
    else:
        flatten_subtree(root_key, root_block)

    # sort by full path to item
    if sort:
        flattened.sort(key=itemgetter(0))

    return flattened
