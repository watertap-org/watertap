#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import pandas as pd
from pyomo.network import Arc, Port
from pyomo.environ import Block
from idaes.core import UnitModelBlockData, FlowsheetBlockData

import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


def list_ports(block, descend_into=False):
    """
    Lists all the inlet and outlet ports on a unit model or flowsheet
    Args:
        block : a unit model or flowsheet model block
        descend_into : whether or not to consider arcs in nested flowsheets or sub-blocks

    Returns:
        df : DataFrame with the unit model name, port name, and port path
    """
    if not isinstance(block, (UnitModelBlockData, FlowsheetBlockData)):
        raise TypeError(
            f"Expected a UnitModelBlockData or FlowsheetBlockData instance, but "
            f"got {type(block).__name__!r}."
        )

    # Finds the flowsheet whether a flowsheet or unit model was passed
    if isinstance(block, FlowsheetBlockData):
        flowsheet = block
    else:
        flowsheet = block.parent_block()

    port_to_arc = {}
    for arc in flowsheet.component_objects(Arc, descend_into=descend_into):
        # Check if the arc's expanded block exists and is deactivated
        arc_expanded = arc.expanded_block
        is_deactivated = arc_expanded is not None and not arc_expanded.active

        source_label = arc.source.name + (" (deactivated)" if is_deactivated else "")
        dest_label = arc.dest.name + (" (deactivated)" if is_deactivated else "")

        # Assigns the destination to source ports
        if arc.source.name not in port_to_arc:
            port_to_arc[arc.source.name] = {"Source": [], "Destination": []}
        port_to_arc[arc.source.name]["Destination"].append(dest_label)

        # Assigns the source to destination ports
        if arc.dest.name not in port_to_arc:
            port_to_arc[arc.dest.name] = {"Source": [], "Destination": []}
        port_to_arc[arc.dest.name]["Source"].append(source_label)

    rows = []

    # If a flowsheet was passed, collect all unit models
    if isinstance(block, FlowsheetBlockData):
        units = [
            u
            for u in block.component_objects(Block)
            if isinstance(u, UnitModelBlockData)
        ]
    else:
        units = [block]

    # For each unit, identify its name and all the port information
    for unit in units:
        ports = dict(unit.component_map(Port))
        for name, port in ports.items():
            connected = port_to_arc.get(
                port.name, {"Source": None, "Destination": None}
            )
            if not connected["Source"] and not connected["Destination"]:
                _log.warning(f"Port {port.name} is not connected to any stream.")
            rows.append(
                {
                    "Unit Model": type(unit).__name__.removeprefix("_Scalar"),
                    "Port Name": name,
                    "Port": port.name,
                    "Source": connected["Source"] or None,
                    "Destination": connected["Destination"] or None,
                }
            )

    # Display table
    df = pd.DataFrame(rows)
    print(df.to_string(index=False))

    return df
