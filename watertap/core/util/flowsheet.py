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

import os
from collections import defaultdict
import pandas as pd

from pyomo.network import Arc, Port
from pyomo.environ import (
    Var,
    Param,
    Expression,
    Objective,
    Block,
    value,
    units as pyunits,
)
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
        arc_name = arc.name
        # Check if the arc's expanded block exists and is deactivated
        arc_expanded = arc.expanded_block
        is_deactivated = arc_expanded is not None and not arc_expanded.active

        source_label = arc.source.name + (" (deactivated)" if is_deactivated else "")
        dest_label = arc.dest.name + (" (deactivated)" if is_deactivated else "")

        # Assigns the destination to source ports
        if arc.source.name not in port_to_arc:
            port_to_arc[arc.source.name] = {"Arc": [], "Source": [], "Destination": []}
        port_to_arc[arc.source.name]["Arc"].append(arc_name)
        port_to_arc[arc.source.name]["Destination"].append(dest_label)

        # Assigns the source to destination ports
        if arc.dest.name not in port_to_arc:
            port_to_arc[arc.dest.name] = {"Arc": [], "Source": [], "Destination": []}
        port_to_arc[arc.dest.name]["Arc"].append(arc_name)
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
                port.name, {"Arc": None, "Source": None, "Destination": None}
            )
            if not connected["Source"] and not connected["Destination"]:
                _log.warning(f"Port {port.name} is not connected to any stream.")
            rows.append(
                {
                    "Unit Model": type(unit).__name__.removeprefix("_Scalar"),
                    "Port Name": name,
                    "Port": port.name,
                    "Source": connected["Source"] or "None",
                    "Destination": connected["Destination"] or "None",
                    "Arc": connected["Arc"] or "None",
                }
            )

    # Display table
    df = pd.DataFrame(rows)
    print(df.to_string(index=False))

    return df


def get_block_data(
    blk,
    descend_into=True,
    components=[Var, Expression, Param],
    sweep_mode=False,
    sweep_cols=dict(),
):
    """
    Get the data from a block for export and reporting.

    Args:
        blk: The block to extract data from.
        descend_into (bool): Whether to descend into sub-blocks. Default is True.
        components (list): List of component types to extract. Default is [Var, Expression, Param].
        sweep_mode (bool): Whether to return only the "value" entry. Default is False.
        sweep_cols (dict): Optional user-provided sweep column names and corresponding Pyomo components.

    Returns:
        dict: A dictionary containing the extracted data where keys are model
            component names and values are the component values
    """

    if not all(c in [Var, Expression, Param, Objective] for c in components):
        raise ValueError(
            "The only accepted components for export are Var, Expression, Param, and Objective."
        )

    data = defaultdict(dict)

    for c in blk.component_objects(components, descend_into=descend_into):
        if c.is_reference():
            # The object is a Reference and we will get to the referenced
            # object eventually. Notably, if the Reference referent is
            # on a sub-block we will still descend_into it even if descend_into=False.
            continue
        if c.is_indexed():
            for ci in c.values():
                data[ci.name]["component_type"] = (
                    type(ci).__name__.removeprefix("Scalar").removesuffix("Data")
                )
                data[ci.name]["value"] = value(ci)
                data[ci.name]["units"] = "/".join(
                    "year" if x == "a" else x
                    for x in pyunits.get_units(ci).getname().split("/")
                )
        else:
            data[c.name]["component_type"] = (
                type(c).__name__.removeprefix("Scalar").removesuffix("Data")
            )
            data[c.name]["value"] = value(c)
            data[c.name]["units"] = "/".join(
                "year" if x == "a" else x
                for x in pyunits.get_units(c).getname().split("/")
            )

    if sweep_mode:
        # Only return "value" entry
        data = {k: v["value"] for k, v in data.items()}
        # Add user-provided sweep column names for convenience
        for k, v in sweep_cols.items():
            data[k] = value(v)

    return data


def block_data_to_df(blk_data):
    """
    Convert block data dictionary to a pandas DataFrame.

    Args:
        blk_data (dict): The block data dictionary obtained from `get_block_data`.

    Returns:
        df: A DataFrame containing the block data with columns
            for model component, value, units, and component type.
    """

    df = pd.DataFrame(blk_data).T
    df["model_component"] = df.index
    df.reset_index(inplace=True, drop=True)
    if df.empty:
        raise ValueError("Model export failed: no data to export.")
    df = df[["model_component", "value", "units", "component_type"]]

    return df


def export_block_data_to_csv(
    blk,
    save_as=None,
    **kwargs,
):
    """
    Export block data to a CSV file.

    Args:
        blk: The Pyomo block containing the model components.
        save_as (str, optional): The file path to save the csv file.
            Defaults to "cwd/watertap_model_results.csv".
        **kwargs: Additional keyword arguments passed to `get_block_data`.

    Returns:
        pd.DataFrame: A DataFrame containing the exported block data.
    """

    if save_as is None:
        save_as = f"{os.getcwd()}/watertap_model_results"
    save_as = save_as.replace(".csv", "")

    components = kwargs.get("components", [Var, Expression, Param])
    descend_into = kwargs.get("descend_into", True)

    blk_data = get_block_data(blk, components=components, descend_into=descend_into)
    blk_df = block_data_to_df(blk_data)

    blk_df.to_csv(f"{save_as}.csv", index=False)
    _log.info(f"{blk.name} data exported to {save_as}.csv")

    return blk_df
