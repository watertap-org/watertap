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
from pyomo.network import Port
from pyomo.environ import Block
from idaes.core import UnitModelBlockData, FlowsheetBlockData


def list_ports(block):
    """
    Lists all the inlet and outlet ports on a unit model or flowsheet
    Args:
        block : a unit model or flowsheet model block

    Returns:
        df : DataFrame with the unit model name, port name, and port path
    """
    if not isinstance(block, (UnitModelBlockData, FlowsheetBlockData)):
        raise TypeError(
            f"Expected a UnitModelBlockData or FlowsheetBlockData instance, but "
            f"got {type(block).__name__!r}."
        )

    rows = []

    if isinstance(block, FlowsheetBlockData):
        units = [
            u
            for u in block.component_objects(Block)
            if isinstance(u, UnitModelBlockData)
        ]
    else:
        units = [block]

    for unit in units:
        ports = dict(unit.component_map(Port))
        for name, port in ports.items():
            rows.append(
                {
                    "Unit Model": type(unit).__name__.removeprefix("_Scalar"),
                    "Port Name": name,
                    "Port": port.name,
                }
            )

    df = pd.DataFrame(rows)

    pd.set_option("display.max_rows", None)
    print(df.to_string(index=False))

    return df
