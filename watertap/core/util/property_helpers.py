#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import pandas as pd
from idaes.core.util.config import is_physical_parameter_block


def get_property_metadata(prop_pkg):
    """Get all supported properties from a WaterTAP/IDAES property package as a Pandas DataFrame."""
    try:
        assert is_physical_parameter_block(prop_pkg)
    except:
        raise TypeError("get_property_metadata expected a PhysicalParameterBlock.")
    metadata = prop_pkg.get_metadata()
    pd.set_option("display.max_rows", None)
    df = pd.DataFrame(
        {
            "Description": [v._doc for v in metadata.properties],
            "Name": [v._name for v in metadata.properties],
            "Units": [str(v._units) for v in metadata.properties],
        }
    )
    return df
