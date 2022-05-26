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
helper for splitting up the LSRRO testing jobs between the GHA runners
"""

import sys
import platform

_supported_systems = ["Linux", "Windows"]
_number_of_python_versions = 3
_python_version_index = sys.version_info.minor % _number_of_python_versions

_this_platform = platform.system()
try:
    _platform_index = _supported_systems.index(_this_platform)
except ValueError:
    _platform_index = None


def get_test_cases_subset(all_test_cases):
    return [
        test_case
        for idx, test_case in enumerate(all_test_cases)
        if (idx % len(_supported_systems) == _platform_index)
        and (idx % _number_of_python_versions == _python_version_index)
    ]
