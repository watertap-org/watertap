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
This module contains utility functions for WaterTAP reports.
"""

import difflib


def print_report_differences(report_1, report_2):
    """
    Utility to help determine the differences between two reports and provide instructions
    on what changes line by line will make report_1 equivalent to report_2.

    Arguments:
        report_1 - the first report
        report_2 - the second report
    """
    if report_1 != report_2:
        diff_obj = difflib.Differ()
        result = diff_obj.compare(report_1.split('\n'), report_2.split('\n'))
        print('\n'.join(result))
