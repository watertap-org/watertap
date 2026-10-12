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

"""Helpers for selecting representative LSRRO baseline cases."""


def get_test_cases_subset(all_test_cases):
    """Return one representative case per baseline file.

    The test file uses these as the pull-request subset and sends the remaining
    combinations to the nightly-only bucket.
    """

    seen_files = set()
    subset = []
    for csv_file, row_index in all_test_cases:
        if csv_file in seen_files:
            continue
        subset.append((csv_file, row_index))
        seen_files.add(csv_file)
    return subset
