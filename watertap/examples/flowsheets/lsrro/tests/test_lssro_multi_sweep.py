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

# TODO: Load and spot-check parameter sweep analysis corresponding to
#       s*_stage/results_LSRRO.csv's for stages 1-8

# TODO: Maybe we want to randomly generate which of the various test
#       scenarios we run, so we get coverage in expectation.

import os

from watertap.examples.flowsheets.lsrro.multi_sweep import run_case

_this_file_path = os.path.dirname(os.path.abspath(__file__))
