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

import pytest
from pyomo.environ import value, assert_optimal_termination
from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.electrochemical_nutrient_removal.electrochemical_nutrient_removal import (
    main,
)


@pytest.mark.component
def test_treatment():
    m, results = main()
    assert_optimal_termination(results)
