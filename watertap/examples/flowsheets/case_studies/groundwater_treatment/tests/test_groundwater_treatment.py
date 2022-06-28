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

from pyomo.util.check_units import assert_units_consistent
from pyomo.environ import (
    value,
    assert_optimal_termination,
)

from watertap.core.util.initialization import assert_degrees_of_freedom

from watertap.examples.flowsheets.case_studies.groundwater_treatment.groundwater_treatment import (
    build,
    set_operating_conditions,
    initialize_system,
    solve,
    add_costing,
    display_costing,
    display_results,
)


class TestGroundwaterTreatmentFlowsheet:
    @pytest.fixture(scope="class")
    def system_frame(self):
        m = build()
        return m

    @pytest.mark.unit
    def test_build(self, system_frame):
        m = system_frame
        assert_degrees_of_freedom(m, 21)
        assert_units_consistent(m)


# TODO write more tests
