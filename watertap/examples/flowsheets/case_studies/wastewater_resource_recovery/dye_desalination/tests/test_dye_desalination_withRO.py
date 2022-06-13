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
from pyomo.environ import (
    value,
    assert_optimal_termination,
)
from idaes.core.util import get_solver
from watertap.core.util.initialization import assert_degrees_of_freedom
from pyomo.util.check_units import assert_units_consistent
from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.dye_desalination_withRO.dye_desalination_withRO import (
    build,
    set_operating_conditions,
    initialize_system,
    initialize_costing,
    solve,
    add_costing,
    display_results,
    display_costing,
)

solver = get_solver()
