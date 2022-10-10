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

from pyomo.environ import ConcreteModel, assert_optimal_termination, value
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.core.util.scaling as iscale

from watertap.examples.custom_model_demo.demo_simple_filter import main

solver = get_solver()

# -----------------------------------------------------------------------------
@pytest.mark.component
def test_demo_simple_prop_pack():
    m = main()

    assert pytest.approx(0.97, rel=1e-3) == value(
        m.fs.filter.properties_out[0].flow_mass_phase_comp["Liq", "H2O"]
    )
    assert pytest.approx(3.395e-2, rel=1e-3) == value(
        m.fs.filter.properties_out[0].flow_mass_phase_comp["Liq", "NaCl"]
    )
    assert pytest.approx(1.200e-5, rel=1e-3) == value(
        m.fs.filter.properties_out[0].flow_mass_phase_comp["Liq", "TSS"]
    )
