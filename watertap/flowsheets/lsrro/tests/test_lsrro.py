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
import pytest
from pyomo.environ import assert_optimal_termination, value
from watertap.flowsheets.lsrro import lsrro


@pytest.mark.component
def test_lsrro():
    m, results = lsrro.main()
    assert_optimal_termination(results)
    assert pytest.approx(value(m.fs.costing.LCOW), rel=1e-3) == 1.4660
    assert (
        pytest.approx(value(m.fs.costing.specific_energy_consumption), rel=1e-3)
        == 7.77891
    )
    assert (
        pytest.approx(value(sum(m.fs.ROUnits[a].area for a in m.fs.Stages)), rel=1e-3)
        == 195.0059
    )
