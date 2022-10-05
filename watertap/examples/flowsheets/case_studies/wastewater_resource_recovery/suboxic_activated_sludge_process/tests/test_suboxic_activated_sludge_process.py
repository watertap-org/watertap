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
from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.suboxic_activated_sludge_process.suboxic_activated_sludge_process import (
    main,
)

# -----------------------------------------------------------------------------
@pytest.mark.component
def test_suboxicASM():
    m, results = main()

    assert_optimal_termination(results)

    assert value(m.fs.feed.properties[0].flow_mass_comp["H2O"]) == pytest.approx(
        355.3038, rel=1e-3
    )
    assert value(m.fs.feed.properties[0].flow_mass_comp["bod"]) == pytest.approx(
        0.1116, rel=1e-3
    )
    assert value(m.fs.feed.properties[0].flow_mass_comp["tss"]) == pytest.approx(
        0.1195, rel=1e-3
    )
    assert value(m.fs.feed.properties[0].flow_mass_comp["tkn"]) == pytest.approx(
        0.01849, rel=1e-3
    )
    assert value(m.fs.feed.properties[0].flow_mass_comp["phosphorus"]) == pytest.approx(
        0.002133, rel=1e-3
    )

    assert value(
        m.fs.suboxicASM.properties_treated[0].flow_mass_comp["H2O"]
    ) == pytest.approx(355.3038, rel=1e-3)
    assert value(
        m.fs.suboxicASM.properties_treated[0].flow_mass_comp["bod"]
    ) == pytest.approx(0.003550, rel=1e-3)
    assert value(
        m.fs.suboxicASM.properties_treated[0].flow_mass_comp["tss"]
    ) == pytest.approx(0.001780, rel=1e-3)
    assert value(
        m.fs.suboxicASM.properties_treated[0].flow_mass_comp["tkn"]
    ) == pytest.approx(0.002134, rel=1e-3)
    assert value(
        m.fs.suboxicASM.properties_treated[0].flow_mass_comp["phosphorus"]
    ) == pytest.approx(0.0003556, rel=1e-3)

    assert value(m.fs.costing.LCOW) == pytest.approx(0.660649, rel=1e-3)  # in $/m**3
    assert value(m.fs.costing.electricity_intensity) == pytest.approx(0.5003, rel=1e-3)
