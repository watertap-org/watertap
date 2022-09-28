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
from watertap.examples.flowsheets.case_studies.municipal_treatment.municipal_treatment import (
    main,
)


# -----------------------------------------------------------------------------
@pytest.mark.component
def test_municipal_treatment():
    m, results = main()

    assert_optimal_termination(results)

    assert value(m.fs.feed.properties[0].flow_mass_comp["H2O"]) == pytest.approx(
        921.8, rel=1e-3
    )
    assert value(m.fs.feed.properties[0].flow_mass_comp["tds"]) == pytest.approx(
        0.5811, rel=1e-3
    )
    assert value(m.fs.feed.properties[0].flow_mass_comp["toc"]) == pytest.approx(
        3.690e-3, rel=1e-3
    )
    assert value(m.fs.feed.properties[0].flow_mass_comp["tss"]) == pytest.approx(
        6.019e-3, rel=1e-3
    )

    assert value(
        m.fs.backwash_pump.properties[0].flow_mass_comp["H2O"]
    ) == pytest.approx(36.87, rel=1e-3)
    assert value(
        m.fs.backwash_pump.properties[0].flow_mass_comp["tds"]
    ) == pytest.approx(0, abs=1e-6)
    assert value(
        m.fs.backwash_pump.properties[0].flow_mass_comp["toc"]
    ) == pytest.approx(0, abs=1e-6)
    assert value(
        m.fs.backwash_pump.properties[0].flow_mass_comp["tss"]
    ) == pytest.approx(5.356e-5, rel=1e-3)

    assert value(
        m.fs.recharge_pump.properties[0].flow_mass_comp["H2O"]
    ) == pytest.approx(884.8, rel=1e-3)
    assert value(
        m.fs.recharge_pump.properties[0].flow_mass_comp["tds"]
    ) == pytest.approx(5.811e-2, rel=1e-3)
    assert value(
        m.fs.recharge_pump.properties[0].flow_mass_comp["toc"]
    ) == pytest.approx(9.477e-4, rel=1e-3)
    assert value(
        m.fs.recharge_pump.properties[0].flow_mass_comp["tss"]
    ) == pytest.approx(1.656e-6, rel=1e-3)

    assert value(m.fs.costing.LCOW) == pytest.approx(5.0547e-7, rel=1e-3)  # in M$/m**3
    assert value(m.fs.costing.electricity_intensity) == pytest.approx(
        4.4812e-1, rel=1e-3
    )
