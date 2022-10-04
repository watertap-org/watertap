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
from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.biomembrane_filtration.biomembrane_filtration import (
    main,
)

# -----------------------------------------------------------------------------
@pytest.mark.component
def test_biomembrane_filtration():
    m, results = main()

    assert_optimal_termination(results)

    assert value(m.fs.feed.properties[0].flow_mass_comp["H2O"]) == pytest.approx(
        115.807, rel=1e-3
    )
    assert value(m.fs.feed.properties[0].flow_mass_comp["bod"]) == pytest.approx(
        0.01193, rel=1e-3
    )
    assert value(
        m.fs.feed.properties[0].flow_mass_comp["ammonium_as_nitrogen"]
    ) == pytest.approx(3.01167e-3, rel=1e-3)
    assert value(m.fs.feed.properties[0].flow_mass_comp["nitrate"]) == pytest.approx(
        1.50583e-4, rel=1e-3
    )
    assert value(m.fs.feed.properties[0].flow_mass_comp["tss"]) == pytest.approx(
        0.011583, rel=1e-3
    )

    assert value(
        m.fs.dmbr.properties_treated[0].flow_mass_comp["H2O"]
    ) == pytest.approx(110.016, rel=1e-3)
    assert value(
        m.fs.dmbr.properties_treated[0].flow_mass_comp["bod"]
    ) == pytest.approx(1.193e-3, rel=1e-3)
    assert value(
        m.fs.dmbr.properties_treated[0].flow_mass_comp["ammonium_as_nitrogen"]
    ) == pytest.approx(9.035e-4, rel=1e-3)
    assert value(
        m.fs.dmbr.properties_treated[0].flow_mass_comp["nitrate"]
    ) == pytest.approx(1.129e-4, rel=1e-3)
    assert value(
        m.fs.dmbr.properties_treated[0].flow_mass_comp["tss"]
    ) == pytest.approx(2.317e-3, rel=1e-3)

    assert value(m.fs.costing.LCOW) == pytest.approx(0.3118654, rel=1e-3)  # in M$/m**3
    assert value(m.fs.costing.electricity_intensity) == pytest.approx(0.11302, rel=1e-3)
