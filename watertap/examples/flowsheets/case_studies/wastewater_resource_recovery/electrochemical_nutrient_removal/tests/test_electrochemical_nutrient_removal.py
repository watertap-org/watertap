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
def test_electroNP():
    m, results = main()
    assert_optimal_termination(results)

    # test feed water flow
    assert value(m.fs.feed.properties[0].flow_mass_comp["H2O"]) == pytest.approx(
        10.512, rel=1e-3
    )
    assert value(m.fs.feed.properties[0].flow_mass_comp["nitrogen"]) == pytest.approx(
        0.00752, rel=1e-5
    )
    assert value(m.fs.feed.properties[0].flow_mass_comp["phosphorus"]) == pytest.approx(
        0.00752, rel=1e-5
    )
    assert value(m.fs.feed.properties[0].flow_mass_comp["struvite"]) == pytest.approx(
        0, rel=1e-10
    )

    # test struvite product flow
    assert value(
        m.fs.product_struvite.properties[0].flow_mass_comp["H2O"]
    ) == pytest.approx(0, rel=1e-10)
    assert value(
        m.fs.product_struvite.properties[0].flow_mass_comp["phosphorus"]
    ) == pytest.approx(0, rel=1e-10)

    # TODO - resolve discrepency between feed H2O mass flowrate and product_water H2O mass flowrate

    # test costing
    assert value(m.fs.costing.LCOW) == pytest.approx(76.191, rel=1e-3)  # in $/m**3
    assert value(m.fs.costing.LCOS) == pytest.approx(128.387, rel=1e-3)
