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
from pyomo.environ import value
from watertap.examples.flowsheets.full_treatment_train.flowsheet_components.examples import (
    full_example,
)


@pytest.mark.component
def test_flowsheet_mvp_NF_bypass_twostage_1():
    m = full_example.solve_flowsheet_mvp_NF(
        has_bypass=True,
        has_desal_feed=False,
        is_twostage=True,
        has_ERD=False,
        NF_type="ZO",
        NF_base="ion",
        RO_type="0D",
        RO_base="TDS",
        RO_level="detailed",
    )
    assert value(
        m.fs.NF.retentate.flow_mass_phase_comp[0, "Liq", "H2O"]
    ) == pytest.approx(0.1901, rel=1e-3)
    assert value(
        m.fs.NF.retentate.flow_mass_phase_comp[0, "Liq", "Ca"]
    ) == pytest.approx(1.669e-4, rel=1e-3)
    assert value(
        m.fs.RO2.retentate.flow_mass_phase_comp[0, "Liq", "H2O"]
    ) == pytest.approx(0.2408, rel=1e-3)
    assert value(
        m.fs.RO2.retentate.flow_mass_phase_comp[0, "Liq", "TDS"]
    ) == pytest.approx(2.605e-2, rel=1e-3)
    assert value(m.fs.desal_saturation.saturation_index) == pytest.approx(
        0.4073, rel=1e-3
    )


@pytest.mark.component
def test_flowsheet_mvp_cost_optimization():
    kwargs_flowsheet = {
        "has_bypass": True,
        "has_desal_feed": False,
        "is_twostage": True,
        "has_ERD": True,
        "NF_type": "ZO",
        "NF_base": "ion",
        "RO_type": "0D",
        "RO_base": "TDS",
        "RO_level": "detailed",
    }
    m = full_example.solve_optimization(system_recovery=0.68, **kwargs_flowsheet)
    assert value(
        m.fs.mixer_permeate.outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
    ) == pytest.approx(0.6778, rel=1e-3)
    assert value(m.fs.costing.LCOW) == pytest.approx(0.5569, rel=1e-3)
