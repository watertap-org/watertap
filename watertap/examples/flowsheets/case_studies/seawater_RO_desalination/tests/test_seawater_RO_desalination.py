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
from watertap.examples.flowsheets.case_studies.seawater_RO_desalination.seawater_RO_desalination import (
    main,
)


# -----------------------------------------------------------------------------
@pytest.mark.component
def test_seawater_RO_desalination_pressure_exchanger():
    m = main(erd_type="pressure_exchanger")

    f = m.fs.feed
    assert pytest.approx(305.63, rel=1e-4) == value(f.flow_mass_comp[0, "H2O"])
    assert pytest.approx(10.822, rel=1e-4) == value(f.flow_mass_comp[0, "tds"])
    assert pytest.approx(9.2760e-3, rel=1e-4) == value(f.flow_mass_comp[0, "tss"])
    assert pytest.approx(0.3092, rel=1e-4) == value(f.flow_vol[0])

    p1 = m.fs.desalination.P1
    assert pytest.approx(0.8, rel=1e-4) == value(p1.efficiency_pump[0])
    assert pytest.approx(1.1536e6, rel=1e-4) == value(p1.work_mechanical[0])
    assert pytest.approx(6.9e6, rel=1e-4) == value(p1.deltaP[0])
    assert pytest.approx(70.0, rel=1e-4) == value(p1.ratioP[0])

    assert pytest.approx(132.15, rel=1e-4) == value(
        p1.inlet.flow_mass_phase_comp[0.0, "Liq", "H2O"]
    )
    assert pytest.approx(4.6800, rel=1e-4) == value(
        p1.inlet.flow_mass_phase_comp[0.0, "Liq", "TDS"]
    )
    assert pytest.approx(298.00, rel=1e-4) == value(p1.inlet.temperature[0.0])
    assert pytest.approx(1.0e5, rel=1e-4) == value(p1.inlet.pressure[0.0])
    assert pytest.approx(132.15, rel=1e-4) == value(
        p1.outlet.flow_mass_phase_comp[0.0, "Liq", "H2O"]
    )
    assert pytest.approx(4.6800, rel=1e-4) == value(
        p1.outlet.flow_mass_phase_comp[0.0, "Liq", "TDS"]
    )
    assert pytest.approx(298.00, rel=1e-4) == value(p1.outlet.temperature[0.0])
    assert pytest.approx(7.0e6, rel=1e-4) == value(p1.outlet.pressure[0.0])

    ro = m.fs.desalination.RO
    assert pytest.approx(13914.0, rel=1e-4) == value(ro.area)
    assert pytest.approx(0.43681, rel=1e-4) == value(
        ro.recovery_mass_phase_comp[0.0, "Liq", "H2O"]
    )
    assert pytest.approx(0.43293, rel=1e-4) == value(ro.recovery_vol_phase[0.0, "Liq"])
    assert pytest.approx(305.57, rel=1e-4) == value(
        ro.inlet.flow_mass_phase_comp[0.0, "Liq", "H2O"]
    )
    assert pytest.approx(10.822, rel=1e-4) == value(
        ro.inlet.flow_mass_phase_comp[0.0, "Liq", "TDS"]
    )
    assert pytest.approx(298.00, rel=1e-4) == value(ro.inlet.temperature[0.0])
    assert pytest.approx(7.0e6, rel=1e-4) == value(ro.inlet.pressure[0.0])
    assert pytest.approx(172.09, rel=1e-4) == value(
        ro.retentate.flow_mass_phase_comp[0.0, "Liq", "H2O"]
    )
    assert pytest.approx(10.792, rel=1e-4) == value(
        ro.retentate.flow_mass_phase_comp[0.0, "Liq", "TDS"]
    )
    assert pytest.approx(298.02, rel=1e-4) == value(ro.retentate.temperature[0.0])
    assert pytest.approx(6.7759e6, rel=1e-4) == value(ro.retentate.pressure[0.0])
    assert pytest.approx(133.48, rel=1e-4) == value(
        ro.permeate.flow_mass_phase_comp[0.0, "Liq", "H2O"]
    )
    assert pytest.approx(2.9557e-2, rel=1e-4) == value(
        ro.permeate.flow_mass_phase_comp[0.0, "Liq", "TDS"]
    )
    assert pytest.approx(298.02, rel=1e-4) == value(ro.permeate.temperature[0.0])
    assert pytest.approx(1.0132e5, rel=1e-4) == value(ro.permeate.pressure[0.0])

    muni = m.fs.municipal
    assert pytest.approx(147.59, rel=1e-4) == value(muni.electricity[0])
    assert pytest.approx(133.48, rel=1e-4) == value(
        muni.inlet.flow_mass_comp[0.0, "H2O"]
    )
    assert pytest.approx(2.9557e-2, rel=1e-4) == value(
        muni.inlet.flow_mass_comp[0.0, "tds"]
    )
    assert pytest.approx(133.48, rel=1e-4) == value(
        muni.outlet.flow_mass_comp[0.0, "H2O"]
    )
    assert pytest.approx(2.9557e-2, rel=1e-4) == value(
        muni.outlet.flow_mass_comp[0.0, "tds"]
    )

    assert value(m.LCOW) == pytest.approx(0.8411, rel=1e-5)


@pytest.mark.component
def test_seawater_RO_desalination_pump_as_turbine():
    m = main(erd_type="pump_as_turbine")

    f = m.fs.feed
    assert pytest.approx(305.63, rel=1e-4) == value(f.flow_mass_comp[0, "H2O"])
    assert pytest.approx(10.822, rel=1e-4) == value(f.flow_mass_comp[0, "tds"])
    assert pytest.approx(9.2760e-3, rel=1e-4) == value(f.flow_mass_comp[0, "tss"])
    assert pytest.approx(0.3092, rel=1e-4) == value(f.flow_vol[0])

    p1 = m.fs.desalination.P1
    assert pytest.approx(0.8, rel=1e-4) == value(p1.efficiency_pump[0])
    assert pytest.approx(2.6676e6, rel=1e-4) == value(p1.work_mechanical[0])
    assert pytest.approx(6.9e6, rel=1e-4) == value(p1.deltaP[0])
    assert pytest.approx(70.0, rel=1e-4) == value(p1.ratioP[0])

    assert pytest.approx(305.57, rel=1e-4) == value(
        p1.inlet.flow_mass_phase_comp[0.0, "Liq", "H2O"]
    )
    assert pytest.approx(10.822, rel=1e-4) == value(
        p1.inlet.flow_mass_phase_comp[0.0, "Liq", "TDS"]
    )
    assert pytest.approx(298.00, rel=1e-4) == value(p1.inlet.temperature[0.0])
    assert pytest.approx(1.0e5, rel=1e-4) == value(p1.inlet.pressure[0.0])
    assert pytest.approx(305.57, rel=1e-4) == value(
        p1.outlet.flow_mass_phase_comp[0.0, "Liq", "H2O"]
    )
    assert pytest.approx(10.822, rel=1e-4) == value(
        p1.outlet.flow_mass_phase_comp[0.0, "Liq", "TDS"]
    )
    assert pytest.approx(298.00, rel=1e-4) == value(p1.outlet.temperature[0.0])
    assert pytest.approx(7.0e6, rel=1e-4) == value(p1.outlet.pressure[0.0])

    ro = m.fs.desalination.RO
    assert pytest.approx(13914.0, rel=1e-4) == value(ro.area)
    assert pytest.approx(0.43681, rel=1e-4) == value(
        ro.recovery_mass_phase_comp[0.0, "Liq", "H2O"]
    )
    assert pytest.approx(0.43293, rel=1e-4) == value(ro.recovery_vol_phase[0.0, "Liq"])
    assert pytest.approx(305.57, rel=1e-4) == value(
        ro.inlet.flow_mass_phase_comp[0.0, "Liq", "H2O"]
    )
    assert pytest.approx(10.822, rel=1e-4) == value(
        ro.inlet.flow_mass_phase_comp[0.0, "Liq", "TDS"]
    )
    assert pytest.approx(298.00, rel=1e-4) == value(ro.inlet.temperature[0.0])
    assert pytest.approx(7.0e6, rel=1e-4) == value(ro.inlet.pressure[0.0])
    assert pytest.approx(172.09, rel=1e-4) == value(
        ro.retentate.flow_mass_phase_comp[0.0, "Liq", "H2O"]
    )
    assert pytest.approx(10.792, rel=1e-4) == value(
        ro.retentate.flow_mass_phase_comp[0.0, "Liq", "TDS"]
    )
    assert pytest.approx(298.02, rel=1e-4) == value(ro.retentate.temperature[0.0])
    assert pytest.approx(6.7759e6, rel=1e-4) == value(ro.retentate.pressure[0.0])
    assert pytest.approx(133.48, rel=1e-4) == value(
        ro.permeate.flow_mass_phase_comp[0.0, "Liq", "H2O"]
    )
    assert pytest.approx(2.9557e-2, rel=1e-4) == value(
        ro.permeate.flow_mass_phase_comp[0.0, "Liq", "TDS"]
    )
    assert pytest.approx(298.02, rel=1e-4) == value(ro.permeate.temperature[0.0])
    assert pytest.approx(1.0132e5, rel=1e-4) == value(ro.permeate.pressure[0.0])

    muni = m.fs.municipal
    assert pytest.approx(147.59, rel=1e-4) == value(muni.electricity[0])
    assert pytest.approx(133.48, rel=1e-4) == value(
        muni.inlet.flow_mass_comp[0.0, "H2O"]
    )
    assert pytest.approx(2.9557e-2, rel=1e-4) == value(
        muni.inlet.flow_mass_comp[0.0, "tds"]
    )
    assert pytest.approx(133.48, rel=1e-4) == value(
        muni.outlet.flow_mass_comp[0.0, "H2O"]
    )
    assert pytest.approx(2.9557e-2, rel=1e-4) == value(
        muni.outlet.flow_mass_comp[0.0, "tds"]
    )

    assert value(m.LCOW) == pytest.approx(1.11162, rel=1e-5)
