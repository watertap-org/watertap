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
from pyomo.environ import value, units as pyunits
from watertap.flowsheets.seawater_RO_desalination.seawater_RO_desalination import (
    main,
)


@pytest.mark.component
def test_seawater_RO_desalination_pressure_exchanger():
    # NOTE: testing 0D RO by default
    m = main(erd_type="pressure_exchanger")

    f = m.fs.feed
    assert pytest.approx(0.42375, rel=1e-3) == value(
        m.fs.desalination.RO.recovery_vol_phase[0, "Liq"]
    )
    assert pytest.approx(316.30, rel=1e-4) == value(
        f.flow_mass_phase_comp[0, "Liq", "H2O"]
    )
    assert pytest.approx(10.822, rel=1e-4) == value(
        f.flow_mass_phase_comp[0, "Liq", "TDS"]
    )
    assert pytest.approx(9.2760e-3, rel=1e-4) == value(
        f.flow_mass_phase_comp[0, "Liq", "TSS"]
    )
    assert pytest.approx(0.320058, rel=1e-4) == value(f.properties[0].flow_vol)

    p1 = m.fs.desalination.P1
    assert pytest.approx(0.8, rel=1e-4) == value(p1.efficiency_pump[0])
    assert pytest.approx(1168370, rel=1e-4) == value(p1.work_mechanical[0])
    assert pytest.approx(6.9e6, rel=1e-4) == value(p1.deltaP[0])
    assert pytest.approx(70.0, rel=1e-4) == value(p1.ratioP[0])

    assert pytest.approx(133.875, rel=1e-4) == value(
        p1.inlet.flow_mass_phase_comp[0.0, "Liq", "H2O"]
    )
    assert pytest.approx(4.58127, rel=1e-4) == value(
        p1.inlet.flow_mass_phase_comp[0.0, "Liq", "TDS"]
    )
    assert pytest.approx(298.15, rel=1e-4) == value(p1.inlet.temperature[0.0])
    assert pytest.approx(1.0e5, rel=1e-4) == value(p1.inlet.pressure[0.0])
    assert pytest.approx(133.875, rel=1e-4) == value(
        p1.outlet.flow_mass_phase_comp[0.0, "Liq", "H2O"]
    )
    assert pytest.approx(4.5812, rel=1e-4) == value(
        p1.outlet.flow_mass_phase_comp[0.0, "Liq", "TDS"]
    )
    assert pytest.approx(298.15, rel=1e-4) == value(p1.outlet.temperature[0.0])
    assert pytest.approx(7.0e6, rel=1e-4) == value(p1.outlet.pressure[0.0])

    ro = m.fs.desalination.RO
    assert pytest.approx(13914.0, rel=1e-4) == value(ro.area)
    assert pytest.approx(0.4274, rel=1e-4) == value(
        ro.recovery_mass_phase_comp[0.0, "Liq", "H2O"]
    )
    assert pytest.approx(0.42375, rel=1e-4) == value(ro.recovery_vol_phase[0.0, "Liq"])
    assert pytest.approx(316.244, rel=1e-4) == value(
        ro.inlet.flow_mass_phase_comp[0.0, "Liq", "H2O"]
    )
    assert pytest.approx(10.822, rel=1e-4) == value(
        ro.inlet.flow_mass_phase_comp[0.0, "Liq", "TDS"]
    )
    assert pytest.approx(298.15, rel=1e-4) == value(ro.inlet.temperature[0.0])
    assert pytest.approx(7.0e6, rel=1e-4) == value(ro.inlet.pressure[0.0])
    assert pytest.approx(181.071, rel=1e-4) == value(
        ro.retentate.flow_mass_phase_comp[0.0, "Liq", "H2O"]
    )
    assert pytest.approx(10.794, rel=1e-4) == value(
        ro.retentate.flow_mass_phase_comp[0.0, "Liq", "TDS"]
    )
    assert pytest.approx(298.15, rel=1e-4) == value(ro.retentate.temperature[0.0])
    assert pytest.approx(6580200.27, rel=1e-4) == value(ro.retentate.pressure[0.0])
    assert pytest.approx(135.17, rel=1e-4) == value(
        ro.permeate.flow_mass_phase_comp[0.0, "Liq", "H2O"]
    )
    assert pytest.approx(0.0269839, rel=1e-4) == value(
        ro.permeate.flow_mass_phase_comp[0.0, "Liq", "TDS"]
    )
    assert pytest.approx(298.15, rel=1e-4) == value(ro.permeate.temperature[0.0])
    assert pytest.approx(1.0132e5, rel=1e-4) == value(ro.permeate.pressure[0.0])

    muni = m.fs.municipal
    assert pytest.approx(149.900, rel=1e-4) == value(muni.electricity[0])
    assert pytest.approx(135.17, rel=1e-4) == value(
        muni.inlet.flow_mass_phase_comp[0.0, "Liq", "H2O"]
    )
    assert pytest.approx(0.026983, rel=1e-4) == value(
        muni.inlet.flow_mass_phase_comp[0.0, "Liq", "TDS"]
    )
    assert pytest.approx(135.17, rel=1e-4) == value(
        muni.outlet.flow_mass_phase_comp[0.0, "Liq", "H2O"]
    )
    assert pytest.approx(0.0269839, rel=1e-4) == value(
        muni.outlet.flow_mass_phase_comp[0.0, "Liq", "TDS"]
    )

    assert m.fs.costing.base_currency == pyunits.USD_2018
    assert value(m.fs.costing.LCOW) == pytest.approx(0.84878, rel=1e-3)


@pytest.mark.component
def test_seawater_RO_desalination_pump_as_turbine():
    # NOTE: testing 0D RO by default
    m = main(erd_type="pump_as_turbine")

    f = m.fs.feed
    assert pytest.approx(0.42375, rel=1e-3) == value(
        m.fs.desalination.RO.recovery_vol_phase[0, "Liq"]
    )
    assert pytest.approx(316.30, rel=1e-4) == value(
        f.flow_mass_phase_comp[0, "Liq", "H2O"]
    )
    assert pytest.approx(10.822, rel=1e-4) == value(
        f.flow_mass_phase_comp[0, "Liq", "TDS"]
    )
    assert pytest.approx(9.2760e-3, rel=1e-4) == value(
        f.flow_mass_phase_comp[0, "Liq", "TSS"]
    )
    assert pytest.approx(0.320058, rel=1e-4) == value(f.properties[0].flow_vol)

    p1 = m.fs.desalination.P1
    assert pytest.approx(0.8, rel=1e-4) == value(p1.efficiency_pump[0])
    assert pytest.approx(2759955.01, rel=1e-4) == value(p1.work_mechanical[0])
    assert pytest.approx(6.9e6, rel=1e-4) == value(p1.deltaP[0])
    assert pytest.approx(70.0, rel=1e-4) == value(p1.ratioP[0])

    assert pytest.approx(316.24, rel=1e-4) == value(
        p1.inlet.flow_mass_phase_comp[0.0, "Liq", "H2O"]
    )
    assert pytest.approx(10.822, rel=1e-4) == value(
        p1.inlet.flow_mass_phase_comp[0.0, "Liq", "TDS"]
    )
    assert pytest.approx(298.15, rel=1e-4) == value(p1.inlet.temperature[0.0])
    assert pytest.approx(1.0e5, rel=1e-4) == value(p1.inlet.pressure[0.0])
    assert pytest.approx(316.24, rel=1e-4) == value(
        p1.outlet.flow_mass_phase_comp[0.0, "Liq", "H2O"]
    )
    assert pytest.approx(10.822, rel=1e-4) == value(
        p1.outlet.flow_mass_phase_comp[0.0, "Liq", "TDS"]
    )
    assert pytest.approx(298.15, rel=1e-4) == value(p1.outlet.temperature[0.0])
    assert pytest.approx(7.0e6, rel=1e-4) == value(p1.outlet.pressure[0.0])

    ro = m.fs.desalination.RO
    assert pytest.approx(13914.0, rel=1e-4) == value(ro.area)
    assert pytest.approx(0.4274, rel=1e-4) == value(
        ro.recovery_mass_phase_comp[0.0, "Liq", "H2O"]
    )
    assert pytest.approx(0.42375, rel=1e-4) == value(ro.recovery_vol_phase[0.0, "Liq"])
    assert pytest.approx(316.244, rel=1e-4) == value(
        ro.inlet.flow_mass_phase_comp[0.0, "Liq", "H2O"]
    )
    assert pytest.approx(10.822, rel=1e-4) == value(
        ro.inlet.flow_mass_phase_comp[0.0, "Liq", "TDS"]
    )
    assert pytest.approx(298.15, rel=1e-4) == value(ro.inlet.temperature[0.0])
    assert pytest.approx(7.0e6, rel=1e-4) == value(ro.inlet.pressure[0.0])
    assert pytest.approx(181.071, rel=1e-4) == value(
        ro.retentate.flow_mass_phase_comp[0.0, "Liq", "H2O"]
    )
    assert pytest.approx(10.794, rel=1e-4) == value(
        ro.retentate.flow_mass_phase_comp[0.0, "Liq", "TDS"]
    )
    assert pytest.approx(298.15, rel=1e-4) == value(ro.retentate.temperature[0.0])
    assert pytest.approx(6580200.27, rel=1e-4) == value(ro.retentate.pressure[0.0])
    assert pytest.approx(135.17, rel=1e-4) == value(
        ro.permeate.flow_mass_phase_comp[0.0, "Liq", "H2O"]
    )
    assert pytest.approx(0.0269839, rel=1e-4) == value(
        ro.permeate.flow_mass_phase_comp[0.0, "Liq", "TDS"]
    )
    assert pytest.approx(298.15, rel=1e-4) == value(ro.permeate.temperature[0.0])
    assert pytest.approx(1.0132e5, rel=1e-4) == value(ro.permeate.pressure[0.0])

    muni = m.fs.municipal
    assert pytest.approx(149.900, rel=1e-4) == value(muni.electricity[0])
    assert pytest.approx(135.17, rel=1e-4) == value(
        muni.inlet.flow_mass_phase_comp[0.0, "Liq", "H2O"]
    )
    assert pytest.approx(0.026983, rel=1e-4) == value(
        muni.inlet.flow_mass_phase_comp[0.0, "Liq", "TDS"]
    )
    assert pytest.approx(135.17, rel=1e-4) == value(
        muni.outlet.flow_mass_phase_comp[0.0, "Liq", "H2O"]
    )
    assert pytest.approx(0.0269839, rel=1e-4) == value(
        muni.outlet.flow_mass_phase_comp[0.0, "Liq", "TDS"]
    )

    assert m.fs.costing.base_currency == pyunits.USD_2018
    assert value(m.fs.costing.LCOW) == pytest.approx(0.96649, rel=1e-3)


@pytest.mark.component
def test_main_1D_pump_as_turbine():
    main(erd_type="pump_as_turbine", RO_1D=True)


@pytest.mark.component
def test_main_1D_pressure_exchanger():
    main(erd_type="pressure_exchanger", RO_1D=True)


flow_vols = [0.01, 0.1, 1]


@pytest.mark.parametrize("flow_vol", flow_vols)
def test_main_1D_pressure_exchanger_flow_rates(flow_vol):
    main(erd_type="pressure_exchanger", RO_1D=True, flow_vol=flow_vol)


@pytest.mark.parametrize("flow_vol", flow_vols)
def test_main_0D_pressure_exchanger_flow_rates(flow_vol):
    main(erd_type="pressure_exchanger", RO_1D=False, flow_vol=flow_vol)
