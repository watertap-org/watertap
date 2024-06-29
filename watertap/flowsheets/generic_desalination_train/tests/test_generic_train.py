from pyomo.environ import assert_optimal_termination, value
from watertap.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.flowsheets.generic_desalination_train import generic_train as gt

import pytest

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

# -----------------------------------------------------------------------------


@pytest.mark.unit
def test_standard_train():
    # This will test bilding all the units and
    # unit operations as well as opitmization
    # TODO: Might want to separate tests based on unit, and refine but probably not
    # needed as this is a very simple model
    m = gt.build()
    m.fs.feed.properties[0].temperature.fix(298.15)
    m.fs.feed.properties[0].pressure.fix(101325)
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(1)
    m.fs.feed.properties[0].conc_mass_phase_comp["Liq", "TDS"].fix(3.5)  # kg/m3 or g/L
    m.fs.feed.properties[0].conc_mass_phase_comp["Liq", "X"].fix(0.3)  # kg/m3 or g/L
    gt.update_feed(m.fs)

    gt.initialize(m)
    m.fs.Pretreatment.separator.component_removal_percent["X"].fix(50)
    m.fs.Pretreatment.separator.separation_cost["X"].fix(0.5)

    m.fs.Desal_1.desalter.water_recovery.fix(80)
    m.fs.Desal_2.desalter.water_recovery.fix(50)
    m.fs.Desal_2.desalter.recovery_cost.fix(0.01)
    m.fs.Desal_2.desalter.recovery_cost_offset.fix(35)
    # config as crystalizer
    m.fs.Desal_3.desalter.water_recovery.unfix()
    m.fs.Desal_3.desalter.brine_water_mass_percent.fix(80)

    m.fs.Valorizer.separator.product_value["X"].fix(1)
    m.fs.Valorizer.separator.component_removal_percent["X"].fix(50)

    result = gt.solve(m)
    assert_optimal_termination(result)
    assert degrees_of_freedom(m) == 0
    # test overall result costs operation
    assert value(m.fs.water_recovery) == pytest.approx(98.220, rel=1e-3)

    # test costing compute
    assert value(m.fs.costing.LCOW) == pytest.approx(1.2686, rel=1e-3)
    unit_costs = {
        "Product distribution": 0,
        "Waste disposal": 0.0180462225724005,
        "Feed sourcing": 0,
        "Pretreatment": 0.17817144384075734,
        "Desal_1": 0.24431276721281453,
        "Desal_2": 0.06616804112013724,
        "Desal_3": 0.8382712295193712,
        "Valorizer cost": 0,
        "Valorizer revenue": -0.07635918566495165,
    }
    for name, var in m.fs.costing.LCOW_unit.items():
        assert value(var) == pytest.approx(unit_costs[name], rel=1e-5)

    assert value(
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]
    ) == pytest.approx(1, rel=1e-3)
    assert value(
        m.fs.product.properties[0].flow_mass_phase_comp["Liq", "H2O"]
    ) == pytest.approx(0.9859465972001905, rel=1e-3)
    assert value(
        m.fs.disposal.properties[0].flow_mass_phase_comp["Liq", "H2O"]
    ) == pytest.approx(0.014053402799809303, rel=1e-5)

    # test desalter brine water content constraint
    assert value(m.fs.Desal_3.desalter.water_recovery) == pytest.approx(
        82.34735228114, rel=1e-3
    )

    # test pretreatment removal constraints
    assert value(
        m.fs.Pretreatment.separator.product_properties[0].flow_mass_phase_comp[
            "Liq", "X"
        ]
    ) == pytest.approx(0.00015057215854267565, rel=1e-6)

    # test valorizer removal constraints
    assert value(
        m.fs.Valorizer.separator.product_properties[0].flow_mass_phase_comp["Liq", "X"]
    ) == pytest.approx(7.528607927133832e-05, rel=1e-6)
