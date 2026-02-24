from pyomo.environ import (
    RangeSet,
    ConcreteModel,
    value,
    Constraint,
    Objective,
    Param,
    TransformationFactory,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from watertap.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import EnergyBalanceType
from idaes.core.util.initialization import solve_indexed_blocks, propagate_state
from idaes.models.unit_models import (
    Mixer,
    Separator,
    Product,
    Feed,
    StateJunction,
    Mixer,
    Separator,
)
from idaes.models.unit_models.mixer import MomentumMixingType
from idaes.core import UnitModelCostingBlock
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.core.util.misc import StrEnum

from watertap.property_models.NaCl_T_dep_prop_pack import NaClParameterBlock
# from watertap.property_models.NaCl_prop_pack import NaClParameterBlock
from watertap.unit_models.reverse_osmosis_0D import (
    ReverseOsmosis0D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
from watertap.unit_models.reverse_osmosis_1D import (
    ReverseOsmosis1D,
    PressureChangeType,
    MassTransferCoefficient,
    ConcentrationPolarizationType,
)
from watertap.unit_models.pressure_exchanger import PressureExchanger
from watertap.unit_models.pressure_changer import Pump, EnergyRecoveryDevice
from watertap.core.util.initialization import assert_degrees_of_freedom
from watertap.costing import WaterTAPCosting
from watertap.core.util.model_diagnostics.infeasible import *


solver = get_solver()


def calculate_operating_pressure(
    feed_state_block=None,
    over_pressure=0.15,
    water_recovery=0.5,
    NaCl_passage=0.01,
    solver=None,
):
    """
    Estimate operating pressure for RO unit model given the following arguments:

    Arguments:
        feed_state_block:   the state block of the RO feed that has the non-pressure state
                            variables initialized to their values (default=None)
        over_pressure:  the amount of operating pressure above the brine osmotic pressure
                        represented as a fraction (default=0.15)
        water_recovery: the mass-based fraction of inlet H2O that becomes permeate
                        (default=0.5)
        NaCl_passage:   the mass-based fraction of inlet NaCl that becomes permeate
                        (default=0.01)
        solver:     solver object to be used (default=None)
    """
    t = ConcreteModel()  # create temporary model
    prop = feed_state_block.config.parameters
    t.brine = prop.build_state_block([0])

    # specify state block
    t.brine[0].flow_mass_phase_comp["Liq", "H2O"].fix(
        value(feed_state_block.flow_mass_phase_comp["Liq", "H2O"])
        * (1 - water_recovery)
    )
    t.brine[0].flow_mass_phase_comp["Liq", "NaCl"].fix(
        value(feed_state_block.flow_mass_phase_comp["Liq", "NaCl"]) * (1 - NaCl_passage)
    )
    t.brine[0].pressure.fix(
        101325
    )  # valid when osmotic pressure is independent of hydraulic pressure
    t.brine[0].temperature.fix(value(feed_state_block.temperature))

    # calculate osmotic pressure
    # since properties are created on demand, we must touch the property to create it
    t.brine[0].pressure_osm_phase
    # solve state block
    results = solve_indexed_blocks(solver, [t.brine])
    assert_optimal_termination(results)

    return value(t.brine[0].pressure_osm_phase["Liq"]) * (1 + over_pressure)


def build_model(
    flow_vol=1e-3, salt_mass_frac=35e-3, water_recovery=0.5, over_pressure=0.3
):

    m = ConcreteModel()

    m.flow_vol = flow_vol
    m.salt_mass_frac = salt_mass_frac
    m.water_recovery = water_recovery
    m.over_pressure = over_pressure

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = NaClParameterBlock()

    m.fs.feed = Feed(property_package=m.fs.properties)

    m.fs.P1 = Pump(property_package=m.fs.properties, energy_balance_type=EnergyBalanceType.none)
    m.fs.RO = ReverseOsmosis1D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        transformation_scheme="BACKWARD",
        transformation_method="dae.finite_difference",
        module_type="spiral_wound",
        finite_elements=10,
        has_full_reporting=True,
    )

    m.fs.ERD = EnergyRecoveryDevice(property_package=m.fs.properties)

    m.fs.product = Product(property_package=m.fs.properties)
    m.fs.disposal = Product(property_package=m.fs.properties)

    # Connections
    m.fs.feed_to_P1 = Arc(source=m.fs.feed.outlet, destination=m.fs.P1.inlet)
    m.fs.P1_to_RO = Arc(source=m.fs.P1.outlet, destination=m.fs.RO.inlet)
    m.fs.RO_to_ERD = Arc(source=m.fs.RO.retentate, destination=m.fs.ERD.inlet)
    m.fs.ERD_to_disposal = Arc(source=m.fs.ERD.outlet, destination=m.fs.disposal.inlet)
    m.fs.RO_to_product = Arc(source=m.fs.RO.permeate, destination=m.fs.product.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m)

    print(f"Degrees of freedom: {degrees_of_freedom(m)}")

    return m


def scale_model(m):
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1000 * m.flow_vol, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp",
        1e-3 / m.flow_vol / m.salt_mass_frac,
        index=("Liq", "NaCl"),
    )

    iscale.set_scaling_factor(m.fs.P1.control_volume.work, 1e-3)
    iscale.set_scaling_factor(
        m.fs.P1.control_volume.properties_out[0].flow_vol_phase["Liq"], 1
    )
    iscale.set_scaling_factor(m.fs.P1.work_fluid[0], 1)

    iscale.set_scaling_factor(m.fs.RO.area, 1e-2)
    iscale.set_scaling_factor(m.fs.RO.width, 1e-2)
    iscale.set_scaling_factor(m.fs.RO.feed_side.area, 1e2)
    iscale.set_scaling_factor(m.fs.ERD.control_volume.work, 1e-3)

    for e in m.fs.RO.feed_side.difference_elements:
        iscale.set_scaling_factor(
            m.fs.RO.mass_transfer_phase_comp[0, e, "Liq", "NaCl"], 1e4
        )
        iscale.set_scaling_factor(
            m.fs.RO.feed_side.mass_transfer_term[0, e, "Liq", "NaCl"], 1e4
        )

    for u in [m.fs.feed, m.fs.product, m.fs.disposal]:
        u.properties[0].flow_vol_phase["Liq"]
        u.properties[0].mass_frac_phase_comp["Liq", "NaCl"]
        u.properties[0].conc_mass_phase_comp["Liq", "NaCl"]

    iscale.calculate_scaling_factors(m)


def set_operating_conditions(m):

    # FEED
    m.fs.feed.properties[0].pressure.fix(101325)
    m.fs.feed.properties[0].temperature.fix(273.15 + 20)
    # m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(0.965)
    # m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].fix(0.035)

    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): m.flow_vol,
            ("mass_frac_phase_comp", ("Liq", "NaCl")): m.salt_mass_frac,
        },
        hold_state=True,
    )

    ###################################
    # PUMP 1
    operating_pressure = calculate_operating_pressure(
        feed_state_block=m.fs.feed.properties[0],
        over_pressure=m.over_pressure,
        water_recovery=m.water_recovery,
        NaCl_passage=0.01,
        solver=solver,
    )
    print(f"Estimated operating pressure: {operating_pressure*1e-5:.2f} bar")
    m.fs.P1.control_volume.properties_out[0].pressure.fix(operating_pressure)
    # m.fs.P1.control_volume.properties_out[0].pressure.setub(100*pyunits.bar)

    ###################################
    # RO

    m.fs.RO.A_comp.fix(4.2e-12)
    m.fs.RO.B_comp.fix(3.5e-8)
    m.fs.RO.feed_side.channel_height.fix(1e-3)
    m.fs.RO.feed_side.spacer_porosity.fix(0.85)
    m.fs.RO.permeate.pressure[0].fix(101325)
    m.fs.RO.width.fix(5)
    # m.fs.RO.length.fix(5)

    m.fs.RO.area.fix(60)

    ###################################
    # ERD

    m.fs.ERD.efficiency_pump.fix(0.95)
    m.fs.ERD.control_volume.properties_out[0].pressure.fix(101325)


def init_model(m):

    m.fs.RO.feed_side.properties[0, :].flow_mass_phase_comp["Liq", "H2O"] = value(
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]
    )
    m.fs.RO.feed_side.properties[0, :].flow_mass_phase_comp["Liq", "NaCl"] = value(
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"]
    )
    m.fs.RO.feed_side.properties[0, :].temperature = value(
        m.fs.feed.properties[0].temperature
    )
    m.fs.RO.feed_side.properties[0, :].pressure = value(
        m.fs.P1.control_volume.properties_out[0].pressure
    )
    m.fs.RO.initialize()
    m.fs.RO.area.unfix()
    m.fs.RO.recovery_mass_phase_comp[0, "Liq", "H2O"].fix(m.water_recovery)
    m.fs.RO.initialize()

    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_P1)

    m.fs.P1.initialize()
    propagate_state(m.fs.P1_to_RO)

    m.fs.RO.initialize()
    propagate_state(m.fs.RO_to_ERD)

    m.fs.ERD.initialize()
    propagate_state(m.fs.ERD_to_disposal)

    m.fs.disposal.initialize()

    propagate_state(m.fs.RO_to_product)
    m.fs.product.initialize()

    print(f"Degrees of freedom: {degrees_of_freedom(m)}")

    m.fs.RO.area.fix()


if __name__ == "__main__":
    flow_vol = 1e-3
    salt_mass_frac = 35e-3
    water_recovery = 0.5
    over_pressure = 0.3

    m = build_model(
        flow_vol=flow_vol,
        salt_mass_frac=salt_mass_frac,
        water_recovery=water_recovery,
        over_pressure=over_pressure,
    )
    scale_model(m)
    set_operating_conditions(m)
    init_model(m)

    # m.fs.RO.area.fix()
    print(f"Degrees of freedom: {degrees_of_freedom(m)}")
    assert degrees_of_freedom(m) == 0
    results = solver.solve(m)
    assert_optimal_termination(results)
