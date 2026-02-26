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
from watertap.property_models.NaCl_prop_pack import NaClParameterBlock
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
from watertap.flowsheets.ccro.utils.utils import calculate_operating_pressure

solver = get_solver()


def build_model(
    blk,
    m=None,
    add_booster=False,
):
    if m is None:
        m = blk.model()

    blk.add_booster = add_booster
    blk.feed = StateJunction(property_package=m.fs.properties)
    blk.RO = ReverseOsmosis1D(
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

    blk.product = StateJunction(property_package=m.fs.properties)
    blk.disposal = StateJunction(property_package=m.fs.properties)

    if add_booster:
        blk.pump = Pump(property_package=m.fs.properties)
        blk.feed_to_pump = Arc(source=blk.feed.outlet, destination=blk.pump.inlet)
        blk.pump_to_RO = Arc(source=blk.pump.outlet, destination=blk.RO.inlet)
    else:
        blk.feed_to_RO = Arc(source=blk.feed.outlet, destination=blk.RO.inlet)

    blk.RO_to_disposal = Arc(source=blk.RO.retentate, destination=blk.disposal.inlet)
    blk.RO_to_product = Arc(source=blk.RO.permeate, destination=blk.product.inlet)

    TransformationFactory("network.expand_arcs").apply_to(blk)

    print(f"Degrees of freedom {blk.name}: {degrees_of_freedom(blk)}")


def set_stage_op_conditions(blk, m=None):

    if m is None:
        m = blk.model()

    if blk.add_booster:
        operating_pressure = calculate_operating_pressure(
            feed_state_block=m.fs.feed.properties[0],
            over_pressure=m.over_pressure,
            water_recovery=m.water_recovery,
            NaCl_passage=0.01,
            solver=solver,
        )
        print(f"Estimated operating pressure: {operating_pressure*1e-5:.2f} bar")
        blk.pump.control_volume.properties_out[0].pressure.fix(operating_pressure)
        blk.pump.efficiency_pump.fix(0.8)

    ###################################
    # RO

    blk.RO.A_comp.fix(4.2e-12)
    blk.RO.B_comp.fix(3.5e-8)
    blk.RO.feed_side.channel_height.fix(1e-3)
    blk.RO.feed_side.spacer_porosity.fix(0.85)
    blk.RO.permeate.pressure[0].fix(101325)
    blk.RO.width.fix(5)
    blk.RO.area.fix(60)


def scale_stage(blk, m=None):

    if blk.add_booster:
        iscale.set_scaling_factor(blk.pump.control_volume.work, 1e-3)
        iscale.set_scaling_factor(
            blk.pump.control_volume.properties_out[0].flow_vol_phase["Liq"], 1
        )
        iscale.set_scaling_factor(blk.pump.work_fluid[0], 1)

    iscale.set_scaling_factor(blk.RO.area, 1e-2)
    iscale.set_scaling_factor(blk.RO.width, 1e-2)
    iscale.set_scaling_factor(blk.RO.feed_side.area, 1e2)

    for e in blk.RO.feed_side.difference_elements:
        iscale.set_scaling_factor(
            blk.RO.mass_transfer_phase_comp[0, e, "Liq", "NaCl"], 1e4
        )
        iscale.set_scaling_factor(
            blk.RO.feed_side.mass_transfer_term[0, e, "Liq", "NaCl"], 1e4
        )

    for b in [blk.feed, blk.product, blk.disposal]:
        b.properties[0].flow_vol_phase["Liq"]
        b.properties[0].conc_mass_phase_comp["Liq", "NaCl"]

    iscale.calculate_scaling_factors(blk)


def init_stage(blk, m=None):

    if m is None:
        m = blk.model()

    blk.feed.initialize()

    if blk.add_booster:
        propagate_state(blk.feed_to_pump)
        blk.pump.initialize()
        propagate_state(blk.pump_to_RO)
    else:
        propagate_state(blk.feed_to_RO)

    blk.RO.initialize()
    # propagate_state(blk.RO_to_ERD)
    propagate_state(blk.RO_to_disposal)

    # blk.ERD.initialize()
    # propagate_state(blk.ERD_to_disposal)

    blk.disposal.initialize()

    propagate_state(blk.RO_to_product)
    blk.product.initialize()

    print(f"Degrees of {blk.name}: {degrees_of_freedom(blk)}")


if __name__ == "__main__":

    n_stages = 2

    flow_vol = 1e-3
    salt_mass_frac = 35e-3
    water_recovery = 0.5
    over_pressure = 0.5

    m = ConcreteModel()
    m.flow_vol = flow_vol
    m.salt_mass_frac = salt_mass_frac
    m.over_pressure = over_pressure
    m.water_recovery = water_recovery

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = NaClParameterBlock()

    m.fs.stages_set = RangeSet(n_stages)
    m.fs.stage = FlowsheetBlock(m.fs.stages_set, dynamic=False)

    m.fs.feed = Feed(property_package=m.fs.properties)

    for n, stage in m.fs.stage.items():
        if n == m.fs.stages_set.first():
            build_model(stage, m=m, add_booster=True)
        else:
            build_model(stage, m=m)

    m.fs.ERD = EnergyRecoveryDevice(property_package=m.fs.properties)

    m.fs.product = Product(property_package=m.fs.properties)
    m.fs.product_mixer = Mixer(
        property_package=m.fs.properties,
        momentum_mixing_type=MomentumMixingType.minimize,
        inlet_list=[f"stage{i}" for i in range(1, n_stages + 1)],
    )

    m.fs.disposal = Product(property_package=m.fs.properties)

    m.fs.system_recovery = Var(
        initialize=0.5,
        bounds=(0.0001, 0.9999),
        units=pyunits.dimensionless,
        doc="System recovery",
    )

    m.fs.system_recovery_constraint = Constraint(
        expr=m.fs.system_recovery
        == m.fs.product.properties[0].flow_vol_phase["Liq"]
        / m.fs.feed.properties[0].flow_vol_phase["Liq"]
    )

    m.fs.perm_conc = Constraint(
        expr=m.fs.product.properties[0].conc_mass_phase_comp["Liq", "NaCl"] <= 0.25
    )

    ### CONNECTIONS
    m.fs.feed_to_stage1 = Arc(
        source=m.fs.feed.outlet, destination=m.fs.stage[1].feed.inlet
    )

    m.fs.stage1_to_stage2 = Arc(
        source=m.fs.stage[1].disposal.outlet, destination=m.fs.stage[2].feed.inlet
    )
    m.fs.product1_to_mixer = Arc(
        source=m.fs.stage[1].product.outlet, destination=m.fs.product_mixer.stage1
    )

    m.fs.stage2_to_ERD = Arc(
        source=m.fs.stage[2].disposal.outlet, destination=m.fs.ERD.inlet
    )
    m.fs.ERD_to_disposal = Arc(source=m.fs.ERD.outlet, destination=m.fs.disposal.inlet)

    m.fs.product2_to_mixer = Arc(
        source=m.fs.stage[2].product.outlet, destination=m.fs.product_mixer.stage2
    )

    m.fs.product_mixer_to_product = Arc(
        source=m.fs.product_mixer.outlet, destination=m.fs.product.inlet
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    #### SCALING
    print(f"dof after build = {degrees_of_freedom(m)}")

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1000 * flow_vol, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp",
        1e-3 / flow_vol / salt_mass_frac,
        index=("Liq", "NaCl"),
    )

    scale_stage(m.fs.stage[1], m=m)
    scale_stage(m.fs.stage[2], m=m)
    iscale.set_scaling_factor(m.fs.ERD.control_volume.work, 1e-3)
    iscale.calculate_scaling_factors(m)

    #### SET OPERATING CONDITIONS

    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): m.flow_vol,
            ("mass_frac_phase_comp", ("Liq", "NaCl")): m.salt_mass_frac,
            ("pressure", None): 101325,
            ("temperature", None): 298,
        },
        hold_state=True,
    )

    set_stage_op_conditions(m.fs.stage[1], m=m)
    set_stage_op_conditions(m.fs.stage[2], m=m)

    m.fs.ERD.efficiency_pump.fix(0.95)
    m.fs.ERD.control_volume.properties_out[0].pressure.fix(101325)

    m.fs.system_recovery.fix(0.5)

    #### INITIALIZE

    print(f"dof = {degrees_of_freedom(m)}")

    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_stage1)

    init_stage(m.fs.stage[1], m=m)
    m.fs.stage[1].RO.area.unfix()

    propagate_state(m.fs.stage1_to_stage2)
    m.fs.stage[2].RO.area.unfix()

    init_stage(m.fs.stage[2], m=m)

    propagate_state(m.fs.stage2_to_ERD)
    m.fs.ERD.initialize()

    propagate_state(m.fs.ERD_to_disposal)
    m.fs.disposal.initialize()

    propagate_state(m.fs.product1_to_mixer)
    propagate_state(m.fs.product2_to_mixer)
    m.fs.product_mixer.initialize()

    propagate_state(m.fs.product_mixer_to_product)
    m.fs.product.initialize()

    ### SOLVE
    m.fs.system_recovery.unfix()
    m.fs.obj = Objective(expr=m.fs.system_recovery, sense="maximize")
    print(f"dof = {degrees_of_freedom(m)}")

    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)

    m.fs.system_recovery.fix()
    m.fs.stage[1].RO.area.fix()

    assert_degrees_of_freedom(m, 0)
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
