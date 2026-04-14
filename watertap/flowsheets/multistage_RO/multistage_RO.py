from pyomo.environ import (
    RangeSet,
    ConcreteModel,
    value,
    Constraint,
    Objective,
    Expression,
    Param,
    Var,
    TransformationFactory,
    assert_optimal_termination,
    check_optimal_termination,
    units as pyunits,
)
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import propagate_state
from idaes.models.unit_models import (
    Mixer,
    Separator,
    Product,
    Feed,
    StateJunction,
    Mixer,
)
from idaes.models.unit_models.mixer import MomentumMixingType
from idaes.core import UnitModelCostingBlock
import idaes.core.util.scaling as iscale

from watertap.property_models.NaCl_T_dep_prop_pack import NaClParameterBlock
from watertap.unit_models.pressure_changer import Pump, EnergyRecoveryDevice
from watertap.costing import WaterTAPCosting
# from watertap.core.util.model_diagnostics.infeasible import *
from watertap.unit_models.reverse_osmosis_1D import ReverseOsmosis1D
from watertap.core import (  # noqa # pylint: disable=unused-import
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    MembraneChannel1DBlock,
    PressureChangeType,
)
from watertap.core.solvers import get_solver

from watertap.flowsheets.multistage_RO.utils import (
    calculate_operating_pressure,
    overscale_ro,
    solve,
    relax_bounds_for_low_salinity_waters,
)

solver = get_solver()

"""
Build and solve n-stage RO system with optional ERD and booster pumps.
"""


def build_stage(
    blk,
    m=None,
    add_pump=False,
):
    if m is None:
        m = blk.model()

    blk.add_pump = add_pump

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

    # Touch flow vol and concentration for reporting
    for b in [blk.feed, blk.product, blk.disposal]:
        b.properties[0].flow_vol_phase["Liq"]
        b.properties[0].conc_mass_phase_comp["Liq", "NaCl"]

    if add_pump:
        blk.pump = Pump(property_package=m.fs.properties)
        blk.feed_to_pump = Arc(source=blk.feed.outlet, destination=blk.pump.inlet)
        blk.pump_to_RO = Arc(source=blk.pump.outlet, destination=blk.RO.inlet)
    else:
        blk.feed_to_RO = Arc(source=blk.feed.outlet, destination=blk.RO.inlet)

    blk.RO_to_disposal = Arc(source=blk.RO.retentate, destination=blk.disposal.inlet)
    blk.RO_to_product = Arc(source=blk.RO.permeate, destination=blk.product.inlet)

    blk.recovery = Expression(
        expr=blk.product.properties[0].flow_vol_phase["Liq"]
        / blk.feed.properties[0].flow_vol_phase["Liq"]
    )

    blk.flux = Var(
        initialize=25,
        bounds=(0.5, 50),
        units=pyunits.liter / pyunits.m**2 / pyunits.hour,
        doc="Stage flux",
    )

    blk.flux_constr = Constraint(
        expr=blk.flux
        == pyunits.convert(
            blk.product.properties[0].flow_vol_phase["Liq"] / blk.RO.area,
            to_units=pyunits.liter / pyunits.m**2 / pyunits.hour,
        )
    )

    TransformationFactory("network.expand_arcs").apply_to(blk)

    print(f"Degrees of freedom {blk.name}: {degrees_of_freedom(blk)}")


def set_stage_op_conditions(blk, m=None, max_pressure=200e5, ro_op_dict={}):

    if m is None:
        m = blk.model()

    stage_num = blk.index()

    if blk.add_pump:
        operating_pressure = calculate_operating_pressure(
            feed_state_block=m.fs.feed.properties[0],
            over_pressure=m.over_pressure,
            water_recovery=m.water_recovery,
            NaCl_passage=0.01,
            solver=solver,
        )
        # print(f"Estimated operating pressure: {operating_pressure*1e-5:.2f} bar")
        if operating_pressure >= max_pressure:
            operating_pressure = max_pressure

        if stage_num > 1 and m.fs.stage[stage_num - 1].add_pump:
            last_pump = m.fs.stage[stage_num - 1].pump
            c = Constraint(
                expr=blk.pump.control_volume.properties_out[0].pressure
                - last_pump.control_volume.properties_out[0].pressure
                >= 1 * pyunits.bar
            )
            blk.add_component(f"pressure_increase_constr_stage{stage_num}", c)

        print(f"Estimated operating pressure: {operating_pressure*1e-5:.2f} bar")

        blk.pump.control_volume.properties_out[0].pressure.setub(max_pressure)
        blk.pump.control_volume.properties_out[0].pressure.fix(operating_pressure)
        blk.pump.efficiency_pump.fix(0.8)

    ###################################
    # RO
    blk.RO.A_comp.fix(4.2e-12)
    blk.RO.B_comp.fix(3.5e-8)
    blk.RO.feed_side.channel_height.fix(1e-3)
    blk.RO.feed_side.spacer_porosity.fix(0.9)
    blk.RO.permeate.pressure[0].fix(101325)
    blk.RO.width.fix(5)
    blk.RO.area.fix(30)

    for p, val in ro_op_dict.items():
        v = blk.RO.find_component(p)
        if v is not None:
            v.fix(val)
        else:
            msg = f"Component {p} not found in RO unit model"
            raise ValueError(msg)


def scale_stage(blk, m=None):
    overscale_ro(blk.RO, m.fs.properties, full_scaling=True)

    if blk.add_pump:
        iscale.set_scaling_factor(blk.pump.control_volume.work, 1e-3)
        iscale.set_scaling_factor(
            blk.pump.control_volume.properties_out[0].flow_vol_phase["Liq"], 1
        )
        iscale.set_scaling_factor(blk.pump.work_fluid[0], 1)


def init_stage(blk, m=None):
    """
    Initialize one stage
    """

    if m is None:
        m = blk.model()

    blk.feed.initialize()

    if blk.add_pump:
        propagate_state(blk.feed_to_pump)
        blk.pump.initialize()
        propagate_state(blk.pump_to_RO)
    else:
        propagate_state(blk.feed_to_RO)

    blk.RO.initialize()
    propagate_state(blk.RO_to_disposal)

    blk.disposal.initialize()

    propagate_state(blk.RO_to_product)
    blk.product.initialize()

    print(f"DOF for {blk.name}: {degrees_of_freedom(blk)}")


def set_stage_bounds(blk):

    for phase in blk.RO.flux_mass_phase_comp:
        if "H2O" in phase:
            blk.RO.flux_mass_phase_comp[phase].setlb(
                0 * pyunits.kg / (pyunits.m**2 * pyunits.hr)
            )
            blk.RO.flux_mass_phase_comp[phase].setub(
                200 * pyunits.kg / (pyunits.m**2 * pyunits.hr)
            )

    blk.RO.length.setlb(1 * pyunits.meter)
    blk.RO.width.setlb(0.1 * pyunits.meter)
    blk.RO.feed_side.velocity[0, 0].setub(0.3)
    blk.RO.feed_side.velocity[0, 0].setlb(0.1)
    blk.RO.feed_side.cp_modulus.setub(20)
    blk.RO.feed_side.friction_factor_darcy.setub(20)


def add_costing(m, hpro_costing=False):

    cma_pump = {"pump_type": "high_pressure"}
    if hpro_costing:
        cma_ro = {"ro_type": "high_pressure"}

    else:
        cma_ro = {"ro_type": "standard"}

    m.fs.costing = WaterTAPCosting()

    for n, stage in m.fs.stage.items():
        if stage.add_pump:
            stage.pump.costing = UnitModelCostingBlock(
                flowsheet_costing_block=m.fs.costing, costing_method_arguments=cma_pump
            )
        stage.RO.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing, costing_method_arguments=cma_ro
        )
    if m.fs.find_component("ERD") is not None:
        m.fs.ERD.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    m.fs.costing.cost_process()

    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol_phase["Liq"])
    m.fs.costing.add_specific_energy_consumption(
        m.fs.product.properties[0].flow_vol_phase["Liq"], name="SEC"
    )


def run_n_stage_system(
    n_stages=2,
    flow_vol=1,  # L/s
    salinity=35,  # g/L
    water_recovery=0.5, # only for initialization
    over_pressure=0.25,
    perm_conc=0.5,
    pump_dict={1: True},  # first stage always has pump
    hpro_costing=False,
    add_erd=True,
    **kwargs,
):

    print(f"Running {n_stages} stage system with {sum(pump_dict.values())} pumps")

    m = ConcreteModel()
    m.flow_vol = flow_vol * pyunits.liter / pyunits.s
    m.salinity = salinity * pyunits.gram / pyunits.liter
    m.over_pressure = over_pressure
    m.water_recovery = water_recovery

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = NaClParameterBlock()

    m.fs.stages_set = RangeSet(n_stages)
    m.fs.stage = FlowsheetBlock(m.fs.stages_set, dynamic=False)

    m.fs.feed = Feed(property_package=m.fs.properties)

    tot_area_expr = 0
    tot_power_expr = 0

    for n, stage in m.fs.stage.items():
        build_stage(stage, m=m, add_pump=pump_dict.get(n, False))
        tot_area_expr += stage.RO.area
        tot_power_expr += stage.pump.work_mechanical[0] if stage.add_pump else 0

    if add_erd:
        m.fs.ERD = EnergyRecoveryDevice(property_package=m.fs.properties)

    m.fs.product = Product(property_package=m.fs.properties)
    m.fs.product_mixer = Mixer(
        property_package=m.fs.properties,
        momentum_mixing_type=MomentumMixingType.minimize,
        inlet_list=[f"stage{i}" for i in m.fs.stages_set],
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

    m.fs.perm_conc_lb = Constraint(
        expr=m.fs.product.properties[0].conc_mass_phase_comp["Liq", "NaCl"] <= perm_conc
    )

    m.fs.total_area = Expression(
        expr=tot_area_expr, doc="Total membrane area in system"
    )
    m.fs.specific_area = Expression(
        expr=m.fs.total_area / m.fs.product.properties[0].flow_vol_phase["Liq"],
        doc="Specific area of system (membrane area per unit product flow rate)",
    )

    m.fs.system_flux = Expression(
        expr=pyunits.convert(
            m.fs.product.properties[0].flow_vol_phase["Liq"] / m.fs.total_area,
            to_units=pyunits.liter / pyunits.m**2 / pyunits.hour,
        ),
        doc="System flux (based on total membrane area)",
    )

    m.fs.total_power = Expression(
        expr=pyunits.convert(tot_power_expr, to_units=pyunits.kW),
        doc="Total pumping power in system",
    )

    for b in [m.fs.feed, m.fs.product, m.fs.disposal]:
        b.properties[0].flow_vol_phase["Liq"]
        b.properties[0].conc_mass_phase_comp["Liq", "NaCl"]

    ### CONNECTIONS
    m.fs.feed_to_stage1 = Arc(
        source=m.fs.feed.outlet, destination=m.fs.stage[1].feed.inlet
    )

    for n, stage in m.fs.stage.items():
        if not n == m.fs.stages_set.last():
            # connect this stage to the next stage
            a = Arc(
                source=stage.disposal.outlet, destination=m.fs.stage[n + 1].feed.inlet
            )
            m.fs.add_component(f"stage{n}_to_stage{n+1}", a)
        # connect this stage permeate to the mixer
        a = Arc(
            source=stage.product.outlet,
            destination=m.fs.product_mixer.find_component(f"stage{n}"),
        )
        m.fs.add_component(f"stage{n}_product_to_mixer", a)
        if n == m.fs.stages_set.last() and add_erd:
            # connect this stage disposal to the ERD
            a = Arc(source=stage.disposal.outlet, destination=m.fs.ERD.inlet)
            m.fs.add_component(f"stage{n}_to_ERD", a)

    if add_erd:
        m.fs.ERD_to_disposal = Arc(
            source=m.fs.ERD.outlet, destination=m.fs.disposal.inlet
        )
    else:
        a = Arc(source=stage.disposal.outlet, destination=m.fs.disposal.inlet)
        m.fs.add_component(f"stage{n}_to_disposal", a)
    m.fs.product_mixer_to_product = Arc(
        source=m.fs.product_mixer.outlet, destination=m.fs.product.inlet
    )

    add_costing(m, hpro_costing=hpro_costing)

    TransformationFactory("network.expand_arcs").apply_to(m)

    #### SCALING
    rho = 1000 * pyunits.kg / pyunits.m**3  # approximate for scaling
    mass_flow_water = pyunits.convert(m.flow_vol * rho, to_units=pyunits.kg / pyunits.s)
    mass_flow_salt = pyunits.convert(
        m.flow_vol * m.salinity, to_units=pyunits.kg / pyunits.s
    )

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1 / value(mass_flow_water), index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(mass_flow_salt),
        index=("Liq", "NaCl"),
    )

    for n, stage in m.fs.stage.items():
        scale_stage(stage, m=m)

    if m.fs.find_component("ERD") is not None:
        iscale.set_scaling_factor(m.fs.ERD.control_volume.work, 1e-3)

    iscale.calculate_scaling_factors(m)

    #### SET OPERATING CONDITIONS

    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): m.flow_vol,
            ("conc_mass_phase_comp", ("Liq", "NaCl")): m.salinity,
            ("pressure", None): 101325,
            ("temperature", None): 298,
        },
        hold_state=True,
    )

    if salinity <= 10:
        ro_op_dict = {
            "A_comp": 5 * pyunits.liter / pyunits.m**2 / pyunits.hour / pyunits.bar,
            "B_comp": 0.5 * pyunits.liter / pyunits.m**2 / pyunits.hour,
        }
    else:
        ro_op_dict = {
            "A_comp": 1.5 * pyunits.liter / pyunits.m**2 / pyunits.hour / pyunits.bar,
            "B_comp": 0.1 * pyunits.liter / pyunits.m**2 / pyunits.hour,
        }

    for n, stage in m.fs.stage.items():
        set_stage_bounds(stage)
        set_stage_op_conditions(stage, m=m, ro_op_dict=ro_op_dict)

    if add_erd:
        m.fs.ERD.efficiency_pump.fix(0.95)
        m.fs.ERD.control_volume.properties_out[0].pressure.fix(101325)

    #### INITIALIZE

    print(f"DOF before initialization = {degrees_of_freedom(m)}")

    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_stage1)

    for n, stage in m.fs.stage.items():
        init_stage(stage, m=m)
        # stage.RO.area.unfix()
        # stage.RO.width.unfix()
        a = m.fs.find_component(f"stage{n}_product_to_mixer")
        propagate_state(a)
        if n != m.fs.stages_set.last():
            a = m.fs.find_component(f"stage{n}_to_stage{n+1}")
            propagate_state(a)
        if n == m.fs.stages_set.last() and add_erd:
            a = m.fs.find_component(f"stage{n}_to_ERD")
            propagate_state(a)
        elif n == m.fs.stages_set.last():
            a = m.fs.find_component(f"stage{n}_to_disposal")
            propagate_state(a)

    if add_erd:
        m.fs.ERD.initialize()
        propagate_state(m.fs.ERD_to_disposal)

    m.fs.disposal.initialize()

    m.fs.product_mixer.initialize()

    propagate_state(m.fs.product_mixer_to_product)
    m.fs.product.initialize()

    ### SOLVE
    m.fs.obj = Objective(expr=m.fs.costing.LCOW, sense="minimize")

    # Unfix pump pressure 
    for n, stage in m.fs.stage.items():
        stage.RO.area.unfix()
        stage.RO.width.unfix()
        if stage.add_pump:
            stage.pump.control_volume.properties_out[0].pressure.unfix()

    print(f"DOF before optimization = {degrees_of_freedom(m)}")

    results = solve(model=m, tee=True)
    assert_optimal_termination(results)

    # Fix area and pressure, optimize recovery
    for n, stage in m.fs.stage.items():
        if stage.add_pump:
            stage.pump.control_volume.properties_out[0].pressure.fix()
        stage.RO.area.fix()
        stage.RO.width.fix()

    print(f"FINAL DOF = {degrees_of_freedom(m)}")
    results = solve(model=m, tee=True)
    assert_optimal_termination(results)

    return m


def set_system_recovery(m, recovery):
    """
    Fix system recovery and unfix membrane area, width, and pump pressure in each stage
    """

    m.fs.system_recovery.fix(recovery)

    for n, stage in m.fs.stage.items():
        stage.RO.area.unfix()
        stage.RO.width.unfix()
        if stage.add_pump:
            stage.pump.control_volume.properties_out[0].pressure.unfix()

    return m


if __name__ == "__main__":

    sw_ro_op_dict = {
        "A_comp": 1.5 * pyunits.liter / pyunits.m**2 / pyunits.hour / pyunits.bar,
        "B_comp": 0.1 * pyunits.liter / pyunits.m**2 / pyunits.hour,
    }
    bw_ro_op_dict = {
        "A_comp": 5 * pyunits.liter / pyunits.m**2 / pyunits.hour / pyunits.bar,
        "B_comp": 0.5 * pyunits.liter / pyunits.m**2 / pyunits.hour,
    }

    m = run_n_stage_system(
        n_stages=2,
        salinity=35,
        water_recovery=0.5,
        pump_dict={1: True, 2: True},
        ro_op_dict=sw_ro_op_dict,
        add_erd=True,
    )
    m.fs.system_recovery.display()
    for n, stage in m.fs.stage.items():
        stage.RO.area.display()
        stage.RO.length.display()
        stage.RO.width.display()

    m = set_system_recovery(m, 0.5)

    res = solve(model=m, tee=False)

    # report_n_stage_system(m)
    
    m.fs.system_recovery.display()
    for n, stage in m.fs.stage.items():
        print(f"\nSTAGE {n}\n")
        stage.RO.area.display()
        stage.RO.length.display()
        stage.RO.width.display()

    # from collections import defaultdict
    # results = dict()
    # for x in [15, 20, 25, 30, 35, 40, 45]:
    #     results[x] = dict()
    #     for n in [1, 2, 3]:
    #         m = run_n_stage_system(
    #             n_stages=n,
    #             salinity=x,
    #             water_recovery=0.5,
    #             pump_dict={1: True, 2:True, 3:False},
    #             ro_op_dict=sw_ro_op_dict,
    #             add_erd=True,
    #         )
    #         m = set_system_recovery(m, 0.6)

    #         res = solve(model=m, tee=True)

    #         report_n_stage_system(m)
    #         results[x][n] = value(m.fs.costing.LCOW)

    # import pprint
    # pprint.pprint(results)
