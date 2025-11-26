import os
from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    value,
    Var,
    NonNegativeReals,
    assert_optimal_termination,
    Objective,
    units as pyunits,
)
from pyomo.environ import TransformationFactory
from pyomo.network import Arc
from pyomo.util.calc_var_value import calculate_variable_from_constraint

from idaes.core.surrogate.pysmo_surrogate import PysmoRBFTrainer, PysmoSurrogate
from idaes.core.surrogate.surrogate_block import SurrogateBlock
from idaes.core.util.misc import add_object_reference
import idaes.core.util.scaling as iscale
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.scaling import calculate_scaling_factors
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models import Product, Feed
from idaes.models.unit_models.mixer import (
    Mixer,
    MomentumMixingType,
    MaterialBalanceType,
)
from idaes.apps.grid_integration.multiperiod.multiperiod import MultiPeriodModel

from watertap.flowsheets.ccro.utils.ccro_utils import composition_calculator
from watertap.unit_models.pressure_changer import Pump
from watertap.property_models.NaCl_T_dep_prop_pack import NaClParameterBlock
from watertap.unit_models.reverse_osmosis_1D import (
    ReverseOsmosis1D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
from watertap.unit_models.reverse_osmosis_0D import (
    ConcentrationPolarizationType,
    MassTransferCoefficient,
)
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.core.util.initialization import *
from watertap.core.solvers import get_solver

from watertap.unit_models.reverse_osmosis_1D import (
    ReverseOsmosis1D,
    ReverseOsmosis1DData,
)
from watertap.unit_models.pseudo_steady_state import CCRO1D, DeadVolume0D
from watertap.unit_models.pseudo_steady_state.flushing import Flushing

from watertap.costing import (
    WaterTAPCosting,
    PumpType,
    MixerType,
    ROType,
)

from idaes.core.util.scaling import (
    calculate_scaling_factors,
    set_scaling_factor,
    constraint_scaling_transform,
)

from watertap.flowsheets.ccro.utils import *

__all__ = [
    "build_ccro_system",
    "create_ccro_multiperiod",
    "add_ccro_costing",
    "fix_overall_water_recovery",
    "setup_optimization",
]

solver = get_solver()


op_dict_default = dict(
    n_time_points=11,
    rho=1000,
    raw_feed_conc=5,  # g/L
    raw_feed_flowrate=1.8,  # L/min
    recycle_flowrate=49.1,  # L/min
    recycle_conc_start=11.7,
    temperature=25,  # °C
    p1_pressure_start=306,  # psi
    p2_pressure_start=306,
    p1_eff=0.8,
    p2_eff=0.8,
    A_comp=5.96e-12,
    B_comp=3.08e-08,
    membrane_area=7.9,  # m2
    membrane_length=1,  # m
    channel_height=0.0008636,
    spacer_porosity=0.9,
    dead_volume=0.035564,
    accumulation_time=60,
    single_pass_water_recovery=0.063,
    include_costing=True,
    flushing_efficiency=0.9,
)


def build_ccro_system(time_blk=None, configuration=None):
    """
    Build the CCRO steady state flowsheet
    """

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.configuration = configuration

    m.fs.properties = NaClParameterBlock()

    # Low concentration feed to the system
    m.fs.raw_feed = Feed(property_package=m.fs.properties)

    # Dead volume is included in both configurations
    m.fs.dead_volume = DeadVolume0D(property_package=m.fs.properties)

    # Touch relevant parameters
    m.fs.raw_feed.properties[0].flow_vol_phase
    m.fs.raw_feed.properties[0].conc_mass_phase_comp
    m.fs.dead_volume.dead_volume.properties_out[0].flow_vol_phase
    m.fs.dead_volume.dead_volume.properties_out[0].conc_mass_phase_comp
    m.fs.dead_volume.dead_volume.properties_in[0].conc_mass_phase_comp
    m.fs.dead_volume.dead_volume.properties_out[0].conc_mass_phase_comp
    m.fs.dead_volume.dead_volume.properties_in[0].flow_vol_phase

    if m.fs.configuration == "filtration":

        # Feed pump
        m.fs.P1 = Pump(property_package=m.fs.properties)
        # Recirculation pump
        m.fs.P2 = Pump(property_package=m.fs.properties)
        # Mixer to combine feed and recycle streams
        m.fs.M1 = Mixer(
            property_package=m.fs.properties,
            has_holdup=False,
            num_inlets=2,
            momentum_mixing_type=MomentumMixingType.equality,
        )

        m.fs.RO = ReverseOsmosis1D(
            property_package=m.fs.properties,
            has_pressure_change=True,
            pressure_change_type=PressureChangeType.calculated,
            mass_transfer_coefficient=MassTransferCoefficient.calculated,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
            transformation_scheme="BACKWARD",
            transformation_method="dae.finite_difference",
            finite_elements=4,
            module_type="spiral_wound",
            has_full_reporting=True,
        )

        m.fs.product = Product(property_package=m.fs.properties)

        # Add connections
        m.fs.raw_feed_to_P1 = Arc(
            source=m.fs.raw_feed.outlet, destination=m.fs.P1.inlet
        )

        m.fs.P1_to_M1 = Arc(source=m.fs.P1.outlet, destination=m.fs.M1.inlet_1)
        m.fs.P2_to_M1 = Arc(source=m.fs.P2.outlet, destination=m.fs.M1.inlet_2)

        m.fs.M1_to_RO = Arc(source=m.fs.M1.outlet, destination=m.fs.RO.inlet)

        m.fs.RO_permeate_to_product = Arc(
            source=m.fs.RO.permeate, destination=m.fs.product.inlet
        )

        m.fs.RO_retentate_to_dead_volume = Arc(
            source=m.fs.RO.retentate, destination=m.fs.dead_volume.inlet
        )
        m.fs.dead_volume_to_P2 = Arc(
            source=m.fs.dead_volume.outlet, destination=m.fs.P2.inlet
        )

        TransformationFactory("network.expand_arcs").apply_to(m)

        # Touch relevant parameters
        m.fs.M1.inlet_1_state[0].flow_vol_phase
        m.fs.M1.inlet_1_state[0].conc_mass_phase_comp
        m.fs.M1.inlet_2_state[0].flow_vol_phase
        m.fs.M1.inlet_2_state[0].conc_mass_phase_comp
        m.fs.M1.mixed_state[0].flow_vol_phase
        m.fs.M1.mixed_state[0].conc_mass_phase_comp
        
    elif m.fs.configuration == "flushing":

        m.fs.P1 = Pump(property_package=m.fs.properties)
        m.fs.P2 = Pump(property_package=m.fs.properties)

        # Add connections
        m.fs.raw_feed_to_P1 = Arc(
            source=m.fs.raw_feed.outlet, destination=m.fs.P1.inlet
        )
        m.fs.P1_to_dead_volume = Arc(
            source=m.fs.P1.outlet, destination=m.fs.dead_volume.inlet
        )

        m.fs.dead_volume_to_P2 = Arc(
            source=m.fs.dead_volume.outlet, destination=m.fs.P2.inlet
        )

        m.fs.flushing = Flushing()

        TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def create_ccro_multiperiod(n_time_points=10, include_costing=True, op_dict=None):
    """
    Create multiperiod model for CCRO system
    """

    watertap_solver = get_solver()

    mp = MultiPeriodModel(
        n_time_points=n_time_points,
        process_model_func=build_ccro_system,
        linking_variable_func=get_variable_pairs,
        initialization_func=fix_dof_and_initialize,
        unfix_dof_func=unfix_dof,
        solver=watertap_solver,
        outlvl=logging.WARNING,
    )

    mp.op_dict = op_dict
    mp.n_time_points = n_time_points
    mp.include_costing = include_costing

    configuration_list = ["filtration"] * (n_time_points - 1) + ["flushing"]
    flowsheet_options = {
        t: {"time_blk": t, "configuration": configuration_list[t]}
        for t in range(n_time_points)
    }

    # Build instances of the process model for each time period
    mp.build_multi_period_model(model_data_kwargs=flowsheet_options)

    if include_costing:
        mp.costing = WaterTAPCosting()
        for t, m in enumerate(mp.get_active_process_blocks(), 1):
            if t == 1:
                register_costed_unit(mp, m.fs.P1, register_electricity_flow_only=True)
                register_costed_unit(mp, m.fs.RO)
            elif t == n_time_points - 1:
                register_costed_unit(mp, m.fs.P2)
                register_costed_unit(mp, m.fs.P1)
            # Last time period is flushing
            elif t == n_time_points:
                register_costed_unit(mp, m.fs.P2, register_electricity_flow_only=True)
                register_costed_unit(mp, m.fs.P1, register_electricity_flow_only=True)

    for t, m in enumerate(mp.get_active_process_blocks(), 1):
        if t == 1:
            fix_dof_and_initialize(m, op_dict=op_dict)
            unfix_dof(m, unfix_dead_volume_state=False, op_dict=op_dict)
            old_m = m
        # Last time period is flushing
        if t == n_time_points:
            copy_time_period_links(old_m, m)
            fix_dof_and_initialize(m, op_dict=op_dict)
            unfix_dof(m, unfix_dead_volume_state=False, op_dict=op_dict)
        else:
            copy_state_prop_time_period_links(old_m, m)

        print(f"Period {t} DOF:", degrees_of_freedom(m))
        assert degrees_of_freedom(m) == 0

        results = solve(model=m, tee=True)
        assert_optimal_termination(results)
        unfix_dof(m, unfix_dead_volume_state=True, op_dict=op_dict)
        old_m = m

    print("Multi-period DOF:", degrees_of_freedom(mp))
    add_multiperiod_variables(mp)
    add_multiperiod_constraints(mp)

    if include_costing:
        mp.costing.cost_process()
        mp.costing.add_LCOW(mp.avg_product_flow_rate)
        mp.costing.add_specific_energy_consumption(mp.avg_product_flow_rate, name="SEC")

        mp.costing.initialize()

    scale_multiperiod_model(mp)
    calculate_scaling_factors(mp)

    return mp


def fix_overall_water_recovery(mp, overall_water_recovery):

    mp.overall_recovery.fix(overall_water_recovery)

    # Fixed for accumulation time for initialization
    for t, m in enumerate(mp.get_active_process_blocks()):
        m.fs.dead_volume.accumulation_time.unfix()
        # m.fs.dead_volume.accumulation_time.setlb(1)
        # m.fs.dead_volume.accumulation_time.setub(400)

        set_scaling_factor(m.fs.dead_volume.accumulation_time, 1e-2)

    # Equal accumulation time across all filtration periods
    @mp.Constraint(list(range(1, mp.n_time_points - 1)))
    def accumulation_time_cons(mp, t):
        blks = list(mp.get_active_process_blocks())
        return blks[t].fs.dead_volume.accumulation_time[0] == (
            blks[0].fs.dead_volume.accumulation_time[0]
        )


def setup_optimization(
    mp, overall_water_recovery=0.5, max_cycle_time_hr=1, recycle_flow_bounds=(0.1, 20)
):
    """
    Setup the multiperiod model for optimization.
    """

    m0 = mp.get_active_process_blocks()[0]

    fix_overall_water_recovery(mp, overall_water_recovery)
    mp.global_dead_volume_constraint.activate()
    mp.equal_dead_volume_constraint.activate()
    mp.equal_delta_dead_volume_constraint.activate()
    mp.ro_membrane_area_constraint.activate()
    mp.ro_membrane_length_constraint.activate()
    mp.total_cycle_time.setub(max_cycle_time_hr * pyunits.hours)
    # mp.total_cycle_time.setlb(0.98 * max_cycle_time_hr * pyunits.hours)
    mp.equal_recycle_rate.activate()

    for t, m in enumerate(mp.get_active_process_blocks(), 1):
        if m.fs.find_component("RO") is not None:
            m.fs.RO.length.unfix()

            m.fs.RO.width.unfix()
            m.fs.RO.length.setlb(0.1)

            m.fs.RO.width.setlb(0.1)

            # m.fs.RO.mixed_permeate[0].conc_mass_phase_comp.display()
            m.fs.RO.mixed_permeate[0].conc_mass_phase_comp["Liq", "NaCl"].setub(0.5)
            m.fs.RO.area.setlb(1)
            m.fs.RO.area.unfix()
            m.fs.RO.flux_mass_phase_comp[0.0, 1.0, "Liq", "H2O"].setlb(
                0.001 * pyunits.kg / pyunits.m**2 / pyunits.hr
            )

            # m.fs.RO.feed_side.velocity[0, 0].setub(0.3)
            # m.fs.RO.feed_side.velocity[0, 0].setub(0.3)
            m.fs.RO.feed_side.velocity[0, 0].setlb(0.05)

        m.fs.dead_volume.volume.unfix()
        m.fs.dead_volume.delta_state.volume[0, "Liq"].unfix()


        # if m.fs.configuration == "filtration":
        #     m.fs.RO.recovery_vol_phase[0.0, "Liq"].setub(0.1)


        if t != mp.n_time_points:
            print("unfixing P2 outlet flowrate at t =", t)
            m.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"].unfix()
            m.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"].setlb(
                recycle_flow_bounds[0] * pyunits.L / pyunits.min
            )
            m.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"].setub(
                recycle_flow_bounds[1] * pyunits.L / pyunits.min
            )
        # else:
        #     m.fs.raw_feed.properties[0].flow_vol_phase["Liq"].unfix()
        #     m.fs.raw_feed.properties[0].flow_vol_phase["Liq"].setlb(
        #         recycle_flow_bounds[0] * pyunits.L / pyunits.min
        #     )
        #     m.fs.raw_feed.properties[0].flow_vol_phase["Liq"].setub(
        #         recycle_flow_bounds[1] * pyunits.L / pyunits.min
        #     )
        #     m.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"].setlb(
        #         recycle_flow_bounds[0] * pyunits.L / pyunits.min
        #     )
        #     m.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"].setub(
        #         recycle_flow_bounds[1] * pyunits.L / pyunits.min
        #     )

    if mp.include_costing:
        mp.cost_objective = Objective(expr=mp.costing.LCOW)

    print("DOF for optimization:", degrees_of_freedom(mp))


if __name__ == "__main__":
    op_dict = dict(
        rho=1000,
        raw_feed_conc=5,  # g/L
        raw_feed_flowrate=1.8,  # L/min
        recycle_flowrate=49.1,  # L/min
        recycle_conc_start=11.7,
        temperature=298,  # K
        p1_pressure_start=306,  # psi
        p2_pressure_start=306,
        p1_eff=0.8,
        p2_eff=0.8,
        A_comp=5.96e-12,
        B_comp=3.08e-08,
        membrane_area=7.9,  # m2
        membrane_length=1,  # m
        channel_height=0.0008636,
        spacer_porosity=0.9,
        dead_volume=0.035564,
        accumulation_time=60,
        single_pass_water_recovery=0.063,
        include_costing=True,
        flushing_efficiency=0.9,
    )

    op_dict = config_op_dict(op_dict)
    mp = create_ccro_multiperiod(n_time_points=11, include_costing=True, op_dict=op_dict)
    results = solve(mp, tee=False)
    print_results_table(mp, w=15)


    setup_optimization(
        mp,
        overall_water_recovery=0.8,
        max_cycle_time_hr=1,
        recycle_flow_bounds=(0.1, 100),
    )

    results = solve(mp)
    print_results_table(mp, w=16)
