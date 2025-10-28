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

from watertap.unit_models.pseudo_steady_state import CCRO1D, DeadVolume0D
from watertap.unit_models.pseudo_steady_state.flushing import FlushingSurrogate

from idaes.core.surrogate.pysmo_surrogate import PysmoRBFTrainer, PysmoSurrogate
from idaes.core.surrogate.surrogate_block import SurrogateBlock

from watertap.costing import (
    WaterTAPCosting,
    PumpType,
    MixerType,
    ROType,
)

from watertap.flowsheets.ccro.utils import *

__all__ = [
    "build_system",
    "create_multiperiod",
    "add_multiperiod_constraints",
    "add_costing",
    "setup_optimization",
]

solver = get_solver()

# Get filepath for surrogate model
filepath = os.path.dirname(os.path.abspath(__file__))
surrogate_filename = os.path.join(filepath, "data/flushing_surrogate_n_5.json")

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
)


def build_system(time_blk=None, configuration=None):
    """
    Build the CCRO steady state flowsheet
    """

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.configuration = configuration

    m.fs.properties = NaClParameterBlock()

    m.fs.raw_feed = Feed(property_package=m.fs.properties)
    # Dead volume is included in both configurations
    m.fs.dead_volume = DeadVolume0D(property_package=m.fs.properties)

    m.fs.dead_volume.dead_volume.properties_in[0].flow_vol_phase["Liq"]
    m.fs.dead_volume.dead_volume.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"]

    if m.fs.configuration == "filtration":

        m.fs.product = Product(property_package=m.fs.properties)

        # Raw feed to the system

        m.fs.P1 = Pump(property_package=m.fs.properties)
        m.fs.P2 = Pump(property_package=m.fs.properties)

        m.fs.M1 = Mixer(
            property_package=m.fs.properties,
            has_holdup=False,
            num_inlets=2,
            momentum_mixing_type=MomentumMixingType.equality,
        )

        m.fs.RO = CCRO1D(
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

        # m.fs.single_pass_water_recovery = Var(
        #     initialize=single_pass_water_recovery,
        #     bounds=(0, 0.99),
        #     domain=NonNegativeReals,
        #     units=pyunits.dimensionless,
        #     doc="Water recovery after single pass through RO or accumulation time",
        # )

        m.fs.raw_feed.properties[0].flow_vol_phase
        m.fs.raw_feed.properties[0].conc_mass_phase_comp
        m.fs.M1.inlet_1_state[0].flow_vol_phase
        m.fs.M1.inlet_1_state[0].conc_mass_phase_comp
        m.fs.M1.inlet_2_state[0].flow_vol_phase
        m.fs.M1.inlet_2_state[0].conc_mass_phase_comp
        m.fs.M1.mixed_state[0].flow_vol_phase
        m.fs.M1.mixed_state[0].conc_mass_phase_comp
        m.fs.dead_volume.dead_volume.properties_out[0].flow_vol_phase
        m.fs.dead_volume.dead_volume.properties_out[0].conc_mass_phase_comp

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

        m.fs.flushing = FlushingSurrogate(surrogate_model_file=surrogate_filename)

        m.fs.dead_volume.dead_volume.properties_out[0].flow_vol_phase
        m.fs.dead_volume.dead_volume.properties_out[0].conc_mass_phase_comp

        TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def add_multiperiod_constraints(mp):
    """
    Add constraints to the multiperiod model.
    """

    # Get all filtration time blocks
    blks = list(mp.get_active_process_blocks())

    n_time_points = len(blks)
    accumulation_time = mp.op_dict["accumulation_time"]

    # # RO membrane area should be the same across all time periods - except flushing
    @mp.Constraint(list(range(1, n_time_points - 1)))
    def ro_membrane_area_constraint(mp, t):
        return blks[t].fs.RO.area == blks[0].fs.RO.area

    mp.ro_membrane_area_constraint.deactivate()

    # RO membrane length should be the same across all time periods - except flushing
    @mp.Constraint(list(range(1, n_time_points - 1)))
    def ro_membrane_length_constraint(mp, t):
        return blks[t].fs.RO.length == blks[0].fs.RO.length

    mp.ro_membrane_length_constraint.deactivate()

    for t in mp.ro_membrane_area_constraint:
        iscale.constraint_scaling_transform(mp.ro_membrane_area_constraint[t], 1e-1)
        iscale.constraint_scaling_transform(mp.ro_membrane_length_constraint[t], 1e-1)

    # RO membrane length should be the same across all time periods - except flushing
    @mp.Constraint(list(range(1, n_time_points)))
    def equal_dead_volume_constraint(mp, t):
        return (
            blks[t].fs.dead_volume.volume[0, "Liq"]
            == blks[0].fs.dead_volume.volume[0, "Liq"]
        )

    for t in mp.ro_membrane_area_constraint:
        iscale.constraint_scaling_transform(mp.equal_dead_volume_constraint[t], 1e2)

    @mp.Constraint(list(range(n_time_points)))
    def equal_delta_dead_volume_constraint(mp, t):
        return (
            blks[t].fs.dead_volume.volume[0, "Liq"]
            == blks[t].fs.dead_volume.delta_state.volume[0, "Liq"]
        )

    # RO membrane length should be the same across all time periods - except flushing
    @mp.Constraint(list(range(1, n_time_points - 1)))
    def equal_recycle_rate(mp, t):
        return (
            blks[t].fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"]
            == blks[0].fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"]
        )

    for t in mp.ro_membrane_area_constraint:
        iscale.constraint_scaling_transform(mp.equal_recycle_rate[t], 1e2)
    mp.equal_recycle_rate.deactivate()
    for t in mp.ro_membrane_area_constraint:
        iscale.constraint_scaling_transform(
            mp.equal_delta_dead_volume_constraint[t], 1e2
        )
    mp.equal_delta_dead_volume_constraint.deactivate()
    mp.equal_dead_volume_constraint.deactivate()

    mp.dead_volume_to_area_ratio = Var(
        initialize=1,
        domain=NonNegativeReals,
        units=pyunits.m**3 / pyunits.m**2,
        doc="Global dead volume over all time periods",
    )
    iscale.set_scaling_factor(mp.dead_volume_to_area_ratio, 1e2)
    mp.dead_volume_to_area_multiplier = Var(
        initialize=2,
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="Global dead volume over all time periods",
    )
    iscale.set_scaling_factor(mp.dead_volume_to_area_multiplier, 1)
    mp.dead_volume_to_area_multiplier.fix(1)
    mp.dead_volume_to_area_ratio.fix(
        value(1 * pyunits.m * (3.14 * 0.1016**2) * pyunits.m**2 / (7.2 * pyunits.m**2))
    )

    # mp.dead_volume_to_area_ratio.display()
    mp.global_dead_volume_constraint = Constraint(
        expr=blks[0].fs.dead_volume.volume[0, "Liq"]
        == blks[0].fs.RO.area
        * mp.dead_volume_to_area_ratio
        * mp.dead_volume_to_area_multiplier
    )
    calculate_variable_from_constraint(
        blks[0].fs.dead_volume.volume[0, "Liq"], mp.global_dead_volume_constraint
    )
    # blks[0].fs.dead_volume.volume.display()
    mp.global_dead_volume_constraint.deactivate()
    iscale.constraint_scaling_transform(mp.global_dead_volume_constraint, 1e2)

    # Density at the start of cycle should be the same as end of flushing (Initial condition and after flushing)
    @mp.Constraint()
    def cycle_end_density_constraint(mp):
        return (
            blks[0].fs.dead_volume.delta_state.dens_mass_phase[0, "Liq"]
            == blks[-1]
            .fs.dead_volume.dead_volume.properties_out[0]
            .dens_mass_phase["Liq"]
        )

    # Mass fraction at the start of cycle should be the same as end of flushing (Initial condition and after flushing)
    @mp.Constraint()
    def cycle_end_mass_frac_constraint(mp):
        return (
            blks[0].fs.dead_volume.delta_state.mass_frac_phase_comp[0, "Liq", "NaCl"]
            == blks[-1]
            .fs.dead_volume.dead_volume.properties_out[0]
            .mass_frac_phase_comp["Liq", "NaCl"]
        )

    # Permeate and feed
    mp.overall_recovery = Var(
        initialize=0.5,
        # bounds=(0, 1),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="Overall water recovery over all time periods",
    )

    mp.total_permeate = Var(
        initialize=0.5,
        domain=NonNegativeReals,
        units=pyunits.m**3,
        doc="Total permeate over all time periods",
    )

    mp.total_feed = Var(
        initialize=0.5,
        domain=NonNegativeReals,
        units=pyunits.m**3,
        doc="Total feed over all time periods",
    )

    mp.avg_product_flow_rate = Var(
        initialize=1,
        domain=NonNegativeReals,
        units=pyunits.m**3 / pyunits.s,
        doc="Average permeate production over all time periods",
    )
    mp.total_cycle_time = Var(
        initialize=1,
        domain=NonNegativeReals,
        units=pyunits.s,
        doc="Total cycle time including flushing",
    )
    mp.total_cycle_time_constraint = Constraint(
        expr=mp.total_cycle_time
        == sum(
            blks[t].fs.dead_volume.accumulation_time[0]
            for t in range(n_time_points - 1)
        )
        + blks[-1].fs.flushing.flushing_time
    )
    iscale.set_scaling_factor(mp.total_cycle_time, 1e-3)
    iscale.constraint_scaling_transform(mp.total_cycle_time_constraint, 1e-3)
    mp.final_concentration = Var(
        initialize=0.5,
        domain=NonNegativeReals,
        units=pyunits.g / pyunits.L,
        doc="Final concentration of the product stream",
    )
    mp.final_concentration_constraint = Constraint(
        expr=mp.final_concentration == blks[-1].fs.flushing.pre_flushing_concentration
    )
    iscale.set_scaling_factor(mp.final_concentration, 1e-1)
    iscale.constraint_scaling_transform(mp.final_concentration_constraint, 1e-1)

    # Total permeate
    @mp.Constraint()
    def total_permeate_constraint(mp):
        blks = list(mp.get_active_process_blocks())
        return mp.total_permeate == sum(
            blks[t].fs.product.properties[0].flow_vol_phase["Liq"]
            * blks[t].fs.dead_volume.accumulation_time[0]
            for t in range(n_time_points - 1)
        )

    @mp.Constraint()
    def eq_avg_product_flow_rate(mp):
        blks = mp.get_active_process_blocks()
        return mp.avg_product_flow_rate == pyunits.convert(
            mp.total_permeate
            / (
                accumulation_time * (len(blks) - 1) + blks[-1].fs.flushing.flushing_time
            ),
            to_units=pyunits.m**3 / pyunits.s,
        )

    # Total feed
    @mp.Constraint()
    def total_feed_constraint(mp):
        blks = list(mp.get_active_process_blocks())
        return (
            mp.total_feed
            == sum(
                blks[t].fs.raw_feed.properties[0].flow_vol_phase["Liq"]
                * blks[t].fs.dead_volume.accumulation_time[0]
                for t in range(n_time_points - 1)
            )
            + blks[0].fs.dead_volume.volume[0, "Liq"]
        )

    # Overall water recovery
    @mp.Constraint()
    def overall_water_recovery_constraint(mp):
        return mp.total_permeate == mp.overall_recovery * mp.total_feed


def create_multiperiod(n_time_points=10, include_costing=True, op_dict=None):
    """
    Create MultiPeriod model
    """

    watertap_solver = get_solver()

    mp = MultiPeriodModel(
        n_time_points=n_time_points,
        process_model_func=build_system,
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

    for t, m in enumerate(mp.get_active_process_blocks()):
        if t == 0:
            m0 = m
            fix_dof_and_initialize(m=m, op_dict=op_dict)
            unfix_dof(
                m=m, unfix_dead_volume_state=False, op_dict=op_dict
            )  # ensure we do not unfix dead volume stuff
            old_m = m
        # Last time period is flushing
        if t == n_time_points - 1:
            copy_time_period_links(old_m, m)
            fix_dof_and_initialize(m=m, op_dict=op_dict)
            unfix_dof(m=m, unfix_dead_volume_state=False, op_dict=op_dict)
        else:
            copy_state_prop_time_period_links(old_m, m)
        print("DOFs:", degrees_of_freedom(m))
        assert degrees_of_freedom(m) == 0
        results = solve(model=m, tee=False)
        assert_optimal_termination(results)
        unfix_dof(m=m, unfix_dead_volume_state=True, op_dict=op_dict)
        old_m = m

    print("DOF_multi_period:", degrees_of_freedom(mp))
    add_multiperiod_constraints(mp)

    calculate_scaling_factors(mp)

    if include_costing:
        add_costing(m=m0, mp=mp)
        m0.fs.costing.initialize()

    return mp


def add_costing(m=None, mp=None):
    """
    Add costing blocks to steady-state model.
    """

    m.fs.costing = WaterTAPCosting()
    costing_method_arguments = dict(mp=mp)

    m.fs.RO.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method_arguments=costing_method_arguments,
    )

    m.fs.costing.cost_process()

    m.fs.costing.add_LCOW(mp.avg_product_flow_rate)

    m.fs.costing.add_specific_energy_consumption(mp.avg_product_flow_rate, name="SEC")


def setup_optimization(mp, overall_water_recovery=0.5):
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
    mp.total_cycle_time.setub(1 * pyunits.hours)
    mp.equal_recycle_rate.activate()
    for t, m in enumerate(mp.get_active_process_blocks()):
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

            m.fs.RO.feed_side.velocity[0, 0].setub(0.3)
            m.fs.RO.feed_side.velocity[0, 0].setlb(0.05)
        m.fs.dead_volume.volume.unfix()
        m.fs.dead_volume.delta_state.volume[0, "Liq"].unfix()
        if t != len(mp.get_active_process_blocks()) - 1:
            print("unfixing P2 outlet flowrate at t =", t)
            m.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"].unfix()
            m.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"].setub(
                20 * pyunits.L / pyunits.min
            )
            m.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"].setlb(
                0.1 * pyunits.L / pyunits.min
            )
    # m.fs.dead_volume.volume.fix()
    if mp.include_costing:
        mp.cost_objective = Objective(expr=m0.fs.costing.LCOW)
    print("DOF_optimization:", degrees_of_freedom(mp))


if __name__ == "__main__":
    op_dict = dict(
        n_time_points=11,
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
    )
    op_dict = config_op_dict(op_dict)
    mp = create_multiperiod(n_time_points=10, include_costing=False, op_dict=op_dict)
    results = solve(mp)
    print_results_table(mp)
    setup_optimization(mp, overall_water_recovery=0.5)
