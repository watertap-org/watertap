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
from watertap.unit_models.pseudo_steady_state.flushing import FlushingSurrogate

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

# # Get filepath for surrogate model
# filepath = os.path.dirname(os.path.abspath(__file__))
# surrogate_filename = os.path.join(
#     filepath, "data/flushing_surrogate_multiple_tau_n_2.json"
# )

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

        m.fs.M1 = Mixer(
            property_package=m.fs.properties,
            has_holdup=False,
            num_inlets=2,
            momentum_mixing_type=MomentumMixingType.equality,
        )

        # m.fs.RO = CCRO1D(
        #     property_package=m.fs.properties,
        #     has_pressure_change=True,
        #     pressure_change_type=PressureChangeType.calculated,
        #     mass_transfer_coefficient=MassTransferCoefficient.calculated,
        #     concentration_polarization_type=ConcentrationPolarizationType.calculated,
        #     transformation_scheme="BACKWARD",
        #     transformation_method="dae.finite_difference",
        #     finite_elements=4,
        #     module_type="spiral_wound",
        #     has_full_reporting=True,
        # )

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

        m.fs.flushing = FlushingSurrogate()  # surrogate_model_file=surrogate_filename)

        TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def add_multiperiod_variables(mp):
    """
    Add variables to the multiperiod model.
    """

    bs = list(mp.get_active_process_blocks())
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

    # Permeate and feed
    mp.overall_recovery = Var(
        initialize=0.5,
        # bounds=(0, 1),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="Overall water recovery over all time periods",
    )
    mp.apparent_recovery = Var(
        initialize=0.5,
        # bounds=(0, 1),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="Apparent water recovery over all time periods",  # from BLM report on OCWD pilot
    )

    mp.final_concentration = Var(
        initialize=0.5,
        domain=NonNegativeReals,
        units=pyunits.g / pyunits.L,
        doc="Final concentration of the product stream",
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

    # mp.total_brine = Var(
    #     initialize=0.5,
    #     domain=NonNegativeReals,
    #     units=pyunits.m**3,
    #     doc="Total brine over all time periods",
    # )

    mp.avg_product_flow_rate = Var(
        initialize=1,
        domain=NonNegativeReals,
        units=pyunits.m**3 / pyunits.s,
        doc="Average permeate production over all time periods",
    )

    mp.total_filtration_time = Var(
        initialize=1,
        domain=NonNegativeReals,
        units=pyunits.s,
        doc="Total filtration time excluding flushing",
    )

    mp.total_cycle_time = Var(
        initialize=1,
        domain=NonNegativeReals,
        units=pyunits.s,
        doc="Total cycle time including flushing",
    )

    mp.cycle_time_ratio = Var(
        initialize=1,
        domain=NonNegativeReals,
        bounds=(0, 1.0001),
        units=pyunits.dimensionless,
        doc="Ratio of total cycle time to filtration time",
    )


def scale_multiperiod_model(mp):

    iscale.set_scaling_factor(mp.dead_volume_to_area_ratio, 1e2)
    iscale.set_scaling_factor(mp.dead_volume_to_area_multiplier, 1)
    iscale.set_scaling_factor(mp.total_cycle_time, 1e-3)
    iscale.set_scaling_factor(mp.total_filtration_time, 1e-3)
    iscale.set_scaling_factor(mp.final_concentration, 1e-1)

    iscale.constraint_scaling_transform(mp.total_filtration_time_constraint, 1e-3)
    iscale.constraint_scaling_transform(mp.total_cycle_time_constraint, 1e-3)
    iscale.constraint_scaling_transform(mp.final_concentration_constraint, 1e-1)
    iscale.constraint_scaling_transform(mp.global_dead_volume_constraint, 1e2)

    for c in mp.ro_membrane_area_constraint.values():
        iscale.constraint_scaling_transform(c, 1e-1)
    for c in mp.ro_membrane_length_constraint.values():
        iscale.constraint_scaling_transform(c, 1e-1)
    for c in mp.equal_dead_volume_constraint.values():
        iscale.constraint_scaling_transform(c, 1e2)
    for c in mp.equal_delta_dead_volume_constraint.values():
        iscale.constraint_scaling_transform(c, 1e2)
    for c in mp.equal_recycle_rate.values():
        iscale.constraint_scaling_transform(c, 1e2)
    pass


def add_multiperiod_constraints(mp):
    """
    Add constraints to the multiperiod model.
    """

    # Get all filtration time blocks
    blks = list(mp.get_active_process_blocks())
    b0 = blks[mp.TIME.first()]

    accumulation_time = mp.op_dict["accumulation_time"]

    # RO membrane area should be the same across all time periods - except flushing
    @mp.Constraint(
        mp.TIME,
        doc="RO membrane area equality through all filtration periods",
    )
    def ro_membrane_area_constraint(b, t):
        if t in [b.TIME.first(), b.TIME.last()]:
            return Constraint.Skip
        return blks[t].fs.RO.area == b0.fs.RO.area

    mp.ro_membrane_area_constraint.deactivate()

    # RO membrane length should be the same across all time periods - except flushing
    @mp.Constraint(
        mp.TIME,
        doc="RO membrane length equality through all filtration periods",
    )
    def ro_membrane_length_constraint(b, t):
        if t in [b.TIME.first(), b.TIME.last()]:
            return Constraint.Skip
        return blks[t].fs.RO.length == b0.fs.RO.length

    mp.ro_membrane_length_constraint.deactivate()

    # Dead volume should be the same across all time periods
    @mp.Constraint(
        mp.TIME,
        doc="Dead volume equality through all filtration periods",
    )
    def equal_dead_volume_constraint(b, t):
        if t == b.TIME.first():
            return Constraint.Skip
        return (
            blks[t].fs.dead_volume.volume[0, "Liq"]
            == b0.fs.dead_volume.volume[0, "Liq"]
        )

    mp.equal_dead_volume_constraint.deactivate()

    @mp.Constraint(
        mp.TIME, doc="Delta dead volume equality through all filtration periods"
    )
    def equal_delta_dead_volume_constraint(b, t):
        return (
            blks[t].fs.dead_volume.volume[0, "Liq"]
            == blks[t].fs.dead_volume.delta_state.volume[0, "Liq"]
        )

    mp.equal_delta_dead_volume_constraint.deactivate()

    # Recycle rate should be the same across all time periods
    @mp.Constraint(
        mp.TIME,
        doc="Recycle rate equality through all filtration periods",
    )
    def equal_recycle_rate(b, t):
        if t in [b.TIME.first(), b.TIME.last()]:
            return Constraint.Skip
        return (
            blks[t].fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"]
            == blks[0].fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"]
        )

    mp.equal_recycle_rate.deactivate()

    @mp.Constraint(doc="Global dead volume constraint")
    def global_dead_volume_constraint(b):
        return (
            b0.fs.dead_volume.volume[0, "Liq"]
            == b0.fs.RO.area
            * b.dead_volume_to_area_ratio
            * b.dead_volume_to_area_multiplier
        )

    calculate_variable_from_constraint(
        b0.fs.dead_volume.volume[0, "Liq"], mp.global_dead_volume_constraint
    )

    mp.global_dead_volume_constraint.deactivate()

    # Density at the start of cycle should be the same as end of flushing
    # (Initial condition and after flushing)
    @mp.Constraint(
        doc="Density equality between end of flushing and start of filtration"
    )
    def cycle_end_density_constraint(b):
        return (
            b0.fs.dead_volume.delta_state.dens_mass_phase[0, "Liq"]
            == blks[-1]
            .fs.dead_volume.dead_volume.properties_out[0]
            .dens_mass_phase["Liq"]
        )

    # Mass fraction at the start of cycle should be the same as end of flushing
    # (Initial condition and after flushing)
    @mp.Constraint(
        doc="Mass fraction equality between end of flushing and start of filtration"
    )
    def cycle_end_mass_frac_constraint(b):
        return (
            b0.fs.dead_volume.delta_state.mass_frac_phase_comp[0, "Liq", "NaCl"]
            == blks[-1]
            .fs.dead_volume.dead_volume.properties_out[0]
            .mass_frac_phase_comp["Liq", "NaCl"]
        )

    @mp.Constraint(doc="Total filtration time constraint")
    def total_filtration_time_constraint(b):
        return b.total_filtration_time == sum(
            blks[t].fs.dead_volume.accumulation_time[0]
            for t in range(b.n_time_points - 1)
        )

    @mp.Constraint(doc="Total cycle time constraint")
    def total_cycle_time_constraint(b):
        return (
            b.total_cycle_time
            == b.total_filtration_time + blks[-1].fs.flushing.flushing_time
        )

    @mp.Constraint(doc="Cycle time ratio constraint")
    def cycle_time_ratio_constraint(b):
        return b.cycle_time_ratio == (b.total_filtration_time / b.total_cycle_time)

    @mp.Constraint(doc="Final concentration constraint")
    def final_concentration_constraint(b):
        return b.final_concentration == blks[-1].fs.flushing.pre_flushing_concentration

    # Total permeate
    @mp.Constraint(doc="Total permeate produced over all time periods")
    def total_permeate_constraint(b):
        return b.total_permeate == sum(
            blks[t].fs.product.properties[0].flow_vol_phase["Liq"]
            * blks[t].fs.dead_volume.accumulation_time[0]
            for t in range(b.n_time_points - 1)
        )

    # # Total brine
    # @mp.Constraint(doc="Total brine produced over all time periods")
    # def total_brine_constraint(mp):
    #     blks = list(mp.get_active_process_blocks())
    #     return mp.total_brine == sum(
    #         blks[t].fs.brine.properties[0].flow_vol_phase["Liq"]
    #         * blks[t].fs.dead_volume.accumulation_time[0]
    #         for t in range(n_time_points - 1)
    #     )

    @mp.Constraint(doc="Average product flow rate over all time periods")
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
    @mp.Constraint(doc="Total feed volume over all time periods")
    def total_feed_constraint(b):
        return (
            b.total_feed
            == sum(
                blks[t].fs.raw_feed.properties[0].flow_vol_phase["Liq"]
                * blks[t].fs.dead_volume.accumulation_time[0]
                for t in range(b.n_time_points - 1)
            )
            + b0.fs.dead_volume.volume[0, "Liq"]
        )

    # Overall water recovery
    @mp.Constraint(doc="Overall water recovery for system")
    def overall_water_recovery_constraint(b):
        return b.total_permeate == b.overall_recovery * b.total_feed


def create_multiperiod(n_time_points=10, include_costing=True, op_dict=None):
    """
    Create multiperiod model for CCRO system
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

    if include_costing:
        # there is probably a more elegant way to do this
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

    # if include_costing:
    #     add_costing(m=m0, mp=mp)
    #     add_object_reference(mp, "costing", m0.fs.costing)
    #     add_object_reference(mp, "RO_costing", m0.fs.RO.costing)
    #     m0.fs.costing.initialize()

    return mp


def register_costed_unit(
    mp, unit, costing_method_arguments={}, register_electricity_flow_only=False
):
    if register_electricity_flow_only:
        lb = unit.work_mechanical[0.0].lb
        # set lower bound to 0 to avoid negative defined flow warning when lb is not >= 0
        unit.work_mechanical.setlb(0)
        mp.costing.cost_flow(
            pyunits.convert(unit.work_mechanical[0.0], to_units=pyunits.kW),
            "electricity",
        )
        # set lower bound back to its original value that was assigned to lb
        unit.work_mechanical.setlb(lb)
    else:
        unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=mp.costing,
            costing_method_arguments=costing_method_arguments,
        )


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
    mp = create_multiperiod(n_time_points=11, include_costing=True, op_dict=op_dict)
    results = solve(mp, tee=False)
    print_results_table(mp, w=16)


    setup_optimization(
        mp,
        overall_water_recovery=0.8,
        max_cycle_time_hr=1,
        recycle_flow_bounds=(0.1, 100),
    )

    results = solve(mp)
    print_results_table(mp, w=16)
    # m = build_system(configuration="flushing")
