import os
import math
import numpy as np
from pyomo.environ import (
    ConcreteModel,
    value,
    Param,
    Var,
    Constraint,
    Objective,
    NonNegativeReals,
    TransformationFactory,
    Block,
    RangeSet,
    check_optimal_termination,
    units as pyunits,
)
from pyomo.network import Arc
from idaes.core import FlowsheetBlock, UnitModelCostingBlock

from watertap.core.solvers import get_solver
from idaes.core.util.initialization import propagate_state as _prop_state
import idaes.core.util.scaling as iscale
from idaes.core.util.scaling import (
    constraint_scaling_transform,
    calculate_scaling_factors,
    set_scaling_factor,
)

from watertap.core.util.model_diagnostics.infeasible import *
from watertap.unit_models.reverse_osmosis_1D import (
    ReverseOsmosis1D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)

from watertap.unit_models.pressure_changer import Pump
from watertap.core.util.initialization import *
from idaes.models.unit_models.mixer import (
    Mixer,
    MomentumMixingType,
    MaterialBalanceType,
)
from idaes.models.unit_models.separator import (
    SplittingType,
    EnergySplittingType,
)
from watertap.costing import (
    WaterTAPCosting,
    PumpType,
    MixerType,
    ROType,
)

# from analysisWaterTAP.utils.flowsheet_utils import *
from watertap.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery import (
    calculate_operating_pressure,
)

# from analysisWaterTAP.utils import flowsheet_utils as fsTool
# from analysisWaterTAP.flowsheets.lssro_oaro.costing.LSRRO_ORARO_costing import *
from idaes.models.unit_models import Product, Feed, StateJunction, Separator
from idaes.core.util.model_statistics import *
from watertap.costing import WaterTAPCosting
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.NaCl_prop_pack import NaClParameterBlock

__all__ = [
    "build_ro",
    "build_ro_stage",
    "init_ro_system",
    "init_ro_stage",
    "set_operating_conditions",
    "set_ro_system_operating_conditions",
    "add_ro_costing",
    "add_ro_scaling",
    "display_ro_system_build",
    "display_dof_breakdown",
    "display_flow_table",
    "report_RO",
    "print_RO_costing_breakdown",
]


def propagate_state(arc):
    _prop_state(arc)
    print(f"Propogation of {arc.source.name} to {arc.destination.name} successful.")
    arc.source.display()
    print(arc.destination.name)
    arc.destination.display()
    print("\n")


def _initialize(blk, verbose=False):
    if verbose:
        print("\n")
        print(
            f"{blk.name:<30s}{f'Degrees of Freedom at Initialization = {degrees_of_freedom(blk):<10.0f}'}"
        )
        print("\n")
    try:
        blk.initialize()
    except:
        print("----------------------------------\n")
        print(f"Initialization of {blk.name} failed.")
        print("\n----------------------------------\n")

        blk.report()
        print_infeasible_bounds(blk)
        print_close_to_bounds(blk)
        assert False


def print_RO_op_pressure_est(blk):
    solver = get_solver()
    operating_pressure = calculate_operating_pressure(
        feed_state_block=blk.feed.properties[0],
        over_pressure=0.15,
        water_recovery=0.8,
        NaCl_passage=0.01,
        solver=solver,
    )

    operating_pressure_psi = pyunits.convert(
        operating_pressure * pyunits.Pa, to_units=pyunits.psi
    )()
    operating_pressure_bar = pyunits.convert(
        operating_pressure * pyunits.Pa, to_units=pyunits.bar
    )()
    print(
        f"\nOperating Pressure Estimate = {round(operating_pressure_bar, 2)} bar = {round(operating_pressure_psi, 2)} psi\n"
    )


_log = idaeslog.getModelLogger("my_model", level=idaeslog.DEBUG, tag="model")


def init_system(m, verbose=True, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    print("\n\n-------------------- INITIALIZING SYSTEM --------------------\n\n")

    assert_no_degrees_of_freedom(m)
    _initialize(m.fs.feed)
    propagate_state(m.fs.feed_to_pump)

    m.fs.pump.initialize()
    propagate_state(m.fs.pump_to_ro)

    m.fs.ro.initialize()
    propagate_state(m.fs.ro_to_product)
    propagate_state(m.fs.ro_to_disposal)

    _initialize(m.fs.product)
    _initialize(m.fs.disposal)


def set_operating_conditions(
    m, Qin=None, Qout=None, Cin=None, water_recovery=None, ro_pressure=25e5
):
    print(
        "\n\n-------------------- SETTING SYSTEM OPERATING CONDITIONS --------------------\n\n"
    )
    if Cin is None:
        Cin = 35

    m.fs.water_recovery = Var(
        initialize=0.5,
        bounds=(0, 0.99),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="System Water Recovery",
    )

    m.fs.feed_salinity = Var(
        initialize=35,
        bounds=(0, 2000),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="System Water Recovery",
    )

    m.fs.feed_flow_mass = Var(
        initialize=1,
        bounds=(0.00001, 1e6),
        domain=NonNegativeReals,
        units=pyunits.kg / pyunits.s,
        doc="System Feed Flowrate",
    )

    m.fs.perm_flow_mass = Var(
        initialize=1,
        bounds=(0.00001, 1e6),
        domain=NonNegativeReals,
        units=pyunits.kg / pyunits.s,
        doc="System Produce Flowrate",
    )

    if water_recovery is not None:
        m.fs.water_recovery.fix(water_recovery)
    else:
        m.fs.water_recovery.unfix()

    feed_temperature = 273.15 + 20
    supply_pressure = 101325

    #     # initialize feed
    m.fs.feed.pressure[0].fix(supply_pressure)
    m.fs.feed.temperature[0].fix(feed_temperature)

    m.fs.eq_water_recovery = Constraint(
        expr=m.fs.feed.properties[0].flow_vol * m.fs.water_recovery
        == m.fs.product.properties[0].flow_vol
    )

    if Qin is not None:
        m.fs.feed_flow_mass.fix(Qin)

    iscale.set_scaling_factor(m.fs.feed_flow_mass, 1)
    m.fs.feed_salinity.fix(Cin)
    iscale.set_scaling_factor(m.fs.feed_salinity, 0.1)

    m.fs.feed_flow_constraint = Constraint(
        expr=m.fs.feed_flow_mass == m.fs.perm_flow_mass / m.fs.water_recovery
    )
    iscale.set_scaling_factor(m.fs.perm_flow_mass, 1)

    m.fs.nacl_mass_constraint = Constraint(
        expr=m.fs.feed.flow_mass_phase_comp[0, "Liq", "NaCl"] * 1000
        == m.fs.feed_flow_mass * m.fs.feed_salinity
    )

    m.fs.h2o_mass_constraint = Constraint(
        expr=m.fs.feed.flow_mass_phase_comp[0, "Liq", "H2O"]
        == m.fs.feed_flow_mass * (1 - m.fs.feed_salinity / 1000)
    )

    m.fs.feed.properties[0].flow_vol_phase["Liq"]
    m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"]

    m.fs.feed.flow_mass_phase_comp[0, "Liq", "NaCl"].value = (
        m.fs.feed_flow_mass.value * m.fs.feed_salinity.value / 1000
    )
    m.fs.feed.flow_mass_phase_comp[0, "Liq", "H2O"].value = (
        m.fs.feed_flow_mass.value * (1 - m.fs.feed_salinity.value / 1000)
    )

    scale_flow = calc_scale(m.fs.feed.flow_mass_phase_comp[0, "Liq", "H2O"].value)
    scale_tds = calc_scale(m.fs.feed.flow_mass_phase_comp[0, "Liq", "NaCl"].value)

    m.fs.RO_properties.set_default_scaling(
        "flow_mass_phase_comp", 10**-scale_flow, index=("Liq", "H2O")
    )
    m.fs.RO_properties.set_default_scaling(
        "flow_mass_phase_comp", 10**-scale_tds, index=("Liq", "NaCl")
    )


def calc_scale(value):
    return math.floor(math.log(value, 10))


def set_ro_system_operating_conditions(m, mem_area=10000, RO_pressure=70e5):
    print(
        "\n\n-------------------- SETTING RO OPERATING CONDITIONS --------------------\n\n"
    )
    mem_A = 8.4e-12  # membrane water permeability coefficient [m/s-Pa]
    mem_B = 3.5e-8  # membrane salt permeability coefficient [m/s]
    height = 1e-3  # channel height in membrane stage [m]
    spacer_porosity = 0.95  # spacer porosity in membrane stage [-]
    area = mem_area  # membrane area [m^2]
    length = 7  # effective membrane width [m]
    pressure_atm = 101325  # atmospheric pressure [Pa]
    pump_efi = 0.8  # pump efficiency [-]

    m.fs.pump.efficiency_pump.fix(pump_efi)
    m.fs.pump.control_volume.properties_out[0].pressure.fix(RO_pressure)

    m.fs.ro.A_comp.fix(mem_A)
    m.fs.ro.B_comp.fix(mem_B)
    m.fs.ro.area.fix(mem_area)  # membrane area (m^2)
    m.fs.ro.length.fix(length)
    # m.fs.ro.feed_side.velocity[0, 0].fix(0.35)
    m.fs.ro.feed_side.channel_height.fix(height)
    m.fs.ro.feed_side.spacer_porosity.fix(spacer_porosity)

    m.fs.ro.width.setub(100000)
    m.fs.ro.feed_side.friction_factor_darcy.setub(10)

    m.fs.ro.mixed_permeate[0].pressure.fix(pressure_atm)

    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")


def set_scaling(m):

    # set_scaling_factor(m.fs.pump.work, 1e5)
    set_scaling_factor(m.fs.ro.area, 1e5)
    set_scaling_factor(m.fs.ro.width, 1e5)
    set_scaling_factor(m.fs.ro.length, 1)
    iscale.calculate_scaling_factors(m)


def add_costing(m, blk=None, costing_blk=None):
    if costing_blk is None:
        costing_blk = m.fs.costing

    m.fs.pump.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    m.fs.ro.costing = UnitModelCostingBlock(flowsheet_costing_block=costing_blk)

    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol)
    m.fs.costing.initialize()


def solve(model, solver=None, tee=True, raise_on_failure=True):
    # ---solving---
    if solver is None:
        solver = get_solver()

    print("\n--------- SOLVING ---------\n")

    results = solver.solve(model, tee=tee)

    if check_optimal_termination(results):
        print("\n--------- OPTIMAL SOLVE!!! ---------\n")
        return results
    msg = (
        "The current configuration is infeasible. Please adjust the decision variables."
    )
    if raise_on_failure:
        print_infeasible_bounds(m)
        print_variables_close_to_bounds(m)
        print_infeasible_constraints(m)
        raise RuntimeError(msg)
    else:
        return results


def display_flow_table(blk):
    print("\n\n")
    print("RO System Flow Table")
    print(
        f'{"NODE":<34s}{"MASS FLOW RATE H2O (KG/S)":<30s}{"PRESSURE (BAR)":<20s}{"MASS FLOW RATE NACL (KG/S)":<30s}{"CONC. (G/L)":<20s}'
    )
    print(
        f'{"Feed":<34s}{blk.feed.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{value(pyunits.convert(blk.feed.properties[0.0].pressure, to_units=pyunits.bar)):<30.1f}{blk.feed.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{blk.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}'
    )
    print(
        f'{"Product":<34s}{blk.product.properties[0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(blk.product.properties[0].pressure, to_units=pyunits.bar)():<30.1f}{blk.product.properties[0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{blk.product.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}'
    )
    print(
        f'{"Disposal":<34s}{blk.disposal.properties[0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(blk.disposal.properties[0].pressure, to_units=pyunits.bar)():<30.1f}{blk.disposal.properties[0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{blk.disposal.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}'
    )

    for idx, stage in blk.stage.items():
        print(
            f'{"RO Stage " + str(idx) + " Feed":<34s}{stage.feed.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(stage.module.feed_side.properties[0, 0].pressure, to_units=pyunits.bar)():<30.1f}{stage.feed.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{stage.module.feed_side.properties[0,0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}'
        )
    for idx, stage in blk.stage.items():
        print(
            f'{"RO Stage " + str(idx) + " Permeate":<34s}{stage.permeate.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(stage.permeate.properties[0.0].pressure, to_units=pyunits.bar)():<30.1f}{stage.permeate.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{stage.module.mixed_permeate[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}'
        )
    for idx, stage in blk.stage.items():
        print(
            f'{"RO Stage " + str(idx) + " Retentate":<34s}{stage.retentate.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(stage.retentate.properties[0.0].pressure, to_units=pyunits.bar)():<30.1f}{stage.retentate.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{stage.module.feed_side.properties[0.0,1.0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}'
        )


def report_RO(m):
    print(f"\n\n-------------------- RO Report --------------------\n")
    print(f'{"Recovery":<30s}{value(100*m.fs.water_recovery):<10.1f}{"%"}')
    print(
        f'{"RO Operating Pressure":<30s}{value(pyunits.convert(m.fs.pump.control_volume.properties_out[0.0].pressure, to_units=pyunits.psi)):<10.1f}{"psi"}'
    )
    print(f'{"RO Membrane Area":<30s}{value(m.fs.ro.area):<10.1f}{"m^2"}')


def build_system():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.RO_properties = NaClParameterBlock()
    m.fs.costing = WaterTAPCosting()

    m.fs.feed = Feed(property_package=m.fs.RO_properties)
    m.fs.product = Product(property_package=m.fs.RO_properties)
    m.fs.disposal = Product(property_package=m.fs.RO_properties)
    m.fs.pump = Pump(property_package=m.fs.RO_properties)

    m.fs.ro = ReverseOsmosis1D(
        property_package=m.fs.RO_properties,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        transformation_scheme="BACKWARD",
        transformation_method="dae.finite_difference",
        finite_elements=10,
        has_full_reporting=True,
    )

    m.fs.feed_to_pump = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.pump.inlet,
    )

    m.fs.pump_to_ro = Arc(
        source=m.fs.pump.outlet,
        destination=m.fs.ro.inlet,
    )

    m.fs.ro_to_product = Arc(
        source=m.fs.ro.permeate,
        destination=m.fs.product.inlet,
    )

    m.fs.ro_to_disposal = Arc(
        source=m.fs.ro.retentate,
        destination=m.fs.disposal.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def optimize(
    m,
    water_recovery=0.5,
    fixed_pressure=None,
    ro_mem_area=None,
):

    m.fs.lcow_objective = Objective(expr=m.fs.costing.LCOW)

    if water_recovery is not None:
        print(f"\n------- Fixed Recovery at {100*water_recovery}% -------")
        m.fs.water_recovery.fix(water_recovery)
    else:
        lower_bound = 0.01
        upper_bound = 0.99
        print(f"\n------- Unfixed Recovery -------")
        print(f"Lower Bound: {lower_bound}")
        print(f"Upper Bound: {upper_bound}")
        m.fs.water_recovery.unfix()
        m.fs.water_recovery.setlb(0.01)
        m.fs.water_recovery.setub(0.99)

    if fixed_pressure is not None:
        print(f"\n------- Fixed RO Pump Pressure at {fixed_pressure} -------\n")
        m.fs.pump.control_volume.properties_out[0].pressure.fix(fixed_pressure)
    else:
        lower_bound = 100 * pyunits.psi
        upper_bound = 900 * pyunits.psi
        print(f"------- Unfixed RO Pump Pressure -------")
        print(f"Lower Bound: {value(lower_bound)} {pyunits.get_units(lower_bound)}")
        print(f"Upper Bound: {value(upper_bound)} {pyunits.get_units(upper_bound)}")
        m.fs.pump.control_volume.properties_out[0].pressure.unfix()
        m.fs.pump.control_volume.properties_out[0].pressure.setlb(lower_bound)
        m.fs.pump.control_volume.properties_out[0].pressure.setub(upper_bound)

    if ro_mem_area is not None:
        print(f"\n------- Fixed RO Membrane Area at {ro_mem_area} -------\n")
        m.fs.ro.area.fix(ro_mem_area)
    else:
        lower_bound = 1e3
        upper_bound = 2e5
        print(f"\n------- Unfixed RO Membrane Area -------")
        print(f"Lower Bound: {lower_bound} m2")
        print(f"Upper Bound: {upper_bound} m2")
        print("\n")
        m.fs.ro.area.unfix()
        m.fs.ro.area.setub(1e6)


if __name__ == "__main__":
    file_dir = os.path.dirname(os.path.abspath(__file__))
    m = build_system()
    set_operating_conditions(m, Qin=100, Cin=35)
    set_ro_system_operating_conditions(m, mem_area=10000, RO_pressure=50e5)
    set_scaling(m)
    init_system(m)
    add_costing(m)
    solve(m)
    report_RO(m)
    optimize(
        m,
        water_recovery=0.5,
    )
    solve(m)
    report_RO(m)
