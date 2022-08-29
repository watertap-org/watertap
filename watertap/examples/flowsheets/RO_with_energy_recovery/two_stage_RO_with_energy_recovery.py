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
from pyomo.environ import (
    ConcreteModel,
    value,
    Constraint,
    Expression,
    Objective,
    Param,
    TransformationFactory,
    units as pyunits,
    assert_optimal_termination,
)
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import solve_indexed_blocks, propagate_state
from idaes.models.unit_models import Mixer, Separator, Product, Feed
from idaes.models.unit_models.mixer import MomentumMixingType
from idaes.core import UnitModelCostingBlock
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

import watertap.property_models.NaCl_prop_pack as props
from watertap.unit_models.reverse_osmosis_0D import (
    ReverseOsmosis0D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
from watertap.unit_models.pressure_exchanger import PressureExchanger
from watertap.unit_models.pressure_changer import Pump
from watertap.core.util.initialization import assert_degrees_of_freedom
from watertap.costing import WaterTAPCosting


def main():
    # set up solver
    solver = get_solver()

    # build, set, and initialize
    m = build()
    set_operating_conditions(m, water_recovery=0.5, over_pressure=0.3, solver=solver)
    initialize_system(m, solver=solver)

    # simulate and display
    solve(m, solver=solver)
    print("\n***---Simulation results---***")
    display_system(m)
    display_design(m)
    display_state(m)

    # optimize and display
    optimize_set_up(m)
    optimize(m, solver=solver)
    print("\n***---Optimization results---***")
    display_system(m)
    display_design(m)
    display_state(m)


def build():
    # flowsheet set up
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic = False)
    m.fs.properties = props.NaClParameterBlock()
    m.fs.costing = WaterTAPCosting()

    # unit models (1st stage)
    m.fs.feed = Feed(property_package = m.fs.properties)
    m.fs.S1 = Separator(
        property_package = m.fs.properties, outlet_list = ["P1", "PXR1"]
    )
    m.fs.P1 = Pump(property_package = m.fs.properties)
    m.fs.PXR1 = PressureExchanger(property_package = m.fs.properties)
    m.fs.P2 = Pump(property_package = m.fs.properties)
    m.fs.M1 = Mixer(
            property_package = m.fs.properties,
            momentum_mixing_type = MomentumMixingType.equality,  # booster pump will match pressure
            inlet_list = ["P1", "P2"],
    )
    m.fs.RO1 = ReverseOsmosis0D(
            property_package = m.fs.properties,
            has_pressure_change = True,
            pressure_change_type = PressureChangeType.calculated,
            mass_transfer_coefficient = MassTransferCoefficient.calculated,
            concentration_polarization_type = ConcentrationPolarizationType.calculated,
    )
    m.fs.disposal1 = Product(property_package = m.fs.properties)

    # unit models (2nd stage)
    m.fs.P3 = Pump(property_package = m.fs.properties)
    m.fs.RO2 = ReverseOsmosis0D(
            property_package = m.fs.properties,
            has_pressure_change = True,
            pressure_change_type = PressureChangeType.calculated,
            mass_transfer_coefficient = MassTransferCoefficient.calculated,
            concentration_polarization_type = ConcentrationPolarizationType.calculated,
    )
    m.fs.product = Product(property_package = m.fs.properties)
    m.fs.disposal2 = Product(property_package = m.fs.properties)


    # costing (1st stage)
    m.fs.P1.costing = UnitModelCostingBlock(
        flowsheet_costing_block = m.fs.costing
    )
    m.fs.P2.costing = UnitModelCostingBlock(
        flowsheet_costing_block = m.fs.costing
    )

    m.fs.RO1.costing = UnitModelCostingBlock(
        flowsheet_costing_block = m.fs.costing
    )
    m.fs.PXR1.costing = UnitModelCostingBlock(
        flowsheet_costing_block = m.fs.costing
    )

    # costing (2nd stage)
    m.fs.P3.costing = UnitModelCostingBlock(
        flowsheet_costing_block = m.fs.costing
    )
    m.fs.RO2.costing = UnitModelCostingBlock(
        flowsheet_costing_block = m.fs.costing
    )
    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(m.fs.product.properties[0].flow_vol)
    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(m.fs.product.properties[0].flow_vol)

    # connections (1st stage)
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.S1.inlet)
    m.fs.s02 = Arc(source=m.fs.S1.P1, destination=m.fs.P1.inlet)
    m.fs.s03 = Arc(source=m.fs.P1.outlet, destination=m.fs.M1.P1)
    m.fs.s04 = Arc(source=m.fs.M1.outlet, destination=m.fs.RO1.inlet)
    m.fs.s05 = Arc(source=m.fs.RO1.permeate, destination=m.fs.P3.inlet)
    m.fs.s06 = Arc(source=m.fs.RO1.retentate, destination=m.fs.PXR1.high_pressure_inlet)
    m.fs.s07 = Arc(source=m.fs.PXR1.high_pressure_outlet, destination=m.fs.disposal1.inlet)
    m.fs.s08 = Arc(source=m.fs.S1.PXR1, destination=m.fs.PXR1.low_pressure_inlet)
    m.fs.s09 = Arc(source=m.fs.PXR1.low_pressure_outlet, destination=m.fs.P2.inlet)
    m.fs.s10 = Arc(source=m.fs.P2.outlet, destination=m.fs.M1.P2)

    # connections (2nd stage)
    m.fs.s11 = Arc(source=m.fs.P3.outlet, destination=m.fs.RO2.inlet)
    m.fs.s12 = Arc(source=m.fs.RO2.permeate, destination=m.fs.product.inlet)
    m.fs.s13 = Arc(source=m.fs.RO2.retentate, destination=m.fs.disposal2.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling
    # set default property values
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1e2, index=("Liq", "NaCl"))

    # set unit model values (1st stage)
    iscale.set_scaling_factor(m.fs.P1.control_volume.work, 1e-3)
    iscale.set_scaling_factor(m.fs.P2.control_volume.work, 1e-3)
    iscale.set_scaling_factor(m.fs.PXR1.low_pressure_side.work, 1e-3)
    iscale.set_scaling_factor(m.fs.PXR1.high_pressure_side.work, 1e-3)
    iscale.set_scaling_factor(m.fs.RO1.area, 1e1)

    # set unit model values (2nd stage)
    iscale.set_scaling_factor(m.fs.P3.control_volume.work, 1e-3)
    iscale.set_scaling_factor(m.fs.RO2.area, 1e1)

    # touch properties used in specifying and initializing the model (1st stage)
    m.fs.feed.properties[0].flow_vol_phase["Liq"]
    m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"]
    m.fs.S1.mixed_state[0].mass_frac_phase_comp
    m.fs.S1.PXR1_state[0].flow_vol_phase["Liq"]

    # should RO1 properties be touched? If so, what's the naming convention? prop_out[0]? permeate vs permeate_side?
    # m.fs.RO1.permeate.properties[0].flow_vol_phase["Liq"]
    # m.fs.RO1.permeate.properties[0].mass_frac_phase_comp["Liq", "NaCl"]

    # unused scaling factors needed by IDAES base costing module
    # calculate and propagate scaling factors
    iscale.calculate_scaling_factors(m)

    return m

def set_operating_conditions(m, water_recovery=0.5, over_pressure=0.3, solver=None):
    if solver is None:
        solver = get_solver()

    # ---specifications---
    # feed
    # state variables
    m.fs.feed.properties[0].pressure.fix(101325)  # stage 1 feed pressure [Pa]
    m.fs.feed.properties[0].temperature.fix(273.15 + 25)  # stage 1 feed temperature [K]
    # properties (cannot be fixed for initialization routines, must calculate the state variables)
    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): 1e-3,  # feed volumetric flow rate [m3/s]
            ("mass_frac_phase_comp", ("Liq", "NaCl")): 0.035,
        },  # feed NaCl mass fraction [-]
        hold_state=True,  # fixes the calculated component mass flow rates
    )

    # separator, no degrees of freedom (i.e. equal flow rates in PXR determines split fraction)

    # pumps 1 & 3, high pressure pumps, 2 degrees of freedom (efficiency and outlet pressure)
    m.fs.P1.efficiency_pump.fix(0.80)  # pump efficiency [-]
    operating_pressure = calculate_operating_pressure(
        feed_state_block=m.fs.feed.properties[0],
        over_pressure=over_pressure,
        water_recovery=water_recovery,
        NaCl_passage=0.01,
        solver=solver,
    )
    m.fs.P1.control_volume.properties_out[0].pressure.fix(operating_pressure)

    m.fs.P3.efficiency_pump.fix(0.80)  # pump efficiency [-]
    operating_pressure = calculate_operating_pressure(      # are all these operating conditions the same as P1?
        feed_state_block=m.fs.feed.properties[0],
        over_pressure=over_pressure,
        water_recovery=water_recovery,
        NaCl_passage=0.01,
        solver=solver,
    )
    m.fs.P3.control_volume.properties_out[0].pressure.fix(operating_pressure)

    # pressure exchanger
    m.fs.PXR1.efficiency_pressure_exchanger.fix(
        0.95
    )  # pressure exchanger efficiency [-]

    # pump 2, booster pumps, 1 degree of freedom (efficiency, pressure must match high pressure pump)
    m.fs.P2.efficiency_pump.fix(0.80)

    # mixer, no degrees of freedom

    # RO units
    m.fs.RO1.A_comp.fix(4.2e-12)  # membrane water permeability coefficient [m/s-Pa]
    m.fs.RO1.B_comp.fix(3.5e-8)  # membrane salt permeability coefficient [m/s]
    m.fs.RO1.channel_height.fix(1e-3)  # channel height in membrane stage [m]
    m.fs.RO1.spacer_porosity.fix(0.97)  # spacer porosity in membrane stage [-]
    m.fs.RO1.permeate.pressure[0].fix(101325)  # atmospheric pressure [Pa]
    m.fs.RO1.width.fix(5)  # stage width [m]

    # Assume these are the same as RO1?
    m.fs.RO2.A_comp.fix(4.2e-12)  # membrane water permeability coefficient [m/s-Pa]
    m.fs.RO2.B_comp.fix(3.5e-8)  # membrane salt permeability coefficient [m/s]
    m.fs.RO2.channel_height.fix(1e-3)  # channel height in membrane stage [m]
    m.fs.RO2.spacer_porosity.fix(0.97)  # spacer porosity in membrane stage [-]
    m.fs.RO2.permeate.pressure[0].fix(101325)  # atmospheric pressure [Pa]
    m.fs.RO2.width.fix(5)  # stage width [m]

    # initialize RO
    m.fs.RO1.feed_side.properties_in[0].flow_mass_phase_comp["Liq", "H2O"] = value(
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]
    )
    m.fs.RO1.feed_side.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"] = value(
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"]
    )
    m.fs.RO1.feed_side.properties_in[0].temperature = value(
        m.fs.feed.properties[0].temperature
    )
    m.fs.RO1.feed_side.properties_in[0].pressure = value(
        m.fs.P1.control_volume.properties_out[0].pressure
    )
    m.fs.RO1.area.fix(50)  # guess area for RO initialization
    m.fs.RO1.initialize(optarg=solver.options)

    # initialize RO2
    m.fs.RO2.feed_side.properties_in[0].flow_mass_phase_comp["Liq", "H2O"] = value(
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]
    )
    m.fs.RO2.feed_side.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"] = value(
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"]
    )
    m.fs.RO2.feed_side.properties_in[0].temperature = value(
        m.fs.feed.properties[0].temperature
    )
    m.fs.RO2.feed_side.properties_in[0].pressure = value(
        m.fs.P3.control_volume.properties_out[0].pressure
    )
    m.fs.RO2.area.fix(50)  # guess area for RO initialization
    m.fs.RO2.initialize(optarg=solver.options)


    # unfix guessed area, and fix water recovery - should water recovery be 0.5 for both stages?
    m.fs.RO1.area.unfix()
    m.fs.RO1.recovery_mass_phase_comp[0, "Liq", "H2O"].fix(water_recovery)

    m.fs.RO2.area.unfix()
    m.fs.RO2.recovery_mass_phase_comp[0, "Liq", "H2O"].fix(water_recovery)

    # check degrees of freedom
    if degrees_of_freedom(m) != 0:
        raise RuntimeError(
            "The set_operating_conditions function resulted in {} "
            "degrees of freedom rather than 0. This error suggests "
            "that too many or not enough variables are fixed for a "
            "simulation.".format(degrees_of_freedom(m))
        )


def calculate_operating_pressure(
    feed_state_block=None,
    over_pressure=0.15,
    water_recovery=0.5,
    NaCl_passage=0.01,
    solver=None,
):
    """
    estimate operating pressure for RO unit model given the following arguments:

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
    t.brine = prop.build_state_block([0], default={})

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


def solve(blk, solver=None, tee=False, check_termination=True):
    if solver is None:
        solver = get_solver()
    results = solver.solve(blk, tee=tee)
    if check_termination:
        assert_optimal_termination(results)
    return results


def initialize_system(m, solver=None):
    if solver is None:
        solver = get_solver()
    optarg = solver.options

    # ---initialize ROs---
    m.fs.RO1.initialize(optarg=optarg)
    m.fs.RO2.initialize(optarg=optarg)


    # ---initialize feed block---
    m.fs.feed.initialize(optarg=optarg)

    # ---initialize splitter and pressure exchanger---
    # pressure exchanger high pressure inlets
    propagate_state(m.fs.s06)  # propagate to PXR high pressure inlet from RO retentate
    m.fs.PXR1.high_pressure_side.properties_in.initialize(optarg=optarg)

    # splitter inlets
    propagate_state(m.fs.s01)  # propagate to splitter inlet from feed
    m.fs.S1.mixed_state.initialize(
        optarg=optarg
    )  # initialize inlet state block to solve for mass fraction
    # splitter outlet to PXR, enforce same volumetric flow as PXR high pressure inlet
    m.fs.S1.PXR1_state.calculate_state(
        var_args={
            (
                "flow_vol_phase",
                "Liq",
            ): value(  # same volumetric flow rate as PXR high pressure inlet
                m.fs.PXR1.high_pressure_side.properties_in[0].flow_vol_phase["Liq"]
            ),
            ("mass_frac_phase_comp", ("Liq", "NaCl")): value(
                m.fs.S1.mixed_state[0].mass_frac_phase_comp["Liq", "NaCl"]
            ),  # same as splitter inlet
            ("pressure", None): value(
                m.fs.S1.mixed_state[0].pressure
            ),  # same as splitter inlet
            ("temperature", None): value(m.fs.S1.mixed_state[0].temperature),
        },  # same as splitter inlet
    )
    # splitter initialization
    m.fs.S1.PXR1_state[0].flow_mass_phase_comp[
        "Liq", "NaCl"
    ].fix()  # fix the single degree of freedom for unit
    m.fs.S1.initialize(optarg=optarg)
    m.fs.S1.PXR1_state[0].flow_mass_phase_comp[
        "Liq", "NaCl"
    ].unfix()  # unfix for flowsheet simulation and optimization

    # pressure exchanger low pressure inlet
    propagate_state(m.fs.s08)

    # pressure exchanger initialization
    m.fs.PXR1.initialize(optarg=optarg)

    # ---initialize pumps 1 & 3---
    propagate_state(m.fs.s02)
    m.fs.P1.initialize(optarg=optarg)

    propagate_state(m.fs.s05)
    m.fs.P3.initialize(optarg=optarg)

    # ---initialize pump 2---
    propagate_state(m.fs.s09)
    m.fs.P2.control_volume.properties_out[0].pressure.fix(
        value(m.fs.P2.control_volume.properties_out[0].pressure)
    )
    m.fs.P2.initialize(optarg=optarg)
    m.fs.P2.control_volume.properties_out[0].pressure.unfix()

    # ---initialize mixers---
    propagate_state(m.fs.s03)
    propagate_state(m.fs.s10)
    m.fs.M1.initialize(optarg=optarg, outlvl=idaeslog.INFO)

    m.fs.costing.initialize()

def optimize_set_up(m):
    # objective
    m.fs.objective = Objective(expr=m.fs.costing.LCOW)

    # unfix decision variables and add bounds
    # pumps
    m.fs.P1.control_volume.properties_out[0].pressure.unfix()
    m.fs.P1.control_volume.properties_out[0].pressure.setlb(10e5)
    m.fs.P1.control_volume.properties_out[0].pressure.setub(80e5)
    m.fs.P1.deltaP.setlb(0)
    m.fs.P2.control_volume.properties_out[0].pressure.setlb(10e5)
    m.fs.P2.control_volume.properties_out[0].pressure.setub(80e5)
    m.fs.P2.deltaP.setlb(0)
    m.fs.P3.control_volume.properties_out[0].pressure.unfix()
    m.fs.P3.control_volume.properties_out[0].pressure.setlb(10e5)
    m.fs.P3.control_volume.properties_out[0].pressure.setub(80e5)
    m.fs.P3.deltaP.setlb(0)

    # RO units
    m.fs.RO1.area.setlb(1)
    m.fs.RO1.area.setub(150)
    m.fs.RO2.area.setlb(1)
    m.fs.RO2.area.setub(150)

    # additional specifications
    m.fs.product_salinity = Param(
        initialize=500e-6, mutable=True
    )  # product NaCl mass fraction [-]
    m.fs.minimum_water_flux = Param(
        initialize=1.0 / 3600.0, mutable=True
    )  # minimum water flux [kg/m2-s]

    # additional constraints
    m.fs.eq_product_quality = Constraint(
        expr=m.fs.product.properties[0].mass_frac_phase_comp["Liq", "NaCl"]
        <= m.fs.product_salinity
    )
    iscale.constraint_scaling_transform(
        m.fs.eq_product_quality, 1e3
    )  # scaling constraint
    m.fs.eq_minimum_water_flux1 = Constraint(
        expr=m.fs.RO1.flux_mass_phase_comp[0, 1, "Liq", "H2O"] >= m.fs.minimum_water_flux
    )
    m.fs.eq_minimum_water_flux2 = Constraint(
        expr=m.fs.RO2.flux_mass_phase_comp[0, 1, "Liq", "H2O"] >= m.fs.minimum_water_flux
    )

    # ---checking model---
    assert_degrees_of_freedom(m, 1)


def optimize(m, solver=None, check_termination=True):
    # --solve---
    return solve(m, solver=solver, check_termination=check_termination)


def display_system(m):
    print("---system metrics---")
    feed_flow_mass = sum(
        m.fs.feed.flow_mass_phase_comp[0, "Liq", j].value for j in ["H2O", "NaCl"]
    )
    feed_mass_frac_NaCl = (
        m.fs.feed.flow_mass_phase_comp[0, "Liq", "NaCl"].value / feed_flow_mass
    )
    print("Feed: %.2f kg/s, %.0f ppm" % (feed_flow_mass, feed_mass_frac_NaCl * 1e6))

    prod_flow_mass = sum(
        m.fs.product.flow_mass_phase_comp[0, "Liq", j].value for j in ["H2O", "NaCl"]
    )
    prod_mass_frac_NaCl = (
        m.fs.product.flow_mass_phase_comp[0, "Liq", "NaCl"].value / prod_flow_mass
    )
    print("Product: %.3f kg/s, %.0f ppm" % (prod_flow_mass, prod_mass_frac_NaCl * 1e6))

    print(
        "Volumetric recovery: %.1f%%"
        % (value(m.fs.RO2.recovery_vol_phase[0, "Liq"]) * 100)
    )
    print(
        "Water recovery: %.1f%%"
        % (value(m.fs.RO2.recovery_mass_phase_comp[0, "Liq", "H2O"]) * 100)
    )
    print(
        "Energy Consumption: %.1f kWh/m3"
        % value(m.fs.costing.specific_energy_consumption)
    )
    print("Levelized cost of water: %.2f $/m3" % value(m.fs.costing.LCOW))


def display_design(m):
    print("---decision variables---")
    print("RO operating pressure %.1f bar" % (m.fs.RO1.inlet.pressure[0].value / 1e5))
    print("RO membrane area %.1f m2" % (m.fs.RO1.area.value))

    print("---design variables---")
    print("Separator 1")
    print("Split fraction %.2f" % (m.fs.S1.split_fraction[0, "PXR1"].value * 100))
    print(
        "Pump 1\noutlet pressure: %.1f bar\npower %.2f kW"
        % (
            m.fs.P1.outlet.pressure[0].value / 1e5,
            m.fs.P1.work_mechanical[0].value / 1e3,
        )
    )
    print(
        "Pump 2\noutlet pressure: %.1f bar\npower %.2f kW"
        % (
            m.fs.P2.outlet.pressure[0].value / 1e5,
            m.fs.P2.work_mechanical[0].value / 1e3,
        )
    )
    print(
        "Pump 3\noutlet pressure: %.1f bar\npower %.2f kW"
        % (
            m.fs.P3.outlet.pressure[0].value / 1e5,
            m.fs.P3.work_mechanical[0].value / 1e3,
        )
    )

def display_state(m):
    print("---state---")

    def print_state(s, b):
        flow_mass = sum(
            b.flow_mass_phase_comp[0, "Liq", j].value for j in ["H2O", "NaCl"]
        )
        mass_frac_ppm = b.flow_mass_phase_comp[0, "Liq", "NaCl"].value / flow_mass * 1e6
        pressure_bar = b.pressure[0].value / 1e5
        print(
            s
            + ": %.3f kg/s, %.0f ppm, %.1f bar"
            % (flow_mass, mass_frac_ppm, pressure_bar)
        )
    print("--1st stage--")
    print_state("Feed      ", m.fs.feed.outlet)
    print_state("Split 1   ", m.fs.S1.P1)
    print_state("P1 out    ", m.fs.P1.outlet)
    print_state("Split 2   ", m.fs.S1.PXR1)
    print_state("PXR LP out", m.fs.PXR1.low_pressure_outlet)
    print_state("P2 out    ", m.fs.P2.outlet)
    print_state("Mix out   ", m.fs.M1.outlet)
    print_state("RO1 perm   ", m.fs.RO1.permeate)
    print_state("RO1 reten  ", m.fs.RO1.retentate)
    print_state("PXR HP out", m.fs.PXR1.high_pressure_outlet)

    print("--2nd stage--")
    print_state("P3 out    ", m.fs.P3.outlet)
    print_state("RO2 perm   ", m.fs.RO2.permeate)
    print_state("RO2 reten  ", m.fs.RO2.retentate)

if __name__ == "__main__":
    main()
