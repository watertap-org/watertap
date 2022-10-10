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
import os
from pyomo.environ import (
    ConcreteModel,
    value,
    Constraint,
    Expression,
    Objective,
    Var,
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
from idaes.core.util.misc import StrEnum

import watertap.property_models.NaCl_prop_pack as props
from watertap.unit_models.reverse_osmosis_0D import (
    ReverseOsmosis0D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
from watertap.unit_models.pressure_exchanger import PressureExchanger
from watertap.unit_models.pressure_changer import Pump, EnergyRecoveryDevice
from watertap.core.util.initialization import assert_degrees_of_freedom
from watertap.costing import WaterTAPCosting


class ERDtype(StrEnum):
    pressure_exchanger = "pressure_exchanger"
    pump_as_turbine = "pump_as_turbine"


def erd_type_not_found(erd_type):
    raise NotImplementedError(
        "erd_type was {}, but can only "
        "be pressure_exchanger or pump_as_turbine"
        "".format(erd_type.value)
    )


def main(erd_type=ERDtype.pressure_exchanger):
    # set up solver
    solver = get_solver()

    # build, set, and initialize
    m = build(erd_type=erd_type)
    set_operating_conditions(m)
    initialize_system(m, solver=solver)

    # optimize and display
    optimize_set_up(m)
    solve(m, solver=solver)

    print("\n***---Simulation results---***")
    display_system(m)
    display_design(m)
    if erd_type == ERDtype.pressure_exchanger:
        display_state(m)
    else:
        pass

    return m


def build(erd_type=ERDtype.pressure_exchanger):
    # flowsheet set up
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.erd_type = erd_type
    m.fs.properties = props.NaClParameterBlock()
    m.fs.costing = WaterTAPCosting()

    # Control volume flow blocks
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)
    m.fs.disposal = Product(property_package=m.fs.properties)

    # --- Main pump ---
    m.fs.P1 = Pump(property_package=m.fs.properties)
    m.fs.P1.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    # --- Reverse Osmosis Block ---
    m.fs.RO = ReverseOsmosis0D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
    )
    m.fs.RO.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    # --- ERD blocks ---
    if erd_type == ERDtype.pressure_exchanger:
        m.fs.S1 = Separator(property_package=m.fs.properties, outlet_list=["P1", "PXR"])

        m.fs.PXR = PressureExchanger(property_package=m.fs.properties)
        m.fs.P2 = Pump(property_package=m.fs.properties)
        m.fs.M1 = Mixer(
            property_package=m.fs.properties,
            momentum_mixing_type=MomentumMixingType.equality,
            inlet_list=["P1", "P2"],
        )

        # add costing for PX and recirculation pump
        m.fs.PXR.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.P2.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        # mixer and separator have no associated costing
    elif erd_type == ERDtype.pump_as_turbine:
        # add energy recovery turbine block
        m.fs.ERD = EnergyRecoveryDevice(property_package=m.fs.properties)
        # add costing for ERD config
        m.fs.ERD.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    else:
        erd_type_not_found(erd_type)

    # process costing and add system level metrics
    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(m.fs.product.properties[0].flow_vol)
    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(m.fs.product.properties[0].flow_vol)
    m.fs.costing.add_specific_electrical_carbon_intensity(
        m.fs.product.properties[0].flow_vol
    )

    # connections
    if erd_type == ERDtype.pressure_exchanger:
        m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.S1.inlet)
        m.fs.s02 = Arc(source=m.fs.S1.P1, destination=m.fs.P1.inlet)
        m.fs.s03 = Arc(source=m.fs.P1.outlet, destination=m.fs.M1.P1)
        m.fs.s04 = Arc(source=m.fs.M1.outlet, destination=m.fs.RO.inlet)
        m.fs.s05 = Arc(source=m.fs.RO.permeate, destination=m.fs.product.inlet)
        m.fs.s06 = Arc(
            source=m.fs.RO.retentate, destination=m.fs.PXR.high_pressure_inlet
        )
        m.fs.s07 = Arc(
            source=m.fs.PXR.high_pressure_outlet, destination=m.fs.disposal.inlet
        )
        m.fs.s08 = Arc(source=m.fs.S1.PXR, destination=m.fs.PXR.low_pressure_inlet)
        m.fs.s09 = Arc(source=m.fs.PXR.low_pressure_outlet, destination=m.fs.P2.inlet)
        m.fs.s10 = Arc(source=m.fs.P2.outlet, destination=m.fs.M1.P2)
    elif erd_type == ERDtype.pump_as_turbine:
        m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.P1.inlet)
        m.fs.s02 = Arc(source=m.fs.P1.outlet, destination=m.fs.RO.inlet)
        m.fs.s03 = Arc(source=m.fs.RO.permeate, destination=m.fs.product.inlet)
        m.fs.s04 = Arc(source=m.fs.RO.retentate, destination=m.fs.ERD.inlet)
        m.fs.s05 = Arc(source=m.fs.ERD.outlet, destination=m.fs.disposal.inlet)
    else:
        # this case should be caught in the previous conditional
        erd_type_not_found(erd_type)

    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling
    # set default property values
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )
    # set unit model values
    iscale.set_scaling_factor(m.fs.P1.control_volume.work, 1e-3)
    iscale.set_scaling_factor(m.fs.RO.area, 1e-2)
    m.fs.feed.properties[0].flow_vol_phase["Liq"]
    m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"]
    if erd_type == ERDtype.pressure_exchanger:
        iscale.set_scaling_factor(m.fs.P2.control_volume.work, 1e-3)
        iscale.set_scaling_factor(m.fs.PXR.low_pressure_side.work, 1e-3)
        iscale.set_scaling_factor(m.fs.PXR.high_pressure_side.work, 1e-3)
        # touch properties used in specifying and initializing the model
        m.fs.S1.mixed_state[0].mass_frac_phase_comp
        m.fs.S1.PXR_state[0].flow_vol_phase["Liq"]
    elif erd_type == ERDtype.pump_as_turbine:
        iscale.set_scaling_factor(m.fs.ERD.control_volume.work, 1e-3)
    else:
        erd_type_not_found(erd_type)
    # unused scaling factors needed by IDAES base costing module
    # calculate and propagate scaling factors
    iscale.calculate_scaling_factors(m)

    return m


def set_operating_conditions(
    m,
    water_recovery=0.5,
    over_pressure=0.3,
    solver=None,
):

    if solver is None:
        solver = get_solver()
    # ---specifications---
    # feed
    # state variables
    m.fs.feed.properties[0].pressure.fix(101325)  # feed pressure [Pa]
    m.fs.feed.properties[0].temperature.fix(273.15 + 25)  # feed temperature [K]
    # properties (cannot be fixed for initialization routines, must calculate the state variables)
    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): 1e-3,  # feed volumetric flow rate [m3/s]
            ("mass_frac_phase_comp", ("Liq", "NaCl")): 0.035,
        },  # feed NaCl mass fraction [-]
        hold_state=True,  # fixes the calculated component mass flow rates
    )

    # pump 1, high pressure pump, 2 degrees of freedom (efficiency and outlet pressure)
    m.fs.P1.efficiency_pump.fix(0.80)  # pump efficiency [-]
    operating_pressure = calculate_operating_pressure(
        feed_state_block=m.fs.feed.properties[0],
        over_pressure=over_pressure,
        water_recovery=water_recovery,
        NaCl_passage=0.01,
        solver=solver,
    )
    m.fs.P1.control_volume.properties_out[0].pressure.fix(operating_pressure)

    # RO unit
    m.fs.RO.A_comp.fix(4.2e-12)  # membrane water permeability coefficient [m/s-Pa]
    m.fs.RO.B_comp.fix(3.5e-8)  # membrane salt permeability coefficient [m/s]
    m.fs.RO.channel_height.fix(1e-3)  # channel height in membrane stage [m]
    m.fs.RO.spacer_porosity.fix(0.97)  # spacer porosity in membrane stage [-]
    m.fs.RO.permeate.pressure[0].fix(101325)  # atmospheric pressure [Pa]
    m.fs.RO.width.fix(5)  # stage width [m]
    # initialize RO
    m.fs.RO.feed_side.properties_in[0].flow_mass_phase_comp["Liq", "H2O"] = value(
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]
    )
    m.fs.RO.feed_side.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"] = value(
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"]
    )
    m.fs.RO.feed_side.properties_in[0].temperature = value(
        m.fs.feed.properties[0].temperature
    )
    m.fs.RO.feed_side.properties_in[0].pressure = value(
        m.fs.P1.control_volume.properties_out[0].pressure
    )

    m.fs.RO.area.fix(50)  # guess area for RO initialization

    if m.fs.erd_type == ERDtype.pressure_exchanger:
        # pressure exchanger
        m.fs.PXR.efficiency_pressure_exchanger.fix(
            0.95
        )  # pressure exchanger efficiency [-]
        # pump 2, booster pump, 1 degree of freedom (efficiency, pressure must match high pressure pump)
        m.fs.P2.efficiency_pump.fix(0.80)

    elif m.fs.erd_type == ERDtype.pump_as_turbine:
        # energy recovery turbine - efficiency and outlet pressure
        m.fs.ERD.efficiency_pump.fix(0.95)
        m.fs.ERD.control_volume.properties_out[0].pressure.fix(101325)
    else:
        erd_type_not_found(m.fs.erd_type)

    m.fs.RO.initialize(optarg=solver.options)

    # unfix guessed area, and fix water recovery
    m.fs.RO.area.unfix()

    m.fs.RO.recovery_mass_phase_comp[0, "Liq", "H2O"].fix(water_recovery)

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

    # ---initialize RO---
    m.fs.RO.initialize(optarg=optarg)

    # ---initialize feed block---
    m.fs.feed.initialize(optarg=optarg)

    # --- initialize ERD ---
    if m.fs.erd_type == ERDtype.pressure_exchanger:
        initialize_pressure_exchanger(m, optarg)

    elif m.fs.erd_type == ERDtype.pump_as_turbine:
        initialize_pump_as_turbine(m, optarg)

    else:
        erd_type_not_found(m.fs.erd_type)

    m.fs.costing.initialize()


def initialize_pressure_exchanger(m, optarg):
    # ---initialize splitter and pressure exchanger---
    # pressure exchanger high pressure inlet
    propagate_state(m.fs.s06)  # propagate to PXR high pressure inlet from RO retentate
    m.fs.PXR.high_pressure_side.properties_in.initialize(optarg=optarg)

    # splitter inlet
    propagate_state(m.fs.s01)  # propagate to splitter inlet from feed
    m.fs.S1.mixed_state.initialize(
        optarg=optarg
    )  # initialize inlet state block to solve for mass fraction
    # splitter outlet to PXR, enforce same volumetric flow as PXR high pressure inlet
    m.fs.S1.PXR_state.calculate_state(
        var_args={
            (
                "flow_vol_phase",
                "Liq",
            ): value(  # same volumetric flow rate as PXR high pressure inlet
                m.fs.PXR.high_pressure_side.properties_in[0].flow_vol_phase["Liq"]
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
    m.fs.S1.PXR_state[0].flow_mass_phase_comp[
        "Liq", "NaCl"
    ].fix()  # fix the single degree of freedom for unit
    m.fs.S1.initialize(optarg=optarg)
    m.fs.S1.PXR_state[0].flow_mass_phase_comp[
        "Liq", "NaCl"
    ].unfix()  # unfix for flowsheet simulation and optimization

    # pressure exchanger low pressure inlet
    propagate_state(m.fs.s08)

    # pressure exchanger initialization
    m.fs.PXR.initialize(optarg=optarg)

    # ---initialize pump 1---
    propagate_state(m.fs.s02)
    m.fs.P1.initialize(optarg=optarg)

    # ---initialize pump 2---
    propagate_state(m.fs.s09)
    m.fs.P2.control_volume.properties_out[0].pressure.fix(
        value(m.fs.P2.control_volume.properties_out[0].pressure)
    )
    m.fs.P2.initialize(optarg=optarg)
    m.fs.P2.control_volume.properties_out[0].pressure.unfix()

    # ---initialize mixer---
    propagate_state(m.fs.s03)
    propagate_state(m.fs.s10)
    m.fs.M1.initialize(optarg=optarg, outlvl=idaeslog.INFO)


def initialize_pump_as_turbine(m, optarg):
    propagate_state(m.fs.s05)
    m.fs.ERD.initialize(optarg=optarg)
    propagate_state(m.fs.s01)
    m.fs.P1.initialize(optarg=optarg)
    propagate_state(m.fs.s02)


def optimize_set_up(m):
    # add objective
    m.fs.objective = Objective(expr=m.fs.costing.LCOW)

    # unfix decision variables and add bounds
    # pump 1 and pump 2
    m.fs.P1.control_volume.properties_out[0].pressure.unfix()
    m.fs.P1.control_volume.properties_out[0].pressure.setlb(10e5)
    m.fs.P1.control_volume.properties_out[0].pressure.setub(80e5)
    m.fs.P1.deltaP.setlb(0)
    if m.fs.erd_type == ERDtype.pressure_exchanger:
        m.fs.P2.control_volume.properties_out[0].pressure.setlb(10e5)
        m.fs.P2.control_volume.properties_out[0].pressure.setub(80e5)
        m.fs.P2.deltaP.setlb(0)
    else:
        # no additional optimization needed for pump_as_turbine configuration
        pass

    # RO
    m.fs.RO.area.unfix()
    m.fs.RO.area.setlb(1)
    m.fs.RO.area.setub(150)

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
    m.fs.eq_minimum_water_flux = Constraint(
        expr=m.fs.RO.flux_mass_phase_comp[0, 1, "Liq", "H2O"] >= m.fs.minimum_water_flux
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
        % (value(m.fs.RO.recovery_vol_phase[0, "Liq"]) * 100)
    )
    print(
        "Water recovery: %.1f%%"
        % (value(m.fs.RO.recovery_mass_phase_comp[0, "Liq", "H2O"]) * 100)
    )
    print(
        "Energy Consumption: %.1f kWh/m3"
        % value(m.fs.costing.specific_energy_consumption)
    )
    print("Levelized cost of water: %.2f $/m3" % value(m.fs.costing.LCOW))


def display_design(m):
    print("---decision variables---")
    print("Operating pressure %.1f bar" % (m.fs.RO.inlet.pressure[0].value / 1e5))
    print("Membrane area %.1f m2" % (m.fs.RO.area.value))

    print("---design variables---")
    print(
        "Pump 1\noutlet pressure: %.1f bar\npower %.2f kW"
        % (
            m.fs.P1.outlet.pressure[0].value / 1e5,
            m.fs.P1.work_mechanical[0].value / 1e3,
        )
    )
    if m.fs.erd_type == ERDtype.pressure_exchanger:
        print("Separator")
        print("Split fraction %.2f" % (m.fs.S1.split_fraction[0, "PXR"].value * 100))
        print(
            "Pump 2\noutlet pressure: %.1f bar\npower %.2f kW"
            % (
                m.fs.P2.outlet.pressure[0].value / 1e5,
                m.fs.P2.work_mechanical[0].value / 1e3,
            )
        )
    elif m.fs.erd_type == ERDtype.pump_as_turbine:
        print(
            "ERD\ninlet pressure: %.1f bar\npower recovered %.2f kW"
            % (
                m.fs.ERD.inlet.pressure[0].value / 1e5,
                -1 * m.fs.ERD.work_mechanical[0].value / 1e3,
            )
        )
    else:
        erd_type_not_found(m.fs.erd_type)


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

    print_state("Feed      ", m.fs.feed.outlet)
    print_state("Split 1   ", m.fs.S1.P1)
    print_state("P1 out    ", m.fs.P1.outlet)
    print_state("Split 2   ", m.fs.S1.PXR)
    print_state("PXR LP out", m.fs.PXR.low_pressure_outlet)
    print_state("P2 out    ", m.fs.P2.outlet)
    print_state("Mix out   ", m.fs.M1.outlet)
    print_state("RO perm   ", m.fs.RO.permeate)
    print_state("RO reten  ", m.fs.RO.retentate)
    print_state("PXR HP out", m.fs.PXR.high_pressure_outlet)


if __name__ == "__main__":
    # m = main(erd_type=ERDtype.pressure_exchanger)
    m = main(erd_type=ERDtype.pump_as_turbine)
