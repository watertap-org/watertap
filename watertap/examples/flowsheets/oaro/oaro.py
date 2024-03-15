#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

from pyomo.environ import (
    ConcreteModel,
    value,
    Constraint,
    Objective,
    Var,
    Param,
    NonNegativeReals,
    TransformationFactory,
    units as pyunits,
    check_optimal_termination,
)
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import (
    propagate_state,
)
from idaes.models.unit_models import Product, Feed
from idaes.core import UnitModelCostingBlock
import idaes.core.util.scaling as iscale
from idaes.core.util.misc import StrEnum

import watertap.property_models.NaCl_prop_pack as props
from watertap.unit_models.reverse_osmosis_0D import (
    ReverseOsmosis0D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
from watertap.unit_models.osmotically_assisted_reverse_osmosis_0D import (
    OsmoticallyAssistedReverseOsmosis0D,
)
from watertap.unit_models.pressure_changer import Pump, EnergyRecoveryDevice
from watertap.core.util.initialization import assert_degrees_of_freedom
from watertap.costing import WaterTAPCosting


class ERDtype(StrEnum):
    pump_as_turbine = "pump_as_turbine"


def erd_type_not_found(erd_type):
    raise NotImplementedError(
        "erd_type was {}, but can only " "be pump_as_turbine" "".format(erd_type.value)
    )


def main(erd_type=ERDtype.pump_as_turbine, raise_on_failure=False):
    # set up solver
    solver = get_solver()

    # build, set, and initialize
    m = build(erd_type=erd_type)
    set_operating_conditions(m)
    initialize_system(m, solver=solver)

    optimize_set_up(m)
    solve(m, solver=solver)

    print("\n***---Simulation results---***")
    display_system(m)
    if erd_type == ERDtype.pump_as_turbine:
        display_state(m)
    else:
        pass

    return m


def build(erd_type=ERDtype.pump_as_turbine):
    # TODO: add costing later (OARO unit model does not have costing method)

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

    m.fs.P2 = Pump(property_package=m.fs.properties)
    m.fs.P2.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    m.fs.P3 = Pump(property_package=m.fs.properties)
    m.fs.P3.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    # --- Reverse Osmosis Block ---
    m.fs.RO = ReverseOsmosis0D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        has_full_reporting=True,
    )
    m.fs.RO.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    # --- Osmotically Assisted Reverse Osmosis Block ---
    m.fs.OARO = OsmoticallyAssistedReverseOsmosis0D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        has_full_reporting=True,
    )
    m.fs.OARO.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    # --- ERD blocks ---
    if erd_type == ERDtype.pump_as_turbine:
        # add energy recovery turbine block
        m.fs.ERD1 = EnergyRecoveryDevice(property_package=m.fs.properties)
        m.fs.ERD2 = EnergyRecoveryDevice(property_package=m.fs.properties)
        # add costing for ERD config
        m.fs.ERD1.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.ERD2.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
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

    # system water recovery
    m.fs.volumetric_recovery = Var(
        initialize=0.5,
        bounds=(0, 1),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="System Volumetric Recovery of Water",
    )
    m.fs.eq_volumetric_recovery = Constraint(
        expr=m.fs.feed.properties[0].flow_vol_phase["Liq"] * m.fs.volumetric_recovery
        == m.fs.product.properties[0].flow_vol_phase["Liq"]
    )

    m.fs.water_recovery = Var(
        initialize=0.5,
        bounds=(0, 1),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="System Water Recovery",
    )
    m.fs.eq_water_recovery = Constraint(
        expr=m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]
        * m.fs.water_recovery
        == m.fs.product.properties[0].flow_mass_phase_comp["Liq", "H2O"]
    )

    # connections
    if erd_type == ERDtype.pump_as_turbine:
        m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.P1.inlet)
        m.fs.s02 = Arc(source=m.fs.P1.outlet, destination=m.fs.OARO.feed_inlet)
        m.fs.s03 = Arc(source=m.fs.OARO.feed_outlet, destination=m.fs.ERD1.inlet)
        m.fs.s04 = Arc(source=m.fs.ERD1.outlet, destination=m.fs.disposal.inlet)
        m.fs.s05 = Arc(source=m.fs.OARO.permeate_outlet, destination=m.fs.P2.inlet)
        m.fs.s06 = Arc(source=m.fs.P2.outlet, destination=m.fs.RO.inlet)
        m.fs.s07 = Arc(source=m.fs.RO.permeate, destination=m.fs.product.inlet)
        m.fs.s08 = Arc(source=m.fs.RO.retentate, destination=m.fs.ERD2.inlet)
        m.fs.s09 = Arc(source=m.fs.ERD2.outlet, destination=m.fs.P3.inlet)
        m.fs.s10 = Arc(source=m.fs.P3.outlet, destination=m.fs.OARO.permeate_inlet)

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
    iscale.set_scaling_factor(m.fs.P2.control_volume.work, 1e-3)
    iscale.set_scaling_factor(m.fs.P3.control_volume.work, 1e-3)
    iscale.set_scaling_factor(m.fs.RO.area, 1e-2)
    iscale.set_scaling_factor(m.fs.OARO.area, 1e-2)
    m.fs.feed.properties[0].flow_vol_phase["Liq"]
    m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"]
    if erd_type == ERDtype.pump_as_turbine:
        iscale.set_scaling_factor(m.fs.ERD1.control_volume.work, 1e-3)
        iscale.set_scaling_factor(m.fs.ERD2.control_volume.work, 1e-3)
    else:
        erd_type_not_found(erd_type)
    # unused scaling factors needed by IDAES base costing module
    # calculate and propagate scaling factors
    iscale.calculate_scaling_factors(m)

    return m


def set_operating_conditions(
    m,
    solver=None,
):

    if solver is None:
        solver = get_solver()

    # ---specifications---
    # feed
    # state variables
    pressure_atmospheric = 101325
    feed_temperature = 273.15 + 25
    m.fs.feed.properties[0].pressure.fix(pressure_atmospheric)  # feed pressure [Pa]
    m.fs.feed.properties[0].temperature.fix(feed_temperature)  # feed temperature [K]

    # properties (cannot be fixed for initialization routines, must calculate the state variables)
    feed_flow_mass = 1
    feed_mass_frac_NaCl = 0.035
    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )

    # pump 1, high pressure pump, 2 degrees of freedom (efficiency and outlet pressure)
    m.fs.P1.efficiency_pump.fix(0.80)  # pump efficiency [-]
    m.fs.P1.control_volume.properties_out[0].pressure.fix(50e5)

    # pump 2, 2 degrees of freedom (efficiency and outlet pressure)
    m.fs.P2.efficiency_pump.fix(0.80)  # pump efficiency [-]
    m.fs.P2.control_volume.properties_out[0].pressure.fix(18e5)

    # pump 3, 2 degrees of freedom (efficiency, outlet pressure)
    m.fs.P3.efficiency_pump.fix(0.80)  # pump efficiency [-]
    m.fs.P3.control_volume.properties_out[0].pressure.fix(4e5)

    # initial guess for states of pump 3 output (temperature and concentrations)
    m.fs.P3.control_volume.properties_out[0].temperature.value = feed_temperature
    permeate_flow_mass = 0.75
    permeate_mass_frac_NaCl = 0.015
    permeate_mass_frac_H2O = 1 - permeate_mass_frac_NaCl
    m.fs.P3.control_volume.properties_out[0].flow_mass_phase_comp[
        "Liq", "H2O"
    ].value = (permeate_flow_mass * permeate_mass_frac_H2O)
    m.fs.P3.control_volume.properties_out[0].flow_mass_phase_comp[
        "Liq", "NaCl"
    ].value = (permeate_flow_mass * permeate_mass_frac_NaCl)

    # Initialize OARO
    membrane_area = 50
    A = 4.2e-12
    B = 1.3e-8

    m.fs.OARO.area.fix(membrane_area)

    m.fs.OARO.A_comp.fix(A)
    m.fs.OARO.B_comp.fix(B)

    m.fs.OARO.structural_parameter.fix(300e-6)

    m.fs.OARO.permeate_side.channel_height.fix(0.001)
    m.fs.OARO.permeate_side.spacer_porosity.fix(0.75)
    m.fs.OARO.feed_side.channel_height.fix(0.002)
    m.fs.OARO.feed_side.spacer_porosity.fix(0.75)
    m.fs.OARO.feed_side.velocity[0, 0].fix(0.1)

    # RO unit
    m.fs.RO.A_comp.fix(4.2e-12)  # membrane water permeability coefficient [m/s-Pa]
    m.fs.RO.B_comp.fix(3.5e-8)  # membrane salt permeability coefficient [m/s]
    m.fs.RO.feed_side.channel_height.fix(1e-3)  # channel height in membrane stage [m]
    m.fs.RO.feed_side.spacer_porosity.fix(0.97)  # spacer porosity in membrane stage [-]
    m.fs.RO.permeate.pressure[0].fix(101325)  # atmospheric pressure [Pa]
    m.fs.RO.width.fix(5)  # stage width [m]
    m.fs.RO.area.fix(50)  # guess area for RO initialization

    if m.fs.erd_type == ERDtype.pump_as_turbine:
        # energy recovery turbine - efficiency and outlet pressure
        m.fs.ERD1.efficiency_pump.fix(0.95)
        m.fs.ERD1.control_volume.properties_out[0].pressure.fix(pressure_atmospheric)

        m.fs.ERD2.efficiency_pump.fix(0.95)
        m.fs.ERD2.control_volume.properties_out[0].pressure.fix(pressure_atmospheric)
    else:
        erd_type_not_found(m.fs.erd_type)

    # check degrees of freedom
    if degrees_of_freedom(m) != 0:
        raise RuntimeError(
            "The set_operating_conditions function resulted in {} "
            "degrees of freedom rather than 0. This error suggests "
            "that too many or not enough variables are fixed for a "
            "simulation.".format(degrees_of_freedom(m))
        )


def solve(blk, solver=None, tee=True):
    if solver is None:
        solver = get_solver()
    results = solver.solve(blk, tee=tee)
    if not check_optimal_termination(results):
        results = solver.solve(blk, tee=tee)
    return results


def initialize_loop(m, solver):

    propagate_state(m.fs.s05)
    m.fs.P2.initialize()

    # ---initialize RO---
    propagate_state(m.fs.s06)
    m.fs.RO.initialize()
    propagate_state(m.fs.s07)

    propagate_state(m.fs.s08)
    m.fs.ERD2.initialize()

    propagate_state(m.fs.s09)
    m.fs.P3.initialize()

    propagate_state(m.fs.s10)
    propagate_state(m.fs.s02)
    m.fs.OARO.initialize()

    propagate_state(m.fs.s03)
    m.fs.ERD1.initialize()
    propagate_state(m.fs.s04)
    m.fs.disposal.initialize()


def initialize_system(m, solver=None, verbose=True):
    if solver is None:
        solver = get_solver()

    # ---initialize feed block---
    m.fs.feed.initialize()

    propagate_state(m.fs.s01)
    m.fs.P1.initialize()

    # ---initialize OARO---
    propagate_state(m.fs.s02)
    propagate_state(m.fs.s10)
    m.fs.OARO.initialize()

    initialize_loop(m, solver)

    print(f"DOF: {degrees_of_freedom(m)}")

    # Now that the units are initialized, we can fix the
    # permeate side outlet pressure and unfix the RO pump
    # (which allows for control over the flow mass composition
    # into the OARO permeate_side).
    m.fs.OARO.permeate_side.properties_out[0].pressure.fix(101325)
    m.fs.P2.control_volume.properties_out[0].pressure.unfix()

    print(f"DOF: {degrees_of_freedom(m)}")

    m.fs.costing.initialize()


def optimize_set_up(m):
    # add objective
    m.fs.objective = Objective(expr=m.fs.costing.LCOW)

    # unfix decision variables and add bounds
    # pump 1
    m.fs.P1.control_volume.properties_out[0].pressure.unfix()
    m.fs.P1.control_volume.properties_out[0].pressure.setlb(10e5)
    m.fs.P1.control_volume.properties_out[0].pressure.setub(80e5)
    m.fs.P1.deltaP.setlb(0)

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
    assert_degrees_of_freedom(m, 2)


def optimize(m, solver=None):
    # --solve---
    return solve(m, solver=solver)


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

    print("Volumetric recovery: %.1f%%" % (value(m.fs.volumetric_recovery) * 100))
    print("Water recovery: %.1f%%" % (value(m.fs.water_recovery) * 100))

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
    if m.fs.erd_type == ERDtype.pump_as_turbine:
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
    print_state("P1 out    ", m.fs.P1.outlet)
    print_state("OARO Feed in", m.fs.OARO.feed_inlet)
    print_state("OARO Feed out", m.fs.OARO.feed_outlet)
    print_state("OARO Perm in", m.fs.OARO.permeate_inlet)
    print_state("OARO Perm out", m.fs.OARO.permeate_outlet)
    print_state("ERD1 out    ", m.fs.ERD1.outlet)
    print_state("P2 out    ", m.fs.P2.outlet)
    print_state("RO reten  ", m.fs.RO.retentate)
    print_state("RO perm   ", m.fs.RO.permeate)
    print_state("ERD2 out    ", m.fs.ERD2.outlet)
    print_state("P3 out    ", m.fs.P3.outlet)


if __name__ == "__main__":
    m = main(erd_type=ERDtype.pump_as_turbine)
