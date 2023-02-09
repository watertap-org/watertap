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
import itertools
from pyomo.environ import (
    ConcreteModel,
    value,
    Param,
    Var,
    Constraint,
    Expression,
    Objective,
    TransformationFactory,
    Block,
    NonNegativeReals,
    RangeSet,
    Set,
    check_optimal_termination,
    units as pyunits,
)
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.exceptions import InitializationError
from idaes.core.util.model_statistics import degrees_of_freedom, fixed_variables_set
from idaes.core.util.initialization import (
    solve_indexed_blocks,
    propagate_state,
    # propagate_state as _pro_state,
)
from idaes.models.unit_models import Mixer, Separator, Product, Feed
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
from watertap.unit_models.osmotically_assisted_reverse_osmosis_0D import (
    OsmoticallyAssistedReverseOsmosis0D,
    MassTransferCoefficient,
    PressureChangeType,
)
from watertap.unit_models.pressure_changer import Pump, EnergyRecoveryDevice
from watertap.core.util.initialization import assert_degrees_of_freedom
from watertap.costing import WaterTAPCosting

# from watertap.core.util.infeasible import *


class ERDtype(StrEnum):
    pump_as_turbine = "pump_as_turbine"


def erd_type_not_found(erd_type):
    raise NotImplementedError(
        "erd_type was {}, but can only " "be pump_as_turbine" "".format(erd_type.value)
    )


# def propagate_state(arc):
#     _pro_state(arc)
#     print(arc.destination.name)
#     arc.destination.display()


def main(number_of_stages, erd_type=ERDtype.pump_as_turbine):
    # set up solver
    solver = get_solver()

    # build, set, and initialize
    m = build(number_of_stages=number_of_stages, erd_type=erd_type)
    set_operating_conditions(m)
    initialize_system(m, solver=solver)

    # print_close_to_bounds(m)
    # print_infeasible_constraints(m)

    # optimize_set_up(m)
    solve(m, solver=solver)

    print("\n***---Simulation results---***")
    display_system(m)
    display_design(m)
    if erd_type == ERDtype.pump_as_turbine:
        display_state(m)
    else:
        pass

    return m


def build(number_of_stages, erd_type=ERDtype.pump_as_turbine):

    # flowsheet set up
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.erd_type = erd_type
    m.fs.properties = props.NaClParameterBlock()
    m.fs.costing = WaterTAPCosting()

    # stage set up
    m.fs.NumberOfStages = Param(initialize=number_of_stages)
    m.fs.Stages = RangeSet(m.fs.NumberOfStages)
    m.fs.NonFirstStages = RangeSet(2, m.fs.NumberOfStages)
    m.fs.NonFinalStages = RangeSet(m.fs.NumberOfStages - 1)
    if number_of_stages > 1:
        m.fs.IntermediateStages = RangeSet(2, m.fs.NumberOfStages - 1)
    else:
        m.fs.IntermediateStages = RangeSet(0)
    m.fs.FirstStage = m.fs.Stages.first()
    m.fs.LastStage = m.fs.Stages.last()

    # Control volume flow blocks
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)
    m.fs.disposal = Product(property_package=m.fs.properties)

    # --- Main pump ---
    m.fs.PrimaryPumps = Pump(m.fs.Stages, property_package=m.fs.properties)
    for pump in m.fs.PrimaryPumps.values():
        pump.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
        )

    # --- Recycle pump ---
    m.fs.RecyclePumps = Pump(m.fs.NonFirstStages, property_package=m.fs.properties)
    for pump in m.fs.RecyclePumps.values():
        pump.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
        )

    m.fs.total_pump_work = Expression(
        expr=sum(
            pyunits.convert(pump.work_mechanical[0], to_units=pyunits.kW)
            for pump in itertools.chain(
                m.fs.PrimaryPumps.values(), m.fs.RecyclePumps.values()
            )
        )
    )

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
    m.fs.OAROUnits = OsmoticallyAssistedReverseOsmosis0D(
        m.fs.NonFinalStages,
        property_package=m.fs.properties,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        has_full_reporting=True,
    )
    for stage in m.fs.OAROUnits.values():
        stage.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method_arguments={"oaro_type": "standard"},
        )

    # --- ERD blocks ---
    if erd_type == ERDtype.pump_as_turbine:
        # add energy recovery turbine block
        m.fs.EnergyRecoveryDevices = EnergyRecoveryDevice(
            m.fs.Stages, property_package=m.fs.properties
        )
        # add costing for ERD config
        for erd in m.fs.EnergyRecoveryDevices.values():
            erd.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    else:
        erd_type_not_found(erd_type)

    m.fs.recovered_pump_work = Expression(
        expr=sum(
            pyunits.convert(erd.work_mechanical[0], to_units=pyunits.kW)
            for erd in m.fs.EnergyRecoveryDevices.values()
        )
    )
    m.fs.net_pump_work = Expression(
        expr=m.fs.total_pump_work + m.fs.recovered_pump_work
    )

    # additional parameters, variables or expressions ---------------------------------------------------------------------------
    m.fs.oaro_min_pressure = Param(initialize=10e5, units=pyunits.Pa, mutable=True)
    m.fs.oaro_max_pressure = Param(initialize=85e5, units=pyunits.Pa, mutable=True)
    m.fs.ro_min_pressure = Param(initialize=10e5, units=pyunits.Pa, mutable=True)
    m.fs.ro_max_pressure = Param(initialize=65e5, units=pyunits.Pa, mutable=True)
    m.fs.recycle_pump_min_pressure = Param(
        initialize=1e5, units=pyunits.Pa, mutable=True
    )
    m.fs.recycle_pump_max_pressure = Param(
        initialize=10e5, units=pyunits.Pa, mutable=True
    )

    # process costing and add system level metrics
    m.fs.costing.utilization_factor.fix(0.9)
    m.fs.costing.factor_total_investment.fix(2)
    m.fs.costing.factor_maintenance_labor_chemical.fix(0.03)
    m.fs.costing.factor_capital_annualization.fix(0.1)
    m.fs.costing.electricity_cost.set_value(0.07)
    m.fs.costing.reverse_osmosis.factor_membrane_replacement.fix(0.15)
    m.fs.costing.reverse_osmosis.membrane_cost.fix(30)
    m.fs.costing.reverse_osmosis.high_pressure_membrane_cost.fix(50)
    m.fs.costing.high_pressure_pump.cost.fix(53 / 1e5 * 3600)
    m.fs.costing.energy_recovery_device.pressure_exchanger_cost.fix(535)

    m.fs.costing.cost_process()
    product_flow_vol_total = m.fs.product.properties[0].flow_vol
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
        # Connect the feed to the first pump
        m.fs.feed_to_pump = Arc(
            source=m.fs.feed.outlet, destination=m.fs.PrimaryPumps[1].inlet
        )

        # Connect first EnergyRecoveryDevice to disposal
        m.fs.ERD_to_disposal = Arc(
            source=m.fs.EnergyRecoveryDevices[1].outlet, destination=m.fs.disposal.inlet
        )

        # Connect the Pump n to OARO n
        m.fs.pump_to_OARO = Arc(
            m.fs.NonFinalStages,
            rule=lambda fs, n: {
                "source": fs.PrimaryPumps[n].outlet,
                "destination": fs.OAROUnits[n].feed_inlet,
            },
        )

        # Connect OARO n permeate outlet to Pump n+1
        m.fs.OARO_to_pump = Arc(
            m.fs.NonFinalStages,
            rule=lambda fs, n: {
                "source": fs.OAROUnits[n].permeate_outlet,
                "destination": fs.PrimaryPumps[n + 1].inlet,
            },
        )

        # Connect OARO n feed outlet to EnergyRecoveryDevice n
        m.fs.OARO_to_ERD = Arc(
            m.fs.NonFinalStages,
            rule=lambda fs, n: {
                "source": fs.OAROUnits[n].feed_outlet,
                "destination": fs.EnergyRecoveryDevices[n].inlet,
            },
        )

        # Connect EnergyRecoveryDevice n to RecyclePumps n
        m.fs.ERD_to_recyclepump = Arc(
            m.fs.NonFirstStages,
            rule=lambda fs, n: {
                "source": fs.EnergyRecoveryDevices[n].outlet,
                "destination": fs.RecyclePumps[n].inlet,
            },
        )

        # Connect RecyclePumps n to OARO n-1 permeate_inlet
        m.fs.recyclepump_to_OARO = Arc(
            m.fs.NonFirstStages,
            rule=lambda fs, n: {
                "source": fs.RecyclePumps[n].outlet,
                "destination": fs.OAROUnits[n - 1].permeate_inlet,
            },
        )

        last_stage = m.fs.LastStage
        # Connect the last Pump to RO
        m.fs.pump_to_ro = Arc(
            source=m.fs.PrimaryPumps[last_stage].outlet, destination=m.fs.RO.inlet
        )

        # Connect the primary RO permeate to the product
        m.fs.ro_to_product = Arc(
            source=m.fs.RO.permeate, destination=m.fs.product.inlet
        )

        # Connect RO retentate to the last EnergyRecoveryDevice
        m.fs.ro_to_ERD = Arc(
            source=m.fs.RO.retentate,
            destination=m.fs.EnergyRecoveryDevices[last_stage].inlet,
        )

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
    for pump in m.fs.PrimaryPumps.values():
        iscale.set_scaling_factor(pump.control_volume.work, 1e-3)

    for pump in m.fs.RecyclePumps.values():
        iscale.set_scaling_factor(pump.control_volume.work, 1e-3)

    iscale.set_scaling_factor(m.fs.RO.area, 1e-2)

    for stage in m.fs.NonFinalStages:
        iscale.set_scaling_factor(m.fs.OAROUnits[stage].area, 1e-2)

    m.fs.feed.properties[0].flow_vol_phase["Liq"]
    m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"]

    if erd_type == ERDtype.pump_as_turbine:
        for erd in m.fs.EnergyRecoveryDevices.values():
            iscale.set_scaling_factor(erd.control_volume.work, 1e-3)
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
    print(f"DOF before set: {degrees_of_freedom(m)}")
    pressure_atmospheric = 101325
    feed_temperature = 273.15 + 25
    m.fs.feed.properties[0].pressure.fix(pressure_atmospheric)  # feed pressure [Pa]
    m.fs.feed.properties[0].temperature.fix(feed_temperature)  # feed temperature [K]

    # properties (cannot be fixed for initialization routines, must calculate the state variables)
    feed_flow_mass = 1
    feed_mass_frac_NaCl = 0.07
    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )

    # primary pumps
    for idx, pump in m.fs.PrimaryPumps.items():
        pump.control_volume.properties_out[0].pressure = 65e5
        pump.efficiency_pump.fix(0.80)
        pump.control_volume.properties_out[0].pressure.fix()

    for pump in m.fs.RecyclePumps.values():
        pump.control_volume.properties_out[0].pressure = 4e5
        pump.efficiency_pump.fix(0.80)
        pump.control_volume.properties_out[0].pressure.fix()

    # initial guess for states of recycle pumps (temperature and concentrations)
    # for stage in m.fs.NonFirstStages:
    #     m.fs.RecyclePumps[stage].control_volume.properties_out[
    #         0
    #     ].temperature.value = feed_temperature
    #     permeate_flow_mass = 1 / stage + 0.3
    #     permeate_mass_frac_NaCl = 0.07 / stage
    #     permeate_mass_frac_H2O = 1 - permeate_mass_frac_NaCl
    #     m.fs.RecyclePumps[stage].control_volume.properties_out[0].flow_mass_phase_comp[
    #         "Liq", "H2O"
    #     ].value = (permeate_flow_mass * permeate_mass_frac_H2O)
    #     m.fs.RecyclePumps[stage].control_volume.properties_out[0].flow_mass_phase_comp[
    #         "Liq", "NaCl"
    #     ].value = (permeate_flow_mass * permeate_mass_frac_NaCl)

    # Initialize OARO
    membrane_area = 50
    A = 4.2e-12
    B = 1.3e-8

    for stage in m.fs.OAROUnits.values():

        stage.area.fix(membrane_area)

        stage.A_comp.fix(A)
        stage.B_comp.fix(B)

        stage.structural_parameter.fix(300e-6)

        stage.permeate_side.channel_height.fix(0.001)
        stage.permeate_side.spacer_porosity.fix(0.75)
        stage.feed_side.channel_height.fix(0.002)
        stage.feed_side.spacer_porosity.fix(0.75)
        stage.feed_side.velocity[0, 0].fix(0.1)

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
        for erd in m.fs.EnergyRecoveryDevices.values():
            erd.efficiency_pump.fix(0.95)
            erd.control_volume.properties_out[0].pressure.fix(pressure_atmospheric)
    else:
        erd_type_not_found(m.fs.erd_type)

    print(f"DOF after set: {degrees_of_freedom(m)}")

    # check degrees of freedom
    if degrees_of_freedom(m) != 0:
        print(
            "The set_operating_conditions function resulted in {} "
            "degrees of freedom rather than 0. This error suggests "
            "that too many or not enough variables are fixed for a "
            "simulation.".format(degrees_of_freedom(m))
        )
        raise


def recycle_pump_initializer(pump, oaro, solvent_multiplier, solute_multiplier):
    # for vname in pump.control_volume.properties_in[0].vars:
    #     if vname == "flow_mass_phase_comp":
    #         for phase, comp in pump.control_volume.properties_in[0].vars[vname]:
    #             if comp in pump.config.property_package.solute_set:
    #                 pump.control_volume.properties_out[0].vars[vname][
    #                     phase, comp
    #                 ].value = (
    #                     solute_multiplier
    #                     * pump.control_volume.properties_out[0]
    #                     .vars[vname][phase, comp]
    #                     .value
    #                 )
    #             elif comp in pump.config.property_package.solvent_set:
    #                 pump.control_volume.properties_out[0].vars[vname][
    #                     phase, comp
    #                 ].value = (
    #                     solvent_multiplier
    #                     * pump.control_volume.properties_in[0]
    #                     .vars[vname][phase, comp]
    #                     .value
    #                 )
    #             else:
    #                 raise RuntimeError(f"Unknown component {comp}")
    #     else:  # copy the state
    #         for idx in pump.control_volume.properties_in[0].vars[vname]:
    #             pump.control_volume.properties_out[0].vars[vname][idx].value = (
    #                 pump.control_volume.properties_in[0].vars[vname][idx].value
    #             )

    feed_temperature = 273.15 + 25
    pump.control_volume.properties_out[0].temperature.value = feed_temperature
    pump.control_volume.properties_out[0].flow_mass_phase_comp["Liq", "H2O"].value = (
        oaro.feed_inlet.flow_mass_phase_comp[0, "Liq", "H2O"] * solvent_multiplier
    )
    pump.control_volume.properties_out[0].flow_mass_phase_comp["Liq", "NaCl"].value = (
        oaro.feed_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"] * solute_multiplier
    )


def solve(blk, solver=None, tee=True):
    if solver is None:
        solver = get_solver()
    results = solver.solve(blk, tee=tee)
    if not check_optimal_termination(results):
        results = solver.solve(blk, tee=tee)
    return results


def initialize_loop(m, solver):

    for stage in m.fs.IntermediateStages:
        propagate_state(m.fs.OARO_to_pump[stage - 1])
        m.fs.PrimaryPumps[stage].initialize()

        # ---initialize loop---
        propagate_state(m.fs.pump_to_OARO[stage])
        # recycle_pump_initializer(
        #     m.fs.RecyclePumps[stage + 1],
        #     m.fs.OAROUnits[stage],
        #     solvent_multiplier=0.8,
        #     solute_multiplier=0.5,
        # )
        feed_temperature = 273.15 + 25
        m.fs.RecyclePumps[stage + 1].control_volume.properties_out[
            0
        ].temperature.value = feed_temperature
        m.fs.RecyclePumps[stage + 1].control_volume.properties_out[
            0
        ].flow_mass_phase_comp["Liq", "H2O"].value = (
            m.fs.OAROUnits[stage].feed_inlet.flow_mass_phase_comp[0, "Liq", "H2O"] * 0.8
        )
        m.fs.RecyclePumps[stage + 1].control_volume.properties_out[
            0
        ].flow_mass_phase_comp["Liq", "NaCl"].value = (
            m.fs.OAROUnits[stage].feed_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"]
            * 0.5
        )
        propagate_state(m.fs.recyclepump_to_OARO[stage + 1])
        m.fs.OAROUnits[stage].initialize()

        propagate_state(m.fs.OARO_to_ERD[stage])
        m.fs.EnergyRecoveryDevices[stage].initialize()

        propagate_state(m.fs.ERD_to_recyclepump[stage])
        m.fs.RecyclePumps[stage].initialize()

        propagate_state(m.fs.recyclepump_to_OARO[stage])
        propagate_state(m.fs.pump_to_OARO[stage - 1])
        m.fs.OAROUnits[stage - 1].initialize()


def initialize_system(m, solver=None, verbose=True):
    if solver is None:
        solver = get_solver()

    last_stage = m.fs.LastStage
    first_stage = m.fs.FirstStage

    # ---initialize feed block---
    m.fs.feed.initialize()

    propagate_state(m.fs.feed_to_pump)
    m.fs.PrimaryPumps[first_stage].initialize()

    # ---initialize first OARO unit---
    propagate_state(m.fs.pump_to_OARO[first_stage])
    recycle_pump_initializer(
        m.fs.RecyclePumps[first_stage + 1],
        m.fs.OAROUnits[first_stage],
        solvent_multiplier=0.8,
        solute_multiplier=0.5,
    )
    # print(m.fs.recyclepump_to_OARO[first_stage + 1].destination.name)
    # m.fs.recyclepump_to_OARO[first_stage + 1].destination.display()
    propagate_state(m.fs.recyclepump_to_OARO[first_stage + 1])
    m.fs.OAROUnits[first_stage].initialize()

    print(f"DOF before loop: {degrees_of_freedom(m)}")
    initialize_loop(m, solver)
    print(f"DOF after loop: {degrees_of_freedom(m)}")

    # ---initialize RO loop---
    propagate_state(m.fs.OARO_to_pump[last_stage - 1])
    m.fs.PrimaryPumps[last_stage].initialize()

    propagate_state(m.fs.pump_to_ro)
    print(f"DOF after prop_state to RO: {degrees_of_freedom(m)}")
    print(f"fixed variables set after prop_state to RO: {fixed_variables_set(m.fs.RO)}")
    m.fs.RO.initialize()
    print(f"fixed variables set after RO: {fixed_variables_set(m.fs.RO)}")
    print(f"DOF after RO: {degrees_of_freedom(m)}")

    propagate_state(m.fs.ro_to_ERD)
    m.fs.EnergyRecoveryDevices[last_stage].initialize()

    propagate_state(m.fs.ERD_to_recyclepump[last_stage])
    m.fs.RecyclePumps[last_stage].initialize()

    propagate_state(m.fs.recyclepump_to_OARO[last_stage])
    propagate_state(m.fs.pump_to_OARO[last_stage - 1])
    m.fs.OAROUnits[last_stage - 1].initialize()

    # ---initialize first ERD---
    propagate_state(m.fs.OARO_to_ERD[first_stage])
    m.fs.EnergyRecoveryDevices[first_stage].initialize()

    print(f"DOF: {degrees_of_freedom(m)}")

    # Now that the units are initialized, we can fix the
    # permeate side outlet pressure and unfix the RO pump
    # (which allows for control over the flow mass composition
    # into the OARO permeate_side).
    for stage in m.fs.NonFinalStages:
        m.fs.OAROUnits[stage].permeate_side.properties_out[0].pressure.fix(101325)
        m.fs.PrimaryPumps[stage + 1].control_volume.properties_out[0].pressure.unfix()

    print(f"DOF: {degrees_of_freedom(m)}")

    m.fs.costing.initialize()


def optimize_set_up(m):
    # add objective
    m.fs.objective = Objective(expr=m.fs.costing.LCOW)

    # Primary pumps
    for idx, pump in m.fs.PrimaryPumps.items():
        pump.control_volume.properties_out[0].pressure.unfix()
        pump.deltaP.setlb(0)
        if idx < m.fs.Stages.last():
            pump.max_oaro_pressure_con = Constraint(
                expr=(
                    m.fs.oaro_min_pressure,
                    pump.control_volume.properties_out[0].pressure,
                    m.fs.oaro_max_pressure,
                )
            )
            iscale.constraint_scaling_transform(
                pump.max_oaro_pressure_con,
                iscale.get_scaling_factor(
                    pump.control_volume.properties_out[0].pressure
                ),
            )
        else:
            pump.max_ro_pressure_con = Constraint(
                expr=(
                    m.fs.ro_min_pressure,
                    pump.control_volume.properties_out[0].pressure,
                    m.fs.ro_max_pressure,
                )
            )
            iscale.constraint_scaling_transform(
                pump.max_ro_pressure_con,
                iscale.get_scaling_factor(
                    pump.control_volume.properties_out[0].pressure
                ),
            )

    # Recycle pumps
    for idx, pump in m.fs.RecyclePumps.items():
        pump.control_volume.properties_out[0].pressure.unfix()
        pump.max_recycle_pump_pressure_con = Constraint(
            expr=(
                m.fs.recycle_pump_min_pressure,
                pump.control_volume.properties_out[0].pressure,
                m.fs.recycle_pump_max_pressure,
            )
        )
        iscale.constraint_scaling_transform(
            pump.max_recycle_pump_pressure_con,
            iscale.get_scaling_factor(pump.control_volume.properties_out[0].pressure),
        )
        pump.deltaP.setlb(0)

    # # OARO Units
    # for stage in m.fs.OAROUnits.values():
    #     stage.area.unfix()
    #     stage.area.setlb(1)
    #     stage.area.setub(2000)
    #
    # # RO
    # m.fs.RO.area.unfix()
    # m.fs.RO.area.setlb(1)
    # m.fs.RO.area.setub(2000)

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
    # assert_degrees_of_freedom(m, 2 * m.fs.NumberOfStages)


def display_system(m):
    print("----system metrics----")
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

    brine_flow_mass = sum(
        m.fs.disposal.flow_mass_phase_comp[0, "Liq", j].value for j in ["H2O", "NaCl"]
    )
    brine_mass_frac_NaCl = (
        m.fs.disposal.flow_mass_phase_comp[0, "Liq", "NaCl"].value / brine_flow_mass
    )
    print("Brine: %.3f kg/s, %.0f ppm" % (brine_flow_mass, brine_mass_frac_NaCl * 1e6))

    print("Volumetric water recovery: %.1f%%" % (value(m.fs.water_recovery) * 100))
    print(f"Number of Stages: {value(m.fs.NumberOfStages)}")
    total_area = value(
        sum(m.fs.OAROUnits[a].area for a in m.fs.NonFinalStages) + m.fs.RO.area
    )
    print(f"Total Membrane Area: {total_area:.2f}")
    print(
        "Energy Consumption: %.1f kWh/m3"
        % value(m.fs.costing.specific_energy_consumption)
    )

    print("Levelized cost of water: %.2f $/m3" % value(m.fs.costing.LCOW))
    print(
        f"Primary Pump Capital Cost ($/m3):"
        f"{value(m.fs.costing.factor_capital_annualization*sum(m.fs.PrimaryPumps[stage].costing.capital_cost for stage in m.fs.Stages)/ m.fs.costing.annual_water_production)}"
    )
    print(
        f"Recycle Pump Capital Cost ($/m3): "
        f"{value(m.fs.costing.factor_capital_annualization*sum(m.fs.RecyclePumps[stage].costing.capital_cost for stage in m.fs.NonFirstStages) / m.fs.costing.annual_water_production)}"
    )
    print(
        f"ERD Capital Cost ($/m3):"
        f"{value(m.fs.costing.factor_capital_annualization*sum(erd.costing.capital_cost for erd in m.fs.EnergyRecoveryDevices.values()) / m.fs.costing.annual_water_production)}"
    )
    print(
        f"Membrane Capital Cost ($/m3): "
        f"{value(m.fs.costing.factor_capital_annualization*(sum(m.fs.OAROUnits[stage].costing.capital_cost for stage in m.fs.NonFinalStages) + m.fs.RO.costing.capital_cost) / m.fs.costing.annual_water_production)}"
    )
    print(
        f"Indirect Capital Cost ($/m3): "
        f"{value(m.fs.costing.factor_capital_annualization*(m.fs.costing.total_capital_cost - m.fs.costing.aggregate_capital_cost) / m.fs.costing.annual_water_production)}"
    )
    electricity_cost = value(
        m.fs.costing.aggregate_flow_costs["electricity"]
        * m.fs.costing.utilization_factor
        / m.fs.costing.annual_water_production
    )
    print(f"Electricity cost ($/m3): {electricity_cost}")

    print("\n")


def display_design(m):
    print("--decision variables--")
    for stage in m.fs.NonFinalStages:
        print(
            "OARO Stage %d feed operating pressure %.1f bar"
            % (stage, m.fs.OAROUnits[stage].feed_inlet.pressure[0].value / 1e5)
        )
        print(
            "OARO Stage %d permeate operating pressure %.1f bar"
            % (stage, m.fs.OAROUnits[stage].permeate_inlet.pressure[0].value / 1e5)
        )
        print(
            "OARO tage %d membrane area      %.1f m2"
            % (stage, m.fs.OAROUnits[stage].area.value)
        )
        print(
            "OARO Stage %d water perm. coeff.  %.1f LMH/bar"
            % (stage, m.fs.OAROUnits[stage].A_comp[0, "H2O"].value * (3.6e11))
        )
        print(
            "OARO Stage %d salt perm. coeff.  %.1f LMH"
            % (stage, m.fs.OAROUnits[stage].B_comp[0, "NaCl"].value * (1000.0 * 3600.0))
        )
    print(
        "RO feed operating pressure %.1f bar" % (m.fs.RO.inlet.pressure[0].value / 1e5)
    )
    print(
        "RO permeate operating pressure %.1f bar"
        % (m.fs.RO.permeate.pressure[0].value / 1e5)
    )
    print("RO membrane area      %.1f m2" % (m.fs.RO.area.value))
    print(
        "RO water perm. coeff.  %.1f LMH/bar"
        % (m.fs.RO.A_comp[0, "H2O"].value * (3.6e11))
    )
    print(
        "RO salt perm. coeff.  %.1f LMH"
        % (m.fs.RO.B_comp[0, "NaCl"].value * (1000.0 * 3600.0))
    )


def display_state(m):
    print("--------state---------")

    def print_state(s, b):
        flow_mass = sum(
            b.flow_mass_phase_comp[0, "Liq", j].value for j in ["H2O", "NaCl"]
        )
        mass_frac_ppm = b.flow_mass_phase_comp[0, "Liq", "NaCl"].value / flow_mass * 1e6
        pressure_bar = b.pressure[0].value / 1e5
        print(
            s.ljust(20)
            + ": %.3f kg/s, %.0f ppm, %.1f bar"
            % (flow_mass, mass_frac_ppm, pressure_bar)
        )

    print_state("Feed", m.fs.feed.outlet)

    for stage in m.fs.Stages:

        print_state(f"Primary Pump {stage} out", m.fs.PrimaryPumps[stage].outlet)
        print_state(
            f"Energy Recovery Device {stage} out",
            m.fs.EnergyRecoveryDevices[stage].outlet,
        )

        if stage == m.fs.LastStage:
            pass
        else:
            print_state(f"OARO {stage} feed outlet", m.fs.OAROUnits[stage].feed_outlet)
            print_state(
                f"OARO {stage} permeate outlet", m.fs.OAROUnits[stage].permeate_outlet
            )

        if stage == m.fs.FirstStage:
            pass
        else:
            print_state(f"Recycle Pump {stage} out", m.fs.RecyclePumps[stage].outlet)

    print_state(f"RO permeate", m.fs.RO.permeate)
    print_state(f"RO retentate", m.fs.RO.retentate)

    print_state(f"Disposal", m.fs.disposal.inlet)
    print_state(f"Product", m.fs.product.inlet)


if __name__ == "__main__":
    m = main(5, erd_type=ERDtype.pump_as_turbine)
