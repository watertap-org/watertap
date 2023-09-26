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
    NonNegativeReals,
    RangeSet,
    check_optimal_termination,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import (
    propagate_state as _pro_state,
)
from idaes.models.unit_models import Mixer, Separator, Product, Feed
from idaes.core import UnitModelCostingBlock
import idaes.core.util.scaling as iscale
from idaes.core.util.misc import StrEnum

import watertap.property_models.NaCl_prop_pack as props
from watertap.core.util.initialization import (
    assert_degrees_of_freedom,
)
from watertap.unit_models.reverse_osmosis_0D import (
    ReverseOsmosis0D,
    ConcentrationPolarizationType,
)
from watertap.unit_models.osmotically_assisted_reverse_osmosis_0D import (
    OsmoticallyAssistedReverseOsmosis0D,
    MassTransferCoefficient,
    PressureChangeType,
)
from watertap.unit_models.pressure_changer import Pump, EnergyRecoveryDevice
from watertap.costing import WaterTAPCosting
from watertap.costing.unit_models.pump import cost_low_pressure_pump


class ERDtype(StrEnum):
    pump_as_turbine = "pump_as_turbine"


def erd_type_not_found(erd_type):
    raise NotImplementedError(
        "erd_type was {}, but can only " "be pump_as_turbine" "".format(erd_type.value)
    )


def propagate_state(arc):
    _pro_state(arc)
    print(arc.destination.name)
    arc.destination.display()


def main(number_of_stages, system_recovery, erd_type=ERDtype.pump_as_turbine):
    # set up solver
    solver = get_solver()

    # build, set, and initialize
    m = build(number_of_stages=number_of_stages, erd_type=erd_type)
    set_operating_conditions(m)
    initialize_system(
        m,
        number_of_stages,
        solvent_multiplier=0.5,
        solute_multiplier=0.7,
        solver=solver,
    )

    optimize_set_up(
        m, number_of_stages=number_of_stages, water_recovery=system_recovery
    )

    results = solve(m, solver=solver)
    assert_optimal_termination(results)

    print("\n***---Optimization results---***")
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
    if number_of_stages > 2:
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
            flowsheet_costing_block=m.fs.costing, costing_method=cost_low_pressure_pump
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
            costing_method_arguments={"oaro_type": "high_pressure"},
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

    # --- Separators ---
    m.fs.Separators = Separator(
        m.fs.NonFirstStages,
        property_package=m.fs.properties,
        outlet_list=["treat", "purge"],
    )
    if number_of_stages > 1:
        separator_list = ["product"]
        for i in range(2, number_of_stages + 1):
            separator_list.append("makeup" + str(i))
        m.fs.ProductSeparator = Separator(
            property_package=m.fs.properties, outlet_list=separator_list
        )

    # --- Mixers ---
    m.fs.Mixers = Mixer(
        m.fs.NonFirstStages,
        property_package=m.fs.properties,
        inlet_list=["treat", "makeup"],
    )
    if number_of_stages > 1:
        mixer_list = []
        for i in range(1, number_of_stages + 1):
            mixer_list.append("purge" + str(i))
        m.fs.WasteMixer = Mixer(property_package=m.fs.properties, inlet_list=mixer_list)

    m.fs.recovered_pump_work = Expression(
        expr=sum(
            pyunits.convert(erd.work_mechanical[0], to_units=pyunits.kW)
            for erd in m.fs.EnergyRecoveryDevices.values()
        )
    )
    m.fs.net_pump_work = Expression(
        expr=m.fs.total_pump_work + m.fs.recovered_pump_work
    )

    # additional parameters, variables or expressions
    m.fs.oaro_min_pressure = Param(initialize=10e5, units=pyunits.Pa, mutable=True)
    m.fs.oaro_max_pressure = Param(initialize=65e5, units=pyunits.Pa, mutable=True)
    m.fs.ro_min_pressure = Param(initialize=10e5, units=pyunits.Pa, mutable=True)
    m.fs.ro_max_pressure = Param(initialize=85e5, units=pyunits.Pa, mutable=True)
    m.fs.recycle_pump_min_pressure = Param(
        initialize=1e5, units=pyunits.Pa, mutable=True
    )
    m.fs.recycle_pump_max_pressure = Param(
        initialize=85e5, units=pyunits.Pa, mutable=True
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
    m.fs.costing.add_annual_water_production(product_flow_vol_total)
    m.fs.costing.add_LCOW(product_flow_vol_total)
    m.fs.costing.add_specific_energy_consumption(product_flow_vol_total)
    m.fs.costing.add_specific_electrical_carbon_intensity(product_flow_vol_total)

    # Expressions for parameter sweep -----------------------------------------
    # Final permeate concentration as mass fraction
    m.fs.product.properties[0].mass_frac_phase_comp

    # Touch feed concentration as mass concentration
    m.fs.feed.properties[0].conc_mass_phase_comp

    # Touch final brine concentration as mass concentration
    m.fs.disposal.properties[0].conc_mass_phase_comp

    # Touch final brine concentration as mass fraction
    m.fs.disposal.properties[0].mass_frac_phase_comp

    m.fs.system_salt_rejection = Expression(
        expr=1
        - m.fs.product.properties[0].conc_mass_phase_comp["Liq", "NaCl"]
        / m.fs.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"]
    )
    m.fs.annual_feed = Expression(
        expr=pyunits.convert(
            m.fs.feed.properties[0].flow_vol, to_units=pyunits.m**3 / pyunits.year
        )
        * m.fs.costing.utilization_factor
    )
    m.fs.final_permeate_concentration = Expression(
        expr=m.fs.product.flow_mass_phase_comp[0, "Liq", "NaCl"]
        / sum(m.fs.product.flow_mass_phase_comp[0, "Liq", j] for j in ["H2O", "NaCl"])
    )
    m.fs.costing.add_LCOW(m.fs.feed.properties[0].flow_vol, name="LCOW_feed")

    m.fs.costing.add_specific_energy_consumption(
        m.fs.feed.properties[0].flow_vol, name="specific_energy_consumption_feed"
    )

    m.fs.costing.primary_pump_capex_lcow = Expression(
        expr=m.fs.costing.factor_capital_annualization
        * sum(m.fs.PrimaryPumps[n].costing.capital_cost for n in m.fs.Stages)
        / m.fs.costing.annual_water_production
    )

    m.fs.total_membrane_area = Expression(
        expr=sum(oaro.area for oaro in m.fs.OAROUnits.values()) + m.fs.RO.area
    )

    m.fs.costing.recycle_pump_capex_lcow = Expression(
        expr=m.fs.costing.factor_capital_annualization
        * (
            sum(m.fs.RecyclePumps[n].costing.capital_cost for n in m.fs.NonFirstStages)
            if number_of_stages > 1
            else 0 * m.fs.costing.base_currency
        )
        / m.fs.costing.annual_water_production,
    )

    m.fs.costing.erd_capex_lcow = Expression(
        expr=m.fs.costing.factor_capital_annualization
        * sum(erd.costing.capital_cost for erd in m.fs.EnergyRecoveryDevices.values())
        / m.fs.costing.annual_water_production
    )

    m.fs.costing.electricity_lcow = Expression(
        expr=m.fs.costing.aggregate_flow_costs["electricity"]
        * m.fs.costing.utilization_factor
        / m.fs.costing.annual_water_production
    )

    m.fs.costing.pumping_energy_aggregate_lcow = Expression(
        expr=m.fs.costing.factor_total_investment
        * (
            m.fs.costing.primary_pump_capex_lcow
            + (
                m.fs.costing.recycle_pump_capex_lcow
                if number_of_stages > 1
                else 0 * m.fs.costing.base_currency
            )
            + m.fs.costing.erd_capex_lcow
        )
        * (
            1
            + m.fs.costing.factor_maintenance_labor_chemical
            / m.fs.costing.factor_capital_annualization
        )
        + m.fs.costing.electricity_lcow
    )

    m.fs.costing.membrane_capex_lcow = Expression(
        expr=m.fs.costing.factor_capital_annualization
        * (
            sum(m.fs.OAROUnits[n].costing.capital_cost for n in m.fs.NonFinalStages)
            + m.fs.RO.costing.capital_cost
        )
        / m.fs.costing.annual_water_production
    )

    m.fs.costing.indirect_capex_lcow = Expression(
        expr=m.fs.costing.factor_capital_annualization
        * (m.fs.costing.total_capital_cost - m.fs.costing.aggregate_capital_cost)
        / m.fs.costing.annual_water_production
    )

    m.fs.costing.membrane_replacement_lcow = Expression(
        expr=(
            sum(
                m.fs.OAROUnits[n].costing.fixed_operating_cost
                for n in m.fs.NonFinalStages
            )
            + m.fs.RO.costing.fixed_operating_cost
        )
        / m.fs.costing.annual_water_production
    )

    m.fs.costing.chemical_labor_maintenance_lcow = Expression(
        expr=m.fs.costing.maintenance_labor_chemical_operating_cost
        / m.fs.costing.annual_water_production
    )

    m.fs.costing.membrane_aggregate_lcow = Expression(
        expr=m.fs.costing.factor_total_investment
        * m.fs.costing.membrane_capex_lcow
        * (
            1
            + m.fs.costing.factor_maintenance_labor_chemical
            / m.fs.costing.factor_capital_annualization
        )
        + m.fs.costing.membrane_replacement_lcow
    )

    # system water recovery
    m.fs.water_recovery = Var(
        initialize=0.5,
        bounds=(0, 1),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="System Volumetric Recovery of Water",
    )
    m.fs.eq_water_recovery = Constraint(
        expr=(
            m.fs.feed.properties[0].flow_vol_phase["Liq"]
            + sum(
                m.fs.Mixers[stage].makeup_state[0].flow_vol_phase["Liq"]
                for stage in m.fs.NonFirstStages
            )
        )
        * m.fs.water_recovery
        == m.fs.product.properties[0].flow_vol_phase["Liq"]
    )

    m.fs.mass_water_recovery = Var(
        initialize=0.5,
        bounds=(0, 1),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="System Mass Recovery of Water",
    )
    m.fs.eq_mass_water_recovery = Constraint(
        expr=(
            m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]
            + sum(
                m.fs.Mixers[stage].makeup_state[0].flow_mass_phase_comp["Liq", "H2O"]
                for stage in m.fs.NonFirstStages
            )
        )
        * m.fs.mass_water_recovery
        == m.fs.product.properties[0].flow_mass_phase_comp["Liq", "H2O"]
    )

    # connections
    if erd_type == ERDtype.pump_as_turbine:
        last_stage = m.fs.LastStage
        first_stage = m.fs.FirstStage

        # Connect the feed to the first pump
        m.fs.feed_to_pump = Arc(
            source=m.fs.feed.outlet, destination=m.fs.PrimaryPumps[first_stage].inlet
        )

        if number_of_stages > 1:
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

            # Connect EnergyRecoveryDevice n to Separator n
            m.fs.ERD_to_separator = Arc(
                m.fs.NonFirstStages,
                rule=lambda fs, n: {
                    "source": fs.EnergyRecoveryDevices[n].outlet,
                    "destination": fs.Separators[n].inlet,
                },
            )

            m.fs.separator_to_mixer = Arc(
                m.fs.NonFirstStages,
                rule=lambda fs, n: {
                    "source": fs.Separators[n].treat,
                    "destination": fs.Mixers[n].treat,
                },
            )

            m.fs.makeup_to_mixer = Arc(
                m.fs.NonFirstStages,
                rule=lambda fs, n: {
                    "source": getattr(fs.ProductSeparator, "makeup" + str(n)),
                    "destination": fs.Mixers[n].makeup,
                },
            )

            m.fs.mixer_to_recyclepump = Arc(
                m.fs.NonFirstStages,
                rule=lambda fs, n: {
                    "source": fs.Mixers[n].outlet,
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

            # Connect purge streams to waste mixer
            m.fs.purge_to_wastemixer = Arc(
                m.fs.NonFirstStages,
                rule=lambda fs, n: {
                    "source": fs.Separators[n].purge,
                    "destination": getattr(fs.WasteMixer, "purge" + str(n)),
                },
            )

        if number_of_stages > 1:
            # Connect first EnergyRecoveryDevice to Mixer
            m.fs.ERD_to_wastemixer = Arc(
                source=m.fs.EnergyRecoveryDevices[first_stage].outlet,
                destination=m.fs.WasteMixer.purge1,
            )

            # Connect waste mixer to disposal
            m.fs.wastemixer_to_disposal = Arc(
                source=m.fs.WasteMixer.outlet,
                destination=m.fs.disposal.inlet,
            )
        else:
            m.fs.ERD_to_disposal = Arc(
                source=m.fs.EnergyRecoveryDevices[first_stage].outlet,
                destination=m.fs.disposal.inlet,
            )

        # Connect the last Pump to RO
        m.fs.pump_to_ro = Arc(
            source=m.fs.PrimaryPumps[last_stage].outlet, destination=m.fs.RO.inlet
        )

        # Connect the primary RO permeate to the product
        # m.fs.ro_to_product = Arc(
        #     source=m.fs.RO.permeate, destination=m.fs.product.inlet
        # )
        # Connect the primary RO permeate to the product
        if number_of_stages > 1:
            m.fs.ro_to_productseparator = Arc(
                source=m.fs.RO.permeate, destination=m.fs.ProductSeparator.inlet
            )
            m.fs.productseparator_to_product = Arc(
                source=m.fs.ProductSeparator.product, destination=m.fs.product.inlet
            )
        else:
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
    pressure_atmospheric = 101325
    feed_temperature = 273.15 + 25
    m.fs.feed.properties[0].pressure.fix(pressure_atmospheric)  # feed pressure [Pa]
    m.fs.feed.properties[0].temperature.fix(feed_temperature)  # feed temperature [K]

    Cin = 75000e-6
    Qin = 5.416667e-3
    m.fs.feed.properties.calculate_state(
        var_args={
            ("mass_frac_phase_comp", ("Liq", "NaCl")): value(
                Cin
            ),  # feed mass concentration
            ("flow_vol_phase", "Liq"): Qin,
        },  # volumetric feed flowrate [-]
        hold_state=True,  # fixes the calculated component mass flow rates
    )

    # primary pumps
    for idx, pump in m.fs.PrimaryPumps.items():
        # pump.control_volume.properties_out[0].pressure = 10e5 + 65e5 / float(idx)
        if idx == m.fs.Stages.last():
            pump.control_volume.properties_out[0].pressure = 85e5
        else:
            pump.control_volume.properties_out[0].pressure = 65e5
        pump.control_volume.properties_out[0].pressure.fix()
        pump.efficiency_pump.fix(0.75)

    for idx, pump in m.fs.RecyclePumps.items():
        pump.control_volume.properties_out[0].pressure = 2e5
        pump.control_volume.properties_out[0].pressure.fix()
        pump.efficiency_pump.fix(0.75)

    # Initialize OARO
    A_OARO = 1.0e-12
    B_OARO = 8.0e-8
    spacer_porosity = 0.75

    for idx, stage in m.fs.OAROUnits.items():
        stage.A_comp.fix(A_OARO)
        stage.B_comp.fix(B_OARO)

        stage.structural_parameter.fix(1200e-6)

        stage.permeate_side.channel_height.fix(2e-3)
        stage.permeate_side.spacer_porosity.fix(spacer_porosity)
        stage.feed_side.channel_height.fix(2e-3)
        stage.feed_side.spacer_porosity.fix(spacer_porosity)
        stage.feed_side.velocity[0, 0].fix(0.1)
        stage.recovery_mass_phase_comp[0, "Liq", "H2O"].fix(0.3)

    # RO unit
    A_RO = 4.2e-12
    B_RO = 3.5e-8
    width = 5 * Qin / 1e-3  # effective membrane width [m]
    m.fs.RO.A_comp.fix(A_RO)  # membrane water permeability coefficient [m/s-Pa]
    m.fs.RO.B_comp.fix(B_RO)  # membrane salt permeability coefficient [m/s]
    m.fs.RO.feed_side.channel_height.fix(2e-3)  # channel height in membrane stage [m]
    m.fs.RO.feed_side.spacer_porosity.fix(
        spacer_porosity
    )  # spacer porosity in membrane stage [-]
    m.fs.RO.permeate.pressure[0].fix(101325)  # atmospheric pressure [Pa]
    m.fs.RO.width.fix(width / float(m.fs.NumberOfStages))  # stage width [m]
    m.fs.RO.recovery_mass_phase_comp[0, "Liq", "H2O"].fix(0.3)

    # ERD unit
    if m.fs.erd_type == ERDtype.pump_as_turbine:
        # energy recovery turbine - efficiency and outlet pressure
        for erd in m.fs.EnergyRecoveryDevices.values():
            erd.efficiency_pump.fix(0.9)
            erd.control_volume.properties_out[0].pressure.fix(pressure_atmospheric)
    else:
        erd_type_not_found(m.fs.erd_type)

    # Makeup
    for i in range(2, m.fs.NumberOfStages + 1):
        m.fs.ProductSeparator.split_fraction[:, "makeup" + str(i)].fix(0)

    # Separator
    for sep in m.fs.Separators.values():
        sep.split_fraction[:, "treat"] = 1
        sep.split_fraction[:, "treat"].fix()

    # check degrees of freedom
    if degrees_of_freedom(m) != 0:
        raise RuntimeError(
            "The set_operating_conditions function resulted in {} "
            "degrees of freedom rather than 0. This error suggests "
            "that too many or not enough variables are fixed for a "
            "simulation.".format(degrees_of_freedom(m))
        )


def purge_initializer(sep, oaro, solute_multiplier):
    split_fraction_treat = value(
        oaro.feed_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"]
        * solute_multiplier
        / sep.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"]
    )
    if split_fraction_treat <= 1:
        sep.split_fraction[:, "treat"].value = (
            oaro.feed_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"]
            * solute_multiplier
            / sep.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"]
        )
    else:
        sep.split_fraction[:, "treat"].value = 1


def makeup_initializer(mixer, oaro, solvent_multiplier):
    makeup_flow_mass = value(
        oaro.feed_inlet.flow_mass_phase_comp[0, "Liq", "H2O"] * solvent_multiplier
        - mixer.treat.flow_mass_phase_comp[0, "Liq", "H2O"]
    )
    if makeup_flow_mass >= 0:
        mixer.makeup_state[0].flow_mass_phase_comp["Liq", "H2O"].value = (
            oaro.feed_inlet.flow_mass_phase_comp[0, "Liq", "H2O"] * solvent_multiplier
            - mixer.treat.flow_mass_phase_comp[0, "Liq", "H2O"]
        )
    else:
        mixer.makeup_state[0].flow_mass_phase_comp["Liq", "H2O"].value = 1e-9
    mixer.makeup_state[0].flow_mass_phase_comp["Liq", "NaCl"].value = 0


def recycle_pump_initializer(pump, oaro, solvent_multiplier, solute_multiplier):

    feed_temperature = 273.15 + 25
    pump.control_volume.properties_out[0].temperature.value = feed_temperature
    pump.control_volume.properties_out[0].flow_mass_phase_comp["Liq", "H2O"].value = (
        oaro.feed_inlet.flow_mass_phase_comp[0, "Liq", "H2O"] * solvent_multiplier
    )
    pump.control_volume.properties_out[0].flow_mass_phase_comp["Liq", "NaCl"].value = (
        oaro.feed_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"] * solute_multiplier
    )


def solve(model, solver=None, tee=False, raise_on_failure=False):
    # ---solving---
    if solver is None:
        solver = get_solver(options={"bound_push": 1e-2})

    results = solver.solve(model, tee=tee)
    if check_optimal_termination(results):
        return results
    msg = (
        "The current configuration is infeasible. Please adjust the decision variables."
    )
    if raise_on_failure:
        raise RuntimeError(msg)
    else:
        print(msg)
        return results


def initialize_loop(m, solvent_multiplier, solute_multiplier):

    for stage in m.fs.IntermediateStages:
        propagate_state(m.fs.OARO_to_pump[stage - 1])
        m.fs.PrimaryPumps[stage].initialize()

        # ---initialize loop---
        propagate_state(m.fs.pump_to_OARO[stage])
        recycle_pump_initializer(
            m.fs.RecyclePumps[stage + 1],
            m.fs.OAROUnits[stage],
            solvent_multiplier=solvent_multiplier,
            solute_multiplier=solute_multiplier,
        )
        propagate_state(m.fs.recyclepump_to_OARO[stage + 1])
        m.fs.OAROUnits[stage].initialize()

        propagate_state(m.fs.OARO_to_ERD[stage])
        m.fs.EnergyRecoveryDevices[stage].initialize()
        propagate_state(m.fs.ERD_to_separator[stage])
        purge_initializer(
            m.fs.Separators[stage],
            m.fs.OAROUnits[stage - 1],
            solute_multiplier=solute_multiplier,
        )
        m.fs.Separators[stage].initialize()

        propagate_state(m.fs.separator_to_mixer[stage])
        makeup_initializer(
            m.fs.Mixers[stage],
            m.fs.OAROUnits[stage - 1],
            solvent_multiplier=solvent_multiplier,
        )
        m.fs.Mixers[stage].initialize()
        propagate_state(m.fs.mixer_to_recyclepump[stage])
        m.fs.RecyclePumps[stage].initialize()

        propagate_state(m.fs.recyclepump_to_OARO[stage])
        propagate_state(m.fs.pump_to_OARO[stage - 1])
        m.fs.OAROUnits[stage - 1].initialize()


def initialize_system(
    m,
    number_of_stages=None,
    solvent_multiplier=None,
    solute_multiplier=None,
    solver=None,
):
    if solver is None:
        solver = get_solver()

    last_stage = m.fs.LastStage
    first_stage = m.fs.FirstStage

    # ---initialize feed block---
    m.fs.feed.initialize()

    propagate_state(m.fs.feed_to_pump)
    m.fs.PrimaryPumps[first_stage].initialize()

    # ---initialize first OARO unit---
    if number_of_stages > 1:
        propagate_state(m.fs.pump_to_OARO[first_stage])
        recycle_pump_initializer(
            m.fs.RecyclePumps[first_stage + 1],
            m.fs.OAROUnits[first_stage],
            solvent_multiplier=solvent_multiplier,
            solute_multiplier=solute_multiplier,
        )
        propagate_state(m.fs.recyclepump_to_OARO[first_stage + 1])
        m.fs.OAROUnits[first_stage].initialize()

        if number_of_stages > 2:
            initialize_loop(m, solvent_multiplier, solute_multiplier)

        # ---initialize RO loop---
        propagate_state(m.fs.OARO_to_pump[last_stage - 1])
        m.fs.PrimaryPumps[last_stage].initialize()

    propagate_state(m.fs.pump_to_ro)
    m.fs.RO.initialize()

    propagate_state(m.fs.ro_to_ERD)
    m.fs.EnergyRecoveryDevices[last_stage].initialize()

    if number_of_stages > 1:
        purge_initializer(
            m.fs.Separators[last_stage],
            m.fs.OAROUnits[last_stage - 1],
            solute_multiplier=solute_multiplier,
        )
        propagate_state(m.fs.ERD_to_separator[last_stage])
        m.fs.Separators[last_stage].initialize()
        propagate_state(m.fs.separator_to_mixer[last_stage])
        makeup_initializer(
            m.fs.Mixers[last_stage],
            m.fs.OAROUnits[last_stage - 1],
            solvent_multiplier=solvent_multiplier,
        )
        m.fs.Mixers[last_stage].initialize()
        propagate_state(m.fs.mixer_to_recyclepump[last_stage])
        m.fs.RecyclePumps[last_stage].initialize()

        propagate_state(m.fs.recyclepump_to_OARO[last_stage])
        propagate_state(m.fs.pump_to_OARO[last_stage - 1])
        m.fs.OAROUnits[last_stage - 1].initialize()

        # ---initialize first ERD---
        propagate_state(m.fs.OARO_to_ERD[first_stage])
        m.fs.EnergyRecoveryDevices[first_stage].initialize()

    m.fs.costing.initialize()


def optimize_set_up(
    m,
    number_of_stages=None,
    water_recovery=None,
):
    # add objective
    m.fs.objective = Objective(expr=m.fs.costing.LCOW)

    # Primary pumps
    for idx, pump in m.fs.PrimaryPumps.items():
        pump.control_volume.properties_out[0].pressure.unfix()
        pump.deltaP.setlb(0)
        if idx == m.fs.Stages.last():
            pump.max_ro_pressure_con = Constraint(
                expr=(
                    m.fs.ro_min_pressure,
                    pump.control_volume.properties_out[0].pressure,
                    m.fs.ro_max_pressure,
                )
            )
        else:
            pump.max_oaro_pressure_con = Constraint(
                expr=(
                    m.fs.oaro_min_pressure,
                    pump.control_volume.properties_out[0].pressure,
                    m.fs.oaro_max_pressure,
                )
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

    # OARO Units
    for stage in m.fs.OAROUnits.values():
        stage.recovery_mass_phase_comp[0, "Liq", "H2O"].unfix()

        stage.area.unfix()
        stage.area.setlb(1)
        stage.area.setub(20000)

        stage.width.unfix()
        stage.width.setlb(0.1)
        stage.width.setub(1000)

        stage.feed_side.velocity[0, 0].unfix()
        stage.feed_side.velocity[0, 0].setlb(0)
        stage.feed_side.velocity[0, 0].setub(1)

        stage.feed_side.N_Re[0, 0].setlb(100)
        stage.feed_side.N_Re[0, 0].setub(2000)

        stage.feed_side.N_Re[0, 1].setlb(100)
        stage.feed_side.N_Re[0, 1].setub(2000)

        stage.permeate_side.N_Re[0, 0].setlb(100)
        stage.permeate_side.N_Re[0, 0].setub(2000)

        stage.permeate_side.N_Re[0, 1].setlb(100)
        stage.permeate_side.N_Re[0, 1].setub(2000)

        stage.permeate_outlet.pressure[0].unfix()
        stage.permeate_outlet.pressure[0].setlb(101325)

        stage.oaro_avg_water_flux_con = Constraint(
            expr=(
                0.1,
                pyunits.convert(
                    stage.flux_mass_phase_comp_avg[0, "Liq", "H2O"]
                    / stage.dens_solvent,
                    to_units=pyunits.L / pyunits.m**2 / pyunits.hr,
                ),
                10,
            )
        )

        stage.oaro_feed_water_flux_con = Constraint(
            expr=(
                0.1 / 5,
                pyunits.convert(
                    stage.flux_mass_phase_comp[0, 0, "Liq", "H2O"] / stage.dens_solvent,
                    to_units=pyunits.L / pyunits.m**2 / pyunits.hr,
                ),
                10 * 1.5,
            )
        )

        stage.oaro_permeate_water_flux_con = Constraint(
            expr=(
                0.1 / 5,
                pyunits.convert(
                    stage.flux_mass_phase_comp[0, 1, "Liq", "H2O"] / stage.dens_solvent,
                    to_units=pyunits.L / pyunits.m**2 / pyunits.hr,
                ),
                10 * 1.5,
            )
        )

        stage.oaro_avg_salt_flux_con = Constraint(
            expr=(
                0,
                pyunits.convert(
                    stage.flux_mass_phase_comp_avg[0, "Liq", "NaCl"],
                    to_units=pyunits.g / pyunits.m**2 / pyunits.hr,
                ),
                50,
            )
        )

        stage.oaro_inlet_salt_flux_con = Constraint(
            expr=(
                0,
                pyunits.convert(
                    stage.flux_mass_phase_comp[0, 0, "Liq", "NaCl"],
                    to_units=pyunits.g / pyunits.m**2 / pyunits.hr,
                ),
                50,
            )
        )

        stage.oaro_permeate_salt_flux_con = Constraint(
            expr=(
                0,
                pyunits.convert(
                    stage.flux_mass_phase_comp[0, 1, "Liq", "NaCl"],
                    to_units=pyunits.g / pyunits.m**2 / pyunits.hr,
                ),
                50,
            )
        )

        stage.min_permeate_flow_rate_con = Constraint(
            expr=sum(
                stage.permeate_inlet.flow_mass_phase_comp[0, "Liq", j]
                for j in ["H2O", "NaCl"]
            )
            >= 0.15
            * sum(
                stage.feed_inlet.flow_mass_phase_comp[0, "Liq", j]
                for j in ["H2O", "NaCl"]
            )
        )

        stage.max_permeate_flow_rate_con = Constraint(
            expr=sum(
                stage.permeate_inlet.flow_mass_phase_comp[0, "Liq", j]
                for j in ["H2O", "NaCl"]
            )
            <= 0.8
            * sum(
                stage.feed_inlet.flow_mass_phase_comp[0, "Liq", j]
                for j in ["H2O", "NaCl"]
            )
        )

    # RO
    m.fs.RO.recovery_mass_phase_comp[0, "Liq", "H2O"].unfix()

    m.fs.RO.area.setlb(1)
    m.fs.RO.area.setub(20000)

    m.fs.RO.width.unfix()
    m.fs.RO.width.setlb(0.1)
    m.fs.RO.width.setub(1000)

    m.fs.RO.feed_side.N_Re[0, 0].setlb(100)
    m.fs.RO.feed_side.N_Re[0, 0].setub(2000)

    m.fs.RO.feed_side.N_Re[0, 1].setlb(100)
    m.fs.RO.feed_side.N_Re[0, 1].setub(2000)

    m.fs.RO.feed_side.properties[0, 0].conc_mass_phase_comp["Liq", "NaCl"].setlb(
        10 * pyunits.g / pyunits.L
    )
    m.fs.RO.feed_side.properties[0, 1].conc_mass_phase_comp["Liq", "NaCl"].setlb(
        10 * pyunits.g / pyunits.L
    )

    m.fs.RO.ro_avg_water_flux_con = Constraint(
        expr=(
            0.4,
            pyunits.convert(
                m.fs.RO.flux_mass_phase_comp_avg[0, "Liq", "H2O"]
                / m.fs.RO.dens_solvent,
                to_units=pyunits.L / pyunits.m**2 / pyunits.hr,
            ),
            40,
        )
    )

    m.fs.RO.ro_feed_water_flux_con = Constraint(
        expr=(
            0.4 / 5,
            pyunits.convert(
                m.fs.RO.flux_mass_phase_comp[0.0, 0.0, "Liq", "H2O"]
                / m.fs.RO.dens_solvent,
                to_units=pyunits.L / pyunits.m**2 / pyunits.hr,
            ),
            40 * 1.5,
        )
    )

    m.fs.RO.ro_permeate_water_flux_con = Constraint(
        expr=(
            0.4 / 5,
            pyunits.convert(
                m.fs.RO.flux_mass_phase_comp[0.0, 1.0, "Liq", "H2O"]
                / m.fs.RO.dens_solvent,
                to_units=pyunits.L / pyunits.m**2 / pyunits.hr,
            ),
            40 * 1.5,
        )
    )

    m.fs.RO.ro_avg_salt_flux_con = Constraint(
        expr=(
            0,
            pyunits.convert(
                m.fs.RO.flux_mass_phase_comp_avg[0, "Liq", "NaCl"],
                to_units=pyunits.g / pyunits.m**2 / pyunits.hr,
            ),
            50,
        )
    )

    m.fs.RO.ro_feed_salt_flux_con = Constraint(
        expr=(
            0,
            pyunits.convert(
                m.fs.RO.flux_mass_phase_comp[0, 0, "Liq", "NaCl"],
                to_units=pyunits.g / pyunits.m**2 / pyunits.hr,
            ),
            50,
        )
    )

    m.fs.RO.ro_retentate_salt_flux_con = Constraint(
        expr=(
            0,
            pyunits.convert(
                m.fs.RO.flux_mass_phase_comp[0, 1, "Liq", "NaCl"],
                to_units=pyunits.g / pyunits.m**2 / pyunits.hr,
            ),
            50,
        )
    )

    # Separator
    for sep in m.fs.Separators.values():
        sep.split_fraction[:, "treat"].unfix()
        sep.split_fraction[:, "treat"].setlb(0)
        sep.split_fraction[:, "treat"].setub(1)

    # Makeup
    if number_of_stages > 1:
        for i in range(2, number_of_stages + 1):
            m.fs.ProductSeparator.split_fraction[:, "makeup" + str(i)].unfix()
            m.fs.ProductSeparator.split_fraction[:, "makeup" + str(i)].setlb(0)
            m.fs.ProductSeparator.split_fraction[:, "makeup" + str(i)].setub(1)

    # additional specifications
    m.fs.product_salinity = Param(
        initialize=500e-6, mutable=True
    )  # product NaCl mass fraction [-]

    # additional constraints
    if water_recovery is not None:
        # product mass flow rate fraction of feed [-]
        m.fs.mass_water_recovery.fix(water_recovery)

    m.fs.eq_product_quality = Constraint(
        expr=m.fs.product.properties[0].mass_frac_phase_comp["Liq", "NaCl"]
        <= m.fs.product_salinity
    )
    iscale.constraint_scaling_transform(
        m.fs.eq_product_quality, 1e3
    )  # scaling constraint

    if number_of_stages > 1:
        m.fs.purge_rate_con = Constraint(
            expr=(
                0,
                sum(
                    m.fs.Separators[stage].purge.flow_mass_phase_comp[0, "Liq", j]
                    for stage in m.fs.NonFirstStages
                    for j in ["H2O", "NaCl"]
                )
                / (
                    sum(
                        m.fs.OAROUnits[1].permeate_outlet.flow_mass_phase_comp[
                            0, "Liq", j
                        ]
                        for j in ["H2O", "NaCl"]
                    )
                    - sum(
                        m.fs.OAROUnits[1].permeate_inlet.flow_mass_phase_comp[
                            0, "Liq", j
                        ]
                        for j in ["H2O", "NaCl"]
                    )
                ),
                0.2,
            )
        )

        assert_degrees_of_freedom(
            m,
            6 * (m.fs.NumberOfStages - 1)
            + 3
            - (1 if (water_recovery is not None) else 0),
        )


def display_system(m):
    print("----system metrics----")
    feed_flow_vol = m.fs.feed.properties[0].flow_vol_phase["Liq"].value * 3600
    feed_flow_mass = sum(
        m.fs.feed.flow_mass_phase_comp[0, "Liq", j].value for j in ["H2O", "NaCl"]
    )
    feed_mass_frac_NaCl = (
        m.fs.feed.flow_mass_phase_comp[0, "Liq", "NaCl"].value / feed_flow_mass
    )
    print(
        "Feed: %.2f m3/h, %.2f kg/s, %.0f ppm"
        % (feed_flow_vol, feed_flow_mass, feed_mass_frac_NaCl * 1e6)
    )

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
    print("Mass water recovery: %.1f%%" % (value(m.fs.mass_water_recovery) * 100))
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
            "OARO Stage %d water recovery: %.3f"
            % (
                stage,
                value(m.fs.OAROUnits[stage].recovery_mass_phase_comp[0, "Liq", "H2O"]),
            )
        )
        print(
            "OARO Stage %d average water flux: %.1f L/m2/h"
            % (
                stage,
                value(m.fs.OAROUnits[stage].flux_mass_phase_comp_avg[0, "Liq", "H2O"])
                / 1e3
                * 1000
                * 3600,
            )
        )
        print(
            "OARO Stage %d average salt flux: %.1f g/m2/h"
            % (
                stage,
                value(m.fs.OAROUnits[stage].flux_mass_phase_comp_avg[0, "Liq", "NaCl"])
                * 1000
                * 3600,
            )
        )
        print(
            "OARO Stage %d feed operating pressure %.1f bar"
            % (stage, m.fs.OAROUnits[stage].feed_inlet.pressure[0].value / 1e5)
        )
        print(
            "OARO Stage %d permeate operating pressure %.1f bar"
            % (stage, m.fs.OAROUnits[stage].permeate_inlet.pressure[0].value / 1e5)
        )
        print(
            "OARO Stage %d membrane area      %.1f m2"
            % (stage, m.fs.OAROUnits[stage].area.value)
        )
        print(
            "OARO Stage %d water perm. coeff.  %.3f LMH/bar"
            % (stage, m.fs.OAROUnits[stage].A_comp[0, "H2O"].value * (3.6e11))
        )
        print(
            "OARO Stage %d salt perm. coeff.  %.3f LMH/bar"
            % (stage, m.fs.OAROUnits[stage].B_comp[0, "NaCl"].value * (1000.0 * 3600.0))
        )
        print(
            "OARO Stage %d feed-side velocity.  %.3f"
            % (stage, m.fs.OAROUnits[stage].feed_side.velocity[0, 0].value)
        )
        print(
            "OARO Stage %d feed-side Reynolds number.  %.3f"
            % (stage, m.fs.OAROUnits[stage].feed_side.N_Re[0, 0].value)
        )
    print(
        "RO water recovery: %.3f"
        % (value(m.fs.RO.recovery_mass_phase_comp[0, "Liq", "H2O"]))
    )
    print(
        "RO feed operating pressure %.1f bar" % (m.fs.RO.inlet.pressure[0].value / 1e5)
    )
    print(
        "RO permeate operating pressure %.1f bar"
        % (m.fs.RO.permeate.pressure[0].value / 1e5)
    )
    print("RO membrane area      %.1f m2" % (m.fs.RO.area.value))
    print("RO membrane width      %.1f m" % (m.fs.RO.width.value))
    print(
        "RO water perm. coeff.  %.3f LMH/bar"
        % (m.fs.RO.A_comp[0, "H2O"].value * (3.6e11))
    )
    print(
        "RO salt perm. coeff.  %.3f LMH"
        % (m.fs.RO.B_comp[0, "NaCl"].value * (1000.0 * 3600.0))
    )


def display_state(m):
    print("--------state---------")

    def print_state(s, b):
        feed_flow_mass = sum(
            m.fs.feed.flow_mass_phase_comp[0, "Liq", j].value for j in ["H2O", "NaCl"]
        )
        flow_mass = sum(
            b.flow_mass_phase_comp[0, "Liq", j].value for j in ["H2O", "NaCl"]
        )
        normalized_flow_mass = flow_mass / feed_flow_mass * 100
        mass_frac_ppm = b.flow_mass_phase_comp[0, "Liq", "NaCl"].value / flow_mass * 1e3
        pressure_bar = b.pressure[0].value / 1e5
        print(
            s.ljust(20)
            + ": %.0f, %.3f g/L, %.1f bar, %.4f kg/s, %.4f kg/s"
            % (
                normalized_flow_mass,
                mass_frac_ppm,
                pressure_bar,
                b.flow_mass_phase_comp[0, "Liq", "H2O"].value,
                b.flow_mass_phase_comp[0, "Liq", "NaCl"].value,
            )
        )

    print_state("Feed", m.fs.feed.outlet)

    for stage in m.fs.Stages:

        print_state(f"Primary Pump {stage} out", m.fs.PrimaryPumps[stage].outlet)
        print_state(
            f"ERD {stage} out",
            m.fs.EnergyRecoveryDevices[stage].outlet,
        )

        if stage == m.fs.LastStage:
            pass
        else:
            print_state(f"OARO {stage} feed inlet", m.fs.OAROUnits[stage].feed_inlet)
            print_state(f"OARO {stage} feed outlet", m.fs.OAROUnits[stage].feed_outlet)
            print_state(
                f"OARO {stage} permeate inlet", m.fs.OAROUnits[stage].permeate_inlet
            )
            print_state(
                f"OARO {stage} permeate outlet", m.fs.OAROUnits[stage].permeate_outlet
            )

        if stage == m.fs.FirstStage:
            pass
        else:
            print_state(f"Recycle Pump {stage} out", m.fs.RecyclePumps[stage].outlet)
            print_state(f"Purge {stage}", m.fs.Separators[stage].purge)
            print_state(f"Make-up {stage}", m.fs.Mixers[stage].makeup)

    print_state(f"RO inlet", m.fs.RO.inlet)
    print_state(f"RO permeate", m.fs.RO.permeate)
    print_state(f"RO retentate", m.fs.RO.retentate)

    print_state(f"Disposal", m.fs.disposal.inlet)
    print_state(f"Product", m.fs.product.inlet)


if __name__ == "__main__":
    m = main(number_of_stages=3, system_recovery=0.5, erd_type=ERDtype.pump_as_turbine)
