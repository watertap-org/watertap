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
    Block,
    NonNegativeReals,
    RangeSet,
    Set,
    check_optimal_termination,
    units as pyunits,
)
from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.solvers import get_solver
from idaes.core.util.exceptions import InitializationError
from idaes.core.util.initialization import propagate_state
from idaes.core.util.misc import StrEnum
from idaes.models.unit_models import Feed, Product, Mixer
from idaes.models.unit_models.mixer import MomentumMixingType
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslogger

from watertap.unit_models.reverse_osmosis_1D import (
    ReverseOsmosis1D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
from watertap.unit_models.pressure_changer import Pump, EnergyRecoveryDevice
from watertap.core.util.initialization import (
    assert_no_degrees_of_freedom,
    assert_degrees_of_freedom,
)
from watertap.costing.watertap_costing_package import (
    WaterTAPCosting,
    make_capital_cost_var,
)
import watertap.property_models.NaCl_prop_pack as props


class ACase(StrEnum):
    fixed = "fixed"
    optimize = "optimize"
    single_optimum = "single_optimum"


class BCase(StrEnum):
    single_optimum = "single_optimum"
    optimize = "optimize"


class ABTradeoff(StrEnum):
    inequality_constraint = "inequality_constraint"
    equality_constraint = "equality_constraint"
    none = "none"


def run_lsrro_case(
    number_of_stages,
    water_recovery=None,
    Cin=None,
    Cbrine=None,
    A_case=ACase.fixed,
    B_case=BCase.optimize,
    AB_tradeoff=ABTradeoff.none,
    A_value=None,
    has_NaCl_solubility_limit=None,
    has_calculated_concentration_polarization=None,
    has_calculated_ro_pressure_drop=None,
    permeate_quality_limit=None,
    AB_gamma_factor=None,
    B_max=None,
    number_of_RO_finite_elements=10,
):
    m = build(
        number_of_stages,
        has_NaCl_solubility_limit,
        has_calculated_concentration_polarization,
        has_calculated_ro_pressure_drop,
        number_of_RO_finite_elements,
        B_max,
    )
    set_operating_conditions(m, Cin)

    initialize(m)
    solve(m)
    print("\n***---Simulation results---***")
    display_system(m)
    display_design(m)
    display_state(m)

    optimize_set_up(
        m,
        water_recovery,
        Cbrine,
        A_case,
        B_case,
        AB_tradeoff,
        A_value,
        permeate_quality_limit,
        AB_gamma_factor,
        B_max,
    )
    res = solve(m, raise_on_failure=False, tee=False)
    print("\n***---Optimization results---***")
    if check_optimal_termination(res):
        display_system(m)
        display_design(m)
        display_state(m)
        display_RO_reports(m)

    return m, res


def build(
    number_of_stages=2,
    has_NaCl_solubility_limit=True,
    has_calculated_concentration_polarization=True,
    has_calculated_ro_pressure_drop=True,
    number_of_RO_finite_elements=10,
    B_max=None,
):
    # ---building model---
    m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.NaClParameterBlock()
    m.fs.costing = WaterTAPCosting()

    # explicitly set the costing parameters used
    m.fs.costing.utilization_factor.fix(0.9)
    m.fs.costing.factor_total_investment.fix(2)
    m.fs.costing.factor_maintenance_labor_chemical.fix(0.03)
    m.fs.costing.factor_capital_annualization.fix(0.1)
    m.fs.costing.electricity_base_cost.set_value(0.07)
    m.fs.costing.reverse_osmosis.factor_membrane_replacement.fix(0.15)
    m.fs.costing.reverse_osmosis.membrane_cost.fix(30)
    m.fs.costing.reverse_osmosis.high_pressure_membrane_cost.fix(50)
    m.fs.costing.high_pressure_pump.cost.fix(53 / 1e5 * 3600)
    m.fs.costing.energy_recovery_device.pressure_exchanger_cost.fix(535)

    m.fs.NumberOfStages = Param(initialize=number_of_stages)
    m.fs.Stages = RangeSet(m.fs.NumberOfStages)
    m.fs.LSRRO_Stages = RangeSet(2, m.fs.NumberOfStages)
    m.fs.NonFinalStages = RangeSet(m.fs.NumberOfStages - 1)
    if number_of_stages > 1:
        m.fs.IntermediateStages = RangeSet(2, m.fs.NumberOfStages - 1)
    else:
        m.fs.IntermediateStages = RangeSet(0)
    m.fs.FirstStage = m.fs.Stages.first()
    m.fs.LastStage = m.fs.Stages.last()

    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)
    m.fs.disposal = Product(property_package=m.fs.properties)

    # Add the mixers
    m.fs.Mixers = Mixer(
        m.fs.NonFinalStages,
        property_package=m.fs.properties,
        momentum_mixing_type=MomentumMixingType.equality,
        inlet_list=["upstream", "downstream"],
    )

    # Add the pumps
    m.fs.PrimaryPumps = Pump(m.fs.Stages, property_package=m.fs.properties)
    for pump in m.fs.PrimaryPumps.values():
        pump.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=cost_high_pressure_pump_lsrro,
        )

    # Add the equalizer pumps
    m.fs.BoosterPumps = Pump(m.fs.LSRRO_Stages, property_package=m.fs.properties)
    for pump in m.fs.BoosterPumps.values():
        pump.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=cost_high_pressure_pump_lsrro,
        )

    m.fs.total_pump_work = Expression(
        expr=sum(
            pyunits.convert(pump.work_mechanical[0], to_units=pyunits.kW)
            for pump in itertools.chain(
                m.fs.PrimaryPumps.values(), m.fs.BoosterPumps.values()
            )
        )
    )

    # Add the stages ROs
    if has_calculated_ro_pressure_drop:
        pressure_change_type = PressureChangeType.calculated
    else:
        pressure_change_type = PressureChangeType.fixed_per_stage

    if has_calculated_concentration_polarization:
        cp_type = ConcentrationPolarizationType.calculated
        kf_type = MassTransferCoefficient.calculated
    else:
        cp_type = ConcentrationPolarizationType.none
        kf_type = MassTransferCoefficient.none

    m.fs.ROUnits = ReverseOsmosis1D(
        m.fs.Stages,
        property_package=m.fs.properties,
        has_pressure_change=has_calculated_ro_pressure_drop,
        pressure_change_type=pressure_change_type,
        mass_transfer_coefficient=kf_type,
        concentration_polarization_type=cp_type,
        transformation_scheme="BACKWARD",
        transformation_method="dae.finite_difference",
        finite_elements=number_of_RO_finite_elements,
        has_full_reporting=True,
    )

    for idx, ro_stage in m.fs.ROUnits.items():
        if idx == m.fs.FirstStage:
            ro_stage.costing = UnitModelCostingBlock(
                flowsheet_costing_block=m.fs.costing,
                costing_method_arguments={"ro_type": "standard"},
            )
        else:
            ro_stage.costing = UnitModelCostingBlock(
                flowsheet_costing_block=m.fs.costing,
                costing_method_arguments={"ro_type": "high_pressure"},
            )

    # Add EnergyRecoveryDevices
    m.fs.EnergyRecoveryDeviceSet = Set(
        initialize=[m.fs.FirstStage, m.fs.LastStage]
        if m.fs.FirstStage < m.fs.LastStage
        else [m.fs.LastStage]
    )
    m.fs.EnergyRecoveryDevices = EnergyRecoveryDevice(
        m.fs.EnergyRecoveryDeviceSet, property_package=m.fs.properties
    )
    for erd in m.fs.EnergyRecoveryDevices.values():
        erd.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method_arguments={
                "energy_recovery_device_type": "pressure_exchanger"
            },
        )

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
    m.fs.ro_min_pressure = Param(initialize=10e5, units=pyunits.Pa, mutable=True)
    m.fs.lsrro_min_pressure = Param(initialize=10e5, units=pyunits.Pa, mutable=True)
    m.fs.ro_max_pressure = Param(initialize=85e5, units=pyunits.Pa, mutable=True)
    m.fs.lsrro_max_pressure = Param(initialize=65e5, units=pyunits.Pa, mutable=True)

    if B_max is not None:
        m.fs.B_max = Param(initialize=B_max, mutable=True, units=pyunits.m / pyunits.s)

    # system water recovery
    m.fs.water_recovery = Var(
        initialize=0.5,
        bounds=(0, 1),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="System Water Recovery",
    )
    m.fs.eq_water_recovery = Constraint(
        expr=m.fs.feed.properties[0].flow_vol * m.fs.water_recovery
        == m.fs.product.properties[0].flow_vol
    )

    # costing and summary quantities
    m.fs.costing.cost_process()

    product_flow_vol_total = m.fs.product.properties[0].flow_vol
    m.fs.costing.add_annual_water_production(product_flow_vol_total)
    m.fs.costing.add_specific_energy_consumption(product_flow_vol_total)
    m.fs.costing.add_LCOW(product_flow_vol_total)

    # objective
    m.fs.objective = Objective(expr=m.fs.costing.LCOW)

    # Expressions for parameter sweep -----------------------------------------
    # Final permeate concentration as mass fraction
    m.fs.product.properties[0].mass_frac_phase_comp

    # Touch feed concentration as mass concentration
    m.fs.feed.properties[0].conc_mass_phase_comp

    # Touch final brine concentration as mass concentration
    m.fs.disposal.properties[0].conc_mass_phase_comp

    # Touch final brine concentration as mass fraction
    m.fs.disposal.properties[0].mass_frac_phase_comp

    m.fs.mass_water_recovery = Expression(
        expr=m.fs.product.flow_mass_phase_comp[0, "Liq", "H2O"]
        / m.fs.feed.flow_mass_phase_comp[0, "Liq", "H2O"]
    )

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

    @m.fs.Expression(m.fs.Stages)
    def stage_recovery_vol(fs, stage):
        return (
            fs.ROUnits[stage].mixed_permeate[0].flow_vol
            / fs.PrimaryPumps[stage].control_volume.properties_in[0].flow_vol
        )

    @m.fs.Expression(m.fs.Stages)
    def stage_recovery_mass_H2O(fs, stage):
        return (
            fs.ROUnits[stage].mixed_permeate[0].flow_mass_phase_comp["Liq", "H2O"]
            / m.fs.PrimaryPumps[stage]
            .control_volume.properties_in[0]
            .flow_mass_phase_comp["Liq", "H2O"]
        )

    m.fs.costing.primary_pump_capex_lcow = Expression(
        expr=m.fs.costing.factor_capital_annualization
        * sum(m.fs.PrimaryPumps[n].costing.capital_cost for n in m.fs.Stages)
        / m.fs.costing.annual_water_production
    )

    m.fs.total_membrane_area = Expression(
        expr=sum(ro.area for ro in m.fs.ROUnits.values())
    )

    m.fs.costing.booster_pump_capex_lcow = Expression(
        expr=m.fs.costing.factor_capital_annualization
        * (
            sum(m.fs.BoosterPumps[n].costing.capital_cost for n in m.fs.LSRRO_Stages)
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
                m.fs.costing.booster_pump_capex_lcow
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
        * sum(m.fs.ROUnits[n].costing.capital_cost for n in m.fs.Stages)
        / m.fs.costing.annual_water_production
    )

    m.fs.costing.indirect_capex_lcow = Expression(
        expr=m.fs.costing.factor_capital_annualization
        * (m.fs.costing.total_investment_cost - m.fs.costing.total_capital_cost)
        / m.fs.costing.annual_water_production
    )

    m.fs.costing.membrane_replacement_lcow = Expression(
        expr=sum(m.fs.ROUnits[n].costing.fixed_operating_cost for n in m.fs.Stages)
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

    # Connections ------------------------------------------------------------
    # Connect the feed to the first pump
    m.fs.feed_to_pump = Arc(
        source=m.fs.feed.outlet, destination=m.fs.PrimaryPumps[1].inlet
    )

    # Connect the primary RO permeate to the product
    m.fs.primary_RO_to_product = Arc(
        source=m.fs.ROUnits[1].permeate, destination=m.fs.product.inlet
    )

    # Connect the Pump n to the Mixer n
    m.fs.pump_to_mixer = Arc(
        m.fs.NonFinalStages,
        rule=lambda fs, n: {
            "source": fs.PrimaryPumps[n].outlet,
            "destination": fs.Mixers[n].upstream,
        },
    )

    # Connect the Mixer n to the Stage n
    m.fs.mixer_to_stage = Arc(
        m.fs.NonFinalStages,
        rule=lambda fs, n: {
            "source": fs.Mixers[n].outlet,
            "destination": fs.ROUnits[n].inlet,
        },
    )

    # Connect the Stage n to the Eq Pump n
    m.fs.stage_permeate_to_booster_pump = Arc(
        m.fs.LSRRO_Stages,
        rule=lambda fs, n: {
            "source": fs.ROUnits[n].permeate,
            "destination": fs.BoosterPumps[n].inlet,
        },
    )

    # Connect the Eq Pump n to the Mixer n-1
    m.fs.booster_pump_to_mixer = Arc(
        m.fs.LSRRO_Stages,
        rule=lambda fs, n: {
            "source": fs.BoosterPumps[n].outlet,
            "destination": fs.Mixers[n - 1].downstream,
        },
    )

    last_stage = m.fs.LastStage
    if number_of_stages > 1:
        # Connect the primary RO permeate to the product
        m.fs.primary_RO_to_erd = Arc(
            source=m.fs.ROUnits[1].retentate,
            destination=m.fs.EnergyRecoveryDevices[1].inlet,
        )
        # Connect 1st stage ERD to primary pump
        m.fs.primary_ERD_to_pump = Arc(
            source=m.fs.EnergyRecoveryDevices[1].outlet,
            destination=m.fs.PrimaryPumps[2].inlet,
        )

    # Connect the Stage n to the Pump n+1
    m.fs.stage_retentate_to_pump = Arc(
        m.fs.IntermediateStages,
        rule=lambda fs, n: {
            "source": fs.ROUnits[n].retentate,
            "destination": fs.PrimaryPumps[n + 1].inlet,
        },
    )
    # Connect the Pump N to the Stage N
    m.fs.pumpN_to_stageN = Arc(
        source=m.fs.PrimaryPumps[last_stage].outlet,
        destination=m.fs.ROUnits[last_stage].inlet,
    )
    # Connect Final Stage to EnergyRecoveryDevice Pump
    m.fs.stage_to_erd = Arc(
        source=m.fs.ROUnits[last_stage].retentate,
        destination=m.fs.EnergyRecoveryDevices[last_stage].inlet,
    )
    # Connect the EnergyRecoveryDevice to the disposal
    m.fs.erd_to_disposal = Arc(
        source=m.fs.EnergyRecoveryDevices[last_stage].outlet,
        destination=m.fs.disposal.inlet,
    )

    # additional bounding
    if has_NaCl_solubility_limit:
        for b in m.component_data_objects(Block, descend_into=True):
            # NaCl solubility limit
            if hasattr(b, "is_property_constructed") and b.is_property_constructed(
                "mass_frac_phase_comp"
            ):
                b.mass_frac_phase_comp["Liq", "NaCl"].setub(0.2614)

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def cost_high_pressure_pump_lsrro(blk, cost_electricity_flow=True):
    t0 = blk.flowsheet().time.first()
    make_capital_cost_var(blk)
    blk.capital_cost_constraint = Constraint(
        expr=blk.capital_cost
        == blk.costing_package.high_pressure_pump.cost
        * pyunits.watt
        / (pyunits.m**3 * pyunits.pascal / pyunits.s)
        * blk.unit_model.outlet.pressure[t0]
        * blk.unit_model.control_volume.properties_out[t0].flow_vol
    )

    if cost_electricity_flow:
        blk.costing_package.cost_flow(
            pyunits.convert(blk.unit_model.work_mechanical[t0], to_units=pyunits.kW),
            "electricity",
        )


def set_operating_conditions(m, Cin=None):
    # ---specifications---
    # parameters
    pump_efi = 0.75  # pump efficiency [-]
    erd_efi = 0.8  # energy recovery device efficiency [-]
    mem_A = 4.2e-12  # membrane water permeability coefficient [m/s-Pa]
    mem_B = 3.5e-8  # membrane salt permeability coefficient [m/s]
    height = 1e-3  # channel height in membrane stage [m]
    spacer_porosity = 0.85  # spacer porosity in membrane stage [-]
    width = 5  # effective membrane width [m]
    area = 100  # membrane area [m^2]
    pressure_atm = 101325  # atmospheric pressure [Pa]

    # feed
    # feed_flow_mass = 1*pyunits.kg/pyunits.s
    if Cin is None:
        Cin = 70
    feed_temperature = 273.15 + 20

    # initialize feed
    m.fs.feed.pressure[0].fix(pressure_atm)
    m.fs.feed.temperature[0].fix(feed_temperature)
    m.fs.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"] = Cin
    m.fs.feed.properties.calculate_state(
        var_args={
            ("conc_mass_phase_comp", ("Liq", "NaCl")): value(
                m.fs.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"]
            ),  # feed mass concentration
            ("flow_vol_phase", "Liq"): 1e-3,
        },  # volumetric feed flowrate [-]
        hold_state=True,  # fixes the calculated component mass flow rates
    )

    # initialize pumps
    for pump in m.fs.PrimaryPumps.values():
        pump.control_volume.properties_out[0].pressure = 75e5
        pump.efficiency_pump.fix(pump_efi)
        pump.control_volume.properties_out[0].pressure.fix()
        iscale.set_scaling_factor(pump.control_volume.work, 1e-3)

    # initialize eq pumps
    for pump in m.fs.BoosterPumps.values():
        pump.efficiency_pump.fix(pump_efi)
        iscale.set_scaling_factor(pump.control_volume.work, 1e-3)

    # initialize stages
    for idx, stage in m.fs.ROUnits.items():
        if idx == m.fs.FirstStage:
            B_scale = 1.0
        else:
            B_scale = 100.0

        stage.A_comp.fix(mem_A)
        stage.B_comp.fix(mem_B * B_scale)
        stage.area.fix(area / float(idx))
        stage.width.fix(width)
        stage.mixed_permeate[0].pressure.fix(pressure_atm)
        if (
            stage.config.mass_transfer_coefficient == MassTransferCoefficient.calculated
        ) or stage.config.pressure_change_type == PressureChangeType.calculated:
            stage.channel_height.fix(height)
            stage.spacer_porosity.fix(spacer_porosity)

    # energy recovery devices
    for erd in m.fs.EnergyRecoveryDevices.values():
        erd.efficiency_pump.fix(erd_efi)
        iscale.set_scaling_factor(erd.control_volume.work, 1e-3)

    # if FirstStage *is* LastStage, we'll just overwrite the pressure
    m.fs.EnergyRecoveryDevices[m.fs.FirstStage].control_volume.properties_out[
        0
    ].pressure.fix(70e5)
    m.fs.EnergyRecoveryDevices[m.fs.LastStage].control_volume.properties_out[
        0
    ].pressure.fix(pressure_atm)

    # ---scaling---
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )
    iscale.calculate_scaling_factors(m)

    # ---checking model---
    assert_units_consistent(m)
    assert_no_degrees_of_freedom(m)

    print(
        "Feed Concentration = %.1f ppt"
        % (value(m.fs.feed.flow_mass_phase_comp[0, "Liq", "NaCl"]) * 1000)
    )


def _lsrro_mixer_guess_initializer(
    mixer, solvent_multiplier, solute_multiplier, optarg
):

    for vname in mixer.upstream.vars:
        if vname == "flow_mass_phase_comp":
            for time, phase, comp in mixer.upstream.vars[vname]:
                if comp in mixer.config.property_package.solute_set:
                    mixer.downstream.vars[vname][time, phase, comp].value = (
                        solute_multiplier
                        * mixer.upstream.vars[vname][time, phase, comp].value
                    )
                elif comp in mixer.config.property_package.solvent_set:
                    mixer.downstream.vars[vname][time, phase, comp].value = (
                        solvent_multiplier
                        * mixer.upstream.vars[vname][time, phase, comp].value
                    )
                else:
                    raise RuntimeError(f"Unknown component {comp}")
        else:  # copy the state
            for idx in mixer.upstream.vars[vname]:
                mixer.downstream.vars[vname][idx].value = mixer.upstream.vars[vname][
                    idx
                ].value

    mixer.initialize(optarg=optarg)


def do_forward_initialization_pass(m, optarg, guess_mixers):
    print("--------------------START FORWARD INITIALIZATION PASS--------------------")
    # start with the feed
    m.fs.feed.initialize(optarg=optarg)

    propagate_state(m.fs.feed_to_pump)

    last_stage = m.fs.LastStage
    first_stage = m.fs.FirstStage

    for stage in m.fs.Stages:
        m.fs.PrimaryPumps[stage].initialize(optarg=optarg)

        if stage == last_stage:
            propagate_state(m.fs.pumpN_to_stageN)
        else:
            propagate_state(m.fs.pump_to_mixer[stage])
            if guess_mixers:
                _lsrro_mixer_guess_initializer(
                    m.fs.Mixers[stage],
                    solvent_multiplier=0.5,
                    solute_multiplier=0.2,
                    optarg=optarg,
                )
            else:
                m.fs.Mixers[stage].initialize(optarg=optarg)
            propagate_state(m.fs.mixer_to_stage[stage])

        m.fs.ROUnits[stage].initialize(optarg=optarg, raise_on_failure=False)

        if stage == first_stage:
            propagate_state(m.fs.primary_RO_to_product)
            m.fs.product.initialize(optarg=optarg)
            if value(m.fs.NumberOfStages) > 1:
                propagate_state(m.fs.primary_RO_to_erd)
                m.fs.EnergyRecoveryDevices[first_stage].initialize(optarg=optarg)
                propagate_state(m.fs.primary_ERD_to_pump)
        else:
            propagate_state(m.fs.stage_permeate_to_booster_pump[stage])
            m.fs.BoosterPumps[stage].initialize(optarg=optarg)
            propagate_state(m.fs.booster_pump_to_mixer[stage])

        if stage in m.fs.IntermediateStages:
            propagate_state(m.fs.stage_retentate_to_pump[stage])

    # for the end stage
    propagate_state(m.fs.stage_to_erd)
    m.fs.EnergyRecoveryDevices[last_stage].initialize(optarg=optarg)
    propagate_state(m.fs.erd_to_disposal)
    m.fs.disposal.initialize(optarg=optarg)


def do_backward_initialization_pass(m, optarg):
    print("--------------------START BACKWARD INITIALIZATION PASS--------------------")

    first_stage = m.fs.FirstStage
    for stage in reversed(m.fs.NonFinalStages):
        m.fs.Mixers[stage].initialize(optarg=optarg)
        propagate_state(m.fs.mixer_to_stage[stage])
        m.fs.ROUnits[stage].initialize(optarg=optarg, raise_on_failure=False)
        if stage == first_stage:
            if value(m.fs.NumberOfStages) > 1:
                propagate_state(m.fs.primary_ERD_to_pump)
                m.fs.EnergyRecoveryDevices[first_stage].initialize(optarg=optarg)
                propagate_state(m.fs.primary_RO_to_erd)
            propagate_state(m.fs.primary_RO_to_product)
            m.fs.product.initialize(optarg=optarg)
        else:
            propagate_state(m.fs.stage_retentate_to_pump[stage])
            propagate_state(m.fs.stage_permeate_to_booster_pump[stage])
            m.fs.BoosterPumps[stage].initialize(optarg=optarg)
            propagate_state(m.fs.booster_pump_to_mixer[stage])


def initialize(m, verbose=True, solver=None):

    # ---initializing---
    # set up solvers
    if solver is None:
        solver = get_solver()

    optarg = solver.options
    do_forward_initialization_pass(m, optarg=optarg, guess_mixers=True)
    for _ in range(m.fs.NumberOfStages.value // 2):
        do_backward_initialization_pass(m, optarg=optarg)
        do_forward_initialization_pass(m, optarg=optarg, guess_mixers=False)

    # # set up SD tool
    seq = SequentialDecomposition()
    seq.options.tear_method = "Direct"
    seq.options.iterLim = m.fs.NumberOfStages
    seq.options.tear_set = list(m.fs.booster_pump_to_mixer.values())
    seq.options.log_info = True

    # run SD tool
    def func_initialize(unit):
        outlvl = idaeslogger.INFO if verbose else idaeslogger.CRITICAL
        if "ROUnits" in unit.name:
            unit.initialize(
                optarg=solver.options, outlvl=outlvl, raise_on_failure=False
            )
        else:
            unit.initialize(optarg=solver.options, outlvl=outlvl)

    seq.run(m, func_initialize)

    m.fs.costing.initialize()


def solve(model, solver=None, tee=False, raise_on_failure=False):
    # ---solving---
    if solver is None:
        solver = get_solver()

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


def optimize_set_up(
    m,
    water_recovery=None,
    Cbrine=None,
    A_case=ACase.fixed,
    B_case=BCase.optimize,
    AB_tradeoff=ABTradeoff.none,
    A_value=None,
    permeate_quality_limit=None,
    AB_gamma_factor=None,
    B_max=None,
):
    """
    Get the LSRRO flowsheet ready to optimize

    Attributes
    ----------
    B_case: 'single_optimum' or anything else to optimize B value at every LSR stage

    A_case: 'fixed' or 'optimize' or 'single_optimum' A at every LSR stage

    AB_tradeoff: 'inequality_constraint' B >= function of A
                 'equality_constraint' B = function of A
                 'none' no constraint relating B value to A value

    A_value: if A_case='fixed', then provide a value to fix A with

    Returns
    -------
    model (Pyomo ConcreteModel) : The LSRRO flowsheet.

    """

    for idx, pump in m.fs.PrimaryPumps.items():
        pump.control_volume.properties_out[0].pressure.unfix()
        pump.deltaP.setlb(0)
        if idx > m.fs.Stages.first():
            pump.max_lsrro_pressure_con = Constraint(
                expr=(
                    m.fs.lsrro_min_pressure,
                    pump.control_volume.properties_out[0].pressure,
                    m.fs.lsrro_max_pressure,
                )
            )
            iscale.constraint_scaling_transform(
                pump.max_lsrro_pressure_con,
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

    if value(m.fs.NumberOfStages) > 1:
        m.fs.EnergyRecoveryDevices[1].control_volume.properties_out[0].pressure.unfix()
        m.fs.EnergyRecoveryDevices[1].deltaP.setub(0)

    # unfix eq pumps
    for idx, pump in m.fs.BoosterPumps.items():
        pump.control_volume.properties_out[0].pressure.unfix()
        pump.max_ro_pressure_con = Constraint(
            expr=(
                m.fs.ro_min_pressure,
                pump.control_volume.properties_out[0].pressure,
                m.fs.ro_max_pressure,
            )
        )
        iscale.constraint_scaling_transform(
            pump.max_ro_pressure_con,
            iscale.get_scaling_factor(pump.control_volume.properties_out[0].pressure),
        )
        pump.deltaP.setlb(0)

    if B_case == BCase.single_optimum:
        m.fs.B_comp_system = Var(
            domain=NonNegativeReals,
            units=pyunits.m * pyunits.s**-1,
            doc="Solute permeability coeff. constant in all LSR stages",
        )
        m.fs.B_comp_system.set_value(
            m.fs.ROUnits[m.fs.LSRRO_Stages.first()].B_comp[0, "NaCl"]
        )
        m.fs.B_comp_system.setlb(3.5e-8)
        if B_max is None:
            m.fs.B_comp_system.setub(3.5e-8 * 1e2)
        else:
            m.fs.B_comp_system.setub(m.fs.B_max)
    if A_case == ACase.single_optimum:
        m.fs.A_comp_system = Var(
            domain=NonNegativeReals,
            units=pyunits.m * pyunits.s**-1 * pyunits.Pa**-1,
            doc="Water permeability coeff. constant in all LSR stages",
        )
        m.fs.A_comp_system.set_value(m.fs.ROUnits[2].A_comp[0, "H2O"])
        m.fs.A_comp_system.setlb(2.78e-12)
        m.fs.A_comp_system.setub(4.2e-11)
    if (
        AB_tradeoff == ABTradeoff.equality_constraint
        or AB_tradeoff == ABTradeoff.inequality_constraint
    ):
        m.fs.AB_tradeoff_coeff = Param(initialize=0.01333, mutable=True)
        m.fs.AB_tradeoff_coeff.set_value(
            AB_gamma_factor * value(m.fs.AB_tradeoff_coeff)
        )

    # unfix stages
    for idx, stage in m.fs.ROUnits.items():
        stage.area.unfix()
        stage.width.unfix()
        stage.area.setlb(1)
        stage.area.setub(20000)
        stage.width.setlb(0.1)
        stage.width.setub(1000)

        if (
            stage.config.mass_transfer_coefficient == MassTransferCoefficient.calculated
        ) or (stage.config.pressure_change_type == PressureChangeType.calculated):
            stage.N_Re[0, 0].unfix()

        if idx > m.fs.Stages.first():
            stage.B_comp.unfix()
            stage.B_comp.setlb(3.5e-8)
            if B_max is not None:
                stage.B_comp.setub(m.fs.B_max)
            else:
                stage.B_comp.setub(None)
            if B_case == BCase.single_optimum:
                stage.B_comp_equal = Constraint(
                    expr=stage.B_comp[0, "NaCl"] == m.fs.B_comp_system
                )
            if A_case == ACase.single_optimum:
                stage.A_comp_equal = Constraint(
                    expr=stage.A_comp[0, "H2O"] == m.fs.A_comp_system
                )

            stage.A_min = Param(
                initialize=2.78e-12, units=pyunits.m**2 / pyunits.kg * pyunits.s
            )
            stage.A_max = Param(
                initialize=4.2e-11, units=pyunits.m**2 / pyunits.kg * pyunits.s
            )
            stage._A_comp_con = Constraint(
                expr=(stage.A_min, stage.A_comp[0, "H2O"], stage.A_max)
            )
            iscale.constraint_scaling_transform(
                stage._A_comp_con, iscale.get_scaling_factor(stage.A_comp[0, "H2O"])
            )
            if A_case == ACase.optimize or A_case == ACase.single_optimum:
                stage.A_comp.unfix()
                stage.A_comp.setlb(2.78e-12)
                stage.A_comp.setub(4.2e-11)
            elif A_case == ACase.fixed:
                if not isinstance(A_value, (int, float)):
                    raise TypeError("A_value must be a numeric value")
                stage.A_comp.unfix()
                stage.A_comp.fix(A_value)
            else:
                raise TypeError(
                    f'A_case must be set to "fix", "single_optimum", "optimize" or None.'
                    f" A_case was set to {A_case}"
                )

            if AB_tradeoff == ABTradeoff.equality_constraint:
                stage.ABtradeoff = Constraint(
                    expr=pyunits.convert(
                        stage.B_comp[0, "NaCl"],
                        to_units=pyunits.L * pyunits.m**-2 * pyunits.hour**-1,
                    )
                    == m.fs.AB_tradeoff_coeff
                    * pyunits.convert(
                        stage.A_comp[0, "H2O"],
                        to_units=pyunits.L
                        * pyunits.m**-2
                        * pyunits.bar**-1
                        * pyunits.hour**-1,
                    )
                    ** 3
                    * pyunits.L**-2
                    * pyunits.m**4
                    * pyunits.hour**2
                    * pyunits.bar**3
                )
            elif AB_tradeoff == ABTradeoff.inequality_constraint:
                stage.ABtradeoff = Constraint(
                    expr=pyunits.convert(
                        stage.B_comp[0, "NaCl"],
                        to_units=pyunits.L * pyunits.m**-2 * pyunits.hour**-1,
                    )
                    >= m.fs.AB_tradeoff_coeff
                    * pyunits.convert(
                        stage.A_comp[0, "H2O"],
                        to_units=pyunits.L
                        * pyunits.m**-2
                        * pyunits.bar**-1
                        * pyunits.hour**-1,
                    )
                    ** 3
                    * pyunits.L**-2
                    * pyunits.m**4
                    * pyunits.hour**2
                    * pyunits.bar**3
                )
            else:
                pass

    min_avg_flux = 1  # minimum average water flux [kg/m2-h]
    # [kg/m2-s]
    min_avg_flux = min_avg_flux / 3600 * pyunits.kg / pyunits.m**2 / pyunits.s

    # additional constraints
    if water_recovery is not None:
        # product mass flow rate fraction of feed [-]
        m.fs.water_recovery.fix(water_recovery)
    if Cbrine is not None:
        # Final brine concentration
        m.fs.ROUnits[m.fs.Stages.last()].feed_side.properties[
            0, 1
        ].conc_mass_phase_comp["Liq", "NaCl"].fix(Cbrine)

    # add upper bound for permeate concentration
    if permeate_quality_limit is not None:
        if isinstance(permeate_quality_limit, (int, float)):
            m.fs.ROUnits[1].mixed_permeate[0].mass_frac_phase_comp["Liq", "NaCl"].setub(
                permeate_quality_limit
            )
        else:
            raise TypeError("permeate_quality_limit must be None, integer, or float")
    # ---checking model---
    assert_units_consistent(m)

    assert_degrees_of_freedom(
        m,
        4 * m.fs.NumberOfStages
        - (1 if (water_recovery is not None) else 0)
        - (1 if value(m.fs.NumberOfStages) == 1 else 0),
    )

    return m


def display_design(m):
    print("--decision variables--")
    for stage in m.fs.Stages:
        print(
            "Stage %d operating pressure %.1f bar"
            % (stage, m.fs.ROUnits[stage].inlet.pressure[0].value / 1e5)
        )
        print(
            "Stage %d membrane area      %.1f m2"
            % (stage, m.fs.ROUnits[stage].area.value)
        )
        print(
            "Stage %d water perm. coeff.  %.1f LMH/bar"
            % (stage, m.fs.ROUnits[stage].A_comp[0, "H2O"].value * (3.6e11))
        )
        print(
            "Stage %d salt perm. coeff.  %.1f LMH"
            % (stage, m.fs.ROUnits[stage].B_comp[0, "NaCl"].value * (1000.0 * 3600.0))
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
        if stage == m.fs.LastStage:
            pass
        else:
            print_state(f"Mixer {stage} recycle", m.fs.Mixers[stage].downstream)
            print_state(f"Mixer {stage} out", m.fs.Mixers[stage].outlet)

        print_state(f"RO {stage} permeate", m.fs.ROUnits[stage].permeate)
        print_state(f"RO {stage} retentate", m.fs.ROUnits[stage].retentate)
        wr = m.fs.ROUnits[stage].recovery_vol_phase[0, "Liq"].value

        if stage == m.fs.FirstStage:
            pass
        else:
            print_state(f"Booster Pump {stage} out", m.fs.BoosterPumps[stage].outlet)

    print_state(f"Disposal", m.fs.disposal.inlet)
    print_state(f"Product", m.fs.product.inlet)


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
    total_area = value(sum(m.fs.ROUnits[a].area for a in m.fs.Stages))
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
        f"Booster Pump Capital Cost ($/m3): "
        f"{value(m.fs.costing.factor_capital_annualization*sum(m.fs.BoosterPumps[stage].costing.capital_cost for stage in m.fs.LSRRO_Stages) / m.fs.costing.annual_water_production)}"
    )
    print(
        f"ERD Capital Cost ($/m3):"
        f"{value(m.fs.costing.factor_capital_annualization*sum(erd.costing.capital_cost for erd in m.fs.EnergyRecoveryDevices.values()) / m.fs.costing.annual_water_production)}"
    )
    print(
        f"Membrane Capital Cost ($/m3): "
        f"{value(m.fs.costing.factor_capital_annualization*sum(m.fs.ROUnits[stage].costing.capital_cost for stage in m.fs.Stages) / m.fs.costing.annual_water_production)}"
    )
    print(
        f"Indirect Capital Cost ($/m3): "
        f"{value(m.fs.costing.factor_capital_annualization*(m.fs.costing.total_investment_cost - m.fs.costing.total_capital_cost) / m.fs.costing.annual_water_production)}"
    )
    electricity_cost = value(
        m.fs.costing.aggregate_flow_costs["electricity"]
        * m.fs.costing.utilization_factor
        / m.fs.costing.annual_water_production
    )
    print(f"Electricity cost ($/m3): {electricity_cost}")

    print("\n")


def display_RO_reports(m):
    for stage in m.fs.ROUnits.values():
        stage.report()
