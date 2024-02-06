#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
import os
import idaes.logger as idaeslog
from pyomo.environ import (
    ConcreteModel,
    Block,
    Expression,
    Objective,
    value,
    TransformationFactory,
    units as pyunits,
)
from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock, MomentumBalanceType
from idaes.core.solvers import get_solver
from idaes.core.util.initialization import propagate_state

import idaes.core.util.scaling as iscale
from idaes.models.unit_models import (
    Mixer,
    Separator,
    Product,
    Translator,
    MomentumMixingType,
)
from idaes.models.unit_models.separator import (
    SplittingType,
    EnergySplittingType,
)

from idaes.core import UnitModelCostingBlock

from watertap.unit_models.pressure_exchanger import PressureExchanger
from watertap.unit_models.pressure_changer import Pump
from watertap.core.util.initialization import assert_degrees_of_freedom, check_solve

import watertap.property_models.seawater_prop_pack as prop_SW
from watertap.unit_models.reverse_osmosis_0D import (
    ReverseOsmosis0D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
from watertap.costing.unit_models.dewatering import (
    cost_centrifuge,
    cost_filter_belt_press,
    cost_filter_plate_press,
)
from watertap.costing import MultiUnitModelCostingBlock

from watertap.core.wt_database import Database
import watertap.core.zero_order_properties as prop_ZO
from watertap.unit_models.zero_order import (
    FeedZO,
    PumpElectricityZO,
    NanofiltrationZO,
    SecondaryTreatmentWWTPZO,
)
from watertap.costing.zero_order_costing import ZeroOrderCosting
from watertap.costing import WaterTAPCosting

# Set up logger
_log = idaeslog.getLogger(__name__)


def main():
    m = build(include_pretreatment=False, include_dewatering=True)
    set_operating_conditions(m)

    assert_degrees_of_freedom(m, 0)
    assert_units_consistent(m)

    initialize_system(m)
    assert_degrees_of_freedom(m, 0)

    results = solve(m, checkpoint="solve flowsheet after initializing system")

    add_costing(m)
    initialize_costing(m)
    assert_degrees_of_freedom(m, 0)  # ensures problem is square

    optimize_operation(m)  # unfixes specific variables for cost optimization

    solve(m, checkpoint="solve flowsheet after costing")
    display_results(m)
    display_costing(m)

    return m, results


def build(include_pretreatment=False, include_dewatering=False):
    # flowsheet set up
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)

    # define property packages
    m.fs.prop_nf = prop_ZO.WaterParameterBlock(solute_list=["dye", "tds"])
    m.fs.prop_ro = prop_SW.SeawaterParameterBlock()

    # define blocks
    if include_pretreatment == True:
        prtrt = m.fs.pretreatment = Block()
        m.fs.wwt_retentate = Product(property_package=m.fs.prop_nf)
    else:
        pass

    dye_sep = m.fs.dye_separation = Block()
    desal = m.fs.desalination = Block()

    # define flowsheet inlets and outlets
    m.fs.feed = FeedZO(property_package=m.fs.prop_nf)

    if include_dewatering == True:
        m.fs.dewaterer = Separator(
            property_package=m.fs.prop_nf,
            outlet_list=["centrate", "precipitant"],
            split_basis=SplittingType.componentFlow,
            energy_split_basis=EnergySplittingType.none,
            momentum_balance_type=MomentumBalanceType.none,
        )
        m.fs.centrate = Product(property_package=m.fs.prop_nf)
        m.fs.precipitant = Product(property_package=m.fs.prop_nf)
    else:
        m.fs.dye_retentate = Product(property_package=m.fs.prop_nf)

    m.fs.permeate = Product(property_package=m.fs.prop_ro)
    m.fs.brine = Product(property_package=m.fs.prop_ro)

    # pretreatment
    if hasattr(m.fs, "pretreatment"):
        prtrt.wwtp = SecondaryTreatmentWWTPZO(
            property_package=m.fs.prop_nf, database=m.db, process_subtype="default"
        )
    else:
        pass

    # nanofiltration components
    dye_sep.P1 = PumpElectricityZO(
        property_package=m.fs.prop_nf, database=m.db, process_subtype="default"
    )

    dye_sep.nanofiltration = NanofiltrationZO(
        property_package=m.fs.prop_nf,
        database=m.db,
        process_subtype="rHGO_dye_rejection",
    )

    # reverse osmosis components

    desal.P2 = Pump(property_package=m.fs.prop_ro)
    desal.RO = ReverseOsmosis0D(
        property_package=m.fs.prop_ro,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
    )

    desal.RO.width.setub(2000)
    desal.RO.area.setub(20000)

    desal.S1 = Separator(property_package=m.fs.prop_ro, outlet_list=["P2", "PXR"])
    desal.M1 = Mixer(
        property_package=m.fs.prop_ro,
        momentum_mixing_type=MomentumMixingType.equality,
        inlet_list=["P2", "P3"],
    )
    desal.PXR = PressureExchanger(property_package=m.fs.prop_ro)
    desal.P3 = Pump(property_package=m.fs.prop_ro)

    # translator blocks
    m.fs.tb_nf_ro = Translator(
        inlet_property_package=m.fs.prop_nf, outlet_property_package=m.fs.prop_ro
    )

    # since the dye << tds: Assume RO_TDS = NF_tds + NF_dye
    @m.fs.tb_nf_ro.Constraint(["H2O", "dye"])
    def eq_flow_mass_comp(blk, j):
        if j == "dye":
            return (
                blk.properties_in[0].flow_mass_comp["dye"]
                + blk.properties_in[0].flow_mass_comp["tds"]
                == blk.properties_out[0].flow_mass_phase_comp["Liq", "TDS"]
            )
        else:
            return (
                blk.properties_in[0].flow_mass_comp["H2O"]
                == blk.properties_out[0].flow_mass_phase_comp["Liq", "H2O"]
            )

    # connections
    if hasattr(m.fs, "pretreatment"):
        m.fs.s_feed = Arc(source=m.fs.feed.outlet, destination=prtrt.wwtp.inlet)
        prtrt.s01 = Arc(source=prtrt.wwtp.treated, destination=dye_sep.P1.inlet)
        prtrt.s02 = Arc(
            source=prtrt.wwtp.byproduct, destination=m.fs.wwt_retentate.inlet
        )
    else:
        m.fs.s_feed = Arc(source=m.fs.feed.outlet, destination=dye_sep.P1.inlet)

    dye_sep.s01 = Arc(
        source=dye_sep.P1.outlet, destination=dye_sep.nanofiltration.inlet
    )
    if hasattr(m.fs, "dewaterer"):
        m.fs.s01 = Arc(
            source=dye_sep.nanofiltration.byproduct, destination=m.fs.dewaterer.inlet
        )
        # TODO: Recycle centrate stream back to the feed via a mixer
        m.fs.s02 = Arc(source=m.fs.dewaterer.centrate, destination=m.fs.centrate.inlet)
        m.fs.s03 = Arc(
            source=m.fs.dewaterer.precipitant, destination=m.fs.precipitant.inlet
        )
    else:
        dye_sep.s02 = Arc(
            source=dye_sep.nanofiltration.byproduct,
            destination=m.fs.dye_retentate.inlet,
        )
    m.fs.s_nf = Arc(
        source=dye_sep.nanofiltration.treated, destination=m.fs.tb_nf_ro.inlet
    )

    m.fs.s_ro = Arc(source=m.fs.tb_nf_ro.outlet, destination=desal.S1.inlet)
    desal.s01 = Arc(source=desal.S1.P2, destination=desal.P2.inlet)
    desal.s02 = Arc(source=desal.P2.outlet, destination=desal.M1.P2)
    desal.s03 = Arc(source=desal.M1.outlet, destination=desal.RO.inlet)
    desal.s04 = Arc(
        source=desal.RO.retentate, destination=desal.PXR.high_pressure_inlet
    )
    desal.s05 = Arc(source=desal.S1.PXR, destination=desal.PXR.low_pressure_inlet)
    desal.s06 = Arc(source=desal.PXR.low_pressure_outlet, destination=desal.P3.inlet)
    desal.s07 = Arc(source=desal.P3.outlet, destination=desal.M1.P3)
    m.fs.s_disposal = Arc(
        source=desal.PXR.high_pressure_outlet, destination=m.fs.brine.inlet
    )

    m.fs.s_permeate = Arc(source=desal.RO.permeate, destination=m.fs.permeate.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling
    m.fs.prop_ro.set_default_scaling("flow_mass_phase_comp", 1e-3, index=("Liq", "H2O"))
    m.fs.prop_ro.set_default_scaling("flow_mass_phase_comp", 1e-1, index=("Liq", "TDS"))

    # set unit model values
    iscale.set_scaling_factor(desal.P2.control_volume.work, 1e-6)
    iscale.set_scaling_factor(desal.P2.control_volume.deltaP, 1e-7)
    iscale.set_scaling_factor(
        desal.P2.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "H2O"],
        1e-1,
    )
    iscale.set_scaling_factor(desal.RO.area, 1e-4)
    iscale.set_scaling_factor(desal.P3.control_volume.work, 1e-5)
    iscale.set_scaling_factor(desal.PXR.low_pressure_side.work, 1e-5)
    iscale.set_scaling_factor(desal.PXR.high_pressure_side.work, 1e-5)

    # calculate and propagate scaling factors
    iscale.calculate_scaling_factors(m)
    return m


def set_operating_conditions(m):
    if hasattr(m.fs, "pretreatment"):
        prtrt = m.fs.pretreatment
    else:
        pass

    dye_sep = m.fs.dye_separation
    desal = m.fs.desalination

    # feed
    flow_vol = 280 / 3600 * pyunits.m**3 / pyunits.s
    conc_mass_dye = 0.2 * pyunits.kg / pyunits.m**3
    conc_mass_tds = 2 * pyunits.kg / pyunits.m**3
    temperature = 298 * pyunits.K
    pressure = 101325 * pyunits.Pa

    m.fs.feed.flow_vol[0].fix(flow_vol)
    m.fs.feed.conc_mass_comp[0, "dye"].fix(conc_mass_dye)
    m.fs.feed.conc_mass_comp[0, "tds"].fix(conc_mass_tds)
    solve(m.fs.feed, checkpoint="solve feed block")

    # pretreatment
    if hasattr(m.fs, "pretreatment"):
        prtrt.wwtp.load_parameters_from_database(use_default_removal=True)
    else:
        pass

    # nanofiltration
    dye_sep.nanofiltration.load_parameters_from_database(use_default_removal=True)

    if hasattr(m.fs, "dewaterer"):
        m.fs.dewaterer.split_fraction[0, "precipitant", "H2O"].fix(0.01)
        m.fs.dewaterer.split_fraction[0, "precipitant", "tds"].fix(0.01)
        m.fs.dewaterer.split_fraction[0, "precipitant", "dye"].fix(0.99)
    else:
        pass

    # nf pump
    dye_sep.P1.load_parameters_from_database(use_default_removal=True)
    dye_sep.P1.applied_pressure.fix(
        dye_sep.nanofiltration.applied_pressure.get_values()[0]
    )
    dye_sep.P1.eta_pump.fix(0.75)  # pump efficiency [-]
    dye_sep.P1.lift_height.unfix()

    # desalination
    desal.P2.efficiency_pump.fix(0.80)
    operating_pressure = 70e5 * pyunits.Pa
    desal.P2.control_volume.properties_out[0].pressure.fix(operating_pressure)
    desal.RO.A_comp.fix(4.2e-12)  # membrane water permeability
    desal.RO.B_comp.fix(3.5e-8)  # membrane salt permeability
    desal.RO.feed_side.channel_height.fix(1e-3)  # channel height in membrane stage [m]
    desal.RO.feed_side.spacer_porosity.fix(
        0.97
    )  # spacer porosity in membrane stage [-]
    desal.RO.permeate.pressure[0].fix(pressure)  # atmospheric pressure [Pa]
    desal.RO.feed_side.velocity[0, 0].fix(0.25)
    desal.RO.recovery_vol_phase[0, "Liq"].fix(0.5)
    m.fs.tb_nf_ro.properties_out[0].temperature.fix(temperature)
    m.fs.tb_nf_ro.properties_out[0].pressure.fix(pressure)

    # pressure exchanger
    desal.PXR.efficiency_pressure_exchanger.fix(0.95)
    # booster pump
    desal.P3.efficiency_pump.fix(0.80)

    return


def initialize_system(m):
    if hasattr(m.fs, "pretreatment"):
        prtrt = m.fs.pretreatment
    else:
        pass

    dye_sep = m.fs.dye_separation
    desal = m.fs.desalination

    # initialize feed
    solve(m.fs.feed, checkpoint="solve flowsheet after initializing feed")

    # initialize pretreatment
    propagate_state(m.fs.s_feed)
    seq = SequentialDecomposition()
    seq.options.tear_set = []
    seq.options.iterLim = 1

    if hasattr(m.fs, "pretreatment"):
        seq.run(prtrt, lambda u: u.initialize())
        propagate_state(prtrt.s01)
    else:
        pass

    # initialize nf
    seq.run(dye_sep, lambda u: u.initialize())

    if hasattr(m.fs, "dewater"):
        seq.run(m.fs.dewaterer, lambda u: u.initialize())
        propagate_state(m.fs.s01)
    else:
        pass

    # initialize ro
    propagate_state(m.fs.s_nf)
    propagate_state(m.fs.s_ro)
    propagate_state(desal.s01)
    propagate_state(m.fs.s_disposal)
    propagate_state(m.fs.s_permeate)

    m.fs.tb_nf_ro.properties_out[0].flow_mass_phase_comp["Liq", "H2O"] = value(
        m.fs.tb_nf_ro.properties_in[0].flow_mass_comp["H2O"]
    )
    m.fs.tb_nf_ro.properties_out[0].flow_mass_phase_comp["Liq", "TDS"] = value(
        m.fs.tb_nf_ro.properties_in[0].flow_mass_comp["tds"]
    ) + value(m.fs.tb_nf_ro.properties_in[0].flow_mass_comp["dye"])

    desal.RO.feed_side.properties_in[0].flow_mass_phase_comp["Liq", "H2O"] = value(
        m.fs.tb_nf_ro.properties_out[0].flow_mass_phase_comp["Liq", "H2O"]
    )
    desal.RO.feed_side.properties_in[0].temperature = value(
        m.fs.tb_nf_ro.properties_out[0].temperature
    )
    desal.RO.feed_side.properties_in[0].pressure = value(
        desal.P2.control_volume.properties_out[0].pressure
    )
    desal.RO.initialize()
    solve(desal, checkpoint="solve flowsheet after initializing desalination")
    return


def optimize_operation(m):
    """
    Unfixes RO operating conditions and sets solver objective
        - Operating pressure: 1 - 83 bar
        - Crossflow velocity: 10 - 30 cm/s
        - Membrane area: 50 - 5000 m2
        - Volumetric recovery: 10 - 75 %
    """
    desal = m.fs.desalination

    # RO operating pressure
    desal.P2.control_volume.properties_out[0].pressure.unfix()
    desal.P2.control_volume.properties_out[0].pressure.setub(
        8300000
    )  # pressure vessel burst pressure
    desal.P2.control_volume.properties_out[0].pressure.setlb(100000)

    # RO inlet velocity
    desal.RO.feed_side.velocity[0, 0].unfix()
    desal.RO.feed_side.velocity[0, 0].setub(0.3)
    desal.RO.feed_side.velocity[0, 0].setlb(0.1)

    # RO membrane area
    desal.RO.area.unfix()
    desal.RO.area.setub(5000)
    desal.RO.area.setlb(50)

    # RO recovery - likely limited by operating pressure
    desal.RO.recovery_vol_phase[0, "Liq"].unfix()
    desal.RO.recovery_vol_phase[0, "Liq"].setub(0.99)
    desal.RO.recovery_vol_phase[0, "Liq"].setlb(0.1)

    # Permeate salt concentration constraint
    m.fs.permeate.properties[0].conc_mass_phase_comp["Liq", "TDS"].setub(0.5)
    m.fs.brine.properties[0].conc_mass_phase_comp
    m.fs.objective = Objective(expr=m.fs.LCOT)
    return


def solve(blk, solver=None, checkpoint=None, tee=False, fail_flag=True):
    if solver is None:
        solver = get_solver()
    results = solver.solve(blk, tee=tee)
    check_solve(results, checkpoint=checkpoint, logger=_log, fail_flag=fail_flag)
    return results


def add_costing(m):
    if hasattr(m.fs, "pretreatment"):
        prtrt = m.fs.pretreatment
    else:
        pass

    dye_sep = m.fs.dye_separation
    desal = m.fs.desalination

    # Zero order costing
    source_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "dye_desalination_global_costing.yaml",
    )

    m.fs.zo_costing = ZeroOrderCosting(case_study_definition=source_file)
    m.fs.ro_costing = WaterTAPCosting()

    # cost nanofiltration module and pump
    if hasattr(m.fs, "pretreatment"):
        prtrt.wwtp.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.zo_costing
        )
    else:
        pass

    dye_sep.nanofiltration.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.zo_costing
    )

    if hasattr(m.fs, "dewaterer"):
        m.fs.dewaterer.costing = MultiUnitModelCostingBlock(
            flowsheet_costing_block=m.fs.ro_costing,
            costing_blocks={
                "centrifuge": {
                    "costing_method": cost_centrifuge,
                    "costing_method_arguments": {"cost_electricity_flow": False},
                },
                "filter_belt_press": {
                    "costing_method": cost_filter_belt_press,
                    "costing_method_arguments": {"cost_electricity_flow": False},
                },
                "filter_plate_press": {
                    "costing_method": cost_filter_plate_press,
                    "costing_method_arguments": {"cost_electricity_flow": False},
                },
            },
            initial_costing_block="centrifuge",
        )
    else:
        pass

    dye_sep.P1.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.zo_costing)
    m.fs.zo_costing.pump_electricity.pump_cost["default"].fix(76)

    # RO Train
    # RO equipment is costed using more detailed costing package

    desal.P2.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.ro_costing,
        costing_method_arguments={"cost_electricity_flow": True},
    )
    desal.RO.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.ro_costing)

    desal.M1.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.ro_costing)
    desal.PXR.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.ro_costing)
    desal.P3.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.ro_costing,
        costing_method_arguments={"cost_electricity_flow": True},
    )
    m.fs.ro_costing.electricity_cost = value(m.fs.zo_costing.electricity_cost)
    m.fs.ro_costing.base_currency = pyunits.USD_2020
    # Aggregate unit level costs and calculate overall process costs
    m.fs.zo_costing.cost_process()
    m.fs.ro_costing.cost_process()

    # Add specific energy consumption
    feed_flowrate = m.fs.feed.flow_vol[0]
    m.fs.zo_costing.add_electricity_intensity(feed_flowrate)
    m.fs.ro_costing.add_specific_energy_consumption(feed_flowrate)

    m.fs.specific_energy_intensity = Expression(
        expr=(
            m.fs.zo_costing.electricity_intensity
            + m.fs.ro_costing.specific_energy_consumption
        ),
        doc="Specific energy consumption of the treatment train on a feed flowrate basis [kWh/m3]",
    )

    # Annual disposal of brine
    m.fs.brine_disposal_cost = Expression(
        expr=(
            m.fs.zo_costing.utilization_factor
            * (
                m.fs.zo_costing.waste_disposal_cost
                * pyunits.convert(
                    m.fs.brine.properties[0].flow_vol,
                    to_units=pyunits.m**3 / m.fs.zo_costing.base_period,
                )
            )
        ),
        doc="Cost of disposing of brine waste",
    )

    if hasattr(m.fs, "pretreatment"):
        m.fs.sludge_disposal_cost = Expression(
            expr=(
                m.fs.zo_costing.utilization_factor
                * (
                    m.fs.zo_costing.waste_disposal_cost
                    * pyunits.convert(
                        m.fs.wwt_retentate.properties[0].flow_vol,
                        to_units=pyunits.m**3 / m.fs.zo_costing.base_period,
                    )
                )
            ),
            doc="Cost of disposing of waste water treatment plant sludge",
        )
    else:
        pass
    if hasattr(m.fs, "dewaterer"):
        m.fs.dye_disposal_cost = Expression(
            expr=(
                m.fs.zo_costing.utilization_factor
                * m.fs.zo_costing.dewatered_dye_disposal_cost
                * pyunits.convert(
                    m.fs.precipitant.properties[0].flow_vol,
                    to_units=pyunits.m**3 / m.fs.zo_costing.base_period,
                )
            ),
            doc="Cost of disposing of dye waste",
        )
    else:
        m.fs.dye_disposal_cost = Expression(
            expr=(
                m.fs.zo_costing.utilization_factor
                * m.fs.zo_costing.dye_disposal_cost
                * pyunits.convert(
                    m.fs.dye_retentate.properties[0].flow_vol,
                    to_units=pyunits.m**3 / m.fs.zo_costing.base_period,
                )
            ),
            doc="Cost of disposing of dye waste",
        )

    if hasattr(m.fs, "dewaterer"):
        # TODO: Remove centrate stream from this calculation after implementing recycle
        m.fs.water_recovery_revenue = Expression(
            expr=(
                m.fs.zo_costing.utilization_factor
                * m.fs.zo_costing.recovered_water_cost
                * pyunits.convert(
                    m.fs.permeate.properties[0].flow_vol
                    + m.fs.centrate.properties[0].flow_vol,
                    to_units=pyunits.m**3 / m.fs.zo_costing.base_period,
                )
            ),
            doc="Savings from water recovered back to the plant",
        )
    else:
        m.fs.water_recovery_revenue = Expression(
            expr=(
                m.fs.zo_costing.utilization_factor
                * m.fs.zo_costing.recovered_water_cost
                * pyunits.convert(
                    m.fs.permeate.properties[0].flow_vol,
                    to_units=pyunits.m**3 / m.fs.zo_costing.base_period,
                )
            ),
            doc="Savings from water recovered back to the plant",
        )

    # Combine results from costing packages and calculate overall metrics
    @m.fs.Expression(doc="Total capital cost of the treatment train")
    def total_capital_cost(b):
        return pyunits.convert(
            m.fs.zo_costing.total_capital_cost, to_units=pyunits.USD_2020
        ) + pyunits.convert(
            m.fs.ro_costing.total_capital_cost, to_units=pyunits.USD_2020
        )

    @m.fs.Expression(doc="Total operating cost of the treatment train")
    def total_operating_cost(b):
        return (
            pyunits.convert(
                m.fs.zo_costing.total_fixed_operating_cost,
                to_units=pyunits.USD_2020 / pyunits.year,
            )
            + pyunits.convert(
                m.fs.zo_costing.total_variable_operating_cost,
                to_units=pyunits.USD_2020 / pyunits.year,
            )
            + pyunits.convert(
                m.fs.ro_costing.total_operating_cost,
                to_units=pyunits.USD_2020 / pyunits.year,
            )
        )

    @m.fs.Expression(doc="Total cost of water recovery and brine/sludge/dye disposed")
    def total_externalities(b):
        if hasattr(m.fs, "pretreatment"):
            return pyunits.convert(
                m.fs.water_recovery_revenue
                - m.fs.dye_disposal_cost
                - m.fs.brine_disposal_cost
                - m.fs.sludge_disposal_cost,
                to_units=pyunits.USD_2020 / pyunits.year,
            )
        else:
            return pyunits.convert(
                m.fs.water_recovery_revenue
                - m.fs.dye_disposal_cost
                - m.fs.brine_disposal_cost,
                to_units=pyunits.USD_2020 / pyunits.year,
            )

    @m.fs.Expression(
        doc="Levelized cost of treatment with respect to volumetric feed flow"
    )
    def LCOT(b):
        return (
            b.total_capital_cost * b.zo_costing.capital_recovery_factor
            + b.total_operating_cost
            - b.total_externalities
        ) / (
            pyunits.convert(
                b.feed.properties[0].flow_vol,
                to_units=pyunits.m**3 / pyunits.year,
            )
            * b.zo_costing.utilization_factor
        )

    @m.fs.Expression(
        doc="Levelized cost of treatment with respect to volumetric feed flow, not including externalities"
    )
    def LCOT_wo_revenue(b):
        return (
            b.total_capital_cost * b.zo_costing.capital_recovery_factor
            + b.total_operating_cost
        ) / (
            pyunits.convert(
                b.feed.properties[0].flow_vol,
                to_units=pyunits.m**3 / pyunits.year,
            )
            * b.zo_costing.utilization_factor
        )

    @m.fs.Expression(
        doc="Levelized cost of water with respect to volumetric permeate flow"
    )
    def LCOW(b):
        if hasattr(m.fs, "pretreatment"):
            return (
                b.total_capital_cost * b.zo_costing.capital_recovery_factor
                + b.total_operating_cost
                - pyunits.convert(
                    -m.fs.dye_disposal_cost
                    - m.fs.brine_disposal_cost
                    - m.fs.sludge_disposal_cost,
                    to_units=pyunits.USD_2020 / pyunits.year,
                )
            ) / (
                pyunits.convert(
                    b.permeate.properties[0].flow_vol,
                    to_units=pyunits.m**3 / pyunits.year,
                )
                * b.zo_costing.utilization_factor
            )
        else:
            return (
                b.total_capital_cost * b.zo_costing.capital_recovery_factor
                + b.total_operating_cost
                - pyunits.convert(
                    -m.fs.dye_disposal_cost - m.fs.brine_disposal_cost,
                    to_units=pyunits.USD_2020 / pyunits.year,
                )
            ) / (
                pyunits.convert(
                    b.permeate.properties[0].flow_vol,
                    to_units=pyunits.m**3 / pyunits.year,
                )
                * b.zo_costing.utilization_factor
            )

    @m.fs.Expression(
        doc="Levelized cost of water with respect to volumetric permeate flow"
    )
    def LCOW_wo_revenue(b):
        return (
            b.total_capital_cost * b.zo_costing.capital_recovery_factor
            + b.total_operating_cost
        ) / (
            pyunits.convert(
                b.permeate.properties[0].flow_vol,
                to_units=pyunits.m**3 / pyunits.year,
            )
            * b.zo_costing.utilization_factor
        )

    assert_units_consistent(m)
    return


def initialize_costing(m):
    m.fs.zo_costing.initialize()
    m.fs.ro_costing.initialize()
    return


def display_results(m):
    print("\nUnit models:")
    if hasattr(m.fs, "pretreatment"):
        m.fs.pretreatment.wwtp.report()
    else:
        pass

    if hasattr(m.fs, "dewaterer"):
        m.fs.dewaterer.report()
    else:
        pass

    m.fs.dye_separation.P1.report()
    m.fs.dye_separation.nanofiltration.report()
    m.fs.desalination.RO.report()

    print("\nStreams:")
    if hasattr(m.fs, "pretreatment") and hasattr(m.fs, "dewater"):
        flow_list = ["feed", "wwt_retentate", "precipitant", "centrate"]
    elif hasattr(m.fs, "pretreatment"):
        flow_list = ["feed", "wwt_retentate", "dye_retentate"]
    elif hasattr(m.fs, "dewaterer"):
        flow_list = ["feed", "precipitant", "centrate"]
    else:
        flow_list = ["feed", "dye_retentate"]

    for f in flow_list:
        m.fs.component(f).report()

    if hasattr(m.fs, "dewaterer"):
        precipitant_vol_flowrate = value(
            pyunits.convert(
                m.fs.precipitant.properties[0].flow_vol,
                to_units=pyunits.m**3 / pyunits.hr,
            )
        )
        precipitant_tds_concentration = m.fs.precipitant.flow_mass_comp[0, "tds"].value
        precipitant_dye_concentration = m.fs.precipitant.flow_mass_comp[0, "dye"].value

        centrate_vol_flowrate = value(
            pyunits.convert(
                m.fs.centrate.properties[0].flow_vol,
                to_units=pyunits.m**3 / pyunits.hr,
            )
        )
        centrate_tds_concentration = m.fs.centrate.flow_mass_comp[0, "tds"].value
        centrate_dye_concentration = m.fs.centrate.flow_mass_comp[0, "dye"].value

        print(
            f"\nPrecipitant volumetric flowrate: {precipitant_vol_flowrate : .3f} m3/hr"
        )
        print(
            f"Precipitant tds concentration: {precipitant_tds_concentration : .3f} g/l"
        )
        print(
            f"Precipitant dye concentration: {precipitant_dye_concentration : .3f} g/l"
        )

        print(f"\nCentrate volumetric flowrate: {centrate_vol_flowrate : .3f} m3/hr")
        print(f"Centrate tds concentration: {centrate_tds_concentration : .3f} g/l")
        print(f"Centrate dye concentration: {centrate_dye_concentration : .3f} g/l")
    else:
        dye_retentate_vol_flowrate = value(
            pyunits.convert(
                m.fs.dye_retentate.properties[0].flow_vol,
                to_units=pyunits.m**3 / pyunits.hr,
            )
        )
        dye_retentate_tds_concentration = m.fs.dye_retentate.flow_mass_comp[
            0, "tds"
        ].value
        dye_retentate_dye_concentration = m.fs.dye_retentate.flow_mass_comp[
            0, "dye"
        ].value

        print(
            f"\nDye retentate volumetric flowrate: {dye_retentate_vol_flowrate : .3f} m3/hr"
        )
        print(
            f"Dye retentate tds concentration: {dye_retentate_tds_concentration : .3f} g/l"
        )
        print(
            f"Dye retentate dye concentration: {dye_retentate_dye_concentration : .3f} g/l"
        )

    if hasattr(m.fs, "pretreatment"):
        wwt_retentate_vol_flowrate = value(
            pyunits.convert(
                m.fs.wwt_retentate.properties[0].flow_vol,
                to_units=pyunits.m**3 / pyunits.hr,
            )
        )
    else:
        pass

    permeate_salt_concentration = (
        m.fs.permeate.properties[0].conc_mass_phase_comp["Liq", "TDS"].value
    )
    permeate_vol_flowrate = value(
        pyunits.convert(
            m.fs.permeate.properties[0].flow_vol, to_units=pyunits.m**3 / pyunits.hr
        )
    )
    brine_salt_concentration = (
        m.fs.brine.properties[0].conc_mass_phase_comp["Liq", "TDS"].value
    )
    brine_vol_flowrate = value(
        pyunits.convert(
            m.fs.brine.properties[0].flow_vol, to_units=pyunits.m**3 / pyunits.hr
        )
    )

    print(f"\nPermeate volumetric flowrate: {permeate_vol_flowrate : .3f} m3/hr")
    print(f"Permeate salt concentration: {permeate_salt_concentration : .3f} g/l")
    print(f"\nBrine volumetric flowrate: {brine_vol_flowrate : .3f} m3/hr")
    print(f"Brine salt concentration: {brine_salt_concentration : .3f} g/l")

    if hasattr(m.fs, "pretreatment"):
        print(
            f"\nWastewater volumetric flowrate: {wwt_retentate_vol_flowrate : .3f} m3/hr"
        )
    else:
        pass

    print("\nSystem Recovery:")
    if hasattr(m.fs, "dewaterer"):
        sys_dye_recovery = (
            m.fs.precipitant.flow_mass_comp[0, "dye"]()
            / m.fs.feed.flow_mass_comp[0, "dye"]()
        )
    else:
        sys_dye_recovery = (
            m.fs.dye_retentate.flow_mass_comp[0, "dye"]()
            / m.fs.feed.flow_mass_comp[0, "dye"]()
        )
    sys_water_recovery = (
        m.fs.permeate.flow_mass_phase_comp[0, "Liq", "H2O"]()
        / m.fs.feed.flow_mass_comp[0, "H2O"]()
    )

    print(f"System dye recovery: {sys_dye_recovery*100 : .3f} %")
    print(f"System water recovery: {sys_water_recovery*100 : .3f} %")

    return


def display_costing(m):
    capex = value(pyunits.convert(m.fs.total_capital_cost, to_units=pyunits.MUSD_2020))
    if hasattr(m.fs, "pretreatment"):
        wwtp_capex = value(
            pyunits.convert(
                m.fs.pretreatment.wwtp.costing.capital_cost, to_units=pyunits.USD_2020
            )
        )

        wwtp_opex = value(
            m.fs.pretreatment.wwtp.energy_electric_flow_vol_inlet
            * m.fs.zo_costing.electricity_cost
            * m.fs.zo_costing.utilization_factor
            * pyunits.convert(
                m.fs.feed.flow_vol[0], to_units=pyunits.m**3 / pyunits.year
            )
        )
    else:
        pass

    if hasattr(m.fs, "dewaterer"):
        dewater_capex = value(
            pyunits.convert(
                m.fs.dewaterer.costing.capital_cost, to_units=pyunits.USD_2020
            )
        )
    else:
        pass

    nf_capex = value(
        pyunits.convert(
            m.fs.dye_separation.nanofiltration.costing.capital_cost
            + m.fs.dye_separation.P1.costing.capital_cost,
            to_units=pyunits.USD_2020,
        )
    )

    ro_capex = value(
        pyunits.convert(m.fs.ro_costing.total_capital_cost, to_units=pyunits.USD_2020)
    )

    opex = value(
        pyunits.convert(
            m.fs.total_operating_cost, to_units=pyunits.MUSD_2020 / pyunits.year
        )
    )

    if hasattr(m.fs, "pretreatment"):
        nf_opex = (
            value(
                pyunits.convert(
                    m.fs.zo_costing.total_operating_cost,
                    to_units=pyunits.USD_2020 / pyunits.year,
                )
            )
            - wwtp_opex
        )
    else:
        nf_opex = value(
            pyunits.convert(
                m.fs.zo_costing.total_operating_cost,
                to_units=pyunits.USD_2020 / pyunits.year,
            )
        )

    ro_opex = value(
        pyunits.convert(
            m.fs.ro_costing.total_operating_cost,
            to_units=pyunits.USD_2020 / pyunits.year,
        )
    )

    externalities = value(
        pyunits.convert(
            m.fs.total_externalities, to_units=pyunits.MUSD_2020 / pyunits.year
        )
    )
    wrr = value(
        pyunits.convert(
            m.fs.water_recovery_revenue, to_units=pyunits.USD_2020 / pyunits.year
        )
    )
    ddc = value(
        pyunits.convert(
            m.fs.dye_disposal_cost, to_units=pyunits.USD_2020 / pyunits.year
        )
    )
    bdc = value(
        pyunits.convert(
            m.fs.brine_disposal_cost, to_units=pyunits.USD_2020 / pyunits.year
        )
    )
    if hasattr(m.fs, "pretreatment"):
        sdc = value(
            pyunits.convert(
                m.fs.sludge_disposal_cost, to_units=pyunits.USD_2020 / pyunits.year
            )
        )
    else:
        pass

    # normalized costs
    feed_flowrate = value(
        pyunits.convert(
            m.fs.feed.properties[0].flow_vol, to_units=pyunits.m**3 / pyunits.hr
        )
    )
    capex_norm = (
        value(pyunits.convert(m.fs.total_capital_cost, to_units=pyunits.USD_2020))
        / feed_flowrate
    )

    annual_investment = value(
        pyunits.convert(
            m.fs.total_capital_cost * m.fs.zo_costing.capital_recovery_factor
            + m.fs.total_operating_cost,
            to_units=pyunits.USD_2020 / pyunits.year,
        )
    )
    opex_fraction = (
        100
        * value(
            pyunits.convert(
                m.fs.total_operating_cost, to_units=pyunits.USD_2020 / pyunits.year
            )
        )
        / annual_investment
    )

    lcot = value(pyunits.convert(m.fs.LCOT, to_units=pyunits.USD_2020 / pyunits.m**3))
    lcot_wo_rev = value(
        pyunits.convert(
            m.fs.LCOT_wo_revenue, to_units=pyunits.USD_2020 / pyunits.m**3
        )
    )

    lcow = value(pyunits.convert(m.fs.LCOW, to_units=pyunits.USD_2020 / pyunits.m**3))
    lcow_wo_rev = value(
        pyunits.convert(
            m.fs.LCOW_wo_revenue, to_units=pyunits.USD_2020 / pyunits.m**3
        )
    )

    sec = m.fs.specific_energy_intensity()

    print("\n System costing metrics:")
    print(f"\nTotal Capital Cost: {capex:.4f} M$")
    if hasattr(m.fs, "pretreatment"):
        print(f"Wastewater Treatment Capital Cost: {wwtp_capex:.4f} $")
    else:
        pass
    if hasattr(m.fs, "dewaterer"):
        print(f"Dewatering Unit Capital Cost: {dewater_capex:.4f} $")
    else:
        pass
    print(f"Nanofiltration (r-HGO) Capital Cost: {nf_capex:.4f} $")
    print(f"Reverse Osmosis Capital Cost: {ro_capex:.4f} $")

    print("\n----------Unit Capital Costs----------\n")
    for u in m.fs.zo_costing._registered_unit_costing:
        print(
            u.name,
            " : {price:0.3f} $".format(
                price=value(pyunits.convert(u.capital_cost, to_units=pyunits.USD_2020))
            ),
        )
    for z in m.fs.ro_costing._registered_unit_costing:
        print(
            z.name,
            " : {price:0.3f} $".format(
                price=value(pyunits.convert(z.capital_cost, to_units=pyunits.USD_2020))
            ),
        )

    print(f"\nTotal Operating Cost: {opex:.4f} M$/year")
    if hasattr(m.fs, "pretreatment"):
        print(f"Wastewater Treatment Operating Cost: {wwtp_opex:.4f} $/yr")
    else:
        pass
    print(f"Nanofiltration (r-HGO) Operating Cost: {nf_opex:.4f} $/yr")
    print(f"Reverse Osmosis Operating Cost: {ro_opex:.4f} $/yr")

    print(f"\nTotal Externalities: {externalities:.4f} M$/year")
    print(f"Water recovery revenue: {wrr: .4f} USD/year")
    print(f"Dye disposal cost: {-1*ddc: .4f} USD/year")
    print(f"Brine disposal cost: {-1*bdc: .4f} USD/year")
    if hasattr(m.fs, "pretreatment"):
        print(f"Sludge disposal cost: {-1*sdc: .4f} USD/year")
    else:
        pass

    print(f"\nTotal Annual Cost: {annual_investment : .4f} $/year")
    print(f"Normalized Capital Cost: {capex_norm:.4f} $/m3feed/hr")
    print(f"Opex Fraction of Annual Cost:{opex_fraction : .4f} %")

    print(f"Levelized cost of treatment with revenue: {lcot:.4f} $/m3 feed")
    print(f"Levelized cost of water with revenue: {lcow:.4f} $/m3 permeate")
    print(f"Levelized cost of treatment without revenue: {lcot_wo_rev:.4f} $/m3 feed")
    print(f"Levelized cost of water without revenue: {lcow_wo_rev:.4f} $/m3 permeate")

    print(f"Specific energy intensity: {sec:.3f} kWh/m3 feed")


if __name__ == "__main__":
    model, results = main()
