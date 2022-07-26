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
    Constraint,
    ConcreteModel,
    Block,
    Expression,
    Objective,
    value,
    TransformationFactory,
    units as pyunits,
    assert_optimal_termination,
)
from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver
from idaes.core.util.initialization import propagate_state

import idaes.core.util.scaling as iscale
from idaes.generic_models.unit_models import Mixer, Separator, Product, Feed
from idaes.generic_models.unit_models.mixer import MomentumMixingType
from idaes.generic_models.unit_models.translator import Translator
from idaes.generic_models.costing import UnitModelCostingBlock
from idaes.core.util.exceptions import ConfigurationError

from watertap.unit_models.pressure_exchanger import PressureExchanger
from watertap.unit_models.pressure_changer import Pump, EnergyRecoveryDevice
from watertap.core.util.initialization import assert_degrees_of_freedom

import watertap.property_models.seawater_prop_pack as prop_SW
from watertap.unit_models.reverse_osmosis_0D import (
    ReverseOsmosis0D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
from watertap.core.wt_database import Database
import watertap.core.zero_order_properties as prop_ZO
from watertap.unit_models.zero_order import (
    FeedZO,
    PumpElectricityZO,
    NanofiltrationZO,
    SecondaryTreatmentWWTPZO,
)
from watertap.core.zero_order_costing import ZeroOrderCosting
from watertap.costing import WaterTAPCosting


def main():
    m = build()
    set_operating_conditions(m)

    assert_degrees_of_freedom(m, 0)
    assert_units_consistent(m)

    initialize_system(m)
    assert_degrees_of_freedom(m, 0)

    results = solve(m)

    add_costing(m)
    initialize_costing(m)
    assert_degrees_of_freedom(m, 0)  # ensures problem is square

    optimize_operation(m)  # unfixes specific variables for cost optimization

    solve(m)
    display_results(m)
    display_costing(m)

    return m, results


def build():
    # flowsheet set up
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(default={"dynamic": False})

    # define property packages
    m.fs.prop_nf = prop_ZO.WaterParameterBlock(default={"solute_list": ["dye", "tds"]})
    m.fs.prop_ro = prop_SW.SeawaterParameterBlock()

    # define blocks
    prtrt = m.fs.pretreatment = Block()
    dye_sep = m.fs.dye_separation = Block()
    desal = m.fs.desalination = Block()

    # define flowsheet inlets and outlets
    m.fs.feed = FeedZO(default={"property_package": m.fs.prop_nf})
    m.fs.wwt_retentate = Product(default={"property_package": m.fs.prop_nf})
    m.fs.dye_retentate = Product(default={"property_package": m.fs.prop_nf})
    m.fs.permeate = Product(default={"property_package": m.fs.prop_ro})
    m.fs.brine = Product(default={"property_package": m.fs.prop_ro})

    # pretreatment
    prtrt.wwtp = SecondaryTreatmentWWTPZO(
        default={
            "property_package": m.fs.prop_nf,
            "database": m.db,
            "process_subtype": "default",
        }
    )

    # nanofiltration components
    dye_sep.P1 = PumpElectricityZO(
        default={
            "property_package": m.fs.prop_nf,
            "database": m.db,
            "process_subtype": "default",
        }
    )

    dye_sep.nanofiltration = NanofiltrationZO(
        default={
            "property_package": m.fs.prop_nf,
            "database": m.db,
            "process_subtype": "rHGO_dye_rejection",
        }
    )

    # reverse osmosis components

    desal.P2 = Pump(default={"property_package": m.fs.prop_ro})
    desal.RO = ReverseOsmosis0D(
        default={
            "property_package": m.fs.prop_ro,
            "has_pressure_change": True,
            "pressure_change_type": PressureChangeType.calculated,
            "mass_transfer_coefficient": MassTransferCoefficient.calculated,
            "concentration_polarization_type": ConcentrationPolarizationType.calculated,
        }
    )

    desal.RO.width.setub(2000)
    desal.RO.area.setub(20000)

    desal.S1 = Separator(
        default={"property_package": m.fs.prop_ro, "outlet_list": ["P2", "PXR"]}
    )
    desal.M1 = Mixer(
        default={
            "property_package": m.fs.prop_ro,
            "momentum_mixing_type": MomentumMixingType.equality,  # booster pump will match pressure
            "inlet_list": ["P2", "P3"],
        }
    )
    desal.PXR = PressureExchanger(default={"property_package": m.fs.prop_ro})
    desal.P3 = Pump(default={"property_package": m.fs.prop_ro})

    # translator blocks
    m.fs.tb_nf_ro = Translator(
        default={
            "inlet_property_package": m.fs.prop_nf,
            "outlet_property_package": m.fs.prop_ro,
        }
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
    m.fs.s_feed = Arc(source=m.fs.feed.outlet, destination=prtrt.wwtp.inlet)
    prtrt.s01 = Arc(source=prtrt.wwtp.treated, destination=dye_sep.P1.inlet)
    prtrt.s02 = Arc(source=prtrt.wwtp.byproduct, destination=m.fs.wwt_retentate.inlet)
    dye_sep.s01 = Arc(
        source=dye_sep.P1.outlet, destination=dye_sep.nanofiltration.inlet
    )
    dye_sep.s02 = Arc(
        source=dye_sep.nanofiltration.byproduct, destination=m.fs.dye_retentate.inlet
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
    iscale.set_scaling_factor(desal.P2.control_volume.work, 1e-5)
    iscale.set_scaling_factor(desal.RO.area, 1e-4)
    iscale.set_scaling_factor(desal.P3.control_volume.work, 1e-5)
    iscale.set_scaling_factor(desal.PXR.low_pressure_side.work, 1e-5)
    iscale.set_scaling_factor(desal.PXR.high_pressure_side.work, 1e-5)

    # calculate and propagate scaling factors
    iscale.calculate_scaling_factors(m)
    return m


def set_operating_conditions(m):
    prtrt = m.fs.pretreatment
    dye_sep = m.fs.dye_separation
    desal = m.fs.desalination

    # feed
    flow_vol = 120 / 3600 * pyunits.m**3 / pyunits.s
    conc_mass_dye = 2.5 * pyunits.kg / pyunits.m**3
    conc_mass_tds = 50.0 * pyunits.kg / pyunits.m**3
    temperature = 298 * pyunits.K
    pressure = 101325 * pyunits.Pa

    m.fs.feed.flow_vol[0].fix(flow_vol)
    m.fs.feed.conc_mass_comp[0, "dye"].fix(conc_mass_dye)
    m.fs.feed.conc_mass_comp[0, "tds"].fix(conc_mass_tds)
    solve(m.fs.feed)

    # pretreatment
    prtrt.wwtp.load_parameters_from_database(use_default_removal=True)

    # nanofiltration
    dye_sep.nanofiltration.load_parameters_from_database(use_default_removal=True)

    # nf pump
    dye_sep.P1.load_parameters_from_database(use_default_removal=True)
    dye_sep.P1.applied_pressure.fix(
        dye_sep.nanofiltration.applied_pressure.get_values()[0]
    )
    dye_sep.P1.lift_height.unfix()

    # desalination
    desal.P2.efficiency_pump.fix(0.80)
    operating_pressure = 70e5 * pyunits.Pa
    desal.P2.control_volume.properties_out[0].pressure.fix(operating_pressure)
    desal.RO.A_comp.fix(4.2e-12)  # membrane water permeability
    desal.RO.B_comp.fix(3.5e-8)  # membrane salt permeability
    desal.RO.channel_height.fix(1e-3)  # channel height in membrane stage [m]
    desal.RO.spacer_porosity.fix(0.97)  # spacer porosity in membrane stage [-]
    desal.RO.permeate.pressure[0].fix(pressure)  # atmospheric pressure [Pa]
    desal.RO.velocity[0, 0].fix(0.25)
    desal.RO.recovery_vol_phase[0, "Liq"].fix(0.5)
    m.fs.tb_nf_ro.properties_out[0].temperature.fix(temperature)
    m.fs.tb_nf_ro.properties_out[0].pressure.fix(pressure)

    # pressure exchanger
    desal.PXR.efficiency_pressure_exchanger.fix(0.95)
    # booster pump
    desal.P3.efficiency_pump.fix(0.80)

    return


def initialize_system(m):
    prtrt = m.fs.pretreatment
    dye_sep = m.fs.dye_separation
    desal = m.fs.desalination

    # initialize feed
    solve(m.fs.feed)

    # pretreatment
    propagate_state(m.fs.s_feed)
    s = SequentialDecomposition()
    s.options.tear_set = []
    s.options.iterLim = 1
    s.run(prtrt, lambda u: u.initialize())

    # initialized nf
    propagate_state(prtrt.s01)
    seq = SequentialDecomposition()
    seq.options.tear_set = []
    seq.options.iterLim = 1
    seq.run(dye_sep, lambda u: u.initialize())

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
    solve(desal)
    desal.RO.initialize()
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
    desal.RO.velocity[0, 0].unfix()
    desal.RO.velocity[0, 0].setub(0.3)
    desal.RO.velocity[0, 0].setlb(0.1)

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

    m.fs.objective = Objective(expr=m.LCOT)
    return


def solve(blk, solver=None, tee=False, check_termination=True):
    if solver is None:
        solver = get_solver()
    results = solver.solve(blk, tee=tee)
    if check_termination:
        assert_optimal_termination(results)
    return results


def add_costing(m):
    prtrt = m.fs.pretreatment
    dye_sep = m.fs.dye_separation
    desal = m.fs.desalination

    # Zero order costing
    source_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "dye_desalination_global_costing.yaml",
    )

    m.fs.zo_costing = ZeroOrderCosting(default={"case_study_definition": source_file})
    m.fs.ro_costing = WaterTAPCosting()

    # cost nanofiltration module and pump
    prtrt.wwtp.costing = UnitModelCostingBlock(
        default={"flowsheet_costing_block": m.fs.zo_costing}
    )
    dye_sep.nanofiltration.costing = UnitModelCostingBlock(
        default={"flowsheet_costing_block": m.fs.zo_costing}
    )
    dye_sep.P1.costing = UnitModelCostingBlock(
        default={"flowsheet_costing_block": m.fs.zo_costing}
    )

    # RO Train
    # RO equipment is costed using more detailed costing package

    desal.P2.costing = UnitModelCostingBlock(
        default={
            "flowsheet_costing_block": m.fs.ro_costing,
            "costing_method_arguments": {"cost_electricity_flow": True},
        }
    )
    desal.RO.costing = UnitModelCostingBlock(
        default={"flowsheet_costing_block": m.fs.ro_costing}
    )

    desal.M1.costing = UnitModelCostingBlock(
        default={"flowsheet_costing_block": m.fs.ro_costing}
    )
    desal.PXR.costing = UnitModelCostingBlock(
        default={"flowsheet_costing_block": m.fs.ro_costing}
    )
    desal.P3.costing = UnitModelCostingBlock(
        default={
            "flowsheet_costing_block": m.fs.ro_costing,
            "costing_method_arguments": {"cost_electricity_flow": True},
        }
    )

    # Aggregate unit level costs and calculate overall process costs
    m.fs.zo_costing.cost_process()
    m.fs.ro_costing.cost_process()

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

    m.fs.dye_recovery_revenue = Expression(
        expr=(
            m.fs.zo_costing.utilization_factor
            * m.fs.zo_costing.dye_mass_cost
            * pyunits.convert(
                m.fs.dye_retentate.flow_mass_comp[0, "dye"],
                to_units=pyunits.kg / m.fs.zo_costing.base_period,
            )
        ),
        doc="Savings from dye recovered back to the plant",
    )

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
    @m.Expression()
    def total_capital_cost(b, doc="Total capital cost of the treatment train"):
        return pyunits.convert(
            m.fs.zo_costing.total_capital_cost, to_units=pyunits.USD_2020
        ) + pyunits.convert(
            m.fs.ro_costing.total_investment_cost, to_units=pyunits.USD_2020
        )

    @m.Expression()
    def total_operating_cost(b, doc="Total operating cost of the treatment train"):
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

    @m.Expression()
    def total_externalities(
        b, doc="Total cost of water/dye recovered and brine/sludge disposed"
    ):
        return pyunits.convert(
            m.fs.water_recovery_revenue
            + m.fs.dye_recovery_revenue
            - m.fs.brine_disposal_cost
            - m.fs.sludge_disposal_cost,
            to_units=pyunits.USD_2020 / pyunits.year,
        )

    @m.Expression()
    def LCOT(b, doc="Levelized cost of treatment with respect to volumetric feed flow"):
        return (
            b.total_capital_cost * b.fs.zo_costing.capital_recovery_factor
            + b.total_operating_cost
            - b.total_externalities
        ) / (
            pyunits.convert(
                b.fs.feed.properties[0].flow_vol,
                to_units=pyunits.m**3 / pyunits.year,
            )
            * b.fs.zo_costing.utilization_factor
        )

    @m.Expression()
    def LCOW(b, doc="Levelized cost of water with respect to volumetric permeate flow"):
        return (
            b.total_capital_cost * b.fs.zo_costing.capital_recovery_factor
            + b.total_operating_cost
            - pyunits.convert(
                m.fs.dye_recovery_revenue
                - m.fs.brine_disposal_cost
                - m.fs.sludge_disposal_cost,
                to_units=pyunits.USD_2020 / pyunits.year,
            )
        ) / (
            pyunits.convert(
                b.fs.permeate.properties[0].flow_vol,
                to_units=pyunits.m**3 / pyunits.year,
            )
            * b.fs.zo_costing.utilization_factor
        )

    assert_units_consistent(m)
    return


def initialize_costing(m):
    m.fs.zo_costing.initialize()
    m.fs.ro_costing.initialize()
    return


def display_results(m):
    print("Unit models:")
    m.fs.pretreatment.wwtp.report()
    m.fs.dye_separation.P1.report()
    m.fs.dye_separation.nanofiltration.report()
    m.fs.desalination.S1.report()
    m.fs.desalination.P2.report()
    m.fs.desalination.P3.report()
    m.fs.desalination.M1.report()
    m.fs.desalination.RO.report()
    m.fs.desalination.PXR.report()

    print("Streams:")
    flow_list = ["feed", "wwt_retentate", "dye_retentate", "permeate", "brine"]
    for f in flow_list:
        m.fs.component(f).report()
    return


def display_costing(m):
    print("\n System costing metrics:")
    capex = value(pyunits.convert(m.total_capital_cost, to_units=pyunits.MUSD_2020))

    opex = value(
        pyunits.convert(
            m.total_operating_cost, to_units=pyunits.MUSD_2020 / pyunits.year
        )
    )

    externalities = value(
        pyunits.convert(
            m.total_externalities, to_units=pyunits.MUSD_2020 / pyunits.year
        )
    )

    lcot = value(pyunits.convert(m.LCOT, to_units=pyunits.USD_2020 / pyunits.m**3))

    lcow = value(pyunits.convert(m.LCOW, to_units=pyunits.USD_2020 / pyunits.m**3))

    print(f"Total Capital Cost: {capex:.4f} M$")

    print(f"Total Operating Cost: {opex:.4f} M$/year")

    print(f"Total Externalities: {externalities:.4f} M$/year")

    print(f"Levelized cost of treatment: {lcot:.4f} $/m3 feed")

    print(f"Levelized cost of water: {lcow:.4f} $/m3 permeate")


if __name__ == "__main__":
    model, results = main()
