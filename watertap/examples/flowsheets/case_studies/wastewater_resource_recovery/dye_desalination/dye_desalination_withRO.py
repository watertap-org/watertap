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
    Block,
    Set,
    Expression,
    value,
    TransformationFactory,
    units as pyunits,
    assert_optimal_termination,
)
from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver
from idaes.core.util.initialization import (
    propagate_state,
    fix_state_vars,
    revert_state_vars,
)
import idaes.core.util.scaling as iscale
from idaes.generic_models.unit_models import Mixer, Separator, Product, Feed
from idaes.generic_models.unit_models.mixer import MomentumMixingType
from idaes.generic_models.unit_models.translator import Translator
from idaes.generic_models.costing import UnitModelCostingBlock
from idaes.core.util.exceptions import ConfigurationError

from watertap.core.util.initialization import assert_degrees_of_freedom

from watertap.unit_models.pressure_exchanger import PressureExchanger
from watertap.unit_models.pressure_changer import Pump, EnergyRecoveryDevice
from watertap.core.util.initialization import check_dof, assert_degrees_of_freedom

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
)
from watertap.core.zero_order_costing import ZeroOrderCosting


def main():
    m = build()
    set_operating_conditions(m)

    assert_degrees_of_freedom(m, 0)
    assert_units_consistent(m)

    initialize_system(m)
    # results = solve(m)
    return m


def build(erd_type="pressure_exchanger"):
    # flowsheet set up
    m = ConcreteModel()
    m.db = Database()
    m.erd_type = erd_type

    m.fs = FlowsheetBlock(default={"dynamic": False})

    # define property packages
    m.fs.prop_nf = prop_ZO.WaterParameterBlock(default={"solute_list": ["dye", "tds"]})
    m.fs.prop_ro = prop_SW.SeawaterParameterBlock()

    # define blocks
    dye_sep = m.fs.dye_separation = Block()
    desal = m.fs.desalination = Block()
    print("\ndefine blocks:")
    check_dof(m)

    # define flowsheet inlets and outlets
    m.fs.feed = FeedZO(default={"property_package": m.fs.prop_nf})
    m.fs.dye_retentate = Product(default={"property_package": m.fs.prop_nf})
    m.fs.permeate = Product(default={"property_package": m.fs.prop_ro})
    m.fs.brine = Product(default={"property_package": m.fs.prop_ro})
    print("\ndefine fs i/o:")
    check_dof(m)

    # nanofiltration components
    dye_sep.P1 = PumpElectricityZO(
        default={
            "property_package": m.fs.prop_nf,
            "database": m.db,
            "process_subtype": "default",
        }
    )
    print("\ndefine pump electricity ZO")
    check_dof(m)

    dye_sep.nanofiltration = NanofiltrationZO(
        default={
            "property_package": m.fs.prop_nf,
            "database": m.db,
            "process_subtype": "rHGO_dye_rejection",
        }
    )
    print("\ndefine NF ZO")
    check_dof(m)
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
    print("\ndefine RO Od")
    check_dof(m)

    desal.RO.width.setub(2000)
    desal.RO.area.setub(20000)

    if erd_type == "pressure_exchanger":
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
    elif erd_type == "pump_as_turbine":
        desal.ERD = EnergyRecoveryDevice(default={"property_package": m.fs.prop_ro})
    else:
        raise ConfigurationError(
            "erd_type was {}, but can only "
            "be pressure_exchanger or pump_as_turbine"
            "".format(erd_type)
        )
    print("\nset PX/erd conditions")
    check_dof(m)

    # translator blocks
    m.fs.tb_nf_ro = Translator(
        default={
            "inlet_property_package": m.fs.prop_nf,
            "outlet_property_package": m.fs.prop_ro,
        }
    )
    print("\ntranslator block")
    check_dof(m)

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

    print("\nAdd constraint between nf/ro")
    check_dof(m)

    # connections
    m.fs.s_feed = Arc(source=m.fs.feed.outlet, destination=dye_sep.P1.inlet)
    dye_sep.s01 = Arc(
        source=dye_sep.P1.outlet, destination=dye_sep.nanofiltration.inlet
    )
    dye_sep.s02 = Arc(
        source=dye_sep.nanofiltration.byproduct, destination=m.fs.dye_retentate.inlet
    )
    m.fs.s_nf = Arc(
        source=dye_sep.nanofiltration.treated, destination=m.fs.tb_nf_ro.inlet
    )

    if erd_type == "pressure_exchanger":
        m.fs.s_ro = Arc(source=m.fs.tb_nf_ro.outlet, destination=desal.S1.inlet)
        desal.s01 = Arc(source=desal.S1.P2, destination=desal.P2.inlet)
        desal.s02 = Arc(source=desal.P2.outlet, destination=desal.M1.P2)
        desal.s03 = Arc(source=desal.M1.outlet, destination=desal.RO.inlet)
        desal.s04 = Arc(
            source=desal.RO.retentate, destination=desal.PXR.high_pressure_inlet
        )
        desal.s05 = Arc(source=desal.S1.PXR, destination=desal.PXR.low_pressure_inlet)
        desal.s06 = Arc(
            source=desal.PXR.low_pressure_outlet, destination=desal.P3.inlet
        )
        desal.s07 = Arc(source=desal.P3.outlet, destination=desal.M1.P3)
        m.fs.s_disposal = Arc(
            source=desal.PXR.high_pressure_outlet, destination=m.fs.brine.inlet
        )
    elif erd_type == "pump_as_turbine":
        m.fs.s_ro = Arc(source=m.fs.tb_nf_ro.outlet, destination=desal.P2.inlet)
        m.fs.s01 = Arc(source=desal.P1.outlet, destination=desal.RO.inlet)
        m.fs.s_disposal = Arc(source=desal.ERD.outlet, destination=m.fs.disposal.inlet)

    m.fs.s_permeate = Arc(source=desal.RO.permeate, destination=m.fs.permeate.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m)
    print("\nExpand network arcs")
    check_dof(m)

    # scaling
    m.fs.prop_ro.set_default_scaling("flow_mass_phase_comp", 1e-3, index=("Liq", "H2O"))
    m.fs.prop_ro.set_default_scaling("flow_mass_phase_comp", 1e-1, index=("Liq", "TDS"))

    # set unit model values
    iscale.set_scaling_factor(desal.P2.control_volume.work, 1e-5)
    iscale.set_scaling_factor(desal.RO.area, 1e-4)
    if erd_type == "pressure_exchanger":
        iscale.set_scaling_factor(desal.P3.control_volume.work, 1e-5)
        iscale.set_scaling_factor(desal.PXR.low_pressure_side.work, 1e-5)
        iscale.set_scaling_factor(desal.PXR.high_pressure_side.work, 1e-5)
    elif erd_type == "pump_as_turbine":
        iscale.set_scaling_factor(desal.ERD.control_volume.work, 1e-5)

    # calculate and propagate scaling factors
    iscale.calculate_scaling_factors(m)
    return m


def set_operating_conditions(m):
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

    # nanofiltration
    dye_sep.nanofiltration.load_parameters_from_database(use_default_removal=True)
    print("\nLoad NF parameters from database")
    check_dof(m)

    # nf pump
    dye_sep.P1.load_parameters_from_database(use_default_removal=True)
    dye_sep.P1.applied_pressure.fix(
        dye_sep.nanofiltration.applied_pressure.get_values()[0]
    )
    dye_sep.P1.lift_height.unfix()
    print("\nLoad/fix low pressure pump parameters")
    check_dof(m)

    # desalination
    desal.P2.efficiency_pump.fix(0.80)
    operating_pressure = 70e5 * pyunits.Pa
    desal.P2.control_volume.properties_out[0].pressure.fix(operating_pressure)
    desal.RO.A_comp.fix(4.2e-12)  # membrane water permeability
    desal.RO.B_comp.fix(3.5e-8)  # membrane salt permeability
    desal.RO.channel_height.fix(1e-3)  # channel height in membrane stage [m]
    desal.RO.spacer_porosity.fix(0.97)  # spacer porosity in membrane stage [-]
    desal.RO.permeate.pressure[0].fix(pressure)  # atmospheric pressure [Pa]
    desal.RO.width.fix(1000)  # stage width [m]
    desal.RO.area.fix(
        flow_vol * 4.5e4 * pyunits.s / pyunits.m
    )  # TODO - verify this value and change as needed based on CP effects
    desal.RO.feed_side.properties_in[0].temperature.fix(temperature)
    desal.RO.feed_side.properties_in[0].pressure = value(
        operating_pressure
    )  # TODO - check if this is the right variable to fix

    print("\nLoad/fix RO parameters")
    check_dof(m)

    if m.erd_type == "pressure_exchanger":
        # pressure exchanger
        desal.PXR.efficiency_pressure_exchanger.fix(0.95)
        # booster pump
        desal.P3.efficiency_pump.fix(0.80)

    elif m.erd_type == "pump_as_turbine":
        desal.ERD.efficiency_pump.fix(0.95)
        desal.ERD.control_volume.properties_out[0].pressure.fix(
            pressure
        )  # atmospheric pressure [Pa]
    print("\nFix PX/ERD parameters")
    check_dof(m)


def initialize_system(m):
    dye_sep = m.fs.dye_separation
    desal = m.fs.desalination

    # initialize feed
    solve(m.fs.feed)

    # initialized nf
    propagate_state(m.fs.s_feed)
    solve(dye_sep)

    # seq = SequentialDecomposition()
    # seq.options.tear_set = []
    # seq.options.iterLim = 1
    # seq.run(m, lambda u: u.initialize())


def solve(blk, solver=None, tee=False, check_termination=True):
    if solver is None:
        solver = get_solver()
    results = solver.solve(blk, tee=tee)
    if check_termination:
        assert_optimal_termination(results)
    return results


if __name__ == "__main__":
    m = main()
