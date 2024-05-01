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
)

from idaes.core import FlowsheetBlock
from watertap.core.solvers import get_solver

import watertap.property_models.seawater_prop_pack as props
import watertap.property_models.multicomp_aq_sol_prop_pack as props2
from watertap.unit_models.pressure_changer import (
    Pump,
    EnergyRecoveryDevice,
    VariableEfficiency,
)

import idaes.core.util.scaling as iscale
from watertap.unit_models.tests.unit_test_harness import UnitTestHarness

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------


def build_pump():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.SeawaterParameterBlock()

    m.fs.unit = Pump(property_package=m.fs.properties)

    # fully specify system
    feed_flow_mass = 1
    feed_mass_frac_TDS = 0.035
    feed_pressure_in = 1e5
    feed_pressure_out = 5e5
    feed_temperature = 273.15 + 25
    efi_pump = 0.75

    feed_mass_frac_H2O = 1 - feed_mass_frac_TDS
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(
        feed_flow_mass * feed_mass_frac_TDS
    )
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )
    m.fs.unit.inlet.pressure[0].fix(feed_pressure_in)
    m.fs.unit.inlet.temperature[0].fix(feed_temperature)
    m.fs.unit.outlet.pressure[0].fix(feed_pressure_out)
    m.fs.unit.efficiency_pump.fix(efi_pump)

    iscale.set_scaling_factor(m.fs.unit.control_volume.work[0], 1e-2)
    iscale.set_scaling_factor(m.fs.unit.work_fluid[0], 1e-2)
    iscale.calculate_scaling_factors(m.fs.unit)

    return m


def build_energy_recovery_device():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.SeawaterParameterBlock()

    m.fs.unit = EnergyRecoveryDevice(property_package=m.fs.properties)

    # fully specify system
    feed_flow_mass = 1
    feed_mass_frac_TDS = 0.035
    feed_pressure_in = 5e5
    feed_pressure_out = 1e5
    feed_temperature = 273.15 + 25
    efi_pump = 0.75

    feed_mass_frac_H2O = 1 - feed_mass_frac_TDS
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(
        feed_flow_mass * feed_mass_frac_TDS
    )
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )
    m.fs.unit.inlet.pressure[0].fix(feed_pressure_in)
    m.fs.unit.inlet.temperature[0].fix(feed_temperature)
    m.fs.unit.outlet.pressure[0].fix(feed_pressure_out)
    m.fs.unit.efficiency_pump.fix(efi_pump)

    iscale.set_scaling_factor(m.fs.unit.control_volume.work[0], 1e-2)
    iscale.set_scaling_factor(m.fs.unit.work_fluid[0], 1e-2)
    iscale.calculate_scaling_factors(m.fs.unit)

    return m


def build_pump_variable_flow():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.SeawaterParameterBlock()

    m.fs.unit = Pump(
        property_package=m.fs.properties,
        variable_efficiency=VariableEfficiency.flow,
    )
    # fully specify system
    feed_flow_mass = 1
    feed_mass_frac_TDS = 0.035
    feed_pressure_in = 1e5
    feed_pressure_out = 5e5
    feed_temperature = 273.15 + 25
    efi_pump = 0.75

    feed_mass_frac_H2O = 1 - feed_mass_frac_TDS
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(
        feed_flow_mass * feed_mass_frac_TDS
    )
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )
    m.fs.unit.inlet.pressure[0].fix(feed_pressure_in)
    m.fs.unit.inlet.temperature[0].fix(feed_temperature)
    m.fs.unit.outlet.pressure[0].fix(feed_pressure_out)

    m.fs.unit.bep_eta.fix(efi_pump)
    m.fs.unit.flow_ratio.fix(1)

    iscale.set_scaling_factor(m.fs.unit.control_volume.work[0], 1e-2)
    iscale.set_scaling_factor(m.fs.unit.work_fluid[0], 1e-2)
    iscale.set_scaling_factor(m.fs.unit.bep_flow, 1e4)
    iscale.calculate_scaling_factors(m.fs.unit)

    return m


def build_pump_isothermal_energybalancetype_none():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props2.MCASParameterBlock(
        solute_list=["TDS"],
        mw_data={"TDS": 58e-3},
        ignore_neutral_charge=True,
    )

    m.fs.unit = Pump(property_package=m.fs.properties)

    # fully specify system
    feed_flow_mol = 1
    feed_mol_frac_TDS = 0.035
    feed_pressure_in = 1e5
    feed_pressure_out = 5e5
    feed_temperature = 273.15 + 25
    efi_pump = 0.75

    feed_mass_frac_H2O = 1 - feed_mol_frac_TDS
    m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "TDS"].fix(
        feed_flow_mol * feed_mol_frac_TDS
    )
    m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mol * feed_mass_frac_H2O
    )
    m.fs.unit.inlet.pressure[0].fix(feed_pressure_in)
    m.fs.unit.inlet.temperature[0].fix(feed_temperature)
    m.fs.unit.outlet.pressure[0].fix(feed_pressure_out)
    m.fs.unit.efficiency_pump.fix(efi_pump)

    iscale.set_scaling_factor(m.fs.unit.control_volume.work[0], 1e-2)
    iscale.set_scaling_factor(m.fs.unit.work_fluid[0], 1e-2)
    iscale.calculate_scaling_factors(m.fs.unit)

    return m


class TestPump(UnitTestHarness):
    def configure(self):
        m = build_pump()

        self.unit_solutions[m.fs.unit.work_mechanical[0]] = 521.05643
        self.unit_solutions[m.fs.unit.deltaP[0]] = 4e5
        self.unit_solutions[m.fs.unit.ratioP[0]] = 5
        self.unit_solutions[m.fs.unit.work_fluid[0]] = 390.7923225
        self.unit_solutions[m.fs.unit.control_volume.properties_out[0].temperature] = (
            298.15
        )
        self.unit_solutions[m.fs.unit.control_volume.properties_out[0].enth_flow] = (
            9.931834173e4
        )
        self.unit_solutions[
            m.fs.unit.control_volume.properties_out[0].flow_mass_phase_comp[
                "Liq", "H2O"
            ]
        ] = 0.965
        self.unit_solutions[
            m.fs.unit.control_volume.properties_out[0].flow_mass_phase_comp[
                "Liq", "TDS"
            ]
        ] = 0.035

        # Conservation check
        comp_lst = ["TDS", "H2O"]

        flow_mass_inlet = sum(
            m.fs.unit.control_volume.properties_in[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )

        flow_mass_outlet = sum(
            m.fs.unit.control_volume.properties_out[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )

        self.unit_solutions[flow_mass_outlet] = value(flow_mass_inlet)

        return m


class TestEnergyRecoveryDevice(UnitTestHarness):
    def configure(self):
        m = build_energy_recovery_device()

        self.unit_solutions[m.fs.unit.work_mechanical[0]] = -293.094242
        self.unit_solutions[m.fs.unit.deltaP[0]] = -4e5
        self.unit_solutions[m.fs.unit.ratioP[0]] = 0.2
        self.unit_solutions[m.fs.unit.work_fluid[0]] = -390.7923225
        self.unit_solutions[m.fs.unit.control_volume.properties_out[0].temperature] = (
            298.15
        )
        self.unit_solutions[m.fs.unit.control_volume.properties_out[0].enth_flow] = (
            9.89617297e4
        )
        self.unit_solutions[
            m.fs.unit.control_volume.properties_out[0].flow_mass_phase_comp[
                "Liq", "H2O"
            ]
        ] = 0.965
        self.unit_solutions[
            m.fs.unit.control_volume.properties_out[0].flow_mass_phase_comp[
                "Liq", "TDS"
            ]
        ] = 0.035

        # Conservation check
        comp_lst = ["TDS", "H2O"]

        flow_mass_inlet = sum(
            m.fs.unit.control_volume.properties_in[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )

        flow_mass_outlet = sum(
            m.fs.unit.control_volume.properties_out[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )

        self.unit_solutions[flow_mass_outlet] = value(flow_mass_inlet)

        return m


class TestPumpVariableFlow(UnitTestHarness):
    def configure(self):
        m = build_pump_variable_flow()

        self.unit_solutions[m.fs.unit.work_mechanical[0]] = 521.05643
        self.unit_solutions[m.fs.unit.deltaP[0]] = 4e5
        self.unit_solutions[m.fs.unit.ratioP[0]] = 5
        self.unit_solutions[m.fs.unit.work_fluid[0]] = 390.7923225
        self.unit_solutions[m.fs.unit.control_volume.properties_out[0].temperature] = (
            298.15
        )
        self.unit_solutions[m.fs.unit.control_volume.properties_out[0].enth_flow] = (
            9.93183417e4
        )
        self.unit_solutions[
            m.fs.unit.control_volume.properties_out[0].flow_mass_phase_comp[
                "Liq", "H2O"
            ]
        ] = 0.965
        self.unit_solutions[
            m.fs.unit.control_volume.properties_out[0].flow_mass_phase_comp[
                "Liq", "TDS"
            ]
        ] = 0.035

        # Conservation check
        comp_lst = ["TDS", "H2O"]

        flow_mass_inlet = sum(
            m.fs.unit.control_volume.properties_in[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )

        flow_mass_outlet = sum(
            m.fs.unit.control_volume.properties_out[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )

        self.unit_solutions[flow_mass_outlet] = value(flow_mass_inlet)

        return m


class TestPumpIsothermal(UnitTestHarness):
    def configure(self):
        m = build_pump_isothermal_energybalancetype_none()

        self.unit_solutions[m.fs.unit.work_mechanical[0]] = 10.3466667
        self.unit_solutions[m.fs.unit.deltaP[0]] = 4e5
        self.unit_solutions[m.fs.unit.ratioP[0]] = 5
        self.unit_solutions[m.fs.unit.work_fluid[0]] = 7.76
        self.unit_solutions[m.fs.unit.control_volume.properties_out[0].temperature] = (
            298.15
        )
        self.unit_solutions[
            m.fs.unit.control_volume.properties_out[0].flow_mass_phase_comp[
                "Liq", "H2O"
            ]
        ] = 0.01737
        self.unit_solutions[
            m.fs.unit.control_volume.properties_out[0].flow_mass_phase_comp[
                "Liq", "TDS"
            ]
        ] = 0.00203

        # Conservation check
        comp_lst = ["TDS", "H2O"]

        flow_mass_inlet = sum(
            m.fs.unit.control_volume.properties_in[0].flow_mol_phase_comp["Liq", j]
            for j in comp_lst
        )

        flow_mass_outlet = sum(
            m.fs.unit.control_volume.properties_out[0].flow_mol_phase_comp["Liq", j]
            for j in comp_lst
        )

        self.unit_solutions[flow_mass_outlet] = value(flow_mass_inlet)

        return m
