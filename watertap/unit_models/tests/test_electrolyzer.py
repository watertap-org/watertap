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

import pytest
import pyomo.environ as pyo

from pyomo.network import Port
from pyomo.util.check_units import assert_units_consistent
from idaes.core import (
    FlowsheetBlock,
    EnergyBalanceType,
    MaterialBalanceType,
    MomentumBalanceType,
    UnitModelCostingBlock,
)
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    badly_scaled_var_generator,
)
from idaes.core.util.testing import initialization_tester
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from watertap.unit_models.electrolyzer import Electrolyzer
from watertap.costing import WaterTAPCosting

__author__ = "Hunter Barber"

solver = get_solver()

# -----------------------------------------------------------------------------
class TestElectrolyzer:
    @pytest.fixture(scope="class")
    def chlor_alkali_elec(self):

        m = pyo.ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = MCASParameterBlock(
            solute_list=[
                "NA+",
                "CL-",
                "CL2-v",
                "H2-v",
                "OH-",
            ],
            mw_data={
                "H2O": 0.018015,
                "NA+": 0.022989,
                "CL-": 0.03545,
                "CL2-v": 0.0709,
                "H2-v": 0.002016,
                "OH-": 0.017007,
            },
            charge={"NA+": 1, "CL-": -1, "OH-": -1},
            ignore_neutral_charge=True,
        )

        m.fs.unit = Electrolyzer(
            property_package=m.fs.properties,
        )

        # fix property parameters
        m.fs.properties.dens_mass_const = 1200

        # feed specifications
        anolyte_blk = m.fs.unit.anolyte
        catholyte_blk = m.fs.unit.catholyte
        anolyte_blk.properties_in[0].pressure.fix(101325)
        anolyte_blk.properties_in[0].temperature.fix(273.15 + 90)
        anolyte_blk.properties_in[0].flow_mol_phase_comp["Liq", "H2O"].fix(5.551)
        anolyte_blk.properties_in[0].flow_mol_phase_comp["Liq", "NA+"].fix(0.3422)
        anolyte_blk.properties_in[0].flow_mol_phase_comp["Liq", "CL-"].fix(0.3422)
        anolyte_blk.properties_in[0].flow_mol_phase_comp["Liq", "CL2-v"].fix(0)
        anolyte_blk.properties_in[0].flow_mol_phase_comp["Liq", "H2-v"].fix(0)
        anolyte_blk.properties_in[0].flow_mol_phase_comp["Liq", "OH-"].fix(0)
        catholyte_blk.properties_in[0].pressure.fix(101325)
        catholyte_blk.properties_in[0].temperature.fix(273.15 + 90)
        catholyte_blk.properties_in[0].flow_mol_phase_comp["Liq", "H2O"].fix(5.551)
        catholyte_blk.properties_in[0].flow_mol_phase_comp["Liq", "NA+"].fix(1.288)
        catholyte_blk.properties_in[0].flow_mol_phase_comp["Liq", "CL-"].fix(0)
        catholyte_blk.properties_in[0].flow_mol_phase_comp["Liq", "CL2-v"].fix(0)
        catholyte_blk.properties_in[0].flow_mol_phase_comp["Liq", "H2-v"].fix(0)
        catholyte_blk.properties_in[0].flow_mol_phase_comp["Liq", "OH-"].fix(1.288)

        # touch properties
        anolyte_blk.properties_in[0].flow_vol_phase
        anolyte_blk.properties_in[0].conc_mass_phase_comp
        anolyte_blk.properties_in[0].conc_mol_phase_comp
        catholyte_blk.properties_in[0].flow_vol_phase
        catholyte_blk.properties_in[0].conc_mass_phase_comp
        catholyte_blk.properties_in[0].conc_mol_phase_comp
        anolyte_blk.properties_out[0].flow_vol_phase
        anolyte_blk.properties_out[0].conc_mass_phase_comp
        anolyte_blk.properties_out[0].conc_mol_phase_comp
        catholyte_blk.properties_out[0].flow_vol_phase
        catholyte_blk.properties_out[0].conc_mass_phase_comp
        catholyte_blk.properties_out[0].conc_mol_phase_comp

        return m

    @pytest.mark.unit
    def test_config(self, chlor_alkali_elec):

        m = chlor_alkali_elec
        u_config = m.fs.unit.config

        # check unit config arguments
        assert len(u_config) == 8
        assert not u_config.dynamic
        assert not u_config.has_holdup
        assert u_config.material_balance_type == MaterialBalanceType.useDefault
        assert u_config.energy_balance_type == EnergyBalanceType.none
        assert u_config.momentum_balance_type == MomentumBalanceType.pressureTotal

        # check properties
        assert u_config.property_package is m.fs.properties
        assert len(u_config.property_package.solute_set) == 5
        assert len(u_config.property_package.solvent_set) == 1

    @pytest.mark.unit
    def test_build(self, chlor_alkali_elec):

        m = chlor_alkali_elec

        # test units
        assert assert_units_consistent(m) is None

        # test ports
        port_lst = [
            "anolyte_inlet",
            "anolyte_outlet",
            "catholyte_inlet",
            "catholyte_outlet",
        ]
        for port_str in port_lst:
            port = getattr(m.fs.unit, port_str)
            assert len(port.vars) == 3
            assert isinstance(port, Port)

        # test statistics
        assert number_variables(m) == 212
        assert number_total_constraints(m) == 153
        assert number_unused_variables(m) == 15

    @pytest.mark.unit
    def test_dof(self, chlor_alkali_elec):

        m = chlor_alkali_elec

        # test initial degrees of freedom
        assert degrees_of_freedom(m) == 10

        # fix electrolysis reaction variables, TODO: transfer the following variables to generic importable blocks
        # membrane properties
        m.fs.unit.membrane_ion_transport_number["Liq", "NA+"].fix(1)
        # anode properties, Cl- --> 0.5 Cl2 + e-
        m.fs.unit.anode_electrochem_potential.fix(1.21)
        m.fs.unit.anode_stoich["Liq", "CL-"].fix(-1)
        m.fs.unit.anode_stoich["Liq", "CL2-v"].fix(0.5)
        # cathode properties, H20 + e- --> 0.5 H2 + OH-
        m.fs.unit.cathode_electrochem_potential.fix(-0.99)
        m.fs.unit.cathode_stoich["Liq", "H2O"].fix(-1)
        m.fs.unit.cathode_stoich["Liq", "H2-v"].fix(0.5)
        m.fs.unit.cathode_stoich["Liq", "OH-"].fix(1)

        # fix design and performance variables
        # membrane properties
        m.fs.unit.membrane_current_density.fix(4000)
        # anode properties
        m.fs.unit.anode_current_density.fix(3000)
        m.fs.unit.anode_overpotential.fix(0.1)  # assumed
        # cathode properties
        m.fs.unit.cathode_current_density.fix(3000)
        m.fs.unit.cathode_overpotential.fix(0.1)  # assumed
        # electrolyzer cell design
        m.fs.unit.current.fix(30000)
        # performance variables
        m.fs.unit.efficiency_current.fix(0.9)
        m.fs.unit.efficiency_voltage.fix(0.8)

        # test degrees of freedom satisfied
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_init(self, chlor_alkali_elec):

        m = chlor_alkali_elec

        # set default scaling factors
        prop = m.fs.properties
        prop.set_default_scaling("flow_mol_phase_comp", 1, index=("Liq", "H2O"))
        prop.set_default_scaling("flow_mol_phase_comp", 1, index=("Liq", "NA+"))
        prop.set_default_scaling("flow_mol_phase_comp", 1, index=("Liq", "CL-"))
        prop.set_default_scaling("flow_mol_phase_comp", 1, index=("Liq", "CL2-v"))
        prop.set_default_scaling("flow_mol_phase_comp", 1, index=("Liq", "H2-v"))
        prop.set_default_scaling("flow_mol_phase_comp", 1, index=("Liq", "OH-"))

        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        assert len(list(unscaled_variables_generator(m))) == 0

        # test initialization
        initialization_tester(m)

        # check variable scaling
        assert len(list(badly_scaled_var_generator(m, zero=1e-6))) == 0

    @pytest.mark.component
    def test_solve(self, chlor_alkali_elec):

        m = chlor_alkali_elec
        results = solver.solve(m)

        # check for optimal solution
        assert pyo.check_optimal_termination(results)

        # re-check variable scaling post solve
        assert len(list(badly_scaled_var_generator(m, zero=1e-6))) == 0

    @pytest.mark.component
    def test_solution(self, chlor_alkali_elec):

        m = chlor_alkali_elec

        # test report
        m.fs.unit.report()

        # check solution values
        assert pytest.approx(7.500, rel=1e-3) == pyo.value(m.fs.unit.membrane_area)
        assert pytest.approx(10.00, rel=1e-3) == pyo.value(m.fs.unit.anode_area)
        assert pytest.approx(10.00, rel=1e-3) == pyo.value(m.fs.unit.cathode_area)
        assert pytest.approx(2.750, rel=1e-3) == pyo.value(m.fs.unit.voltage_cell)
        assert pytest.approx(1.167e-5, rel=1e-3) == pyo.value(m.fs.unit.resistance)
        assert pytest.approx(82510, rel=1e-3) == pyo.value(m.fs.unit.power)
        assert pytest.approx(2.200, rel=1e-3) == pyo.value(m.fs.unit.voltage_reversible)
        assert pytest.approx(0.2798, rel=1e-3) == pyo.value(m.fs.unit.electron_flow)
        assert pytest.approx(0.7200, rel=1e-3) == pyo.value(m.fs.unit.efficiency_power)

        # check flow at outlet
        assert pytest.approx(0.1399, rel=1e-3) == pyo.value(
            m.fs.unit.anolyte.properties_out[0].flow_mol_phase_comp["Liq", "CL2-v"]
        )
        assert pytest.approx(1.568, rel=1e-3) == pyo.value(
            m.fs.unit.catholyte.properties_out[0].flow_mol_phase_comp["Liq", "NA+"]
        )
        assert pytest.approx(1.568, rel=1e-3) == pyo.value(
            m.fs.unit.catholyte.properties_out[0].flow_mol_phase_comp["Liq", "OH-"]
        )

        # check charge balance
        m.fs.unit.anolyte.properties_in[0].assert_electroneutrality(
            tee=True,
            tol=1e-5,
            solve=False,
        )
        m.fs.unit.anolyte.properties_out[0].assert_electroneutrality(
            tee=True,
            tol=1e-5,
            solve=False,
            defined_state=False,
        )
        m.fs.unit.catholyte.properties_in[0].assert_electroneutrality(
            tee=True,
            tol=1e-5,
            solve=False,
        )
        m.fs.unit.catholyte.properties_out[0].assert_electroneutrality(
            tee=True,
            tol=1e-5,
            solve=False,
            defined_state=False,
        )

    @pytest.mark.unit
    def test_costing_build(self, chlor_alkali_elec):

        m = chlor_alkali_elec

        # build costing model block
        m.fs.costing = WaterTAPCosting()
        m.fs.costing.base_currency = pyo.units.USD_2020
        m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()

        # check general config after building costing block
        assert assert_units_consistent(m) is None
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_costing_results(self, chlor_alkali_elec):

        m = chlor_alkali_elec

        # scale, initialize, and solve with costing
        calculate_scaling_factors(m)
        m.fs.unit.costing.initialize()
        results = solver.solve(m)

        # check for optimal solution
        assert pyo.check_optimal_termination(results)

        # check solution values
        assert pytest.approx(2.0 * 17930, rel=1e-3) == pyo.value(
            m.fs.unit.costing.capital_cost
        )
        assert pytest.approx(82.50, rel=1e-3) == pyo.value(
            m.fs.costing.aggregate_flow_electricity
        )
        assert pytest.approx(50040, rel=1e-3) == pyo.value(
            m.fs.costing.aggregate_flow_costs["electricity"]
        )
