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
"""
Tests for CSTR unit model with injection.
Authors: Andrew Lee, Adam Atia, Vibhav Dabadghao
"""

import pytest
from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    units,
    value,
    assert_optimal_termination,
)
from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)
from idaes.models.properties.examples.saponification_thermo import (
    SaponificationParameterBlock,
)
from idaes.models.properties.examples.saponification_reactions import (
    SaponificationReactionParameterBlock,
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import (
    PhysicalParameterTestBlock,
    ReactionParameterTestBlock,
    initialization_tester,
)
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.solvers import get_solver
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent

from watertap.unit_models.cstr_injection import CSTR_Injection, ElectricityConsumption
from idaes.core import UnitModelCostingBlock
from watertap.costing import WaterTAPCosting
from watertap.property_models.activated_sludge.asm1_properties import ASM1ParameterBlock
from watertap.property_models.activated_sludge.asm1_reactions import (
    ASM1ReactionParameterBlock,
)
from watertap.property_models.activated_sludge.asm2d_properties import (
    ASM2dParameterBlock,
)
from watertap.property_models.activated_sludge.asm2d_reactions import (
    ASM2dReactionParameterBlock,
)
from watertap.property_models.activated_sludge.modified_asm2d_properties import (
    ModifiedASM2dParameterBlock,
)
from watertap.property_models.activated_sludge.modified_asm2d_reactions import (
    ModifiedASM2dReactionParameterBlock,
)

from watertap.property_models.anaerobic_digestion.adm1_properties import (
    ADM1ParameterBlock,
)
from watertap.property_models.anaerobic_digestion.adm1_reactions import (
    ADM1ReactionParameterBlock,
)

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = PhysicalParameterTestBlock()
    m.fs.reactions = ReactionParameterTestBlock(property_package=m.fs.properties)

    m.fs.unit = CSTR_Injection(
        property_package=m.fs.properties, reaction_package=m.fs.reactions
    )

    # Check unit config arguments
    assert len(m.fs.unit.config) == 16

    assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
    assert m.fs.unit.config.energy_balance_type == EnergyBalanceType.useDefault
    assert m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
    assert not m.fs.unit.config.has_heat_transfer
    assert not m.fs.unit.config.has_pressure_change
    assert not m.fs.unit.config.has_equilibrium_reactions
    assert not m.fs.unit.config.has_phase_equilibrium
    assert not m.fs.unit.config.has_heat_of_reaction
    assert m.fs.unit.config.property_package is m.fs.properties
    assert m.fs.unit.config.reaction_package is m.fs.reactions
    assert m.fs.unit.config.electricity_consumption == ElectricityConsumption.fixed
    assert not m.fs.unit.config.has_aeration


# -----------------------------------------------------------------------------
class TestSaponification(object):
    @pytest.fixture(scope="class")
    def sapon(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = SaponificationParameterBlock()
        m.fs.reactions = SaponificationReactionParameterBlock(
            property_package=m.fs.properties
        )

        m.fs.unit = CSTR_Injection(
            property_package=m.fs.properties,
            reaction_package=m.fs.reactions,
            has_equilibrium_reactions=False,
            has_heat_transfer=True,
            has_heat_of_reaction=True,
            has_pressure_change=True,
        )

        m.fs.unit.inlet.flow_vol.fix(1.0e-03)
        m.fs.unit.inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
        m.fs.unit.inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
        m.fs.unit.inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
        m.fs.unit.inlet.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
        m.fs.unit.inlet.conc_mol_comp[0, "Ethanol"].fix(0.0)

        m.fs.unit.inlet.temperature.fix(303.15)
        m.fs.unit.inlet.pressure.fix(101325.0)

        m.fs.unit.volume.fix(1.5e-03)
        m.fs.unit.heat_duty.fix(0)
        m.fs.unit.deltaP.fix(0)
        m.fs.unit.injection.fix(0)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, sapon):
        assert hasattr(sapon.fs.unit, "inlet")
        assert len(sapon.fs.unit.inlet.vars) == 4
        assert hasattr(sapon.fs.unit.inlet, "flow_vol")
        assert hasattr(sapon.fs.unit.inlet, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.inlet, "temperature")
        assert hasattr(sapon.fs.unit.inlet, "pressure")

        assert hasattr(sapon.fs.unit, "outlet")
        assert len(sapon.fs.unit.outlet.vars) == 4
        assert hasattr(sapon.fs.unit.outlet, "flow_vol")
        assert hasattr(sapon.fs.unit.outlet, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.outlet, "temperature")
        assert hasattr(sapon.fs.unit.outlet, "pressure")

        assert hasattr(sapon.fs.unit, "cstr_performance_eqn")
        assert hasattr(sapon.fs.unit, "volume")
        assert hasattr(sapon.fs.unit, "hydraulic_retention_time")
        assert hasattr(sapon.fs.unit, "heat_duty")
        assert hasattr(sapon.fs.unit, "deltaP")

        assert hasattr(sapon.fs.unit, "injection")
        for k in sapon.fs.unit.injection:
            assert (
                sapon.fs.unit.injection[k]
                == sapon.fs.unit.control_volume.mass_transfer_term[k]
            )

        assert number_variables(sapon) == 34
        assert number_total_constraints(sapon) == 18
        assert number_unused_variables(sapon) == 0

    @pytest.mark.component
    def test_units(self, sapon):
        assert_units_consistent(sapon)
        assert_units_equivalent(sapon.fs.unit.volume[0], units.m**3)
        assert_units_equivalent(sapon.fs.unit.heat_duty[0], units.W)
        assert_units_equivalent(sapon.fs.unit.deltaP[0], units.Pa)

    @pytest.mark.unit
    def test_dof(self, sapon):
        assert degrees_of_freedom(sapon) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, sapon):
        initialization_tester(sapon)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, sapon):
        results = solver.solve(sapon)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, sapon):
        assert pytest.approx(101325.0, abs=1e-2) == value(
            sapon.fs.unit.outlet.pressure[0]
        )
        assert pytest.approx(304.09, abs=1e-2) == value(
            sapon.fs.unit.outlet.temperature[0]
        )
        assert pytest.approx(20.32, abs=1e-2) == value(
            sapon.fs.unit.outlet.conc_mol_comp[0, "EthylAcetate"]
        )
        assert pytest.approx(0.03888, abs=1e-3) == value(
            sapon.fs.unit.electricity_consumption[0]
        )
        assert pytest.approx(1.5, abs=1e-3) == value(
            sapon.fs.unit.hydraulic_retention_time[0]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, sapon):
        assert (
            abs(
                value(
                    sapon.fs.unit.inlet.flow_vol[0] - sapon.fs.unit.outlet.flow_vol[0]
                )
            )
            <= 1e-6
        )
        assert (
            abs(
                value(
                    sapon.fs.unit.inlet.flow_vol[0]
                    * sum(
                        sapon.fs.unit.inlet.conc_mol_comp[0, j]
                        for j in sapon.fs.properties.component_list
                    )
                    - sapon.fs.unit.outlet.flow_vol[0]
                    * sum(
                        sapon.fs.unit.outlet.conc_mol_comp[0, j]
                        for j in sapon.fs.properties.component_list
                    )
                )
            )
            <= 1e-6
        )

        assert pytest.approx(3904.51, abs=1e-2) == value(
            sapon.fs.unit.control_volume.heat_of_reaction[0]
        )
        assert (
            abs(
                value(
                    (
                        sapon.fs.unit.inlet.flow_vol[0]
                        * sapon.fs.properties.dens_mol
                        * sapon.fs.properties.cp_mol
                        * (
                            sapon.fs.unit.inlet.temperature[0]
                            - sapon.fs.properties.temperature_ref
                        )
                    )
                    - (
                        sapon.fs.unit.outlet.flow_vol[0]
                        * sapon.fs.properties.dens_mol
                        * sapon.fs.properties.cp_mol
                        * (
                            sapon.fs.unit.outlet.temperature[0]
                            - sapon.fs.properties.temperature_ref
                        )
                    )
                    + sapon.fs.unit.control_volume.heat_of_reaction[0]
                )
            )
            <= 1e-3
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_costing(self, sapon):
        m = sapon

        m.fs.costing = WaterTAPCosting()

        m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(m.fs.unit.control_volume.properties_out[0].flow_vol)
        solver = get_solver()
        results = solver.solve(m)

        assert_optimal_termination(results)

        # Check solutions
        assert pytest.approx(7.75429 * 2, rel=1e-5) == value(
            m.fs.unit.costing.capital_cost
        )
        assert pytest.approx(0.00082698, rel=1e-5) == value(m.fs.costing.LCOW)

    @pytest.mark.unit
    def test_get_performance_contents(self, sapon):
        perf_dict = sapon.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "Volume": sapon.fs.unit.volume[0],
                "Injection [('Liq', 'H2O')]": sapon.fs.unit.injection[0, "Liq", "H2O"],
                "Injection [('Liq', 'NaOH')]": sapon.fs.unit.injection[
                    0, "Liq", "NaOH"
                ],
                "Injection [('Liq', 'EthylAcetate')]": sapon.fs.unit.injection[
                    0, "Liq", "EthylAcetate"
                ],
                "Injection [('Liq', 'SodiumAcetate')]": sapon.fs.unit.injection[
                    0, "Liq", "SodiumAcetate"
                ],
                "Injection [('Liq', 'Ethanol')]": sapon.fs.unit.injection[
                    0, "Liq", "Ethanol"
                ],
                "Heat Duty": sapon.fs.unit.heat_duty[0],
                "Pressure Change": sapon.fs.unit.deltaP[0],
                "Electricity Consumption": sapon.fs.unit.electricity_consumption[0],
            }
        }


class TestCSTR_withASM1(object):
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = ASM1ParameterBlock()
        m.fs.reactions = ASM1ReactionParameterBlock(property_package=m.fs.properties)

        m.fs.unit = CSTR_Injection(
            property_package=m.fs.properties,
            reaction_package=m.fs.reactions,
            has_aeration=True,
            electricity_consumption=ElectricityConsumption.calculated,
        )

        m.fs.unit.inlet.flow_vol.fix(20648 * units.m**3 / units.day)
        m.fs.unit.inlet.temperature.fix(308.15 * units.K)
        m.fs.unit.inlet.pressure.fix(1 * units.atm)
        m.fs.unit.inlet.conc_mass_comp[0, "S_I"].fix(27 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_S"].fix(58 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_I"].fix(92 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_S"].fix(363 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_BH"].fix(50 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_BA"].fix(0 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_P"].fix(0 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_O"].fix(0 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_NO"].fix(0 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_NH"].fix(23 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_ND"].fix(5 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_ND"].fix(16 * units.g / units.m**3)
        m.fs.unit.inlet.alkalinity.fix(7 * units.mol / units.m**3)

        m.fs.unit.volume.fix(500)
        m.fs.unit.injection.fix(0)
        m.fs.unit.injection[0, "Liq", "S_O"].fix(2e-3)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, model):
        assert hasattr(model.fs.unit, "inlet")
        assert len(model.fs.unit.inlet.vars) == 5
        assert hasattr(model.fs.unit.inlet, "flow_vol")
        assert hasattr(model.fs.unit.inlet, "conc_mass_comp")
        assert hasattr(model.fs.unit.inlet, "temperature")
        assert hasattr(model.fs.unit.inlet, "pressure")
        assert hasattr(model.fs.unit.inlet, "alkalinity")

        assert hasattr(model.fs.unit, "outlet")
        assert len(model.fs.unit.inlet.vars) == 5
        assert hasattr(model.fs.unit.inlet, "flow_vol")
        assert hasattr(model.fs.unit.inlet, "conc_mass_comp")
        assert hasattr(model.fs.unit.inlet, "temperature")
        assert hasattr(model.fs.unit.inlet, "pressure")
        assert hasattr(model.fs.unit.inlet, "alkalinity")

        assert hasattr(model.fs.unit, "cstr_performance_eqn")
        assert hasattr(model.fs.unit, "volume")
        assert hasattr(model.fs.unit, "hydraulic_retention_time")
        assert hasattr(model.fs.unit, "KLa")
        assert hasattr(model.fs.unit, "S_O_eq")

        assert hasattr(model.fs.unit, "injection")
        for k in model.fs.unit.injection:
            assert (
                model.fs.unit.injection[k]
                == model.fs.unit.control_volume.mass_transfer_term[k]
            )

        assert number_variables(model) == 102
        assert number_total_constraints(model) == 49
        assert number_unused_variables(model) == 3

    @pytest.mark.component
    def test_units(self, model):
        assert_units_consistent(model)
        assert_units_equivalent(model.fs.unit.volume[0], units.m**3)
        assert_units_equivalent(model.fs.unit.electricity_consumption[0], units.kW)

    @pytest.mark.unit
    def test_dof(self, model):
        assert degrees_of_freedom(model) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, model):
        initialization_tester(model)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        assert pytest.approx(101325.0, abs=1e-2) == value(
            model.fs.unit.outlet.pressure[0]
        )
        assert pytest.approx(308.15, abs=1e-2) == value(
            model.fs.unit.outlet.temperature[0]
        )
        assert pytest.approx(2.197e-3, abs=1e-2) == value(
            model.fs.unit.outlet.conc_mass_comp[0, "S_O"]
        )
        assert pytest.approx(18.3765, abs=1e-3) == value(
            model.fs.unit.electricity_consumption[0]
        )
        assert pytest.approx(2092.2123, abs=1e-3) == value(
            model.fs.unit.hydraulic_retention_time[0]
        )
        assert pytest.approx(8.2694, abs=1e-3) == value(model.fs.unit.KLa)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, model):
        assert (
            abs(
                value(
                    model.fs.unit.inlet.flow_vol[0] - model.fs.unit.outlet.flow_vol[0]
                )
            )
            <= 1e-6
        )
        assert abs(
            value(
                model.fs.unit.inlet.flow_vol[0]
                * sum(
                    model.fs.unit.inlet.conc_mass_comp[0, j]
                    for j in model.fs.properties.solute_set
                )
                - model.fs.unit.outlet.flow_vol[0]
                * sum(
                    model.fs.unit.outlet.conc_mass_comp[0, j]
                    for j in model.fs.properties.solute_set
                )
                + sum(
                    model.fs.unit.control_volume.rate_reaction_generation[0, "Liq", j]
                    for j in model.fs.properties.solute_set
                )
                <= 1e-6
            )
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_costing(self, model):
        m = model

        m.fs.costing = WaterTAPCosting()

        m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(m.fs.unit.control_volume.properties_out[0].flow_vol)
        solver = get_solver()
        results = solver.solve(m)

        assert_optimal_termination(results)

        # Check solutions
        assert pytest.approx(613502.91662, rel=1e-5) == value(
            m.fs.unit.costing.capital_cost
        )
        assert pytest.approx(0.0132455, rel=1e-5) == value(m.fs.costing.LCOW)

    @pytest.mark.unit
    def test_report(self, model):
        model.fs.unit.report()


@pytest.mark.build
@pytest.mark.unit
def test_with_asm2d():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = ASM2dParameterBlock()
    m.fs.reactions = ASM2dReactionParameterBlock(property_package=m.fs.properties)

    m.fs.unit = CSTR_Injection(
        property_package=m.fs.properties,
        reaction_package=m.fs.reactions,
        has_aeration=True,
        electricity_consumption=ElectricityConsumption.calculated,
    )


@pytest.mark.build
@pytest.mark.unit
def test_with_mod_asm2d():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = ModifiedASM2dParameterBlock()
    m.fs.reactions = ModifiedASM2dReactionParameterBlock(
        property_package=m.fs.properties
    )

    m.fs.unit = CSTR_Injection(
        property_package=m.fs.properties,
        reaction_package=m.fs.reactions,
        has_aeration=True,
        electricity_consumption=ElectricityConsumption.calculated,
    )


@pytest.mark.unit
def test_error_without_oxygen():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = ADM1ParameterBlock()
    m.fs.reactions = ADM1ReactionParameterBlock(property_package=m.fs.properties)

    # Expect exception if has_aeration=True but S_O or S_O2 not listed in component_list of prop package.
    with pytest.raises(
        ConfigurationError,
        match="has_aeration was set to True, but the property package has neither 'S_O' nor 'S_O2' in its list of components.",
    ):
        m.fs.unit = CSTR_Injection(
            property_package=m.fs.properties,
            reaction_package=m.fs.reactions,
            has_aeration=True,
        )
