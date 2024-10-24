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
"""
Tests for CSTR unit model.
Authors: Marcus Holly
"""

import pytest

from pyomo.environ import (
    ConcreteModel,
    units,
    value,
    Objective,
)

from watertap.unit_models.tests.unit_test_harness import UnitTestHarness
import idaes.core.util.scaling as iscale

from idaes.core import (
    FlowsheetBlock,
    UnitModelCostingBlock,
)
from watertap.unit_models.cstr import CSTR
from watertap.costing import WaterTAPCosting

from watertap.property_models.unit_specific.activated_sludge.asm1_properties import (
    ASM1ParameterBlock,
)
from watertap.property_models.unit_specific.activated_sludge.asm1_reactions import (
    ASM1ReactionParameterBlock,
)

from idaes.models.properties.examples.saponification_thermo import (
    SaponificationParameterBlock,
)
from idaes.models.properties.examples.saponification_reactions import (
    SaponificationReactionParameterBlock,
)
from watertap.core.solvers import get_solver
from idaes.core.initialization import (
    BlockTriangularizationInitializer,
    SingleControlVolumeUnitInitializer,
    InitializationStatus,
)

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
def build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = SaponificationParameterBlock()
    m.fs.reactions = SaponificationReactionParameterBlock(
        property_package=m.fs.properties
    )

    m.fs.unit = CSTR(
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

    iscale.set_scaling_factor(
        m.fs.unit.control_volume.properties_out[0.0].conc_mol_comp["H2O"], 1e-5
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.properties_out[0.0].pressure, 1e-6
    )

    iscale.calculate_scaling_factors(m.fs.unit)

    return m


def build_ASM1():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.props_ASM1 = ASM1ParameterBlock()
    m.fs.ASM1_rxn_props = ASM1ReactionParameterBlock(property_package=m.fs.props_ASM1)

    m.fs.unit = CSTR(
        property_package=m.fs.props_ASM1, reaction_package=m.fs.ASM1_rxn_props
    )

    m.fs.unit.inlet.flow_vol[0].fix(1.2199 * units.m**3 / units.s)
    m.fs.unit.inlet.alkalinity[0].fix(4.5102 * units.mole / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "S_I"].fix(0.061909 * units.kg / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "S_S"].fix(0.012366 * units.kg / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "X_I"].fix(1.4258 * units.kg / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "X_S"].fix(0.090508 * units.kg / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "X_BH"].fix(2.8404 * units.kg / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "X_BA"].fix(0.20512 * units.kg / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "X_P"].fix(0.58681 * units.kg / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "S_O"].fix(0.00036092 * units.kg / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "S_NO"].fix(0.012424 * units.kg / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "S_NH"].fix(0.0076936 * units.kg / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "S_ND"].fix(0.0019068 * units.kg / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "X_ND"].fix(0.0053166 * units.kg / units.m**3)

    m.fs.unit.inlet.temperature[0].fix(308.15 * units.K)
    m.fs.unit.inlet.pressure[0].fix(84790.0 * units.Pa)

    m.fs.unit.volume[0].fix(1000 * units.m**3)

    # Set scaling factors for badly scaled variables
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.properties_out[0.0].conc_mass_comp["S_O"], 1e5
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_generation[0.0, "Liq", "X_BA"], 1e5
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_generation[0.0, "Liq", "X_P"], 1e5
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_generation[0.0, "Liq", "S_O"], 1e5
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_generation[0.0, "Liq", "S_NH"], 1e5
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_generation[0.0, "Liq", "S_ND"], 1e5
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_generation[0.0, "Liq", "X_ND"], 1e5
    )

    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_extent[0.0, "R1"], 1e5
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_extent[0.0, "R3"], 1e5
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_extent[0.0, "R5"], 1e5
    )

    iscale.set_scaling_factor(
        m.fs.unit.control_volume.reactions[0.0].reaction_rate["R1"], 1e7
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.reactions[0.0].reaction_rate["R2"], 1e7
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.reactions[0.0].reaction_rate["R3"], 1e7
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.reactions[0.0].reaction_rate["R4"], 1e7
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.reactions[0.0].reaction_rate["R5"], 1e7
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.reactions[0.0].reaction_rate["R6"], 1e7
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.reactions[0.0].reaction_rate["R7"], 1e7
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.reactions[0.0].reaction_rate["R8"], 1e7
    )

    iscale.calculate_scaling_factors(m.fs.unit)

    return m


class TestCSTR(UnitTestHarness):
    def configure(self):
        m = build()

        self.unit_solutions[m.fs.unit.outlet.pressure[0]] = 101325.0
        self.unit_solutions[m.fs.unit.outlet.temperature[0]] = 304.09
        self.unit_solutions[m.fs.unit.outlet.conc_mol_comp[0, "EthylAcetate"]] = 20.32

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.inlet.flow_vol[0]
                * sum(
                    m.fs.unit.inlet.conc_mol_comp[0, j]
                    for j in m.fs.properties.component_list
                ),
                "out": m.fs.unit.outlet.flow_vol[0]
                * sum(
                    m.fs.unit.inlet.conc_mol_comp[0, j]
                    for j in m.fs.properties.component_list
                ),
            },
            "Check 2": {
                "in": (
                    m.fs.unit.inlet.flow_vol[0]
                    * m.fs.properties.dens_mol
                    * m.fs.properties.cp_mol
                    * (m.fs.unit.inlet.temperature[0] - m.fs.properties.temperature_ref)
                ),
                "out": m.fs.unit.outlet.flow_vol[0]
                * m.fs.properties.dens_mol
                * m.fs.properties.cp_mol
                * (m.fs.unit.outlet.temperature[0] - m.fs.properties.temperature_ref)
                - m.fs.unit.control_volume.heat_of_reaction[0],
            },
        }

        return m


class TestInitializers:
    @pytest.fixture
    def model(self):
        m = build()
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = SaponificationParameterBlock()
        m.fs.reactions = SaponificationReactionParameterBlock(
            property_package=m.fs.properties
        )

        m.fs.unit = CSTR(
            property_package=m.fs.properties,
            reaction_package=m.fs.reactions,
            has_equilibrium_reactions=False,
            has_heat_transfer=True,
            has_heat_of_reaction=True,
            has_pressure_change=True,
        )

        m.fs.unit.inlet.flow_vol[0].set_value(1.0e-03)
        m.fs.unit.inlet.conc_mol_comp[0, "H2O"].set_value(55388.0)
        m.fs.unit.inlet.conc_mol_comp[0, "NaOH"].set_value(100.0)
        m.fs.unit.inlet.conc_mol_comp[0, "EthylAcetate"].set_value(100.0)
        m.fs.unit.inlet.conc_mol_comp[0, "SodiumAcetate"].set_value(0.0)
        m.fs.unit.inlet.conc_mol_comp[0, "Ethanol"].set_value(0.0)

        m.fs.unit.inlet.temperature[0].set_value(303.15)
        m.fs.unit.inlet.pressure[0].set_value(101325.0)

        m.fs.unit.volume[0].fix(1.5e-03)
        m.fs.unit.heat_duty[0].fix(0)
        m.fs.unit.deltaP[0].fix(0)

        return m

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_general_hierarchical(self, model):
        initializer = SingleControlVolumeUnitInitializer()
        initializer.initialize(model.fs.unit)

        assert initializer.summary[model.fs.unit]["status"] == InitializationStatus.Ok

        assert value(model.fs.unit.outlet.flow_vol[0]) == pytest.approx(1e-3, rel=1e-5)
        assert value(model.fs.unit.outlet.conc_mol_comp[0, "H2O"]) == pytest.approx(
            55388, rel=1e-5
        )
        assert value(model.fs.unit.outlet.conc_mol_comp[0, "NaOH"]) == pytest.approx(
            20.31609, rel=1e-5
        )
        assert value(
            model.fs.unit.outlet.conc_mol_comp[0, "EthylAcetate"]
        ) == pytest.approx(20.31609, rel=1e-5)
        assert value(
            model.fs.unit.outlet.conc_mol_comp[0, "SodiumAcetate"]
        ) == pytest.approx(79.683910, rel=1e-5)
        assert value(model.fs.unit.outlet.conc_mol_comp[0, "Ethanol"]) == pytest.approx(
            79.683910, rel=1e-5
        )
        assert value(model.fs.unit.outlet.temperature[0]) == pytest.approx(
            304.0856, rel=1e-5
        )
        assert value(model.fs.unit.outlet.pressure[0]) == pytest.approx(
            101325, rel=1e-5
        )

        assert not model.fs.unit.inlet.flow_vol[0].fixed
        assert not model.fs.unit.inlet.conc_mol_comp[0, "H2O"].fixed
        assert not model.fs.unit.inlet.conc_mol_comp[0, "NaOH"].fixed
        assert not model.fs.unit.inlet.conc_mol_comp[0, "EthylAcetate"].fixed
        assert not model.fs.unit.inlet.conc_mol_comp[0, "SodiumAcetate"].fixed
        assert not model.fs.unit.inlet.conc_mol_comp[0, "Ethanol"].fixed

        assert not model.fs.unit.inlet.temperature[0].fixed
        assert not model.fs.unit.inlet.pressure[0].fixed

    @pytest.mark.component
    def test_block_triangularization(self, model):
        initializer = BlockTriangularizationInitializer(constraint_tolerance=2e-5)
        initializer.initialize(model.fs.unit)

        assert initializer.summary[model.fs.unit]["status"] == InitializationStatus.Ok

        assert value(model.fs.unit.outlet.flow_vol[0]) == pytest.approx(1e-3, rel=1e-5)
        assert value(model.fs.unit.outlet.conc_mol_comp[0, "H2O"]) == pytest.approx(
            55388, rel=1e-5
        )
        assert value(model.fs.unit.outlet.conc_mol_comp[0, "NaOH"]) == pytest.approx(
            20.31609, rel=1e-5
        )
        assert value(
            model.fs.unit.outlet.conc_mol_comp[0, "EthylAcetate"]
        ) == pytest.approx(20.31609, rel=1e-5)
        assert value(
            model.fs.unit.outlet.conc_mol_comp[0, "SodiumAcetate"]
        ) == pytest.approx(79.683910, rel=1e-5)
        assert value(model.fs.unit.outlet.conc_mol_comp[0, "Ethanol"]) == pytest.approx(
            79.683910, rel=1e-5
        )
        assert value(model.fs.unit.outlet.temperature[0]) == pytest.approx(
            304.0856, rel=1e-5
        )
        assert value(model.fs.unit.outlet.pressure[0]) == pytest.approx(
            101325, rel=1e-5
        )

        assert not model.fs.unit.inlet.flow_vol[0].fixed
        assert not model.fs.unit.inlet.conc_mol_comp[0, "H2O"].fixed
        assert not model.fs.unit.inlet.conc_mol_comp[0, "NaOH"].fixed
        assert not model.fs.unit.inlet.conc_mol_comp[0, "EthylAcetate"].fixed
        assert not model.fs.unit.inlet.conc_mol_comp[0, "SodiumAcetate"].fixed
        assert not model.fs.unit.inlet.conc_mol_comp[0, "Ethanol"].fixed

        assert not model.fs.unit.inlet.temperature[0].fixed
        assert not model.fs.unit.inlet.pressure[0].fixed


class TestCosting(UnitTestHarness):
    def configure(self):
        m = build_ASM1()
        m.fs.unit.initialize()

        solver.solve(m)

        # Add unit model costing
        m.fs.costing = WaterTAPCosting()

        m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(m.fs.unit.control_volume.properties_out[0].flow_vol)
        m.objective = Objective(expr=m.fs.costing.LCOW)

        iscale.set_scaling_factor(m.fs.unit.costing.capital_cost, 1e-7)

        iscale.calculate_scaling_factors(m.fs.unit)

        self.unit_solutions[m.fs.unit.costing.capital_cost] = 566989.10
        self.unit_solutions[m.fs.costing.LCOW] = 0.002127

        solver.solve(m, tee=True)

        assert pytest.approx(0.002127, rel=1e-3) == value(m.fs.costing.LCOW)

        component_list = [
            "S_I",
            "S_S",
            "X_I",
            "X_S",
            "X_BH",
            "X_BA",
            "X_P",
            "S_O",
            "S_NO",
            "S_NH",
            "S_ND",
            "X_ND",
        ]

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.inlet.flow_vol[0],
                "out": m.fs.unit.outlet.flow_vol[0],
            },
        }

        return m
