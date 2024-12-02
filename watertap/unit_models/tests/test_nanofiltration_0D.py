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

import re
import pytest

from pyomo.environ import ConcreteModel, Constraint, Set, Suffix, TransformationFactory, units as pyunits, value, Var
from pyomo.network import Port

from idaes.core import FlowsheetBlock
from idaes.core.initialization import InitializationStatus
from idaes.core.scaling import set_scaling_factor
from idaes.core.solvers import get_solver
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.util.testing import PhysicalParameterTestBlock

from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
)
from watertap.unit_models.nanofiltration0D import (
    Nanofiltration0D,
    Nanofiltration0DInitializer,
    Nanofiltration0DScaler,
)


class TestBuild:
    @pytest.mark.unit
    def test_config(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = MCASParameterBlock(
            solute_list=["Na_+", "Cl_-"],
            charge={"Na_+": 1, "Cl_-": -1},
        )
        m.fs.unit = Nanofiltration0D(property_package=m.fs.properties)

        assert len(m.fs.unit.config) == 8
        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert m.fs.unit.config.property_package is m.fs.properties
        assert m.fs.unit.config.property_package_args == {}
        assert m.fs.unit.config.include_pressure_balance
        assert not m.fs.unit.config.has_retentate_pressure_drop
        assert m.fs.unit.config.include_temperature_equality
        assert m.fs.unit.config.electroneutrality_ion == "Cl_-"

    @pytest.mark.unit
    def test_default_initializer_and_scaler(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = MCASParameterBlock(
            solute_list=["Na_+", "Cl_-"],
            charge={"Na_+": 1, "Cl_-": -1},
        )
        m.fs.unit = Nanofiltration0D(property_package=m.fs.properties)

        assert m.fs.unit.default_initializer is Nanofiltration0DInitializer
        assert m.fs.unit.default_scaler is Nanofiltration0DScaler

    @pytest.mark.unit
    def test_multiphase(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = PhysicalParameterTestBlock()

        with pytest.raises(
            ConfigurationError,
            match="Nanofiltration0D model only supports single phase "
            "property packages.",
        ):
            m.fs.unit = Nanofiltration0D(property_package=m.fs.properties)

    @pytest.mark.unit
    def test_electroneutrality_invalid(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = MCASParameterBlock(
            solute_list=["Na_+", "Cl_-"],
            charge={"Na_+": 1, "Cl_-": -1},
        )

        with pytest.raises(
                ConfigurationError,
                match=re.escape("electroneutrality_ion (foo) "
                "is not a valid component in property package."),
        ):
            m.fs.unit = Nanofiltration0D(
                property_package=m.fs.properties,
                electroneutrality_ion="foo",
            )

    @pytest.mark.unit
    def test_electroneutrality_divalent(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = MCASParameterBlock(
            solute_list=["Na_+", "Cl_-"],
            charge={"Na_+": 1, "Cl_-": -2},
        )

        with pytest.raises(
                ConfigurationError,
                match=re.escape(
                    "electroneutrality_ion (Cl_-) must be a monovalent ion."
                ),
        ):
            m.fs.unit = Nanofiltration0D(
                property_package=m.fs.properties,
            )

    @pytest.mark.unit
    def test_invalid_pressure_drop(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = MCASParameterBlock(
            solute_list=["Na_+", "Cl_-"],
            charge={"Na_+": 1, "Cl_-": -1},
        )

        with pytest.raises(
                ConfigurationError,
                match="Inconsistent configuration arguments. has_retentate_pressure_drop=True "
                "requires that include_pressure_balance=True.",
        ):
            m.fs.unit = Nanofiltration0D(
                property_package=m.fs.properties,
                include_pressure_balance=False,
                has_retentate_pressure_drop=True,
            )

    @pytest.mark.unit
    def test_default_build(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = MCASParameterBlock(
            solute_list=["Ca_2+", "SO4_2-", "Mg_2+", "Na_+", "Cl_-"],
            charge={"Ca_2+": 2, "Mg_2+": 2, "SO4_2-": -2, "Na_+": 1, "Cl_-": -1},
        )
        m.fs.unit = Nanofiltration0D(property_package=m.fs.properties)

        assert hasattr(m.fs.unit, "properties_in")
        assert hasattr(m.fs.unit, "properties_retentate")
        assert hasattr(m.fs.unit, "properties_permeate")

        assert isinstance(m.fs.unit.inlet, Port)
        assert isinstance(m.fs.unit.retentate, Port)
        assert isinstance(m.fs.unit.permeate, Port)

        assert isinstance(m.fs.unit.split_species, Set)
        assert len(m.fs.unit.split_species) == 1
        assert "Na_+" in m.fs.unit.split_species

        assert isinstance(m.fs.unit.solvent_recovery, Var)
        assert value(m.fs.unit.solvent_recovery) == pytest.approx(0.8, rel=1e-8)

        assert isinstance(m.fs.unit.multivalent_recovery, Var)
        assert value(m.fs.unit.multivalent_recovery) == pytest.approx(1e-10, abs=1e-12)
        assert m.fs.unit.multivalent_recovery.fixed

        assert isinstance(m.fs.unit.solute_recovery, Var)
        assert len(m.fs.unit.solute_recovery) == 1
        assert value(m.fs.unit.solute_recovery["Na_+"]) == pytest.approx(0.9, rel=1e-8)

        assert not hasattr(m.fs.unit, "deltaP")

        assert isinstance(m.fs.unit.material_balances, Constraint)
        assert len(m.fs.unit.material_balances) == 6  # 1 time point, 6 components
        for (t, j), cd in m.fs.unit.material_balances.items():
            assert str(cd.expr) == str(
                m.fs.unit.properties_in[t].flow_mol_phase_comp["Liq", j]
                == m.fs.unit.properties_retentate[t].flow_mol_phase_comp["Liq", j]
                + m.fs.unit.properties_permeate[t].flow_mol_phase_comp["Liq", j]
            )

        assert isinstance(m.fs.unit.recovery_constraint, Constraint)
        assert len(m.fs.unit.recovery_constraint) == 5  # 1 time point, 5 components (Cl- is skipped)
        for (t, j), cd in m.fs.unit.recovery_constraint.items():
            assert j != "Cl_-"

            if j == "H2O":
                split = m.fs.unit.solvent_recovery
            elif j == "Na_+":
                split = m.fs.unit.solute_recovery[j]
            else:
                split = m.fs.unit.multivalent_recovery

            assert str(cd.expr) == str(
                split * m.fs.unit.properties_in[t].flow_mol_phase_comp["Liq", j]
                == m.fs.unit.properties_permeate[t].flow_mol_phase_comp["Liq", j]
            )

        assert isinstance(m.fs.unit.permeate_electronegativity, Constraint)
        assert len(m.fs.unit.permeate_electronegativity) == 1
        assert str(m.fs.unit.permeate_electronegativity[0].expr) == str(
            0 == sum(
                m.fs.properties.get_component(j).config.charge
                * m.fs.unit.properties_permeate[0].flow_mol_phase_comp["Liq", j]
                for j in m.fs.properties.ion_set
            )
        )

        assert isinstance(m.fs.unit.retentate_pressure_balance, Constraint)
        assert len(m.fs.unit.retentate_pressure_balance) == 1
        assert str(m.fs.unit.retentate_pressure_balance[0].expr) == str(
            m.fs.unit.properties_in[0].pressure
            == m.fs.unit.properties_retentate[0].pressure
        )

        assert isinstance(m.fs.unit.retentate_temperature_equality, Constraint)
        assert len(m.fs.unit.retentate_temperature_equality) == 1
        assert str(m.fs.unit.retentate_temperature_equality[0].expr) == str(
            m.fs.unit.properties_in[0].temperature
            == m.fs.unit.properties_retentate[0].temperature
        )

        assert isinstance(m.fs.unit.permeate_temperature_equality, Constraint)
        assert len(m.fs.unit.permeate_temperature_equality) == 1
        assert str(m.fs.unit.permeate_temperature_equality[0].expr) == str(
            m.fs.unit.properties_in[0].temperature
            == m.fs.unit.properties_permeate[0].temperature
        )

    @pytest.mark.unit
    def test_build_w_pressure_drop(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = MCASParameterBlock(
            solute_list=["Ca_2+", "SO4_2-", "Mg_2+", "Na_+", "Cl_-"],
            charge={"Ca_2+": 2, "Mg_2+": 2, "SO4_2-": -2, "Na_+": 1, "Cl_-": -1},
        )
        m.fs.unit = Nanofiltration0D(
            property_package=m.fs.properties,
            has_retentate_pressure_drop=True,
        )

        assert hasattr(m.fs.unit, "properties_in")
        assert hasattr(m.fs.unit, "properties_retentate")
        assert hasattr(m.fs.unit, "properties_permeate")

        assert isinstance(m.fs.unit.inlet, Port)
        assert isinstance(m.fs.unit.retentate, Port)
        assert isinstance(m.fs.unit.permeate, Port)

        assert isinstance(m.fs.unit.split_species, Set)
        assert len(m.fs.unit.split_species) == 1
        assert "Na_+" in m.fs.unit.split_species

        assert isinstance(m.fs.unit.solvent_recovery, Var)
        assert value(m.fs.unit.solvent_recovery) == pytest.approx(0.8, rel=1e-8)

        assert isinstance(m.fs.unit.multivalent_recovery, Var)
        assert value(m.fs.unit.multivalent_recovery) == pytest.approx(1e-10, abs=1e-12)
        assert m.fs.unit.multivalent_recovery.fixed

        assert isinstance(m.fs.unit.solute_recovery, Var)
        assert len(m.fs.unit.solute_recovery) == 1
        assert value(m.fs.unit.solute_recovery["Na_+"]) == pytest.approx(0.9, rel=1e-8)

        assert isinstance(m.fs.unit.deltaP, Var)

        assert isinstance(m.fs.unit.material_balances, Constraint)
        assert len(m.fs.unit.material_balances) == 6  # 1 time point, 6 components
        for (t, j), cd in m.fs.unit.material_balances.items():
            assert str(cd.expr) == str(
                m.fs.unit.properties_in[t].flow_mol_phase_comp["Liq", j]
                == m.fs.unit.properties_retentate[t].flow_mol_phase_comp["Liq", j]
                + m.fs.unit.properties_permeate[t].flow_mol_phase_comp["Liq", j]
            )

        assert isinstance(m.fs.unit.recovery_constraint, Constraint)
        assert len(m.fs.unit.recovery_constraint) == 5  # 1 time point, 5 components (Cl- is skipped)
        for (t, j), cd in m.fs.unit.recovery_constraint.items():
            assert j != "Cl_-"

            if j == "H2O":
                split = m.fs.unit.solvent_recovery
            elif j == "Na_+":
                split = m.fs.unit.solute_recovery[j]
            else:
                split = m.fs.unit.multivalent_recovery

            assert str(cd.expr) == str(
                split * m.fs.unit.properties_in[t].flow_mol_phase_comp["Liq", j]
                == m.fs.unit.properties_permeate[t].flow_mol_phase_comp["Liq", j]
            )

        assert isinstance(m.fs.unit.permeate_electronegativity, Constraint)
        assert len(m.fs.unit.permeate_electronegativity) == 1
        assert str(m.fs.unit.permeate_electronegativity[0].expr) == str(
            0 == sum(
                m.fs.properties.get_component(j).config.charge
                * m.fs.unit.properties_permeate[0].flow_mol_phase_comp["Liq", j]
                for j in m.fs.properties.ion_set
            )
        )

        assert isinstance(m.fs.unit.retentate_pressure_balance, Constraint)
        assert len(m.fs.unit.retentate_pressure_balance) == 1
        assert str(m.fs.unit.retentate_pressure_balance[0].expr) == str(
            m.fs.unit.properties_in[0].pressure
            == m.fs.unit.properties_retentate[0].pressure
            - m.fs.unit.deltaP[0]
        )

        assert isinstance(m.fs.unit.retentate_temperature_equality, Constraint)
        assert len(m.fs.unit.retentate_temperature_equality) == 1
        assert str(m.fs.unit.retentate_temperature_equality[0].expr) == str(
            m.fs.unit.properties_in[0].temperature
            == m.fs.unit.properties_retentate[0].temperature
        )

        assert isinstance(m.fs.unit.permeate_temperature_equality, Constraint)
        assert len(m.fs.unit.permeate_temperature_equality) == 1
        assert str(m.fs.unit.permeate_temperature_equality[0].expr) == str(
            m.fs.unit.properties_in[0].temperature
            == m.fs.unit.properties_permeate[0].temperature
        )

    @pytest.mark.unit
    def test_build_no_pt_electroneutrality(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = MCASParameterBlock(
            solute_list=["Ca_2+", "SO4_2-", "Mg_2+", "Na_+", "Cl_-"],
            charge={"Ca_2+": 2, "Mg_2+": 2, "SO4_2-": -2, "Na_+": 1, "Cl_-": -1},
        )
        m.fs.unit = Nanofiltration0D(
            property_package=m.fs.properties,
            electroneutrality_ion=None,
            include_pressure_balance=False,
            include_temperature_equality=False,
        )

        assert hasattr(m.fs.unit, "properties_in")
        assert hasattr(m.fs.unit, "properties_retentate")
        assert hasattr(m.fs.unit, "properties_permeate")

        assert isinstance(m.fs.unit.inlet, Port)
        assert isinstance(m.fs.unit.retentate, Port)
        assert isinstance(m.fs.unit.permeate, Port)

        assert isinstance(m.fs.unit.split_species, Set)
        assert len(m.fs.unit.split_species) == 2
        for j in m.fs.unit.split_species:
            assert j in ["Na_+", "Cl_-"]

        assert isinstance(m.fs.unit.solvent_recovery, Var)
        assert value(m.fs.unit.solvent_recovery) == pytest.approx(0.8, rel=1e-8)

        assert isinstance(m.fs.unit.multivalent_recovery, Var)
        assert value(m.fs.unit.multivalent_recovery) == pytest.approx(1e-10, abs=1e-12)
        assert m.fs.unit.multivalent_recovery.fixed

        assert isinstance(m.fs.unit.solute_recovery, Var)
        assert len(m.fs.unit.solute_recovery) == 2
        assert value(m.fs.unit.solute_recovery["Na_+"]) == pytest.approx(0.9, rel=1e-8)
        assert value(m.fs.unit.solute_recovery["Cl_-"]) == pytest.approx(0.9, rel=1e-8)

        assert not hasattr(m.fs.unit, "deltaP")

        assert isinstance(m.fs.unit.material_balances, Constraint)
        assert len(m.fs.unit.material_balances) == 6  # 1 time point, 6 components
        for (t, j), cd in m.fs.unit.material_balances.items():
            assert str(cd.expr) == str(
                m.fs.unit.properties_in[t].flow_mol_phase_comp["Liq", j]
                == m.fs.unit.properties_retentate[t].flow_mol_phase_comp["Liq", j]
                + m.fs.unit.properties_permeate[t].flow_mol_phase_comp["Liq", j]
            )

        assert isinstance(m.fs.unit.recovery_constraint, Constraint)
        assert len(m.fs.unit.recovery_constraint) == 6  # 1 time point, 6 components
        for (t, j), cd in m.fs.unit.recovery_constraint.items():
            if j == "H2O":
                split = m.fs.unit.solvent_recovery
            elif j in ["Na_+", "Cl_-"]:
                split = m.fs.unit.solute_recovery[j]
            else:
                split = m.fs.unit.multivalent_recovery

            assert str(cd.expr) == str(
                split * m.fs.unit.properties_in[t].flow_mol_phase_comp["Liq", j]
                == m.fs.unit.properties_permeate[t].flow_mol_phase_comp["Liq", j]
            )

        assert not hasattr(m.fs.unit, "permeate_electronegativity")
        assert not hasattr(m.fs.unit, "retentate_pressure_balance")
        assert not hasattr(m.fs.unit, "retentate_temperature_equality")
        assert not hasattr(m.fs.unit, "permeate_temperature_equality")


@pytest.fixture
def mcas_model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(
        solute_list=["Ca_2+", "SO4_2-", "Mg_2+", "Na_+", "Cl_-"],
        diffusivity_data={
            ("Liq", "Ca_2+"): 9.2e-10,
            ("Liq", "SO4_2-"): 1.06e-09,
            ("Liq", "Mg_2+"): 7.06e-10,
            ("Liq", "Na_+"): 1.33e-09,
            ("Liq", "Cl_-"): 2.03e-09,
        },
        mw_data={
            "H2O": 0.018,
            "Ca_2+": 0.04,
            "Mg_2+": 0.024,
            "SO4_2-": 0.096,
            "Na_+": 0.023,
            "Cl_-": 0.035,
        },
        stokes_radius_data={
            "Ca_2+": 3.09e-10,
            "Mg_2+": 3.47e-10,
            "SO4_2-": 2.3e-10,
            "Cl_-": 1.21e-10,
            "Na_+": 1.84e-10,
        },
        charge={"Ca_2+": 2, "Mg_2+": 2, "SO4_2-": -2, "Na_+": 1, "Cl_-": -1},
    )

    m.fs.unit = Nanofiltration0D(
        property_package=m.fs.properties,
        has_retentate_pressure_drop=True,
    )

    return m


class TestNanofiltration0DInitializer:
    @pytest.mark.component
    def test_initialize(self, mcas_model):
        # Fix DoF
        mcas_model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix(53.6036)
        mcas_model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Ca_2+"].fix(0.00955)
        mcas_model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Mg_2+"].fix(0.5808)
        mcas_model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "SO4_2-"].fix(0.02225)
        mcas_model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(1.61977)
        mcas_model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(0.48357)
        mcas_model.fs.unit.inlet.temperature[0].fix(298.15)
        mcas_model.fs.unit.inlet.pressure[0].fix(4e5)

        mcas_model.fs.unit.solvent_recovery.fix()
        mcas_model.fs.unit.solute_recovery.fix()
        mcas_model.fs.unit.permeate.pressure[0].fix(1e5)
        mcas_model.fs.unit.deltaP.fix(-1e4)

        # Initialize
        initializer = Nanofiltration0DInitializer()

        initializer.initialize(mcas_model.fs.unit)

        assert initializer.summary[mcas_model.fs.unit]["status"] == InitializationStatus.Ok

        for c in mcas_model.component_data_objects(Constraint, descend_into=True):
            print(c.name)
            assert value(c.body) <= 1e-5


class TestNanofiltration0DScaler:
    @pytest.mark.unit
    def test_variable_scaling_routine(self, mcas_model):
        scaler = Nanofiltration0DScaler()

        scaler.variable_scaling_routine(mcas_model.fs.unit)

        # Check that submodels have scaling factors
        assert isinstance(mcas_model.fs.unit.properties_in[0].scaling_factor, Suffix)
        assert len(mcas_model.fs.unit.properties_in[0].scaling_factor) > 1

        assert isinstance(mcas_model.fs.unit.properties_retentate[0].scaling_factor, Suffix)
        assert len(mcas_model.fs.unit.properties_retentate[0].scaling_factor) > 1

        assert isinstance(mcas_model.fs.unit.properties_permeate[0].scaling_factor, Suffix)
        assert len(mcas_model.fs.unit.properties_permeate[0].scaling_factor) > 1

        # Check scaling factors for unit variables
        assert isinstance(mcas_model.fs.unit.scaling_factor, Suffix)
        assert len(mcas_model.fs.unit.scaling_factor) == 3
        assert mcas_model.fs.unit.scaling_factor[mcas_model.fs.unit.multivalent_recovery] == 1e2
        assert mcas_model.fs.unit.scaling_factor[mcas_model.fs.unit.solvent_recovery] == 10
        assert mcas_model.fs.unit.scaling_factor[mcas_model.fs.unit.solute_recovery["Na_+"]] == 10

    @pytest.mark.unit
    def test_constraint_scaling_routine(self, mcas_model):
        scaler = Nanofiltration0DScaler()

        scaler.constraint_scaling_routine(mcas_model.fs.unit)

        # Check that submodels have scaling factors
        # No constraint scaling factors applied ot sub-models, so no check for length of Suffix
        assert isinstance(mcas_model.fs.unit.properties_in[0].scaling_factor, Suffix)
        assert isinstance(mcas_model.fs.unit.properties_retentate[0].scaling_factor, Suffix)
        assert isinstance(mcas_model.fs.unit.properties_permeate[0].scaling_factor, Suffix)

        # Check scaling factors for unit variables
        assert isinstance(mcas_model.fs.unit.scaling_factor, Suffix)
        assert len(mcas_model.fs.unit.scaling_factor) == 15

        assert mcas_model.fs.unit.scaling_factor[mcas_model.fs.unit.material_balances[0.0, "Ca_2+"]] == 10.0
        assert mcas_model.fs.unit.scaling_factor[mcas_model.fs.unit.material_balances[0.0, "Cl_-"]] == 10.0
        assert mcas_model.fs.unit.scaling_factor[mcas_model.fs.unit.material_balances[0.0, "H2O"]] == 10.0
        assert mcas_model.fs.unit.scaling_factor[mcas_model.fs.unit.material_balances[0.0, "Mg_2+"]] == 10.0
        assert mcas_model.fs.unit.scaling_factor[mcas_model.fs.unit.material_balances[0.0, "Na_+"]] == 10.0
        assert mcas_model.fs.unit.scaling_factor[mcas_model.fs.unit.material_balances[0.0, "SO4_2-"]] == 10.0
        assert mcas_model.fs.unit.scaling_factor[mcas_model.fs.unit.permeate_electronegativity[0.0]] == 5.0
        assert mcas_model.fs.unit.scaling_factor[mcas_model.fs.unit.permeate_temperature_equality[0.0]] == pytest.approx(3.354016e-3, rel=1e-5)
        assert mcas_model.fs.unit.scaling_factor[mcas_model.fs.unit.recovery_constraint[0.0, "Ca_2+"]] == 10.0
        assert mcas_model.fs.unit.scaling_factor[mcas_model.fs.unit.recovery_constraint[0.0, "H2O"]] == 10.0
        assert mcas_model.fs.unit.scaling_factor[mcas_model.fs.unit.recovery_constraint[0.0, "Mg_2+"]] == 10.0
        assert mcas_model.fs.unit.scaling_factor[mcas_model.fs.unit.recovery_constraint[0.0, "Na_+"]] == 10.0
        assert mcas_model.fs.unit.scaling_factor[mcas_model.fs.unit.recovery_constraint[0.0, "SO4_2-"]] == 10.0
        assert mcas_model.fs.unit.scaling_factor[mcas_model.fs.unit.retentate_pressure_balance[0.0]] == pytest.approx(9.869233e-06, rel=1e-5)
        assert mcas_model.fs.unit.scaling_factor[mcas_model.fs.unit.retentate_temperature_equality[0.0]] == pytest.approx(3.354016e-3, rel=1e-5)


# class TestMCAS:
#     @pytest.mark.component
#     def test_nf0d(self):
#         m = ConcreteModel()
#         m.fs = FlowsheetBlock(dynamic=False)
#         m.fs.properties = MCASParameterBlock(
#             solute_list=["Ca_2+", "SO4_2-", "Mg_2+", "Na_+", "Cl_-"],
#             diffusivity_data={
#                 ("Liq", "Ca_2+"): 9.2e-10,
#                 ("Liq", "SO4_2-"): 1.06e-09,
#                 ("Liq", "Mg_2+"): 7.06e-10,
#                 ("Liq", "Na_+"): 1.33e-09,
#                 ("Liq", "Cl_-"): 2.03e-09,
#             },
#             mw_data={
#                 "H2O": 0.018,
#                 "Ca_2+": 0.04,
#                 "Mg_2+": 0.024,
#                 "SO4_2-": 0.096,
#                 "Na_+": 0.023,
#                 "Cl_-": 0.035,
#             },
#             stokes_radius_data={
#                 "Ca_2+": 3.09e-10,
#                 "Mg_2+": 3.47e-10,
#                 "SO4_2-": 2.3e-10,
#                 "Cl_-": 1.21e-10,
#                 "Na_+": 1.84e-10,
#             },
#             charge={"Ca_2+": 2, "Mg_2+": 2, "SO4_2-": -2, "Na_+": 1, "Cl_-": -1},
#         )
#
#         m.fs.unit = Nanofiltration0D(property_package=m.fs.properties)
#
#         # Fix other inlet state variables
#         m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix(53.6036)
#         m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Ca_2+"].fix(0.00955)
#         m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Mg_2+"].fix(0.5808)
#         m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "SO4_2-"].fix(0.02225)
#         m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(1.61977)
#         m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(0.48357)
#         m.fs.unit.inlet.temperature[0].fix(298.15)
#         m.fs.unit.inlet.pressure[0].fix(4e5)
#
#         m.fs.unit.solvent_recovery.fix()
#         m.fs.unit.solute_recovery.fix()
#         m.fs.unit.permeate.pressure[0].fix(1e5)
#
#         # Scale model
#         set_scaling_factor(m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H2O"], 1e-1)
#         set_scaling_factor(m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Ca_2+"], 1e3)
#         set_scaling_factor(m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Mg_2+"], 10)
#         set_scaling_factor(m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "SO4_2-"], 1e2)
#         set_scaling_factor(m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Cl_-"], 1)
#         set_scaling_factor(m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Na_+"], 10)
#
#         scaler = m.fs.unit.default_scaler()
#         scaler.scale_model(m.fs.unit)
#
#         initializer = m.fs.unit.default_initializer()
#         try:
#             initializer.initialize(m.fs.unit)
#         except:
#             pass
#
#         solver = get_solver("ipopt_v2", writer_config={"scale_model": True, "linear_presolve": True})
#         solver.solve(m)
#
#         sm = TransformationFactory("core.scale_model").create_using(m, rename=False)
#
#         dt = DiagnosticsToolbox(sm.fs.unit)
#         dt.report_structural_issues()
#         dt.report_numerical_issues()
#
#         m.fs.unit.inlet.display()
#         m.fs.unit.retentate.display()
#         m.fs.unit.permeate.display()
#
#         assert False