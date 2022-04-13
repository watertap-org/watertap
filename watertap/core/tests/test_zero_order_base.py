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
"""
Tests for general zero-order base class
"""
import pytest

from idaes.core import declare_process_block_class, FlowsheetBlock
from idaes.core.util.exceptions import ConfigurationError
from pyomo.environ import ConcreteModel, Var

from watertap.core import ZeroOrderBaseData
from watertap.core import WaterParameterBlock
import idaes.logger as idaeslog


@declare_process_block_class("DerivedZOBase")
class DerivedZOBaseData(ZeroOrderBaseData):
    def build(self):
        super().build()


@pytest.mark.unit
def test_private_attributes():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.params = WaterParameterBlock(default={"solute_list": ["A", "B", "C"]})

    m.fs.unit = DerivedZOBase(default={"property_package": m.fs.params})

    assert m.fs.unit._tech_type is None
    assert m.fs.unit._has_recovery_removal is False
    assert m.fs.unit._fixed_perf_vars == []
    assert m.fs.unit._initialize is None
    assert m.fs.unit._scaling is None
    assert m.fs.unit._get_Q is None
    assert m.fs.unit._stream_table_dict == {}
    assert m.fs.unit._perf_var_dict == {}


class TestZOBase:
    @pytest.fixture
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.params = WaterParameterBlock(default={"solute_list": ["A", "B", "C"]})

        m.fs.params.del_component(m.fs.params.phase_list)
        m.fs.params.del_component(m.fs.params.solvent_set)
        m.fs.params.del_component(m.fs.params.solute_set)
        m.fs.params.del_component(m.fs.params.component_list)

        return m

    @pytest.mark.unit
    def test_phase_list(self, model):
        model.fs.params.phase_list = ["foo"]

        with pytest.raises(
            ConfigurationError,
            match="fs.unit configured with invalid property "
            "package. Zero-order models only support property "
            "packages with a single phase named 'Liq'.",
        ):
            model.fs.unit = DerivedZOBase(default={"property_package": model.fs.params})

    @pytest.mark.unit
    def test_no_solvent_set(self, model):
        model.fs.params.phase_list = ["Liq"]

        with pytest.raises(
            ConfigurationError,
            match="fs.unit configured with invalid property "
            "package. Zero-order models only support property "
            "packages which include 'H2O' as the only Solvent.",
        ):
            model.fs.unit = DerivedZOBase(default={"property_package": model.fs.params})

    @pytest.mark.unit
    def test_invalid_solvent_set(self, model):
        model.fs.params.phase_list = ["Liq"]
        model.fs.params.solvent_set = ["foo"]

        with pytest.raises(
            ConfigurationError,
            match="fs.unit configured with invalid property "
            "package. Zero-order models only support property "
            "packages which include 'H2O' as the only Solvent.",
        ):
            model.fs.unit = DerivedZOBase(default={"property_package": model.fs.params})

    @pytest.mark.unit
    def test_no_solute_set(self, model):
        model.fs.params.phase_list = ["Liq"]
        model.fs.params.solvent_set = ["H2O"]

        with pytest.raises(
            ConfigurationError,
            match="fs.unit configured with invalid property "
            "package. Zero-order models require property "
            "packages to declare all dissolved species as "
            "Solutes.",
        ):
            model.fs.unit = DerivedZOBase(default={"property_package": model.fs.params})

    @pytest.mark.unit
    def test_non_solvent_or_solute(self, model):
        model.fs.params.phase_list = ["Liq"]
        model.fs.params.solvent_set = ["H2O"]
        model.fs.params.solute_set = ["A", "B", "C"]
        model.fs.params.component_list = ["H2O", "A", "B", "C", "foo"]

        with pytest.raises(
            ConfigurationError,
            match="fs.unit configured with invalid property "
            "package. Zero-order models only support `H2O` as "
            "a solvent and all other species as Solutes.",
        ):
            model.fs.unit = DerivedZOBase(default={"property_package": model.fs.params})

    @pytest.mark.unit
    def test_load_parameters_from_database(self, model):
        model.fs.params.phase_list = ["Liq"]
        model.fs.params.solvent_set = ["H2O"]
        model.fs.params.solute_set = ["A", "B", "C"]
        model.fs.params.component_list = ["H2O", "A", "B", "C"]

        model.fs.unit = DerivedZOBase(default={"property_package": model.fs.params})

        with pytest.raises(NotImplementedError):
            model.fs.unit.load_parameters_from_database()

    @pytest.mark.unit
    def test_set_param_from_data(self, model, caplog):
        caplog.set_level(idaeslog.DEBUG, logger="watertap")
        log = idaeslog.getLogger("idaes.watertap.core.zero_order_base")
        log.setLevel(idaeslog.DEBUG)

        model.fs.params.phase_list = ["Liq"]
        model.fs.params.solvent_set = ["H2O"]
        model.fs.params.solute_set = ["A", "B", "C"]
        model.fs.params.component_list = ["H2O", "A", "B", "C"]

        model.fs.unit = DerivedZOBase(default={"property_package": model.fs.params})

        model.fs.unit.recovery_vol = Var(model.fs.time)

        model.fs.unit.set_param_from_data(
            model.fs.unit.recovery_vol,
            {"recovery_vol": {"value": 0.42, "units": "m^3/m^3"}},
        )

        assert model.fs.unit.recovery_vol[0].value == 0.42
        assert model.fs.unit.recovery_vol[0].fixed

        assert "fs.unit.recovery_vol fixed to value 0.42 dimensionless" in caplog.text

    @pytest.mark.unit
    def test_set_param_from_data_no_entry(self, model):
        model.fs.params.phase_list = ["Liq"]
        model.fs.params.solvent_set = ["H2O"]
        model.fs.params.solute_set = ["A", "B", "C"]
        model.fs.params.component_list = ["H2O", "A", "B", "C"]

        model.fs.unit = DerivedZOBase(default={"property_package": model.fs.params})

        model.fs.unit.recovery_vol = Var(model.fs.time)

        with pytest.raises(
            KeyError,
            match="fs.unit - database provided does not "
            "contain an entry for recovery_vol for technology.",
        ):
            model.fs.unit.set_param_from_data(model.fs.unit.recovery_vol, {})

    @pytest.mark.unit
    def test_set_param_from_data_no_value(self, model):
        model.fs.params.phase_list = ["Liq"]
        model.fs.params.solvent_set = ["H2O"]
        model.fs.params.solute_set = ["A", "B", "C"]
        model.fs.params.component_list = ["H2O", "A", "B", "C"]

        model.fs.unit = DerivedZOBase(default={"property_package": model.fs.params})

        model.fs.unit.recovery_vol = Var(model.fs.time)

        with pytest.raises(
            KeyError,
            match="fs.unit - no value provided for recovery_vol"
            " \(index: None\) in database.",
        ):
            model.fs.unit.set_param_from_data(
                model.fs.unit.recovery_vol, {"recovery_vol": {}}
            )

    @pytest.mark.unit
    def test_set_param_from_data_no_units(self, model):
        model.fs.params.phase_list = ["Liq"]
        model.fs.params.solvent_set = ["H2O"]
        model.fs.params.solute_set = ["A", "B", "C"]
        model.fs.params.component_list = ["H2O", "A", "B", "C"]

        model.fs.unit = DerivedZOBase(default={"property_package": model.fs.params})

        model.fs.unit.recovery_vol = Var(model.fs.time)

        with pytest.raises(
            KeyError,
            match="fs.unit - no units provided for recovery_vol"
            " \(index: None\) in database.",
        ):
            model.fs.unit.set_param_from_data(
                model.fs.unit.recovery_vol, {"recovery_vol": {"value": 0.42}}
            )

    @pytest.mark.unit
    def test_set_param_from_data_indexed(self, model):
        model.fs.params.phase_list = ["Liq"]
        model.fs.params.solvent_set = ["H2O"]
        model.fs.params.solute_set = ["A", "B", "C"]
        model.fs.params.component_list = ["H2O", "A", "B", "C"]

        model.fs.unit = DerivedZOBase(default={"property_package": model.fs.params})

        model.fs.unit.removal_frac_mass_solute = Var(
            model.fs.time, model.fs.params.solute_set
        )

        model.fs.unit.set_param_from_data(
            model.fs.unit.removal_frac_mass_solute[0, "A"],
            {"removal_frac_mass_solute": {"A": {"value": 0.42, "units": "m^3/m^3"}}},
            index="A",
        )

        assert model.fs.unit.removal_frac_mass_solute[0, "A"].value == 0.42
        assert model.fs.unit.removal_frac_mass_solute[0, "A"].fixed

    @pytest.mark.unit
    def test_set_param_from_data_indexed_no_entry(self, model):
        model.fs.params.phase_list = ["Liq"]
        model.fs.params.solvent_set = ["H2O"]
        model.fs.params.solute_set = ["A", "B", "C"]
        model.fs.params.component_list = ["H2O", "A", "B", "C"]

        model.fs.unit = DerivedZOBase(default={"property_package": model.fs.params})

        model.fs.unit.removal_frac_mass_solute = Var(
            model.fs.time, model.fs.params.solute_set
        )

        with pytest.raises(
            KeyError,
            match="fs.unit - database provided does not "
            "contain an entry for removal_frac_mass_solute with "
            "index A for technology.",
        ):
            model.fs.unit.set_param_from_data(
                model.fs.unit.removal_frac_mass_solute[0, "A"],
                {"removal_frac_mass_solute": {}},
                index="A",
            )

    @pytest.mark.unit
    def test_set_param_from_data_indexed_default_not_removal(self, model):
        model.fs.params.phase_list = ["Liq"]
        model.fs.params.solvent_set = ["H2O"]
        model.fs.params.solute_set = ["A", "B", "C"]
        model.fs.params.component_list = ["H2O", "A", "B", "C"]

        model.fs.unit = DerivedZOBase(default={"property_package": model.fs.params})

        model.fs.unit.test = Var()

        model.fs.unit.removal_frac_mass_solute = Var(
            model.fs.time, model.fs.params.solute_set
        )

        with pytest.raises(
            KeyError,
            match="fs.unit - database provided does not "
            "contain an entry for test with "
            "index A for technology.",
        ):
            model.fs.unit.set_param_from_data(
                model.fs.unit.test, {"test": {}}, index="A", use_default_removal=True
            )

    @pytest.mark.unit
    def test_set_param_from_data_indexed_use_default(self, model):
        model.fs.params.phase_list = ["Liq"]
        model.fs.params.solvent_set = ["H2O"]
        model.fs.params.solute_set = ["A", "B", "C"]
        model.fs.params.component_list = ["H2O", "A", "B", "C"]

        model.fs.unit = DerivedZOBase(default={"property_package": model.fs.params})

        model.fs.unit.removal_frac_mass_solute = Var(
            model.fs.time, model.fs.params.solute_set
        )

        model.fs.unit.set_param_from_data(
            model.fs.unit.removal_frac_mass_solute[0, "A"],
            {
                "removal_frac_mass_solute": {"A": {"value": 0.42, "units": "m^3/m^3"}},
                "default_removal_frac_mass_solute": {"value": 0.70, "units": "kg/kg"},
            },
            index="D",
            use_default_removal=True,
        )

        assert model.fs.unit.removal_frac_mass_solute[0, "A"].value == 0.70
        assert model.fs.unit.removal_frac_mass_solute[0, "A"].fixed

    @pytest.mark.unit
    def test_set_param_from_data_indexed_use_default_undefined(self, model):
        model.fs.params.phase_list = ["Liq"]
        model.fs.params.solvent_set = ["H2O"]
        model.fs.params.solute_set = ["A", "B", "C"]
        model.fs.params.component_list = ["H2O", "A", "B", "C"]

        model.fs.unit = DerivedZOBase(default={"property_package": model.fs.params})

        model.fs.unit.removal_frac_mass_solute = Var(
            model.fs.time, model.fs.params.solute_set
        )

        with pytest.raises(
            KeyError,
            match="fs.unit - database provided does not "
            "contain an entry for removal_frac_mass_solute with "
            "index D for technology and no default removal was "
            "specified.",
        ):
            model.fs.unit.set_param_from_data(
                model.fs.unit.removal_frac_mass_solute[0, "A"],
                {
                    "removal_frac_mass_solute": {
                        "A": {"value": 0.42, "units": "m^3/m^3"}
                    }
                },
                index="D",
                use_default_removal=True,
            )
