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
Tests for general zero-order property package
"""
import pytest
import os

from idaes.core import declare_process_block_class, FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale
from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    Param,
    Set,
    value,
    Var,
)
from pyomo.network import Port
from pyomo.util.check_units import assert_units_consistent


from watertap.core import (
    Database,
    WaterParameterBlock,
    WaterStateBlock,
    ZeroOrderBaseData,
)
from watertap.core.zero_order_sido_reactive import (
    build_sido_reactive,
    initialize_sidor,
    calculate_scaling_factors_sidor,
    _get_Q_sidor,
)

solver = get_solver()

local_path = os.path.dirname(os.path.abspath(__file__))


@declare_process_block_class("DerivedSIDOR")
class DerivedSIDORData(ZeroOrderBaseData):
    def build(self):
        super().build()

        self._tech_type = "test_sidor_data"

        build_sido_reactive(self)


class TestSIDOR:
    @pytest.fixture(scope="module")
    def model(self):
        m = ConcreteModel()

        m.db = Database(dbpath=local_path)

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.water_props = WaterParameterBlock(solute_list=["A", "B", "C"])

        m.fs.unit = DerivedSIDOR(property_package=m.fs.water_props, database=m.db)

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(1000)
        m.fs.unit.inlet.flow_mass_comp[0, "A"].fix(10)
        m.fs.unit.inlet.flow_mass_comp[0, "B"].fix(20)
        m.fs.unit.inlet.flow_mass_comp[0, "C"].fix(30)

        m.fs.unit.load_parameters_from_database(use_default_removal=True)

        return m

    @pytest.mark.unit
    def test_private_attributes(self, model):
        assert model.fs.unit._has_recovery_removal is True
        assert model.fs.unit._fixed_perf_vars == []
        assert model.fs.unit._initialize is initialize_sidor
        assert model.fs.unit._scaling is calculate_scaling_factors_sidor
        assert model.fs.unit._get_Q is _get_Q_sidor
        assert model.fs.unit._stream_table_dict == {
            "Inlet": model.fs.unit.inlet,
            "Treated": model.fs.unit.treated,
            "Byproduct": model.fs.unit.byproduct,
        }
        assert model.fs.unit._perf_var_dict == {
            "Water Recovery": model.fs.unit.recovery_frac_mass_H2O,
            "Solute Removal": model.fs.unit.removal_frac_mass_comp,
            "Reaction Extent": model.fs.unit.extent_of_reaction,
        }

    @pytest.mark.unit
    def test_build(self, model):
        assert isinstance(model.fs.unit.properties_in, WaterStateBlock)
        assert isinstance(model.fs.unit.properties_treated, WaterStateBlock)
        assert isinstance(model.fs.unit.properties_byproduct, WaterStateBlock)

        assert isinstance(model.fs.unit.inlet, Port)
        assert isinstance(model.fs.unit.treated, Port)
        assert isinstance(model.fs.unit.byproduct, Port)

        assert isinstance(model.fs.unit.recovery_frac_mass_H2O, Var)
        assert len(model.fs.unit.recovery_frac_mass_H2O) == 1
        assert isinstance(model.fs.unit.removal_frac_mass_comp, Var)
        assert len(model.fs.unit.removal_frac_mass_comp) == 3

        assert isinstance(model.fs.unit.water_recovery_equation, Constraint)
        assert len(model.fs.unit.water_recovery_equation) == 1
        assert isinstance(model.fs.unit.water_balance, Constraint)
        assert len(model.fs.unit.water_balance) == 1
        assert isinstance(model.fs.unit.solute_removal_equation, Constraint)
        assert len(model.fs.unit.solute_removal_equation) == 3
        assert isinstance(model.fs.unit.solute_treated_equation, Constraint)
        assert len(model.fs.unit.solute_treated_equation) == 3

        assert isinstance(model.fs.unit.reaction_set, Set)
        assert isinstance(model.fs.unit.generation_ratio, Var)

        assert isinstance(model.fs.unit.reaction_conversion, Var)
        assert isinstance(model.fs.unit.extent_of_reaction, Var)

        assert isinstance(model.fs.unit.reaction_extent_equation, Constraint)

        for r in model.fs.unit.reaction_set:
            assert r in ["Rxn1", "Rxn2"]
            assert (0, r) in model.fs.unit.reaction_conversion
            assert (0, r) in model.fs.unit.extent_of_reaction
            assert (0, r) in model.fs.unit.reaction_extent_equation
            for j in model.fs.water_props.component_list:
                assert (r, j) in model.fs.unit.generation_ratio

    @pytest.mark.unit
    def test_key_components(self, model):
        assert str(model.fs.unit.properties_in[0].flow_mass_comp["A"]) in str(
            model.fs.unit.reaction_extent_equation[0, "Rxn1"].body
        )
        assert str(model.fs.unit.properties_in[0].flow_mass_comp["B"]) in str(
            model.fs.unit.reaction_extent_equation[0, "Rxn2"].body
        )

    @pytest.mark.unit
    def test_loading_data(self, model):
        # RXn 1 is all conversion ratios
        # Rxn 2 is all stocihiometry (incl. water)
        cfactor = {
            "Rxn1": {"A": -1, "B": 1, "C": 0, "H2O": 0},
            "Rxn2": {"A": 0, "B": -1, "C": 22 * 2 / 20, "H2O": -1 * 18 / 20},
        }

        for (r, j), p in model.fs.unit.generation_ratio.items():
            assert value(p) == cfactor[r][j]

        assert model.fs.unit.reaction_conversion[0, "Rxn1"].value == 0.8
        assert model.fs.unit.reaction_conversion[0, "Rxn2"].value == 0.1

        assert model.fs.unit.recovery_frac_mass_H2O[0].value == 0.85
        assert model.fs.unit.removal_frac_mass_comp[0, "A"].value == 0.5
        assert model.fs.unit.removal_frac_mass_comp[0, "B"].value == 0.4
        assert model.fs.unit.removal_frac_mass_comp[0, "C"].value == 0.0

    @pytest.mark.unit
    def test_degrees_of_freedom(self, model):
        assert degrees_of_freedom(model) == 0

    @pytest.mark.component
    def test_unit_consistency(self, model):
        assert_units_consistent(model)

    @pytest.mark.component
    def test_scaling(self, model):
        iscale.calculate_scaling_factors(model)

        assert (
            iscale.get_constraint_transform_applied_scaling_factor(
                model.fs.unit.water_recovery_equation[0]
            )
            == 1e5
        )
        assert (
            iscale.get_constraint_transform_applied_scaling_factor(
                model.fs.unit.water_balance[0]
            )
            == 1e5
        )
        assert (
            iscale.get_constraint_transform_applied_scaling_factor(
                model.fs.unit.solute_removal_equation[0, "A"]
            )
            == 1e5
        )
        assert (
            iscale.get_constraint_transform_applied_scaling_factor(
                model.fs.unit.solute_removal_equation[0, "B"]
            )
            == 1e5
        )
        assert (
            iscale.get_constraint_transform_applied_scaling_factor(
                model.fs.unit.solute_removal_equation[0, "C"]
            )
            == 1e5
        )
        assert (
            iscale.get_constraint_transform_applied_scaling_factor(
                model.fs.unit.solute_treated_equation[0, "A"]
            )
            == 1e5
        )
        assert (
            iscale.get_constraint_transform_applied_scaling_factor(
                model.fs.unit.solute_treated_equation[0, "B"]
            )
            == 1e5
        )
        assert (
            iscale.get_constraint_transform_applied_scaling_factor(
                model.fs.unit.solute_treated_equation[0, "C"]
            )
            == 1e5
        )

        assert (
            iscale.get_constraint_transform_applied_scaling_factor(
                model.fs.unit.reaction_extent_equation[0, "Rxn1"]
            )
            == 1e5
        )
        assert (
            iscale.get_constraint_transform_applied_scaling_factor(
                model.fs.unit.reaction_extent_equation[0, "Rxn2"]
            )
            == 1e5
        )

    @pytest.mark.component
    def test_initialization(self, model):
        initialization_tester(model)

    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, model):
        model.fs.unit.treated.display()
        model.fs.unit.byproduct.display()
        assert pytest.approx((1000 - 0.1 * 20 * 18 / 20) * 0.85, rel=1e-5) == value(
            model.fs.unit.treated.flow_mass_comp[0, "H2O"]
        )
        assert pytest.approx((1000 - 0.1 * 20 * 18 / 20) * 0.15, rel=1e-5) == value(
            model.fs.unit.byproduct.flow_mass_comp[0, "H2O"]
        )

        assert pytest.approx((10 * 0.2) * (1 - 0.5), rel=1e-5) == value(
            model.fs.unit.treated.flow_mass_comp[0, "A"]
        )
        assert pytest.approx((10 * 0.2) * 0.5, rel=1e-5) == value(
            model.fs.unit.byproduct.flow_mass_comp[0, "A"]
        )
        assert pytest.approx((20 + 0.8 * 10 - 0.1 * 20) * (1 - 0.4), rel=1e-5) == value(
            model.fs.unit.treated.flow_mass_comp[0, "B"]
        )
        assert pytest.approx((20 + 0.8 * 10 - 0.1 * 20) * 0.4, rel=1e-5) == value(
            model.fs.unit.byproduct.flow_mass_comp[0, "B"]
        )
        assert pytest.approx(30 + 0.1 * 20 * 2 * 22 / 20, rel=1e-5) == value(
            model.fs.unit.treated.flow_mass_comp[0, "C"]
        )
        assert pytest.approx(0, rel=1e-5) == value(
            model.fs.unit.byproduct.flow_mass_comp[0, "C"]
        )

    # Conservation would be similar to above, as need to calculate generation

    @pytest.mark.component
    def test_report(self, model):
        model.fs.unit.report()


class TestSIDORErrors:
    @pytest.mark.unit
    def test_reaction_list(self):
        m = ConcreteModel()

        m.db = Database(dbpath=local_path)

        # Load data from YAML so that we can modify it
        m.db._get_technology("test_sidor_data")
        m.db._cached_files["test_sidor_data"]["default"] = {}

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.water_props = WaterParameterBlock(solute_list=["A", "B", "C"])

        with pytest.raises(
            KeyError,
            match="fs.unit - database provided does not contain a list of "
            "reactions for this technology.",
        ):
            m.fs.unit = DerivedSIDOR(property_package=m.fs.water_props, database=m.db)

    @pytest.mark.unit
    def test_missing_conversion(self):
        m = ConcreteModel()

        m.db = Database(dbpath=local_path)

        # Load data from YAML so that we can modify it
        m.db._get_technology("test_sidor_data")
        m.db._cached_files["test_sidor_data"]["default"]["reactions"]["Rxn3"] = {}

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.water_props = WaterParameterBlock(solute_list=["A", "B", "C"])

        with pytest.raises(
            KeyError,
            match="fs.unit - database provided does not "
            "contain an entry for conversion for reaction Rxn3.",
        ):
            m.fs.unit = DerivedSIDOR(property_package=m.fs.water_props, database=m.db)

    @pytest.mark.unit
    def test_missing_key_reactant(self):
        m = ConcreteModel()

        m.db = Database(dbpath=local_path)

        # Load data from YAML so that we can modify it
        m.db._get_technology("test_sidor_data")
        m.db._cached_files["test_sidor_data"]["default"]["reactions"]["Rxn3"] = {}
        R3 = m.db._cached_files["test_sidor_data"]["default"]["reactions"]["Rxn3"]
        R3["conversion"] = 0.5

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.water_props = WaterParameterBlock(solute_list=["A", "B", "C"])

        with pytest.raises(
            KeyError,
            match="fs.unit - database provided does not "
            "contain an entry for key_reactant for reaction Rxn3.",
        ):
            m.fs.unit = DerivedSIDOR(property_package=m.fs.water_props, database=m.db)

    @pytest.mark.unit
    def test_invlaid_key_reactant(self):
        m = ConcreteModel()

        m.db = Database(dbpath=local_path)

        # Load data from YAML so that we can modify it
        m.db._get_technology("test_sidor_data")
        m.db._cached_files["test_sidor_data"]["default"]["reactions"]["Rxn3"] = {}
        R3 = m.db._cached_files["test_sidor_data"]["default"]["reactions"]["Rxn3"]
        R3["conversion"] = 0.5
        R3["key_reactant"] = "foo"

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.water_props = WaterParameterBlock(solute_list=["A", "B", "C"])

        with pytest.raises(
            ValueError,
            match="fs.unit - key_reactant foo for reaction Rxn3 "
            "is not in the component list used by the assigned property "
            "package.",
        ):
            m.fs.unit = DerivedSIDOR(property_package=m.fs.water_props, database=m.db)

    @pytest.mark.unit
    def test_missing_stoichiometry(self):
        m = ConcreteModel()

        m.db = Database(dbpath=local_path)

        # Load data from YAML so that we can modify it
        m.db._get_technology("test_sidor_data")
        m.db._cached_files["test_sidor_data"]["default"]["reactions"]["Rxn3"] = {}
        R3 = m.db._cached_files["test_sidor_data"]["default"]["reactions"]["Rxn3"]
        R3["conversion"] = 0.5
        R3["key_reactant"] = "A"

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.water_props = WaterParameterBlock(solute_list=["A", "B", "C"])

        with pytest.raises(
            KeyError,
            match="fs.unit - database provided does not "
            "contain an entry for stoichiometry for reaction Rxn3.",
        ):
            m.fs.unit = DerivedSIDOR(property_package=m.fs.water_props, database=m.db)

    @pytest.mark.unit
    def test_ratio_and_order(self):
        m = ConcreteModel()

        m.db = Database(dbpath=local_path)

        # Load data from YAML so that we can modify it
        m.db._get_technology("test_sidor_data")
        m.db._cached_files["test_sidor_data"]["default"]["reactions"]["Rxn3"] = {}
        R3 = m.db._cached_files["test_sidor_data"]["default"]["reactions"]["Rxn3"]
        R3["conversion"] = 0.5
        R3["key_reactant"] = "A"
        R3["stoichiometry"] = {"A": {"order": 1, "conversion_ratio": 1}}

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.water_props = WaterParameterBlock(solute_list=["A", "B", "C"])

        with pytest.raises(
            RuntimeError,
            match="fs.unit - database provides entries for both "
            "conversion_ratio and reaction order in reaction Rxn3. "
            "Please provide only one or the other.",
        ):
            m.fs.unit = DerivedSIDOR(property_package=m.fs.water_props, database=m.db)

    @pytest.mark.unit
    def test_no_ratio_or_order(self):
        m = ConcreteModel()

        m.db = Database(dbpath=local_path)

        # Load data from YAML so that we can modify it
        m.db._get_technology("test_sidor_data")
        m.db._cached_files["test_sidor_data"]["default"]["reactions"]["Rxn3"] = {}
        R3 = m.db._cached_files["test_sidor_data"]["default"]["reactions"]["Rxn3"]
        R3["conversion"] = 0.5
        R3["key_reactant"] = "A"
        R3["stoichiometry"] = {"A": {}}

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.water_props = WaterParameterBlock(solute_list=["A", "B", "C"])

        with pytest.raises(
            RuntimeError,
            match="fs.unit - database provided does not "
            "contain any information for conversion_ratio or reaction "
            "order w.r.t. species A in reaction Rxn3.",
        ):
            m.fs.unit = DerivedSIDOR(property_package=m.fs.water_props, database=m.db)

    @pytest.mark.unit
    def test_order_no_mw(self):
        m = ConcreteModel()

        m.db = Database(dbpath=local_path)

        # Load data from YAML so that we can modify it
        m.db._get_technology("test_sidor_data")
        m.db._cached_files["test_sidor_data"]["default"]["reactions"]["Rxn3"] = {}
        R3 = m.db._cached_files["test_sidor_data"]["default"]["reactions"]["Rxn3"]
        R3["conversion"] = 0.5
        R3["key_reactant"] = "B"
        R3["stoichiometry"] = {
            "A": {"order": 1},
            "B": {"order": 1, "molecular_weight": 1},
        }

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.water_props = WaterParameterBlock(solute_list=["A", "B", "C"])

        with pytest.raises(
            KeyError,
            match="fs.unit - database provided does not "
            "contain an entry for molecular_weight w.r.t. "
            "species A in reaction Rxn3.",
        ):
            m.fs.unit = DerivedSIDOR(property_package=m.fs.water_props, database=m.db)

    @pytest.mark.unit
    def test_order_no_key_order(self):
        m = ConcreteModel()

        m.db = Database(dbpath=local_path)

        # Load data from YAML so that we can modify it
        m.db._get_technology("test_sidor_data")
        m.db._cached_files["test_sidor_data"]["default"]["reactions"]["Rxn3"] = {}
        R3 = m.db._cached_files["test_sidor_data"]["default"]["reactions"]["Rxn3"]
        R3["conversion"] = 0.5
        R3["key_reactant"] = "B"
        R3["stoichiometry"] = {"A": {"order": 1, "molecular_weight": 1}}

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.water_props = WaterParameterBlock(solute_list=["A", "B", "C"])

        with pytest.raises(
            KeyError,
            match="fs.unit - database provided does not "
            "contain an entry for order w.r.t. species "
            "B in reaction Rxn3.",
        ):
            m.fs.unit = DerivedSIDOR(property_package=m.fs.water_props, database=m.db)

    @pytest.mark.unit
    def test_order_no_key_mw(self):
        m = ConcreteModel()

        m.db = Database(dbpath=local_path)

        # Load data from YAML so that we can modify it
        m.db._get_technology("test_sidor_data")
        m.db._cached_files["test_sidor_data"]["default"]["reactions"]["Rxn3"] = {}
        R3 = m.db._cached_files["test_sidor_data"]["default"]["reactions"]["Rxn3"]
        R3["conversion"] = 0.5
        R3["key_reactant"] = "B"
        R3["stoichiometry"] = {
            "A": {"order": 1, "molecular_weight": 1},
            "B": {"order": 1},
        }

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.water_props = WaterParameterBlock(solute_list=["A", "B", "C"])

        with pytest.raises(
            KeyError,
            match="fs.unit - database provided does not "
            "contain an entry for molecular_weight w.r.t. "
            "species B in reaction Rxn3.",
        ):
            m.fs.unit = DerivedSIDOR(property_package=m.fs.water_props, database=m.db)
