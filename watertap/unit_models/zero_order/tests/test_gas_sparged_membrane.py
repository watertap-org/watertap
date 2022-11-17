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
Tests for zero-order gas-sparged membrane unit
"""
import pytest


from idaes.core import declare_process_block_class, FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale
from idaes.config import bin_directory as idaes_bin_directory

from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    value,
    Var,
)
from pyomo.network import Port
from pyomo.util.check_units import assert_units_consistent


from watertap.core import WaterParameterBlock, WaterStateBlock
from watertap.core.wt_database import Database
from watertap.unit_models.zero_order import GasSpargedMembraneZO
from watertap.unit_models.zero_order.gas_sparged_membrane_zo import (
    initialize_sido,
    calculate_scaling_factors_gas_extraction,
    _get_Q_gas_extraction,
)

solver = get_solver()


class TestGasSpargedMembraneZO:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.water_props = WaterParameterBlock(solute_list=["cod"])

        m.fs.unit = GasSpargedMembraneZO(
            property_package=m.fs.water_props, database=m.db
        )

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(1000)
        m.fs.unit.inlet.flow_mass_comp[0, "cod"].fix(1)

        return m

    @pytest.mark.unit
    def test_private_attributes(self, model):
        assert model.fs.unit.config.database is model.db

        assert model.fs.unit._has_recovery_removal is True
        assert model.fs.unit._initialize is initialize_sido
        assert model.fs.unit._scaling is calculate_scaling_factors_gas_extraction
        assert model.fs.unit._get_Q is _get_Q_gas_extraction
        assert model.fs.unit._stream_table_dict == {
            "Inlet": model.fs.unit.inlet,
            "Treated": model.fs.unit.treated,
            "Byproduct": model.fs.unit.byproduct,
        }
        assert model.fs.unit._perf_var_dict == {
            "Water Recovery": model.fs.unit.recovery_frac_mass_H2O,
            "Solute Removal": model.fs.unit.removal_frac_mass_comp,
            "Mass of gas extracted per mass flow of influent(kg/d/(kg/d)": model.fs.unit.gas_mass_influent_ratio,
            "Mass flow of gas extracted (kg/s))": model.fs.unit.flow_mass_gas_extraction,
            "Electricity Demand": model.fs.unit.electricity,
        }

    @pytest.mark.unit
    def test_build(self, model):
        assert isinstance(model.fs.unit.properties_in, WaterStateBlock)
        assert isinstance(model.fs.unit.properties_treated, WaterStateBlock)
        assert isinstance(model.fs.unit.properties_byproduct, WaterStateBlock)

        assert isinstance(model.fs.unit.inlet, Port)
        assert isinstance(model.fs.unit.treated, Port)
        assert isinstance(model.fs.unit.byproduct, Port)

        assert isinstance(model.fs.unit.gas_mass_influent_ratio, Var)
        assert len(model.fs.unit.gas_mass_influent_ratio) == 1
        assert isinstance(model.fs.unit.recovery_frac_mass_H2O, Var)
        assert len(model.fs.unit.recovery_frac_mass_H2O) == 1
        assert isinstance(model.fs.unit.removal_frac_mass_comp, Var)
        assert len(model.fs.unit.removal_frac_mass_comp) == 1

        assert isinstance(model.fs.unit.water_recovery_equation, Constraint)
        assert len(model.fs.unit.water_recovery_equation) == 1
        assert isinstance(model.fs.unit.mass_balance, Constraint)
        assert len(model.fs.unit.mass_balance) == 1
        assert isinstance(model.fs.unit.solute_removal_equation, Constraint)
        assert len(model.fs.unit.solute_removal_equation) == 1
        assert isinstance(model.fs.unit.solute_treated_equation, Constraint)
        assert len(model.fs.unit.solute_treated_equation) == 1
        assert isinstance(model.fs.unit.mass_gas_extraction_equation, Constraint)
        assert len(model.fs.unit.mass_gas_extraction_equation) == 1

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("gas_sparged_membrane")

        model.fs.unit.load_parameters_from_database()

        assert model.fs.unit.recovery_frac_mass_H2O[0].fixed
        assert (
            model.fs.unit.recovery_frac_mass_H2O[0].value
            == data["recovery_frac_mass_H2O"]["value"]
        )

        for (t, j), v in model.fs.unit.removal_frac_mass_comp.items():
            assert v.fixed
            assert v.value == data["removal_frac_mass_comp"][j]["value"]

        assert model.fs.unit.gas_mass_influent_ratio[0].fixed
        assert (
            model.fs.unit.gas_mass_influent_ratio[0].value
            == data["gas_mass_influent_ratio"]["value"]
        )

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
                model.fs.unit.mass_balance[0]
            )
            == 1e5
        )
        assert (
            iscale.get_constraint_transform_applied_scaling_factor(
                model.fs.unit.solute_removal_equation[0, "cod"]
            )
            == 1e5
        )
        assert (
            iscale.get_constraint_transform_applied_scaling_factor(
                model.fs.unit.solute_treated_equation[0, "cod"]
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
        assert pytest.approx(979, rel=1e-5) == value(
            model.fs.unit.treated.flow_mass_comp[0, "H2O"]
        )
        assert pytest.approx(20.58784, rel=1e-5) == value(
            model.fs.unit.byproduct.flow_mass_comp[0, "H2O"]
        )
        assert pytest.approx(0, abs=1e-5) == value(
            model.fs.unit.treated.flow_mass_comp[0, "cod"]
        )
        assert pytest.approx(1, rel=1e-5) == value(
            model.fs.unit.byproduct.flow_mass_comp[0, "cod"]
        )

        assert pytest.approx(1.001, rel=1e-5) == value(
            model.fs.unit.properties_in[0].flow_vol
        )
        assert pytest.approx(0.9990001, rel=1e-5) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["cod"]
        )
        assert pytest.approx(0.979, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].flow_vol
        )
        assert pytest.approx(0, abs=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["cod"]
        )
        assert pytest.approx(0.021588, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].flow_vol
        )
        assert pytest.approx(46.32238, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["cod"]
        )
        assert pytest.approx(368.855328, abs=1e-5) == value(
            model.fs.unit.electricity[0]
        )
        assert pytest.approx(
            5 * value(model.fs.unit.properties_in[0].flow_vol) * 0.08235, rel=1e-5
        ) == value(model.fs.unit.flow_mass_gas_extraction[0])

    @pytest.mark.component
    def test_conservation(self, model):
        for (t, j) in model.fs.unit.inlet.flow_mass_comp.keys():
            if j != "H2O":
                assert (
                    abs(
                        value(
                            model.fs.unit.inlet.flow_mass_comp[t, j]
                            - model.fs.unit.treated.flow_mass_comp[t, j]
                            - model.fs.unit.byproduct.flow_mass_comp[t, j]
                        )
                    )
                    <= 1e-6
                )
            else:
                assert (
                    abs(
                        value(
                            model.fs.unit.inlet.flow_mass_comp[t, j]
                            - model.fs.unit.treated.flow_mass_comp[t, j]
                            - model.fs.unit.byproduct.flow_mass_comp[t, j]
                            - model.fs.unit.flow_mass_gas_extraction[t]
                        )
                    )
                    <= 1e-6
                )

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_report(self, model):

        model.fs.unit.report()
