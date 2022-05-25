###############################################################################
# ProteusLib Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/nawi-hub/proteuslib/"
#
###############################################################################
import pytest
from watertap.property_models.NDMA_prop_pack import (
    NDMAParameterBlock,
    NDMAStateBlock,
)
from watertap.property_models.tests.property_test_harness import PropertyAttributeError
from pyomo.environ import (
    ConcreteModel,
    assert_optimal_termination,
    value,
    Set,
    Param,
    Var,
    units as pyunits,
    Suffix,
    Constraint,
    SolverFactory,
    SolverStatus,
    TerminationCondition,
)
from idaes.core import (
    FlowsheetBlock,
    MaterialFlowBasis,
    MaterialBalanceType,
    EnergyBalanceType,
)
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.util.check_units import assert_units_consistent
import idaes.core.util.scaling as iscale
from idaes.core.util.scaling import badly_scaled_var_generator
from idaes.core.util import get_solver

solver = get_solver()

# Start test class
class TestNDMAPropPack:
    @pytest.fixture(scope="class")
    def NDMA_obj(self):
        model = ConcreteModel()
        model.fs = FlowsheetBlock(default={"dynamic": False})
        model.fs.properties = NDMAParameterBlock()

        return model

    @pytest.mark.unit
    def test_build_properties(self, NDMA_obj):
        model = NDMA_obj

        # Check to make sure we have the correct components and phases
        assert isinstance(model.fs.properties.component_list, Set)
        for j in model.fs.properties.component_list:
            assert j in ["H2O", "NDMA"]
        assert isinstance(model.fs.properties.phase_list, Set)
        for p in model.fs.properties.phase_list:
            assert p in ["Liq"]

        # Check for existance of the state block object
        assert model.fs.properties._state_block_class is NDMAStateBlock

        # Assert that we have the expected set of parameters
        assert isinstance(model.fs.properties.cp, Param)

    @pytest.mark.unit
    def test_metadata(self, NDMA_obj):
        model = NDMA_obj

        # Create the state block and pull in the default metadata
        model.fs.stream = model.fs.properties.build_state_block([0], default={})
        metadata = model.fs.properties.get_metadata().properties

        # check that properties are not built if not demanded
        for v_name in metadata:
            if metadata[v_name]["method"] is not None:
                if model.fs.stream[0].is_property_constructed(v_name):
                    raise PropertyAttributeError(
                        "Property {v_name} is an on-demand property, but was found "
                        "on the stateblock without being demanded".format(v_name=v_name)
                    )

        # check that properties are built if demanded
        for v_name in metadata:
            if metadata[v_name]["method"] is not None:
                if not hasattr(model.fs.stream[0], v_name):
                    raise PropertyAttributeError(
                        "Property {v_name} is an on-demand property, but was not built "
                        "when demanded".format(v_name=v_name)
                    )

        # check the other stateblock functions
        res = model.fs.stream[0].define_state_vars()
        assert res["flow_mass_phase_comp"] == model.fs.stream[0].flow_mass_phase_comp
        assert res["pressure"] == model.fs.stream[0].pressure
        assert res["temperature"] == model.fs.stream[0].temperature
        assert (
            model.fs.stream[0].get_material_flow_terms("Liq", "H2O")
            == model.fs.stream[0].flow_mass_phase_comp["Liq", "H2O"]
        )
        # assert (
        #     model.fs.stream[0].default_energy_balance_type()
        #     == EnergyBalanceType.enthalpyTotal
        # )
        # assert model.fs.stream[0].get_material_flow_basis() == MaterialFlowBasis.mass
        # assert (
        #     model.fs.stream[0].get_enthalpy_flow_terms("Liq")
        #     == model.fs.stream[0].enth_flow["Liq"]
        # )