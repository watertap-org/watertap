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

from pyomo.environ import (ConcreteModel,
                           Set,
                           Param,
                           Var,
                           Expression,
                           Constraint,
                           Block,
                           units as pyunits)
from pyomo.util.check_units import assert_units_consistent, check_units_equivalent
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_variables,
                                              number_total_constraints,
                                              number_unused_variables)
from idaes.core.util.scaling import (calculate_scaling_factors,
                                     unscaled_variables_generator,
                                     unscaled_constraints_generator,
                                     badly_scaled_var_generator)
from idaes.core.util import get_solver

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()
solver.options["nlp_scaling_method"] = "user-scaling"


# -----------------------------------------------------------------------------
class PropertyUnitTestHarness(object):
    def configure_class(self, m):
        self.prop_pack = None
        self.param_args = {}
        self.set_default_scaling_dict = {}  # keys = var string, value = list of (scaling factor, tuple index)

        self.component_dict = {}  # keys = component string, value = component type
        self.phase_dict = {}  # keys = phase string, value = phase type
        self.param_obj_dict = {}  # keys = pyomo object string, value = (type, units)
        self.default_scaling_factor_dict = {}  # keys = (var string, index), value = scaling factor

        self.state_vars_list = None
        self.stateblock_obj_dict = {}  # keys = pyomo object string, value = (type, units)

        self.statistics_dict = {}  # keys = statistic string, value = number
        self.general_methods_dict = {}  # keys = method name string, value = (var, list of tuple indices)
        self.default_methods_dict = {}  # keys = method name string, value = class attribute

        self.configure()

        # attaching objects to model
        m.component_dict = self.component_dict
        m.phase_dict = self.phase_dict
        m.param_obj_dict = self.param_obj_dict
        m.default_scaling_factor_dict = self.default_scaling_factor_dict

        m.stateblock_obj_dict = self.stateblock_obj_dict
        m.state_vars_list = self.state_vars_list

        m.statistics_dict = self.statistics_dict
        m.general_methods_dict = self.general_methods_dict
        m.default_methods_dict = self.default_methods_dict


    def configure(self):
        # Placeholder method to allow user to setup test harness
        pass

    @staticmethod
    def check_block_obj(blk, obj_dict):
        """check that all objects are their expected type and have their expected units"""
        for (obj_str, (obj_type, obj_units)) in obj_dict.items():
            obj = getattr(blk, obj_str)
            if not isinstance(obj, obj_type):
                raise Exception(
                    "{obj_str} is in the testing dictionary, but is not on {blk} or is not"
                    "the expected type. Update the appropriate dictionary in the configure "
                    "function that sets up the PropertyUnitTestHarness".format(
                        obj_str=obj_str, blk=blk))
            if obj_type in [Var, Expression]:
                if not check_units_equivalent(obj, obj_units):
                    raise Exception(
                        "{obj_str} on {blk} has units of {obj_units_1}, but it was assigned units of "
                        "{obj_units} in the testing dictionary. Update the appropriate dictionary with "
                        "the {obj_str} key in the configure function that set up the "
                        "PropertyUnitTestHarness".format(obj_str=obj_str,
                                                         blk=blk,
                                                         obj_units=obj_units,
                                                         obj_units_1=pyunits.get_units(obj)))

    @staticmethod
    def check_block_obj_coverage(blk, obj_dict,
                         type_tpl=(Param, Var, Expression, Constraint)):
        """check that all pyomo objects on block are in the testing dictionary"""
        for obj in blk.component_objects(type_tpl, descend_into=False):
            obj_str = obj.local_name
            if obj_str not in obj_dict:
                raise Exception(
                    "{obj_str} is on {blk}, but is not included in the testing dictionary. "
                    "Add {obj_str} to the appropriate dictionary in the configure function"
                    "that sets up the PropertyUnitTestHarness".format(obj_str=obj_str, blk=blk))

    @pytest.fixture(scope='class')
    def frame(self):
        m = ConcreteModel()
        self.configure_class(m)

        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.properties = self.prop_pack()
        m.fs.stream = m.fs.properties.build_state_block(
            [0], default=self.param_args)

        for v_str, temp_list in self.set_default_scaling_dict.items():
            for temp in temp_list:
                sf = temp[0]
                ind = temp[1]
                m.fs.properties.set_default_scaling(v_str, sf, index=ind)

        return m

    @pytest.mark.unit
    def test_components_phases(self, frame):
        m = frame
        # components
        assert isinstance(m.fs.properties.component_list, Set)
        assert len(m.fs.properties.component_list) == len(m.component_dict)
        for c_str in m.component_dict:
            assert c_str in m.fs.properties.component_list
            c = getattr(m.fs.properties, c_str)
            c_type = m.component_dict[c_str]
            assert isinstance(c, c_type)

        # phases
        assert isinstance(m.fs.properties.phase_list, Set)
        assert len(m.fs.properties.phase_list) == len(m.phase_dict)
        for p_str in m.phase_dict:
            assert p_str in m.fs.properties.phase_list
            p = getattr(m.fs.properties, p_str)
            p_type = m.phase_dict[p_str]
            assert isinstance(p, p_type)

    @pytest.mark.unit
    def test_parameters(self, frame):
        m = frame

        self.check_block_obj(m.fs.properties, m.param_obj_dict)
        self.check_block_obj_coverage(m.fs.properties, m.param_obj_dict)

        # test that the parameter variables are fixed
        for v in m.fs.properties.component_objects(Var):
            assert v.is_fixed()

    @pytest.mark.unit
    def test_default_scaling(self, frame):
        m = frame

        assert isinstance(m.fs.properties.default_scaling_factor, dict)
        assert len(m.default_scaling_factor_dict) == len(m.fs.properties.default_scaling_factor)
        for t, sf in m.fs.properties.default_scaling_factor.items():
            assert t in m.default_scaling_factor_dict.keys()
            assert m.default_scaling_factor_dict[t] == sf

    @pytest.mark.unit
    def test_state_variables(self, frame):
        m = frame
        # state variables
        state_vars_dict = m.fs.stream[0].define_state_vars()
        assert len(state_vars_dict) == len(m.state_vars_list)
        for v_str in state_vars_dict:
            assert v_str in m.state_vars_list

    @pytest.mark.unit
    def test_metadata(self, frame):
        m = frame
        assert hasattr(m.fs.properties, 'define_metadata')
        metadata = m.fs.properties.get_metadata().properties

        # check method naming convention
        for v_str in metadata.keys():
            if v_str in m.state_vars_list:
                assert metadata[v_str]['method'] is None
            else:
                assert metadata[v_str]['method'] == '_' + v_str

    @pytest.mark.unit
    def test_build(self, frame):
        m = frame

        # properties and state variables
        prop_str_list = []
        for obj_str, (obj_type, obj_units) in m.stateblock_obj_dict.items():
            if obj_type is Expression or obj_type is Var:
                prop_str_list.append(obj_str)

        # check that the testing property list matches metadata
        metadata = m.fs.properties.get_metadata().properties
        assert len(metadata) == len(prop_str_list)
        for v_str in metadata.keys():
            assert v_str in prop_str_list

        # check that properties are not built if not demanded
        for obj_str in prop_str_list:
            if obj_str not in m.state_vars_list:
                assert not m.fs.stream[0].is_property_constructed(obj_str)

        # check that properties are built on demand
        for obj_str in prop_str_list:
            assert hasattr(m.fs.stream[0], obj_str)

        # add constraints for on demand variables to testing dictionary
        for obj_str in prop_str_list:
            if (m.stateblock_obj_dict[obj_str][0] is Var
                    and obj_str not in m.state_vars_list):
                m.stateblock_obj_dict['eq_' + obj_str] = (Constraint, None)

        self.check_block_obj(m.fs.stream[0], m.stateblock_obj_dict)
        self.check_block_obj_coverage(m.fs.stream[0], m.stateblock_obj_dict)

    @pytest.mark.unit
    def test_statistics(self, frame):
        m = frame

        assert number_variables(m) == m.statistics_dict['number_variables']
        assert number_total_constraints(m) == m.statistics_dict['number_total_constraints']
        assert number_unused_variables(m) == m.statistics_dict['number_unused_variables']

    @pytest.mark.unit
    def test_general_methods(self, frame):
        m = frame

        # TODO: update for dynamic compatible property packages
        # general methods requirements
        required_general_methods_list = ['get_material_flow_terms',
                                         'get_enthalpy_flow_terms']
        for i in required_general_methods_list:
            assert i in m.general_methods_dict

        for k, (v_str, ind_list) in m.general_methods_dict.items():
            assert hasattr(m.fs.stream[0], k)
            general_method = getattr(m.fs.stream[0], k)
            var = getattr(m.fs.stream[0], v_str)
            for ind in ind_list:
                if ind is None:
                    assert var == general_method(None)
                elif len(ind) == 1:
                    assert var[ind] == general_method(ind[0])
                elif len(ind) == 2:
                    assert var[ind] == general_method(ind[0], ind[1])
                else:
                    raise Exception(
                        "general method {k} was assigned an index of {ind} with a length of "
                        "{ind_length}, but a length of 2 or less is expected.".format(
                            k=k, ind=ind, ind_length=len(ind)))

    @pytest.mark.unit
    def test_default_methods(self, frame):
        m = frame

        # default methods requirements
        required_default_methods_list = ['default_material_balance_type',
                                         'default_energy_balance_type',
                                         'get_material_flow_basis']
        for i in required_default_methods_list:
            assert i in m.default_methods_dict

        for k, v in m.default_methods_dict.items():
            assert hasattr(m.fs.stream[0], k)
            assert getattr(m.fs.stream[0], k)() == v

    @pytest.mark.unit
    def test_scaling(self, frame):
        m = frame

        calculate_scaling_factors(m.fs)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0
        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m))
        assert len(unscaled_constraint_list) == 0



    # @pytest.mark.unit
    # def test_default_scaling(self, frame):
    #     pass
    #
    # @pytest.mark.unit
    # def test_fail(self):
    #     assert False
