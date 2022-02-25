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
from watertap.property_models.coagulation_prop_pack import CoagulationParameterBlock
from watertap.unit_models.coag_floc_zo_jartest import CoagulationFlocculationZO_JarTestModel
from pyomo.environ import (ConcreteModel,
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
                           TerminationCondition)
from idaes.core import (FlowsheetBlock,
                        MaterialFlowBasis,
                        MaterialBalanceType,
                        EnergyBalanceType)
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.util.check_units import assert_units_consistent
import idaes.core.util.scaling as iscale
from idaes.core.util.scaling import badly_scaled_var_generator
from idaes.core.util import get_solver

__author__ = "Austin Ladshaw"

solver = get_solver()

# -----------------------------------------------------------------------------
# Start test class
class TestCoagulationZOJarTest_withChemicals():
    @pytest.fixture(scope="class")
    def coag_obj_w_chems(self):
        model = ConcreteModel()
        model.fs = FlowsheetBlock(default={"dynamic": False})
        model.fs.properties = CoagulationParameterBlock()
        chem_dict = {'Alum':
                        {"parameter_data":
                            {"mw_additive": (200, pyunits.g/pyunits.mol),
                            "moles_salt_per_mole_additive": 3,
                            "mw_salt": (100, pyunits.g/pyunits.mol)}
                        },
                      'Poly':
                        {"parameter_data":
                            {"mw_additive": (25, pyunits.g/pyunits.mol),
                            "moles_salt_per_mole_additive": 0,
                            "mw_salt": (23, pyunits.g/pyunits.mol)}
                        }
                     }
        model.fs.unit = CoagulationFlocculationZO_JarTestModel(default={
            "property_package": model.fs.properties,
            "chemical_additives": chem_dict })

        return model

    @pytest.mark.unit
    def test_build_model(self, coag_obj_w_chems):
        model = coag_obj_w_chems

        assert len(model.fs.unit.config.chemical_additives) == 2
        assert isinstance(model.fs.unit.slope, Var)
        assert isinstance(model.fs.unit.intercept, Var)
        assert isinstance(model.fs.unit.initial_turbidity_ntu, Var)
        assert isinstance(model.fs.unit.final_turbidity_ntu, Var)
        assert isinstance(model.fs.unit.chemical_doses, Var)
        assert len(model.fs.unit.chemical_doses) == 2
        assert isinstance(model.fs.unit.chemical_mw, Param)
        assert len(model.fs.unit.chemical_mw) == 2
        assert isinstance(model.fs.unit.salt_mw, Param)
        assert len(model.fs.unit.salt_mw) == 2
        assert isinstance(model.fs.unit.salt_from_additive_mole_ratio, Param)
        assert len(model.fs.unit.salt_from_additive_mole_ratio) == 2


# -----------------------------------------------------------------------------
# Start test class without chemicals added
class TestCoagulationZOJarTest_withNoChemicals():
    @pytest.fixture(scope="class")
    def coag_obj_wo_chems(self):
        model = ConcreteModel()
        model.fs = FlowsheetBlock(default={"dynamic": False})
        model.fs.properties = CoagulationParameterBlock()
        model.fs.unit = CoagulationFlocculationZO_JarTestModel(default={
            "property_package": model.fs.properties})

        return model

    @pytest.mark.unit
    def test_build_model(self, coag_obj_wo_chems):
        model = coag_obj_wo_chems

        assert len(model.fs.unit.config.chemical_additives) == 0
        assert isinstance(model.fs.unit.slope, Var)
        assert isinstance(model.fs.unit.intercept, Var)
        assert isinstance(model.fs.unit.initial_turbidity_ntu, Var)
        assert isinstance(model.fs.unit.final_turbidity_ntu, Var)
        assert isinstance(model.fs.unit.chemical_doses, Var)
        assert len(model.fs.unit.chemical_doses) == 0
        assert isinstance(model.fs.unit.chemical_mw, Param)
        assert len(model.fs.unit.chemical_mw) == 0
        assert isinstance(model.fs.unit.salt_mw, Param)
        assert len(model.fs.unit.salt_mw) == 0
        assert isinstance(model.fs.unit.salt_from_additive_mole_ratio, Param)
        assert len(model.fs.unit.salt_from_additive_mole_ratio) == 0

# -----------------------------------------------------------------------------
# Start test class with bad config
class TestCoagulationZOJarTest_withBadConfig():
    @pytest.fixture(scope="class")
    def coag_obj_bad_config(self):
        model = ConcreteModel()
        model.fs = FlowsheetBlock(default={"dynamic": False})
        model.fs.properties = CoagulationParameterBlock()

        return model

    @pytest.mark.unit
    def test_build_model_catch_errors(self, coag_obj_bad_config):
        model = coag_obj_bad_config

        bad_dict1 = {'Alum':
                        {"foo_bar":
                            {"mw_additive": (200, pyunits.g/pyunits.mol),
                            "moles_salt_per_mole_additive": 3,
                            "mw_salt": (100, pyunits.g/pyunits.mol)}
                        }
                     }
        with pytest.raises(ConfigurationError):
            model.fs.unit = CoagulationFlocculationZO_JarTestModel(default={
                "property_package": model.fs.properties,
                "chemical_additives": bad_dict1 })

        bad_dict2 = {'Alum':
                        {"parameter_data":
                            {"foo_bar": (200, pyunits.g/pyunits.mol),
                            "moles_salt_per_mole_additive": 3,
                            "mw_salt": (100, pyunits.g/pyunits.mol)}
                        }
                     }
        with pytest.raises(ConfigurationError):
            model.fs.unit = CoagulationFlocculationZO_JarTestModel(default={
                "property_package": model.fs.properties,
                "chemical_additives": bad_dict2 })

        bad_dict3 = {'Alum':
                        {"parameter_data":
                            {"mw_additive": (200, pyunits.g/pyunits.mol),
                            "moles_salt_per_mole_additive": "foo-bar",
                            "mw_salt": (100, pyunits.g/pyunits.mol)}
                        }
                     }
        with pytest.raises(ConfigurationError):
            model.fs.unit = CoagulationFlocculationZO_JarTestModel(default={
                "property_package": model.fs.properties,
                "chemical_additives": bad_dict3 })

        bad_dict4 = {'Alum':
                        {"parameter_data":
                            {"mw_additive": (200, pyunits.g/pyunits.mol),
                            "foo-bar": 3,
                            "mw_salt": (100, pyunits.g/pyunits.mol)}
                        }
                     }
        with pytest.raises(ConfigurationError):
            model.fs.unit = CoagulationFlocculationZO_JarTestModel(default={
                "property_package": model.fs.properties,
                "chemical_additives": bad_dict4 })

        bad_dict5 = {'Alum':
                        {"parameter_data":
                            {"mw_additive": (200, pyunits.g/pyunits.mol),
                            "moles_salt_per_mole_additive": 3,
                            "foo-bar": (100, pyunits.g/pyunits.mol)}
                        }
                     }
        with pytest.raises(ConfigurationError):
            model.fs.unit = CoagulationFlocculationZO_JarTestModel(default={
                "property_package": model.fs.properties,
                "chemical_additives": bad_dict5 })

        bad_dict6 = {'Alum':
                        {"parameter_data":
                            {"mw_additive": "not a tuple",
                            "moles_salt_per_mole_additive": 3,
                            "mw_salt": (100, pyunits.g/pyunits.mol)}
                        }
                     }
        with pytest.raises(ConfigurationError):
            model.fs.unit = CoagulationFlocculationZO_JarTestModel(default={
                "property_package": model.fs.properties,
                "chemical_additives": bad_dict6 })
