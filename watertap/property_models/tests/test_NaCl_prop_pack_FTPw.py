#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import pytest
import watertap.property_models.NaCl_prop_pack_FTPw as props
from idaes.models.properties.tests.test_harness import (
    PropertyTestHarness as PropertyTestHarness_idaes,
)
from watertap.property_models.tests.property_test_harness import (
    PropertyTestHarness,
    PropertyRegressionTest,
    PropertyCalculateStateTest,
)


# -----------------------------------------------------------------------------
@pytest.mark.unit
class TestNaClProperty_idaes(PropertyTestHarness_idaes):
    def configure(self):
        self.prop_pack = props.NaClFTPwParameterBlock
        self.param_args = {}
        self.prop_args = {}
        self.has_density_terms = False


class TestNaClProperty(PropertyTestHarness):
    defined_state = True

    def configure(self):
        self.prop_pack = props.NaClFTPwParameterBlock
        self.param_args = {"defined_state": self.defined_state}
        self.scaling_args = {
            ("flow_mass_phase", "Liq"): 1,
            ("mass_frac_phase_comp", ("Liq", "H2O")): 1,
            ("mass_frac_phase_comp", ("Liq", "NaCl")): 1e2,
        }
        if self.defined_state:
            self.stateblock_statistics = {
                "number_variables": 21,
                "number_total_constraints": 16,
                "number_unused_variables": 0,
                "default_degrees_of_freedom": 5,
            }
        else:
            self.stateblock_statistics = {
                "number_variables": 21,
                "number_total_constraints": 17,
                "number_unused_variables": 0,
                "default_degrees_of_freedom": 4,
            }

        self.default_solution = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 0.965,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 0.035,
            ("dens_mass_phase", "Liq"): 1021.5,
            ("flow_vol_phase", "Liq"): 9.790e-4,
            ("conc_mass_phase_comp", ("Liq", "H2O")): 985.7,
            ("conc_mass_phase_comp", ("Liq", "NaCl")): 35.75,
            ("flow_mol_phase_comp", ("Liq", "H2O")): 53.57,
            ("flow_mol_phase_comp", ("Liq", "NaCl")): 0.5989,
            ("mole_frac_phase_comp", ("Liq", "H2O")): 0.9889,
            ("mole_frac_phase_comp", ("Liq", "NaCl")): 1.106e-2,
            ("molality_phase_comp", ("Liq", "NaCl")): 0.6206,
            ("diffus_phase_comp", ("Liq", "NaCl")): 1.472e-9,
            ("visc_d_phase", "Liq"): 1.055e-3,
            ("osm_coeff", None): 0.9271,
            ("pressure_osm_phase", "Liq"): 2.853e6,
            ("enth_mass_phase", "Liq"): 1.045e5,
        }


class TestNaClPropertyDefinedState(TestNaClProperty):
    defined_state = True


@pytest.mark.component
class TestNaClPropertySolution_1(PropertyRegressionTest):
    defined_state = False

    def configure(self):
        self.prop_pack = props.NaClFTPwParameterBlock
        self.param_args = {"defined_state": self.defined_state}

        self.solver = "ipopt"
        self.optarg = {"nlp_scaling_method": "user-scaling"}

        self.scaling_args = {
            ("flow_mass_phase", "Liq"): 1,
            ("mass_frac_phase_comp", ("Liq", "H2O")): 1,
            ("mass_frac_phase_comp", ("Liq", "NaCl")): 1e2,
        }
        self.state_args = {
            ("flow_mass_phase", "Liq"): 1,
            ("mass_frac_phase_comp", ("Liq", "NaCl")): 0.05,
            ("temperature", None): 273.15 + 25,
            ("pressure", None): 50e5,
        }
        if self.defined_state:
            self.state_args[("mass_frac_phase_comp", ("Liq", "H2O"))] = 0.95
        self.regression_solution = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 0.95,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 0.05,
            ("dens_mass_phase", "Liq"): 1032.8,
            ("flow_vol_phase", "Liq"): 9.682e-4,
            ("conc_mass_phase_comp", ("Liq", "H2O")): 981.1,
            ("conc_mass_phase_comp", ("Liq", "NaCl")): 51.64,
            ("flow_mol_phase_comp", ("Liq", "H2O")): 52.73,
            ("flow_mol_phase_comp", ("Liq", "NaCl")): 0.8556,
            ("mole_frac_phase_comp", ("Liq", "H2O")): 0.9840,
            ("mole_frac_phase_comp", ("Liq", "NaCl")): 1.597e-2,
            ("molality_phase_comp", ("Liq", "NaCl")): 0.9006,
            ("diffus_phase_comp", ("Liq", "NaCl")): 1.471e-9,
            ("visc_d_phase", "Liq"): 1.0875e-3,
            ("osm_coeff", None): 0.9347,
            ("pressure_osm_phase", "Liq"): 4.174e6,
            ("enth_mass_phase", "Liq"): 1.093e5,
        }


class TestNaClPropertySolution_1DefinedState(TestNaClPropertySolution_1):
    defined_state = True


@pytest.mark.component
class TestNaClPropertySolution_2(PropertyRegressionTest):
    defined_state = False

    def configure(self):
        self.prop_pack = props.NaClFTPwParameterBlock
        self.param_args = {"defined_state": self.defined_state}

        self.solver = "ipopt"
        self.optarg = {"nlp_scaling_method": "user-scaling"}

        self.scaling_args = {
            ("flow_mass_phase", "Liq"): 1,
            ("mass_frac_phase_comp", ("Liq", "H2O")): 1,
            ("mass_frac_phase_comp", ("Liq", "NaCl")): 1e2,
        }
        self.state_args = {
            ("flow_mass_phase", "Liq"): 1,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 0.26,
            ("temperature", None): 273.15 + 25,
            ("pressure", None): 50e5,
        }
        if self.defined_state:
            self.state_args[("flow_mass_phase_comp", ("Liq", "H2O"))] = 0.74
        self.regression_solution = {
            ("mass_frac_phase_comp", ("Liq", "H2O")): 0.74,
            ("mass_frac_phase_comp", ("Liq", "NaCl")): 0.26,
            ("dens_mass_phase", "Liq"): 1192,
            ("flow_vol_phase", "Liq"): 8.392e-4,
            ("conc_mass_phase_comp", ("Liq", "H2O")): 881.8,
            ("conc_mass_phase_comp", ("Liq", "NaCl")): 309.8,
            ("flow_mol_phase_comp", ("Liq", "H2O")): 41.08,
            ("flow_mol_phase_comp", ("Liq", "NaCl")): 4.449,
            ("mole_frac_phase_comp", ("Liq", "H2O")): 0.9022,
            ("mole_frac_phase_comp", ("Liq", "NaCl")): 9.773e-2,
            ("molality_phase_comp", ("Liq", "NaCl")): 6.012,
            ("diffus_phase_comp", ("Liq", "NaCl")): 1.580e-9,
            ("visc_d_phase", "Liq"): 1.539e-3,
            ("osm_coeff", None): 1.274,
            ("pressure_osm_phase", "Liq"): 3.796e7,
            ("enth_mass_phase", "Liq"): 1.087e5,
        }


class TestNaClPropertySolution_2DefinedState(TestNaClPropertySolution_2):
    defined_state = True


@pytest.mark.component
class TestNaClPropertySolution_3(PropertyRegressionTest):
    defined_state = False

    def configure(self):
        self.prop_pack = props.NaClFTPwParameterBlock
        self.param_args = {"defined_state": self.defined_state}

        self.solver = "ipopt"
        self.optarg = {"nlp_scaling_method": "user-scaling"}

        self.scaling_args = {
            ("flow_mass_phase", "Liq"): 1,
            ("mass_frac_phase_comp", ("Liq", "H2O")): 1,
            ("mass_frac_phase_comp", ("Liq", "NaCl")): 1e2,
        }
        self.state_args = {
            ("flow_mass_phase", "Liq"): 1,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 0.001,
            ("temperature", None): 273.15 + 25,
            ("pressure", None): 10e5,
        }
        if self.defined_state:
            self.state_args[("flow_mass_phase_comp", ("Liq", "H2O"))] = 0.999
        self.regression_solution = {
            ("mass_frac_phase_comp", ("Liq", "H2O")): 0.999,
            ("mass_frac_phase_comp", ("Liq", "NaCl")): 0.001,
            ("dens_mass_phase", "Liq"): 995.8,
            ("flow_vol_phase", "Liq"): 1.004e-3,
            ("conc_mass_phase_comp", ("Liq", "H2O")): 994.8,
            ("conc_mass_phase_comp", ("Liq", "NaCl")): 0.9958,
            ("flow_mol_phase_comp", ("Liq", "H2O")): 55.45,
            ("flow_mol_phase_comp", ("Liq", "NaCl")): 1.711e-2,
            ("mole_frac_phase_comp", ("Liq", "H2O")): 0.9997,
            ("mole_frac_phase_comp", ("Liq", "NaCl")): 3.084e-4,
            ("molality_phase_comp", ("Liq", "NaCl")): 1.713e-2,
            ("diffus_phase_comp", ("Liq", "NaCl")): 1.508e-9,
            ("visc_d_phase", "Liq"): 9.822e-4,
            ("osm_coeff", None): 0.918,
            ("pressure_osm_phase", "Liq"): 7.797e4,
            ("enth_mass_phase", "Liq"): 1.055e5,
        }


class TestNaClPropertySolution_3DefinedState(TestNaClPropertySolution_3):
    defined_state = True


@pytest.mark.component
class TestNaClCalculateState_1(PropertyCalculateStateTest):
    def configure(self):
        self.prop_pack = props.NaClFTPwParameterBlock
        self.param_args = {}

        self.solver = "ipopt_v2"

        self.scaling_args = {
            ("flow_mass_phase", "Liq"): 1e-1,
            ("mass_frac_phase_comp", ("Liq", "H2O")): 1,
            ("mass_frac_phase_comp", ("Liq", "NaCl")): 1e2,
        }
        self.var_args = {
            ("flow_vol_phase", "Liq"): 2e-2,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 1.033,
            ("temperature", None): 273.15 + 25,
            ("pressure", None): 5e5,
        }
        self.state_solution = {
            ("flow_mass_phase", "Liq"): 20.653,
            ("mass_frac_phase_comp", ("Liq", "NaCl")): 0.05,
            ("mass_frac_phase_comp", ("Liq", "H2O")): 0.95,
        }


@pytest.mark.component
class TestNaClCalculateState_2(PropertyCalculateStateTest):
    def configure(self):
        self.prop_pack = props.NaClFTPwParameterBlock
        self.param_args = {}

        self.solver = "ipopt_v2"

        self.scaling_args = {
            ("flow_mass_phase", "Liq"): 1e1,
            ("mass_frac_phase_comp", ("Liq", "H2O")): 1,
            ("mass_frac_phase_comp", ("Liq", "NaCl")): 1e2,
        }
        self.var_args = {
            ("flow_vol_phase", "Liq"): 2e-4,
            ("mole_frac_phase_comp", ("Liq", "NaCl")): 0.05,
            ("temperature", None): 273.15 + 25,
            ("pressure", None): 5e5,
        }
        self.state_solution = {
            ("flow_mass_phase", "Liq"): 0.22104,
            ("mass_frac_phase_comp", ("Liq", "NaCl")): 0.14586,
            ("mass_frac_phase_comp", ("Liq", "H2O")): 0.85414,
        }


@pytest.mark.component
class TestNaClCalculateState_3(PropertyCalculateStateTest):
    def configure(self):
        self.prop_pack = props.NaClFTPwParameterBlock
        self.param_args = {}

        self.solver = "ipopt_v2"

        self.scaling_args = {
            ("flow_mass_phase", "Liq"): 1,
            ("mass_frac_phase_comp", ("Liq", "H2O")): 1,
            ("mass_frac_phase_comp", ("Liq", "NaCl")): 1e2,
        }
        self.var_args = {
            ("flow_vol_phase", "Liq"): 1e-3,
            ("pressure_osm_phase", "Liq"): 100e5,
            ("temperature", None): 273.15 + 25,
            ("pressure", None): 5e5,
        }
        self.state_solution = {
            ("flow_mass_phase", "Liq"): 1.0759,
            ("mass_frac_phase_comp", ("Liq", "NaCl")): 0.10698,
            ("mass_frac_phase_comp", ("Liq", "H2O")): 0.89302,
        }
