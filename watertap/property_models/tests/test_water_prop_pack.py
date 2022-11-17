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
import pytest
import watertap.property_models.water_prop_pack as props
from idaes.models.properties.tests.test_harness import (
    PropertyTestHarness as PropertyTestHarness_idaes,
)
from watertap.property_models.tests.property_test_harness import (
    PropertyTestHarness,
    PropertyRegressionTest,
    PropertyCalculateStateTest,
)

# -----------------------------------------------------------------------------
class TestSeawaterProperty_idaes(PropertyTestHarness_idaes):
    def configure(self):
        self.prop_pack = props.WaterParameterBlock
        self.param_args = {}
        self.prop_args = {}
        self.has_density_terms = False


class TestSeawaterProperty(PropertyTestHarness):
    def configure(self):
        self.prop_pack = props.WaterParameterBlock
        self.param_args = {}
        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1,
            ("flow_mass_phase_comp", ("Vap", "H2O")): 1,
        }
        self.stateblock_statistics = {
            "number_variables": 20,
            "number_total_constraints": 16,
            "number_unused_variables": 0,
            "default_degrees_of_freedom": 4,
        }  # 4 state vars
        self.default_solution = {
            ("dens_mass_phase", "Liq"): 996.9,
            ("dens_mass_phase", "Vap"): 0.7363,
            ("flow_vol_phase", "Liq"): 5.016e-4,
            ("flow_vol_phase", "Vap"): 0.6790,
            ("flow_mol_phase_comp", ("Liq", "H2O")): 27.75,
            ("flow_mol_phase_comp", ("Vap", "H2O")): 27.75,
            ("mole_frac_phase_comp", ("Liq", "H2O")): 0.5,
            ("mole_frac_phase_comp", ("Vap", "H2O")): 0.5,
            ("enth_mass_phase", "Liq"): 1.049e5,
            ("enth_mass_phase", "Vap"): 2.547e6,
            ("dh_vap_mass", None): 2.442e6,
            ("cp_mass_phase", "Liq"): 4.187e3,
            ("cp_mass_phase", "Vap"): 1.865e3,
        }


@pytest.mark.component
class TestSeawaterPropertySolution_1(PropertyRegressionTest):
    def configure(self):
        self.prop_pack = props.WaterParameterBlock
        self.param_args = {}

        self.solver = "ipopt"
        self.optarg = {"nlp_scaling_method": "user-scaling"}

        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1,
            ("flow_mass_phase_comp", ("Vap", "H2O")): 1e8,
        }
        self.state_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1,
            ("flow_mass_phase_comp", ("Vap", "H2O")): 1e-8,
            ("temperature", None): 273.15 + 50,
            ("pressure", None): 2e5,
        }
        self.regression_solution = {
            ("dens_mass_phase", "Liq"): 988.05,
            ("dens_mass_phase", "Vap"): 1.341,
            ("flow_vol_phase", "Liq"): 1.012e-3,
            ("flow_vol_phase", "Vap"): 7.457e-9,
            ("flow_mol_phase_comp", ("Liq", "H2O")): 55.51,
            ("flow_mol_phase_comp", ("Vap", "H2O")): 5.551e-7,
            ("mole_frac_phase_comp", ("Liq", "H2O")): 1,
            ("mole_frac_phase_comp", ("Vap", "H2O")): 1e-8,
            ("enth_mass_phase", "Liq"): 2.094e5,
            ("enth_mass_phase", "Vap"): 2.591e6,
            ("dh_vap_mass", None): 2.382e6,
            ("cp_mass_phase", "Liq"): 4.181e3,
            ("cp_mass_phase", "Vap"): 1.871e3,
        }


@pytest.mark.component
class TestSeawaterPropertySolution_2(PropertyRegressionTest):
    def configure(self):
        self.prop_pack = props.WaterParameterBlock
        self.param_args = {}

        self.solver = "ipopt"
        self.optarg = {"nlp_scaling_method": "user-scaling"}

        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1e8,
            ("flow_mass_phase_comp", ("Vap", "H2O")): 1,
        }
        self.state_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1e-8,
            ("flow_mass_phase_comp", ("Vap", "H2O")): 1,
            ("temperature", None): 273.15 + 100,
            ("pressure", None): 1e5,
        }
        self.regression_solution = {
            ("dens_mass_phase", "Liq"): 958.3,
            ("dens_mass_phase", "Vap"): 0.5807,
            ("flow_vol_phase", "Liq"): 1.044e-11,
            ("flow_vol_phase", "Vap"): 1.722,
            ("flow_mol_phase_comp", ("Liq", "H2O")): 5.551e-7,
            ("flow_mol_phase_comp", ("Vap", "H2O")): 55.51,
            ("mole_frac_phase_comp", ("Liq", "H2O")): 1e-8,
            ("mole_frac_phase_comp", ("Vap", "H2O")): 1,
            ("enth_mass_phase", "Liq"): 4.190e5,
            ("enth_mass_phase", "Vap"): 2.676e6,
            ("dh_vap_mass", None): 2.257e6,
            ("cp_mass_phase", "Liq"): 4.215e3,
            ("cp_mass_phase", "Vap"): 1.890e3,
        }
