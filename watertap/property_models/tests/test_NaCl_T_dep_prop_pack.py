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

import pytest
import watertap.property_models.NaCl_T_dep_prop_pack as props
from idaes.models.properties.tests.test_harness import (
    PropertyTestHarness as PropertyTestHarness_idaes,
)
from watertap.property_models.tests.property_test_harness import (
    PropertyTestHarness,
    PropertyRegressionTest,
)

# -----------------------------------------------------------------------------
@pytest.mark.unit
class TestNaClProperty_idaes(PropertyTestHarness_idaes):
    def configure(self):
        self.prop_pack = props.NaClParameterBlock
        self.param_args = {}
        self.prop_args = {}
        self.has_density_terms = False


class TestNaClProperty(PropertyTestHarness):
    def configure(self):
        self.prop_pack = props.NaClParameterBlock
        self.param_args = {}
        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 1e2,
        }
        self.stateblock_statistics = {
            "number_variables": 24,
            "number_total_constraints": 20,
            "number_unused_variables": 1,
            "default_degrees_of_freedom": 3,
        }  # 4 state vars, but pressure is not active
        self.default_solution = {
            ("mass_frac_phase_comp", ("Liq", "H2O")): 0.965,
            ("mass_frac_phase_comp", ("Liq", "NaCl")): 0.035,
            ("dens_mass_phase", "Liq"): 1021.2,
            ("flow_vol_phase", "Liq"): 9.790e-4,
            ("conc_mass_phase_comp", ("Liq", "H2O")): 985.5,
            ("conc_mass_phase_comp", ("Liq", "NaCl")): 35.74,
            ("flow_mol_phase_comp", ("Liq", "H2O")): 53.57,
            ("flow_mol_phase_comp", ("Liq", "NaCl")): 0.5989,
            ("mole_frac_phase_comp", ("Liq", "H2O")): 0.9889,
            ("mole_frac_phase_comp", ("Liq", "NaCl")): 1.106e-2,
            ("molality_phase_comp", ("Liq", "NaCl")): 0.6206,
            ("osm_coeff", None): 0.9294,
            ("pressure_osm_phase", "Liq"): 2.851e6,
            ("diffus_phase_comp", ("Liq", "NaCl")): 1.517e-9,
            ("visc_d_phase", "Liq"): 0.000967,
            ("enth_mass_phase", "Liq"): 101.09e3,
            ("cp_mass_phase", "Liq"): 4000,
            ("solubility_comp", "NaCl"): 0.265,
            ("therm_cond_phase", "Liq"): 0.6015,
            ("pressure_sat", None): 0.00293e6,
        }


@pytest.mark.component
class TestNaClPropertySolution_1(PropertyRegressionTest):
    def configure(self):
        self.prop_pack = props.NaClParameterBlock
        self.param_args = {}

        self.solver = "ipopt"
        self.optarg = {"nlp_scaling_method": "user-scaling"}

        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 1e2,
        }
        self.state_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 0.95,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 0.05,
            ("temperature", None): 273.15 + 50,
            ("pressure", None): 5e5,
        }
        self.regression_solution = {
            ("mass_frac_phase_comp", ("Liq", "H2O")): 0.95,
            ("mass_frac_phase_comp", ("Liq", "NaCl")): 0.05,
            ("dens_mass_phase", "Liq"): 1021.7,
            ("flow_vol_phase", "Liq"): 9.787e-4,
            ("conc_mass_phase_comp", ("Liq", "H2O")): 970.6,
            ("conc_mass_phase_comp", ("Liq", "NaCl")): 51.08,
            ("flow_mol_phase_comp", ("Liq", "H2O")): 52.73,
            ("flow_mol_phase_comp", ("Liq", "NaCl")): 0.8556,
            ("mole_frac_phase_comp", ("Liq", "H2O")): 0.9840,
            ("mole_frac_phase_comp", ("Liq", "NaCl")): 1.597e-2,
            ("molality_phase_comp", ("Liq", "NaCl")): 0.9006,
            ("diffus_phase_comp", ("Liq", "NaCl")): 2.3e-9,
            ("visc_d_phase", "Liq"): 0.0005802,
            ("osm_coeff", None): 0.9322,
            ("pressure_osm_phase", "Liq"): 4.4579e6,
            ("enth_mass_phase", "Liq"): 197.9e3,
            ("cp_mass_phase", "Liq"): 3941.9,
            ("solubility_comp", "NaCl"): 0.2686,
            ("therm_cond_phase", "Liq"): 0.6392,
            ("pressure_sat", None): 0.01224e6,
        }


@pytest.mark.component
class TestNaClPropertySolution_2(PropertyRegressionTest):
    def configure(self):
        self.prop_pack = props.NaClParameterBlock
        self.param_args = {}

        self.solver = "ipopt"
        self.optarg = {"nlp_scaling_method": "user-scaling"}

        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 1e2,
        }
        self.state_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 0.85,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 0.15,
            ("temperature", None): 273.15 + 100,
            ("pressure", None): 5e5,
        }
        self.regression_solution = {
            ("mass_frac_phase_comp", ("Liq", "H2O")): 0.85,
            ("mass_frac_phase_comp", ("Liq", "NaCl")): 0.15,
            ("dens_mass_phase", "Liq"): 1063.3,
            ("flow_vol_phase", "Liq"): 0.0009404,
            ("conc_mass_phase_comp", ("Liq", "H2O")): 903.8,
            ("conc_mass_phase_comp", ("Liq", "NaCl")): 159.4,
            ("flow_mol_phase_comp", ("Liq", "H2O")): 47.18,
            ("flow_mol_phase_comp", ("Liq", "NaCl")): 2.566,
            ("mole_frac_phase_comp", ("Liq", "H2O")): 0.9484,
            ("mole_frac_phase_comp", ("Liq", "NaCl")): 0.05159,
            ("molality_phase_comp", ("Liq", "NaCl")): 3.019,
            ("diffus_phase_comp", ("Liq", "NaCl")): 4.577e-9,
            ("visc_d_phase", "Liq"): 0.0004139,
            ("osm_coeff", None): 1.045,
            ("pressure_osm_phase", "Liq"): 1.877e7,
            ("enth_mass_phase", "Liq"): 353.9e3,
            ("cp_mass_phase", "Liq"): 3603.7,
            ("solubility_comp", "NaCl"): 0.2799,
            ("therm_cond_phase", "Liq"): 0.6648,
            ("pressure_sat", None): 0.09011e6,
        }


@pytest.mark.component
class TestNaClPropertySolution_3(PropertyRegressionTest):
    def configure(self):
        self.prop_pack = props.NaClParameterBlock
        self.param_args = {}

        self.solver = "ipopt"
        self.optarg = {"nlp_scaling_method": "user-scaling"}

        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 1e2,
        }
        self.state_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 0.74,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 0.26,
            ("temperature", None): 273.15 + 150,
            ("pressure", None): 5e5,
        }
        self.regression_solution = {
            ("mass_frac_phase_comp", ("Liq", "H2O")): 0.74,
            ("mass_frac_phase_comp", ("Liq", "NaCl")): 0.26,
            ("dens_mass_phase", "Liq"): 1113.8,
            ("flow_vol_phase", "Liq"): 0.0008977,
            ("conc_mass_phase_comp", ("Liq", "H2O")): 824.2,
            ("conc_mass_phase_comp", ("Liq", "NaCl")): 289.6,
            ("flow_mol_phase_comp", ("Liq", "H2O")): 41.07,
            ("flow_mol_phase_comp", ("Liq", "NaCl")): 4.449,
            ("mole_frac_phase_comp", ("Liq", "H2O")): 0.9022,
            ("mole_frac_phase_comp", ("Liq", "NaCl")): 0.09772,
            ("molality_phase_comp", ("Liq", "NaCl")): 6.0121,
            ("diffus_phase_comp", ("Liq", "NaCl")): 4.8422e-9,
            ("visc_d_phase", "Liq"): 0.0003578,
            ("osm_coeff", None): 1.134,
            ("pressure_osm_phase", "Liq"): 4.399e7,
            ("enth_mass_phase", "Liq"): 475.3e3,
            ("cp_mass_phase", "Liq"): 3422.6,
            ("solubility_comp", "NaCl"): 0.2966,
            ("therm_cond_phase", "Liq"): 0.6627,
            ("pressure_sat", None): 0.3688e6,
        }
