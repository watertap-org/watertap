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
import watertap.property_models.NDMA_prop_pack as props
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
class TestNDMAroperty_idaes(PropertyTestHarness_idaes):
    def configure(self):
        self.prop_pack = props.NDMAParameterBlock
        self.param_args = {}
        self.prop_args = {}
        self.has_density_terms = False


class TestNDMAProperty(PropertyTestHarness):
    def configure(self):
        self.prop_pack = props.NDMAParameterBlock
        self.param_args = {}
        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1,
            ("flow_mass_phase_comp", ("Liq", "NDMA")): 1e8,
        }
        self.stateblock_statistics = {
            "number_variables": 15,
            "number_total_constraints": 11,
            "number_unused_variables": 2,  # pressure is unused
            "default_degrees_of_freedom": 2,
        }  # 3 state vars, but pressure is not active
        self.default_solution = {
            ("mass_frac_phase_comp", ("Liq", "H2O")): 0.999999926,
            ("mass_frac_phase_comp", ("Liq", "NDMA")): 74e-9,
            ("dens_mass_phase", "Liq"): 997,
            ("flow_vol_phase", "Liq"): 1.003e-3,
            ("conc_mass_phase_comp", ("Liq", "H2O")): 997,
            ("conc_mass_phase_comp", ("Liq", "NDMA")): 7.3778e-5,
            ("flow_mol_phase_comp", ("Liq", "H2O")): 55.5084,
            ("flow_mol_phase_comp", ("Liq", "NDMA")): 9.9889e-7,
            ("mole_frac_phase_comp", ("Liq", "H2O")): 0.999999982,
            ("mole_frac_phase_comp", ("Liq", "NDMA")): 1.7995e-8,
            ("molality_phase_comp", ("Liq", "NDMA")): 9.9889e-7,
        }


@pytest.mark.component
class TestNDMAPropertySolution_1(PropertyRegressionTest):
    def configure(self):
        self.prop_pack = props.NDMAParameterBlock
        self.param_args = {}

        self.solver = "ipopt"
        self.optarg = {"nlp_scaling_method": "user-scaling"}

        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1,
            ("flow_mass_phase_comp", ("Liq", "NDMA")): 1e2,
        }
        self.state_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 0.95,
            ("flow_mass_phase_comp", ("Liq", "NDMA")): 0.05,
            ("temperature", None): 273.15 + 25,
            ("pressure", None): 1e5,
        }
        self.regression_solution = {
            ("mass_frac_phase_comp", ("Liq", "H2O")): 0.95,
            ("mass_frac_phase_comp", ("Liq", "NDMA")): 0.05,
            ("dens_mass_phase", "Liq"): 997,
            ("flow_vol_phase", "Liq"): 1.003e-3,
            ("conc_mass_phase_comp", ("Liq", "H2O")): 947.15,
            ("conc_mass_phase_comp", ("Liq", "NDMA")): 49.85,
            ("flow_mol_phase_comp", ("Liq", "H2O")): 52.73,
            ("flow_mol_phase_comp", ("Liq", "NDMA")): 0.6749,
            ("mole_frac_phase_comp", ("Liq", "H2O")): 0.9874,
            ("mole_frac_phase_comp", ("Liq", "NDMA")): 1.2637e-2,
            ("molality_phase_comp", ("Liq", "NDMA")): 0.7105,
        }


@pytest.mark.component
class TestNDMACalculateState_1(PropertyCalculateStateTest):
    def configure(self):
        self.prop_pack = props.NDMAParameterBlock
        self.param_args = {}

        self.solver = "ipopt"
        self.optarg = {"nlp_scaling_method": "user-scaling"}

        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1e-1,
            ("flow_mass_phase_comp", ("Liq", "NDMA")): 1e1,
        }
        self.var_args = {
            ("flow_vol_phase", "Liq"): 2e-2,
            ("mass_frac_phase_comp", ("Liq", "NDMA")): 0.05,
            ("temperature", None): 273.15 + 25,
            ("pressure", None): 1e5,
        }
        self.state_solution = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 18.943,
            ("flow_mass_phase_comp", ("Liq", "NDMA")): 0.997,
        }
