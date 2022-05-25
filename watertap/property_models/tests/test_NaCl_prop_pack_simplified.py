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
import watertap.property_models.NaCl_prop_pack_simplified as props
from idaes.generic_models.properties.tests.test_harness import (
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
            "number_variables": 15,
            "number_total_constraints": 11,
            "number_unused_variables": 2,  # pressure is unused
            "default_degrees_of_freedom": 2,
        }  # 4 state vars, but pressure is not active
        self.default_solution = {
            ("mass_frac_phase_comp", ("Liq", "H2O")): 0.9999,
            ("mass_frac_phase_comp", ("Liq", "NaCl")): 1e-4,
        }