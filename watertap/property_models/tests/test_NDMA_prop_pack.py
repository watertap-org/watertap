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
            ("mass_frac_phase_comp", ("Liq", "NDMA")): 0.035,
            ("dens_mass_phase", "Liq"): 1021.5,
            ("flow_vol_phase", "Liq"): 9.790e-4,
            ("conc_mass_phase_comp", ("Liq", "H2O")): 985.7,
            ("conc_mass_phase_comp", ("Liq", "NDMA")): 35.75,
            ("flow_mol_phase_comp", ("Liq", "H2O")): 100,
            ("flow_mol_phase_comp", ("Liq", "NDMA")): 0.5989,
            ("mole_frac_phase_comp", ("Liq", "H2O")): 0.9889,
            ("mole_frac_phase_comp", ("Liq", "NDMA")): 1.106e-2,
            ("molality_comp", "NDMA"): 0.6206,
        }