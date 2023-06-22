#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import watertap.property_models.selective_oil_permeation_prop_pack as props
from watertap.property_models.tests.property_test_harness import (
    PropertyTestHarness,
    PropertyRegressionTest,
)

# -----------------------------------------------------------------------------
class TestSOPProperty(PropertyTestHarness):
    def configure(self):
        self.prop_pack = props.SopParameterBlock
        self.param_args = {}
        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1,
            ("flow_mass_phase_comp", ("Liq", "oil")): 1e2,
        }
        self.stateblock_statistics = {
            "number_variables": 16,
            "number_total_constraints": 12,
            "number_unused_variables": 2,
            "default_degrees_of_freedom": 2,
        }
        self.default_solution = {
            ("visc_d_phase_comp", ("Liq", "H2O")): 1e-3,
            ("visc_d_phase_comp", ("Liq", "oil")): 3.56e-3,
            ("dens_mass_phase_comp", ("Liq", "H2O")): 1e3,
            ("dens_mass_phase_comp", ("Liq", "oil")): 780,
            ("mass_frac_phase_comp", ("Liq", "H2O")): 0.990099,
            ("mass_frac_phase_comp", ("Liq", "oil")): 9.90099e-3,
            ("flow_vol_phase_comp", ("Liq", "H2O")): 1e-3,
            ("flow_vol_phase_comp", ("Liq", "oil")): 1.28205e-5,
            ("flow_vol_phase", "Liq"): 1.01282e-3,
            ("flow_vol", None): 1.01282e-3,
            ("flow_mass_phase", "Liq"): 1.01,
            ("dens_mass_phase", "Liq"): 997.215,
            ("vol_frac_phase_comp", ("Liq", "H2O")): 0.98734,
            ("vol_frac_phase_comp", ("Liq", "oil")): 0.012658,
            ("conc_mass_phase_comp", ("Liq", "H2O")): 987.342,
            ("conc_mass_phase_comp", ("Liq", "oil")): 9.87342,
        }


class TestSOPPropertySolution(PropertyRegressionTest):
    def configure(self):
        self.prop_pack = props.SopParameterBlock
        self.param_args = {}
        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1,
            ("flow_mass_phase_comp", ("Liq", "oil")): 1,
        }
        self.state_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 0.5,
            ("flow_mass_phase_comp", ("Liq", "oil")): 0.5,
            ("temperature", None): 298.15,
            ("pressure", None): 101325,
        }
        self.regression_solution = {
            ("visc_d_phase_comp", ("Liq", "H2O")): 1e-3,
            ("visc_d_phase_comp", ("Liq", "oil")): 3.56e-3,
            ("dens_mass_phase_comp", ("Liq", "H2O")): 1e3,
            ("dens_mass_phase_comp", ("Liq", "oil")): 780,
            ("mass_frac_phase_comp", ("Liq", "H2O")): 0.5,
            ("mass_frac_phase_comp", ("Liq", "oil")): 0.5,
            ("flow_vol_phase_comp", ("Liq", "H2O")): 5e-4,
            ("flow_vol_phase_comp", ("Liq", "oil")): 6.41026e-4,
            ("flow_vol_phase", "Liq"): 1.141025641025641e-3,
            ("flow_vol", None): 1.141025641025641e-3,
            ("flow_mass_phase", "Liq"): 1,
            ("dens_mass_phase", "Liq"): 876.404,
            ("vol_frac_phase_comp", ("Liq", "H2O")): 0.4382,
            ("vol_frac_phase_comp", ("Liq", "oil")): 0.5618,
            ("conc_mass_phase_comp", ("Liq", "H2O")): 438.2022,
            ("conc_mass_phase_comp", ("Liq", "oil")): 438.2022,
        }
