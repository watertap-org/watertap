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
from pyomo.environ import ConcreteModel
from idaes.core import FlowsheetBlock
import watertap.property_models.seawater_prop_pack as props
from idaes.models.properties.tests.test_harness import (
    PropertyTestHarness as PropertyTestHarness_idaes,
)
from watertap.property_models.tests.property_test_harness import (
    PropertyTestHarness,
    PropertyRegressionTest,
    PropertyCalculateStateTest,
)


class TestSeawaterProperty_idaes(PropertyTestHarness_idaes):
    def configure(self):
        self.prop_pack = props.SeawaterParameterBlock
        self.param_args = {}
        self.prop_args = {}
        self.has_density_terms = True


class TestSeawaterProperty(PropertyTestHarness):
    def configure(self):
        self.prop_pack = props.SeawaterParameterBlock
        self.param_args = {}
        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1,
            ("flow_mass_phase_comp", ("Liq", "TDS")): 1e2,
        }
        self.stateblock_statistics = {
            "number_variables": 26,
            "number_total_constraints": 22,
            "number_unused_variables": 0,
            "default_degrees_of_freedom": 4,
        }  # 4 state vars, but pressure is not active
        self.default_solution = {
            ("mass_frac_phase_comp", ("Liq", "H2O")): 0.965,
            ("mass_frac_phase_comp", ("Liq", "TDS")): 0.035,
            ("dens_mass_phase", "Liq"): 1023.5,
            ("dens_mass_solvent", None): 996.9,
            ("flow_vol_phase", "Liq"): 9.770e-4,
            ("conc_mass_phase_comp", ("Liq", "H2O")): 987.7,
            ("conc_mass_phase_comp", ("Liq", "TDS")): 35.82,
            ("flow_mol_phase_comp", ("Liq", "H2O")): 53.57,
            ("flow_mol_phase_comp", ("Liq", "TDS")): 1.1145,
            ("mole_frac_phase_comp", ("Liq", "H2O")): 0.9796,
            ("mole_frac_phase_comp", ("Liq", "TDS")): 2.038e-2,
            ("molality_phase_comp", ("Liq", "TDS")): 1.155,
            ("visc_d_phase", "Liq"): 9.588e-4,
            ("osm_coeff", None): 0.9068,
            ("pressure_osm_phase", "Liq"): 2.588e6,
            ("enth_mass_phase", "Liq"): 9.9766e4,
            ("pressure_sat", None): 3111,
            ("cp_mass_phase", "Liq"): 4001,
            ("therm_cond_phase", "Liq"): 0.6086,
            ("dh_vap_mass", None): 2.356e6,
            ("diffus_phase_comp", ("Liq", "TDS")): 1.471e-9,
            ("boiling_point_elevation_phase", "Liq"): 0.3093,
        }


@pytest.mark.component
class TestSeawaterPropertySolution_1(PropertyRegressionTest):
    def configure(self):
        self.prop_pack = props.SeawaterParameterBlock
        self.param_args = {}

        self.solver = "ipopt"
        self.optarg = {"nlp_scaling_method": "user-scaling"}

        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1,
            ("flow_mass_phase_comp", ("Liq", "TDS")): 1e2,
        }
        self.state_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 0.99,
            ("flow_mass_phase_comp", ("Liq", "TDS")): 0.01,
            ("temperature", None): 273.15 + 50,
            ("pressure", None): 2e5,
        }
        self.regression_solution = {
            ("mass_frac_phase_comp", ("Liq", "H2O")): 0.99,
            ("mass_frac_phase_comp", ("Liq", "TDS")): 0.01,
            ("dens_mass_phase", "Liq"): 995.4,
            ("dens_mass_solvent", None): 988.0,
            ("flow_vol_phase", "Liq"): 1.005e-3,
            ("conc_mass_phase_comp", ("Liq", "H2O")): 985.5,
            ("conc_mass_phase_comp", ("Liq", "TDS")): 9.954,
            ("flow_mol_phase_comp", ("Liq", "H2O")): 54.95,
            ("flow_mol_phase_comp", ("Liq", "TDS")): 0.3184,
            ("mole_frac_phase_comp", ("Liq", "H2O")): 0.9942,
            ("mole_frac_phase_comp", ("Liq", "TDS")): 5.761e-3,
            ("molality_phase_comp", ("Liq", "TDS")): 0.3216,
            ("visc_d_phase", "Liq"): 5.596e-4,
            ("osm_coeff", None): 0.9029,
            ("pressure_osm_phase", "Liq"): 7.710e5,
            ("enth_mass_phase", "Liq"): 2.06690e5,
            ("energy_density_phase", "Liq"): 2.0555e8,
            ("pressure_sat", None): 1.229e4,
            ("cp_mass_phase", "Liq"): 4.130e3,
            ("therm_cond_phase", "Liq"): 0.6400,
            ("dh_vap_mass", None): 2.358e6,
            ("diffus_phase_comp", ("Liq", "TDS")): 1.493e-9,
            ("boiling_point_elevation_phase", "Liq"): 0.0989,
        }


@pytest.mark.component
class TestSeawaterPropertySolution_2(PropertyRegressionTest):
    def configure(self):
        self.prop_pack = props.SeawaterParameterBlock
        self.param_args = {}

        self.solver = "ipopt"
        self.optarg = {"nlp_scaling_method": "user-scaling"}

        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1,
            ("flow_mass_phase_comp", ("Liq", "TDS")): 1e2,
        }
        self.state_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 0.95,
            ("flow_mass_phase_comp", ("Liq", "TDS")): 0.05,
            ("temperature", None): 273.15 + 10,
            ("pressure", None): 100e5,
        }
        self.regression_solution = {
            ("mass_frac_phase_comp", ("Liq", "H2O")): 0.95,
            ("mass_frac_phase_comp", ("Liq", "TDS")): 0.05,
            ("dens_mass_phase", "Liq"): 1.039e3,
            ("dens_mass_solvent", None): 999.5,
            ("flow_vol_phase", "Liq"): 9.628e-4,
            ("conc_mass_phase_comp", ("Liq", "H2O")): 986.8,
            ("conc_mass_phase_comp", ("Liq", "TDS")): 51.93,
            ("flow_mol_phase_comp", ("Liq", "H2O")): 52.73,
            ("flow_mol_phase_comp", ("Liq", "TDS")): 1.592,
            ("mole_frac_phase_comp", ("Liq", "H2O")): 0.9707,
            ("mole_frac_phase_comp", ("Liq", "TDS")): 2.931e-2,
            ("molality_phase_comp", ("Liq", "TDS")): 1.676,
            ("visc_d_phase", "Liq"): 1.443e-3,
            ("osm_coeff", None): 0.9106,
            ("pressure_osm_phase", "Liq"): 3.591e6,
            ("enth_mass_phase", "Liq"): 4.8008e4,
            ("energy_density_phase", "Liq"): 3.9865e7,
            ("pressure_sat", None): 1.194e3,
            ("cp_mass_phase", "Liq"): 3.916e3,
            ("therm_cond_phase", "Liq"): 0.5854,
            ("dh_vap_mass", None): 2.353e6,
            ("diffus_phase_comp", ("Liq", "TDS")): 1.471e-9,
            ("boiling_point_elevation_phase", "Liq"): 0.4069,
        }


# Parameter values from:
# https://web.mit.edu/seawater/2017_MIT_Seawater_Property_Tables_r2b_2023c.pdf
# T = 10C, salinity = 10 g/kg, Q = 1 L/s
@pytest.mark.component
class TestSeawaterPropertyTemp10Salinity10(PropertyRegressionTest):
    def configure(self):
        self.prop_pack = props.SeawaterParameterBlock
        self.param_args = {}

        self.solver = "ipopt"
        self.optarg = {"nlp_scaling_method": "user-scaling"}

        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1,
            ("flow_mass_phase_comp", ("Liq", "TDS")): 10,
        }
        self.state_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 0.9972722,
            ("flow_mass_phase_comp", ("Liq", "TDS")): 0.010073,
            ("temperature", None): 273.15 + 10,
            ("pressure", None): 1e5,
        }
        self.regression_solution = {
            ("mass_frac_phase_comp", ("Liq", "H2O")): 0.99,
            ("mass_frac_phase_comp", ("Liq", "TDS")): 0.01,
            ("dens_mass_phase", "Liq"): 1007.3457,
            ("dens_mass_solvent", None): 999.50934,
            ("flow_vol_phase", "Liq"): 0.001,
            ("conc_mass_phase_comp", ("Liq", "H2O")): 997.27225,
            ("conc_mass_phase_comp", ("Liq", "TDS")): 10.07345,
            ("flow_mol_phase_comp", ("Liq", "H2O")): 55.35702,
            ("flow_mol_phase_comp", ("Liq", "TDS")): 0.32077,
            ("mole_frac_phase_comp", ("Liq", "H2O")): 0.99423,
            ("mole_frac_phase_comp", ("Liq", "TDS")): 0.00576,
            ("molality_phase_comp", ("Liq", "TDS")): 0.32164,
            ("visc_d_phase", "Liq"): 0.00132956,  # 0.0133 Pa-s
            ("osm_coeff", None): 0.89768,
            ("pressure_osm_phase", "Liq"): 679429.93808,  # 0.679 MPa
            ("enth_mass_phase", "Liq"): 41590.0253,  # 41.6 kJ/kg
            ("energy_density_phase", "Liq"): 41795533.49664,
            ("pressure_sat", None): 1222.13185,  # 1.222 kPa
            ("cp_mass_phase", "Liq"): 4136.70887,  # 4136.7 J/kg-K
            ("therm_cond_phase", "Liq"): 0.58768,  # 0.588 W/m-K
            ("dh_vap_mass", None): 2452555.1844,  # 2452.5 kJ/kg
            ("diffus_phase_comp", ("Liq", "TDS")): 1.492e-09,
            ("boiling_point_elevation_phase", "Liq"): 0.07309,  # 0.073 K
        }


# T = 30C, salinity = 30 g/kg, Q = 1 L/s
@pytest.mark.component
class TestSeawaterPropertyTemp30Salinity30(PropertyRegressionTest):
    def configure(self):
        self.prop_pack = props.SeawaterParameterBlock
        self.param_args = {}

        self.solver = "ipopt"
        self.optarg = {"nlp_scaling_method": "user-scaling"}

        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1,
            ("flow_mass_phase_comp", ("Liq", "TDS")): 10,
        }
        self.state_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 0.987677512,
            ("flow_mass_phase_comp", ("Liq", "TDS")): 0.030546727,
            ("temperature", None): 273.15 + 30,
            ("pressure", None): 1e5,
        }
        self.regression_solution = {
            ("mass_frac_phase_comp", ("Liq", "H2O")): 0.97,
            ("mass_frac_phase_comp", ("Liq", "TDS")): 0.03,
            ("dens_mass_phase", "Liq"): 1018.22423,
            ("dens_mass_solvent", None): 995.53714,
            ("flow_vol_phase", "Liq"): 0.001,
            ("conc_mass_phase_comp", ("Liq", "H2O")): 987.67751,
            ("conc_mass_phase_comp", ("Liq", "TDS")): 30.54672,
            ("flow_mol_phase_comp", ("Liq", "H2O")): 54.82443,
            ("flow_mol_phase_comp", ("Liq", "TDS")): 0.9727,
            ("mole_frac_phase_comp", ("Liq", "H2O")): 0.98256,
            ("mole_frac_phase_comp", ("Liq", "TDS")): 0.01743,
            ("molality_phase_comp", ("Liq", "TDS")): 0.98484,
            ("visc_d_phase", "Liq"): 0.00085086,  # 0.000851 Pa-s
            ("osm_coeff", None): 0.90569,
            ("pressure_osm_phase", "Liq"): 2238198.31526,  # 2.238 MPa
            ("enth_mass_phase", "Liq"): 120597.84791,  # 120.6 kJ/kg
            ("energy_density_phase", "Liq"): 122695651.94676,
            ("pressure_sat", None): 4180.36772,  # 4.180 kPa
            ("cp_mass_phase", "Liq"): 4027.48595,  # 4027.5 J/kg-K
            ("therm_cond_phase", "Liq"): 0.61563,  # 0.616 W/m-K
            ("dh_vap_mass", None): 2357037.33712,  # 2356.9 kJ/kg
            ("diffus_phase_comp", ("Liq", "TDS")): 1.473e-09,
            ("boiling_point_elevation_phase", "Liq"): 0.27175,  # 0.272 K
        }


# T = 50C, salinity = 50 g/kg, Q = 1 L/s
@pytest.mark.component
class TestSeawaterPropertyTemp50Salinity50(PropertyRegressionTest):
    def configure(self):
        self.prop_pack = props.SeawaterParameterBlock
        self.param_args = {}

        self.solver = "ipopt"
        self.optarg = {"nlp_scaling_method": "user-scaling"}

        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1,
            ("flow_mass_phase_comp", ("Liq", "TDS")): 10,
        }
        self.state_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 0.973797107,
            ("flow_mass_phase_comp", ("Liq", "TDS")): 0.051252479,
            ("temperature", None): 273.15 + 50,
            ("pressure", None): 1e5,
        }
        self.regression_solution = {
            ("mass_frac_phase_comp", ("Liq", "H2O")): 0.95,
            ("mass_frac_phase_comp", ("Liq", "TDS")): 0.05,
            ("dens_mass_phase", "Liq"): 1025.04958,
            ("dens_mass_solvent", None): 988.04718,
            ("flow_vol_phase", "Liq"): 0.001,
            ("conc_mass_phase_comp", ("Liq", "H2O")): 973.7971,
            ("conc_mass_phase_comp", ("Liq", "TDS")): 51.25247,
            ("flow_mol_phase_comp", ("Liq", "H2O")): 54.05395,
            ("flow_mol_phase_comp", ("Liq", "TDS")): 1.63204,
            ("mole_frac_phase_comp", ("Liq", "H2O")): 0.97069,
            ("mole_frac_phase_comp", ("Liq", "TDS")): 0.0293,
            ("molality_phase_comp", ("Liq", "TDS")): 1.67596,
            ("visc_d_phase", "Liq"): 0.00061694,  # 0.000617 Pa-s
            ("osm_coeff", None): 0.91724,
            ("pressure_osm_phase", "Liq"): 4080987.77845,  # 4.081 MPa
            ("enth_mass_phase", "Liq"): 195859.02578,  # 195.9 kJ/kg
            ("energy_density_phase", "Liq"): 200665213.4249,
            ("pressure_sat", None): 12008.61958,  # 12.009 kPa
            ("cp_mass_phase", "Liq"): 3940.51113,  # 3940.5 J/kg-K
            ("therm_cond_phase", "Liq"): 0.63805,  # 0.638 W/m-K
            ("dh_vap_mass", None): 2262972.85312,  # 2262.9 kJ/kg
            ("diffus_phase_comp", ("Liq", "TDS")): 1.47e-09,
            ("boiling_point_elevation_phase", "Liq"): 0.55617,  # 0.556 K
        }


# T = 70C, salinity = 70 g/kg, Q = 1 L/s
@pytest.mark.component
class TestSeawaterPropertyTemp70Salinity70(PropertyRegressionTest):
    def configure(self):
        self.prop_pack = props.SeawaterParameterBlock
        self.param_args = {}

        self.solver = "ipopt"
        self.optarg = {"nlp_scaling_method": "user-scaling"}

        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1,
            ("flow_mass_phase_comp", ("Liq", "TDS")): 10,
        }
        self.state_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 0.95708085,
            ("flow_mass_phase_comp", ("Liq", "TDS")): 0.07203834,
            ("temperature", None): 273.15 + 70,
            ("pressure", None): 1e5,
        }
        self.regression_solution = {
            ("mass_frac_phase_comp", ("Liq", "H2O")): 0.93,
            ("mass_frac_phase_comp", ("Liq", "TDS")): 0.07,
            ("dens_mass_phase", "Liq"): 1029.1192,
            ("dens_mass_solvent", None): 977.76708,
            ("flow_vol_phase", "Liq"): 0.001,
            ("conc_mass_phase_comp", ("Liq", "H2O")): 957.08085,
            ("conc_mass_phase_comp", ("Liq", "TDS")): 72.03834,
            ("flow_mol_phase_comp", ("Liq", "H2O")): 53.12606,
            ("flow_mol_phase_comp", ("Liq", "TDS")): 2.29393,
            ("mole_frac_phase_comp", ("Liq", "H2O")): 0.9586,
            ("mole_frac_phase_comp", ("Liq", "TDS")): 0.04139,
            ("molality_phase_comp", ("Liq", "TDS")): 2.3968,
            ("visc_d_phase", "Liq"): 0.00048369,  # 0.000484 Pa-s
            ("osm_coeff", None): 0.93293,
            ("pressure_osm_phase", "Liq"): 6237901.49812,  # 6.238 MPa
            ("enth_mass_phase", "Liq"): 268014.57187,  # 268.0 kJ/kg
            ("energy_density_phase", "Liq"): 275718942.07197,
            ("pressure_sat", None): 29912.03673,  # 29.912 kPa
            ("cp_mass_phase", "Liq"): 3862.34273,  # 3862.4 J/kg-K
            ("therm_cond_phase", "Liq"): 0.65529,  # 0.655 W/m-K
            ("dh_vap_mass", None): 2169879.46248,  # 2169.8 kJ/kg
            ("diffus_phase_comp", ("Liq", "TDS")): 1.479e-09,
            ("boiling_point_elevation_phase", "Liq"): 0.94374,  # 0.944 K
        }


# T = 90C, salinity = 90 g/kg
@pytest.mark.component
class TestSeawaterPropertyTemp90Salinity90(PropertyRegressionTest):
    def configure(self):
        self.prop_pack = props.SeawaterParameterBlock
        self.param_args = {}

        self.solver = "ipopt"
        self.optarg = {"nlp_scaling_method": "user-scaling"}

        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1,
            ("flow_mass_phase_comp", ("Liq", "TDS")): 10,
        }
        self.state_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 0.93860506,
            ("flow_mass_phase_comp", ("Liq", "TDS")): 0.092829072,
            ("temperature", None): 273.15 + 90,
            ("pressure", None): 1e5,
        }
        self.regression_solution = {
            ("mass_frac_phase_comp", ("Liq", "H2O")): 0.90999,
            ("mass_frac_phase_comp", ("Liq", "TDS")): 0.09,
            ("dens_mass_phase", "Liq"): 1031.43413,
            ("dens_mass_solvent", None): 965.24563,
            ("flow_vol_phase", "Liq"): 0.001,
            ("conc_mass_phase_comp", ("Liq", "H2O")): 938.60506,
            ("conc_mass_phase_comp", ("Liq", "TDS")): 92.82907,
            ("flow_mol_phase_comp", ("Liq", "H2O")): 52.10049,
            ("flow_mol_phase_comp", ("Liq", "TDS")): 2.95598,
            ("mole_frac_phase_comp", ("Liq", "H2O")): 0.94631,
            ("mole_frac_phase_comp", ("Liq", "TDS")): 0.05368,
            ("molality_phase_comp", ("Liq", "TDS")): 3.14933,
            ("visc_d_phase", "Liq"): 0.00039999,  # 0.0004 Pa-s
            ("osm_coeff", None): 0.9523,
            ("pressure_osm_phase", "Liq"): 8740862.18899,  # 8.741 MPa
            ("enth_mass_phase", "Liq"): 336459.24947,  # 336.5 kJ/kg
            ("energy_density_phase", "Liq"): 346935555.96049,
            ("pressure_sat", None): 66238.94125,  # 66.239 kPa
            ("cp_mass_phase", "Liq"): 3788.28386,  # 3788.3 J/kg-K
            ("therm_cond_phase", "Liq"): 0.66767,  # 0.668 W/m-K
            ("dh_vap_mass", None): 2077246.1356,  # 2077.1 kJ/kg
            ("diffus_phase_comp", ("Liq", "TDS")): 1.494e-09,
            ("boiling_point_elevation_phase", "Liq"): 1.45011,  # 1.450 K
        }


# T = 120C, salinity = 120 g/kg, Q = 1 L/s
@pytest.mark.component
class TestSeawaterPropertyTemp120Salinity120(PropertyRegressionTest):
    def configure(self):
        self.prop_pack = props.SeawaterParameterBlock
        self.param_args = {}

        self.solver = "ipopt"
        self.optarg = {"nlp_scaling_method": "user-scaling"}

        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1,
            ("flow_mass_phase_comp", ("Liq", "TDS")): 10,
        }
        self.state_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 0.90910769,
            ("flow_mass_phase_comp", ("Liq", "TDS")): 0.12396923,
            ("temperature", None): 273.15 + 120,
            ("pressure", None): 1e5,
        }
        self.regression_solution = {
            ("mass_frac_phase_comp", ("Liq", "H2O")): 0.88,
            ("mass_frac_phase_comp", ("Liq", "TDS")): 0.12,
            ("dens_mass_phase", "Liq"): 1033.07692,
            ("dens_mass_solvent", None): 943.02132,
            ("flow_vol_phase", "Liq"): 0.001,
            ("conc_mass_phase_comp", ("Liq", "H2O")): 909.10769,
            ("conc_mass_phase_comp", ("Liq", "TDS")): 123.96923,
            ("flow_mol_phase_comp", ("Liq", "H2O")): 50.46314,
            ("flow_mol_phase_comp", ("Liq", "TDS")): 3.94758,
            ("mole_frac_phase_comp", ("Liq", "H2O")): 0.92744,
            ("mole_frac_phase_comp", ("Liq", "TDS")): 0.07255,
            ("molality_phase_comp", ("Liq", "TDS")): 4.34226,
            ("visc_d_phase", "Liq"): 0.00032278,  # 0.000323 Pa-s
            ("osm_coeff", None): 0.98462,
            ("pressure_osm_phase", "Liq"): 13179594.05754,  # 13.180 MPa
            ("enth_mass_phase", "Liq"): 427777.48946,  # 427.8 kJ/kg
            ("energy_density_phase", "Liq"): 441827053.0114,
            ("pressure_sat", None): 182600.72284,  # 182.601 kPa
            ("cp_mass_phase", "Liq"): 3688.73469,  # 3688.8 J/kg-K
            ("therm_cond_phase", "Liq"): 0.67781,  # 0.678 W/m-K
            ("dh_vap_mass", None): 1937991.723,  # 1937.9 kJ/kg
            ("diffus_phase_comp", ("Liq", "TDS")): 1.524e-09,
            ("boiling_point_elevation_phase", "Liq"): 2.4623,  # 2.462 K
        }


@pytest.mark.component
class TestSeawaterCalculateState_1(PropertyCalculateStateTest):
    def configure(self):
        self.prop_pack = props.SeawaterParameterBlock
        self.param_args = {}

        self.solver = "ipopt"
        self.optarg = {"nlp_scaling_method": "user-scaling"}

        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1e-1,
            ("flow_mass_phase_comp", ("Liq", "TDS")): 1e1,
        }
        self.var_args = {
            ("flow_vol_phase", "Liq"): 2e-2,
            ("mass_frac_phase_comp", ("Liq", "TDS")): 0.05,
            ("temperature", None): 273.15 + 25,
            ("pressure", None): 5e5,
        }
        self.state_solution = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 19.66,
            ("flow_mass_phase_comp", ("Liq", "TDS")): 1.035,
        }


@pytest.mark.component
class TestSeawaterCalculateState_2(PropertyCalculateStateTest):
    def configure(self):
        self.prop_pack = props.SeawaterParameterBlock
        self.param_args = {}

        self.solver = "ipopt"
        self.optarg = {"nlp_scaling_method": "user-scaling"}

        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1,
            ("flow_mass_phase_comp", ("Liq", "TDS")): 1e2,
        }
        self.var_args = {
            ("flow_vol_phase", "Liq"): 1e-3,
            ("pressure_osm_phase", "Liq"): 100e5,
            ("temperature", None): 273.15 + 25,
            ("pressure", None): 5e5,
        }
        self.state_solution = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 0.9604,
            ("flow_mass_phase_comp", ("Liq", "TDS")): 0.1231,
        }


@pytest.mark.component
class TestSeawaterCalculateState_3(PropertyCalculateStateTest):
    def configure(self):
        self.prop_pack = props.SeawaterParameterBlock
        self.param_args = {}

        self.solver = "ipopt"
        self.optarg = {"nlp_scaling_method": "user-scaling"}

        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1,
            ("flow_mass_phase_comp", ("Liq", "TDS")): 1e2,
        }
        self.var_args = {
            ("flow_vol_phase", "Liq"): 1e-3,
            ("mass_frac_phase_comp", ("Liq", "TDS")): 0.06,
            ("pressure_sat", None): 3e4,
            ("pressure", None): 5e5,
        }
        self.state_solution = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 0.9605,
            ("flow_mass_phase_comp", ("Liq", "TDS")): 0.0613,
            ("temperature", None): 343.05,
        }


@pytest.mark.component
class TestSeawaterCalculateState_4(PropertyCalculateStateTest):
    def configure(self):
        self.prop_pack = props.SeawaterParameterBlock
        self.param_args = {}

        self.solver = "ipopt"
        self.optarg = {"nlp_scaling_method": "user-scaling"}

        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1e-4,
            ("flow_mass_phase_comp", ("Liq", "TDS")): 1e1,
            ("enth_mass_phase", ("Liq", "TDS")): 1e-4,
        }
        self.var_args = {
            ("flow_vol_phase", "Liq"): 10,
            ("mass_frac_phase_comp", ("Liq", "TDS")): 1e-5,
            ("temperature", None): 273.15 + 10,
            ("pressure", None): 2e6,
        }
        self.state_solution = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 9995.07,
            ("flow_mass_phase_comp", ("Liq", "TDS")): 0.09995,
            ("enth_mass_phase", "Liq"): 4.3945e4,
        }


@pytest.mark.component
class TestSeawaterCalculateState_5(PropertyCalculateStateTest):
    def configure(self):
        self.prop_pack = props.SeawaterParameterBlock
        self.param_args = {}

        self.solver = "ipopt"
        self.optarg = {"nlp_scaling_method": "user-scaling"}

        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1e-3,
            ("flow_mass_phase_comp", ("Liq", "TDS")): 1e-3,
            ("enth_mass_phase", ("Liq", "TDS")): 1e-5,
        }
        self.var_args = {
            ("flow_vol_phase", "Liq"): 10,
            ("mass_frac_phase_comp", ("Liq", "TDS")): 0.12,
            ("temperature", None): 273.15 + 120,
            ("pressure", None): 1.2e7,
        }
        self.state_solution = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 9091.08,
            ("flow_mass_phase_comp", ("Liq", "TDS")): 1239.69,
            ("enth_mass_phase", "Liq"): 4.36167e5,
        }


@pytest.mark.unit
def test_list_and_print_properties():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.props = props.SeawaterParameterBlock()

    m.fs.props.list_properties()
    m.fs.props.print_properties()
