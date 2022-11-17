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
import watertap.property_models.cryst_prop_pack as props
from pyomo.environ import ConcreteModel, SolverFactory, TerminationCondition
from idaes.core import FlowsheetBlock, ControlVolume0DBlock
from idaes.models.properties.tests.test_harness import (
    PropertyTestHarness as PropertyTestHarness_idaes,
)
from watertap.property_models.tests.property_test_harness import (
    PropertyTestHarness,
    PropertyRegressionTest,
    PropertyCalculateStateTest,
)


# -----------------------------------------------------------------------------


class TestNaClProperty_idaes(PropertyTestHarness_idaes):
    def configure(self):
        self.prop_pack = props.NaClParameterBlock
        self.param_args = {}
        self.prop_args = {}
        self.has_density_terms = False


class TestDefaultNaClwaterProperty:

    # Create block and stream for running default tests
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.NaClParameterBlock(
        heat_of_crystallization_model=props.HeatOfCrystallizationModel.constant
    )
    m.fs.stream = m.fs.properties.build_state_block([0], defined_state=True)

    m.fs.cv = ControlVolume0DBlock(
        dynamic=False, has_holdup=False, property_package=m.fs.properties
    )
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_material_balances()
    m.fs.cv.add_energy_balances()
    m.fs.cv.add_momentum_balances()

    # Create instance of PropertyTesthARNESS class and add attributes needed for tests
    xv = PropertyTestHarness()

    xv.stateblock_statistics = {
        "number_variables": 43,
        "number_total_constraints": 37,
        "number_unused_variables": 0,
        "default_degrees_of_freedom": 6,
    }

    xv.scaling_args = {
        ("flow_mass_phase_comp", ("Liq", "H2O")): 1e0,
        ("flow_mass_phase_comp", ("Liq", "NaCl")): 1e1,
        ("flow_mol_phase_comp", ("Sol", "NaCl")): 1e2,
        ("flow_mol_phase_comp", ("Vap", "H2O")): 1e3,
    }

    xv.default_solution = {
        ("solubility_mass_phase_comp", ("Liq", "NaCl")): 359.51,
        ("solubility_mass_frac_phase_comp", ("Liq", "NaCl")): 0.265,
        ("dens_mass_solvent", "Liq"): 996.89,
        ("dens_mass_solvent", "Vap"): 0.7363,
        ("dens_mass_solute", "Liq"): 3199.471,
        ("dens_mass_solute", "Sol"): 2115,
        ("dens_mass_phase", "Liq"): 1021.50,
        ("dh_vap_mass_solvent", None): 2441.80,
        ("cp_mass_solvent", "Liq"): 4186.52,
        ("cp_mass_solvent", "Vap"): 1864.52,
        ("cp_mass_solute", "Sol"): 864.15,  # cp_mass_solute liquid ignored for now
        ("cp_mass_phase", "Liq"): 4008.30,
        ("flow_vol_phase", "Liq"): 9.79e-4,
        ("flow_vol_phase", "Sol"): 0,
        ("flow_vol_phase", "Vap"): 0,
        ("pressure_sat", None): 2932.43,
        ("temperature_sat_solvent", None): 296.79,
        ("conc_mass_phase_comp", ("Liq", "NaCl")): 35.753,
        ("conc_mass_phase_comp", ("Liq", "H2O")): 985.753,
        ("enth_mass_solvent", "Liq"): 104.9212,
        ("enth_mass_solvent", "Vap"): 2546.7289,
        ("enth_mass_solute", "Sol"): 21.4739,
        ("enth_mass_phase", "Liq"): 101.0922,
        ("dh_crystallization_mass_comp", "NaCl"): -520,
        ("mass_frac_phase_comp", ("Liq", "H2O")): 0.965,
        ("mass_frac_phase_comp", ("Liq", "NaCl")): 0.035,
        ("mass_frac_phase_comp", ("Sol", "NaCl")): 1.0,
        ("mass_frac_phase_comp", ("Vap", "H2O")): 1.0,
        ("flow_mol_phase_comp", ("Liq", "H2O")): 53.57,
        ("flow_mol_phase_comp", ("Liq", "NaCl")): 0.5989,
        ("flow_mol_phase_comp", ("Sol", "NaCl")): 1.6318503332497833e-09,
        ("flow_mol_phase_comp", ("Vap", "H2O")): 1.6318424528638692e-07,
        ("mole_frac_phase_comp", ("Liq", "H2O")): 0.9889,
        ("mole_frac_phase_comp", ("Liq", "NaCl")): 0.01106,
        ("mole_frac_phase_comp", ("Sol", "NaCl")): 1.0,
        ("mole_frac_phase_comp", ("Vap", "H2O")): 1.0,
    }

    # Configure class
    xv.configure_class(m)

    @pytest.mark.unit
    def test_components_phases(self):
        self.xv.test_components_phases(self.m)

    @pytest.mark.unit
    def test_parameters(self):
        self.xv.test_parameters(self.m)

    @pytest.mark.unit
    def test_state_variables(self):
        self.xv.test_state_variables(self.m)

    @pytest.mark.unit
    def test_permanent_properties(self):
        self.xv.test_permanent_properties(self.m)

    @pytest.mark.unit
    def test_on_demand_properties(self):
        self.xv.test_on_demand_properties(self.m)

    @pytest.mark.unit
    def test_stateblock_statistics(self):
        self.xv.test_stateblock_statistics(self.m)

    @pytest.mark.unit
    def test_units_consistent(self):
        self.xv.test_units_consistent(self.m)

    @pytest.mark.unit
    def test_scaling(self):
        self.xv.test_scaling(self.m)

    @pytest.mark.component
    def test_default_initialization(self):
        self.xv.test_default_initialization(self.m)

    @pytest.mark.component
    def test_property_control_volume(self):
        self.xv.test_property_control_volume(self.m)


@pytest.mark.component
class TestNaClPropertySolution_1(PropertyRegressionTest):
    # Test pure liquid solution 1 - same solution as NaCl prop pack
    def configure(self):
        self.prop_pack = props.NaClParameterBlock
        self.param_args = {}

        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1e0,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 1e2,
            ("flow_mass_phase_comp", ("Vap", "H2O")): 1e2,
        }

        self.state_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 0.95,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 0.05,
            ("flow_mass_phase_comp", ("Sol", "NaCl")): 0.0,
            ("flow_mass_phase_comp", ("Vap", "H2O")): 0.0,
            ("temperature", None): 273.15 + 25,
            ("pressure", None): 50e5,
        }

        self.regression_solution = {
            ("mass_frac_phase_comp", ("Liq", "H2O")): 0.95,
            ("mass_frac_phase_comp", ("Liq", "NaCl")): 0.05,
            ("solubility_mass_phase_comp", ("Liq", "NaCl")): 359.50,
            ("solubility_mass_frac_phase_comp", ("Liq", "NaCl")): 0.265,
            ("dens_mass_solvent", "Liq"): 996.89,
            ("dens_mass_solvent", "Vap"): 36.335,
            ("dens_mass_phase", "Liq"): 1032.2,
            ("flow_vol_phase", "Liq"): 9.687e-4,
            ("conc_mass_phase_comp", ("Liq", "H2O")): 980.6,
            ("conc_mass_phase_comp", ("Liq", "NaCl")): 51.61,
            ("flow_mol_phase_comp", ("Liq", "H2O")): 52.73,
            ("flow_mol_phase_comp", ("Liq", "NaCl")): 0.8556,
            ("mole_frac_phase_comp", ("Liq", "H2O")): 0.9840,
            ("mole_frac_phase_comp", ("Liq", "NaCl")): 1.597e-2,
            ("cp_mass_solvent", "Liq"): 4186.52,
            ("cp_mass_phase", "Liq"): 3940.03,
            ("enth_mass_solvent", "Liq"): 104.92,
            ("enth_mass_phase", "Liq"): 99.45,
        }


@pytest.mark.component
class TestNaClPropertySolution_2(PropertyRegressionTest):
    # Test pure liquid solution 2 - same solution as NaCl prop pack
    def configure(self):
        self.prop_pack = props.NaClParameterBlock
        self.param_args = {}

        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1e0,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 1e2,
            ("flow_mass_phase_comp", ("Vap", "H2O")): 1e2,
        }

        self.state_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 0.74,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 0.26,
            ("flow_mass_phase_comp", ("Sol", "NaCl")): 0.0,
            ("flow_mass_phase_comp", ("Vap", "H2O")): 0.0,
            ("temperature", None): 273.15 + 25,
            ("pressure", None): 50e5,
        }

        self.regression_solution = {
            ("mass_frac_phase_comp", ("Liq", "H2O")): 0.74,
            ("mass_frac_phase_comp", ("Liq", "NaCl")): 0.26,
            ("solubility_mass_phase_comp", ("Liq", "NaCl")): 359.50,
            ("solubility_mass_frac_phase_comp", ("Liq", "NaCl")): 0.265,
            ("dens_mass_solvent", "Liq"): 996.89,
            ("dens_mass_phase", "Liq"): 1193.4,
            ("flow_vol_phase", "Liq"): 8.379e-4,
            ("conc_mass_phase_comp", ("Liq", "H2O")): 883.14,
            ("conc_mass_phase_comp", ("Liq", "NaCl")): 310.29,
            ("flow_mol_phase_comp", ("Liq", "H2O")): 41.08,
            ("flow_mol_phase_comp", ("Liq", "NaCl")): 4.449,
            ("mole_frac_phase_comp", ("Liq", "H2O")): 0.9022,
            ("mole_frac_phase_comp", ("Liq", "NaCl")): 9.773e-2,
            ("cp_mass_solvent", "Liq"): 4186.52,
            ("cp_mass_phase", "Liq"): 3276.56,
            ("enth_mass_solvent", "Liq"): 104.92,
            ("enth_mass_phase", "Liq"): 68.78,
        }


@pytest.mark.component
class TestNaClPropertySolution_3(PropertyRegressionTest):
    # Test pure liquid solution 3 - same solution as NaCl prop pack
    def configure(self):
        self.prop_pack = props.NaClParameterBlock
        self.param_args = {}

        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1e0,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 1e3,
            ("flow_mass_phase_comp", ("Vap", "H2O")): 1e3,
        }

        self.state_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 0.999,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 0.001,
            ("flow_mass_phase_comp", ("Sol", "NaCl")): 0.0,
            ("flow_mass_phase_comp", ("Vap", "H2O")): 0.0,
            ("temperature", None): 273.15 + 25,
            ("pressure", None): 10e5,
        }

        self.regression_solution = {
            ("mass_frac_phase_comp", ("Liq", "H2O")): 0.999,
            ("mass_frac_phase_comp", ("Liq", "NaCl")): 0.001,
            ("solubility_mass_phase_comp", ("Liq", "NaCl")): 359.50,
            ("solubility_mass_frac_phase_comp", ("Liq", "NaCl")): 0.265,
            ("dens_mass_solvent", "Liq"): 996.89,
            ("dens_mass_solvent", "Vap"): 7.267,
            ("dens_mass_phase", "Liq"): 997.59,
            ("flow_vol_phase", "Liq"): 1.002e-3,
            ("conc_mass_phase_comp", ("Liq", "H2O")): 996.59,
            ("conc_mass_phase_comp", ("Liq", "NaCl")): 0.9976,
            ("flow_mol_phase_comp", ("Liq", "H2O")): 55.45,
            ("flow_mol_phase_comp", ("Liq", "NaCl")): 1.711e-2,
            ("mole_frac_phase_comp", ("Liq", "H2O")): 0.9997,
            ("mole_frac_phase_comp", ("Liq", "NaCl")): 3.084e-4,
            ("cp_mass_solvent", "Liq"): 4186.52,
            ("cp_mass_phase", "Liq"): 4180.98,
            ("enth_mass_solvent", "Liq"): 104.92,
            ("enth_mass_phase", "Liq"): 104.40,
        }


@pytest.mark.component
class TestNaClPropertySolution_4(PropertyRegressionTest):
    # Test pure solid solution 1 - check solid properties
    def configure(self):

        self.prop_pack = props.NaClParameterBlock
        self.param_args = {}

        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1e1,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 1e3,
            ("flow_mass_phase_comp", ("Sol", "NaCl")): 1e0,
            ("flow_mass_phase_comp", ("Vap", "H2O")): 1e2,
        }

        self.state_args = {
            ("flow_vol_phase", "Liq"): 0,
            ("flow_vol_phase", "Vap"): 0,
            ("mass_frac_phase_comp", ("Liq", "NaCl")): 0.05,
            ("flow_mass_phase_comp", ("Sol", "NaCl")): 1.0,
            ("temperature", None): 273.15 + 25,
            ("pressure", None): 5e5,
        }

        self.regression_solution = {
            ("mass_frac_phase_comp", ("Sol", "NaCl")): 1.0,
            ("dens_mass_solute", "Sol"): 2115,
            ("cp_mass_solute", "Sol"): 864.16,
            ("flow_vol_phase", "Sol"): 1 / 2115,  # mass floe / density
            ("cp_mass_solute", "Sol"): 864.16,
            ("enth_mass_solute", "Sol"): 21.474,
            ("dh_crystallization_mass_comp", "NaCl"): -520,
            ("flow_mol_phase_comp", ("Sol", "NaCl")): 1 / 58.44e-3,  # mass flow / mw
            ("mole_frac_phase_comp", ("Sol", "NaCl")): 1.0,
        }


@pytest.mark.component
class TestNaClPropertySolution_5(PropertyRegressionTest):
    # Test pure vapor solution 1 - check vapor properties
    def configure(self):

        self.prop_pack = props.NaClParameterBlock
        self.param_args = {}

        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1e1,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 1e3,
            ("flow_mass_phase_comp", ("Sol", "NaCl")): 1e3,
            ("flow_mass_phase_comp", ("Vap", "H2O")): 1e0,
        }

        self.state_args = {
            ("flow_vol_phase", "Liq"): 0,
            ("flow_vol_phase", "Sol"): 0,
            ("mass_frac_phase_comp", ("Liq", "NaCl")): 0.05,
            ("flow_mass_phase_comp", ("Vap", "H2O")): 1.0,
            ("temperature", None): 273.15 + 25,
            ("pressure", None): 5e5,
        }

        self.regression_solution = {
            ("mass_frac_phase_comp", ("Vap", "H2O")): 1.0,
            ("dens_mass_solvent", "Vap"): 3.633,
            ("dh_vap_mass_solvent", None): 2441.808,
            ("cp_mass_solvent", "Vap"): 1864.52,
            ("flow_vol_phase", "Vap"): 1 / 3.633,  # mass flow / density
            ("pressure_sat", None): 2905.28,
            ("flow_mol_phase_comp", ("Vap", "H2O")): 1 / 18.01528e-3,  # mass flow / mw
            ("mole_frac_phase_comp", ("Vap", "H2O")): 1.0,
        }


@pytest.mark.component
class TestNaClPropertySolution_6(PropertyRegressionTest):
    # Test for S-L-V solution 1 with similar magnitude flowrates in all phases and high liquid salt. conc. - check all properties
    def configure(self):
        self.prop_pack = props.NaClParameterBlock
        self.param_args = {}

        self.scaling_args = {
            ("flow_mass_phase_comp", ("Sol", "NaCl")): 10,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 1e1,
        }

        self.state_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 0.75,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 0.25,
            ("flow_mass_phase_comp", ("Sol", "NaCl")): 0.25,
            ("flow_mass_phase_comp", ("Vap", "H2O")): 0.25,
            ("temperature", None): 273.15 + 50,
            ("pressure", None): 5e5,
        }

        self.regression_solution = {
            ("solubility_mass_phase_comp", ("Liq", "NaCl")): 362.93,
            ("solubility_mass_frac_phase_comp", ("Liq", "NaCl")): 0.2686,
            ("dens_mass_solvent", "Liq"): 988.04,
            ("dens_mass_solvent", "Vap"): 3.352,
            ("dens_mass_solute", "Liq"): 2645.21,
            ("dens_mass_solute", "Sol"): 2115,
            ("dens_mass_phase", "Liq"): 1171.53,
            ("dh_vap_mass_solvent", None): 2382.08,
            ("cp_mass_solvent", "Liq"): 4180.92,
            ("cp_mass_solvent", "Vap"): 1871.21,
            ("cp_mass_solute", "Sol"): 873.34,
            ("cp_mass_phase", "Liq"): 3300.83,
            ("flow_vol_phase", "Liq"): 8.53e-4,
            ("flow_vol_phase", "Sol"): 1.182e-4,
            ("flow_vol_phase", "Vap"): 7.46e-2,
            ("pressure_sat", None): 9799.91,
            ("temperature_sat_solvent", None): 318.52,
            ("conc_mass_phase_comp", ("Liq", "NaCl")): 292.88,
            ("conc_mass_phase_comp", ("Liq", "H2O")): 878.65,
            ("enth_mass_solvent", "Liq"): 209.40,
            ("enth_mass_solvent", "Vap"): 2591.48,
            ("enth_mass_solute", "Sol"): 43.19,
            ("enth_mass_phase", "Liq"): 152.36,
            ("dh_crystallization_mass_comp", "NaCl"): -520,
            ("mass_frac_phase_comp", ("Liq", "H2O")): 0.75,
            ("mass_frac_phase_comp", ("Liq", "NaCl")): 0.25,
            ("mass_frac_phase_comp", ("Sol", "NaCl")): 1,
            ("mass_frac_phase_comp", ("Vap", "H2O")): 1,
            ("flow_mol_phase_comp", ("Liq", "H2O")): 41.63,
            ("flow_mol_phase_comp", ("Liq", "NaCl")): 4.28,
            ("flow_mol_phase_comp", ("Sol", "NaCl")): 4.28,
            ("flow_mol_phase_comp", ("Vap", "H2O")): 13.88,
            ("mole_frac_phase_comp", ("Liq", "H2O")): 0.9068,
            ("mole_frac_phase_comp", ("Liq", "NaCl")): 0.0932,
            ("mole_frac_phase_comp", ("Sol", "NaCl")): 1.0,
            ("mole_frac_phase_comp", ("Vap", "H2O")): 1.0,
        }


@pytest.mark.component
class TestNaClPropertySolution_7(PropertyRegressionTest):
    # Test for S-L-V solution 2 with flowrates in all phases of same magnitude but low liquid salt. conc. - check all properties
    def configure(self):
        self.prop_pack = props.NaClParameterBlock
        self.param_args = {}

        self.scaling_args = {
            ("flow_mass_phase_comp", ("Sol", "NaCl")): 10,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 1e3,
        }

        self.state_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 0.999,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 0.001,
            ("flow_mass_phase_comp", ("Sol", "NaCl")): 0.25,
            ("flow_mass_phase_comp", ("Vap", "H2O")): 0.25,
            ("temperature", None): 273.15 + 50,
            ("pressure", None): 5e5,
        }

        self.regression_solution = {
            ("solubility_mass_phase_comp", ("Liq", "NaCl")): 362.93,
            ("solubility_mass_frac_phase_comp", ("Liq", "NaCl")): 0.2686,
            ("dens_mass_solvent", "Liq"): 988.04,
            ("dens_mass_solvent", "Vap"): 3.352,
            ("dens_mass_solute", "Liq"): 3073.05,
            ("dens_mass_solute", "Sol"): 2115,
            ("dens_mass_phase", "Liq"): 988.72,
            ("dh_vap_mass_solvent", None): 2382.08,
            ("cp_mass_solvent", "Liq"): 4180.92,
            ("cp_mass_solvent", "Vap"): 1871.21,
            ("cp_mass_solute", "Sol"): 873.34,
            ("cp_mass_phase", "Liq"): 4175.95,
            ("flow_vol_phase", "Liq"): 1.011e-3,
            ("flow_vol_phase", "Sol"): 1.182e-4,
            ("flow_vol_phase", "Vap"): 7.46e-2,
            ("pressure_sat", None): 12614.93,
            ("temperature_sat_solvent", None): 323.55,
            ("conc_mass_phase_comp", ("Liq", "NaCl")): 0.9887,
            ("conc_mass_phase_comp", ("Liq", "H2O")): 987.73,
            ("enth_mass_solvent", "Liq"): 209.40,
            ("enth_mass_solvent", "Vap"): 2591.48,
            ("enth_mass_solute", "Sol"): 43.19,
            ("enth_mass_phase", "Liq"): 208.81,
            ("dh_crystallization_mass_comp", "NaCl"): -520,
            ("mass_frac_phase_comp", ("Liq", "H2O")): 0.999,
            ("mass_frac_phase_comp", ("Liq", "NaCl")): 0.001,
            ("mass_frac_phase_comp", ("Sol", "NaCl")): 1,
            ("mass_frac_phase_comp", ("Vap", "H2O")): 1,
            ("flow_mol_phase_comp", ("Liq", "H2O")): 55.45,
            ("flow_mol_phase_comp", ("Liq", "NaCl")): 0.0171,
            ("flow_mol_phase_comp", ("Sol", "NaCl")): 4.28,
            ("flow_mol_phase_comp", ("Vap", "H2O")): 13.88,
            ("mole_frac_phase_comp", ("Liq", "H2O")): 0.999,
            ("mole_frac_phase_comp", ("Liq", "NaCl")): 3.084e-4,
            ("mole_frac_phase_comp", ("Sol", "NaCl")): 1.0,
            ("mole_frac_phase_comp", ("Vap", "H2O")): 1.0,
        }


@pytest.mark.component
class TestNaClPropertySolution_8(PropertyRegressionTest):
    # Test for S-L-V solution 3 with outlet data from Dutta et al. - proper crystallization system
    # # Dutta recorded properties for solids and liquids at crystallizer temperature:
    # # Solubility @ 55C: 0.27 kg/kg
    # # Heat of vaporization @ 55C: 2400 kJ/kg
    # # Vapor density @ 55C: 68.7e-3 kg/m3
    # # Solid heat capacity: 877 J/kgK
    # # Liquid heat capacity @ 20 C: 3290 kJ/kgK
    # # Liquid density @ 20 C : 1185 kg/m3

    def configure(self):
        self.prop_pack = props.NaClParameterBlock
        self.param_args = {}

        self.scaling_args = {
            ("flow_mass_phase_comp", ("Sol", "NaCl")): 1e0,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 1e0,
            ("flow_mass_phase_comp", ("Vap", "H2O")): 1e-1,
        }

        self.state_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 18.37,  # 84t/h
            ("mass_frac_phase_comp", ("Liq", "NaCl")): 0.2126,
            ("flow_mass_phase_comp", ("Sol", "NaCl")): 5.55,  # 20t/h
            ("flow_mass_phase_comp", ("Vap", "H2O")): 20.55,  # 74t/h
            ("temperature", None): 273.15 + 20,
            ("pressure", None): 10e3,
        }

        self.regression_solution = {
            ("solubility_mass_phase_comp", ("Liq", "NaCl")): 358.88,
            ("solubility_mass_frac_phase_comp", ("Liq", "NaCl")): 0.2644,
            ("dens_mass_solvent", "Liq"): 998.02,
            ("dens_mass_solvent", "Vap"): 73.91e-3,
            ("dens_mass_solute", "Liq"): 2847.63,
            ("dens_mass_solute", "Sol"): 2115,
            ("dens_mass_phase", "Liq"): 1157.91,
            ("dh_vap_mass_solvent", None): 2453.66,
            ("cp_mass_solvent", "Liq"): 4189.43,
            ("cp_mass_solvent", "Vap"): 1863.46,
            ("cp_mass_solute", "Sol"): 862.16,
            ("cp_mass_phase", "Liq"): 3380.07,
            ("flow_vol_phase", "Liq"): 0.02015,
            ("flow_vol_phase", "Sol"): 2.624e-3,
            ("flow_vol_phase", "Vap"): 278.04,
            ("pressure_sat", None): 1679.64,
            ("temperature_sat_solvent", None): 287.88,
            ("conc_mass_phase_comp", ("Liq", "NaCl")): 246.17,
            ("conc_mass_phase_comp", ("Liq", "H2O")): 911.74,
            ("enth_mass_solvent", "Liq"): 84.00,
            ("enth_mass_solvent", "Vap"): 2537.66,
            ("enth_mass_solute", "Sol"): 17.16,
            ("enth_mass_phase", "Liq"): 59.03,
            ("dh_crystallization_mass_comp", "NaCl"): -520,
            ("mass_frac_phase_comp", ("Liq", "H2O")): 0.7874,
            ("mass_frac_phase_comp", ("Sol", "NaCl")): 1,
            ("mass_frac_phase_comp", ("Vap", "H2O")): 1,
            ("flow_mol_phase_comp", ("Liq", "H2O")): 1019.69,
            ("flow_mol_phase_comp", ("Liq", "NaCl")): 84.88,
            ("flow_mol_phase_comp", ("Sol", "NaCl")): 94.97,
            ("flow_mol_phase_comp", ("Vap", "H2O")): 1140.698,
            ("mole_frac_phase_comp", ("Liq", "H2O")): 0.9232,
            ("mole_frac_phase_comp", ("Liq", "NaCl")): 0.0768,
            ("mole_frac_phase_comp", ("Sol", "NaCl")): 1.0,
            ("mole_frac_phase_comp", ("Vap", "H2O")): 1.0,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 4.96,
        }


@pytest.mark.component
class TestNaClPropertySolution_9(PropertyRegressionTest):
    # Test for S-L-V solution 4 with outlet data from Dutta et al. - proper crystallization system
    # # Dutta recorded properties for solids and liquids at crystallizer temperature:
    # # Solubility @ 55C: 0.27 kg/kg
    # # Heat of vaporization @ 55C: 2400 kJ/kg
    # # Vapor density @ 55C: 68.7e-3 kg/m3
    # # Solid heat capacity: 877 J/kgK
    # # Liquid heat capacity @ 20 C: 3290 kJ/kgK
    # # Liquid density @ 20 C : 1185 kg/m3

    def configure(self):
        self.prop_pack = props.NaClParameterBlock
        self.param_args = {}

        self.scaling_args = {
            ("flow_mass_phase_comp", ("Sol", "NaCl")): 1e0,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 1e0,
            ("flow_mass_phase_comp", ("Vap", "H2O")): 1e-1,
        }

        self.state_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 18.37,  # 84t/h
            ("mass_frac_phase_comp", ("Liq", "NaCl")): 0.2126,
            ("flow_mass_phase_comp", ("Sol", "NaCl")): 5.55,  # 20t/h
            ("flow_mass_phase_comp", ("Vap", "H2O")): 20.55,  # 74t/h
            ("temperature", None): 273.15 + 55,
            ("pressure", None): 10e3,
        }

        self.regression_solution = {
            ("solubility_mass_phase_comp", ("Liq", "NaCl")): 363.71,
            ("solubility_mass_frac_phase_comp", ("Liq", "NaCl")): 0.2695,
            ("dens_mass_solvent", "Liq"): 985.71,
            ("dens_mass_solvent", "Vap"): 66.03e-3,
            ("dens_mass_solute", "Liq"): 2694.61,
            ("dens_mass_solute", "Sol"): 2115,
            ("dens_mass_phase", "Liq"): 1139.33,
            ("dh_vap_mass_solvent", None): 2369.98,
            ("cp_mass_solvent", "Liq"): 4181.59,
            ("cp_mass_solvent", "Vap"): 1872.79,
            ("cp_mass_solute", "Sol"): 875.04,
            ("cp_mass_phase", "Liq"): 3390.43,
            ("flow_vol_phase", "Liq"): 0.02047,
            ("flow_vol_phase", "Sol"): 2.624e-3,
            ("flow_vol_phase", "Vap"): 311.233,
            ("pressure_sat", None): 13201.77,
            ("temperature_sat_solvent", None): 324.47,
            ("conc_mass_phase_comp", ("Liq", "NaCl")): 242.22,
            ("conc_mass_phase_comp", ("Liq", "H2O")): 897.107,
            ("enth_mass_solvent", "Liq"): 230.30,
            ("enth_mass_solvent", "Vap"): 2600.279,
            ("enth_mass_solute", "Sol"): 47.57,
            ("enth_mass_phase", "Liq"): 177.39,
            ("dh_crystallization_mass_comp", "NaCl"): -520,
            ("mass_frac_phase_comp", ("Liq", "H2O")): 0.7874,
            ("mass_frac_phase_comp", ("Sol", "NaCl")): 1,
            ("mass_frac_phase_comp", ("Vap", "H2O")): 1,
            ("flow_mol_phase_comp", ("Liq", "H2O")): 1019.69,
            ("flow_mol_phase_comp", ("Liq", "NaCl")): 84.88,
            ("flow_mol_phase_comp", ("Sol", "NaCl")): 94.97,
            ("flow_mol_phase_comp", ("Vap", "H2O")): 1140.698,
            ("mole_frac_phase_comp", ("Liq", "H2O")): 0.9232,
            ("mole_frac_phase_comp", ("Liq", "NaCl")): 0.0768,
            ("mole_frac_phase_comp", ("Sol", "NaCl")): 1.0,
            ("mole_frac_phase_comp", ("Vap", "H2O")): 1.0,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 4.96,
        }


@pytest.mark.component
class TestNaClCalculateState_1(PropertyCalculateStateTest):
    # Test pure liquid solution with mass fractions
    def configure(self):
        self.prop_pack = props.NaClParameterBlock
        self.param_args = {}

        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1e-1,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 1e0,
            ("flow_mass_phase_comp", ("Vap", "H2O")): 1,
            ("flow_mass_phase_comp", ("Sol", "NaCl")): 1,
        }

        self.var_args = {
            ("flow_vol_phase", "Liq"): 2e-2,
            ("flow_vol_phase", "Sol"): 0,
            ("flow_vol_phase", "Vap"): 0,
            ("mass_frac_phase_comp", ("Liq", "NaCl")): 0.05,
            ("temperature", None): 273.15 + 25,
            ("pressure", None): 5e5,
        }

        self.state_solution = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 19.6,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 1.032,
            ("flow_mass_phase_comp", ("Sol", "NaCl")): 0,
            ("flow_mass_phase_comp", ("Vap", "H2O")): 0,
        }


@pytest.mark.component
class TestNaClCalculateState_2(PropertyCalculateStateTest):
    # Test pure liquid with mole fractions
    def configure(self):
        self.prop_pack = props.NaClParameterBlock
        self.param_args = {}

        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1e1,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 1e2,
            # The rest are expected to be zero
            ("flow_mass_phase_comp", ("Sol", "NaCl")): 1,
            ("flow_mass_phase_comp", ("Vap", "H2O")): 1,
        }
        self.var_args = {
            ("flow_vol_phase", "Liq"): 2e-4,
            ("flow_vol_phase", "Sol"): 0,
            ("flow_vol_phase", "Vap"): 0,
            ("mole_frac_phase_comp", ("Liq", "NaCl")): 0.05,
            ("temperature", None): 273.15 + 25,
            ("pressure", None): 5e5,
        }
        self.state_solution = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 18.84e-2,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 3.215e-2,
            ("flow_mass_phase_comp", ("Sol", "NaCl")): 0,
            ("flow_mass_phase_comp", ("Vap", "H2O")): 0,
        }


@pytest.mark.component
class TestNaClCalculateState_3(PropertyCalculateStateTest):
    # Test pure liquid solution with pressure_sat defined instead of temperature
    def configure(self):
        self.prop_pack = props.NaClParameterBlock
        self.param_args = {}

        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1e0,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 1e1,
            ("flow_mass_phase_comp", ("Sol", "NaCl")): 1e-1,
            ("flow_mass_phase_comp", ("Vap", "H2O")): 1e2,
        }

        self.var_args = {
            ("flow_vol_phase", "Liq"): 2e-2,
            ("mass_frac_phase_comp", ("Liq", "NaCl")): 0.05,
            ("flow_vol_phase", "Sol"): 0,
            ("flow_vol_phase", "Vap"): 0,
            ("pressure_sat", None): 2905,
            ("pressure", None): 5e5,
        }

        self.state_solution = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 19.6,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 1.032,
            ("flow_mass_phase_comp", ("Sol", "NaCl")): 0,
            ("flow_mass_phase_comp", ("Vap", "H2O")): 0,
            ("temperature", None): 273.15 + 25,
        }


@pytest.mark.component
class TestNaClCalculateState_4(PropertyCalculateStateTest):
    # Test pure solid solution with mass fractions
    def configure(self):

        self.prop_pack = props.NaClParameterBlock
        self.param_args = {}

        self.scaling_args = {
            ("flow_mass_phase_comp", ("Sol", "NaCl")): 1e-1,
            ("flow_mass_phase_comp", ("Vap", "H2O")): 1e-1,
        }

        self.var_args = {
            ("flow_vol_phase", "Liq"): 0,
            ("mass_frac_phase_comp", ("Liq", "NaCl")): 0,
            ("flow_vol_phase", "Sol"): 2e-2,
            ("flow_vol_phase", "Vap"): 0,
            ("temperature", None): 273.15 + 25,
            ("pressure", None): 5e5,
        }

        self.state_solution = {
            ("flow_mass_phase_comp", ("Sol", "NaCl")): 2115
            * 2e-2,  # solid density is constant
            ("flow_mass_phase_comp", ("Liq", "H2O")): 0,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 0,
            ("flow_mass_phase_comp", ("Vap", "H2O")): 0,
        }


@pytest.mark.component
class TestNaClCalculateState_5(PropertyCalculateStateTest):
    # Test solid-liquid-vapor mixture solution with mass fractions
    def configure(self):
        self.prop_pack = props.NaClParameterBlock
        self.param_args = {}

        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1e-1,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 1e1,
            ("flow_mass_phase_comp", ("Sol", "NaCl")): 1e-2,
        }

        self.var_args = {
            ("flow_vol_phase", "Liq"): 2e-2,
            ("mass_frac_phase_comp", ("Liq", "NaCl")): 0.05,
            ("flow_vol_phase", "Sol"): 2e-2,
            ("flow_vol_phase", "Vap"): 2e-2,
            ("temperature", None): 273.15 + 25,
            ("pressure", None): 5e5,
        }

        self.state_solution = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 19.6,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 1.032,
            ("flow_mass_phase_comp", ("Vap", "H2O")): 0.07265,
            ("flow_mass_phase_comp", ("Sol", "NaCl")): 42.3,
        }


@pytest.mark.component
class TestNaClCalculateState_6(PropertyCalculateStateTest):
    # Test liquid-solid-vapor mixture with mole fractions
    def configure(self):
        self.prop_pack = props.NaClParameterBlock
        self.param_args = {}

        self.scaling_args = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1e1,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 1e2,
            ("flow_mass_phase_comp", ("Sol", "NaCl")): 1e2,
            ("flow_mass_phase_comp", ("Vap", "H2O")): 1e3,
        }

        self.var_args = {
            ("flow_vol_phase", "Liq"): 2e-4,
            ("flow_vol_phase", "Sol"): 2e-4,
            ("flow_vol_phase", "Vap"): 2e-4,
            ("mole_frac_phase_comp", ("Liq", "NaCl")): 0.05,
            ("temperature", None): 273.15 + 25,
            ("pressure", None): 5e5,
        }

        self.state_solution = {
            ("flow_mass_phase_comp", ("Liq", "H2O")): 18.84e-2,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 3.215e-2,
            ("flow_mass_phase_comp", ("Sol", "NaCl")): 0.423,
            ("flow_mass_phase_comp", ("Vap", "H2O")): 3.632
            * 2e-4,  # Density from ideal gas law * vol. flow
        }
