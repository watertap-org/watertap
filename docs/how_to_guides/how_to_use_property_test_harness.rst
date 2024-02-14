.. _how_to_use_property_test_harness:

How to use the property test harness
====================================

Overview
--------

This guide shows you how to use the property test harness to generate tests for WaterTAP property models. The purpose of this
tool is to standardize testing so developers don't need to write tests from the ground-up for each property model.

How To
------

Begin by importing the following functions - note that the property model import will differ for each test file.
The following example assumes a test file is being created for the NaCl property package.

.. testsetup::

    # quiet idaes logs
    import idaes.logger as idaeslogger
    idaeslogger.getLogger('ideas.core').setLevel('CRITICAL')
    idaeslogger.getLogger('ideas.core.util.scaling').setLevel('CRITICAL')
    idaeslogger.getLogger('idaes.init').setLevel('CRITICAL')

.. testcode::

    import pytest

    # Import the property test harness functions
    from watertap.property_models.tests.property_test_harness import PropertyTestHarness
    from watertap.property_models.tests.property_test_harness import PropertyRegressionTest
    from idaes.models.properties.tests.test_harness import PropertyTestHarness as PropertyTestHarness_idaes

    # Import the property model to be tested
    import watertap.property_models.NaCl_prop_pack as props

Next, test the configuration of the property package against the IDAES property test harness.

.. testcode::

    @pytest.mark.unit
    class TestNaClProperty_idaes(PropertyTestHarness_idaes):
        def configure(self):
            self.prop_pack = props.NaClParameterBlock
            self.param_args = {}
            self.prop_args = {}
            self.has_density_terms = False

Then, test the outputs of the property package against the WaterTAP property test harness by specifying the property package, scaling,
stateblock statistics, and the expected solutions for the model's variables.

.. testcode::

    class TestNaClProperty(PropertyTestHarness):
        def configure(self):
            self.prop_pack = props.NaClParameterBlock
            self.param_args = {}
            self.scaling_args = {
                ("flow_mass_phase_comp", ("Liq", "H2O")): 1,
                ("flow_mass_phase_comp", ("Liq", "NaCl")): 1e2,
            }
            self.stateblock_statistics = {
                "number_variables": 20,
                "number_total_constraints": 16,
                "number_unused_variables": 0,
                "default_degrees_of_freedom": 4,
            }  # 4 state vars, but pressure is not active
            self.default_solution = {
                ("mass_frac_phase_comp", ("Liq", "H2O")): 0.965,
                ("mass_frac_phase_comp", ("Liq", "NaCl")): 0.035,
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

Finally, test the regression outputs of the property model by specifying the property package, solver, state arguments, and the expected solutions.

.. testcode::

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
                ("temperature", None): 273.15 + 25,
                ("pressure", None): 50e5,
            }
            self.regression_solution = {
                ("mass_frac_phase_comp", ("Liq", "H2O")): 0.95,
                ("mass_frac_phase_comp", ("Liq", "NaCl")): 0.05,
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

