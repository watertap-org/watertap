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
import numpy as np
import idaes.core.util.scaling as iscale
import pytest
from idaes.core import (
    FlowsheetBlock,
)

from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from pyomo.environ import (
    ConcreteModel,
    assert_optimal_termination,
    value,
)

from watertap.unit_models.electrodialysis_bipolar_1D import (
    Electrodialysis_Bipolar_1D,
    ElectricalOperationMode,
    LimitingCurrentDensitybpemMethod,
    PressureDropMethod,
    FrictionFactorMethod,
    HydraulicDiameterMethod,
)
from watertap.core.solvers import get_solver
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock


__author__ = "Johnson Dhanasekaran"

solver = get_solver()


# -----------------------------------------------------------------------------
# Start test class


class Test_membrane_characteristics:
    @pytest.fixture(scope="class")
    def bped(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        ion_dict = {
            "solute_list": ["Na_+", "Cl_-", "H_+", "OH_-"],
            "mw_data": {
                "Na_+": 23e-3,
                "Cl_-": 35.5e-3,
                "H_+": 1e-3,
                "OH_-": 17.0e-3,
            },
            "elec_mobility_data": {
                ("Liq", "Na_+"): 5.19e-8,
                ("Liq", "Cl_-"): 7.92e-8,
                ("Liq", "H_+"): 36.23e-8,
                ("Liq", "OH_-"): 20.64e-8,
            },
            "charge": {"Na_+": 1, "Cl_-": -1, "H_+": 1, "OH_-": -1},
            "diffusivity_data": {
                ("Liq", "Na_+"): 1.33e-9,
                ("Liq", "Cl_-"): 2.03e-9,
                ("Liq", "H_+"): 9.31e-9,
                ("Liq", "OH_-"): 5.27e-9,
            },
        }
        m.fs.properties = MCASParameterBlock(**ion_dict)
        m.fs.unit = Electrodialysis_Bipolar_1D(
            property_package=m.fs.properties,
            # operation_mode=ElectricalOperationMode.Constant_Voltage,
            limiting_current_density_method_bpem=LimitingCurrentDensitybpemMethod.Empirical,
        )

        m.fs.unit.diffus_mass.fix((2.03 + 1.96) * 10**-9 / 2)
        m.fs.unit.conc_water.fix(50 * 1e3)
        m.fs.unit.k2_zero.fix(2 * 10**-6)
        m.fs.unit.relative_permittivity.fix(30)
        m.fs.unit.membrane_fixed_catalyst_cel.fix(5e3)
        m.fs.unit.membrane_fixed_catalyst_ael.fix(5e3)
        m.fs.unit.k_a.fix(447)
        m.fs.unit.k_b.fix(5e4)

        m.fs.unit.ion_trans_number_membrane["bpem", "Na_+"].fix(0)
        m.fs.unit.ion_trans_number_membrane["bpem", "Cl_-"].fix(0)
        m.fs.unit.ion_trans_number_membrane["bpem", "H_+"].fix(1)
        m.fs.unit.ion_trans_number_membrane["bpem", "OH_-"].fix(1)

        m.fs.unit.ion_trans_number_membrane["cem", "Na_+"].fix(0.94)
        m.fs.unit.ion_trans_number_membrane["cem", "Cl_-"].fix(0)
        m.fs.unit.ion_trans_number_membrane["cem", "H_+"].fix(0.03)
        m.fs.unit.ion_trans_number_membrane["cem", "OH_-"].fix(0.03)

        m.fs.unit.ion_trans_number_membrane["aem", "Na_+"].fix(0)
        m.fs.unit.ion_trans_number_membrane["aem", "Cl_-"].fix(0.94)
        m.fs.unit.ion_trans_number_membrane["aem", "H_+"].fix(0.03)
        m.fs.unit.ion_trans_number_membrane["aem", "OH_-"].fix(0.03)

        m.fs.unit.solute_diffusivity_membrane["bpem", "Na_+"].fix(0)
        m.fs.unit.solute_diffusivity_membrane["bpem", "Cl_-"].fix(0)
        m.fs.unit.solute_diffusivity_membrane["bpem", "H_+"].fix(0)
        m.fs.unit.solute_diffusivity_membrane["bpem", "OH_-"].fix(0)

        m.fs.unit.solute_diffusivity_membrane["cem", "Na_+"].fix(
            (1.8e-10 + 1.25e-10) / 2
        )
        m.fs.unit.solute_diffusivity_membrane["cem", "Cl_-"].fix(
            (1.8e-10 + 1.25e-10) / 2
        )
        m.fs.unit.solute_diffusivity_membrane["cem", "H_+"].fix(0)
        m.fs.unit.solute_diffusivity_membrane["cem", "OH_-"].fix(0)

        m.fs.unit.solute_diffusivity_membrane["aem", "Na_+"].fix(
            (1.8e-10 + 1.25e-10) / 2
        )
        m.fs.unit.solute_diffusivity_membrane["aem", "Cl_-"].fix(
            (1.8e-10 + 1.25e-10) / 2
        )
        m.fs.unit.solute_diffusivity_membrane["aem", "H_+"].fix(0)
        m.fs.unit.solute_diffusivity_membrane["aem", "OH_-"].fix(0)

        # Set inlet streams.
        m.fs.unit.inlet_basic.pressure.fix(101325)
        m.fs.unit.inlet_basic.temperature.fix(298.15)
        m.fs.unit.inlet_acidic.pressure.fix(101325)
        m.fs.unit.inlet_acidic.temperature.fix(298.15)
        m.fs.unit.inlet_diluate.pressure.fix(101325)
        m.fs.unit.inlet_diluate.temperature.fix(298.15)
        m.fs.unit.spacer_porosity.fix(1)

        m.fs.unit.inlet_basic.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(7.38e-3)
        m.fs.unit.inlet_basic.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(7.38e-3)
        m.fs.unit.inlet_basic.flow_mol_phase_comp[0, "Liq", "H_+"].fix(7.38e-3)
        m.fs.unit.inlet_basic.flow_mol_phase_comp[0, "Liq", "OH_-"].fix(7.38e-3)
        m.fs.unit.inlet_acidic.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(7.38e-3)
        m.fs.unit.inlet_acidic.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(7.38e-3)
        m.fs.unit.inlet_acidic.flow_mol_phase_comp[0, "Liq", "H_+"].fix(7.38e-3)
        m.fs.unit.inlet_acidic.flow_mol_phase_comp[0, "Liq", "OH_-"].fix(7.38e-3)
        m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(7.38e-3)
        m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(7.38e-3)
        m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "H_+"].fix(7.38e-3)
        m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "OH_-"].fix(7.38e-3)

        m.fs.unit.inlet_basic.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.40e-1)
        m.fs.unit.inlet_acidic.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.40e-1)
        m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.40e-1)

        m.fs.unit.shadow_factor.fix(1)

        m.fs.unit.water_trans_number_membrane["cem"].fix((5.8 + 4.3) / 2)
        m.fs.unit.water_permeability_membrane["cem"].fix((2.16e-14 + 1.75e-14) / 2)

        m.fs.unit.water_trans_number_membrane["aem"].fix((5.8 + 4.3) / 2)
        m.fs.unit.water_permeability_membrane["aem"].fix((2.16e-14 + 1.75e-14) / 2)

        m.fs.unit.water_trans_number_membrane["bpem"].fix((5.8 + 4.3) / 2)
        m.fs.unit.water_permeability_membrane["bpem"].fix((2.16e-14 + 1.75e-14) / 2)
        m.fs.unit.electrodes_resistance.fix(0)
        m.fs.unit.current_utilization.fix(1)
        m.fs.unit.channel_height.fix(2.7e-4)
        m.fs.unit.membrane_areal_resistance_coef_0.fix((1.89e-4 + 1.77e-4) / 2)
        m.fs.unit.membrane_areal_resistance_coef_1.fix(0)
        m.fs.unit.cell_width.fix(0.1)
        m.fs.unit.cell_length.fix(0.79)
        m.fs.unit.membrane_thickness["aem"].fix(4e-4)
        m.fs.unit.membrane_thickness["cem"].fix(4e-4)

        # Set scaling of critical quantities.
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e3, index=("Liq", "Na_+")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e3, index=("Liq", "Cl_-")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e3, index=("Liq", "H_+")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e3, index=("Liq", "OH_-")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e2, index=("Liq", "H2O")
        )
        m.fs.unit.cell_triplet_num.fix(1)
        iscale.set_scaling_factor(m.fs.unit.k_a, 1e-2)
        iscale.set_scaling_factor(m.fs.unit.k_b, 1e-4)
        iscale.set_scaling_factor(m.fs.unit.voltage_x, 1e-1)
        iscale.set_scaling_factor(m.fs.unit.flux_splitting, 1e4)
        iscale.set_scaling_factor(m.fs.unit.current_density_x, 1e-3)
        return m

    @pytest.mark.unit
    def test_IV_curve(self, bped):

        m = bped

        m.fs.unit.salt_conc_ael_ref.fix(2000)
        m.fs.unit.salt_conc_cel_ref.fix(2000)
        m.fs.unit.membrane_fixed_charge.fix(5e3)
        m.fs.unit.membrane_thickness["bpem"].fix(8e-4)

        # Experimental data from Wilhelm et al. (2002) with additional inputs from Mareev et al. (2020)
        expt_membrane_potential = np.array(
            [0.088, 0.184, 0.4045, 0.690, 0.858, 0.95, 1.3]
        )  # in volts

        # experimental current density: 11.7, 15.5, 20.8, 24.4, 30.6, 40.0, 100.6]  # in mA/cm2
        expected_current_density = np.array(
            [19.250, 19.269, 19.606, 22.764, 29.009, 35.278, 99.223]
        )  # in mA/cm2 (computed numerically)

        for indx, v in enumerate(expt_membrane_potential):

            m.fs.unit.voltage_membrane_drop[0, 1].fix(v)

            iscale.calculate_scaling_factors(m.fs)
            assert degrees_of_freedom(m) == 0
            initialization_tester(m)
            results = solver.solve(m)
            assert_optimal_termination(results)

            current_density_computed = (
                0.1 * m.fs.unit.current_density_x[0, 1]
            )  # convert to mA/cm2
            assert value(current_density_computed) == pytest.approx(
                expected_current_density[indx], rel=1e-3
            )
            m.fs.unit.voltage_membrane_drop.unfix()

    @pytest.mark.unit
    def test_limiting_current(self, bped):

        # Data on limiting current and potential barrier to water splitting have been obtained from:
        # Fumatech, Technical Data Sheet for Fumasep FBM, 2020. With additional modelling parameters obtained from
        # Ionescu, Viorel. Advanced Topics in Optoelectronics, Microelectronics, and Nanotechnologies (2023)
        # Expected current density is 1000 A/m2

        m = bped

        m.fs.unit.salt_conc_ael_ref.fix(500 + 2 * 250)
        m.fs.unit.salt_conc_cel_ref.fix(500 + 2 * 250)
        m.fs.unit.membrane_fixed_charge.fix(1.5e3)
        m.fs.unit.membrane_thickness["bpem"].fix(1.3e-4)
        m.fs.unit.current_applied.fix(1e2)

        # Test computing limiting current in  bipolar membrane
        iscale.calculate_scaling_factors(m)
        assert degrees_of_freedom(m) == 0
        initialization_tester(m)
        results = solver.solve(m)
        assert_optimal_termination(results)
        assert value(m.fs.unit.current_dens_lim_bpem[0, 1]) == pytest.approx(
            987.120, rel=1e-3
        )


class Test_Operation:

    @pytest.fixture(scope="class")
    def bped_const_v_no_salt_calc(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        ion_dict = {
            "solute_list": ["Na_+", "Cl_-", "H_+", "OH_-"],
            "mw_data": {
                "Na_+": 23e-3,
                "Cl_-": 35.5e-3,
                "H_+": 1e-3,
                "OH_-": 17.0e-3,
            },
            "elec_mobility_data": {
                ("Liq", "Na_+"): 5.19e-8,
                ("Liq", "Cl_-"): 7.92e-8,
                ("Liq", "H_+"): 36.23e-8,
                ("Liq", "OH_-"): 20.64e-8,
            },
            "charge": {"Na_+": 1, "Cl_-": -1, "H_+": 1, "OH_-": -1},
            "diffusivity_data": {
                ("Liq", "Na_+"): 1.33e-9,
                ("Liq", "Cl_-"): 2.03e-9,
                ("Liq", "H_+"): 9.31e-9,
                ("Liq", "OH_-"): 5.27e-9,
            },
        }
        m.fs.properties = MCASParameterBlock(**ion_dict)
        m.fs.unit = Electrodialysis_Bipolar_1D(
            property_package=m.fs.properties,
            operation_mode=ElectricalOperationMode.Constant_Voltage,
            limiting_current_density_method_bpem=LimitingCurrentDensitybpemMethod.Empirical,
            salt_calculation=False,
        )
        return m

    @pytest.fixture(scope="class")
    def bped_const_current_no_salt_calc(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        ion_dict = {
            "solute_list": ["Na_+", "Cl_-", "H_+", "OH_-"],
            "mw_data": {
                "Na_+": 23e-3,
                "Cl_-": 35.5e-3,
                "H_+": 1e-3,
                "OH_-": 17.0e-3,
            },
            "elec_mobility_data": {
                ("Liq", "Na_+"): 5.19e-8,
                ("Liq", "Cl_-"): 7.92e-8,
                ("Liq", "H_+"): 36.23e-8,
                ("Liq", "OH_-"): 20.64e-8,
            },
            "charge": {"Na_+": 1, "Cl_-": -1, "H_+": 1, "OH_-": -1},
            "diffusivity_data": {
                ("Liq", "Na_+"): 1.33e-9,
                ("Liq", "Cl_-"): 2.03e-9,
                ("Liq", "H_+"): 9.31e-9,
                ("Liq", "OH_-"): 5.27e-9,
            },
        }
        m.fs.properties = MCASParameterBlock(**ion_dict)
        m.fs.unit = Electrodialysis_Bipolar_1D(
            property_package=m.fs.properties,
            operation_mode=ElectricalOperationMode.Constant_Current,
            limiting_current_density_method_bpem=LimitingCurrentDensitybpemMethod.Empirical,
            salt_calculation=False,
        )
        return m

    @pytest.fixture(scope="class")
    def bped_const_v_salt_calc(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        ion_dict = {
            "solute_list": ["Na_+", "Cl_-", "H_+", "OH_-"],
            "mw_data": {
                "Na_+": 23e-3,
                "Cl_-": 35.5e-3,
                "H_+": 1e-3,
                "OH_-": 17.0e-3,
            },
            "elec_mobility_data": {
                ("Liq", "Na_+"): 5.19e-8,
                ("Liq", "Cl_-"): 7.92e-8,
                ("Liq", "H_+"): 36.23e-8,
                ("Liq", "OH_-"): 20.64e-8,
            },
            "charge": {"Na_+": 1, "Cl_-": -1, "H_+": 1, "OH_-": -1},
            "diffusivity_data": {
                ("Liq", "Na_+"): 1.33e-9,
                ("Liq", "Cl_-"): 2.03e-9,
                ("Liq", "H_+"): 9.31e-9,
                ("Liq", "OH_-"): 5.27e-9,
            },
        }
        m.fs.properties = MCASParameterBlock(**ion_dict)
        m.fs.unit = Electrodialysis_Bipolar_1D(
            property_package=m.fs.properties,
            operation_mode=ElectricalOperationMode.Constant_Voltage,
            limiting_current_density_method_bpem=LimitingCurrentDensitybpemMethod.Empirical,
            salt_calculation=True,
        )
        return m

    @pytest.fixture(scope="class")
    def bped_const_current_salt_calc(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        ion_dict = {
            "solute_list": ["Na_+", "Cl_-", "H_+", "OH_-"],
            "mw_data": {
                "Na_+": 23e-3,
                "Cl_-": 35.5e-3,
                "H_+": 1e-3,
                "OH_-": 17.0e-3,
            },
            "elec_mobility_data": {
                ("Liq", "Na_+"): 5.19e-8,
                ("Liq", "Cl_-"): 7.92e-8,
                ("Liq", "H_+"): 36.23e-8,
                ("Liq", "OH_-"): 20.64e-8,
            },
            "charge": {"Na_+": 1, "Cl_-": -1, "H_+": 1, "OH_-": -1},
            "diffusivity_data": {
                ("Liq", "Na_+"): 1.33e-9,
                ("Liq", "Cl_-"): 2.03e-9,
                ("Liq", "H_+"): 9.31e-9,
                ("Liq", "OH_-"): 5.27e-9,
            },
        }
        m.fs.properties = MCASParameterBlock(**ion_dict)
        m.fs.unit = Electrodialysis_Bipolar_1D(
            property_package=m.fs.properties,
            operation_mode=ElectricalOperationMode.Constant_Current,
            limiting_current_density_method_bpem=LimitingCurrentDensitybpemMethod.Empirical,
            salt_calculation=True,
        )
        return m

    @pytest.mark.unit
    def test_assign_evaluate(
        self,
        bped_const_v_no_salt_calc,
        bped_const_current_no_salt_calc,
        bped_const_v_salt_calc,
        bped_const_current_salt_calc,
    ):
        bped_m = (
            bped_const_v_no_salt_calc,
            bped_const_current_no_salt_calc,
            bped_const_v_salt_calc,
            bped_const_current_salt_calc,
        )
        for m in bped_m:

            m.fs.unit.diffus_mass.fix((2.03 + 1.96) * 10**-9 / 2)
            m.fs.unit.membrane_fixed_charge.fix(5e3)
            m.fs.unit.salt_conc_ael_ref.fix(2000)
            m.fs.unit.salt_conc_cel_ref.fix(2000)
            m.fs.unit.conc_water.fix(50 * 1e3)
            m.fs.unit.k2_zero.fix(2 * 10**-6)
            m.fs.unit.relative_permittivity.fix(30)
            m.fs.unit.membrane_fixed_catalyst_cel.fix(5e3)
            m.fs.unit.membrane_fixed_catalyst_ael.fix(5e3)
            m.fs.unit.k_a.fix(447)
            m.fs.unit.k_b.fix(5e4)

            m.fs.unit.ion_trans_number_membrane["bpem", "Na_+"].fix(0)
            m.fs.unit.ion_trans_number_membrane["bpem", "Cl_-"].fix(0)
            m.fs.unit.ion_trans_number_membrane["bpem", "H_+"].fix(1)
            m.fs.unit.ion_trans_number_membrane["bpem", "OH_-"].fix(1)

            m.fs.unit.ion_trans_number_membrane["cem", "Na_+"].fix(0.94)
            m.fs.unit.ion_trans_number_membrane["cem", "Cl_-"].fix(0)
            m.fs.unit.ion_trans_number_membrane["cem", "H_+"].fix(0.03)
            m.fs.unit.ion_trans_number_membrane["cem", "OH_-"].fix(0.03)

            m.fs.unit.ion_trans_number_membrane["aem", "Na_+"].fix(0)
            m.fs.unit.ion_trans_number_membrane["aem", "Cl_-"].fix(0.94)
            m.fs.unit.ion_trans_number_membrane["aem", "H_+"].fix(0.03)
            m.fs.unit.ion_trans_number_membrane["aem", "OH_-"].fix(0.03)

            m.fs.unit.solute_diffusivity_membrane["bpem", "Na_+"].fix(0)
            m.fs.unit.solute_diffusivity_membrane["bpem", "Cl_-"].fix(0)
            m.fs.unit.solute_diffusivity_membrane["bpem", "H_+"].fix(0)
            m.fs.unit.solute_diffusivity_membrane["bpem", "OH_-"].fix(0)

            m.fs.unit.solute_diffusivity_membrane["cem", "Na_+"].fix(
                (1.8e-10 + 1.25e-10) / 2
            )
            m.fs.unit.solute_diffusivity_membrane["cem", "Cl_-"].fix(
                (1.8e-10 + 1.25e-10) / 2
            )
            m.fs.unit.solute_diffusivity_membrane["cem", "H_+"].fix(0)
            m.fs.unit.solute_diffusivity_membrane["cem", "OH_-"].fix(0)

            m.fs.unit.solute_diffusivity_membrane["aem", "Na_+"].fix(
                (1.8e-10 + 1.25e-10) / 2
            )
            m.fs.unit.solute_diffusivity_membrane["aem", "Cl_-"].fix(
                (1.8e-10 + 1.25e-10) / 2
            )
            m.fs.unit.solute_diffusivity_membrane["aem", "H_+"].fix(0)
            m.fs.unit.solute_diffusivity_membrane["aem", "OH_-"].fix(0)

            # Set inlet streams.
            m.fs.unit.inlet_basic.pressure.fix(101325)
            m.fs.unit.inlet_basic.temperature.fix(298.15)
            m.fs.unit.inlet_acidic.pressure.fix(101325)
            m.fs.unit.inlet_acidic.temperature.fix(298.15)
            m.fs.unit.inlet_diluate.pressure.fix(101325)
            m.fs.unit.inlet_diluate.temperature.fix(298.15)
            m.fs.unit.spacer_porosity.fix(1)

            m.fs.unit.inlet_basic.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(7.38e-2)
            m.fs.unit.inlet_basic.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(0)
            m.fs.unit.inlet_basic.flow_mol_phase_comp[0, "Liq", "H_+"].fix(0)
            m.fs.unit.inlet_basic.flow_mol_phase_comp[0, "Liq", "OH_-"].fix(7.38e-2)
            m.fs.unit.inlet_acidic.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(0)
            m.fs.unit.inlet_acidic.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(7.38e-2)
            m.fs.unit.inlet_acidic.flow_mol_phase_comp[0, "Liq", "H_+"].fix(7.38e-2)
            m.fs.unit.inlet_acidic.flow_mol_phase_comp[0, "Liq", "OH_-"].fix(0)
            m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(7.38e-2)
            m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(7.38e-2)
            m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "H_+"].fix(0)
            m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "OH_-"].fix(0)

            m.fs.unit.inlet_basic.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.40e-1)
            m.fs.unit.inlet_acidic.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.40e-1)
            m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.40e-1)

            m.fs.unit.shadow_factor.fix(1)

            m.fs.unit.water_trans_number_membrane["cem"].fix((5.8 + 4.3) / 2)
            m.fs.unit.water_permeability_membrane["cem"].fix((2.16e-14 + 1.75e-14) / 2)

            m.fs.unit.water_trans_number_membrane["aem"].fix((5.8 + 4.3) / 2)
            m.fs.unit.water_permeability_membrane["aem"].fix((2.16e-14 + 1.75e-14) / 2)

            m.fs.unit.water_trans_number_membrane["bpem"].fix((5.8 + 4.3) / 2)
            m.fs.unit.water_permeability_membrane["bpem"].fix((2.16e-14 + 1.75e-14) / 2)
            m.fs.unit.electrodes_resistance.fix(0)
            m.fs.unit.current_utilization.fix(1)
            m.fs.unit.channel_height.fix(2.7e-4)
            m.fs.unit.membrane_areal_resistance_coef_0.fix((1.89e-4 + 1.77e-4) / 2)
            m.fs.unit.membrane_areal_resistance_coef_1.fix(0)
            m.fs.unit.cell_width.fix(0.1)
            m.fs.unit.cell_length.fix(0.79)
            m.fs.unit.membrane_thickness["bpem"].fix(8e-4)
            m.fs.unit.membrane_thickness["aem"].fix(4e-4)
            m.fs.unit.membrane_thickness["cem"].fix(4e-4)

            # Set scaling of critical quantities.
            m.fs.properties.set_default_scaling(
                "flow_mol_phase_comp", 1e3, index=("Liq", "Na_+")
            )
            m.fs.properties.set_default_scaling(
                "flow_mol_phase_comp", 1e3, index=("Liq", "Cl_-")
            )
            m.fs.properties.set_default_scaling(
                "flow_mol_phase_comp", 1e3, index=("Liq", "H_+")
            )
            m.fs.properties.set_default_scaling(
                "flow_mol_phase_comp", 1e3, index=("Liq", "OH_-")
            )
            m.fs.properties.set_default_scaling(
                "flow_mol_phase_comp", 1e2, index=("Liq", "H2O")
            )
            m.fs.unit.cell_triplet_num.fix(10)
            iscale.set_scaling_factor(m.fs.unit.k_a, 1e-2)
            iscale.set_scaling_factor(m.fs.unit.k_b, 1e-4)
            iscale.set_scaling_factor(m.fs.unit.voltage_x, 1e-1)
            iscale.set_scaling_factor(m.fs.unit.flux_splitting, 1e4)
            iscale.set_scaling_factor(m.fs.unit.current_density_x, 1e-3)

        # Test constant Voltage operation with no salt calculation

        bped_m[0].fs.unit.voltage_applied.fix(1e1)
        iscale.calculate_scaling_factors(bped_m[0])
        assert degrees_of_freedom(bped_m[0]) == 0
        initialization_tester(bped_m[0])
        assert value(
            bped_m[0].fs.unit.outlet_diluate.flow_mol_phase_comp[0, "Liq", "Na_+"]
        ) == pytest.approx(0.06896, rel=1e-3)
        assert value(
            bped_m[0].fs.unit.outlet_diluate.flow_mol_phase_comp[0, "Liq", "Cl_-"]
        ) == pytest.approx(0.06896, rel=1e-3)
        assert value(
            bped_m[0].fs.unit.outlet_diluate.flow_mol_phase_comp[0, "Liq", "H2O"]
        ) == pytest.approx(0.2002, rel=1e-3)

        assert value(
            bped_m[0].fs.unit.outlet_acidic.flow_mol_phase_comp[0, "Liq", "H_+"]
        ) == pytest.approx(0.074934, rel=1e-3)
        assert value(
            bped_m[0].fs.unit.outlet_acidic.flow_mol_phase_comp[0, "Liq", "Na_+"]
        ) == pytest.approx(0.002597, rel=1e-3)
        assert value(
            bped_m[0].fs.unit.outlet_acidic.flow_mol_phase_comp[0, "Liq", "Cl_-"]
        ) == pytest.approx(0.07604, rel=1e-3)
        assert value(
            bped_m[0].fs.unit.outlet_acidic.flow_mol_phase_comp[0, "Liq", "H2O"]
        ) == pytest.approx(0.2734, rel=1e-3)

        assert value(
            bped_m[0].fs.unit.outlet_basic.flow_mol_phase_comp[0, "Liq", "OH_-"]
        ) == pytest.approx(0.074934, rel=1e-3)
        assert value(
            bped_m[0].fs.unit.outlet_basic.flow_mol_phase_comp[0, "Liq", "Na_+"]
        ) == pytest.approx(0.07604, rel=1e-3)
        assert value(
            bped_m[0].fs.unit.outlet_basic.flow_mol_phase_comp[0, "Liq", "Cl_-"]
        ) == pytest.approx(0.002597, rel=1e-3)
        assert value(
            bped_m[0].fs.unit.outlet_basic.flow_mol_phase_comp[0, "Liq", "H2O"]
        ) == pytest.approx(0.24515, rel=1e-3)

        # Test constant Current operation with no salt calculation

        bped_m[1].fs.unit.current_applied.fix(2e2)
        iscale.calculate_scaling_factors(bped_m[1])
        assert degrees_of_freedom(bped_m[1]) == 0
        initialization_tester(bped_m[1])
        assert value(
            bped_m[1].fs.unit.outlet_diluate.flow_mol_phase_comp[0, "Liq", "Na_+"]
        ) == pytest.approx(0.05111, rel=1e-3)
        assert value(
            bped_m[1].fs.unit.outlet_diluate.flow_mol_phase_comp[0, "Liq", "Cl_-"]
        ) == pytest.approx(0.05054, rel=1e-3)
        assert value(
            bped_m[1].fs.unit.outlet_diluate.flow_mol_phase_comp[0, "Liq", "H2O"]
        ) == pytest.approx(0.03734, rel=1e-3)

        assert value(
            bped_m[1].fs.unit.outlet_acidic.flow_mol_phase_comp[0, "Liq", "H_+"]
        ) == pytest.approx(0.09233, rel=1e-3)
        assert value(
            bped_m[1].fs.unit.outlet_acidic.flow_mol_phase_comp[0, "Liq", "OH_-"]
        ) == pytest.approx(0.0006218, rel=1e-3)
        assert value(
            bped_m[1].fs.unit.outlet_acidic.flow_mol_phase_comp[0, "Liq", "Na_+"]
        ) == pytest.approx(0.003246, rel=1e-3)
        assert value(
            bped_m[1].fs.unit.outlet_acidic.flow_mol_phase_comp[0, "Liq", "Cl_-"]
        ) == pytest.approx(0.09385, rel=1e-3)
        assert value(
            bped_m[1].fs.unit.outlet_acidic.flow_mol_phase_comp[0, "Liq", "H2O"]
        ) == pytest.approx(0.4408, rel=1e-3)

        assert value(
            bped_m[1].fs.unit.outlet_basic.flow_mol_phase_comp[0, "Liq", "H_+"]
        ) == pytest.approx(0.0006218, rel=1e-3)
        assert value(
            bped_m[1].fs.unit.outlet_basic.flow_mol_phase_comp[0, "Liq", "OH_-"]
        ) == pytest.approx(0.092331, rel=1e-3)
        assert value(
            bped_m[1].fs.unit.outlet_basic.flow_mol_phase_comp[0, "Liq", "Na_+"]
        ) == pytest.approx(0.09325, rel=1e-3)
        assert value(
            bped_m[1].fs.unit.outlet_basic.flow_mol_phase_comp[0, "Liq", "Cl_-"]
        ) == pytest.approx(0.003216, rel=1e-3)
        assert value(
            bped_m[1].fs.unit.outlet_basic.flow_mol_phase_comp[0, "Liq", "H2O"]
        ) == pytest.approx(0.22267, rel=1e-3)

        # Test constant Voltage operation with salt calculation

        bped_m[2].fs.unit.voltage_applied.fix(1e1)
        bped_m[2].fs.unit.cell_triplet_num.fix(50)
        iscale.calculate_scaling_factors(bped_m[2])
        assert degrees_of_freedom(bped_m[2]) == 0
        initialization_tester(bped_m[2])
        assert value(
            bped_m[2].fs.unit.outlet_diluate.flow_mol_phase_comp[0, "Liq", "Na_+"]
        ) == pytest.approx(0.06089, rel=1e-3)
        assert value(
            bped_m[2].fs.unit.outlet_diluate.flow_mol_phase_comp[0, "Liq", "Cl_-"]
        ) == pytest.approx(0.06131, rel=1e-3)
        assert value(
            bped_m[2].fs.unit.outlet_diluate.flow_mol_phase_comp[0, "Liq", "H2O"]
        ) == pytest.approx(0.17156, rel=1e-3)

        assert value(
            bped_m[2].fs.unit.outlet_acidic.flow_mol_phase_comp[0, "Liq", "H_+"]
        ) == pytest.approx(0.07375, rel=1e-3)
        assert value(
            bped_m[2].fs.unit.outlet_acidic.flow_mol_phase_comp[0, "Liq", "Na_+"]
        ) == pytest.approx(0.012023, rel=1e-3)
        assert value(
            bped_m[2].fs.unit.outlet_acidic.flow_mol_phase_comp[0, "Liq", "Cl_-"]
        ) == pytest.approx(0.07416, rel=1e-3)
        assert value(
            bped_m[2].fs.unit.outlet_acidic.flow_mol_phase_comp[0, "Liq", "H2O"]
        ) == pytest.approx(0.27893, rel=1e-3)

        assert value(
            bped_m[2].fs.unit.outlet_basic.flow_mol_phase_comp[0, "Liq", "OH_-"]
        ) == pytest.approx(0.073755, rel=1e-3)
        assert value(
            bped_m[2].fs.unit.outlet_basic.flow_mol_phase_comp[0, "Liq", "Na_+"]
        ) == pytest.approx(0.07469, rel=1e-3)
        assert value(
            bped_m[2].fs.unit.outlet_basic.flow_mol_phase_comp[0, "Liq", "Cl_-"]
        ) == pytest.approx(0.01213, rel=1e-3)
        assert value(
            bped_m[2].fs.unit.outlet_basic.flow_mol_phase_comp[0, "Liq", "H2O"]
        ) == pytest.approx(0.2695, rel=1e-3)

        # Test constant Current operation with salt calculation

        bped_m[3].fs.unit.current_applied.fix(2e1)
        bped_m[3].fs.unit.cell_triplet_num.fix(50)
        iscale.set_scaling_factor(bped_m[3].fs.unit.voltage_x, 1e-0)
        iscale.calculate_scaling_factors(bped_m[3])
        assert degrees_of_freedom(bped_m[3]) == 0
        initialization_tester(bped_m[3])
        assert value(
            bped_m[3].fs.unit.outlet_diluate.flow_mol_phase_comp[0, "Liq", "Na_+"]
        ) == pytest.approx(0.05166, rel=1e-3)
        assert value(
            bped_m[3].fs.unit.outlet_diluate.flow_mol_phase_comp[0, "Liq", "Cl_-"]
        ) == pytest.approx(0.050728, rel=1e-3)
        assert value(
            bped_m[3].fs.unit.outlet_diluate.flow_mol_phase_comp[0, "Liq", "H2O"]
        ) == pytest.approx(0.111034, rel=1e-3)

        assert value(
            bped_m[3].fs.unit.outlet_acidic.flow_mol_phase_comp[0, "Liq", "H_+"]
        ) == pytest.approx(0.08209, rel=1e-3)
        assert value(
            bped_m[3].fs.unit.outlet_acidic.flow_mol_phase_comp[0, "Liq", "Na_+"]
        ) == pytest.approx(0.01318, rel=1e-3)
        assert value(
            bped_m[3].fs.unit.outlet_acidic.flow_mol_phase_comp[0, "Liq", "Cl_-"]
        ) == pytest.approx(0.083934, rel=1e-3)
        assert value(
            bped_m[3].fs.unit.outlet_acidic.flow_mol_phase_comp[0, "Liq", "H2O"]
        ) == pytest.approx(0.36229, rel=1e-3)

        assert value(
            bped_m[3].fs.unit.outlet_basic.flow_mol_phase_comp[0, "Liq", "H_+"]
        ) == pytest.approx(0.0003109, rel=1e-3)
        assert value(
            bped_m[3].fs.unit.outlet_basic.flow_mol_phase_comp[0, "Liq", "OH_-"]
        ) == pytest.approx(0.082093, rel=1e-3)
        assert value(
            bped_m[3].fs.unit.outlet_basic.flow_mol_phase_comp[0, "Liq", "Na_+"]
        ) == pytest.approx(0.0827538, rel=1e-3)
        assert value(
            bped_m[3].fs.unit.outlet_basic.flow_mol_phase_comp[0, "Liq", "Cl_-"]
        ) == pytest.approx(0.012938, rel=1e-3)
        assert value(
            bped_m[3].fs.unit.outlet_basic.flow_mol_phase_comp[0, "Liq", "H2O"]
        ) == pytest.approx(0.23807, rel=1e-3)


class Test_BPED_pressure_drop_components:

    @pytest.fixture(scope="class")
    def bped_m0(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        ion_dict = {
            "solute_list": ["Na_+", "Cl_-", "H_+", "OH_-"],
            "mw_data": {
                "Na_+": 23e-3,
                "Cl_-": 35.5e-3,
                "H_+": 1e-3,
                "OH_-": 17.0e-3,
            },
            "elec_mobility_data": {
                ("Liq", "Na_+"): 5.19e-8,
                ("Liq", "Cl_-"): 7.92e-8,
                ("Liq", "H_+"): 36.23e-8,
                ("Liq", "OH_-"): 20.64e-8,
            },
            "charge": {"Na_+": 1, "Cl_-": -1, "H_+": 1, "OH_-": -1},
            "diffusivity_data": {
                ("Liq", "Na_+"): 1.33e-9,
                ("Liq", "Cl_-"): 2.03e-9,
                ("Liq", "H_+"): 9.31e-9,
                ("Liq", "OH_-"): 5.27e-9,
            },
        }
        m.fs.properties = MCASParameterBlock(**ion_dict)
        m.fs.unit = Electrodialysis_Bipolar_1D(
            property_package=m.fs.properties,
            salt_calculation=True,
            limiting_current_density_method_bpem=LimitingCurrentDensitybpemMethod.Empirical,
            pressure_drop_method=PressureDropMethod.experimental,
            has_pressure_change=True,
        )
        return m

    @pytest.fixture(scope="class")
    def bped_m1(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        ion_dict = {
            "solute_list": ["Na_+", "Cl_-", "H_+", "OH_-"],
            "mw_data": {
                "Na_+": 23e-3,
                "Cl_-": 35.5e-3,
                "H_+": 1e-3,
                "OH_-": 17.0e-3,
            },
            "elec_mobility_data": {
                ("Liq", "Na_+"): 5.19e-8,
                ("Liq", "Cl_-"): 7.92e-8,
                ("Liq", "H_+"): 36.23e-8,
                ("Liq", "OH_-"): 20.64e-8,
            },
            "charge": {"Na_+": 1, "Cl_-": -1, "H_+": 1, "OH_-": -1},
            "diffusivity_data": {
                ("Liq", "Na_+"): 1.33e-9,
                ("Liq", "Cl_-"): 2.03e-9,
                ("Liq", "H_+"): 9.31e-9,
                ("Liq", "OH_-"): 5.27e-9,
            },
        }
        m.fs.properties = MCASParameterBlock(**ion_dict)
        m.fs.unit = Electrodialysis_Bipolar_1D(
            property_package=m.fs.properties,
            salt_calculation=True,
            limiting_current_density_method_bpem=LimitingCurrentDensitybpemMethod.Empirical,
            pressure_drop_method=PressureDropMethod.Darcy_Weisbach,
            friction_factor_method=FrictionFactorMethod.fixed,
            hydraulic_diameter_method=HydraulicDiameterMethod.conventional,
            has_pressure_change=True,
        )
        return m

    @pytest.fixture(scope="class")
    def bped_m2(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        ion_dict = {
            "solute_list": ["Na_+", "Cl_-", "H_+", "OH_-"],
            "mw_data": {
                "Na_+": 23e-3,
                "Cl_-": 35.5e-3,
                "H_+": 1e-3,
                "OH_-": 17.0e-3,
            },
            "elec_mobility_data": {
                ("Liq", "Na_+"): 5.19e-8,
                ("Liq", "Cl_-"): 7.92e-8,
                ("Liq", "H_+"): 36.23e-8,
                ("Liq", "OH_-"): 20.64e-8,
            },
            "charge": {"Na_+": 1, "Cl_-": -1, "H_+": 1, "OH_-": -1},
            "diffusivity_data": {
                ("Liq", "Na_+"): 1.33e-9,
                ("Liq", "Cl_-"): 2.03e-9,
                ("Liq", "H_+"): 9.31e-9,
                ("Liq", "OH_-"): 5.27e-9,
            },
        }
        m.fs.properties = MCASParameterBlock(**ion_dict)
        m.fs.unit = Electrodialysis_Bipolar_1D(
            property_package=m.fs.properties,
            salt_calculation=True,
            limiting_current_density_method_bpem=LimitingCurrentDensitybpemMethod.Empirical,
            pressure_drop_method=PressureDropMethod.Darcy_Weisbach,
            friction_factor_method=FrictionFactorMethod.Gurreri,
            hydraulic_diameter_method=HydraulicDiameterMethod.conventional,
            has_pressure_change=True,
        )
        return m

    @pytest.fixture(scope="class")
    def bped_m3(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        ion_dict = {
            "solute_list": ["Na_+", "Cl_-", "H_+", "OH_-"],
            "mw_data": {
                "Na_+": 23e-3,
                "Cl_-": 35.5e-3,
                "H_+": 1e-3,
                "OH_-": 17.0e-3,
            },
            "elec_mobility_data": {
                ("Liq", "Na_+"): 5.19e-8,
                ("Liq", "Cl_-"): 7.92e-8,
                ("Liq", "H_+"): 36.23e-8,
                ("Liq", "OH_-"): 20.64e-8,
            },
            "charge": {"Na_+": 1, "Cl_-": -1, "H_+": 1, "OH_-": -1},
            "diffusivity_data": {
                ("Liq", "Na_+"): 1.33e-9,
                ("Liq", "Cl_-"): 2.03e-9,
                ("Liq", "H_+"): 9.31e-9,
                ("Liq", "OH_-"): 5.27e-9,
            },
        }
        m.fs.properties = MCASParameterBlock(**ion_dict)
        m.fs.unit = Electrodialysis_Bipolar_1D(
            property_package=m.fs.properties,
            salt_calculation=True,
            limiting_current_density_method_bpem=LimitingCurrentDensitybpemMethod.Empirical,
            pressure_drop_method=PressureDropMethod.Darcy_Weisbach,
            friction_factor_method=FrictionFactorMethod.Kuroda,
            hydraulic_diameter_method=HydraulicDiameterMethod.conventional,
            has_pressure_change=True,
        )
        return m

    @pytest.fixture(scope="class")
    def bped_m4(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        ion_dict = {
            "solute_list": ["Na_+", "Cl_-", "H_+", "OH_-"],
            "mw_data": {
                "Na_+": 23e-3,
                "Cl_-": 35.5e-3,
                "H_+": 1e-3,
                "OH_-": 17.0e-3,
            },
            "elec_mobility_data": {
                ("Liq", "Na_+"): 5.19e-8,
                ("Liq", "Cl_-"): 7.92e-8,
                ("Liq", "H_+"): 36.23e-8,
                ("Liq", "OH_-"): 20.64e-8,
            },
            "charge": {"Na_+": 1, "Cl_-": -1, "H_+": 1, "OH_-": -1},
            "diffusivity_data": {
                ("Liq", "Na_+"): 1.33e-9,
                ("Liq", "Cl_-"): 2.03e-9,
                ("Liq", "H_+"): 9.31e-9,
                ("Liq", "OH_-"): 5.27e-9,
            },
        }
        m.fs.properties = MCASParameterBlock(**ion_dict)
        m.fs.unit = Electrodialysis_Bipolar_1D(
            property_package=m.fs.properties,
            salt_calculation=True,
            limiting_current_density_method_bpem=LimitingCurrentDensitybpemMethod.Empirical,
            pressure_drop_method=PressureDropMethod.Darcy_Weisbach,
            friction_factor_method=FrictionFactorMethod.Kuroda,
            hydraulic_diameter_method=HydraulicDiameterMethod.fixed,
            has_pressure_change=True,
        )
        return m

    @pytest.fixture(scope="class")
    def bped_m5(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        ion_dict = {
            "solute_list": ["Na_+", "Cl_-", "H_+", "OH_-"],
            "mw_data": {
                "Na_+": 23e-3,
                "Cl_-": 35.5e-3,
                "H_+": 1e-3,
                "OH_-": 17.0e-3,
            },
            "elec_mobility_data": {
                ("Liq", "Na_+"): 5.19e-8,
                ("Liq", "Cl_-"): 7.92e-8,
                ("Liq", "H_+"): 36.23e-8,
                ("Liq", "OH_-"): 20.64e-8,
            },
            "charge": {"Na_+": 1, "Cl_-": -1, "H_+": 1, "OH_-": -1},
            "diffusivity_data": {
                ("Liq", "Na_+"): 1.33e-9,
                ("Liq", "Cl_-"): 2.03e-9,
                ("Liq", "H_+"): 9.31e-9,
                ("Liq", "OH_-"): 5.27e-9,
            },
        }
        m.fs.properties = MCASParameterBlock(**ion_dict)
        m.fs.unit = Electrodialysis_Bipolar_1D(
            property_package=m.fs.properties,
            salt_calculation=True,
            limiting_current_density_method_bpem=LimitingCurrentDensitybpemMethod.Empirical,
            pressure_drop_method=PressureDropMethod.Darcy_Weisbach,
            friction_factor_method=FrictionFactorMethod.Kuroda,
            hydraulic_diameter_method=HydraulicDiameterMethod.spacer_specific_area_known,
            has_pressure_change=True,
        )
        return m

    @pytest.mark.unit
    def test_deltaP_various_methods(
        self, bped_m0, bped_m1, bped_m2, bped_m3, bped_m4, bped_m5
    ):
        bped_m = (bped_m0, bped_m1, bped_m2, bped_m3, bped_m4, bped_m5)
        for m in bped_m:
            m.fs.unit.diffus_mass.fix((2.03 + 1.96) * 10**-9 / 2)
            m.fs.unit.membrane_fixed_charge.fix(5e3)
            m.fs.unit.salt_conc_ael_ref.fix(2000)
            m.fs.unit.salt_conc_cel_ref.fix(2000)
            m.fs.unit.conc_water.fix(50 * 1e3)
            m.fs.unit.k2_zero.fix(2 * 10**-6)
            m.fs.unit.relative_permittivity.fix(30)
            m.fs.unit.membrane_fixed_catalyst_cel.fix(5e3)
            m.fs.unit.membrane_fixed_catalyst_ael.fix(5e3)
            m.fs.unit.k_a.fix(447)
            m.fs.unit.k_b.fix(5e4)

            m.fs.unit.ion_trans_number_membrane["bpem", "Na_+"].fix(0)
            m.fs.unit.ion_trans_number_membrane["bpem", "Cl_-"].fix(0)
            m.fs.unit.ion_trans_number_membrane["bpem", "H_+"].fix(1)
            m.fs.unit.ion_trans_number_membrane["bpem", "OH_-"].fix(1)

            m.fs.unit.ion_trans_number_membrane["cem", "Na_+"].fix(0.94)
            m.fs.unit.ion_trans_number_membrane["cem", "Cl_-"].fix(0)
            m.fs.unit.ion_trans_number_membrane["cem", "H_+"].fix(0.03)
            m.fs.unit.ion_trans_number_membrane["cem", "OH_-"].fix(0.03)

            m.fs.unit.ion_trans_number_membrane["aem", "Na_+"].fix(0)
            m.fs.unit.ion_trans_number_membrane["aem", "Cl_-"].fix(0.94)
            m.fs.unit.ion_trans_number_membrane["aem", "H_+"].fix(0.03)
            m.fs.unit.ion_trans_number_membrane["aem", "OH_-"].fix(0.03)

            m.fs.unit.solute_diffusivity_membrane["bpem", "Na_+"].fix(0)
            m.fs.unit.solute_diffusivity_membrane["bpem", "Cl_-"].fix(0)
            m.fs.unit.solute_diffusivity_membrane["bpem", "H_+"].fix(0)
            m.fs.unit.solute_diffusivity_membrane["bpem", "OH_-"].fix(0)

            m.fs.unit.solute_diffusivity_membrane["cem", "Na_+"].fix(
                (1.8e-10 + 1.25e-10) / 2
            )
            m.fs.unit.solute_diffusivity_membrane["cem", "Cl_-"].fix(
                (1.8e-10 + 1.25e-10) / 2
            )
            m.fs.unit.solute_diffusivity_membrane["cem", "H_+"].fix(0)
            m.fs.unit.solute_diffusivity_membrane["cem", "OH_-"].fix(0)

            m.fs.unit.solute_diffusivity_membrane["aem", "Na_+"].fix(
                (1.8e-10 + 1.25e-10) / 2
            )
            m.fs.unit.solute_diffusivity_membrane["aem", "Cl_-"].fix(
                (1.8e-10 + 1.25e-10) / 2
            )
            m.fs.unit.solute_diffusivity_membrane["aem", "H_+"].fix(0)
            m.fs.unit.solute_diffusivity_membrane["aem", "OH_-"].fix(0)

            # Set inlet streams.
            m.fs.unit.inlet_basic.pressure.fix(101325 * 3)
            m.fs.unit.inlet_basic.temperature.fix(298.15)
            m.fs.unit.inlet_acidic.pressure.fix(101325 * 3)
            m.fs.unit.inlet_acidic.temperature.fix(298.15)
            m.fs.unit.inlet_diluate.pressure.fix(101325 * 3)
            m.fs.unit.inlet_diluate.temperature.fix(298.15)
            m.fs.unit.spacer_porosity.fix(1)

            m.fs.unit.inlet_basic.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(7.38e-2)
            m.fs.unit.inlet_basic.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(0)
            m.fs.unit.inlet_basic.flow_mol_phase_comp[0, "Liq", "H_+"].fix(0)
            m.fs.unit.inlet_basic.flow_mol_phase_comp[0, "Liq", "OH_-"].fix(7.38e-2)
            m.fs.unit.inlet_acidic.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(0)
            m.fs.unit.inlet_acidic.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(7.38e-2)
            m.fs.unit.inlet_acidic.flow_mol_phase_comp[0, "Liq", "H_+"].fix(7.38e-2)
            m.fs.unit.inlet_acidic.flow_mol_phase_comp[0, "Liq", "OH_-"].fix(0)
            m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(7.38e-2)
            m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(7.38e-2)
            m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "H_+"].fix(0)
            m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "OH_-"].fix(0)

            m.fs.unit.inlet_basic.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.40e-1)
            m.fs.unit.inlet_acidic.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.40e-1)
            m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.40e-1)

            m.fs.unit.shadow_factor.fix(1)

            m.fs.unit.water_trans_number_membrane["cem"].fix((5.8 + 4.3) / 2)
            m.fs.unit.water_permeability_membrane["cem"].fix((2.16e-14 + 1.75e-14) / 2)

            m.fs.unit.water_trans_number_membrane["aem"].fix((5.8 + 4.3) / 2)
            m.fs.unit.water_permeability_membrane["aem"].fix((2.16e-14 + 1.75e-14) / 2)

            m.fs.unit.water_trans_number_membrane["bpem"].fix((5.8 + 4.3) / 2)
            m.fs.unit.water_permeability_membrane["bpem"].fix((2.16e-14 + 1.75e-14) / 2)
            m.fs.unit.electrodes_resistance.fix(0)
            m.fs.unit.current_utilization.fix(1)
            m.fs.unit.channel_height.fix(2.7e-4)
            m.fs.unit.membrane_areal_resistance_coef_0.fix((1.89e-4 + 1.77e-4) / 2)
            m.fs.unit.membrane_areal_resistance_coef_1.fix(0)
            m.fs.unit.cell_width.fix(0.1)
            m.fs.unit.cell_length.fix(0.79)
            m.fs.unit.membrane_thickness["bpem"].fix(8e-4)
            m.fs.unit.membrane_thickness["aem"].fix(4e-4)
            m.fs.unit.membrane_thickness["cem"].fix(4e-4)

            # Set scaling of critical quantities.
            m.fs.properties.set_default_scaling(
                "flow_mol_phase_comp", 1e3, index=("Liq", "Na_+")
            )
            m.fs.properties.set_default_scaling(
                "flow_mol_phase_comp", 1e3, index=("Liq", "Cl_-")
            )
            m.fs.properties.set_default_scaling(
                "flow_mol_phase_comp", 1e3, index=("Liq", "H_+")
            )
            m.fs.properties.set_default_scaling(
                "flow_mol_phase_comp", 1e3, index=("Liq", "OH_-")
            )
            m.fs.properties.set_default_scaling(
                "flow_mol_phase_comp", 1e2, index=("Liq", "H2O")
            )
            m.fs.unit.current_applied.fix(2e1)
            m.fs.unit.cell_triplet_num.fix(50)
            iscale.set_scaling_factor(m.fs.unit.k_a, 1e-2)
            iscale.set_scaling_factor(m.fs.unit.k_b, 1e-4)
            iscale.set_scaling_factor(m.fs.unit.voltage_x, 1e-0)
            iscale.set_scaling_factor(m.fs.unit.flux_splitting, 1e4)

        # Test bped_m0
        bped_m[0].fs.unit.pressure_drop.fix(4e4)
        iscale.calculate_scaling_factors(bped_m[0])
        assert degrees_of_freedom(bped_m[0]) == 0
        initialization_tester(bped_m[0])
        results = solver.solve(bped_m[0])
        assert_optimal_termination(results)
        assert value(bped_m[0].fs.unit.pressure_drop_total[0]) == pytest.approx(
            31600, rel=1e-3
        )

        # Test bped_m1
        bped_m[1].fs.unit.friction_factor.fix(20)
        iscale.calculate_scaling_factors(bped_m[1])
        assert degrees_of_freedom(bped_m[1]) == 0
        initialization_tester(bped_m[1])
        results = solver.solve(bped_m[1])
        assert_optimal_termination(results)
        assert value(bped_m[1].fs.unit.N_Re) == pytest.approx(3.4456, rel=1e-3)

        assert value(bped_m[1].fs.unit.pressure_drop[0]) == pytest.approx(
            760.09, rel=1e-3
        )

        assert value(bped_m[1].fs.unit.pressure_drop_total[0]) == pytest.approx(
            600.472, rel=1e-3
        )

        # Test bped_m2
        iscale.calculate_scaling_factors(bped_m[2])
        assert degrees_of_freedom(bped_m[2]) == 0
        initialization_tester(bped_m[2])
        results = solver.solve(bped_m[2])
        assert_optimal_termination(results)
        assert value(bped_m[2].fs.unit.N_Re) == pytest.approx(3.446, rel=1e-3)

        assert value(bped_m[2].fs.unit.pressure_drop[0]) == pytest.approx(
            2232.437, rel=1e-3
        )

        assert value(bped_m[2].fs.unit.pressure_drop_total[0]) == pytest.approx(
            1763.625, rel=1e-3
        )

        # Test bped_m3
        iscale.calculate_scaling_factors(bped_m[3])
        assert degrees_of_freedom(bped_m[3]) == 0
        initialization_tester(bped_m[3])
        results = solver.solve(bped_m[3])
        assert_optimal_termination(results)
        assert value(bped_m[3].fs.unit.N_Re) == pytest.approx(3.446, rel=1e-3)

        assert value(bped_m[3].fs.unit.pressure_drop[0]) == pytest.approx(
            786.201, rel=1e-3
        )

        assert value(bped_m[3].fs.unit.pressure_drop_total[0]) == pytest.approx(
            621.099, rel=1e-3
        )

        # Test bped_m4
        bped_m[4].fs.unit.hydraulic_diameter.fix(1e-3)
        iscale.calculate_scaling_factors(bped_m[4])
        assert degrees_of_freedom(bped_m[4]) == 0
        initialization_tester(bped_m[4])
        results = solver.solve(bped_m[4])
        assert_optimal_termination(results)
        assert value(bped_m[4].fs.unit.N_Re) == pytest.approx(6.398, rel=1e-3)

        assert value(bped_m[4].fs.unit.pressure_drop[0]) == pytest.approx(
            310.719, rel=1e-3
        )

        assert value(bped_m[4].fs.unit.pressure_drop_total[0]) == pytest.approx(
            245.468, rel=1e-3
        )

        # Test bped_m5
        bped_m[5].fs.unit.spacer_specific_area.fix(10700)
        iscale.calculate_scaling_factors(bped_m[5])
        assert degrees_of_freedom(bped_m[5]) == 0
        initialization_tester(bped_m[5])
        iscale.calculate_scaling_factors(bped_m[5])
        results = solver.solve(bped_m[5])
        assert_optimal_termination(results)
        assert value(bped_m[5].fs.unit.N_Re) == pytest.approx(3.455, rel=1e-3)

        assert value(bped_m[5].fs.unit.pressure_drop[0]) == pytest.approx(
            783.027, rel=1e-3
        )

        assert value(bped_m[5].fs.unit.pressure_drop_total[0]) == pytest.approx(
            618.592, rel=1e-3
        )
