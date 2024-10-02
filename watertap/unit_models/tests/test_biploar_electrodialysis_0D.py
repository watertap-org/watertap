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
import idaes.logger as idaeslog
import pytest
from idaes.core import (
    FlowsheetBlock,
)
from idaes.core.util.constants import Constants
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from pyomo.environ import (
    ConcreteModel,
    assert_optimal_termination,
    value,
)

from watertap.unit_models.Bipolar_Electrodialysis_0D import (
    Bipolar_Electrodialysis_0D,
    LimitingCurrentDensitybpemMethod,
    LimitingpotentialMethod,
    PressureDropMethod,
    FrictionFactorMethod,
    HydraulicDiameterMethod,
)
from watertap.core.solvers import get_solver
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock

# from watertap.unit_models.electrodialysis_0D import LimitingCurrentDensitybpemMethod

__author__ = "Johnson Dhanasekaran"

solver = get_solver()


# -----------------------------------------------------------------------------
# Start test class


class Test_catalyst:
    @pytest.fixture(scope="class")
    def bped(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        ion_dict = {
            "solute_list": ["Na_+", "Cl_-", "H_+", "OH_-"],
            "mw_data": {
                "H2O": 18e-3,
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
        m.fs.unit = Bipolar_Electrodialysis_0D(
            property_package=m.fs.properties,
            has_catalyst=True,
            limiting_current_density_method_bpem=LimitingCurrentDensitybpemMethod.Empirical,
        )

        m.fs.unit.diffus_mass.fix((2.03 + 1.96) * 10**-9 / 2)
        m.fs.unit.membrane_fixed_charge["bpem"].fix(5e3)
        m.fs.unit.salt_conc_aem["bpem"].fix(2000)
        m.fs.unit.salt_conc_cem["bpem"].fix(2000)
        m.fs.unit.conc_water["bpem"].fix(50 * 1e3)
        m.fs.unit.kr["bpem"].fix(1.3 * 10**10)
        m.fs.unit.k2_zero["bpem"].fix(2 * 10**-6)
        m.fs.unit.relative_permittivity["bpem"].fix(30)
        m.fs.unit.membrane_fixed_catalyst_cem["bpem"].fix(5e3)
        m.fs.unit.membrane_fixed_catalyst_aem["bpem"].fix(5e3)
        m.fs.unit.k_a.fix(447)
        m.fs.unit.k_b.fix(5e4)

        m.fs.unit.ion_trans_number_membrane["bpem", "Na_+"].fix(0.5)
        m.fs.unit.ion_trans_number_membrane["bpem", "Cl_-"].fix(0.5)
        m.fs.unit.ion_trans_number_membrane["bpem", "H_+"].fix(0.1)
        m.fs.unit.ion_trans_number_membrane["bpem", "OH_-"].fix(0.1)

        # Set inlet streams.
        m.fs.unit.inlet_basate.pressure.fix(101325)
        m.fs.unit.inlet_basate.temperature.fix(298.15)
        m.fs.unit.inlet_acidate.pressure.fix(101325)
        m.fs.unit.inlet_acidate.temperature.fix(298.15)
        m.fs.unit.spacer_porosity.fix(1)

        m.fs.unit.inlet_basate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(7.38e-4)
        m.fs.unit.inlet_basate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(7.38e-4)
        m.fs.unit.inlet_basate.flow_mol_phase_comp[0, "Liq", "H_+"].fix(7.38e-4)
        m.fs.unit.inlet_basate.flow_mol_phase_comp[0, "Liq", "OH_-"].fix(7.38e-4)
        m.fs.unit.inlet_acidate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(7.38e-4)
        m.fs.unit.inlet_acidate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(7.38e-4)
        m.fs.unit.inlet_acidate.flow_mol_phase_comp[0, "Liq", "H_+"].fix(7.38e-4)
        m.fs.unit.inlet_acidate.flow_mol_phase_comp[0, "Liq", "OH_-"].fix(7.38e-4)

        m.fs.unit.inlet_basate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.40e-1)
        m.fs.unit.inlet_acidate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.40e-1)

        m.fs.unit.shadow_factor.fix(1)
        m.fs.unit.water_trans_number_membrane["bpem"].fix((5.8 + 4.3) / 2)
        m.fs.unit.water_permeability_membrane["bpem"].fix((2.16e-14 + 1.75e-14) / 2)
        m.fs.unit.electrodes_resistance.fix(0)
        m.fs.unit.current_utilization.fix(1)
        m.fs.unit.channel_height.fix(2.7e-4)
        m.fs.unit.membrane_areal_resistance.fix((1.89e-4 + 1.77e-4) / 2)
        m.fs.unit.cell_width.fix(0.1)
        m.fs.unit.cell_length.fix(0.79)
        m.fs.unit.membrane_thickness["bpem"].fix(8e-4)

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
        m.fs.unit.cell_num.fix(1)
        iscale.set_scaling_factor(m.fs.unit.k_a, 1e-2)
        iscale.set_scaling_factor(m.fs.unit.k_b, 1e-4)
        iscale.set_scaling_factor(m.fs.unit.voltage, 1e-1)
        iscale.set_scaling_factor(m.fs.unit.flux_splitting, 1e3)
        return m

    @pytest.mark.unit
    def test_assign(self, bped):

        # Experimental data from Wilhelm et al. (2002) with additional inputs from Mareev et al. (2020)
        expt_current_density = np.array(
            [11.7, 15.5, 20.8, 24.4, 30.6, 40.0, 100.6]
        )  # in mA/cm2
        expt_membrane_potential = np.array(
            [0.088, 0.184, 0.4045, 0.690, 0.858, 0.95, 1.3]
        )  # in volts

        for indx, v in enumerate(expt_membrane_potential):
            m = bped
            m.fs.unit.potential_membrane_bpem[0].fix(v)

            iscale.calculate_scaling_factors(m.fs)
            assert degrees_of_freedom(m) == 0
            initialization_tester(m, outlvl=idaeslog.DEBUG)
            results = solver.solve(m)
            assert_optimal_termination(results)
            if indx < 2:
                rel_tol = 1e0
            else:
                rel_tol = 1.5e-1

            current_density_computed = (
                0.1
                * m.fs.unit.current[0]
                / (m.fs.unit.cell_width * m.fs.unit.cell_length)
            )  # convert to mA/cm2
            assert value(current_density_computed) == pytest.approx(
                expt_current_density[indx], rel=rel_tol
            )


class Test_limiting_parameters:
    @pytest.fixture(scope="class")
    def limiting_current_check(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        ion_dict = {
            "solute_list": ["Na_+", "Cl_-", "H_+", "OH_-"],
            "mw_data": {
                "H2O": 18e-3,
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
        m.fs.unit = Bipolar_Electrodialysis_0D(
            property_package=m.fs.properties,
            limiting_current_density_method_bpem=LimitingCurrentDensitybpemMethod.Empirical,
            limiting_potential_method_bpem=LimitingpotentialMethod.InitialValue,
            limiting_potential_data=0.5,
            has_catalyst=False,
        )
        return m

    @pytest.fixture(scope="class")
    def potential_barrier_check(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        ion_dict = {
            "solute_list": ["Na_+", "Cl_-", "H_+", "OH_-"],
            "mw_data": {
                "H2O": 18e-3,
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
        m.fs.unit = Bipolar_Electrodialysis_0D(
            property_package=m.fs.properties,
            limiting_current_density_method_bpem=LimitingCurrentDensitybpemMethod.Empirical,
            limiting_potential_method_bpem=LimitingpotentialMethod.Empirical,
            has_catalyst=False,
        )

        m.fs.unit.kr["bpem"].fix(1.33 * 10**11)
        m.fs.unit.k2_zero["bpem"].fix(2 * 10**-5)
        m.fs.unit.relative_permittivity["bpem"].fix(20)

        return m

    @pytest.mark.unit
    def test_assign(self, limiting_current_check, potential_barrier_check):
        check_m = (limiting_current_check, potential_barrier_check)
        for m in check_m:
            m.fs.unit.shadow_factor.fix(1)
            m.fs.unit.current.fix(1e1)
            m.fs.unit.water_trans_number_membrane["bpem"].fix((5.8 + 4.3) / 2)
            m.fs.unit.water_permeability_membrane["bpem"].fix((2.16e-14 + 1.75e-14) / 2)
            m.fs.unit.electrodes_resistance.fix(0)
            m.fs.unit.cell_num.fix(1)
            m.fs.unit.current_utilization.fix(1)
            m.fs.unit.channel_height.fix(2.7e-4)
            m.fs.unit.membrane_areal_resistance.fix((1.89e-4 + 1.77e-4) / 2)
            m.fs.unit.cell_width.fix(0.1)
            m.fs.unit.cell_length.fix(0.79)
            m.fs.unit.membrane_thickness["bpem"].fix(1.3e-4)
            m.fs.unit.diffus_mass.fix((2.03 + 1.96) * 10**-9 / 2)
            m.fs.unit.ion_trans_number_membrane["bpem", "Na_+"].fix(0.5)
            m.fs.unit.ion_trans_number_membrane["bpem", "Cl_-"].fix(0.5)
            m.fs.unit.ion_trans_number_membrane["bpem", "H_+"].fix(0.1)
            m.fs.unit.ion_trans_number_membrane["bpem", "OH_-"].fix(0.1)
            m.fs.unit.conc_water["bpem"].fix(55 * 1e3)
            m.fs.unit.kr["bpem"].fix(1.33 * 10**11)
            m.fs.unit.k2_zero["bpem"].fix(2 * 10**-5)
            m.fs.unit.relative_permittivity["bpem"].fix(20)
            m.fs.unit.diffus_mass.fix((2.03 + 1.96) * 10**-9 / 2)
            m.fs.unit.salt_conc_aem["bpem"].fix(500 + 2 * 250)
            m.fs.unit.salt_conc_cem["bpem"].fix(500 + 2 * 250)
            m.fs.unit.membrane_fixed_charge["bpem"].fix(1.5e3)

            # Set inlet streams.
            m.fs.unit.inlet_basate.pressure.fix(101325)
            m.fs.unit.inlet_basate.temperature.fix(298.15)
            m.fs.unit.inlet_basate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.40e-1)
            m.fs.unit.inlet_basate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(7.38e-1)
            m.fs.unit.inlet_basate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(7.38e-1)
            m.fs.unit.inlet_basate.flow_mol_phase_comp[0, "Liq", "H_+"].fix(7.38e-1)
            m.fs.unit.inlet_basate.flow_mol_phase_comp[0, "Liq", "OH_-"].fix(7.38e-1)
            m.fs.unit.inlet_acidate.pressure.fix(101325)
            m.fs.unit.inlet_acidate.temperature.fix(298.15)
            m.fs.unit.inlet_acidate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.40e-1)
            m.fs.unit.inlet_acidate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(7.38e-1)
            m.fs.unit.inlet_acidate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(7.38e-1)
            m.fs.unit.inlet_acidate.flow_mol_phase_comp[0, "Liq", "H_+"].fix(7.38e-1)
            m.fs.unit.inlet_acidate.flow_mol_phase_comp[0, "Liq", "OH_-"].fix(7.38e-1)
            m.fs.unit.spacer_porosity.fix(1)

            # Set scaling of critical quantities.
            m.fs.properties.set_default_scaling(
                "flow_mol_phase_comp", 1e1, index=("Liq", "H2O")
            )
            m.fs.properties.set_default_scaling(
                "flow_mol_phase_comp", 1e0, index=("Liq", "Na_+")
            )
            m.fs.properties.set_default_scaling(
                "flow_mol_phase_comp", 1e0, index=("Liq", "Cl_-")
            )
            m.fs.properties.set_default_scaling(
                "flow_mol_phase_comp", 1e0, index=("Liq", "H_+")
            )
            m.fs.properties.set_default_scaling(
                "flow_mol_phase_comp", 1e0, index=("Liq", "OH_-")
            )

        # Data on limiting current and potential barrier to water splitting have been obtained from:
        # Fumatech, Technical Data Sheet for Fumasep FBM, 2020. With additional modelling parameters obtained from
        # Ionescu, Viorel. Advanced Topics in Optoelectronics, Microelectronics, and Nanotechnologies (2023)
        iscale.calculate_scaling_factors(check_m[1])

        # Test computing limiting current in  bipolar membrane
        iscale.calculate_scaling_factors(check_m[0])
        assert degrees_of_freedom(check_m[0]) == 0
        initialization_tester(check_m[0], outlvl=idaeslog.DEBUG)
        results = solver.solve(check_m[0])
        assert_optimal_termination(results)
        assert value(check_m[0].fs.unit.current_dens_lim_bpem[0]) == pytest.approx(
            1000, rel=1e-1
        )

        # Test computing limiting current and potential barrier to water splitting in  bipolar membrane
        iscale.calculate_scaling_factors(check_m[1])
        assert degrees_of_freedom(check_m[1]) == 0
        initialization_tester(check_m[1], outlvl=idaeslog.DEBUG)
        results = solver.solve(check_m[1])
        assert_optimal_termination(results)
        assert value(check_m[0].fs.unit.current_dens_lim_bpem[0]) == pytest.approx(
            1000, rel=1e-1
        )
        assert value(check_m[1].fs.unit.potential_barrier_bpem[0]) == pytest.approx(
            0.8, rel=1e-1
        )


class Test_BPED_pressure_drop_components:
    @pytest.fixture(scope="class")
    def bped_m0(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        ion_dict = {
            "solute_list": ["Na_+", "Cl_-", "H_+", "OH_-"],
            "mw_data": {
                "H2O": 18e-3,
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
        m.fs.unit = Bipolar_Electrodialysis_0D(
            property_package=m.fs.properties,
            has_catalyst=False,
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
                "H2O": 18e-3,
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
        m.fs.unit = Bipolar_Electrodialysis_0D(
            property_package=m.fs.properties,
            has_catalyst=False,
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
                "H2O": 18e-3,
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
        m.fs.unit = Bipolar_Electrodialysis_0D(
            property_package=m.fs.properties,
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
                "H2O": 18e-3,
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
        m.fs.unit = Bipolar_Electrodialysis_0D(
            property_package=m.fs.properties,
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
                "H2O": 18e-3,
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
        m.fs.unit = Bipolar_Electrodialysis_0D(
            property_package=m.fs.properties,
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
                "H2O": 18e-3,
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
        m.fs.unit = Bipolar_Electrodialysis_0D(
            property_package=m.fs.properties,
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
            m.fs.unit.shadow_factor.fix(1)
            m.fs.unit.current.fix(1e2)
            m.fs.unit.water_trans_number_membrane["bpem"].fix((5.8 + 4.3) / 2)
            m.fs.unit.water_permeability_membrane["bpem"].fix((2.16e-14 + 1.75e-14) / 2)
            m.fs.unit.electrodes_resistance.fix(0)
            m.fs.unit.cell_num.fix(10)
            m.fs.unit.current_utilization.fix(1)
            m.fs.unit.channel_height.fix(2.7e-4)
            m.fs.unit.membrane_areal_resistance.fix((1.89e-4 + 1.77e-4) / 2)
            m.fs.unit.cell_width.fix(0.1)
            m.fs.unit.cell_length.fix(0.79)
            m.fs.unit.membrane_thickness["bpem"].fix(1.3e-4)
            m.fs.unit.diffus_mass.fix((2.03 + 1.96) * 10**-9 / 2)
            m.fs.unit.ion_trans_number_membrane["bpem", "Na_+"].fix(0.5)
            m.fs.unit.ion_trans_number_membrane["bpem", "Cl_-"].fix(0.5)
            m.fs.unit.ion_trans_number_membrane["bpem", "H_+"].fix(0.1)
            m.fs.unit.ion_trans_number_membrane["bpem", "OH_-"].fix(0.1)
            m.fs.unit.conc_water["bpem"].fix(55 * 1e3)
            m.fs.unit.kr["bpem"].fix(1.33 * 10**11)
            m.fs.unit.k2_zero["bpem"].fix(2 * 10**-5)
            m.fs.unit.relative_permittivity["bpem"].fix(20)
            m.fs.unit.diffus_mass.fix((2.03 + 1.96) * 10**-9 / 2)
            m.fs.unit.salt_conc_aem["bpem"].fix(500 + 2 * 250)
            m.fs.unit.salt_conc_cem["bpem"].fix(500 + 2 * 250)
            m.fs.unit.membrane_fixed_charge["bpem"].fix(1.5e3)

            # Set inlet streams.
            m.fs.unit.inlet_basate.pressure.fix(201035)
            m.fs.unit.inlet_basate.temperature.fix(298.15)
            m.fs.unit.inlet_basate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.40e-1)
            m.fs.unit.inlet_basate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(7.38e-4)
            m.fs.unit.inlet_basate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(7.38e-4)
            m.fs.unit.inlet_basate.flow_mol_phase_comp[0, "Liq", "H_+"].fix(7.38e-4)
            m.fs.unit.inlet_basate.flow_mol_phase_comp[0, "Liq", "OH_-"].fix(7.38e-4)
            m.fs.unit.inlet_acidate.pressure.fix(201035)
            m.fs.unit.inlet_acidate.temperature.fix(298.15)
            m.fs.unit.inlet_acidate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.40e-1)
            m.fs.unit.inlet_acidate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(7.38e-4)
            m.fs.unit.inlet_acidate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(7.38e-4)
            m.fs.unit.inlet_acidate.flow_mol_phase_comp[0, "Liq", "H_+"].fix(7.38e-4)
            m.fs.unit.inlet_acidate.flow_mol_phase_comp[0, "Liq", "OH_-"].fix(7.38e-4)
            m.fs.unit.spacer_porosity.fix(0.83)

            # Set scaling of critical quantities.
            m.fs.properties.set_default_scaling(
                "flow_mol_phase_comp", 1e1, index=("Liq", "H2O")
            )
            m.fs.properties.set_default_scaling(
                "flow_mol_phase_comp", 1e0, index=("Liq", "Na_+")
            )
            m.fs.properties.set_default_scaling(
                "flow_mol_phase_comp", 1e0, index=("Liq", "Cl_-")
            )
            m.fs.properties.set_default_scaling(
                "flow_mol_phase_comp", 1e0, index=("Liq", "H_+")
            )
            m.fs.properties.set_default_scaling(
                "flow_mol_phase_comp", 1e0, index=("Liq", "OH_-")
            )

            iscale.set_scaling_factor(m.fs.unit.cell_length, 1)
            iscale.set_scaling_factor(m.fs.unit.cell_num, 0.1)
            iscale.set_scaling_factor(m.fs.unit.voltage, 1e-1)

        # Test bped_m0
        bped_m[0].fs.unit.pressure_drop.fix(1e4)
        iscale.calculate_scaling_factors(bped_m[0])
        assert degrees_of_freedom(bped_m[0]) == 0
        initialization_tester(bped_m[0], outlvl=idaeslog.DEBUG)
        results = solver.solve(bped_m[0])
        assert_optimal_termination(results)
        assert value(bped_m[0].fs.unit.pressure_drop_total[0]) == pytest.approx(
            7900, rel=1e-3
        )

        # Test bped_m1
        bped_m[1].fs.unit.current.fix(1e2)
        iscale.set_scaling_factor(bped_m[1].fs.unit.current, 1e-2)
        bped_m[1].fs.unit.diffus_mass.fix(1.6e-9)
        bped_m[1].fs.unit.friction_factor.fix(20)
        iscale.calculate_scaling_factors(bped_m[1])
        assert degrees_of_freedom(bped_m[1]) == 0
        initialization_tester(bped_m[1], outlvl=idaeslog.DEBUG)
        results = solver.solve(bped_m[1])
        assert_optimal_termination(results)
        assert value(bped_m[1].fs.unit.N_Re) == pytest.approx(9.491, rel=1e-3)

        assert value(bped_m[1].fs.unit.pressure_drop[0]) == pytest.approx(
            10087.548, rel=1e-3
        )

        assert value(bped_m[1].fs.unit.pressure_drop_total[0]) == pytest.approx(
            7969.163, rel=1e-3
        )

        # Test bped_m2
        bped_m[2].fs.unit.diffus_mass.fix(1.6e-9)
        iscale.calculate_scaling_factors(bped_m[2])
        assert degrees_of_freedom(bped_m[2]) == 0
        initialization_tester(bped_m[2], outlvl=idaeslog.DEBUG)
        results = solver.solve(bped_m[2])
        assert_optimal_termination(results)
        assert value(bped_m[2].fs.unit.N_Re) == pytest.approx(9.491, rel=1e-3)

        assert value(bped_m[2].fs.unit.pressure_drop[0]) == pytest.approx(
            40080.279, rel=1e-3
        )

        assert value(bped_m[2].fs.unit.pressure_drop_total[0]) == pytest.approx(
            31663.421, rel=1e-3
        )

        # Test bped_m3
        bped_m[3].fs.unit.diffus_mass.fix(1.6e-9)
        iscale.calculate_scaling_factors(bped_m[3])
        assert degrees_of_freedom(bped_m[3]) == 0
        initialization_tester(bped_m[3], outlvl=idaeslog.DEBUG)
        results = solver.solve(bped_m[3])
        assert_optimal_termination(results)
        assert value(bped_m[3].fs.unit.N_Re) == pytest.approx(9.491, rel=1e-3)

        assert value(bped_m[3].fs.unit.pressure_drop[0]) == pytest.approx(
            7574.199, rel=1e-3
        )

        assert value(bped_m[3].fs.unit.pressure_drop_total[0]) == pytest.approx(
            5983.618, rel=1e-3
        )

        # Test bped_m4
        bped_m[4].fs.unit.diffus_mass.fix(1.6e-9)
        bped_m[4].fs.unit.hydraulic_diameter.fix(1e-3)
        iscale.calculate_scaling_factors(bped_m[4])
        assert degrees_of_freedom(bped_m[4]) == 0
        initialization_tester(bped_m[4], outlvl=idaeslog.DEBUG)
        results = solver.solve(bped_m[4])
        assert_optimal_termination(results)
        assert value(bped_m[4].fs.unit.N_Re) == pytest.approx(21.235, rel=1e-3)

        assert value(bped_m[4].fs.unit.pressure_drop[0]) == pytest.approx(
            2263.540, rel=1e-3
        )

        assert value(bped_m[4].fs.unit.pressure_drop_total[0]) == pytest.approx(
            1788.196, rel=1e-3
        )

        # Test bped_m5
        bped_m[5].fs.unit.diffus_mass.fix(1.6e-9)
        bped_m[5].fs.unit.spacer_specific_area.fix(10700)
        iscale.calculate_scaling_factors(bped_m[5])
        assert degrees_of_freedom(bped_m[5]) == 0
        initialization_tester(bped_m[5], outlvl=idaeslog.DEBUG)
        iscale.calculate_scaling_factors(bped_m[5])
        results = solver.solve(bped_m[5])
        assert_optimal_termination(results)
        assert value(bped_m[5].fs.unit.N_Re) == pytest.approx(7.641, rel=1e-3)

        assert value(bped_m[5].fs.unit.pressure_drop[0]) == pytest.approx(
            10486.482, rel=1e-3
        )

        assert value(bped_m[5].fs.unit.pressure_drop_total[0]) == pytest.approx(
            8284.321, rel=1e-3
        )
