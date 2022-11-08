###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# 'https://github.com/watertap-org/watertap/'
#
###############################################################################
import pytest
import pyomo.environ as pyo
from pyomo.environ import (
    ConcreteModel,
    check_optimal_termination,
    value,
)
from pyomo.network import Port
from pyomo.util.check_units import assert_units_consistent

from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    MomentumBalanceType,
)
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    badly_scaled_var_generator,
)
from idaes.core import UnitModelCostingBlock

from watertap.property_models.ion_DSPMDE_prop_pack import (
    DSPMDEParameterBlock,
)
from watertap.unit_models.gac import (
    GAC,
    FilmTransferCoefficientType,
    SurfaceDiffusionCoefficientType,
)
from watertap.costing import WaterTAPCosting

__author__ = "Hunter Barber"

solver = get_solver()


# -----------------------------------------------------------------------------
class TestGACSimplified:
    @pytest.fixture(scope="class")
    def gac_frame_simplified(self):
        ms = ConcreteModel()
        ms.fs = FlowsheetBlock(dynamic=False)

        ms.fs.properties = DSPMDEParameterBlock(
            solute_list=["DCE"], mw_data={"H2O": 0.018, "DCE": 0.09896}
        )

        ms.fs.unit = GAC(
            property_package=ms.fs.properties,
            film_transfer_coefficient_type="fixed",
            surface_diffusion_coefficient_type="fixed",
        )

        # feed specifications
        ms.fs.unit.process_flow.properties_in[0].pressure.fix(
            101325
        )  # feed pressure [Pa]
        ms.fs.unit.process_flow.properties_in[0].temperature.fix(
            273.15 + 25
        )  # feed temperature [K]
        ms.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp["Liq", "H2O"].fix(
            55555.55426666667
        )
        ms.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp["Liq", "DCE"].fix(
            0.0002344381568310428
        )

        # trial problem from Hand, 1984 for removal of trace DCE
        ms.fs.unit.conc_ratio_replace.fix(0.50)
        ms.fs.unit.freund_k.fix(37.9e-6 * (1e6**0.8316))
        ms.fs.unit.freund_ninv.fix(0.8316)
        ms.fs.unit.ebct.fix(300)  # seconds
        ms.fs.unit.bed_voidage.fix(0.449)
        ms.fs.unit.bed_length.fix(6)  # assumed
        ms.fs.unit.particle_porosity.fix(0.5)
        ms.fs.unit.particle_dens_app.fix(722)
        ms.fs.unit.particle_dia.fix(0.00106)
        ms.fs.unit.kf.fix(3.29e-5)
        ms.fs.unit.ds.fix(1.77e-13)
        ms.fs.unit.a0.fix(3.68421)
        ms.fs.unit.a1.fix(13.1579)
        ms.fs.unit.b0.fix(0.784576)
        ms.fs.unit.b1.fix(0.239663)
        ms.fs.unit.b2.fix(0.484422)
        ms.fs.unit.b3.fix(0.003206)
        ms.fs.unit.b4.fix(0.134987)

        return ms

    @pytest.mark.unit
    def test_config_simplified(self, gac_frame_simplified):
        ms = gac_frame_simplified
        # check unit config arguments
        assert len(ms.fs.unit.config) == 9

        assert not ms.fs.unit.config.dynamic
        assert not ms.fs.unit.config.has_holdup
        assert ms.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
        assert (
            ms.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
        )
        assert (
            ms.fs.unit.config.film_transfer_coefficient_type
            == FilmTransferCoefficientType.fixed
        )
        assert (
            ms.fs.unit.config.surface_diffusion_coefficient_type
            == SurfaceDiffusionCoefficientType.fixed
        )

        assert ms.fs.unit.config.property_package is ms.fs.properties
        assert len(ms.fs.unit.config.property_package.solute_set) == 1
        assert len(ms.fs.unit.config.property_package.solvent_set) == 1

    @pytest.mark.unit
    def test_build_simplified(self, gac_frame_simplified):
        ms = gac_frame_simplified

        # test units
        assert assert_units_consistent(ms) is None

        # test ports
        port_lst = ["inlet", "outlet", "adsorbed"]
        for port_str in port_lst:
            port = getattr(ms.fs.unit, port_str)
            assert len(port.vars) == 3  # number of state variables for property package
            assert isinstance(port, Port)

        # test statistics
        assert number_variables(ms) == 77
        assert number_total_constraints(ms) == 44
        assert number_unused_variables(ms) == 10  # dens parameters from properties

    @pytest.mark.unit
    def test_dof_simplified(self, gac_frame_simplified):
        ms = gac_frame_simplified
        assert degrees_of_freedom(ms) == 0

    @pytest.mark.unit
    def test_calculate_scaling_simplified(self, gac_frame_simplified):
        ms = gac_frame_simplified

        ms.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e-4, index=("Liq", "H2O")
        )
        ms.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e4, index=("Liq", "DCE")
        )
        calculate_scaling_factors(ms)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(ms))
        assert len(unscaled_var_list) == 0

    @pytest.mark.component
    def test_initialize_simplified(self, gac_frame_simplified):
        initialization_tester(gac_frame_simplified)

    @pytest.mark.component
    def test_var_scaling_init_simplified(self, gac_frame_simplified):
        ms = gac_frame_simplified
        badly_scaled_var_lst = list(
            badly_scaled_var_generator(ms, large=1e2, small=1e-2, zero=1e-8)
        )
        assert badly_scaled_var_lst == []

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_solve_simplified(self, gac_frame_simplified):
        ms = gac_frame_simplified
        results = solver.solve(ms)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.component
    def test_var_scaling_solve_simplified(self, gac_frame_simplified):
        ms = gac_frame_simplified
        badly_scaled_var_lst = list(
            badly_scaled_var_generator(ms, large=1e2, small=1e-2, zero=1e-8)
        )
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_solution_simplified(self, gac_frame_simplified):
        ms = gac_frame_simplified

        # Approx data pulled from graph in Hand, 1984 at ~30 days
        # 30 days adjusted to actual solution to account for web plot data extraction error within reason
        assert (
            pytest.approx(29.563, rel=1e-3) == value(ms.fs.unit.elap_time) / 24 / 3600
        )


# -----------------------------------------------------------------------------
class TestGACRobust:
    @pytest.fixture(scope="class")
    def gac_frame_robust(self):
        mr = ConcreteModel()
        mr.fs = FlowsheetBlock(dynamic=False)

        mr.fs.properties = DSPMDEParameterBlock(
            solute_list=["TCE"], mw_data={"H2O": 0.018, "TCE": 0.1314}
        )

        mr.fs.unit = GAC(
            property_package=mr.fs.properties,
            film_transfer_coefficient_type="calculated",
            surface_diffusion_coefficient_type="calculated",
        )

        # feed specifications
        mr.fs.unit.process_flow.properties_in[0].pressure.fix(
            101325
        )  # feed pressure [Pa]
        mr.fs.unit.process_flow.properties_in[0].temperature.fix(
            273.15 + 25
        )  # feed temperature [K]
        mr.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp["Liq", "H2O"].fix(
            824.0736620370348
        )
        mr.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp["Liq", "TCE"].fix(
            5.644342973110135e-05
        )
        # touch variables relevant to design, may be fixed as opposed to ftp in flowsheets
        mr.fs.unit.process_flow.properties_in[0].flow_vol_phase["Liq"]
        mr.fs.unit.process_flow.properties_in[0].conc_mass_phase_comp
        # test construction of other variables in .report()
        mr.fs.unit.process_flow.properties_out[0].flow_vol_phase["Liq"]
        mr.fs.unit.adsorbed_contam[0].flow_vol_phase["Liq"]

        # trial problem from Crittenden, 2012 for removal of TCE
        mr.fs.unit.conc_ratio_replace.fix(0.80)
        mr.fs.unit.freund_k.fix(1062e-6 * (1e6**0.48))
        mr.fs.unit.freund_ninv.fix(0.48)
        mr.fs.unit.ebct.fix(10 * 60)
        mr.fs.unit.bed_voidage.fix(0.44)
        mr.fs.unit.particle_porosity.fix(0.641)
        mr.fs.unit.particle_dens_app.fix(803.4)
        mr.fs.unit.particle_dia.fix(0.001026)
        mr.fs.unit.velocity_sup.fix(5 / 3600)
        mr.fs.unit.molal_volume.fix(9.81e-5)
        mr.fs.unit.tort.fix(1)
        mr.fs.unit.spdfr.fix(1)
        mr.fs.unit.sphericity.fix(1.5)
        mr.fs.unit.a0.fix(0.8)
        mr.fs.unit.a1.fix(0)
        mr.fs.unit.b0.fix(0.023)
        mr.fs.unit.b1.fix(0.793673)
        mr.fs.unit.b2.fix(0.039324)
        mr.fs.unit.b3.fix(0.009326)
        mr.fs.unit.b4.fix(0.08275)

        return mr

    @pytest.mark.unit
    def test_config_robust(self, gac_frame_robust):
        mr = gac_frame_robust
        # check unit config arguments
        assert len(mr.fs.unit.config) == 9

        assert not mr.fs.unit.config.dynamic
        assert not mr.fs.unit.config.has_holdup
        assert mr.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
        assert (
            mr.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
        )
        assert (
            mr.fs.unit.config.film_transfer_coefficient_type
            == FilmTransferCoefficientType.calculated
        )
        assert (
            mr.fs.unit.config.surface_diffusion_coefficient_type
            == SurfaceDiffusionCoefficientType.calculated
        )

        assert mr.fs.unit.config.property_package is mr.fs.properties
        assert len(mr.fs.unit.config.property_package.solute_set) == 1
        assert len(mr.fs.unit.config.property_package.solvent_set) == 1

    @pytest.mark.unit
    def test_build_robust(self, gac_frame_robust):
        mr = gac_frame_robust

        # test units
        assert assert_units_consistent(mr) is None

        # test ports
        port_lst = ["inlet", "outlet", "adsorbed"]
        for port_str in port_lst:
            port = getattr(mr.fs.unit, port_str)
            assert len(port.vars) == 3  # number of state variables for property package
            assert isinstance(port, Port)

        # test statistics
        assert number_variables(mr) == 88
        assert number_total_constraints(mr) == 53
        assert number_unused_variables(mr) == 10  # dens parameters from properties

    @pytest.mark.unit
    def test_dof_robust(self, gac_frame_robust):
        mr = gac_frame_robust
        assert degrees_of_freedom(mr) == 0

    @pytest.mark.unit
    def test_calculate_scaling_robust(self, gac_frame_robust):
        mr = gac_frame_robust

        mr.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e-2, index=("Liq", "H2O")
        )
        mr.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e5, index=("Liq", "TCE")
        )
        calculate_scaling_factors(mr)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(mr))
        assert len(unscaled_var_list) == 0

    @pytest.mark.component
    def test_initialize_robust(self, gac_frame_robust):
        initialization_tester(gac_frame_robust)

    @pytest.mark.component
    def test_var_scaling_init_robust(self, gac_frame_robust):
        mr = gac_frame_robust
        badly_scaled_var_lst = list(
            badly_scaled_var_generator(mr, large=1e2, small=1e-2, zero=1e-8)
        )
        assert badly_scaled_var_lst == []

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_solve_robust(self, gac_frame_robust):
        mr = gac_frame_robust
        results = solver.solve(mr)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.component
    def test_var_scaling_solve_robust(self, gac_frame_robust):
        mr = gac_frame_robust
        badly_scaled_var_lst = list(
            badly_scaled_var_generator(mr, large=1e2, small=1e-2, zero=1e-8)
        )
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_solution_robust(self, gac_frame_robust):
        mr = gac_frame_robust

        # values calculated independently and near to those reported in Crittenden, 2012
        assert pytest.approx(1.139, rel=1e-3) == value(mr.fs.unit.mass_throughput)
        assert pytest.approx(12830000, rel=1e-3) == value(mr.fs.unit.elap_time)
        assert pytest.approx(10.68, rel=1e-3) == value(mr.fs.unit.bed_area)

    @pytest.mark.component
    def test_reporting_robust(self, gac_frame_robust):
        mr = gac_frame_robust
        mr.fs.unit.report()

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_costing_robust(self, gac_frame_robust):
        mr = gac_frame_robust

        mr.fs.costing = WaterTAPCosting()
        mr.fs.costing.base_currency = pyo.units.USD_2020

        mr.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=mr.fs.costing
        )
        mr.fs.costing.cost_process()
        results = solver.solve(mr)

        # Check for optimal solution
        assert check_optimal_termination(results)

        # Check for known cost solution of default twin alternating contactors
        assert value(mr.fs.costing.gac.num_contactors_op) == 1
        assert value(mr.fs.costing.gac.num_contactors_redundant) == 1
        assert pytest.approx(56900.93523, rel=1e-5) == value(
            mr.fs.unit.costing.contactor_cost
        )
        assert pytest.approx(4.359114384, rel=1e-5) == value(
            mr.fs.unit.costing.adsorbent_unit_cost
        )
        assert pytest.approx(17454.52868, rel=1e-5) == value(
            mr.fs.unit.costing.adsorbent_cost
        )
        assert pytest.approx(81692.69369, rel=1e-5) == value(
            mr.fs.unit.costing.other_process_cost
        )
        assert pytest.approx(156048.1576, rel=1e-5) == value(
            mr.fs.unit.costing.capital_cost
        )
        assert pytest.approx(13535.92023, rel=1e-5) == value(
            mr.fs.unit.costing.gac_makeup_cost
        )
        assert pytest.approx(29524.89977, rel=1e-5) == value(
            mr.fs.unit.costing.gac_regen_cost
        )
        assert pytest.approx(43060.81999, rel=1e-5) == value(
            mr.fs.unit.costing.fixed_operating_cost
        )

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_costing_modular_contactors_robust(self, gac_frame_robust):
        mr = gac_frame_robust

        mr.fs.costing = WaterTAPCosting()
        mr.fs.costing.base_currency = pyo.units.USD_2020

        mr.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=mr.fs.costing
        )
        mr.fs.costing.cost_process()

        mr.fs.costing.gac.num_contactors_op.fix(4)
        mr.fs.costing.gac.num_contactors_redundant.fix(2)

        results = solver.solve(mr)

        # Check for known cost solution when changing volume scale of vessels in parallel
        assert value(mr.fs.costing.gac.num_contactors_op) == 4
        assert value(mr.fs.costing.gac.num_contactors_redundant) == 2
        assert pytest.approx(89035.16691, rel=1e-5) == value(
            mr.fs.unit.costing.contactor_cost
        )
        assert pytest.approx(69693.33132, rel=1e-5) == value(
            mr.fs.unit.costing.other_process_cost
        )
        assert pytest.approx(176183.0269, rel=1e-5) == value(
            mr.fs.unit.costing.capital_cost
        )

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_costing_max_gac_ref_robust(self, gac_frame_robust):
        mr = gac_frame_robust

        # scale flow up 10x
        mr.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp["Liq", "H2O"].fix(
            10 * 824.0736620370348
        )
        mr.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp["Liq", "TCE"].fix(
            10 * 5.644342973110135e-05
        )

        mr.fs.costing = WaterTAPCosting()
        mr.fs.costing.base_currency = pyo.units.USD_2020

        mr.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=mr.fs.costing
        )
        mr.fs.costing.cost_process()
        # not necessarily an optimum solution because poor scaling but just checking the conditional
        results = solver.solve(mr)

        # Check for bed_mass_gac_cost_ref to be overwritten if bed_mass_gac is greater than bed_mass_gac_cost_max_ref
        assert value(mr.fs.unit.bed_mass_gac) > value(
            mr.fs.costing.gac.bed_mass_max_ref
        )
        assert value(mr.fs.unit.costing.bed_mass_gac_ref) == (
            pytest.approx(value(mr.fs.costing.gac.bed_mass_max_ref), 1e-5)
        )


# -----------------------------------------------------------------------------
class TestGACMulti:
    @pytest.fixture(scope="class")
    def gac_frame_multi(self):
        mm = ConcreteModel()
        mm.fs = FlowsheetBlock(dynamic=False)

        # inserting arbitrary BackGround Solutes, Cations, and Anions to check handling
        mm.fs.properties = DSPMDEParameterBlock(
            solute_list=["TCE", "BGSOL", "BGCAT", "BGAN"],
            mw_data={
                "H2O": 0.018,
                "TCE": 0.1314,
                "BGSOL": 0.1,
                "BGCAT": 0.1,
                "BGAN": 0.1,
            },
            charge={"BGCAT": 1, "BGAN": -2},
        )

        # testing target_species arg
        mm.fs.unit = GAC(
            property_package=mm.fs.properties,
            film_transfer_coefficient_type="calculated",
            surface_diffusion_coefficient_type="calculated",
            target_species={"TCE"},
        )

        # feed specifications
        mm.fs.unit.process_flow.properties_in[0].pressure.fix(
            101325
        )  # feed pressure [Pa]
        mm.fs.unit.process_flow.properties_in[0].temperature.fix(
            273.15 + 25
        )  # feed temperature [K]
        mm.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp["Liq", "H2O"].fix(
            824.0736620370348
        )
        mm.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp["Liq", "TCE"].fix(
            5.644342973110135e-05
        )
        mm.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp[
            "Liq", "BGSOL"
        ].fix(5e-05)
        mm.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp[
            "Liq", "BGCAT"
        ].fix(2e-05)
        mm.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp["Liq", "BGAN"].fix(
            1e-05
        )
        # touch variables relevant to design, may be fixed as opposed to ftp in flowsheets
        mm.fs.unit.process_flow.properties_in[0].flow_vol_phase["Liq"]
        mm.fs.unit.process_flow.properties_in[0].conc_mass_phase_comp
        # test construction of other variables in .report()
        mm.fs.unit.process_flow.properties_out[0].flow_vol_phase["Liq"]
        mm.fs.unit.adsorbed_contam[0].flow_vol_phase["Liq"]

        # trial problem from Crittenden, 2012 for removal of TCE
        mm.fs.unit.conc_ratio_replace.fix(0.80)
        mm.fs.unit.freund_k.fix(1062e-6 * (1e6**0.48))
        mm.fs.unit.freund_ninv.fix(0.48)
        mm.fs.unit.ebct.fix(10 * 60)
        mm.fs.unit.bed_voidage.fix(0.44)
        mm.fs.unit.particle_porosity.fix(0.641)
        mm.fs.unit.particle_dens_app.fix(803.4)
        mm.fs.unit.particle_dia.fix(0.001026)
        mm.fs.unit.velocity_sup.fix(5 / 3600)
        mm.fs.unit.molal_volume.fix(9.81e-5)
        mm.fs.unit.tort.fix(1)
        mm.fs.unit.spdfr.fix(1)
        mm.fs.unit.sphericity.fix(1.5)
        mm.fs.unit.a0.fix(0.8)
        mm.fs.unit.a1.fix(0)
        mm.fs.unit.b0.fix(0.023)
        mm.fs.unit.b1.fix(0.793673)
        mm.fs.unit.b2.fix(0.039324)
        mm.fs.unit.b3.fix(0.009326)
        mm.fs.unit.b4.fix(0.08275)

        return mm

    @pytest.mark.unit
    def test_config_multi(self, gac_frame_multi):
        mm = gac_frame_multi

        # checking non-unity solute set and nonzero ion set handling
        assert len(mm.fs.unit.config.property_package.solute_set) == 2
        assert len(mm.fs.unit.config.property_package.solvent_set) == 1
        assert len(mm.fs.unit.config.property_package.ion_set) == 2

        assert degrees_of_freedom(mm) == 0

    @pytest.mark.unit
    def test_calculate_scaling_multi(self, gac_frame_multi):
        mm = gac_frame_multi

        mm.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e-2, index=("Liq", "H2O")
        )
        for j in mm.fs.properties.ion_set | mm.fs.properties.solute_set:
            mm.fs.properties.set_default_scaling(
                "flow_mol_phase_comp", 1e5, index=("Liq", j)
            )

        calculate_scaling_factors(mm)
        initialization_tester(gac_frame_multi)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(mm))
        assert len(unscaled_var_list) == 0

    @pytest.mark.unit
    def test_var_scaling_init_multi(self, gac_frame_multi):
        mm = gac_frame_multi
        badly_scaled_var_lst = list(
            badly_scaled_var_generator(mm, large=1e2, small=1e-2, zero=1e-8)
        )
        assert badly_scaled_var_lst == []

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_solve_multi(self, gac_frame_multi):
        mm = gac_frame_multi
        results = solver.solve(mm)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.unit
    def test_var_scaling_solve_multi(self, gac_frame_multi):
        mm = gac_frame_multi
        badly_scaled_var_lst = list(
            badly_scaled_var_generator(mm, large=1e2, small=1e-2, zero=1e-8)
        )
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_solution_multi(self, gac_frame_multi):
        mm = gac_frame_multi

        # values calculated independently and near to those reported in Crittenden, 2012
        assert pytest.approx(1.139, rel=1e-3) == value(mm.fs.unit.mass_throughput)
        assert pytest.approx(12830000, rel=1e-3) == value(mm.fs.unit.elap_time)
        assert pytest.approx(10.68, rel=1e-3) == value(mm.fs.unit.bed_area)

    @pytest.mark.component
    def test_reporting_multi(self, gac_frame_multi):
        mm = gac_frame_multi
        mm.fs.unit.report()
