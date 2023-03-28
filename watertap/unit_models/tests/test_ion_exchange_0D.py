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

import pytest
from pyomo.environ import (
    ConcreteModel,
    value,
    Var,
    Param,
    Objective,
    units as pyunits,
    assert_optimal_termination,
)
from pyomo.network import Port
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
    MCASStateBlock,
)
from watertap.unit_models.ion_exchange_0D import (
    IonExchange0D,
    IonExchangeType,
    RegenerantChem,
    IsothermType,
    DiffusionControlType,
)
from watertap.costing import WaterTAPCosting
from watertap.core.util.initialization import check_dof

from idaes.core import (
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.solvers.get_solver import get_solver
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
    set_scaling_factor,
)
from pyomo.util.check_units import assert_units_consistent
import idaes.logger as idaeslog

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


def get_comp_lst(ix, comp=Var):
    cs = []
    for c in ix.component_objects(comp):
        if "ref" in c.name or "process_flow" in c.name or "regeneration_" in c.name:
            continue
        cs.append(c.name.split("fs.ix.")[1])
    return cs


def ix_scaling(m, sf=1e4, est_removal=0.99, est_recov=0.99, sf_mass=1e4):

    ix = m.fs.ix
    prop_in = ix.properties_in[0]
    prop_out = ix.properties_out[0]
    prop_regen = ix.properties_regen[0]
    target_ion = ix.config.target_ion
    ions = ix.config.property_package.ion_set
    sf_flow = 1 / prop_in.flow_mol_phase_comp["Liq", "H2O"].value
    set_scaling_factor(prop_in.flow_mol_phase_comp["Liq", "H2O"], sf_flow)
    set_scaling_factor(prop_out.flow_mol_phase_comp["Liq", "H2O"], sf_flow)
    set_scaling_factor(prop_regen.flow_mol_phase_comp["Liq", "H2O"], sf)
    for ion in ions:
        set_scaling_factor(
            prop_in.flow_mol_phase_comp["Liq", ion],
            1 / prop_in.flow_mol_phase_comp["Liq", ion].value,
        )
        if ion == target_ion:
            sf_targ = 1 / (
                prop_in.flow_mol_phase_comp["Liq", ion].value * (1 - est_removal)
            )
            set_scaling_factor(prop_out.flow_mol_phase_comp["Liq", ion], sf_targ)
            set_scaling_factor(prop_out.conc_equiv_phase_comp["Liq", ion], sf_targ)
            set_scaling_factor(prop_out.conc_mol_phase_comp["Liq", ion], sf_targ)
            set_scaling_factor(prop_out.conc_mass_phase_comp["Liq", ion], sf_targ)
            set_scaling_factor(prop_regen.flow_mol_phase_comp["Liq", ion], sf_targ)
        else:
            set_scaling_factor(
                prop_out.flow_mol_phase_comp["Liq", ion],
                1 / (prop_in.flow_mol_phase_comp["Liq", ion].value),
            )
            set_scaling_factor(
                prop_regen.flow_mol_phase_comp["Liq", ion],
                1 / (prop_in.flow_mol_phase_comp["Liq", ion].value),
            )

    calculate_scaling_factors(m)
    set_scaling_factor(ix.mass_in[target_ion], 1 / sf_mass)
    set_scaling_factor(ix.mass_removed[target_ion], 1 / sf_mass)
    set_scaling_factor(ix.mass_out[target_ion], sf_mass)
    return m


class TestIonExchangeDefault:
    @pytest.fixture(scope="class")
    def IX_default(self):
        target_ion = "Ca_2+"
        ion_props = {
            "solute_list": ["Ca_2+"],
            "diffusivity_data": {("Liq", "Ca_2+"): 9.2e-10},
            "mw_data": {"H2O": 0.018, "Ca_2+": 0.04},
            "charge": {"Ca_2+": 2},
        }
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = MCASParameterBlock(**ion_props)
        ix_config = {
            "property_package": m.fs.properties,
            "target_ion": target_ion,
        }
        ix = m.fs.ix = IonExchange0D(**ix_config)
        ix.process_flow.properties_in.calculate_state(
            var_args={
                ("flow_vol_phase", "Liq"): 0.5,
                ("conc_mass_phase_comp", ("Liq", target_ion)): 0.1,
                ("pressure", None): 101325,
                ("temperature", None): 298,
            },
            hold_state=True,
        )

        ix.service_flow_rate.fix(15)
        ix.langmuir[target_ion].fix(0.9)
        ix.resin_max_capacity.fix(3)
        ix.bed_depth.fix(1.7)
        ix.dimensionless_time.fix()
        ix.number_columns.fix(8)
        ix.resin_diam.fix()
        ix.resin_bulk_dens.fix()
        ix.bed_porosity.fix()
        ix.regen_dose.fix()

        return m

    @pytest.mark.unit
    def test_config(self, IX_default):
        m = IX_default

        assert len(m.fs.ix.config) == 12

        assert not m.fs.ix.config.dynamic
        assert not m.fs.ix.config.has_holdup
        assert m.fs.ix.config.property_package is m.fs.properties
        assert not m.fs.ix.config.hazardous_waste
        assert isinstance(m.fs.ix.regen_chem, RegenerantChem)
        assert m.fs.ix.regen_chem is RegenerantChem.NaCl
        assert isinstance(m.fs.ix.ion_exchange_type, IonExchangeType)
        assert m.fs.ix.ion_exchange_type is IonExchangeType.cation
        assert isinstance(m.fs.ix.config.isotherm, IsothermType)
        assert m.fs.ix.config.isotherm is IsothermType.langmuir
        assert isinstance(m.fs.ix.config.diffusion_control, DiffusionControlType)
        assert m.fs.ix.config.diffusion_control is DiffusionControlType.liquid
        assert isinstance(m.fs.ix.config.energy_balance_type, EnergyBalanceType)
        assert m.fs.ix.config.energy_balance_type is EnergyBalanceType.none
        assert isinstance(m.fs.ix.config.momentum_balance_type, MomentumBalanceType)
        assert m.fs.ix.config.momentum_balance_type is MomentumBalanceType.pressureTotal

    @pytest.mark.unit
    def test_default_build(self, IX_default):
        m = IX_default
        ix = m.fs.ix
        # test ports and variables
        port_lst = ["inlet", "outlet", "regen"]
        port_vars_lst = ["flow_mol_phase_comp", "pressure", "temperature"]
        for port_str in port_lst:
            assert hasattr(ix, port_str)
            port = getattr(ix, port_str)
            assert len(port.vars) == 3
            assert isinstance(port, Port)
            for var_str in port_vars_lst:
                assert hasattr(port, var_str)
                var = getattr(port, var_str)
                assert isinstance(var, Var)

        #         # test unit objects (including parameters, variables, and constraints)
        ix_params = [
            "diff_ion_resin",
            "underdrain_h",
            "distributor_h",
            "p_drop_psi_to_m",
            "holdup_A",
            "holdup_B",
            "holdup_exp",
            "Pe_p_A",
            "Pe_p_exp",
            "Sh_A",
            "Sh_exp_A",
            "Sh_exp_B",
            "bed_expansion_frac_A",
            "bed_expansion_frac_B",
            "bed_expansion_frac_C",
            "p_drop_A",
            "p_drop_B",
            "p_drop_C",
            "pump_efficiency",
            "t_regen",
            "rinse_bv",
            "bw_rate",
            "t_bw",
            "service_to_regen_flow_ratio",
            "number_columns_redund",
        ]

        for p in ix_params:
            assert hasattr(ix, p)
            param = getattr(ix, p)
            assert isinstance(param, Param)

        ix_vars = [
            "resin_max_capacity",
            "resin_eq_capacity",
            "resin_unused_capacity",
            "resin_diam",
            "resin_bulk_dens",
            "resin_particle_dens",
            "langmuir",
            "resin_surf_per_vol",
            "col_height_to_diam_ratio",
            "bed_vol_tot",
            "bed_depth",
            "bed_porosity",
            "col_height",
            "col_diam",
            "number_columns",
            "partition_ratio",
            "fluid_mass_transfer_coeff",
            "t_breakthru",
            "t_contact",
            "num_transfer_units",
            "dimensionless_time",
            "mass_in",
            "mass_removed",
            "mass_out",
            "vel_bed",
            "vel_inter",
            "holdup",
            "service_flow_rate",
            "Re",
            "Sc",
            "Sh",
            "Pe_p",
            "Pe_bed",
            "c_norm",
            "regen_dose",
        ]

        for v in ix_vars:
            assert hasattr(ix, v)
            var = getattr(ix, v)
            assert isinstance(var, Var)

        #         # test statistics
        assert number_variables(m) == 74
        assert number_total_constraints(m) == 46
        assert number_unused_variables(m) == 13

    @pytest.mark.unit
    def test_dof(self, IX_default):
        m = IX_default
        check_dof(m, fail_flag=True)

    @pytest.mark.unit
    def test_calculate_scaling(self, IX_default):
        m = IX_default
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e-4, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e-2, index=("Liq", "Ca_2+")
        )

        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

    @pytest.mark.component
    def test_initialize(self, IX_default):
        m = IX_default
        initialization_tester(m, unit=m.fs.ix, outlvl=idaeslog.DEBUG)

    @pytest.mark.component
    def test_solve(self, IX_default):
        m = IX_default
        results = solver.solve(m, tee=True)
        assert_units_consistent(m)
        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_conservation(self, IX_default):
        m = IX_default
        ix = m.fs.ix
        target_ion = ix.config.target_ion
        assert (
            abs(
                value(
                    ix.mass_in[target_ion]
                    - ix.mass_out[target_ion]
                    - ix.mass_removed[target_ion]
                )
            )
            <= 1e-6
        )

    @pytest.mark.component
    def test_solution(self, IX_default):
        m = IX_default
        ix = m.fs.ix
        target_ion = ix.config.target_ion

        results_dict = {
            "resin_max_capacity": 3,
            "resin_eq_capacity": 1.5547810762853225,
            "resin_unused_capacity": 1.4452189237146775,
            "resin_diam": 0.0007,
            "resin_bulk_dens": 0.7,
            "resin_particle_dens": 1.4,
            "langmuir": 0.9,
            "resin_surf_per_vol": 4285.714285714286,
            "col_height_to_diam_ratio": 1.0408526790314099,
            "bed_vol_tot": 120.0,
            "bed_depth": 1.7,
            "bed_porosity": 0.5,
            "col_height": 3.488715,
            "col_diam": 3.351785579537064,
            "number_columns": 8,
            "partition_ratio": 217.66935067994388,
            "fluid_mass_transfer_coeff": 0.00033082475202981344,
            "t_breakthru": 52360.64416318653,
            "t_contact": 120.0,
            "num_transfer_units": 340.2768878020938,
            "dimensionless_time": 1,
            "mass_in": 130901.61040796709,
            "mass_removed": 130601.61040796706,
            "mass_out": 300.00000000000443,
            "vel_bed": 0.007083333333333333,
            "vel_inter": 0.014166666666666666,
            "holdup": 111.54173661784432,
            "service_flow_rate": 15,
            "Re": 4.958333333333333,
            "Sc": 1086.9565217391305,
            "Sh": 251.7144852400754,
            "Pe_p": 0.10782790064157834,
            "Pe_bed": 261.86775870097597,
            "c_norm": 0.4919290557789296,
            "regen_dose": 300,
        }

        for v, val in results_dict.items():
            var = getattr(ix, v)
            if var.is_indexed():
                assert pytest.approx(val, rel=1e-3) == value(var[target_ion])
            else:
                assert pytest.approx(val, rel=1e-3) == value(var)

    @pytest.mark.component
    def test_costing(self, IX_default):
        m = IX_default
        ix = m.fs.ix

        m.fs.costing = WaterTAPCosting()
        ix.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(ix.process_flow.properties_out[0].flow_vol_phase["Liq"])
        m.fs.costing.add_specific_energy_consumption(
            ix.process_flow.properties_out[0].flow_vol_phase["Liq"]
        )

        results = solver.solve(m, tee=True)
        assert_optimal_termination(results)

        assert pytest.approx(3296998.040376, rel=1e-5) == value(
            m.fs.costing.aggregate_capital_cost
        )
        assert pytest.approx(2150466.47169, rel=1e-5) == value(
            m.fs.costing.total_operating_cost
        )
        assert pytest.approx(6593996.0807, rel=1e-5) == value(
            m.fs.costing.total_capital_cost
        )
        assert pytest.approx(0.197884, rel=1e-5) == value(m.fs.costing.LCOW)
        assert pytest.approx(0.0572452, rel=1e-5) == value(
            m.fs.costing.specific_energy_consumption
        )


# class TestIonExchangeWithInert:
#     @pytest.fixture(scope="class")
#     def IX_frame_with_inert(self):
#         target_ion = "Ca_2+"
#         inert_ions = {"Cl_-": 1e-4, "SO4_2-": 2e-3}
#         mass_frac = 1.5e-4
#         ions = [target_ion] + [x for x in inert_ions.keys()]
#         ix_in = get_ix_in(ions)
#         m = ConcreteModel()
#         m.fs = FlowsheetBlock(dynamic=False)
#         m.fs.properties = MCASParameterBlock(**ix_in)
#         ix_unit_in = {
#             "property_package": m.fs.properties,
#             "target_ion": target_ion,
#             "hazardous_waste": False,
#         }

#         ix = m.fs.ix = IonExchange0D(**ix_unit_in)
#         mass_flow_in = 5 * pyunits.kg / pyunits.s
#         ix_prop_in = inert_ions
#         ix_prop_in[target_ion] = mass_frac

#         # Fix mole flow rates of each ion and water
#         for ion, x in ix_prop_in.items():
#             mol_comp_flow = (
#                 x
#                 * pyunits.kg
#                 / pyunits.kg
#                 * mass_flow_in
#                 / m.fs.ix.config.property_package.mw_comp[ion]
#             )
#             ix.inlet.flow_mol_phase_comp[0, "Liq", ion].fix(mol_comp_flow)
#         H2O_mass_frac = 1 - sum(x for x in ix_prop_in.values())
#         H2O_mol_comp_flow = (
#             H2O_mass_frac
#             * pyunits.kg
#             / pyunits.kg
#             * mass_flow_in
#             / m.fs.ix.config.property_package.mw_comp["H2O"]
#         )

#         ix.inlet.pressure[0].fix(101325)
#         ix.inlet.temperature[0].fix(298.15)
#         ix.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix(H2O_mol_comp_flow)
#         ix.langmuir[target_ion].fix(0.8)
#         ix.resin_max_capacity.fix(1.5)
#         ix.service_flow_rate.fix(10)
#         ix.number_columns.fix(5)
#         ix.bed_depth.fix(1.7634969859421077)
#         ix.resin_diam.fix()
#         ix.resin_bulk_dens.fix()
#         ix.bed_porosity.fix()
#         ix.dimensionless_time.fix()
#         ix.regen_dose.fix()
#         ix.regen_recycle.fix()
#         ix.t_regen.fix()
#         ix.t_bw.fix()
#         ix.bw_rate.fix()
#         ix.rinse_bv.fix()
#         ix.pump_efficiency.fix()
#         ix.p_drop_A.fix()
#         ix.p_drop_B.fix()
#         ix.p_drop_C.fix()
#         ix.bed_expansion_frac_A.fix()
#         ix.bed_expansion_frac_B.fix()
#         ix.bed_expansion_frac_C.fix()
#         ix.service_to_regen_flow_ratio.fix()
#         ix.number_columns_redund.fix()

#         return m

#     @pytest.mark.unit
#     def test_config(self, IX_frame_with_inert):
#         m = IX_frame_with_inert

#         assert len(m.fs.ix.config) == 7

#         assert not m.fs.ix.config.dynamic
#         assert not m.fs.ix.config.has_holdup
#         assert m.fs.ix.config.property_package is m.fs.properties
#         assert not m.fs.ix.config.hazardous_waste
#         assert isinstance(m.fs.ix.regen_chem, RegenerantChem)
#         assert m.fs.ix.regen_chem is RegenerantChem.HCl
#         assert isinstance(m.fs.ix.ion_exchange_type, IonExchangeType)
#         assert m.fs.ix.ion_exchange_type is IonExchangeType.cation

#     @pytest.mark.unit
#     def test_build(self, IX_frame_with_inert):
#         m = IX_frame_with_inert
#         ix = m.fs.ix
#         # test ports and variables
#         port_lst = ["inlet", "outlet", "regen"]
#         port_vars_lst = ["flow_mol_phase_comp", "pressure", "temperature"]
#         for port_str in port_lst:
#             assert hasattr(ix, port_str)
#             port = getattr(ix, port_str)
#             assert len(port.vars) == 3
#             assert isinstance(port, Port)
#             for var_str in port_vars_lst:
#                 assert hasattr(port, var_str)
#                 var = getattr(port, var_str)
#                 assert isinstance(var, Var)

#         # test unit objects (including parameters, variables, and constraints)
#         ix_params = [
#             "diff_ion_resin",
#             "underdrain_h",
#             "distributor_h",
#             "holdup_A",
#             "holdup_B",
#             "holdup_exp",
#             "Pe_p_A",
#             "Pe_p_exp",
#             "Sh_A",
#             "Sh_exp",
#         ]

#         for p in ix_params:
#             assert hasattr(ix, p)
#             param = getattr(ix, p)
#             assert isinstance(param, Param)

#         ix_vars = [
#             "bed_depth_to_diam_ratio",
#             "bed_expansion_frac_A",
#             "bed_expansion_frac_B",
#             "bed_expansion_frac_C",
#             "p_drop_A",
#             "p_drop_B",
#             "p_drop_C",
#             "resin_max_capacity",
#             "resin_eq_capacity",
#             "resin_unused_capacity",
#             "resin_diam",
#             "resin_bulk_dens",
#             "resin_particle_dens",
#             "langmuir",
#             "separation_factor",
#             "resin_surf_per_vol",
#             "bed_vol",
#             "bed_vol_tot",
#             "bed_depth",
#             "bed_porosity",
#             "col_height",
#             "col_diam",
#             "col_vol_per",
#             "number_columns",
#             "number_columns_redund",
#             "partition_ratio",
#             "fluid_mass_transfer_coeff",
#             "rate_coeff",
#             "t_breakthru",
#             "t_cycle",
#             "t_contact",
#             "t_waste",
#             "num_transfer_units",
#             "HTU",
#             "dimensionless_time",
#             "lh",
#             "mass_in",
#             "mass_removed",
#             "mass_out",
#             "vel_bed",
#             "vel_inter",
#             "holdup",
#             "service_flow_rate",
#             "pressure_drop",
#             "Re",
#             "Sc",
#             "Sh",
#             "Pe_p",
#             "Pe_bed",
#             "c_norm",
#             "service_to_regen_flow_ratio",
#             "regen_recycle",
#             "regen_dose",
#             "t_regen",
#             "bw_rate",
#             "bw_flow",
#             "t_bw",
#             "bed_expansion_frac",
#             "bed_expansion_h",
#             "rinse_bv",
#             "rinse_flow",
#             "t_rinse",
#             "main_pump_power",
#             "regen_pump_power",
#             "bw_pump_power",
#             "rinse_pump_power",
#             "pump_efficiency",
#         ]

#         for v in ix_vars:
#             assert hasattr(ix, v)
#             var = getattr(ix, v)
#             assert isinstance(var, Var)

#         # # test state block objects
#         stateblock_lst = [
#             "properties_in",
#             "properties_out",
#             "properties_regen",
#         ]

#         for sb_str in stateblock_lst:
#             sb = getattr(ix, sb_str)
#             assert isinstance(sb, MCASStateBlock)

#         # test statistics
#         assert number_variables(m) == 164
#         assert number_total_constraints(m) == 116
#         assert number_unused_variables(m) == 19

#     @pytest.mark.unit
#     def test_dof(self, IX_frame_with_inert):
#         m = IX_frame_with_inert
#         check_dof(m, fail_flag=True)

#     @pytest.mark.unit
#     def test_calculate_scaling(self, IX_frame_with_inert):
#         m = IX_frame_with_inert
#         m = ix_scaling(m, sf=1, est_recov=0.95, est_removal=0.95, sf_mass=1)

#         # check that all variables have scaling factors
#         unscaled_var_list = list(unscaled_variables_generator(m.fs.ix))
#         assert len(unscaled_var_list) == 0

#     @pytest.mark.component
#     def test_initialize(self, IX_frame_with_inert):
#         m = IX_frame_with_inert
#         initialization_tester(m, outlvl=idaeslog.DEBUG)

#     @pytest.mark.component
#     def test_solve(self, IX_frame_with_inert):
#         m = IX_frame_with_inert
#         results = solver.solve(m, tee=True)

#         # Check for optimal solution
#         assert_optimal_termination(results)

#     @pytest.mark.component
#     def test_conservation(self, IX_frame_with_inert):
#         m = IX_frame_with_inert
#         ix = m.fs.ix
#         target_ion = ix.config.target_ion
#         inert_ions = ix.config.property_package.ion_set
#         assert (
#             abs(
#                 value(
#                     ix.mass_in[target_ion]
#                     - ix.mass_out[target_ion]
#                     - ix.mass_removed[target_ion]
#                 )
#             )
#             <= 1e-6
#         )
#         for ion in inert_ions:
#             if ion == target_ion:
#                 continue
#             assert (
#                 value(
#                     ix.properties_in[0].conc_equiv_phase_comp["Liq", ion]
#                     - ix.properties_out[0].conc_equiv_phase_comp["Liq", ion]
#                 )
#                 <= 1e-6
#             )

#     @pytest.mark.component
#     def test_solution(self, IX_frame_with_inert):
#         m = IX_frame_with_inert
#         ix = m.fs.ix
#         target_ion = ix.config.target_ion

#         results_dict = {
#             "bed_depth_to_diam_ratio": 7.026492764446821,
#             "bed_expansion_frac_A": -0.0123,
#             "bed_expansion_frac_B": 0.102,
#             "bed_expansion_frac_C": -0.00135,
#             "p_drop_A": 0.609,
#             "p_drop_B": 0.173,
#             "p_drop_C": 0.000828,
#             "resin_max_capacity": 1.5,
#             "resin_eq_capacity": 0.8080397292960807,
#             "resin_unused_capacity": 0.6919602707039192,
#             "resin_diam": 0.0007,
#             "resin_bulk_dens": 0.7,
#             "resin_particle_dens": 1.4,
#             "langmuir": 0.8,
#             "separation_factor": 1.2500000000000009,
#             "resin_surf_per_vol": 4285.714285714286,
#             "bed_vol": 0.3600000000001751,
#             "bed_vol_tot": 1.8000000000008756,
#             "bed_depth": 1.7637358582688671,
#             "bed_porosity": 0.5,
#             "col_height": 3.5820211097127084,
#             "col_diam": 0.5097879165308988,
#             "col_vol_per": 0.73113419619473,
#             "number_columns": 5,
#             "number_columns_redund": 1,
#             "partition_ratio": 75.41704142414085,
#             "fluid_mass_transfer_coeff": 4.322369555862131e-05,
#             "rate_coeff": 0.00026463487076706926,
#             "t_breakthru": 27330.134911686862,
#             "t_cycle": 30630.134911686488,
#             "t_contact": 179.99999999992596,
#             "t_waste": 3299.99999999963,
#             "num_transfer_units": 66.68798743329957,
#             "HTU": 0.02644757964143663,
#             "dimensionless_time": 1,
#             "lh": 0,
#             "mass_in": 1024.8800588731324,
#             "mass_removed": 1018.1300589135858,
#             "mass_out": 6.749999959546603,
#             "vel_bed": 0.0048992662729692125,
#             "vel_inter": 0.009798532545942183,
#             "holdup": 102.66180257369356,
#             "service_flow_rate": 10,
#             "pressure_drop": 6.910024848393857,
#             "Re": 3.4294863910767805,
#             "Sc": 1086.9565217396594,
#             "Sh": 32.88759444677707,
#             "Pe_p": 0.09033997470879263,
#             "Pe_bed": 227.62264675835323,
#             "c_norm": 0.48299134878783817,
#             "service_to_regen_flow_ratio": 3,
#             "regen_recycle": 1,
#             "regen_dose": 300,
#             "t_regen": 1800,
#             "bw_rate": 5,
#             "bw_flow": 0.0014174458069060663,
#             "t_bw": 600,
#             "bed_expansion_frac": 0.46395000000000003,
#             "bed_expansion_h": 0.8182852514438409,
#             "rinse_bv": 5,
#             "rinse_flow": 0.005000000000002573,
#             "t_rinse": 899.9999999996298,
#             "main_pump_power": 0.29784481412490327,
#             "regen_pump_power": 0.0992816047082775,
#             "bw_pump_power": 0.08443577678496997,
#             "rinse_pump_power": 0.2978448141249117,
#             "pump_efficiency": 0.8,
#         }

#         for v, val in results_dict.items():
#             var = getattr(ix, v)
#             if var.is_indexed():
#                 assert pytest.approx(val, rel=1e-3) == value(var[target_ion])
#             else:
#                 assert pytest.approx(val, rel=1e-3) == value(var)


# class TestIonExchangeCosting:
#     @pytest.fixture(scope="class")
#     def IX_frame_costing(self):
#         target_ion = "Ca_2+"
#         mass_frac = 1e-4
#         ions = [target_ion]
#         ix_in = get_ix_in(ions)
#         m = ConcreteModel()
#         m.fs = FlowsheetBlock(dynamic=False)
#         m.fs.properties = MCASParameterBlock(**ix_in)
#         ix_unit_in = {
#             "property_package": m.fs.properties,
#             "target_ion": target_ion,
#             "hazardous_waste": False,
#         }
#         ix = m.fs.ix = IonExchange0D(**ix_unit_in)
#         mass_flow_in = 50 * pyunits.kg / pyunits.s
#         ix_prop_in = {target_ion: mass_frac}
#         for ion, x in ix_prop_in.items():
#             mol_comp_flow = (
#                 x
#                 * pyunits.kg
#                 / pyunits.kg
#                 * mass_flow_in
#                 / m.fs.ix.config.property_package.mw_comp[ion]
#             )
#             ix.inlet.flow_mol_phase_comp[0, "Liq", ion].fix(mol_comp_flow)
#         H2O_mass_frac = 1 - sum(x for x in ix_prop_in.values())
#         H2O_mol_comp_flow = (
#             H2O_mass_frac
#             * pyunits.kg
#             / pyunits.kg
#             * mass_flow_in
#             / m.fs.ix.config.property_package.mw_comp["H2O"]
#         )

#         ix.inlet.pressure[0].fix(101325)
#         ix.inlet.temperature[0].fix(298.15)
#         ix.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix(H2O_mol_comp_flow)
#         ix.langmuir[target_ion].fix(0.7)
#         ix.resin_max_capacity.fix(3)
#         ix.service_flow_rate.fix(15)
#         ix.number_columns.fix(4)
#         ix.bed_depth.fix(1.6975735806794794)
#         ix.resin_diam.fix()
#         ix.resin_bulk_dens.fix()
#         ix.bed_porosity.fix()
#         ix.dimensionless_time.fix()
#         ix.regen_dose.fix()
#         ix.regen_recycle.fix()
#         ix.t_regen.fix()
#         ix.t_bw.fix()
#         ix.bw_rate.fix()
#         ix.rinse_bv.fix()
#         ix.pump_efficiency.fix()
#         ix.p_drop_A.fix()
#         ix.p_drop_B.fix()
#         ix.p_drop_C.fix()
#         ix.bed_expansion_frac_A.fix()
#         ix.bed_expansion_frac_B.fix()
#         ix.bed_expansion_frac_C.fix()
#         ix.service_to_regen_flow_ratio.fix()
#         ix.number_columns_redund.fix()

#         m.fs.costing = WaterTAPCosting()
#         ix.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
#         m.fs.costing.cost_process()
#         m.fs.costing.add_LCOW(ix.properties_out[0].flow_vol_phase["Liq"])

#         m = ix_scaling(m, sf=1, est_recov=0.99, est_removal=0.99, sf_mass=1e2)

#         def obj_rule(m):
#             return m.fs.costing.LCOW

#         m.obj = Objective(rule=obj_rule)

#         return m

#     @pytest.mark.component
#     def test_costing(self, IX_frame_costing):
#         m = IX_frame_costing

#         m.fs.ix.initialize(outlvl=idaeslog.DEBUG)
#         m.fs.costing.initialize()

#         assert degrees_of_freedom(m) == 0

#         results = solver.solve(m, tee=True)
#         assert_optimal_termination(results)

#         assert pytest.approx(418913.095, rel=1e-5) == value(
#             m.fs.costing.aggregate_capital_cost
#         )
#         assert pytest.approx(1025648.186, rel=1e-5) == value(
#             m.fs.costing.total_operating_cost
#         )
#         assert pytest.approx(837826.191, rel=1e-5) == value(
#             m.fs.costing.total_capital_cost
#         )
#         assert pytest.approx(0.7827819, rel=1e-5) == value(m.fs.costing.LCOW)
