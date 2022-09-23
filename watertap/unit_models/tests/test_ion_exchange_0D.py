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
from pyomo.environ import (
    ConcreteModel,
    Constraint,
    TerminationCondition,
    SolverStatus,
    value,
    Var,
    Param,
    Constraint,
    units as pyunits,
    assert_optimal_termination,
)
from pyomo.network import Port
from idaes.core import (
    FlowsheetBlock,
)
from watertap.property_models.ion_DSPMDE_prop_pack import (
    DSPMDEParameterBlock,
    DSPMDEStateBlock,
)
from watertap.unit_models.ion_exchange_0D import (
    IonExchange0D,
    IonExchangeType,
    RegenerantChem,
)

from watertap.core.util.initialization import check_dof

from idaes.core.solvers.get_solver import get_solver
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
    unused_variables_set,
)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    get_scaling_factor,
    unscaled_variables_generator,
    unscaled_constraints_generator,
    constraints_with_scale_factor_generator,
    badly_scaled_var_generator,
    set_scaling_factor,
)
import idaes.logger as idaeslog

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


def ix_scaling(m, sf=1e4, est_removal=0.99, est_recov=0.95, sf_mass=1e4, sf_inert=1e6):

    ix = m.fs.unit
    prop_in = ix.properties_in[0]
    prop_out = ix.properties_out[0]
    prop_regen = ix.properties_regen[0]
    target_ion = ix.config.target_ion
    ions = ix.config.property_package.ion_set
    sf_flow_in = 1 / prop_in.flow_mol_phase_comp["Liq", "H2O"].value
    sf_flow_out = 1 / (prop_in.flow_mol_phase_comp["Liq", "H2O"].value)
    sf_flow_regen = 1 / (
        prop_in.flow_mol_phase_comp["Liq", "H2O"].value * (1 - est_recov)
    )
    set_scaling_factor(prop_in.flow_mol_phase_comp["Liq", "H2O"], sf_flow_in)
    set_scaling_factor(prop_out.flow_mol_phase_comp["Liq", "H2O"], sf_flow_out)
    set_scaling_factor(prop_regen.flow_mol_phase_comp["Liq", "H2O"], sf_flow_regen)
    for ion in ions:
        set_scaling_factor(
            prop_in.flow_mol_phase_comp["Liq", ion],
            1 / prop_in.flow_mol_phase_comp["Liq", ion].value,
        )
        if ion == target_ion:
            sf_targ = 1 / (
                prop_in.flow_mol_phase_comp["Liq", ion].value * (1 - est_removal)
            )
            sf_targ_regen = 1 / (
                prop_in.flow_mol_phase_comp["Liq", ion].value * (est_removal)
            )
            set_scaling_factor(prop_out.flow_mol_phase_comp["Liq", ion], sf_targ)
            set_scaling_factor(prop_out.conc_equiv_phase_comp["Liq", ion], sf_targ)
            set_scaling_factor(prop_out.conc_mol_phase_comp["Liq", ion], sf_targ)
            set_scaling_factor(prop_out.conc_mass_phase_comp["Liq", ion], sf_targ)
            set_scaling_factor(
                prop_regen.flow_mol_phase_comp["Liq", ion], sf_targ_regen
            )
        else:
            set_scaling_factor(
                prop_out.flow_mol_phase_comp["Liq", ion],
                1 / (prop_in.flow_mol_phase_comp["Liq", ion].value),
            )
            set_scaling_factor(prop_regen.flow_mol_phase_comp["Liq", ion], sf_inert)
    calculate_scaling_factors(m)
    set_scaling_factor(ix.mass_in[target_ion], 1 / sf_mass)
    set_scaling_factor(ix.mass_removed[target_ion], 1 / sf_mass)
    set_scaling_factor(ix.mass_out[target_ion], sf_mass)
    return m


def get_ix_in(ions):
    diff_data = {
        "Na_+": 1.33e-9,
        "Ca_2+": 9.2e-10,
        "Cl_-": 2.03e-9,
        "Mg_2+": 0.706e-9,
        "SO4_2-": 1.06e-9,
        "PFAS_-": 0.49e-9,
        "Hardness_2+": 0.706e-9,
    }
    mw_data = {
        "Na_+": 23e-3,
        "Ca_2+": 40e-3,
        "Cl_-": 35e-3,
        "Mg_2+": 24e-3,
        "SO4_2-": 96e-3,
        "PFAS_-": 414.1e-3,
        "Hardness_2+": 100.0869e-3,
    }
    charge_data = {
        "Na_+": 1,
        "Ca_2+": 2,
        "Cl_-": -1,
        "Mg_2+": 2,
        "SO4_2-": -2,
        "PFAS_-": -1,
        "Hardness_2+": 2,
    }
    ix_in = {
        "solute_list": [],
        "diffusivity_data": {},
        "mw_data": {"H2O": 18e-3},
        "charge": {},
    }
    for ion in ions:
        ix_in["solute_list"].append(ion)
        ix_in["diffusivity_data"][("Liq", ion)] = diff_data[ion]
        ix_in["mw_data"][ion] = mw_data[ion]
        ix_in["charge"][ion] = charge_data[ion]
    return ix_in


class TestIonExchangeNoInert:
    @pytest.fixture(scope="class")
    def IX_frame_no_inert(self):
        target_ion = "Ca_2+"
        mass_frac = 1.5e-4
        ions = [target_ion]
        ix_in = get_ix_in(ions)
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.properties = DSPMDEParameterBlock(default=ix_in)
        ix_unit_in = {
            "property_package": m.fs.properties,
            "target_ion": target_ion,
            "hazardous_waste": False,
        }
        ix = m.fs.unit = IonExchange0D(default=ix_unit_in)
        mass_flow_in = 5 * pyunits.kg / pyunits.s
        ix_prop_in = {target_ion: mass_frac}
        for ion, x in ix_prop_in.items():
            mol_comp_flow = (
                x
                * pyunits.kg
                / pyunits.kg
                * mass_flow_in
                / m.fs.unit.config.property_package.mw_comp[ion]
            )
            ix.inlet.flow_mol_phase_comp[0, "Liq", ion].fix(mol_comp_flow)
        H2O_mass_frac = 1 - sum(x for x in ix_prop_in.values())
        H2O_mol_comp_flow = (
            H2O_mass_frac
            * pyunits.kg
            / pyunits.kg
            * mass_flow_in
            / m.fs.unit.config.property_package.mw_comp["H2O"]
        )

        ix.inlet.pressure[0].fix(101325)
        ix.inlet.temperature[0].fix(298.15)
        ix.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix(H2O_mol_comp_flow)
        ix.langmuir[target_ion].fix(0.8)
        ix.resin_max_capacity.fix(1.5)
        ix.service_flow_rate.fix(10)
        ix.number_columns.fix(5)
        ix.resin_diam.fix()
        ix.resin_bulk_dens.fix()
        ix.bed_porosity.fix()
        ix.dimensionless_time.fix()
        ix.lh.fix()
        ix.regen_dose.fix()
        ix.regen_recycle.fix()
        ix.t_regen.fix()
        ix.t_bw.fix()
        ix.bw_rate.fix()
        ix.rinse_bv.fix()
        ix.pump_efficiency.fix()
        ix.p_drop_A.fix()
        ix.p_drop_B.fix()
        ix.p_drop_C.fix()
        ix.bed_expansion_frac_A.fix()
        ix.bed_expansion_frac_B.fix()
        ix.bed_expansion_frac_C.fix()
        ix.service_to_regen_flow_ratio.fix()
        ix.number_columns_redund.fix()

        m = ix_scaling(m, sf=1, est_recov=0.99, est_removal=0.99, sf_mass=1e2)

        return m

    @pytest.mark.unit
    def test_config(self, IX_frame_no_inert):
        m = IX_frame_no_inert

        assert len(m.fs.unit.config) == 7

        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert m.fs.unit.config.property_package is m.fs.properties
        assert not m.fs.unit.config.hazardous_waste
        assert isinstance(m.fs.unit.regen_chem, RegenerantChem)
        assert m.fs.unit.regen_chem is RegenerantChem.HCl
        assert isinstance(m.fs.unit.ion_exchange_type, IonExchangeType)
        assert m.fs.unit.ion_exchange_type is IonExchangeType.cation

    @pytest.mark.unit
    def test_build(self, IX_frame_no_inert):
        m = IX_frame_no_inert
        ix = m.fs.unit
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

        # test unit objects (including parameters, variables, and constraints)
        ix_params = [
            "diff_ion_resin",
            "underdrain_h",
            "distributor_h",
            "holdup_A",
            "holdup_B",
            "holdup_exp",
            "Pe_p_A",
            "Pe_p_exp",
            "Sh_A",
            "Sh_exp",
        ]

        for p in ix_params:
            assert hasattr(ix, p)
            param = getattr(ix, p)
            assert isinstance(param, Param)

        ix_vars = [
            "bed_depth_to_diam_ratio",
            "bed_expansion_frac_A",
            "bed_expansion_frac_B",
            "bed_expansion_frac_C",
            "p_drop_A",
            "p_drop_B",
            "p_drop_C",
            "resin_max_capacity",
            "resin_eq_capacity",
            "resin_unused_capacity",
            "resin_diam",
            "resin_bulk_dens",
            "resin_particle_dens",
            "langmuir",
            "separation_factor",
            "resin_surf_per_vol",
            "bed_vol",
            "bed_vol_tot",
            "bed_depth",
            "bed_porosity",
            "col_height",
            "col_diam",
            "col_vol_per",
            "number_columns",
            "number_columns_redund",
            "partition_ratio",
            "fluid_mass_transfer_coeff",
            "rate_coeff",
            "t_breakthru",
            "t_cycle",
            "t_contact",
            "t_waste",
            "num_transfer_units",
            "HTU",
            "dimensionless_time",
            "lh",
            "mass_in",
            "mass_removed",
            "mass_out",
            "vel_bed",
            "vel_inter",
            "holdup",
            "service_flow_rate",
            "pressure_drop",
            "Re",
            "Sc",
            "Sh",
            "Pe_p",
            "Pe_bed",
            "c_norm",
            "service_to_regen_flow_ratio",
            "regen_recycle",
            "regen_dose",
            "t_regen",
            "bw_rate",
            "bw_flow",
            "t_bw",
            "bed_expansion_frac",
            "bed_expansion_h",
            "rinse_bv",
            "rinse_flow",
            "t_rinse",
            "main_pump_power",
            "regen_pump_power",
            "bw_pump_power",
            "rinse_pump_power",
            "pump_efficiency",
        ]

        for v in ix_vars:
            assert hasattr(ix, v)
            var = getattr(ix, v)
            assert isinstance(var, Var)

        # # test state block objects
        stateblock_lst = [
            "properties_in",
            "properties_out",
            "properties_regen",
        ]

        for sb_str in stateblock_lst:
            sb = getattr(ix, sb_str)
            assert isinstance(sb, DSPMDEStateBlock)

        # test statistics
        assert number_variables(m) == 123
        assert number_total_constraints(m) == 82
        assert number_unused_variables(m) == 15

    @pytest.mark.unit
    def test_dof(self, IX_frame_no_inert):
        m = IX_frame_no_inert
        check_dof(m, fail_flag=True)

    @pytest.mark.unit
    def test_calculate_scaling(self, IX_frame_no_inert):
        m = IX_frame_no_inert
        m = ix_scaling(m, sf=1, est_recov=0.99, est_removal=0.99, sf_mass=1e2)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m.fs.unit))
        assert len(unscaled_var_list) == 0

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_initialize(self, IX_frame_no_inert):
        m = IX_frame_no_inert
        initialization_tester(m)

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_solve(self, IX_frame_no_inert):
        m = IX_frame_no_inert
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_conservation(self, IX_frame_no_inert):
        m = IX_frame_no_inert
        ix = m.fs.unit
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

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_solution(self, IX_frame_no_inert):
        m = IX_frame_no_inert
        ix = m.fs.unit
        target_ion = ix.config.target_ion

        results_dict = {
            "bed_depth_to_diam_ratio": 7.0246377079257325,
            "bed_expansion_frac_A": -0.0123,
            "bed_expansion_frac_B": 0.102,
            "bed_expansion_frac_C": -0.00135,
            "p_drop_A": 0.609,
            "p_drop_B": 0.173,
            "p_drop_C": 0.000828,
            "resin_max_capacity": 1.5,
            "resin_eq_capacity": 0.8080397292955411,
            "resin_unused_capacity": 0.691960270704459,
            "resin_diam": 0.0007,
            "resin_bulk_dens": 0.7,
            "resin_particle_dens": 1.4,
            "langmuir": 0.8,
            "separation_factor": 1.25,
            "resin_surf_per_vol": 4285.714285714286,
            "bed_vol": 0.36,
            "bed_vol_tot": 1.8,
            "bed_depth": 1.7633544270675023,
            "bed_porosity": 0.5,
            "col_height": 3.58146271350547,
            "col_diam": 0.509843049907697,
            "col_vol_per": 0.731178347966125,
            "number_columns": 5,
            "number_columns_redund": 1,
            "partition_ratio": 75.4170414009172,
            "fluid_mass_transfer_coeff": 4.32206105893926e-05,
            "rate_coeff": 0.0002646159832003629,
            "t_breakthru": 27330.134904330193,
            "t_cycle": 30630.134904330196,
            "t_contact": 180.0,
            "t_waste": 3300.0,
            "num_transfer_units": 66.68322776649146,
            "HTU": 0.02644374734292022,
            "dimensionless_time": 1,
            "lh": 0,
            "mass_in": 1024.8800589123816,
            "mass_removed": 1018.1300589123816,
            "mass_out": 6.749999999999946,
            "vel_bed": 0.004898206741854173,
            "vel_inter": 0.009796413483708346,
            "holdup": 102.65685725360079,
            "service_flow_rate": 10,
            "pressure_drop": 6.907170449869965,
            "Re": 3.428744719297921,
            "Sc": 1086.9565217391305,
            "Sh": 32.88524718758132,
            "Pe_p": 0.0903305963042154,
            "Pe_bed": 227.54979556097942,
            "c_norm": 0.48299134878783684,
            "service_to_regen_flow_ratio": 3,
            "regen_recycle": 1,
            "regen_dose": 300,
            "t_regen": 1800,
            "bw_rate": 5,
            "bw_flow": 0.0014177524164314234,
            "t_bw": 600,
            "bed_expansion_frac": 0.46395000000000003,
            "bed_expansion_h": 0.8181082864379677,
            "rinse_bv": 5,
            "rinse_flow": 0.005000000000000002,
            "t_rinse": 900.0,
            "main_pump_power": 0.2977217801537614,
            "regen_pump_power": 0.09924059338462904,
            "bw_pump_power": 0.08441915464745203,
            "rinse_pump_power": 0.29772178015376144,
            "pump_efficiency": 0.8,
        }

        for v, val in results_dict.items():
            var = getattr(ix, v)
            if var.is_indexed():
                assert pytest.approx(val, rel=1e-3) == value(var[target_ion])
            else:
                assert pytest.approx(val, rel=1e-3) == value(var)


class TestIonExchangeWithInert:
    @pytest.fixture(scope="class")
    def IX_frame_with_inert(self):
        target_ion = "Ca_2+"
        inert_ions = {"Cl_-": 1e-4, "SO4_2-": 2e-3}
        mass_frac = 1.5e-4
        ions = [target_ion] + [x for x in inert_ions.keys()]
        ix_in = get_ix_in(ions)
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.properties = DSPMDEParameterBlock(default=ix_in)
        ix_unit_in = {
            "property_package": m.fs.properties,
            "target_ion": target_ion,
            "hazardous_waste": False,
        }

        ix = m.fs.unit = IonExchange0D(default=ix_unit_in)
        mass_flow_in = 5 * pyunits.kg / pyunits.s
        ix_prop_in = inert_ions
        ix_prop_in[target_ion] = mass_frac

        # Fix mole flow rates of each ion and water
        for ion, x in ix_prop_in.items():
            mol_comp_flow = (
                x
                * pyunits.kg
                / pyunits.kg
                * mass_flow_in
                / m.fs.unit.config.property_package.mw_comp[ion]
            )
            ix.inlet.flow_mol_phase_comp[0, "Liq", ion].fix(mol_comp_flow)
        H2O_mass_frac = 1 - sum(x for x in ix_prop_in.values())
        H2O_mol_comp_flow = (
            H2O_mass_frac
            * pyunits.kg
            / pyunits.kg
            * mass_flow_in
            / m.fs.unit.config.property_package.mw_comp["H2O"]
        )

        ix.inlet.pressure[0].fix(101325)
        ix.inlet.temperature[0].fix(298.15)
        ix.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix(H2O_mol_comp_flow)
        ix.langmuir[target_ion].fix(0.8)
        ix.resin_max_capacity.fix(1.5)
        ix.service_flow_rate.fix(10)
        ix.number_columns.fix(5)
        ix.resin_diam.fix()
        ix.resin_bulk_dens.fix()
        ix.bed_porosity.fix()
        ix.dimensionless_time.fix()
        ix.lh.fix()
        ix.regen_dose.fix()
        ix.regen_recycle.fix()
        ix.t_regen.fix()
        ix.t_bw.fix()
        ix.bw_rate.fix()
        ix.rinse_bv.fix()
        ix.pump_efficiency.fix()
        ix.p_drop_A.fix()
        ix.p_drop_B.fix()
        ix.p_drop_C.fix()
        ix.bed_expansion_frac_A.fix()
        ix.bed_expansion_frac_B.fix()
        ix.bed_expansion_frac_C.fix()
        ix.service_to_regen_flow_ratio.fix()
        ix.number_columns_redund.fix()

        m = ix_scaling(m, sf=1, est_recov=0.95, est_removal=0.95, sf_mass=1)

        return m

    @pytest.mark.unit
    def test_config(self, IX_frame_with_inert):
        m = IX_frame_with_inert

        assert len(m.fs.unit.config) == 7

        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert m.fs.unit.config.property_package is m.fs.properties
        assert not m.fs.unit.config.hazardous_waste
        assert isinstance(m.fs.unit.regen_chem, RegenerantChem)
        assert m.fs.unit.regen_chem is RegenerantChem.HCl
        assert isinstance(m.fs.unit.ion_exchange_type, IonExchangeType)
        assert m.fs.unit.ion_exchange_type is IonExchangeType.cation

    @pytest.mark.unit
    def test_build(self, IX_frame_with_inert):
        m = IX_frame_with_inert
        ix = m.fs.unit
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

        # test unit objects (including parameters, variables, and constraints)
        ix_params = [
            "diff_ion_resin",
            "underdrain_h",
            "distributor_h",
            "holdup_A",
            "holdup_B",
            "holdup_exp",
            "Pe_p_A",
            "Pe_p_exp",
            "Sh_A",
            "Sh_exp",
        ]

        for p in ix_params:
            assert hasattr(ix, p)
            param = getattr(ix, p)
            assert isinstance(param, Param)

        ix_vars = [
            "bed_depth_to_diam_ratio",
            "bed_expansion_frac_A",
            "bed_expansion_frac_B",
            "bed_expansion_frac_C",
            "p_drop_A",
            "p_drop_B",
            "p_drop_C",
            "resin_max_capacity",
            "resin_eq_capacity",
            "resin_unused_capacity",
            "resin_diam",
            "resin_bulk_dens",
            "resin_particle_dens",
            "langmuir",
            "separation_factor",
            "resin_surf_per_vol",
            "bed_vol",
            "bed_vol_tot",
            "bed_depth",
            "bed_porosity",
            "col_height",
            "col_diam",
            "col_vol_per",
            "number_columns",
            "number_columns_redund",
            "partition_ratio",
            "fluid_mass_transfer_coeff",
            "rate_coeff",
            "t_breakthru",
            "t_cycle",
            "t_contact",
            "t_waste",
            "num_transfer_units",
            "HTU",
            "dimensionless_time",
            "lh",
            "mass_in",
            "mass_removed",
            "mass_out",
            "vel_bed",
            "vel_inter",
            "holdup",
            "service_flow_rate",
            "pressure_drop",
            "Re",
            "Sc",
            "Sh",
            "Pe_p",
            "Pe_bed",
            "c_norm",
            "service_to_regen_flow_ratio",
            "regen_recycle",
            "regen_dose",
            "t_regen",
            "bw_rate",
            "bw_flow",
            "t_bw",
            "bed_expansion_frac",
            "bed_expansion_h",
            "rinse_bv",
            "rinse_flow",
            "t_rinse",
            "main_pump_power",
            "regen_pump_power",
            "bw_pump_power",
            "rinse_pump_power",
            "pump_efficiency",
        ]

        for v in ix_vars:
            assert hasattr(ix, v)
            var = getattr(ix, v)
            assert isinstance(var, Var)

        # # test state block objects
        stateblock_lst = [
            "properties_in",
            "properties_out",
            "properties_regen",
        ]

        for sb_str in stateblock_lst:
            sb = getattr(ix, sb_str)
            assert isinstance(sb, DSPMDEStateBlock)

        # test statistics
        assert number_variables(m) == 161
        assert number_total_constraints(m) == 116
        assert number_unused_variables(m) == 17

    @pytest.mark.unit
    def test_dof(self, IX_frame_with_inert):
        m = IX_frame_with_inert
        check_dof(m, fail_flag=True)

    @pytest.mark.unit
    def test_calculate_scaling(self, IX_frame_with_inert):
        m = IX_frame_with_inert
        m = ix_scaling(m, sf=1, est_recov=0.95, est_removal=0.95, sf_mass=1)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m.fs.unit))
        assert len(unscaled_var_list) == 0

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_initialize(self, IX_frame_with_inert):
        m = IX_frame_with_inert
        initialization_tester(m)

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_solve(self, IX_frame_with_inert):
        m = IX_frame_with_inert
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_conservation(self, IX_frame_with_inert):
        m = IX_frame_with_inert
        ix = m.fs.unit
        target_ion = ix.config.target_ion
        inert_ions = ix.config.property_package.ion_set
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
        for ion in inert_ions:
            if ion == target_ion:
                continue
            assert (
                value(
                    ix.properties_in[0].conc_equiv_phase_comp["Liq", ion]
                    - ix.properties_out[0].conc_equiv_phase_comp["Liq", ion]
                )
                <= 1e-6
            )

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_solution(self, IX_frame_with_inert):
        m = IX_frame_with_inert
        ix = m.fs.unit
        target_ion = ix.config.target_ion

        results_dict = {
            "bed_depth_to_diam_ratio": 7.026492764446821,
            "bed_expansion_frac_A": -0.0123,
            "bed_expansion_frac_B": 0.102,
            "bed_expansion_frac_C": -0.00135,
            "p_drop_A": 0.609,
            "p_drop_B": 0.173,
            "p_drop_C": 0.000828,
            "resin_max_capacity": 1.5,
            "resin_eq_capacity": 0.8080397292960807,
            "resin_unused_capacity": 0.6919602707039192,
            "resin_diam": 0.0007,
            "resin_bulk_dens": 0.7,
            "resin_particle_dens": 1.4,
            "langmuir": 0.8,
            "separation_factor": 1.2500000000000009,
            "resin_surf_per_vol": 4285.714285714286,
            "bed_vol": 0.3600000000001751,
            "bed_vol_tot": 1.8000000000008756,
            "bed_depth": 1.7637358582688671,
            "bed_porosity": 0.5,
            "col_height": 3.5820211097127084,
            "col_diam": 0.5097879165308988,
            "col_vol_per": 0.73113419619473,
            "number_columns": 5,
            "number_columns_redund": 1,
            "partition_ratio": 75.41704142414085,
            "fluid_mass_transfer_coeff": 4.322369555862131e-05,
            "rate_coeff": 0.00026463487076706926,
            "t_breakthru": 27330.134911686862,
            "t_cycle": 30630.134911686488,
            "t_contact": 179.99999999992596,
            "t_waste": 3299.99999999963,
            "num_transfer_units": 66.68798743329957,
            "HTU": 0.02644757964143663,
            "dimensionless_time": 1,
            "lh": 0,
            "mass_in": 1024.8800588731324,
            "mass_removed": 1018.1300589135858,
            "mass_out": 6.749999959546603,
            "vel_bed": 0.0048992662729692125,
            "vel_inter": 0.009798532545942183,
            "holdup": 102.66180257369356,
            "service_flow_rate": 10,
            "pressure_drop": 6.910024848393857,
            "Re": 3.4294863910767805,
            "Sc": 1086.9565217396594,
            "Sh": 32.88759444677707,
            "Pe_p": 0.09033997470879263,
            "Pe_bed": 227.62264675835323,
            "c_norm": 0.48299134878783817,
            "service_to_regen_flow_ratio": 3,
            "regen_recycle": 1,
            "regen_dose": 300,
            "t_regen": 1800,
            "bw_rate": 5,
            "bw_flow": 0.0014174458069060663,
            "t_bw": 600,
            "bed_expansion_frac": 0.46395000000000003,
            "bed_expansion_h": 0.8182852514438409,
            "rinse_bv": 5,
            "rinse_flow": 0.005000000000002573,
            "t_rinse": 899.9999999996298,
            "main_pump_power": 0.29784481412490327,
            "regen_pump_power": 0.0992816047082775,
            "bw_pump_power": 0.08443577678496997,
            "rinse_pump_power": 0.2978448141249117,
            "pump_efficiency": 0.8,
        }

        for v, val in results_dict.items():
            var = getattr(ix, v)
            if var.is_indexed():
                assert pytest.approx(val, rel=1e-3) == value(var[target_ion])
            else:
                assert pytest.approx(val, rel=1e-3) == value(var)
