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
"""
Tests for ASM2d-ADM1 Translator unit model.

Verified against approximated results results from:

Flores-Alsina, X., Solon, K., Mbamba, C.K., Tait, S., Gernaey, K.V., Jeppsson, U. and Batstone, D.J., 2016.
Modelling phosphorus (P), sulfur (S) and iron (Fe) interactions for dynamic simulations of anaerobic digestion processes.
Water Research, 95, pp.370-382.
"""

import pytest
from pyomo.environ import (
    ConcreteModel,
    value,
    assert_optimal_termination,
    Param,
    Objective,
    SolverFactory,
)

from idaes.core import FlowsheetBlock
from idaes.core.util.model_diagnostics import DegeneracyHunter
import idaes.core.util.scaling as iscale

from pyomo.environ import units

from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
    unused_variables_set,
    variables_set,
    fixed_variables_set,
    fixed_variables_generator,
    unfixed_variables_set,
    unfixed_variables_generator,
)

import idaes.logger as idaeslog
from idaes.core.util.testing import initialization_tester

from watertap.unit_models.translators.translator_asm2d_adm1 import Translator_ASM2d_ADM1
from watertap.property_models.anaerobic_digestion.modified_adm1_properties import (
    ModifiedADM1ParameterBlock,
)

from watertap.property_models.activated_sludge.modified_asm2d_properties import (
    ModifiedASM2dParameterBlock,
)

from watertap.property_models.anaerobic_digestion.modified_adm1_reactions import (
    ModifiedADM1ReactionParameterBlock,
)

from watertap.property_models.activated_sludge.modified_asm2d_reactions import (
    ModifiedASM2dReactionParameterBlock,
    DecaySwitch,
)

from watertap.core.util.model_diagnostics.infeasible import *

# from watertap.core.util.model_diagnostics.infeasible import print_close_to_bounds, print_infeasible_constraints

from pyomo.util.check_units import assert_units_consistent

from idaes.core.util.scaling import (
    calculate_scaling_factors,
    set_scaling_factor,
    get_scaling_factor,
    constraint_scaling_transform,
)

import idaes.core.util.scaling as iscale

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_config():
    m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.props_ASM2d = ModifiedASM2dParameterBlock()
    m.fs.props_ADM1 = ModifiedADM1ParameterBlock()
    m.fs.ADM1_rxn_props = ModifiedADM1ReactionParameterBlock(
        property_package=m.fs.props_ADM1
    )
    m.fs.ASM2d_rxn_props = ModifiedASM2dReactionParameterBlock(
        property_package=m.fs.props_ASM2d
    )

    m.fs.unit = Translator_ASM2d_ADM1(
        inlet_property_package=m.fs.props_ASM2d,
        outlet_property_package=m.fs.props_ADM1,
        inlet_reaction_package=m.fs.ASM2d_rxn_props,
        outlet_reaction_package=m.fs.ADM1_rxn_props,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )

    assert len(m.fs.unit.config) == 12

    assert m.fs.unit.config.outlet_state_defined == True
    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert not m.fs.unit.config.has_phase_equilibrium
    assert m.fs.unit.config.inlet_property_package is m.fs.props_ASM2d
    assert m.fs.unit.config.outlet_property_package is m.fs.props_ADM1
    assert m.fs.unit.config.inlet_reaction_package is m.fs.ASM2d_rxn_props
    assert m.fs.unit.config.outlet_reaction_package is m.fs.ADM1_rxn_props


# -----------------------------------------------------------------------------
class TestAsm2dAdm1_decay_on(object):
    @pytest.fixture(scope="class")
    def asmadm(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props_ASM2d = ModifiedASM2dParameterBlock()
        m.fs.props_ADM1 = ModifiedADM1ParameterBlock()
        m.fs.ADM1_rxn_props = ModifiedADM1ReactionParameterBlock(
            property_package=m.fs.props_ADM1
        )
        m.fs.ASM2d_rxn_props = ModifiedASM2dReactionParameterBlock(
            property_package=m.fs.props_ASM2d,
            decay_switch=DecaySwitch.on,
        )

        m.fs.unit = Translator_ASM2d_ADM1(
            inlet_property_package=m.fs.props_ASM2d,
            outlet_property_package=m.fs.props_ADM1,
            inlet_reaction_package=m.fs.ASM2d_rxn_props,
            outlet_reaction_package=m.fs.ADM1_rxn_props,
            has_phase_equilibrium=False,
            outlet_state_defined=True,
        )

        # TODO: Check influent flow_vol
        m.fs.unit.inlet.flow_vol.fix(178.4674 * units.m**3 / units.day)
        m.fs.unit.inlet.temperature.fix(308.15 * units.K)
        m.fs.unit.inlet.pressure.fix(1 * units.atm)
        eps = 1e-9 * units.kg / units.m**3

        m.fs.unit.inlet.conc_mass_comp[0, "S_O2"].fix(eps)
        m.fs.unit.inlet.conc_mass_comp[0, "S_F"].fix(0.02644 * units.kg / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_A"].fix(0.01766 * units.kg / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_I"].fix(0.02723 * units.kg / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_NH4"].fix(
            0.01858 * units.kg / units.m**3
        )
        m.fs.unit.inlet.conc_mass_comp[0, "S_N2"].fix(0.00507 * units.kg / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_NO3"].fix(
            0.00002 * units.kg / units.m**3
        )
        m.fs.unit.inlet.conc_mass_comp[0, "S_PO4"].fix(
            0.00469 * units.kg / units.m**3
        )
        m.fs.unit.inlet.conc_mass_comp[0, "S_IC"].fix(0.07899 * units.kg / units.m**3)

        m.fs.unit.inlet.conc_mass_comp[0, "X_I"].fix(10.96441 * units.kg / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_S"].fix(19.08476 * units.kg / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_H"].fix(9.47939 * units.kg / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_PAO"].fix(3.8622 * units.kg / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_PP"].fix(0.45087 * units.kg / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_PHA"].fix(
            0.02464 * units.kg / units.m**3
        )
        m.fs.unit.inlet.conc_mass_comp[0, "X_AUT"].fix(
            0.33379 * units.kg / units.m**3
        )
        m.fs.unit.inlet.conc_mass_comp[0, "S_K"].fix(0.01979 * units.kg / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_Mg"].fix(0.18987 * units.kg / units.m**3)

        # constraint_scaling_transform(m.fs.unit.SIC_output[0], 1e-3)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, asmadm):
        assert isinstance(asmadm.fs.unit.f_sI_xc, Param)
        assert value(asmadm.fs.unit.f_sI_xc) == 1e-9
        assert isinstance(asmadm.fs.unit.f_xI_xc, Param)
        assert value(asmadm.fs.unit.f_xI_xc) == 0.1
        assert isinstance(asmadm.fs.unit.f_ch_xc, Param)
        assert value(asmadm.fs.unit.f_ch_xc) == 0.275
        assert isinstance(asmadm.fs.unit.f_pr_xc, Param)
        assert value(asmadm.fs.unit.f_pr_xc) == 0.275
        assert isinstance(asmadm.fs.unit.f_li_xc, Param)
        assert value(asmadm.fs.unit.f_li_xc) == 0.35

        assert isinstance(asmadm.fs.unit.f_XPHA_Sva, Param)
        assert value(asmadm.fs.unit.f_XPHA_Sva) == 0.1
        assert isinstance(asmadm.fs.unit.f_XPHA_Sbu, Param)
        assert value(asmadm.fs.unit.f_XPHA_Sbu) == 0.1
        assert isinstance(asmadm.fs.unit.f_XPHA_Spro, Param)
        assert value(asmadm.fs.unit.f_XPHA_Spro) == 0.4
        assert isinstance(asmadm.fs.unit.f_XPHA_Sac, Param)
        assert value(asmadm.fs.unit.f_XPHA_Sac) == 0.4

        assert isinstance(asmadm.fs.unit.C_PHA, Param)
        assert value(asmadm.fs.unit.C_PHA) == 0.3 / 12

        assert hasattr(asmadm.fs.unit, "inlet")
        assert len(asmadm.fs.unit.inlet.vars) == 4
        assert hasattr(asmadm.fs.unit.inlet, "flow_vol")
        assert hasattr(asmadm.fs.unit.inlet, "conc_mass_comp")
        assert hasattr(asmadm.fs.unit.inlet, "temperature")
        assert hasattr(asmadm.fs.unit.inlet, "pressure")

        assert hasattr(asmadm.fs.unit, "outlet")
        assert len(asmadm.fs.unit.outlet.vars) == 6
        assert hasattr(asmadm.fs.unit.outlet, "flow_vol")
        assert hasattr(asmadm.fs.unit.outlet, "conc_mass_comp")
        assert hasattr(asmadm.fs.unit.outlet, "temperature")
        assert hasattr(asmadm.fs.unit.outlet, "pressure")
        assert hasattr(asmadm.fs.unit.outlet, "anions")
        assert hasattr(asmadm.fs.unit.outlet, "cations")

        assert number_variables(asmadm) == 265
        assert number_total_constraints(asmadm) == 34

        # TODO: Result of SN2_AS2 being unused. Remove? It's also unused in the c-code
        assert number_unused_variables(asmadm.fs.unit) == 1

    @pytest.mark.component
    def test_units(self, asmadm):
        assert_units_consistent(asmadm)

    @pytest.mark.unit
    def test_dof(self, asmadm):
        assert degrees_of_freedom(asmadm) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, asmadm):
        def check_jac(model):
            jac, jac_scaled, nlp = iscale.constraint_autoscale_large_jac(
                model, min_scale=1e-8
            )
            # cond_number = iscale.jacobian_cond(model, jac=jac_scaled)  # / 1e10
            # print("--------------------------")
            print("Extreme Jacobian entries:")
            extreme_entries = iscale.extreme_jacobian_entries(
                model, jac=jac_scaled, zero=1e-20, large=10
            )
            extreme_entries = sorted(extreme_entries, key=lambda x: x[0], reverse=True)

            print("EXTREME_ENTRIES")
            print(f"\nThere are {len(extreme_entries)} extreme Jacobian entries")
            for i in extreme_entries:
                print(i[0], i[1], i[2])

            print("--------------------------")
            print("Extreme Jacobian columns:")
            extreme_cols = iscale.extreme_jacobian_columns(model, jac=jac_scaled)
            for val, var in extreme_cols:
                print(val, var.name)
            print("------------------------")
            print("Extreme Jacobian rows:")
            extreme_rows = iscale.extreme_jacobian_rows(model, jac=jac_scaled)
            for val, con in extreme_rows:
                print(val, con.name)

        check_jac(asmadm)

        asmadm.obj = Objective(expr=0)

        # initial point
        solver.options["max_iter"] = 0
        solver.solve(asmadm, tee=False)
        dh = DegeneracyHunter(asmadm, solver=SolverFactory("cbc"))
        dh.check_residuals(tol=1e-8)
        dh.check_variable_bounds(tol=1e-8)

        # solved model
        solver.options["max_iter"] = 10000
        solver.solve(asmadm, tee=False)
        badly_scaled_var_list = iscale.badly_scaled_var_generator(
            asmadm, large=1e1, small=1e-1
        )
        for x in badly_scaled_var_list:
            print(f"{x[0].name}\t{x[0].value}\tsf: {iscale.get_scaling_factor(x[0])}")
        dh.check_residuals(tol=1e-8)
        dh.check_variable_bounds(tol=1e-8)
        dh.check_rank_equality_constraints(dense=True)
        ds = dh.find_candidate_equations(verbose=False, tee=False)
        ids = dh.find_irreducible_degenerate_sets(verbose=False)

        print("------------------------")
        print("Infeasible Constraints")
        print_close_to_bounds(asmadm)
        print_infeasible_constraints(asmadm)

        initialization_tester(asmadm)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, asmadm):

        asmadm.fs.unit.initialize(outlvl=idaeslog.INFO_HIGH)
        solver = get_solver()
        results = solver.solve(asmadm, tee=True)
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, asmadm):
        assert pytest.approx(101325.0, rel=1e-3) == value(
            asmadm.fs.unit.outlet.pressure[0]
        )
        assert pytest.approx(308.15, rel=1e-3) == value(
            asmadm.fs.unit.outlet.temperature[0]
        )
        assert pytest.approx(0.01843, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "S_su"]
        )
        assert pytest.approx(0.00801, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "S_aa"]
        )
        assert pytest.approx(1e-9, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "S_fa"]
        )
        assert pytest.approx(1e-10, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "S_va"]
        )
        assert pytest.approx(1e-10, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "S_bu"]
        )
        assert pytest.approx(4e-10, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "S_pro"]
        )
        assert pytest.approx(0.0175076, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "S_ac"]
        )
        assert pytest.approx(1e-9, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "S_h2"]
        )
        assert pytest.approx(1e-9, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "S_ch4"]
        )
        assert pytest.approx(0.541823, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "S_IC"]
        )
        assert pytest.approx(0.6985, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "S_IN"]
        )
        assert pytest.approx(0.71189, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "S_IP"]
        )
        assert pytest.approx(0.02723, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "S_I"]
        )

        assert pytest.approx(9.082, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "X_ch"]
        )
        assert pytest.approx(9.5424, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "X_pr"]
        )
        assert pytest.approx(12.7683, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "X_li"]
        )
        assert pytest.approx(1e-9, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "X_su"]
        )
        assert pytest.approx(1e-9, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "X_aa"]
        )
        assert pytest.approx(1e-9, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "X_fa"]
        )
        assert pytest.approx(1e-9, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "X_c4"]
        )
        assert pytest.approx(1e-9, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "X_pro"]
        )
        assert pytest.approx(1e-9, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "X_ac"]
        )
        assert pytest.approx(1e-9, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "X_h2"]
        )
        assert pytest.approx(12.3320, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "X_I"]
        )
        assert pytest.approx(1e-9, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "X_PHA"]
        )
        assert pytest.approx(1e-9, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "X_PP"]
        )
        assert pytest.approx(1e-9, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "X_PAO"]
        )
        assert pytest.approx(0.1701, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "S_K"]
        )
        assert pytest.approx(0.34016, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "S_Mg"]
        )
        assert pytest.approx(0.04515, rel=1e-3) == value(
            asmadm.fs.unit.outlet.cations[0]
        )
        assert pytest.approx(0.04989, rel=1e-3) == value(
            asmadm.fs.unit.outlet.anions[0]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, asmadm):
        assert (
            abs(
                value(
                    asmadm.fs.unit.inlet.flow_vol[0] * asmadm.fs.props_ADM1.dens_mass
                    - asmadm.fs.unit.outlet.flow_vol[0]
                    * asmadm.fs.props_ASM2d.dens_mass
                )
            )
            <= 1e-6
        )

        assert (
            abs(
                value(
                    (
                        asmadm.fs.unit.inlet.flow_vol[0]
                        * asmadm.fs.props_ADM1.dens_mass
                        * asmadm.fs.props_ADM1.cp_mass
                        * (
                            asmadm.fs.unit.inlet.temperature[0]
                            - asmadm.fs.props_ADM1.temperature_ref
                        )
                    )
                    - (
                        asmadm.fs.unit.outlet.flow_vol[0]
                        * asmadm.fs.props_ASM2d.dens_mass
                        * asmadm.fs.props_ASM2d.cp_mass
                        * (
                            asmadm.fs.unit.outlet.temperature[0]
                            - asmadm.fs.props_ASM2d.temperature_ref
                        )
                    )
                )
            )
            <= 1e-6
        )


# -----------------------------------------------------------------------------
class TestAsm2dAdm1_decay_off(object):
    @pytest.fixture(scope="class")
    def asmadm(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props_ASM2d = ModifiedASM2dParameterBlock()
        m.fs.props_ADM1 = ModifiedADM1ParameterBlock()
        m.fs.ADM1_rxn_props = ModifiedADM1ReactionParameterBlock(
            property_package=m.fs.props_ADM1
        )
        m.fs.ASM2d_rxn_props = ModifiedASM2dReactionParameterBlock(
            property_package=m.fs.props_ASM2d,
            decay_switch=DecaySwitch.off,
        )

        m.fs.unit = Translator_ASM2d_ADM1(
            inlet_property_package=m.fs.props_ASM2d,
            outlet_property_package=m.fs.props_ADM1,
            inlet_reaction_package=m.fs.ASM2d_rxn_props,
            outlet_reaction_package=m.fs.ADM1_rxn_props,
            has_phase_equilibrium=False,
            outlet_state_defined=True,
        )

        # TODO: Check influent flow_vol
        m.fs.unit.inlet.flow_vol.fix(178.4674 * units.m**3 / units.day)
        m.fs.unit.inlet.temperature.fix(308.15 * units.K)
        m.fs.unit.inlet.pressure.fix(1 * units.atm)
        eps = 1e-9 * units.kg / units.m**3

        m.fs.unit.inlet.conc_mass_comp[0, "S_O2"].fix(eps)
        m.fs.unit.inlet.conc_mass_comp[0, "S_F"].fix(0.02644 * units.kg / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_A"].fix(0.01766 * units.kg / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_I"].fix(0.02723 * units.kg / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_NH4"].fix(
            0.01858 * units.kg / units.m**3
        )
        m.fs.unit.inlet.conc_mass_comp[0, "S_N2"].fix(0.00507 * units.kg / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_NO3"].fix(
            0.00002 * units.kg / units.m**3
        )
        m.fs.unit.inlet.conc_mass_comp[0, "S_PO4"].fix(
            0.00469 * units.kg / units.m**3
        )
        m.fs.unit.inlet.conc_mass_comp[0, "S_IC"].fix(0.07899 * units.kg / units.m**3)

        m.fs.unit.inlet.conc_mass_comp[0, "X_I"].fix(10.96441 * units.kg / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_S"].fix(19.08476 * units.kg / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_H"].fix(9.47939 * units.kg / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_PAO"].fix(3.8622 * units.kg / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_PP"].fix(0.45087 * units.kg / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_PHA"].fix(
            0.02464 * units.kg / units.m**3
        )
        m.fs.unit.inlet.conc_mass_comp[0, "X_AUT"].fix(
            0.33379 * units.kg / units.m**3
        )
        m.fs.unit.inlet.conc_mass_comp[0, "S_K"].fix(0.01979 * units.kg / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_Mg"].fix(0.18987 * units.kg / units.m**3)

        # constraint_scaling_transform(m.fs.unit.SIC_output[0], 1e-3)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, asmadm):
        assert isinstance(asmadm.fs.unit.f_sI_xc, Param)
        assert value(asmadm.fs.unit.f_sI_xc) == 1e-9
        assert isinstance(asmadm.fs.unit.f_xI_xc, Param)
        assert value(asmadm.fs.unit.f_xI_xc) == 0.1
        assert isinstance(asmadm.fs.unit.f_ch_xc, Param)
        assert value(asmadm.fs.unit.f_ch_xc) == 0.275
        assert isinstance(asmadm.fs.unit.f_pr_xc, Param)
        assert value(asmadm.fs.unit.f_pr_xc) == 0.275
        assert isinstance(asmadm.fs.unit.f_li_xc, Param)
        assert value(asmadm.fs.unit.f_li_xc) == 0.35

        assert isinstance(asmadm.fs.unit.f_XPHA_Sva, Param)
        assert value(asmadm.fs.unit.f_XPHA_Sva) == 0.1
        assert isinstance(asmadm.fs.unit.f_XPHA_Sbu, Param)
        assert value(asmadm.fs.unit.f_XPHA_Sbu) == 0.1
        assert isinstance(asmadm.fs.unit.f_XPHA_Spro, Param)
        assert value(asmadm.fs.unit.f_XPHA_Spro) == 0.4
        assert isinstance(asmadm.fs.unit.f_XPHA_Sac, Param)
        assert value(asmadm.fs.unit.f_XPHA_Sac) == 0.4

        assert isinstance(asmadm.fs.unit.C_PHA, Param)
        assert value(asmadm.fs.unit.C_PHA) == 0.3 / 12

        assert hasattr(asmadm.fs.unit, "inlet")
        assert len(asmadm.fs.unit.inlet.vars) == 4
        assert hasattr(asmadm.fs.unit.inlet, "flow_vol")
        assert hasattr(asmadm.fs.unit.inlet, "conc_mass_comp")
        assert hasattr(asmadm.fs.unit.inlet, "temperature")
        assert hasattr(asmadm.fs.unit.inlet, "pressure")

        assert hasattr(asmadm.fs.unit, "outlet")
        assert len(asmadm.fs.unit.outlet.vars) == 6
        assert hasattr(asmadm.fs.unit.outlet, "flow_vol")
        assert hasattr(asmadm.fs.unit.outlet, "conc_mass_comp")
        assert hasattr(asmadm.fs.unit.outlet, "temperature")
        assert hasattr(asmadm.fs.unit.outlet, "pressure")
        assert hasattr(asmadm.fs.unit.outlet, "anions")
        assert hasattr(asmadm.fs.unit.outlet, "cations")

        assert number_variables(asmadm) == 265
        assert number_total_constraints(asmadm) == 34

        # TODO: Result of SN2_AS2 being unused. Remove? It's also unused in the c-code
        assert number_unused_variables(asmadm.fs.unit) == 1

    @pytest.mark.component
    def test_units(self, asmadm):
        assert_units_consistent(asmadm)

    @pytest.mark.unit
    def test_dof(self, asmadm):
        assert degrees_of_freedom(asmadm) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, asmadm):
        def check_jac(model):
            jac, jac_scaled, nlp = iscale.constraint_autoscale_large_jac(
                model, min_scale=1e-8
            )
            # cond_number = iscale.jacobian_cond(model, jac=jac_scaled)  # / 1e10
            # print("--------------------------")
            print("Extreme Jacobian entries:")
            extreme_entries = iscale.extreme_jacobian_entries(
                model, jac=jac_scaled, zero=1e-20, large=10
            )
            extreme_entries = sorted(extreme_entries, key=lambda x: x[0], reverse=True)

            print("EXTREME_ENTRIES")
            print(f"\nThere are {len(extreme_entries)} extreme Jacobian entries")
            for i in extreme_entries:
                print(i[0], i[1], i[2])

            print("--------------------------")
            print("Extreme Jacobian columns:")
            extreme_cols = iscale.extreme_jacobian_columns(model, jac=jac_scaled)
            for val, var in extreme_cols:
                print(val, var.name)
            print("------------------------")
            print("Extreme Jacobian rows:")
            extreme_rows = iscale.extreme_jacobian_rows(model, jac=jac_scaled)
            for val, con in extreme_rows:
                print(val, con.name)

        check_jac(asmadm)

        asmadm.obj = Objective(expr=0)

        # initial point
        solver.options["max_iter"] = 0
        solver.solve(asmadm, tee=False)
        dh = DegeneracyHunter(asmadm, solver=SolverFactory("cbc"))
        dh.check_residuals(tol=1e-8)
        dh.check_variable_bounds(tol=1e-8)

        # solved model
        solver.options["max_iter"] = 10000
        solver.solve(asmadm, tee=False)
        badly_scaled_var_list = iscale.badly_scaled_var_generator(
            asmadm, large=1e1, small=1e-1
        )
        for x in badly_scaled_var_list:
            print(f"{x[0].name}\t{x[0].value}\tsf: {iscale.get_scaling_factor(x[0])}")
        dh.check_residuals(tol=1e-8)
        dh.check_variable_bounds(tol=1e-8)
        dh.check_rank_equality_constraints(dense=True)
        ds = dh.find_candidate_equations(verbose=False, tee=False)
        ids = dh.find_irreducible_degenerate_sets(verbose=False)

        print("------------------------")
        print("Infeasible Constraints")
        print_close_to_bounds(asmadm)
        print_infeasible_constraints(asmadm)

        initialization_tester(asmadm)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, asmadm):

        asmadm.fs.unit.initialize(outlvl=idaeslog.INFO_HIGH)
        solver = get_solver()
        results = solver.solve(asmadm, tee=True)
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, asmadm):
        assert pytest.approx(101325.0, rel=1e-3) == value(
            asmadm.fs.unit.outlet.pressure[0]
        )
        assert pytest.approx(308.15, rel=1e-3) == value(
            asmadm.fs.unit.outlet.temperature[0]
        )
        assert pytest.approx(0.01843, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "S_su"]
        )
        assert pytest.approx(0.00801, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "S_aa"]
        )
        assert pytest.approx(1e-9, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "S_fa"]
        )
        assert pytest.approx(9.9999e-10, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "S_va"]
        )
        assert pytest.approx(9.9999e-10, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "S_bu"]
        )
        assert pytest.approx(9.9999e-10, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "S_pro"]
        )
        assert pytest.approx(0.0175076, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "S_ac"]
        )
        assert pytest.approx(1e-9, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "S_h2"]
        )
        assert pytest.approx(1e-9, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "S_ch4"]
        )
        assert pytest.approx(0.41266, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "S_IC"]
        )
        assert pytest.approx(0.50648, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "S_IN"]
        )
        assert pytest.approx(0.19475, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "S_IP"]
        )
        assert pytest.approx(0.02723, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "S_I"]
        )

        assert pytest.approx(8.0199, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "X_ch"]
        )
        assert pytest.approx(8.48026, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "X_pr"]
        )
        assert pytest.approx(11.4165, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "X_li"]
        )
        assert pytest.approx(1e-9, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "X_su"]
        )
        assert pytest.approx(1e-9, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "X_aa"]
        )
        assert pytest.approx(1e-9, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "X_fa"]
        )
        assert pytest.approx(1e-9, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "X_c4"]
        )
        assert pytest.approx(1e-9, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "X_pro"]
        )
        assert pytest.approx(1e-9, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "X_ac"]
        )
        assert pytest.approx(1e-9, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "X_h2"]
        )
        assert pytest.approx(11.94574, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "X_I"]
        )
        assert pytest.approx(0.02464, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "X_PHA"]
        )
        assert pytest.approx(0.45087, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "X_PP"]
        )
        assert pytest.approx(3.8622, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "X_PAO"]
        )
        assert pytest.approx(0.01979, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "S_K"]
        )
        assert pytest.approx(0.18987, rel=1e-3) == value(
            asmadm.fs.unit.outlet.conc_mass_comp[0, "S_Mg"]
        )
        assert pytest.approx(0.03439, rel=1e-3) == value(
            asmadm.fs.unit.outlet.cations[0]
        )
        assert pytest.approx(0.03618, rel=1e-3) == value(
            asmadm.fs.unit.outlet.anions[0]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, asmadm):
        assert (
            abs(
                value(
                    asmadm.fs.unit.inlet.flow_vol[0] * asmadm.fs.props_ADM1.dens_mass
                    - asmadm.fs.unit.outlet.flow_vol[0]
                    * asmadm.fs.props_ASM2d.dens_mass
                )
            )
            <= 1e-6
        )

        assert (
            abs(
                value(
                    (
                        asmadm.fs.unit.inlet.flow_vol[0]
                        * asmadm.fs.props_ADM1.dens_mass
                        * asmadm.fs.props_ADM1.cp_mass
                        * (
                            asmadm.fs.unit.inlet.temperature[0]
                            - asmadm.fs.props_ADM1.temperature_ref
                        )
                    )
                    - (
                        asmadm.fs.unit.outlet.flow_vol[0]
                        * asmadm.fs.props_ASM2d.dens_mass
                        * asmadm.fs.props_ASM2d.cp_mass
                        * (
                            asmadm.fs.unit.outlet.temperature[0]
                            - asmadm.fs.props_ASM2d.temperature_ref
                        )
                    )
                )
            )
            <= 1e-6
        )
