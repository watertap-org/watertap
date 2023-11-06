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
"""
Based on flowsheet from:

Flores-Alsina X., Gernaey K.V. and Jeppsson, U. "Benchmarking biological
nutrient removal in wastewater treatment plants: influence of mathematical model
assumptions", 2012, Wat. Sci. Tech., Vol. 65 No. 8, pp. 1496-1505
"""

# Some more information about this module
__author__ = "Alejandro Garciadiego, Andrew Lee"

import pyomo.environ as pyo
from pyomo.environ import (
    value,
    units as pyunits,
)
from pyomo.network import Arc, SequentialDecomposition

from idaes.core import (
    FlowsheetBlock,
    UnitModelCostingBlock,
)
from idaes.models.unit_models import (
    CSTR,
    Feed,
    Separator,
    Product,
    Mixer,
    PressureChanger,
)
from idaes.models.unit_models.separator import SplittingType
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale
from idaes.core.util.tables import (
    create_stream_table_dataframe,
    stream_table_dataframe_to_string,
)
from watertap.unit_models.cstr_injection import CSTR_Injection
from watertap.property_models.anaerobic_digestion.modified_adm1_properties import (
    ModifiedADM1ParameterBlock,
)
from watertap.property_models.anaerobic_digestion.adm1_properties_vapor import (
    ADM1_vaporParameterBlock,
)
from watertap.property_models.anaerobic_digestion.modified_adm1_reactions import (
    ModifiedADM1ReactionParameterBlock,
)
from watertap.property_models.activated_sludge.modified_asm2d_properties import (
    ModifiedASM2dParameterBlock,
)
from watertap.property_models.activated_sludge.modified_asm2d_reactions import (
    ModifiedASM2dReactionParameterBlock,
)
from watertap.unit_models.translators.translator_adm1_asm2d import (
    Translator_ADM1_ASM2D,
)
from watertap.unit_models.translators.translator_asm2d_adm1 import Translator_ASM2d_ADM1
from watertap.unit_models.anaerobic_digestor import AD
from watertap.unit_models.electroNP_ZO import ElectroNPZO
from watertap.unit_models.dewatering import (
    DewateringUnit,
    ActivatedSludgeModelType as dewater_type,
)
from watertap.unit_models.thickener import (
    Thickener,
    ActivatedSludgeModelType as thickener_type,
)
from watertap.core.util.initialization import check_solve
from watertap.costing import WaterTAPCosting

from watertap.core.util.model_diagnostics.infeasible import *
from idaes.core.util.model_diagnostics import DegeneracyHunter


# Set up logger
_log = idaeslog.getLogger(__name__)


def build_flowsheet():
    m = pyo.ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.props_ASM2D = ModifiedASM2dParameterBlock()
    m.fs.rxn_props_ASM2D = ModifiedASM2dReactionParameterBlock(
        property_package=m.fs.props_ASM2D
    )
    m.fs.props_ADM1 = ModifiedADM1ParameterBlock()
    m.fs.props_vap_ADM1 = ADM1_vaporParameterBlock()
    m.fs.rxn_props_ADM1 = ModifiedADM1ReactionParameterBlock(
        property_package=m.fs.props_ADM1
    )

    m.fs.costing = WaterTAPCosting()

    # Feed water stream
    m.fs.FeedWater = Feed(property_package=m.fs.props_ASM2D)

    # Translators
    m.fs.translator_asm2d_adm1 = Translator_ASM2d_ADM1(
        inlet_property_package=m.fs.props_ASM2D,
        outlet_property_package=m.fs.props_ADM1,
        inlet_reaction_package=m.fs.rxn_props_ASM2D,
        outlet_reaction_package=m.fs.rxn_props_ADM1,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )

    # Link units
    m.fs.stream1 = Arc(
        source=m.fs.FeedWater.outlet, destination=m.fs.translator_asm2d_adm1.inlet
    )

    pyo.TransformationFactory("network.expand_arcs").apply_to(m)

    # m.fs.FeedWater.flow_vol.fix(20648 * pyo.units.m**3 / pyo.units.day)
    # m.fs.FeedWater.temperature.fix(308.15 * pyo.units.K)
    # m.fs.FeedWater.pressure.fix(1 * pyo.units.atm)
    # m.fs.FeedWater.conc_mass_comp[0, "S_O2"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    # m.fs.FeedWater.conc_mass_comp[0, "S_F"].fix(30 * pyo.units.g / pyo.units.m**3)
    # m.fs.FeedWater.conc_mass_comp[0, "S_A"].fix(20 * pyo.units.g / pyo.units.m**3)
    # m.fs.FeedWater.conc_mass_comp[0, "S_NH4"].fix(16 * pyo.units.g / pyo.units.m**3)
    # m.fs.FeedWater.conc_mass_comp[0, "S_NO3"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    # m.fs.FeedWater.conc_mass_comp[0, "S_PO4"].fix(3.6 * pyo.units.g / pyo.units.m**3)
    # m.fs.FeedWater.conc_mass_comp[0, "S_I"].fix(30 * pyo.units.g / pyo.units.m**3)
    # m.fs.FeedWater.conc_mass_comp[0, "S_N2"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    # m.fs.FeedWater.conc_mass_comp[0, "X_I"].fix(25 * pyo.units.g / pyo.units.m**3)
    # m.fs.FeedWater.conc_mass_comp[0, "X_S"].fix(125 * pyo.units.g / pyo.units.m**3)
    # m.fs.FeedWater.conc_mass_comp[0, "X_H"].fix(30 * pyo.units.g / pyo.units.m**3)
    # m.fs.FeedWater.conc_mass_comp[0, "X_PAO"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    # m.fs.FeedWater.conc_mass_comp[0, "X_PP"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    # m.fs.FeedWater.conc_mass_comp[0, "X_PHA"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    # m.fs.FeedWater.conc_mass_comp[0, "X_AUT"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    # m.fs.FeedWater.properties[0].conc_mass_comp["S_IC"].fix(
    #     0.07899 * pyo.units.kg / pyo.units.m**3
    # )
    # m.fs.FeedWater.properties[0].conc_mass_comp["S_K"].fix(
    #     1e-6 * pyo.units.g / pyo.units.m**3
    # )
    # m.fs.FeedWater.properties[0].conc_mass_comp["S_Mg"].fix(
    #     1e-6 * pyo.units.g / pyo.units.m**3
    # )

    m.fs.FeedWater.flow_vol.fix(0.00018359)
    m.fs.FeedWater.temperature.fix(308.15 * pyo.units.K)
    m.fs.FeedWater.pressure.fix(1 * pyo.units.atm)
    m.fs.FeedWater.conc_mass_comp[0, "S_O2"].fix(0.0079024)
    m.fs.FeedWater.conc_mass_comp[0, "S_F"].fix(0.00030596)
    m.fs.FeedWater.conc_mass_comp[0, "S_A"].fix(4.1385e-05)
    m.fs.FeedWater.conc_mass_comp[0, "S_NH4"].fix(0.0070896)
    m.fs.FeedWater.conc_mass_comp[0, "S_NO3"].fix(1.1750e-07)
    m.fs.FeedWater.conc_mass_comp[0, "S_PO4"].fix(0.0027761)
    m.fs.FeedWater.conc_mass_comp[0, "S_I"].fix(0.03)
    m.fs.FeedWater.conc_mass_comp[0, "S_N2"].fix(5.2387e-07)
    m.fs.FeedWater.conc_mass_comp[0, "X_I"].fix(19.794)
    m.fs.FeedWater.conc_mass_comp[0, "X_S"].fix(1.2980)
    m.fs.FeedWater.conc_mass_comp[0, "X_H"].fix(64.284)
    m.fs.FeedWater.conc_mass_comp[0, "X_PAO"].fix(3.7115e-07)
    m.fs.FeedWater.conc_mass_comp[0, "X_PP"].fix(7.3446e-07)
    m.fs.FeedWater.conc_mass_comp[0, "X_PHA"].fix(1.2442e-06)
    m.fs.FeedWater.conc_mass_comp[0, "X_AUT"].fix(7.6090e-05)
    m.fs.FeedWater.properties[0].conc_mass_comp["S_IC"].fix(0.10470)
    m.fs.FeedWater.properties[0].conc_mass_comp["S_K"].fix(2.4407e-08)
    m.fs.FeedWater.properties[0].conc_mass_comp["S_Mg"].fix(2.4425e-08)

    # Check degrees of freedom
    print(f"DOF after all: {degrees_of_freedom(m)}")
    assert degrees_of_freedom(m) == 0

    def scale_variables(m):
        for var in m.fs.component_data_objects(pyo.Var, descend_into=True):
            if "flow_vol" in var.name:
                iscale.set_scaling_factor(var, 1e1)
            if "translator_asm2d_adm1.properties_in[0.0].flow_vol" in var.name:
                iscale.set_scaling_factor(var, 1e4)
            if "translator_asm2d_adm1.properties_out[0.0].flow_vol" in var.name:
                iscale.set_scaling_factor(var, 1e4)
            if "temperature" in var.name:
                iscale.set_scaling_factor(var, 1e-2)
            if "pressure" in var.name:
                iscale.set_scaling_factor(var, 1e-3)
            # if "conc_mass_comp" in var.name:
            #     iscale.set_scaling_factor(var, 1e2)
            # if "conc_mass_comp[S_A]" in var.name:
            #     iscale.set_scaling_factor(var, 1e6)
            # if "conc_mass_comp[S_F]" in var.name:
            #     iscale.set_scaling_factor(var, 1e4)
            # if "conc_mass_comp[S_I]" in var.name:
            #     iscale.set_scaling_factor(var, 1e2)
            # if "conc_mass_comp[S_N2]" in var.name:
            #     iscale.set_scaling_factor(var, 1e6)
            # if "conc_mass_comp[S_NH4]" in var.name:
            #     iscale.set_scaling_factor(var, 1e3)
            # if "conc_mass_comp[S_NO3]" in var.name:
            #     iscale.set_scaling_factor(var, 1e6)
            # if "conc_mass_comp[S_O2]" in var.name:
            #     iscale.set_scaling_factor(var, 1e6)
            # if "conc_mass_comp[S_PO4]" in var.name:
            #     iscale.set_scaling_factor(var, 1e3)
            # if "conc_mass_comp[S_K]" in var.name:
            #     iscale.set_scaling_factor(var, 1e6)
            # if "conc_mass_comp[S_Mg]" in var.name:
            #     iscale.set_scaling_factor(var, 1e6)
            # if "conc_mass_comp[S_IC]" in var.name:
            #     iscale.set_scaling_factor(var, 1e1)
            # if "conc_mass_comp[X_AUT]" in var.name:
            #     iscale.set_scaling_factor(var, 1e5)
            # if "conc_mass_comp[X_H]" in var.name:
            #     iscale.set_scaling_factor(var, 1e-1)
            # if "conc_mass_comp[X_I]" in var.name:
            #     iscale.set_scaling_factor(var, 1e-1)
            # if "conc_mass_comp[X_PAO]" in var.name:
            #     iscale.set_scaling_factor(var, 1e6)
            # if "conc_mass_comp[X_PHA]" in var.name:
            #     iscale.set_scaling_factor(var, 1e6)
            # if "conc_mass_comp[X_PP]" in var.name:
            #     iscale.set_scaling_factor(var, 1e6)
            # if "conc_mass_comp[X_S]" in var.name:
            #     iscale.set_scaling_factor(var, 1e0)

    # Apply scaling
    scale_variables(m)
    iscale.calculate_scaling_factors(m)

    # Initialize flowsheet
    # Apply sequential decomposition - 1 iteration should suffice
    seq = SequentialDecomposition()
    seq.options.select_tear_method = "heuristic"
    seq.options.tear_method = "Wegstein"
    seq.options.iterLim = 1

    G = seq.create_graph(m)

    tear_guesses2 = {
        "flow_vol": {0: 1.8e-4},
        "conc_mass_comp": {
            (0, "S_A"): 4.1e-5,
            (0, "S_F"): 0.0003,
            (0, "S_I"): 0.03,
            (0, "S_N2"): 1e-9,
            (0, "S_NH4"): 0.007,
            (0, "S_NO3"): 1e-9,
            (0, "S_O2"): 0.008,
            (0, "S_PO4"): 0.003,
            (0, "S_K"): 1e-9,
            (0, "S_Mg"): 1e-9,
            (0, "S_IC"): 0.1,
            (0, "X_AUT"): 1.8e-5,
            (0, "X_H"): 64.3,
            (0, "X_I"): 19.8,
            (0, "X_PAO"): 1e-9,
            (0, "X_PHA"): 1e-9,
            (0, "X_PP"): 1e-9,
            (0, "X_S"): 1.3,
        },
        "temperature": {0: 308.15},
        "pressure": {0: 101325},
    }

    # Pass the tear_guess to the SD tool
    seq.set_guesses_for(m.fs.translator_asm2d_adm1.inlet, tear_guesses2)

    def function(unit):
        unit.initialize(outlvl=idaeslog.INFO, optarg={"bound_push": 1e-2})
        # badly_scaled_vars = list(iscale.badly_scaled_var_generator(unit))
        # if len(badly_scaled_vars) > 0:
        #     automate_rescale_variables(unit)

    # # Model debug
    # solver = get_solver()
    # m, results = model_debug(m, solver)

    # Diagnostics Toolbox
    results = DiagnosticsToolbox(m)
    results.report_structural_issues()
    results.report_numerical_issues()
    results.display_constraints_with_large_residuals()
    results.display_variables_at_or_outside_bounds()
    results.display_constraints_with_extreme_jacobians()
    results.display_extreme_jacobian_entries()
    results.display_variables_near_bounds()

    # seq.run(m, function)

    # solver = get_solver()
    # results = solver.solve(m, tee=True)
    # check_solve(results, checkpoint="closing recycle", logger=_log, fail_flag=True)

    return m, results


def model_debug(model, solver):

    # check_jac(model)

    model.obj = pyo.Objective(expr=0)
    dh = DegeneracyHunter(model, solver=pyo.SolverFactory("cbc"))

    # # initial point
    # solver.options["max_iter"] = 0
    # results = solver.solve(model, tee=False)
    # dh.check_residuals(tol=1e-8)
    # # dh.check_variable_bounds(tol=1e-8)

    # solved model
    solver.options["max_iter"] = 10000
    results = solver.solve(model, tee=False)
    # badly_scaled_var_list = iscale.badly_scaled_var_generator(model, large=1e1, small=1e-1)
    # for x in badly_scaled_var_list:
    #     print(f"{x[0].name}\t{x[0].value}\tsf: {iscale.get_scaling_factor(x[0])}")
    dh.check_residuals(tol=1e-8)
    # dh.check_variable_bounds(tol=1e-8)
    # dh.check_rank_equality_constraints(dense=True)
    # ds = dh.find_candidate_equations(verbose=True, tee=True)
    # ids = dh.find_irreducible_degenerate_sets(verbose=True)

    """
    variables_near_bounds_list = variables_near_bounds_generator(model)
    for x in variables_near_bounds_list:
        print(x, x.value)
    """

    return model, results


def check_jac(model):
    jac, jac_scaled, nlp = iscale.constraint_autoscale_large_jac(model, min_scale=1e-8)
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


if __name__ == "__main__":
    # This method builds and runs a steady state activated sludge
    # flowsheet.
    m, results = build_flowsheet()
    # print_close_to_bounds(m)
    # print_infeasible_constraints(m)

    stream_table = create_stream_table_dataframe(
        {
            # "Feed": m.fs.FeedWater.outlet,
            # "Mix": m.fs.R1.inlet,
            # "R1": m.fs.R1.outlet,
            # "R2": m.fs.R2.outlet,
            # "R3 inlet": m.fs.R3.inlet,
            # "R3": m.fs.R3.outlet,
            # "R4": m.fs.R4.outlet,
            # "R5": m.fs.R5.outlet,
            # "R6": m.fs.R6.outlet,
            # "R7": m.fs.R7.outlet,
            # "thickener outlet": m.fs.thickener.underflow,
            "ASM-ADM translator inlet": m.fs.translator_asm2d_adm1.inlet,
            "ASM-ADM translator outlet": m.fs.translator_asm2d_adm1.outlet,
        },
        time_point=0,
    )
    # print(stream_table_dataframe_to_string(stream_table))
