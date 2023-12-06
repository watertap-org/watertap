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
from idaes.core.util.model_diagnostics import DiagnosticsToolbox

# Set up logger
_log = idaeslog.getLogger(__name__)


def automate_rescale_variables(m):
    for var, sv in iscale.badly_scaled_var_generator(m):
        if iscale.get_scaling_factor(var) is None:
            continue
        sf = iscale.get_scaling_factor(var)
        iscale.set_scaling_factor(var, sf / sv)
        iscale.calculate_scaling_factors(m)


def main():
    m = build_flowsheet()
    set_operating_conditions(m)

    results = initialize_system(m)

    results = solve(m)

    # try:
    #     results = solve(m)
    # except:
    #     pass

    # print_close_to_bounds(m)
    # print_infeasible_constraints(m)

    # # Use of Degeneracy Hunter for troubleshooting model.
    # m.obj = pyo.Objective(expr=0)
    # solver = get_solver()
    # solver.options["max_iter"] = 10000
    # results = solver.solve(m, tee=True)
    # dh = DegeneracyHunter(m, solver=pyo.SolverFactory("cbc"))
    # # badly_scaled_var_list = iscale.badly_scaled_var_generator(
    # #     m, large=1e1, small=1e-1
    # # )
    # # for x in badly_scaled_var_list:
    # #     print(f"{x[0].name}\t{x[0].value}\tsf: {iscale.get_scaling_factor(x[0])}")
    # dh.check_residuals(tol=1e-8)
    # # dh.check_variable_bounds(tol=1e-8)
    # # dh.check_rank_equality_constraints(dense=True)
    # # ds = dh.find_candidate_equations(verbose=True, tee=True)
    # # ids = dh.find_irreducible_degenerate_sets(verbose=True)

    # # Switch to fixed KLa in R3 and R4 (S_O concentration is controlled in R5)
    # m.fs.R5.KLa.fix(240)
    # m.fs.R6.KLa.fix(240)
    # m.fs.R7.KLa.fix(84)
    # m.fs.R5.outlet.conc_mass_comp[:, "S_O2"].unfix()
    # m.fs.R6.outlet.conc_mass_comp[:, "S_O2"].unfix()
    # m.fs.R7.outlet.conc_mass_comp[:, "S_O2"].unfix()
    # # Resolve with controls in place
    # results = solver.solve(m, tee=False)

    # pyo.assert_optimal_termination(results)
    # check_solve(
    #     results,
    #     checkpoint="re-solve with controls in place",
    #     logger=_log,
    #     fail_flag=True,
    # )

    # print("Numerical issues after solving")
    # dt.report_numerical_issues()

    # add_costing(m)
    # # Assert DOF = 0 after adding costing
    # # assert_degrees_of_freedom(m, 0)
    #
    # # TODO: initialize costing after adding to flowsheet
    # # m.fs.costing.initialize()
    #
    # # results = solve(m)

    # display_results(m)

    return m, results


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

    # ====================================================================
    # Primary Clarifier
    m.fs.CL = Separator(
        property_package=m.fs.props_ASM2D,
        outlet_list=["underflow", "effluent"],
        split_basis=SplittingType.componentFlow,
    )

    # ======================================================================
    # Activated Sludge Process
    # Mixer for feed water and recycled sludge
    m.fs.MX1 = Mixer(
        property_package=m.fs.props_ASM2D, inlet_list=["feed_water", "recycle"]
    )
    # First reactor (anoxic) - standard CSTR
    m.fs.R1 = CSTR(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    # First reactor (anoxic) - standard CSTR
    m.fs.R2 = CSTR(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    # Second reactor (anoxic) - standard CSTR
    m.fs.R3 = CSTR(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )

    m.fs.R4 = CSTR(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    # Third reactor (aerobic) - CSTR with injection
    m.fs.R5 = CSTR_Injection(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    # Fourth reactor (aerobic) - CSTR with injection
    m.fs.R6 = CSTR_Injection(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    # Fifth reactor (aerobic) - CSTR with injection
    m.fs.R7 = CSTR_Injection(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    m.fs.SP1 = Separator(
        property_package=m.fs.props_ASM2D, outlet_list=["underflow", "overflow"]
    )
    # Secondary Clarifier
    # TODO: Replace with more detailed model when available
    m.fs.CL1 = Separator(
        property_package=m.fs.props_ASM2D,
        outlet_list=["underflow", "effluent"],
        split_basis=SplittingType.componentFlow,
    )
    # Mixing sludge recycle and R5 underflow
    m.fs.MX2 = Mixer(
        property_package=m.fs.props_ASM2D, inlet_list=["clarifier", "reactor"]
    )
    # Sludge separator
    m.fs.SP2 = Separator(
        property_package=m.fs.props_ASM2D, outlet_list=["waste", "recycle"]
    )

    # ======================================================================
    # Thickener
    m.fs.thickener = Thickener(
        property_package=m.fs.props_ASM2D,
        activated_sludge_model=thickener_type.modified_ASM2D,
    )

    # ======================================================================
    # Anaerobic digester section
    # Translators
    m.fs.translator_asm2d_adm1 = Translator_ASM2d_ADM1(
        inlet_property_package=m.fs.props_ASM2D,
        outlet_property_package=m.fs.props_ADM1,
        inlet_reaction_package=m.fs.rxn_props_ASM2D,
        outlet_reaction_package=m.fs.rxn_props_ADM1,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )

    # Anaerobic digestor
    m.fs.AD = AD(
        liquid_property_package=m.fs.props_ADM1,
        vapor_property_package=m.fs.props_vap_ADM1,
        reaction_package=m.fs.rxn_props_ADM1,
        has_heat_transfer=True,
        has_pressure_change=False,
    )

    m.fs.translator_adm1_asm2d = Translator_ADM1_ASM2D(
        inlet_property_package=m.fs.props_ADM1,
        outlet_property_package=m.fs.props_ASM2D,
        reaction_package=m.fs.rxn_props_ADM1,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )

    # Dewatering Unit
    m.fs.dewater = DewateringUnit(
        property_package=m.fs.props_ASM2D,
        activated_sludge_model=dewater_type.modified_ASM2D,
    )

    # ElectroNP
    m.fs.electroNP = ElectroNPZO(property_package=m.fs.props_ASM2D)

    m.fs.MX3 = Mixer(
        property_package=m.fs.props_ASM2D, inlet_list=["feed_water1", "recycle1"]
    )
    # m.fs.MX4 = Mixer(
    #     property_package=m.fs.props_ASM2D, inlet_list=["feed_water2", "recycle2"]
    # )
    m.fs.MX5 = Mixer(
        property_package=m.fs.props_ASM2D, inlet_list=["thickener", "clarifier"]
    )

    # Product Blocks
    m.fs.Treated = Product(property_package=m.fs.props_ASM2D)
    m.fs.Sludge = Product(property_package=m.fs.props_ASM2D)
    # Recycle pressure changer - use a simple isothermal unit for now
    m.fs.P1 = PressureChanger(property_package=m.fs.props_ASM2D)

    # Link units related to ASM section
    m.fs.stream2 = Arc(source=m.fs.MX1.outlet, destination=m.fs.R1.inlet)
    m.fs.stream3 = Arc(source=m.fs.R1.outlet, destination=m.fs.R2.inlet)
    m.fs.stream4 = Arc(source=m.fs.R2.outlet, destination=m.fs.MX2.reactor)
    m.fs.stream5 = Arc(source=m.fs.MX2.outlet, destination=m.fs.R3.inlet)
    m.fs.stream6 = Arc(source=m.fs.R3.outlet, destination=m.fs.R4.inlet)
    m.fs.stream7 = Arc(source=m.fs.R4.outlet, destination=m.fs.R5.inlet)
    m.fs.stream8 = Arc(source=m.fs.R5.outlet, destination=m.fs.R6.inlet)
    m.fs.stream9 = Arc(source=m.fs.R6.outlet, destination=m.fs.R7.inlet)
    m.fs.stream10 = Arc(source=m.fs.R7.outlet, destination=m.fs.SP1.inlet)
    m.fs.stream11 = Arc(source=m.fs.SP1.overflow, destination=m.fs.CL1.inlet)
    m.fs.stream12 = Arc(source=m.fs.SP1.underflow, destination=m.fs.MX2.clarifier)
    m.fs.stream13 = Arc(source=m.fs.CL1.effluent, destination=m.fs.Treated.inlet)
    m.fs.stream14 = Arc(source=m.fs.CL1.underflow, destination=m.fs.SP2.inlet)
    # m.fs.stream15 = Arc(source=m.fs.SP2.waste, destination=m.fs.Sludge.inlet)
    m.fs.stream16 = Arc(source=m.fs.SP2.recycle, destination=m.fs.P1.inlet)
    m.fs.stream17 = Arc(source=m.fs.P1.outlet, destination=m.fs.MX1.recycle)

    # Link units related to AD section
    m.fs.stream_AD_translator = Arc(
        source=m.fs.AD.liquid_outlet, destination=m.fs.translator_adm1_asm2d.inlet
    )
    m.fs.stream_SP_thickener = Arc(
        source=m.fs.SP2.waste, destination=m.fs.thickener.inlet
    )
    m.fs.stream3adm = Arc(
        source=m.fs.thickener.underflow, destination=m.fs.MX5.thickener
    )
    # m.fs.stream7adm = Arc(source=m.fs.thickener.overflow, destination=m.fs.MX4.recycle2)
    m.fs.stream9adm = Arc(source=m.fs.CL.underflow, destination=m.fs.MX5.clarifier)

    m.fs.stream_translator_dewater = Arc(
        source=m.fs.translator_adm1_asm2d.outlet, destination=m.fs.dewater.inlet
    )
    m.fs.stream_dewater_electroNP = Arc(
        source=m.fs.dewater.overflow, destination=m.fs.electroNP.inlet
    )
    m.fs.stream_electroNP_mixer = Arc(
        source=m.fs.electroNP.treated, destination=m.fs.MX3.recycle1
    )
    # m.fs.stream_dewater_mixer = Arc(
    #     source=m.fs.dewater.overflow, destination=m.fs.MX3.recycle1
    # )

    # m.fs.stream1a = Arc(source=m.fs.FeedWater.outlet, destination=m.fs.CL.inlet)
    m.fs.stream1a = Arc(source=m.fs.FeedWater.outlet, destination=m.fs.MX3.feed_water1)
    m.fs.stream1b = Arc(source=m.fs.MX3.outlet, destination=m.fs.CL.inlet)
    # m.fs.stream1a = Arc(source=m.fs.FeedWater.outlet, destination=m.fs.MX4.feed_water2)
    # m.fs.stream1a = Arc(source=m.fs.FeedWater.outlet, destination=m.fs.MX3.feed_water1)
    # m.fs.stream1b = Arc(source=m.fs.MX3.outlet, destination=m.fs.MX4.feed_water2)
    # m.fs.stream1c = Arc(source=m.fs.MX4.outlet, destination=m.fs.CL.inlet)
    # m.fs.stream1b = Arc(source=m.fs.MX3.outlet, destination=m.fs.CL.inlet)
    m.fs.stream1d = Arc(source=m.fs.CL.effluent, destination=m.fs.MX1.feed_water)

    m.fs.stream10adm = Arc(
        source=m.fs.MX5.outlet, destination=m.fs.translator_asm2d_adm1.inlet
    )
    m.fs.stream_translator_AD = Arc(
        source=m.fs.translator_asm2d_adm1.outlet, destination=m.fs.AD.inlet
    )

    pyo.TransformationFactory("network.expand_arcs").apply_to(m)

    # Oxygen concentration in reactors 3 and 4 is governed by mass transfer
    # Add additional parameter and constraints
    m.fs.R5.KLa = pyo.Var(
        initialize=240,
        units=pyo.units.hour**-1,
        doc="Lumped mass transfer coefficient for oxygen",
    )
    m.fs.R6.KLa = pyo.Var(
        initialize=240,
        units=pyo.units.hour**-1,
        doc="Lumped mass transfer coefficient for oxygen",
    )
    m.fs.R7.KLa = pyo.Var(
        initialize=84,
        units=pyo.units.hour**-1,
        doc="Lumped mass transfer coefficient for oxygen",
    )
    m.fs.S_O_eq = pyo.Param(
        default=8e-3,
        units=pyo.units.kg / pyo.units.m**3,
        mutable=True,
        doc="Dissolved oxygen concentration at equilibrium",
    )

    @m.fs.R5.Constraint(m.fs.time, doc="Mass transfer constraint for R3")
    def mass_transfer_R5(self, t):
        return pyo.units.convert(
            m.fs.R5.injection[t, "Liq", "S_O2"], to_units=pyo.units.kg / pyo.units.hour
        ) == (
            m.fs.R5.KLa
            * m.fs.R5.volume[t]
            * (m.fs.S_O_eq - m.fs.R5.outlet.conc_mass_comp[t, "S_O2"])
        )

    @m.fs.R6.Constraint(m.fs.time, doc="Mass transfer constraint for R4")
    def mass_transfer_R6(self, t):
        return pyo.units.convert(
            m.fs.R6.injection[t, "Liq", "S_O2"], to_units=pyo.units.kg / pyo.units.hour
        ) == (
            m.fs.R6.KLa
            * m.fs.R6.volume[t]
            * (m.fs.S_O_eq - m.fs.R6.outlet.conc_mass_comp[t, "S_O2"])
        )

    @m.fs.R7.Constraint(m.fs.time, doc="Mass transfer constraint for R4")
    def mass_transfer_R7(self, t):
        return pyo.units.convert(
            m.fs.R7.injection[t, "Liq", "S_O2"], to_units=pyo.units.kg / pyo.units.hour
        ) == (
            m.fs.R7.KLa
            * m.fs.R7.volume[t]
            * (m.fs.S_O_eq - m.fs.R7.outlet.conc_mass_comp[t, "S_O2"])
        )

    return m


def set_operating_conditions(m):
    # Mixer settings
    # m.fs.MX3.recycle1_state[0].flow_vol.fix(0 * pyo.units.m**3 / pyo.units.s)
    # m.fs.MX3.recycle1_state[0].temperature.fix(308.15 * pyo.units.K)
    # m.fs.MX3.recycle1_state[0].pressure.fix(1 * pyo.units.atm)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["S_A"].fix(0 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["S_F"].fix(0 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["S_I"].fix(0 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["S_N2"].fix(0 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["S_NH4"].fix(0 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["S_NO3"].fix(0 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["S_O2"].fix(0 * pyo.units.kg / pyo.units.m**3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["S_PO4"].fix(0 * pyo.units.kg / pyo.units.m**3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["S_K"].fix(0 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["S_Mg"].fix(0 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["S_IC"].fix(0 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["X_AUT"].fix(0 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["X_H"].fix(0 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["X_I"].fix(0 * pyo.units.kg / pyo.units.m**3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["X_PAO"].fix(0 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["X_PHA"].fix(0 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["X_PP"].fix(0 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["X_S"].fix(0 * pyo.units.kg / pyo.units.m**3)

    # # Initialize
    # m.fs.MX3.recycle1_state[0].flow_vol.fix(0.0018543 * pyo.units.m**3 / pyo.units.s)
    # m.fs.MX3.recycle1_state[0].temperature.fix(308.15 * pyo.units.K)
    # m.fs.MX3.recycle1_state[0].pressure.fix(1 * pyo.units.atm)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["S_A"].fix(4.1811 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["S_F"].fix(10.964 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["S_I"].fix(0.03 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["S_N2"].fix(1e-9 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["S_NH4"].fix(0.66789 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["S_NO3"].fix(1e-9 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["S_O2"].fix(1e-9 * pyo.units.kg / pyo.units.m**3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["S_PO4"].fix(0.0044715 * pyo.units.kg / pyo.units.m**3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["S_K"].fix(1e-9 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["S_Mg"].fix(1e-9 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["S_IC"].fix(0.39137 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["X_AUT"].fix(1e-9 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["X_H"].fix(1e-9 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["X_I"].fix(0.074180 * pyo.units.kg / pyo.units.m**3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["X_PAO"].fix(1e-9 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["X_PHA"].fix(1e-9 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["X_PP"].fix(1e-9 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["X_S"].fix(0.0015010 * pyo.units.kg / pyo.units.m**3)

    # m.fs.MX3.recycle1_state[0].flow_vol.fix(0.0017727 * pyo.units.m**3 / pyo.units.s)
    # m.fs.MX3.recycle1_state[0].temperature.fix(308.15 * pyo.units.K)
    # m.fs.MX3.recycle1_state[0].pressure.fix(1 * pyo.units.atm)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["S_A"].fix(3.6148 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["S_F"].fix(9.5378 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["S_I"].fix(0.03 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["S_N2"].fix(1e-9 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["S_NH4"].fix(0.53010 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["S_NO3"].fix(1e-9 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["S_O2"].fix(1e-9 * pyo.units.kg / pyo.units.m**3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["S_PO4"].fix(0.0034553 * pyo.units.kg / pyo.units.m**3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["S_K"].fix(1e-9 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["S_Mg"].fix(1e-9 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["S_IC"].fix(0.30863 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["X_AUT"].fix(1e-9 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["X_H"].fix(1e-9 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["X_I"].fix(0.06691 * pyo.units.kg / pyo.units.m**3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["X_PAO"].fix(1e-9 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["X_PHA"].fix(1e-9 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["X_PP"].fix(1e-9 * pyo.units.kg / pyo.units.m ** 3)
    # m.fs.MX3.recycle1_state[0].conc_mass_comp["X_S"].fix(0.0012503 * pyo.units.kg / pyo.units.m**3)

    # Feed Water Conditions
    print(f"DOF before feed: {degrees_of_freedom(m)}")
    m.fs.FeedWater.flow_vol.fix(20648 * pyo.units.m**3 / pyo.units.day)
    m.fs.FeedWater.temperature.fix(308.15 * pyo.units.K)
    m.fs.FeedWater.pressure.fix(1 * pyo.units.atm)
    m.fs.FeedWater.conc_mass_comp[0, "S_O2"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_F"].fix(30 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_A"].fix(20 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_NH4"].fix(16 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_NO3"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_PO4"].fix(3.6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_I"].fix(30 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_N2"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_I"].fix(25 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_S"].fix(125 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_H"].fix(30 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_PAO"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_PP"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_PHA"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_AUT"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_IC"].fix(
        0.07899 * pyo.units.kg / pyo.units.m**3
    )
    m.fs.FeedWater.conc_mass_comp[0, "S_K"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_Mg"].fix(1e-6 * pyo.units.g / pyo.units.m**3)

    # Primary Clarifier
    # TODO: Update primary clarifier once more detailed model available
    m.fs.CL.split_fraction[0, "effluent", "H2O"].fix(0.993)
    m.fs.CL.split_fraction[0, "effluent", "S_A"].fix(0.993)
    m.fs.CL.split_fraction[0, "effluent", "S_F"].fix(0.993)
    m.fs.CL.split_fraction[0, "effluent", "S_I"].fix(0.993)
    m.fs.CL.split_fraction[0, "effluent", "S_N2"].fix(0.993)
    m.fs.CL.split_fraction[0, "effluent", "S_NH4"].fix(0.993)
    m.fs.CL.split_fraction[0, "effluent", "S_NO3"].fix(0.993)
    m.fs.CL.split_fraction[0, "effluent", "S_O2"].fix(0.993)
    m.fs.CL.split_fraction[0, "effluent", "S_PO4"].fix(0.993)
    m.fs.CL.split_fraction[0, "effluent", "S_IC"].fix(0.993)
    m.fs.CL.split_fraction[0, "effluent", "S_K"].fix(0.993)
    m.fs.CL.split_fraction[0, "effluent", "S_Mg"].fix(0.993)
    m.fs.CL.split_fraction[0, "effluent", "X_AUT"].fix(0.5192)
    m.fs.CL.split_fraction[0, "effluent", "X_H"].fix(0.5192)
    m.fs.CL.split_fraction[0, "effluent", "X_I"].fix(0.5192)
    m.fs.CL.split_fraction[0, "effluent", "X_PAO"].fix(0.5192)
    m.fs.CL.split_fraction[0, "effluent", "X_PHA"].fix(0.5192)
    m.fs.CL.split_fraction[0, "effluent", "X_PP"].fix(0.5192)
    m.fs.CL.split_fraction[0, "effluent", "X_S"].fix(0.5192)

    # Reactor sizing
    m.fs.R1.volume.fix(1000 * pyo.units.m**3)
    m.fs.R2.volume.fix(1000 * pyo.units.m**3)
    m.fs.R3.volume.fix(1000 * pyo.units.m**3)
    m.fs.R4.volume.fix(1000 * pyo.units.m**3)
    m.fs.R5.volume.fix(1333 * pyo.units.m**3)
    m.fs.R6.volume.fix(1333 * pyo.units.m**3)
    m.fs.R7.volume.fix(1333 * pyo.units.m**3)

    # Injection rates to Reactions 3, 4 and 5
    for j in m.fs.props_ASM2D.component_list:
        if j != "S_O2":
            # All components except S_O have no injection
            m.fs.R5.injection[:, :, j].fix(0)
            m.fs.R6.injection[:, :, j].fix(0)
            m.fs.R7.injection[:, :, j].fix(0)
    # Then set injections rates for O2
    m.fs.R5.outlet.conc_mass_comp[:, "S_O2"].fix(1.91e-3)
    m.fs.R6.outlet.conc_mass_comp[:, "S_O2"].fix(2.60e-3)
    m.fs.R7.outlet.conc_mass_comp[:, "S_O2"].fix(3.20e-3)

    # Set fraction of outflow from reactor 5 that goes to recycle
    m.fs.SP1.split_fraction[:, "underflow"].fix(0.60)

    # Secondary Clarifier
    # TODO: Update once more detailed model available
    # m.fs.CL1.split_fraction[0, "effluent", "H2O"].fix(0.49986)
    # m.fs.CL1.split_fraction[0, "effluent", "S_A"].fix(0.49986)
    # m.fs.CL1.split_fraction[0, "effluent", "S_F"].fix(0.49986)
    # m.fs.CL1.split_fraction[0, "effluent", "S_I"].fix(0.49986)
    # m.fs.CL1.split_fraction[0, "effluent", "S_N2"].fix(0.49986)
    # m.fs.CL1.split_fraction[0, "effluent", "S_NH4"].fix(0.49986)
    # m.fs.CL1.split_fraction[0, "effluent", "S_NO3"].fix(0.49986)
    # m.fs.CL1.split_fraction[0, "effluent", "S_O2"].fix(0.49986)
    # m.fs.CL1.split_fraction[0, "effluent", "S_PO4"].fix(0.49986)
    # m.fs.CL1.split_fraction[0, "effluent", "S_IC"].fix(0.49986)
    # m.fs.CL1.split_fraction[0, "effluent", "S_K"].fix(0.49986)
    # m.fs.CL1.split_fraction[0, "effluent", "S_Mg"].fix(0.49986)
    # m.fs.CL1.split_fraction[0, "effluent", "X_AUT"].fix(0.022117)
    # m.fs.CL1.split_fraction[0, "effluent", "X_H"].fix(0.021922)
    # m.fs.CL1.split_fraction[0, "effluent", "X_I"].fix(0.021715)
    # m.fs.CL1.split_fraction[0, "effluent", "X_PAO"].fix(0.022)
    # m.fs.CL1.split_fraction[0, "effluent", "X_PHA"].fix(0.02147)
    # m.fs.CL1.split_fraction[0, "effluent", "X_PP"].fix(0.02144)
    # m.fs.CL1.split_fraction[0, "effluent", "X_S"].fix(0.02221)

    m.fs.CL1.split_fraction[0, "effluent", "H2O"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "S_A"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "S_F"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "S_I"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "S_N2"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "S_NH4"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "S_NO3"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "S_O2"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "S_PO4"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "S_IC"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "S_K"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "S_Mg"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "X_AUT"].fix(0.00187)
    m.fs.CL1.split_fraction[0, "effluent", "X_H"].fix(0.00187)
    m.fs.CL1.split_fraction[0, "effluent", "X_I"].fix(0.00187)
    m.fs.CL1.split_fraction[0, "effluent", "X_PAO"].fix(0.00187)
    m.fs.CL1.split_fraction[0, "effluent", "X_PHA"].fix(0.00187)
    m.fs.CL1.split_fraction[0, "effluent", "X_PP"].fix(0.00187)
    m.fs.CL1.split_fraction[0, "effluent", "X_S"].fix(0.00187)

    # Sludge purge separator
    # m.fs.SP2.split_fraction[:, "recycle"].fix(0.97955)
    m.fs.SP2.split_fraction[:, "recycle"].fix(0.985)

    # Outlet pressure from recycle pump
    m.fs.P1.outlet.pressure.fix(101325)

    # AD
    m.fs.AD.volume_liquid.fix(3400)
    m.fs.AD.volume_vapor.fix(300)
    m.fs.AD.liquid_outlet.temperature.fix(308.15)

    # ElectroNP
    m.fs.electroNP.energy_electric_flow_mass.fix(0.044 * pyunits.kWh / pyunits.kg)
    m.fs.electroNP.magnesium_chloride_dosage.fix(0.388)

    # Check degrees of freedom
    print(f"DOF after all: {degrees_of_freedom(m)}")
    assert degrees_of_freedom(m) == 0

    def scale_variables(m):
        for var in m.fs.component_data_objects(pyo.Var, descend_into=True):
            if "flow_vol" in var.name:
                iscale.set_scaling_factor(var, 1e1)
            # if "thickener.properties_in[0.0].flow_vol" in var.name:
            #     iscale.set_scaling_factor(var, 1e3)
            # if "translator_asm2d_adm1.properties_in[0.0].flow_vol" in var.name:
            #     iscale.set_scaling_factor(var, 1e3)
            if "AD.liquid_phase.properties_in[0.0].flow_vol" in var.name:
                iscale.set_scaling_factor(var, 1e3)
            # if "translator_adm1_asm2d.properties_in[0.0].flow_vol" in var.name:
            #     iscale.set_scaling_factor(var, 1e3)
            # if "dewater.properties_in[0.0].flow_vol" in var.name:
            #     iscale.set_scaling_factor(var, 1e3)
            if "electroNP.mixed_state[0.0].flow_vol" in var.name:
                iscale.set_scaling_factor(var, 1e3)
            if "temperature" in var.name:
                iscale.set_scaling_factor(var, 1e-2)
            if "pressure" in var.name:
                iscale.set_scaling_factor(var, 1e-4)
            if "conc_mass_comp" in var.name:
                iscale.set_scaling_factor(var, 1e1)

            # if "conc_mass_comp[S_IN]" in var.name:
            #     iscale.set_scaling_factor(var, 1e1)
            # if "conc_mass_comp[S_IP]" in var.name:
            #     iscale.set_scaling_factor(var, 1e2)

            # if "electroNP.mixed_state[0.0].conc_mass_comp[S_NH4]" in var.name:
            #     iscale.set_scaling_factor(var, 1e2)
            # if "electroNP.mixed_state[0.0].conc_mass_comp[S_PO4]" in var.name:
            #     iscale.set_scaling_factor(var, 1e2)

            # if "conc_mass_comp[S_O2]" in var.name:
            #     iscale.set_scaling_factor(var, 1e4)
            # if "conc_mass_comp[S_F]" in var.name:
            #     iscale.set_scaling_factor(var, 1e2)
            # if "conc_mass_comp[S_A]" in var.name:
            #     iscale.set_scaling_factor(var, 1e2)

            # if "conc_mass_comp[S_NO3]" in var.name:
            #     iscale.set_scaling_factor(var, 1e4)
            # if "conc_mass_comp[S_I]" in var.name:
            #     iscale.set_scaling_factor(var, 1e2)
            # if "conc_mass_comp[S_N2]" in var.name:
            #     iscale.set_scaling_factor(var, 1e4)
            # if "conc_mass_comp[X_I]" in var.name:
            #     iscale.set_scaling_factor(var, 1e0)
            # if "conc_mass_comp[X_S]" in var.name:
            #     iscale.set_scaling_factor(var, 1e1)
            # if "conc_mass_comp[X_H]" in var.name:
            #     iscale.set_scaling_factor(var, 1e0)
            # if "conc_mass_comp[X_PAO]" in var.name:
            #     iscale.set_scaling_factor(var, 1e4)
            # if "conc_mass_comp[X_PP]" in var.name:
            #     iscale.set_scaling_factor(var, 1e4)
            # if "conc_mass_comp[X_PHA]" in var.name:
            #     iscale.set_scaling_factor(var, 1e4)
            # if "conc_mass_comp[X_AUT]" in var.name:
            #     iscale.set_scaling_factor(var, 1e4)
            # if "conc_mass_comp[S_IC]" in var.name:
            #     iscale.set_scaling_factor(var, 1e2)
            # if "conc_mass_comp[S_K]" in var.name:
            #     iscale.set_scaling_factor(var, 1e4)
            # if "conc_mass_comp[S_Mg]" in var.name:
            #     iscale.set_scaling_factor(var, 1e4)

    for unit in ("R1", "R2", "R3", "R4", "R5", "R6", "R7"):
        block = getattr(m.fs, unit)
        iscale.set_scaling_factor(
            block.control_volume.reactions[0.0].rate_expression, 1e3
        )
        iscale.set_scaling_factor(block.cstr_performance_eqn, 1e3)
        iscale.set_scaling_factor(
            block.control_volume.rate_reaction_stoichiometry_constraint, 1e3
        )
        iscale.set_scaling_factor(block.control_volume.material_balances, 1e3)

    # Apply scaling
    scale_variables(m)
    iscale.calculate_scaling_factors(m)


def initialize_system(m):
    # Initialize flowsheet
    # Apply sequential decomposition - 1 iteration should suffice
    seq = SequentialDecomposition()
    # seq.options.select_tear_method = "heuristic"
    seq.options.tear_method = "Direct"
    seq.options.iterLim = 1
    seq.options.tear_set = [m.fs.stream5, m.fs.stream10adm]
    # seq.options.tear_set = [m.fs.stream2, m.fs.stream5, m.fs.stream10adm]

    # seq = SequentialDecomposition()
    # seq.options.select_tear_method = "heuristic"
    # seq.options.tear_method = "Wegstein"
    # seq.options.iterLim = 1

    G = seq.create_graph(m)
    # Uncomment this code to see tear set and initialization order
    order = seq.calculation_order(G)
    print("Initialization Order")
    for o in order:
        print(o[0].name)

    # # Uncomment this code to see tear set and initialization order
    # heuristic_tear_set = seq.tear_set_arcs(G, method="heuristic")
    # order = seq.calculation_order(G)
    # for o in heuristic_tear_set:
    #     print(o.name)
    # for o in order:
    #     print(o[0].name)

    # Initial guesses for flow into first reactor

    tear_guesses1 = {
        "flow_vol": {0: 0.48},
        "conc_mass_comp": {
            (0, "S_A"): 0.01,
            (0, "S_F"): 0.015,
            (0, "S_I"): 0.03,
            (0, "S_N2"): 1e-9,
            (0, "S_NH4"): 0.013,
            (0, "S_NO3"): 1e-9,
            (0, "S_O2"): 0.0016,
            (0, "S_PO4"): 0.0033,
            (0, "S_K"): 1e-9,
            (0, "S_Mg"): 1e-9,
            (0, "S_IC"): 0.087,
            (0, "X_AUT"): 1e-9,
            (0, "X_H"): 0.79,
            (0, "X_I"): 0.21,
            (0, "X_PAO"): 1e-9,
            (0, "X_PHA"): 1e-9,
            (0, "X_PP"): 1e-9,
            (0, "X_S"): 0.0466,
        },
        "temperature": {0: 308.15},
        "pressure": {0: 101325},
    }

    # tear_guesses2 = {
    #     "flow_vol": {0: 1.1861},
    #     "conc_mass_comp": {
    #         (0, "S_A"): 0.013,
    #         (0, "S_F"): 0.00059,
    #         (0, "S_I"): 0.03,
    #         (0, "S_N2"): 1e-9,
    #         (0, "S_NH4"): 0.01,
    #         (0, "S_NO3"): 1e-9,
    #         (0, "S_O2"): 0.0019,
    #         (0, "S_PO4"): 0.0031,
    #         (0, "S_K"): 1e-9,
    #         (0, "S_Mg"): 1e-9,
    #         (0, "S_IC"): 0.092,
    #         (0, "X_AUT"): 1e-9,
    #         (0, "X_H"): 0.81,
    #         (0, "X_I"): 0.21,
    #         (0, "X_PAO"): 1e-9,
    #         (0, "X_PHA"): 1e-9,
    #         (0, "X_PP"): 1e-9,
    #         (0, "X_S"): 0.022,
    #     },
    #     "temperature": {0: 308.15},
    #     "pressure": {0: 101325},
    # }

    # tear_guesses2 = {
    #     "flow_vol": {0: 1.1717},
    #     "conc_mass_comp": {
    #         (0, "S_A"): 0.0339,
    #         (0, "S_F"): 0.00059,
    #         (0, "S_I"): 0.03,
    #         (0, "S_N2"): 1e-9,
    #         (0, "S_NH4"): 0.0079,
    #         (0, "S_NO3"): 1e-9,
    #         (0, "S_O2"): 0.00192,
    #         (0, "S_PO4"): 0.0027,
    #         (0, "S_K"): 1e-9,
    #         (0, "S_Mg"): 1e-9,
    #         (0, "S_IC"): 0.1036,
    #         (0, "X_AUT"): 1e-9,
    #         (0, "X_H"): 1.38,
    #         (0, "X_I"): 0.248,
    #         (0, "X_PAO"): 1e-9,
    #         (0, "X_PHA"): 1e-9,
    #         (0, "X_PP"): 1e-9,
    #         (0, "X_S"): 0.0258,
    #     },
    #     "temperature": {0: 308.15},
    #     "pressure": {0: 101325},
    # }

    # tear_guesses2 = {
    #     "flow_vol": {0: 1.2},
    #     "conc_mass_comp": {
    #         (0, "S_A"): 0.0378,
    #         (0, "S_F"): 0.000717,
    #         (0, "S_I"): 0.03,
    #         (0, "S_N2"): 1e-9,
    #         (0, "S_NH4"): 0.008,
    #         (0, "S_NO3"): 1e-9,
    #         (0, "S_O2"): 0.00192,
    #         (0, "S_PO4"): 0.0027,
    #         (0, "S_K"): 1e-9,
    #         (0, "S_Mg"): 1e-9,
    #         (0, "S_IC"): 0.1036,
    #         (0, "X_AUT"): 1e-9,
    #         (0, "X_H"): 1.38,
    #         (0, "X_I"): 0.248,
    #         (0, "X_PAO"): 1e-9,
    #         (0, "X_PHA"): 1e-9,
    #         (0, "X_PP"): 1e-9,
    #         (0, "X_S"): 0.0258,
    #     },
    #     "temperature": {0: 308.15},
    #     "pressure": {0: 101325},
    # }

    # tear_guesses2 = {
    #     "flow_vol": {0: 1.22},
    #     "conc_mass_comp": {
    #         (0, "S_A"): 0.046,
    #         (0, "S_F"): 0.0002,
    #         (0, "S_I"): 0.03,
    #         (0, "S_N2"): 1e-9,
    #         (0, "S_NH4"): 0.008,
    #         (0, "S_NO3"): 1e-9,
    #         (0, "S_O2"): 0.0019,
    #         (0, "S_PO4"): 0.003,
    #         (0, "S_K"): 1e-9,
    #         (0, "S_Mg"): 1e-9,
    #         (0, "S_IC"): 0.12,
    #         (0, "X_AUT"): 1e-9,
    #         (0, "X_H"): 3.2,
    #         (0, "X_I"): 0.88,
    #         (0, "X_PAO"): 1e-9,
    #         (0, "X_PHA"): 1e-9,
    #         (0, "X_PP"): 1e-9,
    #         (0, "X_S"): 0.04,
    #     },
    #     "temperature": {0: 308.15},
    #     "pressure": {0: 101325},
    # }

    # tear_guesses2 = {
    #     "flow_vol": {0: 1.21},
    #     "conc_mass_comp": {
    #         (0, "S_A"): 0.046,
    #         (0, "S_F"): 0.0002,
    #         (0, "S_I"): 0.03,
    #         (0, "S_N2"): 1e-9,
    #         (0, "S_NH4"): 0.0116,
    #         (0, "S_NO3"): 1e-9,
    #         (0, "S_O2"): 0.00192,
    #         (0, "S_PO4"): 0.0055,
    #         (0, "S_K"): 1e-9,
    #         (0, "S_Mg"): 1e-9,
    #         (0, "S_IC"): 0.12,
    #         (0, "X_AUT"): 1e-9,
    #         (0, "X_H"): 3.21,
    #         (0, "X_I"): 0.88,
    #         (0, "X_PAO"): 1e-9,
    #         (0, "X_PHA"): 1e-9,
    #         (0, "X_PP"): 1e-9,
    #         (0, "X_S"): 0.04,
    #     },
    #     "temperature": {0: 308.15},
    #     "pressure": {0: 101325},
    # }

    # tear_guesses3 = {
    #     "flow_vol": {0: 1.87e-3},
    #     "conc_mass_comp": {
    #         (0, "S_A"): 0.042,
    #         (0, "S_F"): 0.09,
    #         (0, "S_I"): 0.03,
    #         (0, "S_N2"): 1e-9,
    #         (0, "S_NH4"): 0.018,
    #         (0, "S_NO3"): 1e-9,
    #         (0, "S_O2"): 3.15e-4,
    #         (0, "S_PO4"): 0.0035,
    #         (0, "S_K"): 1e-9,
    #         (0, "S_Mg"): 1e-9,
    #         (0, "S_IC"): 0.08,
    #         (0, "X_AUT"): 1e-9,
    #         (0, "X_H"): 8.725,
    #         (0, "X_I"): 2.794,
    #         (0, "X_PAO"): 1e-9,
    #         (0, "X_PHA"): 1e-9,
    #         (0, "X_PP"): 1e-9,
    #         (0, "X_S"): 7.96,
    #     },
    #     "temperature": {0: 308.15},
    #     "pressure": {0: 101325},
    # }

    # tear_guesses3 = {
    #     "flow_vol": {0: 1.87e-3},
    #     "conc_mass_comp": {
    #         (0, "S_A"): 0.042,
    #         (0, "S_F"): 0.09,
    #         (0, "S_I"): 0.03,
    #         (0, "S_N2"): 1e-9,
    #         (0, "S_NH4"): 0.018,
    #         (0, "S_NO3"): 1e-9,
    #         (0, "S_O2"): 3.15e-4,
    #         (0, "S_PO4"): 0.0035,
    #         (0, "S_K"): 1e-9,
    #         (0, "S_Mg"): 1e-9,
    #         (0, "S_IC"): 0.084,
    #         (0, "X_AUT"): 1e-9,
    #         (0, "X_H"): 8.725,
    #         (0, "X_I"): 2.794,
    #         (0, "X_PAO"): 1e-9,
    #         (0, "X_PHA"): 1e-9,
    #         (0, "X_PP"): 1e-9,
    #         (0, "X_S"): 7.78,
    #     },
    #     "temperature": {0: 308.15},
    #     "pressure": {0: 101325},
    # }

    # tear_guesses3 = {
    #     "flow_vol": {0: 2e-4},
    #     "conc_mass_comp": {
    #         (0, "S_A"): 0.047,
    #         (0, "S_F"): 0.103,
    #         (0, "S_I"): 0.03,
    #         (0, "S_N2"): 1e-9,
    #         (0, "S_NH4"): 0.019545,
    #         (0, "S_NO3"): 1e-9,
    #         (0, "S_O2"): 3.15e-4,
    #         (0, "S_PO4"): 0.0035,
    #         (0, "S_K"): 1e-9,
    #         (0, "S_Mg"): 1e-9,
    #         (0, "S_IC"): 0.084,
    #         (0, "X_AUT"): 1e-9,
    #         (0, "X_H"): 8.719,
    #         (0, "X_I"): 2.7975,
    #         (0, "X_PAO"): 1e-9,
    #         (0, "X_PHA"): 1e-9,
    #         (0, "X_PP"): 1e-9,
    #         (0, "X_S"): 7.77,
    #     },
    #     "temperature": {0: 308.15},
    #     "pressure": {0: 101325},
    # }

    # tear_guesses2 = {
    #     "flow_vol": {0: 1.21},
    #     "conc_mass_comp": {
    #         (0, "S_A"): 0.0465,
    #         (0, "S_F"): 0.0002,
    #         (0, "S_I"): 0.03,
    #         (0, "S_N2"): 1e-9,
    #         (0, "S_NH4"): 0.0077,
    #         (0, "S_NO3"): 1e-9,
    #         (0, "S_O2"): 0.0019,
    #         (0, "S_PO4"): 0.003,
    #         (0, "S_K"): 1e-9,
    #         (0, "S_Mg"): 1e-9,
    #         (0, "S_IC"): 0.12,
    #         (0, "X_AUT"): 1e-9,
    #         (0, "X_H"): 2.65,
    #         (0, "X_I"): 0.88,
    #         (0, "X_PAO"): 1e-9,
    #         (0, "X_PHA"): 1e-9,
    #         (0, "X_PP"): 1e-9,
    #         (0, "X_S"): 0.04,
    #     },
    #     "temperature": {0: 308.15},
    #     "pressure": {0: 101325},
    # }
    #
    # tear_guesses3 = {
    #     "flow_vol": {0: 2e-3},
    #     "conc_mass_comp": {
    #         (0, "S_A"): 0.05,
    #         (0, "S_F"): 0.116,
    #         (0, "S_I"): 0.03,
    #         (0, "S_N2"): 1e-9,
    #         (0, "S_NH4"): 0.02,
    #         (0, "S_NO3"): 1e-9,
    #         (0, "S_O2"): 5e-4,
    #         (0, "S_PO4"): 0.0035,
    #         (0, "S_K"): 1e-9,
    #         (0, "S_Mg"): 1e-9,
    #         (0, "S_IC"): 0.091,
    #         (0, "X_AUT"): 1e-9,
    #         (0, "X_H"): 12.918,
    #         (0, "X_I"): 4.53,
    #         (0, "X_PAO"): 1e-9,
    #         (0, "X_PHA"): 1e-9,
    #         (0, "X_PP"): 1e-9,
    #         (0, "X_S"): 7.19,
    #     },
    #     "temperature": {0: 308.15},
    #     "pressure": {0: 101325},
    # }

    # tear_guesses2 = {
    #     "flow_vol": {0: 1.2032},
    #     "conc_mass_comp": {
    #         (0, "S_A"): 0.046474,
    #         (0, "S_F"): 0.00024104,
    #         (0, "S_I"): 0.03,
    #         (0, "S_N2"): 1e-9,
    #         (0, "S_NH4"): 0.0077237,
    #         (0, "S_NO3"): 1e-9,
    #         (0, "S_O2"): 0.00192,
    #         (0, "S_PO4"): 0.0029988,
    #         (0, "S_K"): 1e-9,
    #         (0, "S_Mg"): 1e-9,
    #         (0, "S_IC"): 0.11901,
    #         (0, "X_AUT"): 1e-9,
    #         (0, "X_H"): 3.2121,
    #         (0, "X_I"): 0.87757,
    #         (0, "X_PAO"): 1e-9,
    #         (0, "X_PHA"): 1e-9,
    #         (0, "X_PP"): 1e-9,
    #         (0, "X_S"): 0.038734,
    #     },
    #     "temperature": {0: 308.15},
    #     "pressure": {0: 101325},
    # }
    #
    # tear_guesses3 = {
    #     "flow_vol": {0: 0.0020324},
    #     "conc_mass_comp": {
    #         (0, "S_A"): 0.051250,
    #         (0, "S_F"): 0.11578,
    #         (0, "S_I"): 0.03,
    #         (0, "S_N2"): 1e-9,
    #         (0, "S_NH4"): 0.019672,
    #         (0, "S_NO3"): 1e-9,
    #         (0, "S_O2"): 0.00054398,
    #         (0, "S_PO4"): 0.0034468,
    #         (0, "S_K"): 1e-9,
    #         (0, "S_Mg"): 1e-9,
    #         (0, "S_IC"): 0.090817,
    #         (0, "X_AUT"): 1e-9,
    #         (0, "X_H"): 12.918,
    #         (0, "X_I"): 4.5299,
    #         (0, "X_PAO"): 1e-9,
    #         (0, "X_PHA"): 1e-9,
    #         (0, "X_PP"): 1e-9,
    #         (0, "X_S"): 7.1869,
    #     },
    #     "temperature": {0: 308.15},
    #     "pressure": {0: 101325},
    # }

    # tear_guesses2 = {
    #     "flow_vol": {0: 1.2032},
    #     "conc_mass_comp": {
    #         (0, "S_A"): 0.046474,
    #         (0, "S_F"): 0.00024092,
    #         (0, "S_I"): 0.03,
    #         (0, "S_N2"): 1e-9,
    #         (0, "S_NH4"): 0.0084295,
    #         (0, "S_NO3"): 1e-9,
    #         (0, "S_O2"): 0.00192,
    #         (0, "S_PO4"): 0.0030042,
    #         (0, "S_K"): 1e-9,
    #         (0, "S_Mg"): 1e-9,
    #         (0, "S_IC"): 0.11901,
    #         (0, "X_AUT"): 1e-9,
    #         (0, "X_H"): 3.2121,
    #         (0, "X_I"): 0.87757,
    #         (0, "X_PAO"): 1e-9,
    #         (0, "X_PHA"): 1e-9,
    #         (0, "X_PP"): 1e-9,
    #         (0, "X_S"): 0.038734,
    #     },
    #     "temperature": {0: 308.15},
    #     "pressure": {0: 101325},
    # }
    #
    # tear_guesses3 = {
    #     "flow_vol": {0: 0.0020324},
    #     "conc_mass_comp": {
    #         (0, "S_A"): 0.051,
    #         (0, "S_F"): 0.116,
    #         (0, "S_I"): 0.03,
    #         (0, "S_N2"): 1e-9,
    #         (0, "S_NH4"): 0.0204,
    #         (0, "S_NO3"): 1e-9,
    #         (0, "S_O2"): 0.00054397,
    #         (0, "S_PO4"): 0.0034522,
    #         (0, "S_K"): 1e-9,
    #         (0, "S_Mg"): 1e-9,
    #         (0, "S_IC"): 0.090817,
    #         (0, "X_AUT"): 1e-9,
    #         (0, "X_H"): 12.918,
    #         (0, "X_I"): 4.5299,
    #         (0, "X_PAO"): 1e-9,
    #         (0, "X_PHA"): 1e-9,
    #         (0, "X_PP"): 1e-9,
    #         (0, "X_S"): 7.1869,
    #     },
    #     "temperature": {0: 308.15},
    #     "pressure": {0: 101325},
    # }

    tear_guesses2 = {
        "flow_vol": {0: 1.2032},
        "conc_mass_comp": {
            (0, "S_A"): 0.0465,
            (0, "S_F"): 0.00024,
            (0, "S_I"): 0.03,
            (0, "S_N2"): 1e-9,
            (0, "S_NH4"): 0.00843,
            (0, "S_NO3"): 1e-9,
            (0, "S_O2"): 0.00192,
            (0, "S_PO4"): 0.003,
            (0, "S_K"): 1e-9,
            (0, "S_Mg"): 1e-9,
            (0, "S_IC"): 0.119,
            (0, "X_AUT"): 1e-9,
            (0, "X_H"): 3.21,
            (0, "X_I"): 0.878,
            (0, "X_PAO"): 1e-9,
            (0, "X_PHA"): 1e-9,
            (0, "X_PP"): 1e-9,
            (0, "X_S"): 0.0387,
        },
        "temperature": {0: 308.15},
        "pressure": {0: 101325},
    }

    tear_guesses3 = {
        "flow_vol": {0: 0.0020324},
        "conc_mass_comp": {
            (0, "S_A"): 0.051,
            (0, "S_F"): 0.116,
            (0, "S_I"): 0.03,
            (0, "S_N2"): 1e-9,
            (0, "S_NH4"): 0.0204,
            (0, "S_NO3"): 1e-9,
            (0, "S_O2"): 0.000544,
            (0, "S_PO4"): 0.00345,
            (0, "S_K"): 1e-9,
            (0, "S_Mg"): 1e-9,
            (0, "S_IC"): 0.0908,
            (0, "X_AUT"): 1e-9,
            (0, "X_H"): 12.92,
            (0, "X_I"): 4.53,
            (0, "X_PAO"): 1e-9,
            (0, "X_PHA"): 1e-9,
            (0, "X_PP"): 1e-9,
            (0, "X_S"): 7.19,
        },
        "temperature": {0: 308.15},
        "pressure": {0: 101325},
    }

    # Pass the tear_guess to the SD tool
    # seq.set_guesses_for(m.fs.R1.inlet, tear_guesses1)
    seq.set_guesses_for(m.fs.R3.inlet, tear_guesses2)
    seq.set_guesses_for(m.fs.translator_asm2d_adm1.inlet, tear_guesses3)

    def function(unit):
        unit.initialize(outlvl=idaeslog.INFO, optarg={"bound_push": 1e-2})
        # badly_scaled_vars = list(iscale.badly_scaled_var_generator(unit))
        # if len(badly_scaled_vars) > 0:
        #     automate_rescale_variables(unit)

    # print("Structural issues after setting operating conditions")
    # dt = DiagnosticsToolbox(model=m)
    # dt.report_structural_issues()
    # # dt.report_numerical_issues()
    # # dt.display_variables_with_extreme_jacobians()
    # # dt.display_constraints_with_extreme_jacobians()
    #

    seq.run(m, function)

    # try:
    #     seq.run(m, function)
    # except:
    #     pass

    #
    # print("Numerical issues after initialization")
    # dt.report_numerical_issues()


def solve(m, solver=None):
    if solver is None:
        solver = get_solver()
    results = solver.solve(m, tee=True)
    check_solve(results, checkpoint="closing recycle", logger=_log, fail_flag=True)
    pyo.assert_optimal_termination(results)
    return results
    # results = solver.solve(m, tee=True)


if __name__ == "__main__":
    # This method builds and runs a steady state activated sludge
    # flowsheet.
    m, results = main()

    stream_table = create_stream_table_dataframe(
        {
            "Feed": m.fs.FeedWater.outlet,
            # "Primary clarifier inlet": m.fs.CL.inlet,
            # "R1 inlet": m.fs.R1.inlet,
            # "R1": m.fs.R1.outlet,
            # "R2": m.fs.R2.outlet,
            "R3 inlet": m.fs.R3.inlet,
            # "R3": m.fs.R3.outlet,
            # "R4": m.fs.R4.outlet,
            # "R5": m.fs.R5.outlet,
            # "R6": m.fs.R6.outlet,
            # "R7": m.fs.R7.outlet,
            "ASM outlet": m.fs.SP2.waste,
            # "thickener outlet": m.fs.thickener.underflow,
            # "thickener overflow": m.fs.thickener.overflow,
            "ASM-ADM translator inlet": m.fs.translator_asm2d_adm1.inlet,
            "AD outlet": m.fs.translator_adm1_asm2d.outlet,
            # "AD liquid inlet": m.fs.AD.inlet,
            # "AD liquid outlet": m.fs.AD.liquid_outlet,
            # "AD vapor outlet": m.fs.AD.vapor_outlet,
            # "ADM-ASM translator outlet": m.fs.translator_adm1_asm2d.outlet,
            "dewater outlet": m.fs.dewater.overflow,
            "electroN-P outlet": m.fs.electroNP.treated,
            # "MX3 feed water": m.fs.MX3.feed_water1,
            # "MX3 recycle": m.fs.MX3.recycle1,
            # "MX3 outlet": m.fs.MX3.outlet,
        },
        time_point=0,
    )
    print(stream_table_dataframe_to_string(stream_table))
