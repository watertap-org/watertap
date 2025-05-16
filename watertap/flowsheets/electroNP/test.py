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
"""
Flowsheet example full Water Resource Recovery Facility
(WRRF; a.k.a., wastewater treatment plant) with ASM2d and ADM1 with P extension.

The flowsheet follows the same formulation as benchmark simulation model no.2 (BSM2)
but comprises different specifications for default values than BSM2.
"""

# Some more information about this module
__author__ = "Chenyu Wang, Adam Atia, Alejandro Garciadiego, Marcus Holly, adapted for NH4 Daly Wettermark"

import pyomo.environ as pyo
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
from watertap.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale
from idaes.core.util.tables import (
    create_stream_table_dataframe,
    stream_table_dataframe_to_string,
)
from watertap.unit_models.cstr_injection import CSTR_Injection
from watertap.unit_models.clarifier import Clarifier
from watertap.property_models.unit_specific.anaerobic_digestion.modified_adm1_properties import (
    ModifiedADM1ParameterBlock,
)
from watertap.property_models.unit_specific.anaerobic_digestion.adm1_properties_vapor import (
    ADM1_vaporParameterBlock,
)
from watertap.property_models.unit_specific.anaerobic_digestion.modified_adm1_reactions import (
    ModifiedADM1ReactionParameterBlock,
)
from watertap.property_models.unit_specific.activated_sludge.modified_asm2d_properties import (
    ModifiedASM2dParameterBlock,
)
from watertap.property_models.unit_specific.activated_sludge.modified_asm2d_reactions import (
    ModifiedASM2dReactionParameterBlock,
)
from watertap.unit_models.translators.translator_adm1_asm2d import (
    Translator_ADM1_ASM2D,
)
from idaes.models.unit_models.mixer import MomentumMixingType
from watertap.unit_models.translators.translator_asm2d_adm1 import (
    Translator_ASM2d_ADM1,
)
from watertap.unit_models.anaerobic_digester import AD
from watertap.unit_models.dewatering import (
    DewateringUnit,
    ActivatedSludgeModelType as dewater_type,
)
from watertap.unit_models.thickener import (
    Thickener,
    ActivatedSludgeModelType as thickener_type,
)
from watertap.core.util.initialization import check_solve
from watertap.unit_models.genericNP_ZO import GenericNPZO

from watertap.costing import WaterTAPCosting
from watertap.costing.unit_models.clarifier import (
    cost_circular_clarifier,
    cost_primary_clarifier,
)

from watertap.flowsheets.plot_network import plot_network
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from idaes.core.util.model_diagnostics import DiagnosticsToolbox

import logging

# Set the logging level for the specific logger
idaeslog.getLogger("idaes.core.util.scaling").setLevel(logging.ERROR)

# Alternatively, set the root logger level if you want to suppress all warnings
logging.getLogger().setLevel(logging.ERROR)

import json
import os

this_file_dir = os.path.dirname(os.path.abspath(__file__))

from parameter_sweep.parallel.parallel_manager_factory import create_parallel_manager
from base_values import *

# Set up logger
_log = idaeslog.getLogger(__name__)


def get_stream_dict(m):
    return {
        "FeedWater__MX3": m.fs.FeedWater.outlet,
        "MX3__CL": m.fs.CL.inlet,
        "CL__MX1": m.fs.CL.effluent,
        "CL__MX4": m.fs.CL.underflow,
        "MX4__translator_asm2d_adm1": m.fs.MX4.outlet,
        "translator_asm2d_adm1__AD": m.fs.translator_asm2d_adm1.inlet,
        "AD__translator_adm1_asm2d": m.fs.AD.liquid_outlet,
        "translator_adm1_asm2d__MX3": m.fs.translator_adm1_asm2d.inlet,
        "P1__MX1": m.fs.P1.outlet,
        "MX1__R1": m.fs.R1.inlet,
        "R1__R2": m.fs.R1.outlet,
        "R2__MX2": m.fs.R2.outlet,
        "R3__R4": m.fs.R3.outlet,
        "R4__R5": m.fs.R4.outlet,
        "R5__R6": m.fs.R5.outlet,
        "R6__R7": m.fs.R6.outlet,
        "R7__SP1": m.fs.R7.outlet,
        "SP1__CL2": m.fs.SP1.overflow,
        "SP1__MX2": m.fs.SP1.underflow,
        "MX2__R3": m.fs.MX2.outlet,
        "CL2__SP2": m.fs.CL2.underflow,
        "thickener__MX4": m.fs.thickener.underflow,
        "thickener__MX3": m.fs.thickener.overflow,
        "MX4__translator_asm2d_adm1": m.fs.translator_asm2d_adm1.inlet,
        "translator_asm2d_adm1__AD": m.fs.translator_asm2d_adm1.outlet,
        "AD__translator_adm1_asm2d": m.fs.translator_adm1_asm2d.inlet,
        "translator_adm1_asm2d__dewater": m.fs.translator_adm1_asm2d.outlet,
        "CL2__Treated": m.fs.Treated.inlet,
        "dewater__Sludge": m.fs.Sludge.inlet,
    }


def main(has_genericNP=False):
    m = build_flowsheet(has_genericNP=has_genericNP)
    add_costing(m)

    # deactivate capital cost constraints
    for c in m.fs.component_objects(pyo.Constraint, descend_into=True):
        if "capital_cost" in c.name:
            c.deactivate()

    set_operating_conditions(m)

    for mx in m.fs.mixers:
        mx.pressure_equality_constraints[0.0, 2].deactivate()
    m.fs.MX3.pressure_equality_constraints[0.0, 2].deactivate()
    m.fs.MX3.pressure_equality_constraints[0.0, 3].deactivate()
    print(f"DOF before initialization: {degrees_of_freedom(m)}")

    initialize_system(m, has_genericNP=has_genericNP, outlvl=idaeslog.INFO)
    for mx in m.fs.mixers:
        mx.pressure_equality_constraints[0.0, 2].deactivate()
    m.fs.MX3.pressure_equality_constraints[0.0, 2].deactivate()
    m.fs.MX3.pressure_equality_constraints[0.0, 3].deactivate()
    print(f"DOF after initialization: {degrees_of_freedom(m)}")
    db = DiagnosticsToolbox(m)
    db.report_structural_issues()
    db.display_variables_with_extreme_jacobians()
    db.display_underconstrained_set()
    # db.display_potential_evaluation_errors()
    # db.report_numerical_issues()

    results = solve(m)

    pyo.assert_optimal_termination(results)
    check_solve(
        results,
        checkpoint="re-solve with controls in place",
        logger=_log,
        fail_flag=True,
    )

    results = 0
    return m, results


def build_flowsheet(has_genericNP=False):
    m = pyo.ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.has_genericNP = has_genericNP

    # Properties
    m.fs.props_ASM2D = ModifiedASM2dParameterBlock()
    m.fs.rxn_props_ASM2D = ModifiedASM2dReactionParameterBlock(
        property_package=m.fs.props_ASM2D
    )
    m.fs.props_ADM1 = ModifiedADM1ParameterBlock()
    m.fs.props_vap_ADM1 = ADM1_vaporParameterBlock()
    m.fs.rxn_props_ADM1 = ModifiedADM1ReactionParameterBlock(
        property_package=m.fs.props_ADM1
    )

    # Feed water stream
    m.fs.FeedWater = Feed(property_package=m.fs.props_ASM2D)

    # ====================================================================
    # Primary Clarifier
    m.fs.CL = Clarifier(
        property_package=m.fs.props_ASM2D,
        outlet_list=["underflow", "effluent"],
        split_basis=SplittingType.componentFlow,
    )

    # ======================================================================
    # Activated Sludge Process
    # Mixer for feed water and recycled sludge
    m.fs.MX1 = Mixer(
        property_package=m.fs.props_ASM2D,
        inlet_list=["feed_water", "recycle"],
        momentum_mixing_type=MomentumMixingType.equality,
    )
    # First reactor (anaerobic) - standard CSTR
    m.fs.R1 = CSTR(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    # Second reactor (anaerobic) - standard CSTR
    m.fs.R2 = CSTR(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    # Third reactor (anoxic) - standard CSTR
    m.fs.R3 = CSTR(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    # Fourth reactor (anoxic) - standard CSTR
    m.fs.R4 = CSTR(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    # Fifth reactor (aerobic) - CSTR with injection
    m.fs.R5 = CSTR_Injection(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    # Sixth reactor (aerobic) - CSTR with injection
    m.fs.R6 = CSTR_Injection(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    # Seventh reactor (aerobic) - CSTR with injection
    m.fs.R7 = CSTR_Injection(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    m.fs.SP1 = Separator(
        property_package=m.fs.props_ASM2D, outlet_list=["underflow", "overflow"]
    )
    # Secondary Clarifier
    # TODO: Replace with more detailed model when available
    m.fs.CL2 = Clarifier(
        property_package=m.fs.props_ASM2D,
        outlet_list=["underflow", "effluent"],
        split_basis=SplittingType.componentFlow,
    )
    # Mixing sludge recycle and R5 underflow
    m.fs.MX2 = Mixer(
        property_package=m.fs.props_ASM2D,
        inlet_list=["reactor", "clarifier"],
        momentum_mixing_type=MomentumMixingType.equality,
    )
    # Sludge separator
    m.fs.SP2 = Separator(
        property_package=m.fs.props_ASM2D, outlet_list=["waste", "recycle"]
    )
    # Recycle pressure changer - use a simple isothermal unit for now
    m.fs.P1 = PressureChanger(property_package=m.fs.props_ASM2D)

    # ======================================================================
    # Thickener
    m.fs.thickener = Thickener(
        property_package=m.fs.props_ASM2D,
        activated_sludge_model=thickener_type.modified_ASM2D,
    )
    # Mixing feed and recycle streams from thickener and dewatering unit
    m.fs.MX3 = Mixer(
        property_package=m.fs.props_ASM2D,
        inlet_list=["feed_water", "recycle1", "recycle2"],
        momentum_mixing_type=MomentumMixingType.equality,
    )
    # Mixing sludge from thickener and primary clarifier
    m.fs.MX4 = Mixer(
        property_package=m.fs.props_ASM2D,
        inlet_list=["thickener", "clarifier"],
        momentum_mixing_type=MomentumMixingType.equality,
    )

    # ======================================================================
    # Anaerobic digester section
    # ASM2d-ADM1 translator
    m.fs.translator_asm2d_adm1 = Translator_ASM2d_ADM1(
        inlet_property_package=m.fs.props_ASM2D,
        outlet_property_package=m.fs.props_ADM1,
        inlet_reaction_package=m.fs.rxn_props_ASM2D,
        outlet_reaction_package=m.fs.rxn_props_ADM1,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
        bio_P=False,
    )

    # Anaerobic digester
    m.fs.AD = AD(
        liquid_property_package=m.fs.props_ADM1,
        vapor_property_package=m.fs.props_vap_ADM1,
        reaction_package=m.fs.rxn_props_ADM1,
        has_heat_transfer=True,
        has_pressure_change=False,
    )

    # ADM1-ASM2d translator
    m.fs.translator_adm1_asm2d = Translator_ADM1_ASM2D(
        inlet_property_package=m.fs.props_ADM1,
        outlet_property_package=m.fs.props_ASM2D,
        inlet_reaction_package=m.fs.rxn_props_ADM1,
        outlet_reaction_package=m.fs.rxn_props_ASM2D,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )

    # ======================================================================
    # Dewatering Unit
    m.fs.dewater = DewateringUnit(
        property_package=m.fs.props_ASM2D,
        activated_sludge_model=dewater_type.modified_ASM2D,
    )

    # ======================================================================
    # ElectroN-P
    if has_genericNP is True:
        m.fs.genericNP = GenericNPZO(
            property_package=m.fs.props_ASM2D
        )  # could also add component set, as a dict. or make list dependent on property package

    # ======================================================================
    # Product Blocks
    m.fs.Treated = Product(property_package=m.fs.props_ASM2D)
    m.fs.Sludge = Product(property_package=m.fs.props_ASM2D)
    # Mixers
    m.fs.mixers = (m.fs.MX1, m.fs.MX2, m.fs.MX4)

    # ======================================================================
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
    m.fs.stream11 = Arc(source=m.fs.SP1.overflow, destination=m.fs.CL2.inlet)
    m.fs.stream12 = Arc(source=m.fs.SP1.underflow, destination=m.fs.MX2.clarifier)
    m.fs.stream13 = Arc(source=m.fs.CL2.effluent, destination=m.fs.Treated.inlet)
    m.fs.stream14 = Arc(source=m.fs.CL2.underflow, destination=m.fs.SP2.inlet)
    m.fs.stream15 = Arc(source=m.fs.SP2.recycle, destination=m.fs.P1.inlet)
    m.fs.stream16 = Arc(source=m.fs.P1.outlet, destination=m.fs.MX1.recycle)
    # m.fs.stream17 = Arc(source=m.fs.SP2.waste, destination=m.fs.Sludge.inlet)

    # Link units related to AD section
    m.fs.stream_AD_translator = Arc(
        source=m.fs.AD.liquid_outlet, destination=m.fs.translator_adm1_asm2d.inlet
    )
    m.fs.stream_SP_thickener = Arc(
        source=m.fs.SP2.waste, destination=m.fs.thickener.inlet
    )
    m.fs.stream3adm = Arc(
        source=m.fs.thickener.underflow, destination=m.fs.MX4.thickener
    )
    m.fs.stream7adm = Arc(source=m.fs.thickener.overflow, destination=m.fs.MX3.recycle2)
    m.fs.stream9adm = Arc(source=m.fs.CL.underflow, destination=m.fs.MX4.clarifier)
    m.fs.stream_translator_dewater = Arc(
        source=m.fs.translator_adm1_asm2d.outlet, destination=m.fs.dewater.inlet
    )
    m.fs.stream1a = Arc(source=m.fs.FeedWater.outlet, destination=m.fs.MX3.feed_water)
    m.fs.stream1b = Arc(source=m.fs.MX3.outlet, destination=m.fs.CL.inlet)
    m.fs.stream1c = Arc(source=m.fs.CL.effluent, destination=m.fs.MX1.feed_water)
    m.fs.stream_dewater_sludge = Arc(
        source=m.fs.dewater.underflow, destination=m.fs.Sludge.inlet
    )
    if has_genericNP is True:
        m.fs.stream_dewater_genericNP = Arc(
            source=m.fs.dewater.overflow, destination=m.fs.genericNP.inlet
        )
        m.fs.stream_genericNP_mixer = Arc(
            source=m.fs.genericNP.treated, destination=m.fs.MX3.recycle1
        )
    else:
        m.fs.stream_dewater_mixer = Arc(
            source=m.fs.dewater.overflow, destination=m.fs.MX3.recycle1
        )
    m.fs.stream10adm = Arc(
        source=m.fs.MX4.outlet, destination=m.fs.translator_asm2d_adm1.inlet
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

    m.fs.aerobic_reactors = (m.fs.R5, m.fs.R6, m.fs.R7)
    for R in m.fs.aerobic_reactors:
        iscale.set_scaling_factor(R.KLa, 1e-2)
        iscale.set_scaling_factor(R.hydraulic_retention_time[0], 1e-3)

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
    # eps = 1e-8
    # for R in m.fs.aerobic_reactors:
    #     R.volume[0].setlb(eps)
    #     R.outlet.conc_mass_comp[0, "S_O2"].setlb(eps)

    # Set bounds for AD liquid phase reactions
    # if hasattr(m.fs, 'AD'):
    #     m.fs.AD.liquid_phase.reactions[0.0].S_H.setlb(eps)
    #     m.fs.AD.liquid_phase.reactions[0.0].conc_mol_co2.setlb(eps)
    #     m.fs.AD.liquid_phase.reactions[0.0].conc_mol_hco3.setlb(eps)
    #     m.fs.AD.liquid_phase.reactions[0.0].conc_mol_nh4.setlb(eps)
    #     m.fs.AD.liquid_phase.reactions[0.0].conc_mol_nh3.setlb(eps)
    # m.fs.AD.liquid_phase.properties_in[0.0].flow_vol.setlb(eps)

    # m.fs.AD.liquid_phase.properties_out[0.0].conc_mass_comp["S_IP"].setlb(eps)
    # m.fs.AD.liquid_phase.properties_out[0.0].conc_mass_comp["S_IN"].setlb(eps)
    # for R in [m.fs.R1, m.fs.R3, m.fs.R4, m.fs.R5, m.fs.R6, m.fs.R7]:
    #     for comp in R.control_volume.properties_out[0.0].conc_mass_comp:
    #         R.control_volume.properties_out[0.0].conc_mass_comp[comp].setlb(eps)

    # Feed Water Conditions
    print(f"DOF before feed: {degrees_of_freedom(m)}")
    m.fs.FeedWater.flow_vol.fix(20935.15 * pyo.units.m**3 / pyo.units.day)
    m.fs.FeedWater.temperature.fix(308.15 * pyo.units.K)
    m.fs.FeedWater.pressure.fix(1 * pyo.units.atm)
    m.fs.FeedWater.conc_mass_comp[0, "S_O2"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_F"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_A"].fix(70 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_NH4"].fix(26.6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_NO3"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_PO4"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_I"].fix(57.45 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_N2"].fix(25.19 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_I"].fix(84 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_S"].fix(94.1 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_H"].fix(370 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_PAO"].fix(
        51.5262 * pyo.units.g / pyo.units.m**3
    )
    m.fs.FeedWater.conc_mass_comp[0, "X_PP"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_PHA"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_AUT"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_IC"].fix(5.652 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_K"].fix(374.6925 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_Mg"].fix(20 * pyo.units.g / pyo.units.m**3)

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
    m.fs.R3.volume.fix(1500 * pyo.units.m**3)
    m.fs.R4.volume.fix(1500 * pyo.units.m**3)
    m.fs.R5.volume.fix(3000 * pyo.units.m**3)
    m.fs.R6.volume.fix(3000 * pyo.units.m**3)
    m.fs.R7.volume.fix(3000 * pyo.units.m**3)

    # Injection rates to Reactions 5, 6 and 7
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
    m.fs.CL2.split_fraction[0, "effluent", "H2O"].fix(0.48956)
    m.fs.CL2.split_fraction[0, "effluent", "S_A"].fix(0.48956)
    m.fs.CL2.split_fraction[0, "effluent", "S_F"].fix(0.48956)
    m.fs.CL2.split_fraction[0, "effluent", "S_I"].fix(0.48956)
    m.fs.CL2.split_fraction[0, "effluent", "S_N2"].fix(0.48956)
    m.fs.CL2.split_fraction[0, "effluent", "S_NH4"].fix(0.48956)
    m.fs.CL2.split_fraction[0, "effluent", "S_NO3"].fix(0.48956)
    m.fs.CL2.split_fraction[0, "effluent", "S_O2"].fix(0.48956)
    m.fs.CL2.split_fraction[0, "effluent", "S_PO4"].fix(0.48956)
    m.fs.CL2.split_fraction[0, "effluent", "S_IC"].fix(0.48956)
    m.fs.CL2.split_fraction[0, "effluent", "S_K"].fix(0.48956)
    m.fs.CL2.split_fraction[0, "effluent", "S_Mg"].fix(0.48956)
    m.fs.CL2.split_fraction[0, "effluent", "X_AUT"].fix(0.00187)
    m.fs.CL2.split_fraction[0, "effluent", "X_H"].fix(0.00187)
    m.fs.CL2.split_fraction[0, "effluent", "X_I"].fix(0.00187)
    m.fs.CL2.split_fraction[0, "effluent", "X_PAO"].fix(0.00187)
    m.fs.CL2.split_fraction[0, "effluent", "X_PHA"].fix(0.00187)
    m.fs.CL2.split_fraction[0, "effluent", "X_PP"].fix(0.00187)
    m.fs.CL2.split_fraction[0, "effluent", "X_S"].fix(0.00187)

    m.fs.CL2.surface_area.fix(1500 * pyo.units.m**2)

    # Sludge purge separator
    m.fs.SP2.split_fraction[:, "recycle"].fix(0.985)

    # Outlet pressure from recycle pump
    m.fs.P1.outlet.pressure.fix(101325)

    # AD
    m.fs.AD.volume_liquid.fix(3400)
    m.fs.AD.volume_vapor.fix(300)
    m.fs.AD.liquid_outlet.temperature.fix(308.15)

    # Dewatering Unit - fix either HRT or volume.
    m.fs.dewater.hydraulic_retention_time.fix(1800 * pyo.units.s)

    # Thickener unit
    m.fs.thickener.hydraulic_retention_time.fix(86400 * pyo.units.s)
    m.fs.thickener.diameter.fix(10 * pyo.units.m)

    # genericNP
    if m.fs.has_genericNP is True:
        m.fs.genericNP.magnesium_chloride_dosage.fix(0.388)
        for comp in m.fs.genericNP.energy_electric_flow_mass:
            if comp != "S_NH4":
                m.fs.genericNP.energy_electric_flow_mass[comp].fix(0)
            else:
                m.fs.genericNP.energy_electric_flow_mass[comp].fix(
                    30.0 * pyo.units.kWh / pyo.units.kg
                )

        m.fs.genericNP.P_removal = 0.2
        m.fs.genericNP.NH4_removal = 0.7
        m.fs.genericNP.NO3_removal = 0.0
        m.fs.genericNP.NO2_removal = 0.0
        m.fs.genericNP.frac_mass_H2O_treated[0].fix(0.99)

    def scale_variables(m):
        for var in m.fs.component_data_objects(pyo.Var, descend_into=True):
            if "flow_vol" in var.name:
                iscale.set_scaling_factor(var, 1e-2)
            if "temperature" in var.name:
                iscale.set_scaling_factor(var, 1e-2)
            if "pressure" in var.name:
                iscale.set_scaling_factor(var, 1e-5)
            if "conc_mass_comp" in var.name:
                iscale.set_scaling_factor(var, 1e1)
            if "heat" in var.name:
                iscale.set_scaling_factor(var, 1e7)
            if "rate_reaction_generation" in var.name:
                iscale.set_scaling_factor(var, 1e-5)

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


with open(os.path.join(this_file_dir, "R3_coeffs.json")) as f:
    R3_coeffs = json.load(f)
with open(os.path.join(this_file_dir, "translator_coeffs.json")) as f:
    translator_coeffs = json.load(f)
import random


def get_adjusted_values(base_values, coeffs, NH4_removal, P_removal, EI):
    """Calculate adjusted values based on coefficients and removal rates

    Args:
        base_values (dict): Base case concentrations
        coeffs (dict): Coefficients for adjustments
        NH4_removal (float): Fraction of ammonia removed
        P_removal (float): Fraction of phosphorus removed
        EI (float): Energy intensity value

    Returns:
        dict: Adjusted values
    """
    adj_values = {}
    for key, val in base_values.items():
        coeffs_key = coeffs[key]
        adj_values[key] = (
            val
            + coeffs_key["NH4_removal"] * NH4_removal
            + coeffs_key["P_removal"] * P_removal
            + coeffs_key["EI"] * EI
            + coeffs_key["NH4_removal^2"] * NH4_removal**2
            + coeffs_key["P_removal^2"] * P_removal**2
            + coeffs_key["EI^2"] * EI**2
        )

    return adj_values


def initialize_system(m, has_genericNP=False, outlvl=idaeslog.ERROR):
    # Initialize flowsheet
    # Apply sequential decomposition - 1 iteration should suffice
    seq = SequentialDecomposition()
    seq.options.tear_method = "Direct"
    seq.options.iterLim = 1
    seq.options.tear_set = [m.fs.stream5, m.fs.stream10adm]

    G = seq.create_graph(m)
    # Uncomment this code to see tear set and initialization order
    order = seq.calculation_order(G)
    _log.log(outlvl, "Initialization Order")
    for o in order:
        _log.log(outlvl, o[0].name)

    def generate_tear_guesses(P_removal, NH4_removal):
        """Generate tear guesses based on P and NH4 removal rates

        Args:
            P_removal (float): Fraction of phosphorus removed (no default)
            NH4_removal (float): Fraction of ammonia removed (no default)

        Returns:
            tuple: (tear_guesses, tear_guesses2) dictionaries containing initial guesses
        """

        # Get energy intensity value
        if has_genericNP:
            EI = pyo.value(m.fs.genericNP.energy_electric_flow_mass["S_NH4"])
        else:
            EI = 0

        # # works for intensity 20-50, fraction 0.45-0.7
        #                 # Adjust values based on removal rates
        # adj_values = {}
        # for key, val in base_values_R3.items():
        #     if key == "S_PO4":
        #         adj_values[key] = val * (1 - 0.7*P_removal)
        #     elif key == "S_NH4":
        #         adj_values[key] = val * (1 - 0.2*NH4_removal)
        #     elif key in ["X_PAO", "X_PP"]:
        #         adj_values[key] = val * (1 + 0.1*P_removal)
        #     elif key == "S_NO3":
        #         adj_values[key] = val * (1 - 0.7*NH4_removal)
        #     elif key == "X_AUT":
        #         adj_values[key] = val * (1 - 0.5*NH4_removal)
        #     else:
        #         adj_values[key] = val
        # adj_values2 = {}
        # for key, val in base_values_translator.items():
        #     if key == "S_PO4":
        #         adj_values2[key] = val * (1 - 0.7*P_removal)
        #     elif key == "S_NH4":
        #         adj_values2[key] = val * (1 - 0.2*NH4_removal)
        #     elif key in ["X_PAO", "X_PP"]:
        #         adj_values2[key] = val * (1 + 0.3*P_removal)
        #     elif key == "S_NO3":
        #         adj_values2[key] = val * (1 - 0.7*NH4_removal)
        #     elif key == "X_AUT":
        #         adj_values2[key] = val * (1 - 0.5*NH4_removal)
        #     else:
        #         adj_values2[key] = val

        adj_values = get_adjusted_values(
            base_values_R3, R3_coeffs, NH4_removal, P_removal, EI
        )
        adj_values2 = get_adjusted_values(
            base_values_translator, translator_coeffs, NH4_removal, P_removal, EI
        )

        tear_guesses = {
            "flow_vol": {0: 1.2366},
            "conc_mass_comp": {(0, key): val for key, val in adj_values.items()},
            "temperature": {0: 308.15},
            "pressure": {0: 101325},
        }

        tear_guesses2 = {
            "flow_vol": {0: 0.003},
            "conc_mass_comp": {(0, key): val for key, val in adj_values2.items()},
            "temperature": {0: 308.15},
            "pressure": {0: 101325},
        }

        return tear_guesses, tear_guesses2

    if has_genericNP:
        tear_guesses, tear_guesses2 = generate_tear_guesses(
            P_removal=pyo.value(m.fs.genericNP.P_removal),
            NH4_removal=pyo.value(m.fs.genericNP.NH4_removal),
        )
    else:
        tear_guesses, tear_guesses2 = generate_tear_guesses(P_removal=0, NH4_removal=0)

    # Pass the tear_guess to the SD tool
    seq.set_guesses_for(m.fs.R3.inlet, tear_guesses)
    seq.set_guesses_for(m.fs.translator_asm2d_adm1.inlet, tear_guesses2)

    def function(unit):
        try:
            unit.initialize(outlvl=outlvl, solver="ipopt-watertap")
        except:
            print(f"Failed to initialize {unit.name}")

    # First try initializing all units
    seq.run(m, function)

    # # If MX2 failed, retry with perturbed coefficients
    # if hasattr(m.fs, 'MX2'):
    #     try:
    #         m.fs.MX2.initialize(outlvl=outlvl, solver="ipopt-watertap")
    #         print('MX2 initialized successfully')
    #     except:
    #         print("Retrying MX2 initialization with perturbed R3 coefficients...")
    #         for attempt in range(5):
    #             try:
    #                 # Randomly perturb each coefficient by up to ±2%
    #                 perturbed_R3_coeffs = {
    #                     key: {
    #                         k: v * (1 + (random.random() * 0.04 - 0.02))
    #                         for k, v in val.items()
    #                     }
    #                     for key, val in R3_coeffs.items()
    #                 }

    #                 # Get new adjusted values with perturbed coefficients
    #                 new_adj_values = get_adjusted_values(
    #                     base_values_R3,
    #                     perturbed_R3_coeffs,
    #                     NH4_removal=pyo.value(m.fs.genericNP.NH4_removal),
    #                     P_removal=pyo.value(m.fs.genericNP.P_removal),
    #                     EI=pyo.value(m.fs.genericNP.energy_electric_flow_mass["S_NH4"])
    #                 )

    #                 # Create new tear guesses with perturbed values
    #                 new_tear_guesses = {
    #                     "flow_vol": tear_guesses["flow_vol"],
    #                     "conc_mass_comp": {(0, key): val
    #                                      for key, val in new_adj_values.items()},
    #                     "temperature": tear_guesses["temperature"],
    #                     "pressure": tear_guesses["pressure"]
    #                 }

    #                 seq.set_guesses_for(m.fs.R3.inlet, new_tear_guesses)
    #                 m.fs.MX2.initialize(outlvl=outlvl, solver="ipopt-watertap")
    #                 print(f"Successfully initialized MX2 on attempt {attempt+1}")
    #                 return
    #             except:
    #                 continue
    #         print("All 5 retry attempts for MX2 initialization failed")


def add_costing(m):
    """Add costing block"""
    # m.fs.FeedWater.properties[0].flow_vol.setlb(0)
    # m.fs.Treated.properties[0].flow_vol.setlb(0)

    m.fs.costing = WaterTAPCosting()

    m.fs.costing.base_currency = pyo.units.USD_2020

    # Costing Blocks
    m.fs.R1.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.R2.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.R3.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.R4.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.R5.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.CL.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=cost_primary_clarifier,
    )

    m.fs.CL2.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=cost_circular_clarifier,
    )

    # TODO: Leaving out mixer costs; consider including later

    # process costing and add system level metrics

    if hasattr(m.fs, "genericNP"):
        m.fs.genericNP.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing
        )
        m.fs.costing.genericNP.phosphorus_recovery_value = 0.1
        m.fs.costing.genericNP.ammonia_recovery_value = 0.1

    m.fs.costing.cost_process()
    # m.fs.costing.add_electricity_intensity(m.fs.FeedWater.properties[0].flow_vol)
    # m.fs.costing.add_annual_water_production(m.fs.Treated.properties[0].flow_vol)
    m.fs.costing.add_LCOW(m.fs.FeedWater.properties[0].flow_vol)
    # m.fs.costing.add_specific_energy_consumption(m.fs.FeedWater.properties[0].flow_vol)

    # if hasattr(m.fs, "genericNP"):
    #     m.fs.costing.aggregate_flow_electricity.setlb(0)
    #     # m.fs.costing.aggregate_flow_magnesium_chloride.setlb(0)
    #     # m.fs.costing.aggregate_flow_phosphorus_salt_product.setlb(0)
    #     # m.fs.costing.aggregate_flow_ammonia_product.setlb(0)
    #     for key in m.fs.costing.aggregate_flow_costs:
    #         m.fs.costing.aggregate_flow_costs[key].setlb(0)

    iscale.set_scaling_factor(m.fs.costing.LCOW, 1e2)
    # iscale.set_scaling_factor(m.fs.costing.electricity_intensity, 1e-1)
    iscale.set_scaling_factor(m.fs.costing.total_capital_cost, 1e-7)
    print(
        f"Total operating cost constraint expression: {m.fs.costing.total_operating_cost_constraint.expr}"
    )


def solve(m, solver=None):
    if solver is None:
        solver = get_solver()
    solver.options["tol"] = 1e-1
    solver.options["max_iter"] = 1000
    m.fs.objective = pyo.Objective(expr=m.fs.costing.LCOW)
    results = solver.solve(m, tee=False)
    if m.fs.has_genericNP is True:
        print(f"P removal: {m.fs.genericNP.P_removal.value}")
        print(f"NH4 removal: {m.fs.genericNP.NH4_removal.value}")
    print(f"\nSolver completed with status: {results.solver.status}")
    print(f"Termination condition: {results.solver.termination_condition}\n")
    check_solve(results, checkpoint="closing recycle", logger=_log, fail_flag=True)
    pyo.assert_optimal_termination(results)

    stream_dict = get_stream_dict(m)
    stream_table = create_stream_table_dataframe(stream_dict, time_point=0)
    # save to csv with name f"BSM2_genericNP_stream_table_N{pyo.value(m.fs.genericNP.NH4_removal)}_P{pyo.value(m.fs.genericNP.P_removal)}_E{pyo.value(m.fs.genericNP.energy_electric_flow_mass['S_NH4'])}.csv"
    if m.fs.has_genericNP is False:
        stream_table.to_csv(
            os.path.join(
                this_file_dir,
                "initialization_training/BSM2_genericNP_stream_table_without_genericNP.csv",
            )
        )
    else:
        stream_table.to_csv(
            os.path.join(
                this_file_dir,
                f"initialization_training/BSM2_genericNP_stream_table_N{pyo.value(m.fs.genericNP.NH4_removal):.2f}_P{pyo.value(m.fs.genericNP.P_removal):.2f}_E{pyo.value(m.fs.genericNP.energy_electric_flow_mass['S_NH4']):.2f}.csv",
            )
        )

    return results


run_parameter_sweep = True

from pyomo.environ import Constraint


def print_constraints(model):
    for constraint in model.component_objects(Constraint, active=True):
        print(f"Constraint: {constraint.name}")
        for index in constraint:
            print(f"  Index: {index}, Expression: {constraint[index].expr}")


if __name__ == "__main__" and not run_parameter_sweep:
    # This method builds and runs a steady state activated sludge flowsheet.
    m, results = main(has_genericNP=False)
    stream_dict = get_stream_dict(m)

    if m.fs.has_genericNP is False:
        stream_dict["dewater__MX3"] = m.fs.dewater.overflow
    else:
        stream_dict.update(
            {
                "dewater__MX3": m.fs.dewater.overflow,
                "dewater__genericNP": m.fs.dewater.overflow,
                "genericNP__MX3": m.fs.genericNP.treated,
            }
        )

    stream_table = create_stream_table_dataframe(stream_dict, time_point=0)

    # save to csv with name f"BSM2_genericNP_stream_table_N{pyo.value(m.fs.genericNP.NH4_removal)}_P{pyo.value(m.fs.genericNP.P_removal)}_E{pyo.value(m.fs.genericNP.energy_electric_flow_mass['S_NH4'])}.csv"
    if m.fs.has_genericNP is False:
        stream_table.to_csv(
            os.path.join(
                this_file_dir,
                "initialization_training/BSM2_genericNP_stream_table_without_genericNP.csv",
            )
        )
    else:
        stream_table.to_csv(
            os.path.join(
                this_file_dir,
                f"initialization_training/BSM2_genericNP_stream_table_N{pyo.value(m.fs.genericNP.NH4_removal)}_P{pyo.value(m.fs.genericNP.P_removal)}_E{pyo.value(m.fs.genericNP.energy_electric_flow_mass['S_NH4'])}.csv",
            )
        )

    if m.fs.has_genericNP is False:
        plot_network(
            m,
            stream_table,
            path_to_save=os.path.join(
                this_file_dir,
                "flowsheet_networksBSM2_genericNP_flowsheet_without_genericNP.png",
            ),
            figsize=(12, 10),
        )
    else:
        plot_network(
            m,
            stream_table,
            path_to_save=os.path.join(
                this_file_dir,
                f"flowsheet_networks/BSM2_genericNP_flowsheet_N{pyo.value(m.fs.genericNP.NH4_removal)}_P{pyo.value(m.fs.genericNP.P_removal)}_E{pyo.value(m.fs.genericNP.energy_electric_flow_mass['S_NH4'])}.png",
            ),
            figsize=(12, 10),
        )

    # print the model.fs.costing.electricity_intensity and LCOW
    print(f"Electricity Intensity: {pyo.value(m.fs.costing.electricity_intensity)}")
    print(f"LCOW: {pyo.value(m.fs.costing.LCOW)}")

    print_constraints(m)


from parameter_sweep import (
    LinearSample,
    parameter_sweep,
)


def build_model(**kwargs):
    # return main(has_genericNP=has_genericNP)[0]
    m = build_flowsheet(has_genericNP=True)

    add_costing(m)

    # deactivate capital cost constraints
    for c in m.fs.component_objects(pyo.Constraint, descend_into=True):
        if "capital_cost" in c.name:
            c.deactivate()

    set_operating_conditions(m)
    for mx in m.fs.mixers:
        mx.pressure_equality_constraints[0.0, 2].deactivate()
    m.fs.MX3.pressure_equality_constraints[0.0, 2].deactivate()
    m.fs.MX3.pressure_equality_constraints[0.0, 3].deactivate()
    print(f"DOF before initialization: {degrees_of_freedom(m)}")

    initialize_system(m, has_genericNP=True)
    iscale.calculate_scaling_factors(m)
    for mx in m.fs.mixers:
        mx.pressure_equality_constraints[0.0, 2].deactivate()
    m.fs.MX3.pressure_equality_constraints[0.0, 2].deactivate()
    m.fs.MX3.pressure_equality_constraints[0.0, 3].deactivate()
    print(f"DOF after initialization: {degrees_of_freedom(m)}")

    return m


def build_sweep_params(model, nx=8, **kwargs):
    sweep_params = {}
    sweep_params["NH4 removal %"] = LinearSample(
        model.fs.genericNP.NH4_removal, 0.05, 0.95, nx
    )
    sweep_params["P removal %"] = LinearSample(
        model.fs.genericNP.P_removal, 0.1, 0.9, 5
    )
    sweep_params["NH4 removal energy intensity (kwh/kg)"] = LinearSample(
        model.fs.genericNP.energy_electric_flow_mass["S_NH4"], 1.0, 100.0, nx
    )

    return sweep_params


def build_outputs(model, **kwargs):
    outputs = {}
    try:
        # outputs["Electricity Intensity"] = model.fs.costing.electricity_intensity
        outputs["Treated Water Flow"] = model.fs.Treated.flow_vol[0]
        outputs["Effluent NH4 Concentration"] = model.fs.Treated.conc_mass_comp[
            0, "S_NH4"
        ]

        # MX2 outlet concentrations
        for comp in [
            "S_A",
            "S_F",
            "S_I",
            "S_N2",
            "S_NH4",
            "S_NO3",
            "S_O2",
            "S_PO4",
            "S_K",
            "S_Mg",
            "S_IC",
            "X_AUT",
            "X_H",
            "X_I",
            "X_PAO",
            "X_PHA",
            "X_PP",
            "X_S",
        ]:
            outputs[f"MX2_{comp}"] = model.fs.MX2.outlet.conc_mass_comp[0, comp]

        # MX4 outlet concentrations
        for comp in [
            "S_A",
            "S_F",
            "S_I",
            "S_N2",
            "S_NH4",
            "S_NO3",
            "S_O2",
            "S_PO4",
            "S_K",
            "S_Mg",
            "S_IC",
            "X_AUT",
            "X_H",
            "X_I",
            "X_PAO",
            "X_PHA",
            "X_PP",
            "X_S",
        ]:
            outputs[f"MX4_{comp}"] = model.fs.MX4.outlet.conc_mass_comp[0, comp]

        outputs["P_removal"] = model.fs.genericNP.P_removal
        outputs["NH4_removal"] = model.fs.genericNP.NH4_removal
        outputs["EI"] = model.fs.genericNP.energy_electric_flow_mass["S_NH4"]
        outputs["LCOW"] = model.fs.costing.LCOW
        outputs["Product_Value"] = model.fs.costing.genericNP.ammonia_recovery_value
    except:
        print("Unable to solve")
        # Return empty values if solve fails
        for key in outputs:
            outputs[key] = None
    return outputs


def reinitialize_system(model):
    set_operating_conditions(model)
    initialize_system(model, has_genericNP=True)


# Create a parallel manager with the desired backend
parallel_manager = create_parallel_manager(
    parallel_manager_class=None,  # Automatically selects based on backend
    number_of_subprocesses=4,  # Number of subprocesses for shared memory backends
    parallel_back_end="MPI",  # Choose your backend here
)
# run with mpirun -np 10 python BSM2_electroNH4.py


def run_analysis(case_num=11, interpolate_nan_outputs=True, output_filename=None):
    if output_filename is None:
        output_filename = (
            this_file_dir + f"/initialization_training/sensitivity_{case_num}"
        )

    global_results = parameter_sweep(
        build_model,
        build_sweep_params,
        build_outputs,
        csv_results_file_name=f"{output_filename}.csv",
        # h5_results_file_name=f"{output_filename}.h5",
        optimize_function=solve,
        reinitialize_function=reinitialize_system,
        interpolate_nan_outputs=interpolate_nan_outputs,
    )

    return global_results


if run_parameter_sweep:
    results = run_analysis()
    # create dataframe of results
    df_results = pd.DataFrame()
    df_results["NH4 removal %"] = results[1]["sweep_params"]["NH4 removal %"]["value"]
    df_results["NH4 removal energy intensity (kwh/kg)"] = results[1]["sweep_params"][
        "NH4 removal energy intensity (kwh/kg)"
    ]["value"]
    df_results["P removal %"] = results[1]["sweep_params"]["P removal %"]["value"]
    # df_results["Electricity Intensity"] = results[1]["outputs"][
    #     "Electricity Intensity"
    # ]["value"]
    df_results["LCOW"] = results[1]["outputs"]["LCOW"]["value"]
    df_results["N value"] = results[1]["outputs"]["Product_Value"]["value"]

    df_results_fixed_p = df_results[
        df_results["P removal %"] == list(df_results["P removal %"])[0]
    ]
    # # Create pivot table for electricity intensity
    # pivot_df_ei = df_results_fixed_p.pivot(
    #     index="NH4 removal %",
    #     columns="NH4 removal energy intensity (kwh/kg)",
    #     values="Electricity Intensity",
    # )
    # pivot_df_ei = pivot_df_ei.round(3)
    # pivot_df_ei.index = pivot_df_ei.index.round(3)
    # pivot_df_ei.columns = pivot_df_ei.columns.round(3)

    # Create pivot table for LCOW
    pivot_df_lcow = df_results_fixed_p.pivot(
        index="NH4 removal %",
        columns="NH4 removal energy intensity (kwh/kg)",
        values="LCOW",
    )
    pivot_df_lcow = pivot_df_lcow.round(3)
    pivot_df_lcow.index = pivot_df_lcow.index.round(3)
    pivot_df_lcow.columns = pivot_df_lcow.columns.round(3)

    plt.show()

    # Plot LCOW heatmap
    plt.figure(figsize=(6, 6))
    heatmap2 = sns.heatmap(
        pivot_df_lcow,
        annot=True,
        cmap="YlOrBr",
        cbar_kws={"label": "LCOT ($/m3) including NH4 recovery value"},
    )
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel("NH4 Removal Energy Intensity (kWh/kg)", fontsize=12)
    plt.ylabel("NH4 Removal %", fontsize=12)
    plt.savefig(
        f"sensitivity_heatmap_lcow_nh4_val_0.1.png", dpi=300, bbox_inches="tight"
    )
    plt.show()
