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

[1] J. Alex, L. Benedetti, J. Copp, K.V. Gernaey, U. Jeppsson, I. Nopens, M.N. Pons,
C.Rosen, J.P. Steyer and P. Vanrolleghem, "Benchmark Simulation Model no. 2 (BSM2)", 2018

[2] J. Alex, L. Benedetti, J. Copp, K.V. Gernaey, U. Jeppsson, I. Nopens, M.N. Pons,
J.P. Steyer and P. Vanrolleghem, "Benchmark Simulation Model no. 1 (BSM1)", 2018
"""

# Some more information about this module
__author__ = "Chenyu Wang, Adam Atia, Alejandro Garciadiego, Marcus Holly"

import pyomo.environ as pyo
from pyomo.network import Arc, SequentialDecomposition

from idaes.core import (
    FlowsheetBlock,
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
from watertap.costing.unit_models.genericNP import cost_genericNP
from watertap.costing.unit_models.clarifier import (
    cost_primary_clarifier,
    cost_circular_clarifier,
)

# from watertap.tools.plot_network import plot_network
from watertap.tools.dash_network import (
    create_dash_app,
    find_available_port,
    run_dash_app,
)

from idaes.core import (
    FlowsheetBlock,
    UnitModelCostingBlock,
)

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from parameter_sweep import LinearSample, parameter_sweep
import time

# Set up logger
_log = idaeslog.getLogger(__name__)

this_file_dir = os.path.dirname(os.path.abspath(__file__))


def safe_value(x):
    try:
        return pyo.value(x)
    except Exception:
        return None


def print_constraint_violations(m):
    print("\n--- Constraints with High Violation (after failed solve) ---")
    for c in m.component_data_objects(pyo.Constraint, active=True):
        if c.is_indexed():
            for idx in c:
                con = c[idx]
                try:
                    body = safe_value(con.body)
                    lb = safe_value(con.lower) if con.has_lb() else None
                    ub = safe_value(con.upper) if con.has_ub() else None
                    if lb is not None and body is not None and body < lb - 1e-4:
                        print(f"{c.name}{idx} violated: body={body}, lb={lb}")
                    if ub is not None and body is not None and body > ub + 1e-4:
                        print(f"{c.name}{idx} violated: body={body}, ub={ub}")
                except Exception as e2:
                    print(f"{c.name}{idx}: ERROR ({e2})")
        else:
            con = c
            try:
                body = safe_value(con.body)
                lb = safe_value(con.lower) if con.has_lb() else None
                ub = safe_value(con.upper) if con.has_ub() else None
                if lb is not None and body is not None and body < lb - 1e-4:
                    print(f"{c.name} violated: body={body}, lb={lb}")
                if ub is not None and body is not None and body > ub + 1e-4:
                    print(f"{c.name} violated: body={body}, ub={ub}")
            except Exception as e2:
                print(f"{c.name}: ERROR ({e2})")


def print_costing_scaling(m):
    print("\n--- Costing Variable Scaling Factors ---")
    for v in m.fs.costing.component_data_objects(pyo.Var, descend_into=True):
        sf = iscale.get_scaling_factor(v)
        print(f"{v.name}: scaling_factor={sf}")
    print("\n--- Costing Constraint Scaling Factors ---")
    for c in m.fs.costing.component_data_objects(pyo.Constraint, descend_into=True):
        sf = iscale.get_scaling_factor(c)
        print(f"{c.name}: scaling_factor={sf}")


# DEBUG: Try initializing costing block variables
# This is a simple approach: set all uninitialized costing variables to a small value
# (or 1.0 if appropriate) if they are not fixed and have no value.
def initialize_costing_block(m):
    for v in m.fs.costing.component_data_objects(pyo.Var, descend_into=True):
        if (not v.fixed) and (v.value is None):
            try:
                v.set_value(1.0)
            except Exception:
                pass


def main(
    has_genericNP=False,
    plot_network_before_solve=False,
    basis="mass",
    p_removal=0.95,
    n_to_p_ratio=0.3,
    energy_intensity=0.044,
    mgcl2_dosage=0.388,
    water_recovery=0.99,
    apply_costing=True,
):
    m = build_flowsheet(has_genericNP=has_genericNP, basis=basis)

    # Add costing BEFORE setting operating conditions (if requested)
    if apply_costing:
        add_costing(m)
        # Deactivate capital cost constraints during initialization
        for c in m.fs.component_objects(pyo.Constraint, descend_into=True):
            if "capital_cost" in c.name:
                c.deactivate()

    set_operating_conditions(
        m,
        p_removal=p_removal,
        n_to_p_ratio=n_to_p_ratio,
        energy_intensity=energy_intensity,
        mgcl2_dosage=mgcl2_dosage,
        water_recovery=water_recovery,
    )

    # Now initialize the process (with costing present if apply_costing=True)
    if plot_network_before_solve:
        # Create empty stream table for initial network plot
        stream_dict = {
            "Feed": m.fs.FeedWater.outlet,
            "Treated": m.fs.Treated.inlet,
            "Sludge": m.fs.Sludge.inlet,
        }
        if has_genericNP:
            stream_dict.update(
                {
                    "genericNP_treated": m.fs.genericNP.treated,
                    "genericNP_byproduct": m.fs.genericNP.byproduct,
                }
            )
        stream_table = create_stream_table_dataframe(stream_dict, time_point=0)

        # Create and run Dash app for initial network visualization
        app = create_dash_app(m, stream_table=stream_table, show_labels=True)
        if app is not None:
            port = find_available_port()
            print(f"Starting initial network visualization on port {port}")
            thread = run_dash_app(app, port)
            if thread is not None:
                print(
                    f"Initial network visualization available at http://127.0.0.1:{port}"
                )
                time.sleep(3)  # Give time for the visualization to load
        else:
            print("Failed to create initial network visualization")
    initialize_system(m, has_genericNP=has_genericNP)
    for mx in m.fs.mixers:
        mx.pressure_equality_constraints[0.0, 2].deactivate()
    m.fs.MX3.pressure_equality_constraints[0.0, 2].deactivate()
    m.fs.MX3.pressure_equality_constraints[0.0, 3].deactivate()
    print(f"DOF after initialization: {degrees_of_freedom(m)}")

    # Solve the process (with costing if apply_costing=True)
    if apply_costing:
        initialize_costing_block(m)
        print_costing_scaling(m)

        def safe_cost_value(x):
            try:
                return pyo.value(x)
            except Exception:
                return "UNINITIALIZED"

        # Costing debug output
        for unit in [
            m.fs.genericNP,
            m.fs.R1,
            m.fs.R2,
            m.fs.R3,
            m.fs.R4,
            m.fs.R5,
            m.fs.CL,
            m.fs.CL2,
        ]:
            if hasattr(unit, "costing"):
                print(f"\n--- {unit.name} costing variables ---")
                for v in unit.costing.component_data_objects(
                    pyo.Var, descend_into=True
                ):
                    print(f"{v.name} = {v.value}")

        print("\n--- Costing Constraint Residuals ---")
        for c in m.fs.costing.component_objects(pyo.Constraint, descend_into=True):
            for idx in c:
                try:
                    body = safe_cost_value(c[idx].body)
                    lb = safe_cost_value(c[idx].lower) if c[idx].has_lb() else None
                    ub = safe_cost_value(c[idx].upper) if c[idx].has_ub() else None
                    res = (
                        (body - lb)
                        if lb not in [None, "UNINITIALIZED"] and body != "UNINITIALIZED"
                        else "UNINITIALIZED"
                    )
                    print(f"{c.name}{idx}: {res} (body={body}, lb={lb}, ub={ub})")
                except Exception as e:
                    print(f"{c.name}{idx}: ERROR ({e})")

        print("\n--- Costing Variables (NaN or extreme values) ---")
        for v in m.fs.costing.component_data_objects(pyo.Var, descend_into=True):
            if v.value is not None and (abs(v.value) > 1e6 or abs(v.value) < 1e-12):
                print(f"{v.name} = {v.value}")

        print("\n--- Costing Constraint Expressions and Residuals ---")
        for c in m.fs.costing.component_objects(pyo.Constraint, descend_into=True):
            for idx in c:
                try:
                    expr = c[idx].expr
                    body = safe_cost_value(c[idx].body)
                    lb = safe_cost_value(c[idx].lower) if c[idx].has_lb() else None
                    ub = safe_cost_value(c[idx].upper) if c[idx].has_ub() else None
                    res = (
                        (body - lb)
                        if lb not in [None, "UNINITIALIZED"] and body != "UNINITIALIZED"
                        else "UNINITIALIZED"
                    )
                    print(
                        f"{c.name}{idx}: expr={expr}, body={body}, lb={lb}, ub={ub}, residual={res}"
                    )
                except Exception as e:
                    print(f"{c.name}{idx}: ERROR ({e})")

        print("\n--- Key Costing Variables ---")
        costing_vars = [
            m.fs.genericNP.electricity[0],
            m.fs.genericNP.MgCl2_flowrate[0],
            m.fs.genericNP.byproduct.flow_vol[0],
            m.fs.genericNP.byproduct.conc_mass_comp[0, "S_PO4"],
            m.fs.genericNP.byproduct.conc_mass_comp[0, "S_NH4"],
            m.fs.genericNP.mixed_state[0].flow_vol,
        ]
        for v in costing_vars:
            try:
                print(f"{v.name} = {pyo.value(v)}")
            except Exception as e:
                print(f"{v.name} = ERROR ({e})")

        print("\n--- All Fixed Variables ---")
        for v in m.component_data_objects(pyo.Var, descend_into=True):
            if v.fixed:
                try:
                    print(f"{v.name} (fixed) = {pyo.value(v)}")
                except Exception as e:
                    print(f"{v.name} (fixed) = ERROR ({e})")

        print("\n--- Uninitialized Costing Variables ---")
        for v in m.fs.costing.component_data_objects(pyo.Var, descend_into=True):
            if v.value is None:
                print(f"{v.name} is uninitialized (None)")

        print("--- Costing Constraints ---")
        for c in m.fs.costing.component_objects(pyo.Constraint, descend_into=True):
            print(c.name)
        print("\n--- Costing Variables (unfixed/uninitialized) ---")
        for v in m.fs.costing.component_data_objects(pyo.Var, descend_into=True):
            if not v.fixed and v.value is None:
                print(f"{v.name} is unfixed and uninitialized")
        print("\n--- Degrees of Freedom ---")
        print(degrees_of_freedom(m))

    # Final solve with costing active
    try:
        results = solve(m, max_iter=5000)
        pyo.assert_optimal_termination(results)
        check_solve(
            results,
            checkpoint=(
                "final solve with costing"
                if apply_costing
                else "final solve without costing"
            ),
            logger=_log,
            fail_flag=True,
        )
        print_constraint_violations(m)
    except Exception:
        print_constraint_violations(m)
        raise

    return m, results


def build_flowsheet(has_genericNP=False, basis="mass"):
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
    # GenericNP
    if has_genericNP is True:
        m.fs.genericNP = GenericNPZO(property_package=m.fs.props_ASM2D, basis=basis)

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
    # KLa for R5 and R6 taken from [1], and KLa for R7 taken from [2]
    m.fs.R5.KLa = pyo.Var(
        initialize=240 / 24,
        units=pyo.units.hour**-1,
        doc="Lumped mass transfer coefficient for oxygen",
    )
    m.fs.R6.KLa = pyo.Var(
        initialize=240 / 24,
        units=pyo.units.hour**-1,
        doc="Lumped mass transfer coefficient for oxygen",
    )
    m.fs.R7.KLa = pyo.Var(
        initialize=84 / 24,
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


def set_operating_conditions(
    m,
    p_removal=0.95,
    n_to_p_ratio=0.3,
    energy_intensity=0.044,
    mgcl2_dosage=0.388,
    water_recovery=0.99,
):
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
        nh4_removal = p_removal * n_to_p_ratio
        m.fs.genericNP.removal_factors["S_PO4"].set_value(p_removal)
        m.fs.genericNP.removal_factors["S_NH4"].set_value(nh4_removal)
        if "S_NO3" in m.fs.genericNP.removal_factors:
            m.fs.genericNP.removal_factors["S_NO3"].set_value(0.0)
        if "S_NO2" in m.fs.genericNP.removal_factors:
            m.fs.genericNP.removal_factors["S_NO2"].set_value(0.0)

        if m.fs.genericNP.config.basis == "mass":
            energy_units = pyo.units.kWh / pyo.units.kg
            mg_units = pyo.units.kg / pyo.units.kg
        else:
            energy_units = pyo.units.kWh / pyo.units.mol
            mg_units = pyo.units.kg / pyo.units.mol

        for comp in m.fs.props_ASM2D.component_list:
            if comp in ("S_PO4", "S_NH4"):
                m.fs.genericNP.energy_electric_flow[comp].set_value(
                    energy_intensity * energy_units
                )
            else:
                m.fs.genericNP.energy_electric_flow[comp].set_value(0 * energy_units)

        m.fs.genericNP.magnesium_chloride_dosage.fix(mgcl2_dosage * mg_units)
        m.fs.genericNP.frac_mass_H2O_treated[0].fix(water_recovery)

    def scale_variables(m):
        for var in m.fs.component_data_objects(pyo.Var, descend_into=True):
            if "flow_vol" in var.name:
                iscale.set_scaling_factor(var, 1e2)
            if "temperature" in var.name:
                iscale.set_scaling_factor(var, 1e-2)
            if "pressure" in var.name:
                iscale.set_scaling_factor(var, 1e-5)
            if "conc_mass_comp" in var.name:
                iscale.set_scaling_factor(var, 1e1)

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

    # Initialize and scale genericNP unit
    if m.fs.has_genericNP:
        if m.fs.genericNP.electricity[0].value is None:
            m.fs.genericNP.electricity[0].set_value(1.0)
        if m.fs.genericNP.MgCl2_flowrate[0].value is None:
            m.fs.genericNP.MgCl2_flowrate[0].set_value(1.0)

        # Set scaling factors for genericNP
        iscale.set_scaling_factor(m.fs.genericNP.electricity, 1e-1)
        iscale.set_scaling_factor(m.fs.genericNP.MgCl2_flowrate, 1e0)
        iscale.set_scaling_factor(m.fs.genericNP.magnesium_chloride_dosage, 1e0)
        iscale.set_scaling_factor(m.fs.genericNP.frac_mass_H2O_treated, 1e0)

    # Apply scaling
    scale_variables(m)
    iscale.calculate_scaling_factors(m)


def initialize_system(m, has_genericNP=False):
    # Initialize flowsheet
    # Apply sequential decomposition - 1 iteration should suffice
    seq = SequentialDecomposition()
    seq.options.tear_method = "Direct"
    seq.options.iterLim = 1
    seq.options.tear_set = [m.fs.stream5, m.fs.stream10adm]

    G = seq.create_graph(m)
    # Uncomment this code to see tear set and initialization order
    order = seq.calculation_order(G)
    print("Initialization Order")
    for o in order:
        print(o[0].name)

    if has_genericNP:
        # P_removal = 0.65 - 0.95
        tear_guesses = {
            "flow_vol": {0: 1.2366},
            "conc_mass_comp": {
                (0, "S_A"): 0.0006,
                (0, "S_F"): 0.0004,
                (0, "S_I"): 0.057,
                (0, "S_N2"): 0.04,
                (0, "S_NH4"): 0.006,
                (0, "S_NO3"): 0.002,
                (0, "S_O2"): 0.0019,
                (0, "S_PO4"): 0.09,
                (0, "S_K"): 0.37,
                (0, "S_Mg"): 0.020,
                (0, "S_IC"): 0.13,
                (0, "X_AUT"): 0.085,
                (0, "X_H"): 3.5,
                (0, "X_I"): 3.1,
                (0, "X_PAO"): 3.4,
                (0, "X_PHA"): 0.087,
                (0, "X_PP"): 1.1,
                (0, "X_S"): 0.057,
            },
            "temperature": {0: 308.15},
            "pressure": {0: 101325},
        }

        tear_guesses2 = {
            "flow_vol": {0: 0.003},
            "conc_mass_comp": {
                (0, "S_A"): 0.1,
                (0, "S_F"): 0.15,
                (0, "S_I"): 0.057,
                (0, "S_N2"): 0.034,
                (0, "S_NH4"): 0.025,
                (0, "S_NO3"): 0.0015,
                (0, "S_O2"): 0.0013,
                (0, "S_PO4"): 0.1,
                (0, "S_K"): 0.38,
                (0, "S_Mg"): 0.024,
                (0, "S_IC"): 0.074,
                (0, "X_AUT"): 0.21,
                (0, "X_H"): 23,
                (0, "X_I"): 11,
                (0, "X_PAO"): 10.5,
                (0, "X_PHA"): 0.006,
                (0, "X_PP"): 2.7,
                (0, "X_S"): 3.9,
            },
            "temperature": {0: 308.15},
            "pressure": {0: 101325},
        }

    else:
        tear_guesses = {
            "flow_vol": {0: 1.2368},
            "conc_mass_comp": {
                (0, "S_A"): 0.0006,
                (0, "S_F"): 0.0004,
                (0, "S_I"): 0.057,
                (0, "S_N2"): 0.047,
                (0, "S_NH4"): 0.0075,
                (0, "S_NO3"): 0.003,
                (0, "S_O2"): 0.0019,
                (0, "S_PO4"): 0.73,
                (0, "S_K"): 0.37,
                (0, "S_Mg"): 0.020,
                (0, "S_IC"): 0.13,
                (0, "X_AUT"): 0.11,
                (0, "X_H"): 3.5,
                (0, "X_I"): 3.2,
                (0, "X_PAO"): 3.2,
                (0, "X_PHA"): 0.084,
                (0, "X_PP"): 1.07,
                (0, "X_S"): 0.057,
            },
            "temperature": {0: 308.15},
            "pressure": {0: 101325},
        }

        tear_guesses2 = {
            "flow_vol": {0: 0.003},
            "conc_mass_comp": {
                (0, "S_A"): 0.097,
                (0, "S_F"): 0.15,
                (0, "S_I"): 0.057,
                (0, "S_N2"): 0.036,
                (0, "S_NH4"): 0.03,
                (0, "S_NO3"): 0.002,
                (0, "S_O2"): 0.0013,
                (0, "S_PO4"): 0.74,
                (0, "S_K"): 0.38,
                (0, "S_Mg"): 0.024,
                (0, "S_IC"): 0.075,
                (0, "X_AUT"): 0.28,
                (0, "X_H"): 23.4,
                (0, "X_I"): 11.4,
                (0, "X_PAO"): 10.1,
                (0, "X_PHA"): 0.0044,
                (0, "X_PP"): 2.7,
                (0, "X_S"): 3.9,
            },
            "temperature": {0: 308.15},
            "pressure": {0: 101325},
        }

    # Pass the tear_guess to the SD tool
    seq.set_guesses_for(m.fs.R3.inlet, tear_guesses)
    seq.set_guesses_for(m.fs.translator_asm2d_adm1.inlet, tear_guesses2)

    def function(unit):
        unit.initialize(outlvl=idaeslog.INFO, solver="ipopt-watertap")

    seq.run(m, function)


def solve(m, solver=None, max_iter=3000):
    if solver is None:
        solver = get_solver()
    # Temporarily increase iteration limit for genericNP convergence
    solver.options["max_iter"] = max_iter
    results = solver.solve(m, tee=True)
    check_solve(results, checkpoint="closing recycle", logger=_log, fail_flag=True)
    pyo.assert_optimal_termination(results)
    return results


def add_costing(m):
    """Add costing block"""
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

    if hasattr(m.fs, "genericNP"):
        m.fs.genericNP.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing
        )
        m.fs.costing.genericNP.phosphorus_recovery_value = 0.1
        m.fs.costing.genericNP.ammonia_recovery_value = 0.1

    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(m.fs.FeedWater.properties[0].flow_vol)

    iscale.set_scaling_factor(m.fs.costing.LCOW, 1e2)
    iscale.set_scaling_factor(m.fs.costing.total_capital_cost, 1e-7)


def build_sweep_params(model, nx=2, **kwargs):
    sweep_params = {}
    sweep_params["NH4 removal %"] = LinearSample(
        model.fs.genericNP.removal_factors["S_NH4"], 0.05, 0.95, nx
    )
    sweep_params["P removal %"] = LinearSample(
        model.fs.genericNP.removal_factors["S_PO4"], 0.1, 0.9, 2
    )
    sweep_params["NH4 removal energy intensity (kwh/kg)"] = LinearSample(
        model.fs.genericNP.energy_electric_flow["S_NH4"], 1.0, 100.0, nx
    )
    return sweep_params


def build_model(**kwargs):
    basis = kwargs.get("basis", "mass")
    p_removal = kwargs.get("p_removal", 0.95)
    n_to_p_ratio = kwargs.get("n_to_p_ratio", 0.3)
    energy_intensity = kwargs.get("energy_intensity", 0.044)
    mgcl2_dosage = kwargs.get("mgcl2_dosage", 0.388)
    water_recovery = kwargs.get("water_recovery", 0.99)

    # return main(has_genericNP=has_genericNP)[0]
    m = build_flowsheet(has_genericNP=True, basis=basis)

    add_costing(m)

    # deactivate capital cost constraints
    for c in m.fs.component_objects(pyo.Constraint, descend_into=True):
        if "capital_cost" in c.name:
            c.deactivate()

    set_operating_conditions(
        m,
        p_removal=p_removal,
        n_to_p_ratio=n_to_p_ratio,
        energy_intensity=energy_intensity,
        mgcl2_dosage=mgcl2_dosage,
        water_recovery=water_recovery,
    )
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


def build_outputs(model, **kwargs):
    outputs = {}
    try:
        outputs["Treated Water Flow"] = model.fs.Treated.flow_vol[0]
        outputs["Effluent NH4 Concentration"] = model.fs.Treated.conc_mass_comp[
            0, "S_NH4"
        ]
        outputs["Effluent PO4 Concentration"] = model.fs.Treated.conc_mass_comp[
            0, "S_PO4"
        ]
        outputs["P_removal"] = model.fs.genericNP.removal_factors["S_PO4"]
        outputs["NH4_removal"] = model.fs.genericNP.removal_factors["S_NH4"]
        outputs["EI"] = model.fs.genericNP.energy_electric_flow["S_NH4"]
        outputs["LCOW"] = model.fs.costing.LCOW
        outputs["Product_Value"] = model.fs.costing.genericNP.ammonia_recovery_value
    except:
        print("Unable to solve")
        for key in outputs:
            outputs[key] = None
    return outputs


def reinitialize_system(
    model,
    p_removal=0.95,
    n_to_p_ratio=0.3,
    energy_intensity=0.044,
    mgcl2_dosage=0.388,
    water_recovery=0.99,
):
    set_operating_conditions(
        model,
        p_removal=p_removal,
        n_to_p_ratio=n_to_p_ratio,
        energy_intensity=energy_intensity,
        mgcl2_dosage=mgcl2_dosage,
        water_recovery=water_recovery,
    )
    initialize_system(model, has_genericNP=True)


def run_analysis(case_num=11, interpolate_nan_outputs=True, output_filename=None):
    if output_filename is None:
        output_filename = os.path.join(
            this_file_dir, f"initialization_training/sensitivity_{case_num}"
        )

    global_results = parameter_sweep(
        build_model,
        build_sweep_params,
        build_outputs,
        csv_results_file_name=f"{output_filename}.csv",
        optimize_function=solve,
        reinitialize_function=reinitialize_system,
        interpolate_nan_outputs=interpolate_nan_outputs,
    )
    return global_results


def plot_results(results):
    df_results = pd.DataFrame()
    df_results["NH4 removal %"] = results[1]["sweep_params"]["NH4 removal %"]["value"]
    df_results["NH4 removal energy intensity (kwh/kg)"] = results[1]["sweep_params"][
        "NH4 removal energy intensity (kwh/kg)"
    ]["value"]
    df_results["P removal %"] = results[1]["sweep_params"]["P removal %"]["value"]
    df_results["LCOW"] = results[1]["outputs"]["LCOW"]["value"]
    df_results["N value"] = results[1]["outputs"]["Product_Value"]["value"]

    df_results_fixed_p = df_results[
        df_results["P removal %"] == list(df_results["P removal %"])[0]
    ]

    pivot_df_lcow = df_results_fixed_p.pivot(
        index="NH4 removal %",
        columns="NH4 removal energy intensity (kwh/kg)",
        values="LCOW",
    )
    pivot_df_lcow = pivot_df_lcow.round(3)
    pivot_df_lcow.index = pivot_df_lcow.index.round(3)
    pivot_df_lcow.columns = pivot_df_lcow.columns.round(3)

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
        os.path.join(this_file_dir, "sensitivity_heatmap_lcow_nh4_val_0.1.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.show()


if __name__ == "__main__":
    run_parameter_sweep = False  # Set to True to run parameter sweep
    plot_network_before_solve = False  # Set to True to plot network before solving

    if not run_parameter_sweep:
        m, results = main(
            has_genericNP=True,
            plot_network_before_solve=plot_network_before_solve,
            apply_costing=True,
        )

        # Create stream table
        stream_dict = {
            "Feed": m.fs.FeedWater.outlet,
            "Treated": m.fs.Treated.inlet,
            "Sludge": m.fs.Sludge.inlet,
        }
        if m.fs.has_genericNP:
            stream_dict.update(
                {
                    "genericNP_treated": m.fs.genericNP.treated,
                    "genericNP_byproduct": m.fs.genericNP.byproduct,
                }
            )

        stream_table = create_stream_table_dataframe(stream_dict, time_point=0)
        print(stream_table_dataframe_to_string(stream_table))

        # Create and run Dash app for final visualization with stream data
        app = create_dash_app(m, stream_table=stream_table, show_labels=True)
        if app is not None:
            port = find_available_port()
            print(f"Starting final network visualization on port {port}")
            thread = run_dash_app(app, port)
            if thread is not None:
                print(
                    f"Final network visualization available at http://127.0.0.1:{port}"
                )

                # Keep the main thread alive to allow viewing the visualizations
                try:
                    while True:
                        time.sleep(1)
                except KeyboardInterrupt:
                    print("\nShutting down...")
            else:
                print("Failed to start final network visualization")
        else:
            print("Failed to create final network visualization")
    else:
        results = run_analysis()
        plot_results(results)
