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
Example of activated sludge process model.
Initial layout:
    * 5 reactors
        * R1, R3, R5: anoxic
        * R2 and R4: aerobic
    * 2 Mixers: 
        * M1 mixes feed and R5 split fraction, outlet to R1 inlet
        * M2 mixes R1 outlet, R2 outlet, outlet to R3 inlet
    * 1 Splitter:
        * separates R5 outlet into effluent, M1 inlet, and R2 inlet


Unit operations are modeled as follows:

    * Anoxic reactors as standard CSTRs (CSTR)
    * Aerobic reactors as AerationTanks

Note that a pressure changer is required in the recycle stream to ensure the
pressure inside the recycle loop is bounded. As the inlet Mixer uses a pressure
minimization constraint and there is no pressure drop in the reactors, if pressure
is not specified at some point within the recycle loop then it becomes unbounded.

"""

__author__ = "Adam Atia"

import pyomo.environ as pyo
from pyomo.network import Arc, SequentialDecomposition

from idaes.core import FlowsheetBlock
from idaes.models.unit_models import (
    CSTR,
    Feed,
    Mixer,
    Separator,
    PressureChanger,
    Product,
    MomentumMixingType
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
from idaes.core.util.initialization import propagate_state

from watertap.unit_models import AerationTank

from watertap.property_models.unit_specific.activated_sludge.asm1_properties import (
    ASM1ParameterBlock,
)
from watertap.property_models.unit_specific.activated_sludge.asm1_reactions import (
    ASM1ReactionParameterBlock,
)

from watertap.core.util.initialization import check_solve, interval_initializer

# Set up logger
_log = idaeslog.getLogger(__name__)


def build_flowsheet():
    m = pyo.ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.props = ASM1ParameterBlock()
    m.fs.rxn_props = ASM1ReactionParameterBlock(property_package=m.fs.props)
    # Feed water stream
    m.fs.feed = Feed(property_package=m.fs.props)

    # Mixer for feed water and recycled sludge
    m.fs.M1 = Mixer(property_package=m.fs.props, 
                    inlet_list=["feed", "recycle"],
                    momentum_mixing_type=MomentumMixingType.none
                    )
    m.fs.M1.outlet.pressure.fix()
    # m.fs.M1.pressure_equality_constraints[0,2].deactivate()
    # First reactor (anoxic) 
    m.fs.R1 = CSTR(property_package=m.fs.props, reaction_package=m.fs.rxn_props)
    
    # Second reactor (aerobic) 
    m.fs.R2 = AerationTank(property_package=m.fs.props, reaction_package=m.fs.rxn_props)
    
    
    # Mixer for R1 and R2 outlets
    m.fs.M2= Mixer(property_package=m.fs.props, inlet_list=["R1_outlet", "R2_outlet"],
                   momentum_mixing_type=MomentumMixingType.none
                    )
    m.fs.M2.outlet.pressure.fix()

    # Third reactor (anoxic)
    m.fs.R3 = CSTR(
        property_package=m.fs.props, reaction_package=m.fs.rxn_props
    )
    # Fourth reactor (aerobic) - CSTR with injection
    m.fs.R4 = AerationTank(
        property_package=m.fs.props, reaction_package=m.fs.rxn_props
    )
    # Fifth reactor (anoxic) 
    m.fs.R5 = CSTR(
        property_package=m.fs.props, reaction_package=m.fs.rxn_props
    )
    m.fs.S1 = Separator(
        property_package=m.fs.props, outlet_list=["effluent","M1_inlet", "R2_inlet"]
    )
    # Clarifier
    # # TODO: Replace with more detailed model when available
    # m.fs.CL1 = Separator(
    #     property_package=m.fs.props,
    #     outlet_list=["underflow", "effluent"],
    #     split_basis=SplittingType.componentFlow,
    # )
    # # Sludge purge splitter
    # m.fs.SP6 = Separator(
    #     property_package=m.fs.props,
    #     outlet_list=["recycle", "waste"],
    #     split_basis=SplittingType.totalFlow,
    # )
    # # Mixing sludge recycle and R5 underflow
    # m.fs.MX6 = Mixer(property_package=m.fs.props, inlet_list=["clarifier", "reactor"])
    
    # Product Blocks
    m.fs.Treated = Product(property_package=m.fs.props)
    # Recycle pressure changer - use a simple isothermal unit for now
    # m.fs.P1 = PressureChanger(property_package=m.fs.props)

    # Link units
    m.fs.feed_to_m1 = Arc(source=m.fs.feed.outlet, destination=m.fs.M1.feed)
    m.fs.m1_to_r1 = Arc(source=m.fs.M1.outlet, destination=m.fs.R1.inlet)
    m.fs.r1_to_m2 = Arc(source=m.fs.R1.outlet, destination=m.fs.M2.R1_outlet)
    m.fs.r2_to_m2 = Arc(source=m.fs.R2.outlet, destination=m.fs.M2.R2_outlet)
    m.fs.m2_to_r3 = Arc(source=m.fs.M2.outlet, destination=m.fs.R3.inlet)
    m.fs.r3_to_r4 = Arc(source=m.fs.R3.outlet, destination=m.fs.R4.inlet)
    m.fs.r4_to_r5 = Arc(source=m.fs.R4.outlet, destination=m.fs.R5.inlet)
    m.fs.r5_to_s1 = Arc(source=m.fs.R5.outlet, destination=m.fs.S1.inlet)
    m.fs.s1_to_effluent = Arc(source=m.fs.S1.effluent, destination=m.fs.Treated.inlet)
    m.fs.s1_to_m1 = Arc(source=m.fs.S1.M1_inlet, destination=m.fs.M1.recycle)
    m.fs.s1_to_r2 = Arc(source=m.fs.S1.R2_inlet, destination=m.fs.R2.inlet)
    # m.fs.stream14 = Arc(source=m.fs.MX6.outlet, destination=m.fs.P1.inlet)
    # m.fs.stream15 = Arc(source=m.fs.P1.outlet, destination=m.fs.MX1.recycle)
    pyo.TransformationFactory("network.expand_arcs").apply_to(m)
    
    return m

def set_operating_conditions(m):
    # Feed Water Conditions
    m.fs.feed.flow_vol.fix(18446 * pyo.units.m**3 / pyo.units.day)
    m.fs.feed.temperature.fix(298.15 * pyo.units.K)
    m.fs.feed.pressure.fix(1 * pyo.units.atm)
    m.fs.feed.conc_mass_comp[0, "S_I"].fix(30 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "S_S"].fix(69.5 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "X_I"].fix(51.2 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "X_S"].fix(202.32 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "X_BH"].fix(28.17 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "X_BA"].fix(0 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "X_P"].fix(0 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "S_O"].fix(0 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "S_NO"].fix(0 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "S_NH"].fix(31.56 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "S_ND"].fix(6.95 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "X_ND"].fix(10.59 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.alkalinity.fix(7 * pyo.units.mol / pyo.units.m**3)

    # Reactor sizing
    m.fs.R1.volume.fix(1000 * pyo.units.m**3)
    m.fs.R2.volume.fix(1000 * pyo.units.m**3)
    m.fs.R3.volume.fix(1333 * pyo.units.m**3)
    m.fs.R4.volume.fix(1333 * pyo.units.m**3)
    m.fs.R5.volume.fix(1333 * pyo.units.m**3)

    # Injection rates to Reactors 2 and 4
    for j in m.fs.props.component_list:
        if j != "S_O":
            # All components except S_O have no injection
            m.fs.R2.injection[:, :, j].fix(0)
            m.fs.R4.injection[:, :, j].fix(0)
    # Then set injections rates for O2
    m.fs.R2.outlet.conc_mass_comp[:, "S_O"].fix(1.72e-3)
    m.fs.R4.outlet.conc_mass_comp[:, "S_O"].fix(2.43e-3)

    # Set fraction of outflow from reactor 5 that goes to recycle
    m.fs.S1.split_fraction[:, "M1_inlet"].fix(0.5)

    # Set fraction of outflow from reactor 5 that goes to effluent
    m.fs.S1.split_fraction[:, "effluent"].fix(0.4)

    # # Outlet pressure from recycle pump
    # m.fs.P1.outlet.pressure.fix(101325)

    # Check degrees of freedom
    print("DOF = ", degrees_of_freedom(m))
    assert degrees_of_freedom(m) == 0

def scale_flowsheet(m):
    # Apply scaling
    for var in m.fs.component_data_objects(pyo.Var, descend_into=True):
        if "flow_vol" in var.name:
            iscale.set_scaling_factor(var, 1e1)
        if "temperature" in var.name:
            iscale.set_scaling_factor(var, 1e-1)
        if "pressure" in var.name:
            iscale.set_scaling_factor(var, 1e-4)
        if "conc_mass_comp" in var.name:
            iscale.set_scaling_factor(var, 1e2)
    iscale.calculate_scaling_factors(m.fs)

def init_and_propagate(blk,arc=None,source=None,destination=None):
    blk.initialize()
    if arc is not None:
        propagate_state(arc)
    elif source is not None:
        propagate_state(source=source,destination=destination)


def initialize_flowsheet(m):
    # Initialize flowsheet
    # interval_initializer(m)
    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_m1)
    propagate_state(source=m.fs.feed.outlet, destination=m.fs.M1.recycle)

    m.fs.M1.initialize()
    propagate_state(m.fs.m1_to_r1)

    m.fs.R1.initialize()
    propagate_state(m.fs.r1_to_m2)
    propagate_state(source=m.fs.R1.outlet, destination=m.fs.M2.R2_outlet)

    m.fs.M2.initialize()
    propagate_state(m.fs.m2_to_r3)

    m.fs.R3.initialize() 
    propagate_state(m.fs.r3_to_r4)

    m.fs.R4.initialize()
    propagate_state(m.fs.r4_to_r5)

    m.fs.R5.initialize()
    propagate_state(m.fs.r5_to_s1)

    m.fs.S1.initialize()
    propagate_state(m.fs.s1_to_effluent)
    propagate_state(m.fs.s1_to_m1)
    propagate_state(m.fs.s1_to_r2)

    m.fs.R2.initialize()
    propagate_state(m.fs.r2_to_m2)

    m.fs.M1.initialize()
    propagate_state(m.fs.m1_to_r1)

    m.fs.R1.initialize()
    propagate_state(m.fs.r1_to_m2)

    m.fs.M2.initialize()
    propagate_state(m.fs.m2_to_r3)

    m.fs.R3.initialize() 
    propagate_state(m.fs.r3_to_r4)

    m.fs.R4.initialize()
    propagate_state(m.fs.r4_to_r5)

    m.fs.R5.initialize()
    propagate_state(m.fs.r5_to_s1)

    m.fs.S1.initialize()
    propagate_state(m.fs.s1_to_effluent)
    propagate_state(m.fs.s1_to_m1)
    propagate_state(m.fs.s1_to_r2)
    
    m.fs.Treated.initialize()

    interval_initializer(m)


    # Apply sequential decomposition - 1 iteration should suffice
    seq = SequentialDecomposition()
    seq.options.select_tear_method = "heuristic"
    seq.options.tear_method = "Wegstein"
    seq.options.iterLim = 1

    G = seq.create_graph(m)

    # # Uncomment this code to see tear set and initialization order
    # heuristic_tear_set = seq.tear_set_arcs(G, method="heuristic")
    # order = seq.calculation_order(G)
    # for o in heuristic_tear_set:
    #     print(o.name)
    # for o in order:
    #     print(o[0].name)

    # # Initial guesses for flow into first reactor
    # tear_guesses = {
    #     "flow_vol": {0: 1.0675},
    #     "conc_mass_comp": {
    #         (0, "S_I"): 0.03,
    #         (0, "S_S"): 0.0146,
    #         (0, "X_I"): 1.149,
    #         (0, "X_S"): 0.0893,
    #         (0, "X_BH"): 2.542,
    #         (0, "X_BA"): 0.149,
    #         (0, "X_P"): 0.448,
    #         (0, "S_O"): 3.93e-4,
    #         (0, "S_NO"): 8.32e-3,
    #         (0, "S_NH"): 7.7e-3,
    #         (0, "S_ND"): 1.94e-3,
    #         (0, "X_ND"): 5.62e-3,
    #     },
    #     "alkalinity": {0: 5e-3},
    #     "temperature": {0: 298.15},
    #     "pressure": {0: 101325},
    # }

    # # Pass the tear_guess to the SD tool
    # seq.set_guesses_for(m.fs.R1.inlet, tear_guesses)

    def function(unit):
        unit.initialize(outlvl=idaeslog.DEBUG)

    seq.run(m, function)

def solve_flowsheet(m):
    # Solve overall flowsheet to close recycle loop
    solver = get_solver()
    results = solver.solve(m)
    check_solve(results, checkpoint="closing recycle", logger=_log, fail_flag=True)

    # Switch to fixed KLa in R3 and R4 (S_O concentration is controlled in R5)
    m.fs.R2.KLa.fix(10)
    m.fs.R4.KLa.fix(10)
    m.fs.R2.outlet.conc_mass_comp[:, "S_O"].unfix()
    m.fs.R4.outlet.conc_mass_comp[:, "S_O"].unfix()

    # Resolve with controls in place
    results = solver.solve(m, tee=True)
    check_solve(
        results,
        checkpoint="re-solve with controls in place",
        logger=_log,
        fail_flag=True,
    )

    return m, results


if __name__ == "__main__":
    # This method builds and runs a steady state activated sludge
    # flowsheet.
    # m, results = build_flowsheet()
    m = build_flowsheet()
    set_operating_conditions(m)
    scale_flowsheet(m)

    initialize_flowsheet(m)

    scale_flowsheet(m)

    m, res = solve_flowsheet(m)


    stream_table = create_stream_table_dataframe(
        {
            "Feed": m.fs.feed.outlet,
            "M1": m.fs.M1.outlet,
            "R1": m.fs.R1.outlet,
            "M2": m.fs.M2.outlet,
            "R2": m.fs.R2.outlet,
            "R3": m.fs.R3.outlet,
            "R4": m.fs.R4.outlet,
            "R5": m.fs.R5.outlet,
            "S1 to M1": m.fs.S1.M1_inlet,
            "S1 to R2": m.fs.S1.R2_inlet,
            "Effluent": m.fs.Treated.inlet

        },
        time_point=0,
    )
    print(stream_table_dataframe_to_string(stream_table))
