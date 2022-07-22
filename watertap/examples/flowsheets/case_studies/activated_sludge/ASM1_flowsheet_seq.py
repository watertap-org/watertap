#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
Example of activated sludge process model using ASM1.


Unit opeations are modeled as follows:

    * Anoixic reactors as standard CSTRs (CSTR)
    * Aerobic reactors as CSTRs with additioanl injection terms (CSTR_Injection)
    * Clarifier as a Separator with split fractions by components

Based on example from:

[1] J. Alex, L. Benedetti, J. Copp, K.V. Gernaey, U. Jeppsson, I. Nopens, M.N. Pons,
J.P. Steyer and P. Vanrolleghem, "Benchmark Simulation Model no. 1 (BSM1)", 2018
"""

# Some more inforation about this module
__author__ = "Andrew Lee"

import pyomo.environ as pyo
from pyomo.network import Arc, SequentialDecomposition

from idaes.core import FlowsheetBlock
from idaes.models.unit_models import CSTR, Feed, Mixer, Separator
from idaes.models.unit_models.separator import SplittingType
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale
from idaes.core.util.initialization import propagate_state

from watertap.unit_models.cstr_injection import CSTR_Injection
from watertap.property_models.activated_sludge.asm1_properties import ASM1ParameterBlock
from watertap.property_models.activated_sludge.asm1_reactions import (
    ASM1ReactionParameterBlock,
)


def build_flowsheet():
    m = pyo.ConcreteModel()

    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.props = ASM1ParameterBlock()
    m.fs.rxn_props = ASM1ReactionParameterBlock(
        default={"property_package": m.fs.props}
    )
    # Feed water stream
    m.fs.FeedWater = Feed(
        default={
            "property_package": m.fs.props,
        }
    )
    m.fs.FeedWater.flow_vol.fix(18446 * pyo.units.m**3 / pyo.units.day)
    m.fs.FeedWater.temperature.fix(298.15 * pyo.units.K)
    m.fs.FeedWater.pressure.fix(1 * pyo.units.atm)
    m.fs.FeedWater.conc_mass_comp[0, "S_I"].fix(30 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_S"].fix(69.5 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_I"].fix(51.2 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_S"].fix(202.32 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_BH"].fix(28.17 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_BA"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_P"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_O"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_NO"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_NH"].fix(31.56 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_ND"].fix(6.95 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_ND"].fix(10.59 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.alkalinity.fix(7 * pyo.units.mol / pyo.units.m**3)
    m.fs.FeedWater.initialize()

    # First reactor (anoxic) - standard CSTR
    m.fs.R1 = CSTR(
        default={
            "property_package": m.fs.props,
            "reaction_package": m.fs.rxn_props,
        }
    )
    # Feed conditions based on manual mass balance of inlet and recycle streams
    m.fs.R1.inlet.flow_vol.fix(92230 * pyo.units.m**3 / pyo.units.day)
    m.fs.R1.inlet.temperature.fix(298.15 * pyo.units.K)
    m.fs.R1.inlet.pressure.fix(1 * pyo.units.atm)
    m.fs.R1.inlet.conc_mass_comp[0, "S_I"].fix(30 * pyo.units.g / pyo.units.m**3)
    m.fs.R1.inlet.conc_mass_comp[0, "S_S"].fix(14.6112 * pyo.units.g / pyo.units.m**3)
    m.fs.R1.inlet.conc_mass_comp[0, "X_I"].fix(1149 * pyo.units.g / pyo.units.m**3)
    m.fs.R1.inlet.conc_mass_comp[0, "X_S"].fix(89.324 * pyo.units.g / pyo.units.m**3)
    m.fs.R1.inlet.conc_mass_comp[0, "X_BH"].fix(
        2542.03 * pyo.units.g / pyo.units.m**3
    )
    m.fs.R1.inlet.conc_mass_comp[0, "X_BA"].fix(148.6 * pyo.units.g / pyo.units.m**3)
    m.fs.R1.inlet.conc_mass_comp[0, "X_P"].fix(448 * pyo.units.g / pyo.units.m**3)
    m.fs.R1.inlet.conc_mass_comp[0, "S_O"].fix(0.3928 * pyo.units.g / pyo.units.m**3)
    m.fs.R1.inlet.conc_mass_comp[0, "S_NO"].fix(8.32 * pyo.units.g / pyo.units.m**3)
    m.fs.R1.inlet.conc_mass_comp[0, "S_NH"].fix(7.696 * pyo.units.g / pyo.units.m**3)
    m.fs.R1.inlet.conc_mass_comp[0, "S_ND"].fix(1.9404 * pyo.units.g / pyo.units.m**3)
    m.fs.R1.inlet.conc_mass_comp[0, "X_ND"].fix(5.616 * pyo.units.g / pyo.units.m**3)
    m.fs.R1.inlet.alkalinity.fix(4.704 * pyo.units.mol / pyo.units.m**3)
    m.fs.R1.volume.fix(1000 * pyo.units.m**3)
    m.fs.R1.initialize()

    solve_fs(m)

    # Second reactor (anoxic) - standard CSTR
    m.fs.R2 = CSTR(
        default={
            "property_package": m.fs.props,
            "reaction_package": m.fs.rxn_props,
        }
    )
    m.fs.R2.volume.fix(1000 * pyo.units.m**3)
    m.fs.stream3 = Arc(source=m.fs.R1.outlet, destination=m.fs.R2.inlet)
    pyo.TransformationFactory("network.expand_arcs").apply_to(m)
    propagate_state(m.fs.stream3)
    m.fs.R2.initialize()

    results = solve_fs(m)

    # Third reactor (aerobic) - CSTR with injection
    m.fs.R3 = CSTR_Injection(
        default={
            "property_package": m.fs.props,
            "reaction_package": m.fs.rxn_props,
        }
    )
    m.fs.R3.volume.fix(1333 * pyo.units.m**3)
    m.fs.R3.injection[...].fix(0)
    m.fs.R3.injection[:, "Liq", "S_O"].fix(0.022)
    m.fs.stream4 = Arc(source=m.fs.R2.outlet, destination=m.fs.R3.inlet)
    pyo.TransformationFactory("network.expand_arcs").apply_to(m)
    propagate_state(m.fs.stream4)
    m.fs.R3.initialize()

    results = solve_fs(m)

    # Fourth reactor (aerobic) - CSTR with injection
    m.fs.R4 = CSTR_Injection(
        default={
            "property_package": m.fs.props,
            "reaction_package": m.fs.rxn_props,
        }
    )
    m.fs.R4.volume.fix(1333 * pyo.units.m**3)
    m.fs.R4.injection[...].fix(0)
    m.fs.R4.injection[:, "Liq", "S_O"].fix(0.021)
    m.fs.stream5 = Arc(source=m.fs.R3.outlet, destination=m.fs.R4.inlet)
    pyo.TransformationFactory("network.expand_arcs").apply_to(m)
    propagate_state(m.fs.stream5)
    m.fs.R4.initialize()

    # Fifth reactor (aerobic) - CSTR with injection
    m.fs.R5 = CSTR_Injection(
        default={
            "property_package": m.fs.props,
            "reaction_package": m.fs.rxn_props,
        }
    )
    m.fs.R5.volume.fix(1333 * pyo.units.m**3)
    m.fs.R5.injection[...].fix(0)
    m.fs.R5.injection[:, "Liq", "S_O"].fix(0.01)
    m.fs.stream6 = Arc(source=m.fs.R4.outlet, destination=m.fs.R5.inlet)
    pyo.TransformationFactory("network.expand_arcs").apply_to(m)
    propagate_state(m.fs.stream6)
    m.fs.R5.initialize()

    results = solve_fs(m)

    m.fs.SP5 = Separator(
        default={
            "property_package": m.fs.props,
            "outlet_list": ["underflow", "overflow"],
        }
    )
    m.fs.SP5.split_fraction[:, "underflow"].fix(0.6)
    m.fs.stream7 = Arc(source=m.fs.R5.outlet, destination=m.fs.SP5.inlet)
    pyo.TransformationFactory("network.expand_arcs").apply_to(m)
    propagate_state(m.fs.stream7)
    m.fs.SP5.initialize()

    results = solve_fs(m)

    # Clarifier
    # TODO: Replace with more detailed model when available
    m.fs.CL1 = Separator(
        default={
            "property_package": m.fs.props,
            "outlet_list": ["underflow", "effluent"],
            "split_basis": SplittingType.componentFlow,
        }
    )
    # TODO: Update once more detailed model available
    m.fs.CL1.split_fraction[0, "effluent", "H2O"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "S_I"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "S_S"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "X_I"].fix(0.00187)
    m.fs.CL1.split_fraction[0, "effluent", "X_S"].fix(0.00187)
    m.fs.CL1.split_fraction[0, "effluent", "X_BH"].fix(0.00187)
    m.fs.CL1.split_fraction[0, "effluent", "X_BA"].fix(0.00187)
    m.fs.CL1.split_fraction[0, "effluent", "X_P"].fix(0.00187)
    m.fs.CL1.split_fraction[0, "effluent", "S_O"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "S_NO"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "S_NH"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "S_ND"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "X_ND"].fix(0.00187)
    m.fs.CL1.split_fraction[0, "effluent", "S_ALK"].fix(0.48956)
    m.fs.stream8 = Arc(source=m.fs.SP5.overflow, destination=m.fs.CL1.inlet)
    pyo.TransformationFactory("network.expand_arcs").apply_to(m)
    propagate_state(m.fs.stream8)
    m.fs.CL1.initialize()

    results = solve_fs(m)

    # Sludge purge splitter
    m.fs.SP6 = Separator(
        default={
            "property_package": m.fs.props,
            "outlet_list": ["recycle", "waste"],
            "split_basis": SplittingType.totalFlow,
        }
    )
    m.fs.SP6.split_fraction[:, "recycle"].fix(0.97955)
    m.fs.stream11 = Arc(source=m.fs.CL1.underflow, destination=m.fs.SP6.inlet)
    pyo.TransformationFactory("network.expand_arcs").apply_to(m)
    propagate_state(m.fs.stream11)
    m.fs.SP6.initialize()

    results = solve_fs(m)

    # Mixing sludge recycle and R5 underflow
    m.fs.MX6 = Mixer(
        default={
            "property_package": m.fs.props,
            "inlet_list": ["clarifier", "reactor"],
        }
    )
    m.fs.stream9 = Arc(source=m.fs.SP5.underflow, destination=m.fs.MX6.reactor)
    m.fs.stream13 = Arc(source=m.fs.SP6.recycle, destination=m.fs.MX6.clarifier)
    pyo.TransformationFactory("network.expand_arcs").apply_to(m)
    propagate_state(m.fs.stream9)
    propagate_state(m.fs.stream13)
    m.fs.MX6.initialize()

    results = solve_fs(m)

    # Mixer for feed water and recycled sludge
    m.fs.MX1 = Mixer(
        default={
            "property_package": m.fs.props,
            "inlet_list": ["feed_water", "recycle"],
        }
    )
    m.fs.stream1 = Arc(source=m.fs.FeedWater.outlet, destination=m.fs.MX1.feed_water)
    m.fs.stream14 = Arc(source=m.fs.MX6.outlet, destination=m.fs.MX1.recycle)
    pyo.TransformationFactory("network.expand_arcs").apply_to(m)
    propagate_state(m.fs.stream1)
    propagate_state(m.fs.stream14)
    m.fs.MX1.initialize()

    results = solve_fs(m)

    # Connect recycle
    # m.fs.stream2 = Arc(source=m.fs.MX1.outlet, destination=m.fs.R1.inlet)
    # pyo.TransformationFactory("network.expand_arcs").apply_to(m)
    m.fs.recyc_T = pyo.Constraint(
        expr=m.fs.MX1.outlet.temperature[0] == m.fs.R1.inlet.temperature[0]
    )
    m.fs.R1.inlet.temperature.unfix()

    # m.fs.recyc_P = pyo.Constraint(expr=m.fs.MX1.outlet.pressure[0] == m.fs.R1.inlet.pressure[0])
    # m.fs.R1.inlet.pressure.unfix()

    # m.fs.R1.inlet.flow_vol.unfix()
    # m.fs.R1.inlet.conc_mass_comp.unfix()
    # m.fs.R1.inlet.alkalinity.unfix()

    # Apply scaling
    iscale.calculate_scaling_factors(m.fs)

    results = solve_fs(m)

    return m, results


def solve_fs(m):
    # Check degrees of freedom
    print(degrees_of_freedom(m))
    assert degrees_of_freedom(m) == 0
    solver = get_solver()
    results = solver.solve(m, tee=True)

    pyo.assert_optimal_termination(results)

    return results


def temp():

    # Initialize flowsheet
    # Apply sequential decomposition - 3 iterations should suffice
    seq = SequentialDecomposition()
    seq.options.select_tear_method = "heuristic"
    seq.options.tear_method = "Wegstein"
    seq.options.iterLim = 3

    G = seq.create_graph(m)

    # # Uncomment this code to see tear set and initialization order
    # heuristic_tear_set = seq.tear_set_arcs(G, method="heuristic")
    # order = seq.calculation_order(G)
    # for o in heuristic_tear_set:
    #     print(o.name)
    # for o in order:
    #     print(o[0].name)

    # Initial guesses for flow into first reactor
    tear_guesses = {
        "flow_vol": {0: 1.0675},
        "conc_mass_comp": {
            (0, "S_I"): 0.03,
            (0, "S_S"): 0.0146,
            (0, "X_I"): 1.149,
            (0, "X_S"): 0.0893,
            (0, "X_BH"): 2.542,
            (0, "X_BA"): 0.149,
            (0, "X_P"): 0.448,
            (0, "S_O"): 3.93e-4,
            (0, "S_NO"): 8.32e-3,
            (0, "S_NH"): 7.7e-3,
            (0, "S_ND"): 1.94e-3,
            (0, "X_ND"): 5.62e-3,
        },
        "alkalinity": {0: 5e-3},
        "temperature": {0: 298.15},
        "pressure": {0: 101325},
    }

    # Pass the tear_guess to the SD tool
    seq.set_guesses_for(m.fs.R1.inlet, tear_guesses)

    def function(unit):
        if unit.local_name == "MX1":
            unit.report()
            unit.initialize(outlvl=idaeslog.INFO)
            unit.report()
        else:
            unit.initialize(outlvl=idaeslog.INFO)

    seq.run(m, function)

    # Solve overall flowsheet to close recycle loop
    solver = get_solver()
    results = solver.solve(m, tee=True)

    pyo.assert_optimal_termination(results)

    # Oxygen concentration in reactors 3 and 4 is governed by mass transfer
    # Add additional parameter and constraints
    m.fs.KLa = pyo.Param(
        default=10 / 3600,
        units=pyo.units.s**-1,
        mutable=True,
        doc="Lumped mass transfer coefficient for oxygen",
    )
    m.fs.S_O_eq = pyo.Param(
        default=8e-3,
        units=pyo.units.kg / pyo.units.m**3,
        mutable=True,
        doc="Dissolved oxygen concentration at equilibrium",
    )

    @m.fs.Constraint(m.fs.time, doc="Mass transfer constraint for R3")
    def mass_transfer_R3(self, t):
        return m.fs.R3.injection[t, "Liq", "S_O"] == (
            m.fs.KLa
            * m.fs.R3.volume[t]
            * (m.fs.S_O_eq - m.fs.R3.outlet.conc_mass_comp[t, "S_O"])
        )

    m.fs.R3.injection[:, "Liq", "S_O"].unfix()

    @m.fs.Constraint(m.fs.time, doc="Mass transfer constraint for R4")
    def mass_transfer_R4(self, t):
        return m.fs.R4.injection[t, "Liq", "S_O"] == (
            m.fs.KLa
            * m.fs.R4.volume[t]
            * (m.fs.S_O_eq - m.fs.R4.outlet.conc_mass_comp[t, "S_O"])
        )

    m.fs.R4.injection[:, "Liq", "S_O"].unfix()

    # Oxygen concentration in R5 is controled
    # Fix outlet concentration and unfix injection
    m.fs.R5.outlet.conc_mass_comp[0, "S_O"].fix(0.491 * pyo.units.g / pyo.units.m**3)
    m.fs.R5.injection[:, "Liq", "S_O"].unfix()

    # Resolve with controls in place
    results = solver.solve(m, tee=True)
    pyo.assert_optimal_termination(results)

    return m, results


if __name__ == "__main__":
    # This method builds and runs a steady state activated sludge
    # flowsheet.
    m, results = build_flowsheet()

    # m.fs.R2.report()
    # m.fs.R3.report()
    # m.fs.R4.report()
    m.fs.CL1.report()
