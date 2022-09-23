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

"""

# Some more information about this module
__author__ = "Alejandro Garciadiego, Andrew Lee"

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
)
from idaes.models.unit_models.separator import SplittingType
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom, large_residuals_set
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale
from idaes.core.util.tables import (
    create_stream_table_dataframe,
    stream_table_dataframe_to_string,
)

from watertap.unit_models.cstr_injection import CSTR_Injection
from watertap.property_models.activated_sludge.asm2d_properties import (
    ASM2dParameterBlock,
)
from watertap.property_models.activated_sludge.asm2d_reactions import (
    ASM2dReactionParameterBlock,
)


def build_flowsheet():
    m = pyo.ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.props = ASM2dParameterBlock()
    m.fs.rxn_props = ASM2dReactionParameterBlock(property_package=m.fs.props)

    # Feed water stream
    m.fs.FeedWater = Feed(property_package=m.fs.props)
    # Mixer for feed water and recycled sludge
    m.fs.MX1 = Mixer(property_package=m.fs.props, inlet_list=["feed_water", "recycle"])
    # First reactor (anoxic) - standard CSTR
    m.fs.R1 = CSTR(property_package=m.fs.props, reaction_package=m.fs.rxn_props)
    # Second reactor (anoxic) - standard CSTR
    m.fs.R2 = CSTR(property_package=m.fs.props, reaction_package=m.fs.rxn_props)
    # Third reactor (aerobic) - CSTR with injection
    m.fs.R3 = CSTR_Injection(
        property_package=m.fs.props, reaction_package=m.fs.rxn_props
    )
    # Fourth reactor (aerobic) - CSTR with injection
    m.fs.R4 = CSTR_Injection(
        property_package=m.fs.props, reaction_package=m.fs.rxn_props
    )
    # Fifth reactor (aerobic) - CSTR with injection
    m.fs.R5 = CSTR_Injection(
        property_package=m.fs.props, reaction_package=m.fs.rxn_props
    )
    m.fs.SP5 = Separator(
        property_package=m.fs.props, outlet_list=["underflow", "overflow"]
    )
    # Clarifier
    # TODO: Replace with more detailed model when available
    m.fs.CL1 = Separator(
        property_package=m.fs.props,
        outlet_list=["underflow", "effluent"],
        split_basis=SplittingType.componentFlow,
    )
    # Sludge purge splitter
    m.fs.SP6 = Separator(
        property_package=m.fs.props,
        outlet_list=["recycle", "waste"],
        split_basis=SplittingType.totalFlow,
    )
    # Mixing sludge recycle and R5 underflow
    m.fs.MX6 = Mixer(property_package=m.fs.props, inlet_list=["clarifier", "reactor"])
    # Product Blocks
    m.fs.Treated = Product(property_package=m.fs.props)
    m.fs.Sludge = Product(property_package=m.fs.props)
    # Recycle pressure changer - use a simple isothermal unit for now
    m.fs.P1 = PressureChanger(property_package=m.fs.props)

    # Link units
    m.fs.stream1 = Arc(source=m.fs.FeedWater.outlet, destination=m.fs.MX1.feed_water)
    m.fs.stream2 = Arc(source=m.fs.MX1.outlet, destination=m.fs.R1.inlet)
    m.fs.stream3 = Arc(source=m.fs.R1.outlet, destination=m.fs.R2.inlet)
    m.fs.stream4 = Arc(source=m.fs.R2.outlet, destination=m.fs.R3.inlet)
    m.fs.stream5 = Arc(source=m.fs.R3.outlet, destination=m.fs.R4.inlet)
    m.fs.stream6 = Arc(source=m.fs.R4.outlet, destination=m.fs.R5.inlet)
    m.fs.stream7 = Arc(source=m.fs.R5.outlet, destination=m.fs.SP5.inlet)
    m.fs.stream8 = Arc(source=m.fs.SP5.overflow, destination=m.fs.CL1.inlet)
    m.fs.stream9 = Arc(source=m.fs.SP5.underflow, destination=m.fs.MX6.reactor)
    m.fs.stream10 = Arc(source=m.fs.CL1.effluent, destination=m.fs.Treated.inlet)
    m.fs.stream11 = Arc(source=m.fs.CL1.underflow, destination=m.fs.SP6.inlet)
    m.fs.stream102 = Arc(source=m.fs.SP6.waste, destination=m.fs.Sludge.inlet)
    m.fs.stream13 = Arc(source=m.fs.SP6.recycle, destination=m.fs.MX6.clarifier)
    m.fs.stream14 = Arc(source=m.fs.MX6.outlet, destination=m.fs.P1.inlet)
    m.fs.stream15 = Arc(source=m.fs.P1.outlet, destination=m.fs.MX1.recycle)
    pyo.TransformationFactory("network.expand_arcs").apply_to(m)

    # Oxygen concentration in reactors 3 and 4 is governed by mass transfer
    # Add additional parameter and constraints
    m.fs.R3.KLa = pyo.Var(
        initialize=3.4,
        units=pyo.units.hour**-1,
        doc="Lumped mass transfer coefficient for oxygen",
    )
    m.fs.R4.KLa = pyo.Var(
        initialize=1.8,
        units=pyo.units.hour**-1,
        doc="Lumped mass transfer coefficient for oxygen",
    )
    m.fs.S_O_eq = pyo.Param(
        default=8e-3,
        units=pyo.units.kg / pyo.units.m**3,
        mutable=True,
        doc="Dissolved oxygen concentration at equilibrium",
    )

    m.fs.R3.Constraint(m.fs.time, doc="Mass transfer constraint for R3")

    def mass_transfer_R3(self, t):
        return pyo.units.convert(
            m.fs.R3.injection[t, "Liq", "S_O2"], to_units=pyo.units.kg / pyo.units.hour
        ) == (
            m.fs.R3.KLa
            * m.fs.R3.volume[t]
            * (m.fs.S_O_eq - m.fs.R3.outlet.conc_mass_comp[t, "S_O2"])
        )

    m.fs.R4.Constraint(m.fs.time, doc="Mass transfer constraint for R4")

    def mass_transfer_R4(self, t):
        return pyo.units.convert(
            m.fs.R4.injection[t, "Liq", "S_O2"], to_units=pyo.units.kg / pyo.units.hour
        ) == (
            m.fs.R4.KLa
            * m.fs.R4.volume[t]
            * (m.fs.S_O_eq - m.fs.R4.outlet.conc_mass_comp[t, "S_O2"])
        )

    # Feed Water Conditions
    m.fs.FeedWater.flow_vol.fix(18446 * pyo.units.m**3 / pyo.units.day)
    m.fs.FeedWater.temperature.fix(298.15 * pyo.units.K)
    m.fs.FeedWater.pressure.fix(1 * pyo.units.atm)
    m.fs.FeedWater.conc_mass_comp[0, "S_O2"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_F"].fix(30 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_A"].fix(20 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_NH4"].fix(16 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_NO3"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_PO4"].fix(3.6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_I"].fix(30 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_N2"].fix(15 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_I"].fix(25 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_S"].fix(125 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_H"].fix(30 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_PAO"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_PP"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_PHA"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_AUT"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_MeOH"].fix(
        1e-6 * pyo.units.g / pyo.units.m**3
    )
    m.fs.FeedWater.conc_mass_comp[0, "X_MeP"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_TSS"].fix(180 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.alkalinity.fix(7 * pyo.units.mol / pyo.units.m**3)

    # Reactor sizing
    m.fs.R1.volume.fix(1000 * pyo.units.m**3)
    m.fs.R2.volume.fix(1000 * pyo.units.m**3)
    m.fs.R3.volume.fix(1333 * pyo.units.m**3)
    m.fs.R4.volume.fix(1333 * pyo.units.m**3)
    m.fs.R5.volume.fix(1333 * pyo.units.m**3)

    # Injection rates to Reactions 3, 4 and 5
    for j in m.fs.props.component_list:
        if j != "S_O2":
            # All components except S_O have no injection
            m.fs.R3.injection[:, :, j].fix(0)
            m.fs.R4.injection[:, :, j].fix(0)
            m.fs.R5.injection[:, :, j].fix(0)
    # Then set injections rates for O2
    m.fs.R3.outlet.conc_mass_comp[:, "S_O2"].fix(2.72e-3)
    m.fs.R4.outlet.conc_mass_comp[:, "S_O2"].fix(2.43e-3)
    m.fs.R5.outlet.conc_mass_comp[:, "S_O2"].fix(4.49e-4)

    # Set fraction of outflow from reactor 5 that goes to recycle
    m.fs.SP5.split_fraction[:, "underflow"].fix(0.1)

    # Clarifier
    # TODO: Update once more detailed model available
    m.fs.CL1.split_fraction[0, "effluent", "H2O"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "S_A"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "S_F"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "S_I"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "S_N2"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "S_NH4"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "S_NO3"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "S_O2"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "S_PO4"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "S_ALK"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "X_AUT"].fix(0.00187)
    m.fs.CL1.split_fraction[0, "effluent", "X_H"].fix(0.00187)
    m.fs.CL1.split_fraction[0, "effluent", "X_I"].fix(0.00187)
    m.fs.CL1.split_fraction[0, "effluent", "X_MeOH"].fix(0.00187)
    m.fs.CL1.split_fraction[0, "effluent", "X_MeP"].fix(0.00187)
    m.fs.CL1.split_fraction[0, "effluent", "X_PAO"].fix(0.00187)
    m.fs.CL1.split_fraction[0, "effluent", "X_PHA"].fix(0.00187)
    m.fs.CL1.split_fraction[0, "effluent", "X_PP"].fix(0.00187)
    m.fs.CL1.split_fraction[0, "effluent", "X_S"].fix(0.00187)
    m.fs.CL1.split_fraction[0, "effluent", "X_TSS"].fix(0.00187)

    # Sludge purge separator
    m.fs.SP6.split_fraction[:, "recycle"].fix(0.65)

    # Outlet pressure from recycle pump
    m.fs.P1.outlet.pressure.fix(101325)

    # Check degrees of freedom
    print(degrees_of_freedom(m))
    assert degrees_of_freedom(m) == 0

    def scale_variables(m):
        for var in m.fs.component_data_objects(pyo.Var, descend_into=True):
            if "flow_vol" in var.name:
                iscale.set_scaling_factor(var, 1e2)
            if "temperature" in var.name:
                iscale.set_scaling_factor(var, 1e-1)
            if "pressure" in var.name:
                iscale.set_scaling_factor(var, 1e-4)
            if "enth_mol" in var.name:
                iscale.set_scaling_factor(var, 1e-4)
            if "conc_mass_comp[S_O2]" in var.name:
                iscale.set_scaling_factor(var, 1e5)
            if "conc_mass_comp[S_F]" in var.name:
                iscale.set_scaling_factor(var, 1e2)
            if "conc_mass_comp[S_A]" in var.name:
                iscale.set_scaling_factor(var, 1e2)
            if "conc_mass_comp[S_NH4]" in var.name:
                iscale.set_scaling_factor(var, 1e2)
            if "conc_mass_comp[S_NO3]" in var.name:
                iscale.set_scaling_factor(var, 1e5)
            if "conc_mass_comp[S_PO4]" in var.name:
                iscale.set_scaling_factor(var, 1e3)
            if "conc_mass_comp[S_I]" in var.name:
                iscale.set_scaling_factor(var, 1e2)
            if "conc_mass_comp[S_N2]" in var.name:
                iscale.set_scaling_factor(var, 1e2)
            if "conc_mass_comp[X_I]" in var.name:
                iscale.set_scaling_factor(var, 1e0)
            if "conc_mass_comp[X_S]" in var.name:
                iscale.set_scaling_factor(var, 1e1)
            if "conc_mass_comp[X_H]" in var.name:
                iscale.set_scaling_factor(var, 1e0)
            if "conc_mass_comp[X_PAO]" in var.name:
                iscale.set_scaling_factor(var, 1e1)
            if "conc_mass_comp[X_PP]" in var.name:
                iscale.set_scaling_factor(var, 1e2)
            if "conc_mass_comp[X_PHA]" in var.name:
                iscale.set_scaling_factor(var, 1e2)
            if "conc_mass_comp[X_AUT]" in var.name:
                iscale.set_scaling_factor(var, 1e5)
            if "conc_mass_comp[X_MeOH]" in var.name:
                iscale.set_scaling_factor(var, 1e5)
            if "conc_mass_comp[X_MeP]" in var.name:
                iscale.set_scaling_factor(var, 1e5)
            if "conc_mass_comp[X_TSS]" in var.name:
                iscale.set_scaling_factor(var, 1e0)

    scale_variables(m)
    iscale.calculate_scaling_factors(m.fs)

    # Initialize flowsheet
    # Apply sequential decomposition - 1 iteration should suffice
    seq = SequentialDecomposition()
    seq.options.select_tear_method = "heuristic"
    seq.options.tear_method = "Wegstein"
    seq.options.iterLim = 0

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
        "flow_vol": {0: 0.31},
        "conc_mass_comp": {
            (0, "S_O2"): 0.14,
            (0, "S_F"): 0.021,
            (0, "S_A"): 0.014,
            (0, "S_NH4"): 0.016,
            (0, "S_NO3"): 1e-6,
            (0, "S_PO4"): 0.0035,
            (0, "S_I"): 0.030,
            (0, "S_N2"): 0.015,
            (0, "X_I"): 0.050,
            (0, "X_S"): 0.138,
            (0, "X_H"): 0.014,
            (0, "X_PAO"): 1e-6,
            (0, "X_PP"): 1e-6,
            (0, "X_PHA"): 1e-6,
            (0, "X_AUT"): 1e-6,
            (0, "X_MeOH"): 1e-6,
            (0, "X_MeP"): 1e-6,
            (0, "X_TSS"): 0.366,
        },
        "alkalinity": {0: 7.2e-3},
        "temperature": {0: 298.15},
        "pressure": {0: 101325},
    }

    # Pass the tear_guess to the SD tool
    seq.set_guesses_for(m.fs.R1.inlet, tear_guesses)

    def function(unit):
        unit.initialize(outlvl=idaeslog.INFO, optarg={"bound_push": 1e-8})

    seq.run(m, function)

    solver = get_solver(options={"bound_push": 1e-4})
    results = solver.solve(m, tee=False)
    pyo.assert_optimal_termination(results)

    return m, results


if __name__ == "__main__":
    # This method builds and runs a steady state activated sludge
    # flowsheet.
    m, results = build_flowsheet()

    stream_table = create_stream_table_dataframe(
        {
            "Feed": m.fs.FeedWater.outlet,
            "Mix": m.fs.R1.inlet,
            "R1": m.fs.R1.outlet,
            "R2": m.fs.R2.outlet,
            "R3": m.fs.R3.outlet,
            "R4": m.fs.R4.outlet,
            "R5": m.fs.R5.outlet,
        },
        time_point=0,
    )
    print(stream_table_dataframe_to_string(stream_table))
