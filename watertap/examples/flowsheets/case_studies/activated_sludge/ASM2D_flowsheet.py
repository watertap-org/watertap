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
Based on flowsheet from:

Flores-Alsina X., Gernaey K.V. and Jeppsson, U. "Benchmarking biological
nutrient removal in wastewater treatment plants: influence of mathematical model
assumptions", 2012, Wat. Sci. Tech., Vol. 65 No. 8, pp. 1496-1505
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
from watertap.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
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

from watertap.core.util.initialization import check_solve

# Set up logger
_log = idaeslog.getLogger(__name__)


def automate_rescale_variables(m):
    for var, sv in iscale.badly_scaled_var_generator(m):
        if iscale.get_scaling_factor(var) is None:
            continue
        sf = iscale.get_scaling_factor(var)
        iscale.set_scaling_factor(var, sf / sv)
        iscale.calculate_scaling_factors(m)


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
    # First reactor (anoxic) - standard CSTR
    m.fs.R2 = CSTR(property_package=m.fs.props, reaction_package=m.fs.rxn_props)
    # Second reactor (anoxic) - standard CSTR
    m.fs.R3 = CSTR(property_package=m.fs.props, reaction_package=m.fs.rxn_props)

    m.fs.R4 = CSTR(property_package=m.fs.props, reaction_package=m.fs.rxn_props)
    # Third reactor (aerobic) - CSTR with injection
    m.fs.R5 = CSTR_Injection(
        property_package=m.fs.props, reaction_package=m.fs.rxn_props
    )
    # Fourth reactor (aerobic) - CSTR with injection
    m.fs.R6 = CSTR_Injection(
        property_package=m.fs.props, reaction_package=m.fs.rxn_props
    )
    # Fifth reactor (aerobic) - CSTR with injection
    m.fs.R7 = CSTR_Injection(
        property_package=m.fs.props, reaction_package=m.fs.rxn_props
    )
    m.fs.SP1 = Separator(
        property_package=m.fs.props, outlet_list=["underflow", "overflow"]
    )
    # Clarifier
    # TODO: Replace with more detailed model when available
    m.fs.CL1 = Separator(
        property_package=m.fs.props,
        outlet_list=["underflow", "effluent"],
        split_basis=SplittingType.componentFlow,
    )
    # Mixing sludge recycle and R5 underflow
    m.fs.MX2 = Mixer(property_package=m.fs.props, inlet_list=["clarifier", "reactor"])
    # Sludge separator
    m.fs.SP2 = Separator(property_package=m.fs.props, outlet_list=["waste", "recycle"])
    # Product Blocks
    m.fs.Treated = Product(property_package=m.fs.props)
    m.fs.Sludge = Product(property_package=m.fs.props)
    # Recycle pressure changer - use a simple isothermal unit for now
    m.fs.P1 = PressureChanger(property_package=m.fs.props)

    # Link units
    m.fs.stream1 = Arc(source=m.fs.FeedWater.outlet, destination=m.fs.MX1.feed_water)
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

    m.fs.stream15 = Arc(source=m.fs.SP2.waste, destination=m.fs.Sludge.inlet)
    m.fs.stream16 = Arc(source=m.fs.SP2.recycle, destination=m.fs.P1.inlet)

    m.fs.stream17 = Arc(source=m.fs.P1.outlet, destination=m.fs.MX1.recycle)
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
    m.fs.FeedWater.conc_mass_comp[0, "S_N2"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_I"].fix(25 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_S"].fix(125 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_H"].fix(30 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_PAO"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_PP"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_PHA"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_AUT"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_MeOH"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_MeP"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_TSS"].fix(180 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.alkalinity.fix(7 * pyo.units.mol / pyo.units.m**3)

    # Reactor sizing
    m.fs.R1.volume.fix(1000 * pyo.units.m**3)
    m.fs.R2.volume.fix(1000 * pyo.units.m**3)
    m.fs.R3.volume.fix(1000 * pyo.units.m**3)
    m.fs.R4.volume.fix(1000 * pyo.units.m**3)
    m.fs.R5.volume.fix(1333 * pyo.units.m**3)
    m.fs.R6.volume.fix(1333 * pyo.units.m**3)
    m.fs.R7.volume.fix(1333 * pyo.units.m**3)

    # Injection rates to Reactions 3, 4 and 5
    for j in m.fs.props.component_list:
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

    # Clarifier
    # TODO: Update once more detailed model available
    m.fs.CL1.split_fraction[0, "effluent", "H2O"].fix(0.49986)
    m.fs.CL1.split_fraction[0, "effluent", "S_A"].fix(0.49986)
    m.fs.CL1.split_fraction[0, "effluent", "S_F"].fix(0.49986)
    m.fs.CL1.split_fraction[0, "effluent", "S_I"].fix(0.49986)
    m.fs.CL1.split_fraction[0, "effluent", "S_N2"].fix(0.49986)
    m.fs.CL1.split_fraction[0, "effluent", "S_NH4"].fix(0.49986)
    m.fs.CL1.split_fraction[0, "effluent", "S_NO3"].fix(0.49986)
    m.fs.CL1.split_fraction[0, "effluent", "S_O2"].fix(0.49986)
    m.fs.CL1.split_fraction[0, "effluent", "S_PO4"].fix(0.49986)
    m.fs.CL1.split_fraction[0, "effluent", "S_ALK"].fix(0.49986)
    m.fs.CL1.split_fraction[0, "effluent", "X_AUT"].fix(0.022117)
    m.fs.CL1.split_fraction[0, "effluent", "X_H"].fix(0.021922)
    m.fs.CL1.split_fraction[0, "effluent", "X_I"].fix(0.021715)
    m.fs.CL1.split_fraction[0, "effluent", "X_MeOH"].fix(0.022)
    m.fs.CL1.split_fraction[0, "effluent", "X_MeP"].fix(0.022)
    m.fs.CL1.split_fraction[0, "effluent", "X_PAO"].fix(0.022)
    m.fs.CL1.split_fraction[0, "effluent", "X_PHA"].fix(0.02147)
    m.fs.CL1.split_fraction[0, "effluent", "X_PP"].fix(0.02144)
    m.fs.CL1.split_fraction[0, "effluent", "X_S"].fix(0.02221)
    m.fs.CL1.split_fraction[0, "effluent", "X_TSS"].fix(0.02194)

    # Sludge purge separator
    m.fs.SP2.split_fraction[:, "recycle"].fix(0.97955)

    # Outlet pressure from recycle pump
    m.fs.P1.outlet.pressure.fix(101325)

    # Check degrees of freedom
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
                iscale.set_scaling_factor(var, 1e-3)
            if "alkalinity" in var.name:
                iscale.set_scaling_factor(var, 1e3)

            if "conc_mass_comp[S_O2]" in var.name:
                iscale.set_scaling_factor(var, 1e3)
            if "conc_mass_comp[S_F]" in var.name:
                iscale.set_scaling_factor(var, 1e3)
            if "conc_mass_comp[S_A]" in var.name:
                iscale.set_scaling_factor(var, 1e3)
            if "conc_mass_comp[S_NH4]" in var.name:
                iscale.set_scaling_factor(var, 1e2)
            if "conc_mass_comp[S_NO3]" in var.name:
                iscale.set_scaling_factor(var, 1e2)
            if "conc_mass_comp[S_PO4]" in var.name:
                iscale.set_scaling_factor(var, 1e3)
            if "conc_mass_comp[S_I]" in var.name:
                iscale.set_scaling_factor(var, 1e2)
            if "conc_mass_comp[S_N2]" in var.name:
                iscale.set_scaling_factor(var, 1e2)
            if "conc_mass_comp[X_I]" in var.name:
                iscale.set_scaling_factor(var, 1e1)
            if "conc_mass_comp[X_S]" in var.name:
                iscale.set_scaling_factor(var, 1e2)
            if "conc_mass_comp[X_H]" in var.name:
                iscale.set_scaling_factor(var, 1e1)
            if "conc_mass_comp[X_PAO]" in var.name:
                iscale.set_scaling_factor(var, 1e1)
            if "conc_mass_comp[X_PP]" in var.name:
                iscale.set_scaling_factor(var, 1e2)
            if "conc_mass_comp[X_PHA]" in var.name:
                iscale.set_scaling_factor(var, 1e2)
            if "conc_mass_comp[X_AUT]" in var.name:
                iscale.set_scaling_factor(var, 1e1)
            if "conc_mass_comp[X_MeOH]" in var.name:
                iscale.set_scaling_factor(var, 1e7)
            if "conc_mass_comp[X_MeP]" in var.name:
                iscale.set_scaling_factor(var, 1e7)
            if "conc_mass_comp[X_TSS]" in var.name:
                iscale.set_scaling_factor(var, 1e0)

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

    # Initialize flowsheet
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

    # Initial guesses for flow into first reactor
    tear_guesses = {
        "flow_vol": {0: 1.07},
        "conc_mass_comp": {
            (0, "S_O2"): 0.001,
            (0, "S_F"): 0.056,
            (0, "S_A"): 0.013,
            (0, "S_NH4"): 0.010,
            (0, "S_NO3"): 0.011,
            (0, "S_PO4"): 0.0024,
            (0, "S_I"): 0.03,
            (0, "S_N2"): 0.017,
            (0, "X_I"): 1.17,
            (0, "X_S"): 0.05,
            (0, "X_H"): 1.172,
            (0, "X_PAO"): 0.200,
            (0, "X_PP"): 0.052,
            (0, "X_PHA"): 0.013,
            (0, "X_AUT"): 0.030,
            (0, "X_MeOH"): 1e-9,
            (0, "X_MeP"): 1e-9,
            (0, "X_TSS"): 2.2,
        },
        "alkalinity": {0: 7e-3},
        "temperature": {0: 298.15},
        "pressure": {0: 101325},
    }

    # Pass the tear_guess to the SD tool
    seq.set_guesses_for(m.fs.R3.inlet, tear_guesses)

    def function(unit):
        unit.initialize(outlvl=idaeslog.INFO, optarg={"bound_push": 1e-8})
        badly_scaled_vars = list(iscale.badly_scaled_var_generator(unit))
        if len(badly_scaled_vars) > 0:
            automate_rescale_variables(unit)

    seq.run(m, function)

    solver = get_solver()
    results = solver.solve(m, tee=False)
    check_solve(results, checkpoint="closing recycle", logger=_log, fail_flag=True)

    # Switch to fixed KLa in R3 and R4 (S_O concentration is controlled in R5)
    m.fs.R5.KLa.fix(240)
    m.fs.R6.KLa.fix(240)
    m.fs.R7.KLa.fix(84)
    m.fs.R5.outlet.conc_mass_comp[:, "S_O2"].unfix()
    m.fs.R6.outlet.conc_mass_comp[:, "S_O2"].unfix()
    m.fs.R7.outlet.conc_mass_comp[:, "S_O2"].unfix()

    # Resolve with controls in place
    results = solver.solve(m, tee=False)

    pyo.assert_optimal_termination(results)
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
            "R6": m.fs.R6.outlet,
            "R7": m.fs.R7.outlet,
        },
        time_point=0,
    )
    print(stream_table_dataframe_to_string(stream_table))
