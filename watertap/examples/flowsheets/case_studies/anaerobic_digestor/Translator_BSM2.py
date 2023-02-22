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
__author__ = "Alejandro Garciadiego, Adam Atia"

import pyomo.environ as pyo
from pyomo.environ import (
    units,
)

from pyomo.environ import (
    Reference,
    Var,
    value,
    Constraint,
    Param,
    units as pyunits,
    check_optimal_termination,
    exp,
    Set,
)
from pyomo.network import Arc, SequentialDecomposition
from idaes.core.util.model_statistics import degrees_of_freedom, large_residuals_set
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
import idaes.logger as idaeslog
from idaes.models.unit_models import (
    Feed,
    Product,
    Translator,
)
from watertap.unit_models.translator_adm1_asm1 import Translator_ADM1_ASM1
from watertap.property_models.anaerobic_digestion.adm1_properties import (
    ADM1ParameterBlock,
)

from watertap.property_models.activated_sludge.asm1_properties import (
    ASM1ParameterBlock,
)

from watertap.property_models.anaerobic_digestion.adm1_reactions import (
    ADM1ReactionParameterBlock,
)

m = pyo.ConcreteModel()

m.fs = FlowsheetBlock(dynamic=False)

m.fs.props_ASM1 = ASM1ParameterBlock()
m.fs.props_ADM1 = ADM1ParameterBlock()
m.fs.ADM1_rxn_props = ADM1ReactionParameterBlock(property_package=m.fs.props_ADM1)

m.fs.ADM1feed = Feed(property_package=m.fs.props_ADM1)

m.fs.ASM1Treated = Product(property_package=m.fs.props_ASM1)

m.fs.T101 = Translator_ADM1_ASM1(
    inlet_property_package=m.fs.props_ADM1,
    outlet_property_package=m.fs.props_ASM1,
    has_phase_equilibrium=False,
    outlet_state_defined=True,
)

m.fs.stream1 = Arc(source=m.fs.ADM1feed.outlet, destination=m.fs.T101.inlet)
m.fs.stream2 = Arc(source=m.fs.T101.outlet, destination=m.fs.ASM1Treated.inlet)

pyo.TransformationFactory("network.expand_arcs").apply_to(m)

# Feed conditions based on manual mass balance of inlet and recycle streams
m.fs.ADM1feed.flow_vol.fix(170 * pyo.units.m**3 / pyo.units.day)
m.fs.ADM1feed.temperature.fix(308.15 * pyo.units.K)
m.fs.ADM1feed.pressure.fix(1 * pyo.units.atm)
m.fs.ADM1feed.conc_mass_comp[0, "S_su"].fix(1 * pyo.units.mg / pyo.units.liter)
m.fs.ADM1feed.conc_mass_comp[0, "S_aa"].fix(1 * pyo.units.mg / pyo.units.liter)
m.fs.ADM1feed.conc_mass_comp[0, "S_fa"].fix(1 * pyo.units.mg / pyo.units.liter)
m.fs.ADM1feed.conc_mass_comp[0, "S_va"].fix(1 * pyo.units.mg / pyo.units.liter)
m.fs.ADM1feed.conc_mass_comp[0, "S_bu"].fix(1 * pyo.units.mg / pyo.units.liter)
m.fs.ADM1feed.conc_mass_comp[0, "S_pro"].fix(1 * pyo.units.mg / pyo.units.liter)
m.fs.ADM1feed.conc_mass_comp[0, "S_ac"].fix(1 * pyo.units.mg / pyo.units.liter)
m.fs.ADM1feed.conc_mass_comp[0, "S_h2"].fix(1e-5 * pyo.units.mg / pyo.units.liter)
m.fs.ADM1feed.conc_mass_comp[0, "S_ch4"].fix(1e-2 * pyo.units.mg / pyo.units.liter)
m.fs.ADM1feed.conc_mass_comp[0, "S_IC"].fix(
    40 * units.mmol / units.liter * 12 * units.mg / units.mmol
)
m.fs.ADM1feed.conc_mass_comp[0, "S_IN"].fix(
    10 * units.mmol / units.liter * 14 * units.mg / units.mmol
)
m.fs.ADM1feed.conc_mass_comp[0, "S_I"].fix(20 * pyo.units.mg / pyo.units.liter)

m.fs.ADM1feed.conc_mass_comp[0, "X_c"].fix(2000 * pyo.units.mg / pyo.units.liter)
m.fs.ADM1feed.conc_mass_comp[0, "X_ch"].fix(5000 * pyo.units.mg / pyo.units.liter)
m.fs.ADM1feed.conc_mass_comp[0, "X_pr"].fix(20000 * pyo.units.mg / pyo.units.liter)
m.fs.ADM1feed.conc_mass_comp[0, "X_li"].fix(5000 * pyo.units.mg / pyo.units.liter)
m.fs.ADM1feed.conc_mass_comp[0, "X_su"].fix(1 * pyo.units.mg / pyo.units.liter)
m.fs.ADM1feed.conc_mass_comp[0, "X_aa"].fix(10 * pyo.units.mg / pyo.units.liter)
m.fs.ADM1feed.conc_mass_comp[0, "X_fa"].fix(10 * pyo.units.mg / pyo.units.liter)
m.fs.ADM1feed.conc_mass_comp[0, "X_c4"].fix(10 * pyo.units.mg / pyo.units.liter)
m.fs.ADM1feed.conc_mass_comp[0, "X_pro"].fix(10 * pyo.units.mg / pyo.units.liter)
m.fs.ADM1feed.conc_mass_comp[0, "X_ac"].fix(10 * pyo.units.mg / pyo.units.liter)
m.fs.ADM1feed.conc_mass_comp[0, "X_h2"].fix(10 * pyo.units.mg / pyo.units.liter)
m.fs.ADM1feed.conc_mass_comp[0, "X_I"].fix(25000 * pyo.units.mg / pyo.units.liter)

m.fs.ADM1feed.cations[0].fix(40 * pyo.units.mmol / pyo.units.liter)
m.fs.ADM1feed.anions[0].fix(20 * pyo.units.mmol / pyo.units.liter)


# Check degrees of freedom
assert degrees_of_freedom(m) == 0


# Initialize flowsheet
# Apply sequential decomposition - 1 iteration should suffice
seq = SequentialDecomposition()
seq.options.select_tear_method = "heuristic"
seq.options.tear_method = "Wegstein"
seq.options.iterLim = 1

G = seq.create_graph(m)


def function(unit):
    unit.initialize()


seq.run(m, function)

# Solve overall flowsheet to close recycle loop
solver = get_solver(options={"bound_push": 1e-8})

results = solver.solve(m, tee=True)

pyo.assert_optimal_termination(results)
