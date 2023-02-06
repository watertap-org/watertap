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
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
import idaes.logger as idaeslog
from watertap.unit_models.anaerobic_digestor import AD
from watertap.property_models.anaerobic_digestion.adm1_properties import (
    ADM1ParameterBlock,
)
from watertap.property_models.anaerobic_digestion.adm1_properties_vapor import (
    ADM1_vaporParameterBlock,
)
from watertap.property_models.anaerobic_digestion.adm1_reactions import (
    ADM1ReactionParameterBlock,
)


def build_flowsheet():
    m = pyo.ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.props = ADM1ParameterBlock()
    m.fs.props_vap = ADM1_vaporParameterBlock()
    m.fs.rxn_props = ADM1ReactionParameterBlock(property_package=m.fs.props)

    m.fs.R1 = AD(
        liquid_property_package=m.fs.props,
        vapor_property_package=m.fs.props_vap,
        reaction_package=m.fs.rxn_props,
        has_heat_transfer=True,
        has_pressure_change=False,
    )

    # Feed conditions based on manual mass balance of inlet and recycle streams
    m.fs.R1.inlet.flow_vol.fix(170 * pyo.units.m**3 / pyo.units.day)
    m.fs.R1.inlet.temperature.fix(308.15 * pyo.units.K)
    m.fs.R1.inlet.pressure.fix(1 * pyo.units.atm)
    m.fs.R1.inlet.conc_mass_comp[0, "S_su"].fix(1 * pyo.units.mg / pyo.units.liter)
    m.fs.R1.inlet.conc_mass_comp[0, "S_aa"].fix(1 * pyo.units.mg / pyo.units.liter)
    m.fs.R1.inlet.conc_mass_comp[0, "S_fa"].fix(1 * pyo.units.mg / pyo.units.liter)
    m.fs.R1.inlet.conc_mass_comp[0, "S_va"].fix(1 * pyo.units.mg / pyo.units.liter)
    m.fs.R1.inlet.conc_mass_comp[0, "S_bu"].fix(1 * pyo.units.mg / pyo.units.liter)
    m.fs.R1.inlet.conc_mass_comp[0, "S_pro"].fix(1 * pyo.units.mg / pyo.units.liter)
    m.fs.R1.inlet.conc_mass_comp[0, "S_ac"].fix(1 * pyo.units.mg / pyo.units.liter)
    m.fs.R1.inlet.conc_mass_comp[0, "S_h2"].fix(1e-5 * pyo.units.mg / pyo.units.liter)
    m.fs.R1.inlet.conc_mass_comp[0, "S_ch4"].fix(1e-2 * pyo.units.mg / pyo.units.liter)
    m.fs.R1.inlet.conc_mass_comp[0, "S_IC"].fix(
        40 * units.mmol / units.liter * 12 * units.mg / units.mmol
    )
    m.fs.R1.inlet.conc_mass_comp[0, "S_IN"].fix(
        10 * units.mmol / units.liter * 14 * units.mg / units.mmol
    )
    m.fs.R1.inlet.conc_mass_comp[0, "S_I"].fix(20 * pyo.units.mg / pyo.units.liter)

    m.fs.R1.inlet.conc_mass_comp[0, "X_c"].fix(2000 * pyo.units.mg / pyo.units.liter)
    m.fs.R1.inlet.conc_mass_comp[0, "X_ch"].fix(5000 * pyo.units.mg / pyo.units.liter)
    m.fs.R1.inlet.conc_mass_comp[0, "X_pr"].fix(20000 * pyo.units.mg / pyo.units.liter)
    m.fs.R1.inlet.conc_mass_comp[0, "X_li"].fix(5000 * pyo.units.mg / pyo.units.liter)
    m.fs.R1.inlet.conc_mass_comp[0, "X_su"].fix(1 * pyo.units.mg / pyo.units.liter)
    m.fs.R1.inlet.conc_mass_comp[0, "X_aa"].fix(10 * pyo.units.mg / pyo.units.liter)
    m.fs.R1.inlet.conc_mass_comp[0, "X_fa"].fix(10 * pyo.units.mg / pyo.units.liter)
    m.fs.R1.inlet.conc_mass_comp[0, "X_c4"].fix(10 * pyo.units.mg / pyo.units.liter)
    m.fs.R1.inlet.conc_mass_comp[0, "X_pro"].fix(10 * pyo.units.mg / pyo.units.liter)
    m.fs.R1.inlet.conc_mass_comp[0, "X_ac"].fix(10 * pyo.units.mg / pyo.units.liter)
    m.fs.R1.inlet.conc_mass_comp[0, "X_h2"].fix(10 * pyo.units.mg / pyo.units.liter)
    m.fs.R1.inlet.conc_mass_comp[0, "X_I"].fix(25000 * pyo.units.mg / pyo.units.liter)

    m.fs.R1.inlet.cations[0].fix(40 * pyo.units.mmol / pyo.units.liter)
    m.fs.R1.inlet.anions[0].fix(20 * pyo.units.mmol / pyo.units.liter)

    m.fs.R1.volume_liquid.fix(3400 * pyo.units.m**3)
    m.fs.R1.volume_vapor.fix(300 * pyo.units.m**3)

    m.fs.R1.liquid_outlet.temperature.fix(308.15 * pyo.units.K)

    # TO DO: Fix initialization
    m.fs.R1.initialize(outlvl=idaeslog.INFO_HIGH, optarg={"bound_push": 1e-8})

    solver = get_solver(options={"bound_push": 1e-8})

    results = solver.solve(m, tee=True)

    pyo.assert_optimal_termination(results)

    return m, results
