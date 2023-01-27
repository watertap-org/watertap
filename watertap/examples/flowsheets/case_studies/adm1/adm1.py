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
__author__ = "Alejandro Garciadiego, Adam Atia"

import pyomo.environ as pyo
from pyomo.environ import (
    units,
)
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom, large_residuals_set
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale
from idaes.core.util.tables import (
    create_stream_table_dataframe,
    stream_table_dataframe_to_string,
)
from pyomo.util.check_units import assert_units_consistent
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
from idaes.core.util.model_diagnostics import DegeneracyHunter

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
m.fs.R1.inlet.conc_mass_comp[0, "S_su"].fix(10 * pyo.units.mg / pyo.units.liter)
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
m.fs.R1.inlet.conc_mass_comp[0, "X_su"].fix(1e-8 * pyo.units.mg / pyo.units.liter)
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
m.fs.R1.liquid_phase.enthalpy_transfer.fix(-1327.08315)
m.fs.R1.liquid_outlet.temperature.fix(308.15 * pyo.units.K)
m.fs.R1.liquid_phase.enthalpy_transfer.unfix()

m.fs.R1.liquid_phase.properties_out[0.0].flow_vol.fix(167.25 * units.m**3 / units.day)

m.fs.R1.liquid_phase.properties_out[0.0].pressure.fix(1 * units.atm)
m.fs.R1.liquid_phase.properties_out[0.0].cations.fix(40.597 * units.mmol / units.liter)
m.fs.R1.liquid_phase.properties_out[0.0].anions.fix(20.298 * units.mmol / units.liter)

m.fs.R1.liquid_phase.properties_out[0.0].flow_vol.unfix()

m.fs.R1.liquid_phase.properties_out[0.0].pressure.unfix()
m.fs.R1.liquid_phase.properties_out[0.0].cations.unfix()
m.fs.R1.liquid_phase.properties_out[0.0].anions.unfix()

m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["S_su"].fix(
    12.133 * units.mg / units.liter
)
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["S_aa"].fix(
    5.3940 * units.mg / units.liter
)
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["S_fa"].fix(
    100.09 * units.mg / units.liter
)
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["S_va"].fix(
    11.798 * units.mg / units.liter
)
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["S_bu"].fix(
    13.44 * units.mg / units.liter
)
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["S_pro"].fix(
    16.019243 * units.mg / units.liter
)
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["S_ac"].fix(
    200.57 * units.mg / units.liter
)
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["S_h2"].fix(
    0.00023594 * units.mg / units.liter
)
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["S_ch4"].fix(
    55.088 * units.mg / units.liter
)
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["S_IC"].fix(
    154.9566 * units.mmol / units.liter * 14 * units.mg / units.mmol
)
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["S_IN"].fix(
    132.1735 * units.mmol / units.liter * 12 * units.mg / units.mmol
)
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["S_I"].fix(
    333.6036 * units.mg / units.liter
)

m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["X_c"].fix(
    313.7 * units.mg / units.liter
)
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["X_ch"].fix(
    28.364 * units.mg / units.liter
)
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["X_pr"].fix(
    104.105 * units.mg / units.liter
)
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["X_li"].fix(
    29.923 * units.mg / units.liter
)
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["X_su"].fix(
    426.43711 * units.mg / units.liter
)
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["X_aa"].fix(
    1198.560 * units.mg / units.liter
)
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["X_fa"].fix(
    246.6627 * units.mg / units.liter
)
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["X_c4"].fix(
    438.67689 * units.mg / units.liter
)
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["X_pro"].fix(
    139.355 * units.mg / units.liter
)
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["X_ac"].fix(
    771.91434 * units.mg / units.liter
)
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["X_h2"].fix(
    321.75439 * units.mg / units.liter
)
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["X_I"].fix(
    26038 * units.mg / units.liter
)
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["S_su"].unfix()
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["S_aa"].unfix()
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["S_fa"].unfix()
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["S_va"].unfix()
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["S_bu"].unfix()
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["S_pro"].unfix()
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["S_ac"].unfix()
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["S_h2"].unfix()
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["S_ch4"].unfix()
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["S_IC"].unfix()
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["S_IN"].unfix()
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["S_I"].unfix()

m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["X_c"].unfix()
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["X_ch"].unfix()
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["X_pr"].unfix()
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["X_li"].unfix()
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["X_su"].unfix()
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["X_aa"].unfix()
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["X_fa"].unfix()
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["X_c4"].unfix()
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["X_pro"].unfix()
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["X_ac"].unfix()
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["X_h2"].unfix()
m.fs.R1.liquid_phase.properties_out[0.0].conc_mass_comp["X_I"].unfix()

vol = 1
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R1"].fix(1.786e-6 * vol)
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R2"].fix(3.235e-6 * vol)
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R3"].fix(1.187e-5 * vol)
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R4"].fix(3.412e-6 * vol)
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R5"].fix(3.404e-6 * vol)
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R6"].fix(1.18693e-5 * vol)
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R7"].fix(3.185e-6 * vol)
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R8"].fix(2.505e-6 * vol)
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R9"].fix(3.230e-6 * vol)
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R10"].fix(2.636e-6 * vol)
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R11"].fix(1.220e-5 * vol)
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R12"].fix(4.184e-6 * vol)
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R13"].fix(9.726e-8 * vol)
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R14"].fix(2.730e-7 * vol)
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R15"].fix(5.626e-8 * vol)
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R16"].fix(9.998e-8 * vol)
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R17"].fix(3.178e-8 * vol)
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R18"].fix(1.761e-7 * vol)
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R19"].fix(7.338e-8 * vol)

m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R1"].unfix()
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R2"].unfix()
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R3"].unfix()
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R4"].unfix()
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R5"].unfix()
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R6"].unfix()
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R7"].unfix()
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R8"].unfix()
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R9"].unfix()
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R10"].unfix()
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R11"].unfix()
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R12"].unfix()
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R13"].unfix()
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R14"].unfix()
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R15"].unfix()
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R16"].unfix()
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R17"].unfix()
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R18"].unfix()
m.fs.R1.liquid_phase.reactions[0.0].reaction_rate["R19"].unfix()

vol2 = 3400
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R1"].fix(1.786e-6 * vol2)
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R2"].fix(3.235e-6 * vol2)
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R3"].fix(1.187e-5 * vol2)
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R4"].fix(3.412e-6 * vol2)
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R5"].fix(3.404e-6 * vol2)
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R6"].fix(1.187e-5 * vol2)
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R7"].fix(3.185e-6 * vol2)
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R8"].fix(2.505e-6 * vol2)
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R9"].fix(3.230e-6 * vol2)
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R10"].fix(2.636e-6 * vol2)
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R11"].fix(1.220e-5 * vol2)
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R12"].fix(4.184e-6 * vol2)
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R13"].fix(9.726e-8 * vol2)
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R14"].fix(2.730e-7 * vol2)
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R15"].fix(5.626e-8 * vol2)
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R16"].fix(9.998e-8 * vol2)
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R17"].fix(3.178e-8 * vol2)
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R18"].fix(1.761e-7 * vol2)
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R19"].fix(7.338e-8 * vol2)

m.fs.R1.liquid_phase.rate_reaction_extent[0, "R1"].unfix()
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R2"].unfix()
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R3"].unfix()
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R4"].unfix()
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R5"].unfix()
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R6"].unfix()
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R7"].unfix()
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R8"].unfix()
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R9"].unfix()
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R10"].unfix()
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R11"].unfix()
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R12"].unfix()
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R13"].unfix()
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R14"].unfix()
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R15"].unfix()
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R16"].unfix()
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R17"].unfix()
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R18"].unfix()
m.fs.R1.liquid_phase.rate_reaction_extent[0, "R19"].unfix()

m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "S_su"].fix(3.846e-6)
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "S_aa"].fix(8.489e-6)
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "S_fa"].fix(0.000192078)
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "S_va"].fix(2.0905e-5)
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "S_bu"].fix(2.4104e-5)
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "S_pro"].fix(2.9088e-5)
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "S_ac"].fix(0.000393281)
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "S_h2"].fix(4.445e-10)
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "S_ch4"].fix(0.0557198)
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "S_IC"].fix(0.002660449)
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "S_IN"].fix(0.002838759)
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "S_I"].fix(0.0006073)
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "X_c"].fix(-0.003328)
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "X_ch"].fix(-0.009782)
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "X_pr"].fix(-0.03915)
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "X_li"].fix(-0.00977995)
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "X_su"].fix(0.000826)
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "X_aa"].fix(0.002300)
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "X_fa"].fix(0.0004585)
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "X_c4"].fix(0.0008301)
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "X_pro"].fix(0.00025046)
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "X_ac"].fix(0.001476)
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "X_h2"].fix(0.000604)
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "X_I"].fix(0.00121)

m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "S_su"].unfix()
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "S_aa"].unfix()
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "S_fa"].unfix()
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "S_va"].unfix()
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "S_bu"].unfix()
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "S_pro"].unfix()
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "S_ac"].unfix()
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "S_h2"].unfix()
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "S_ch4"].unfix()
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "S_IC"].unfix()
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "S_IN"].unfix()
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "S_I"].unfix()
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "X_c"].unfix()
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "X_ch"].unfix()
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "X_pr"].unfix()
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "X_li"].unfix()
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "X_su"].unfix()
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "X_aa"].unfix()
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "X_fa"].unfix()
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "X_c4"].unfix()
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "X_pro"].unfix()
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "X_ac"].unfix()
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "X_h2"].unfix()
m.fs.R1.liquid_phase.rate_reaction_generation[0, "Liq", "X_I"].unfix()

m.fs.R1.liquid_phase.mass_transfer_term[0, "Liq", "S_h2"].fix(-3.50341e-07)
m.fs.R1.liquid_phase.mass_transfer_term[0, "Liq", "S_ch4"].fix(-0.05561126)
m.fs.R1.liquid_phase.mass_transfer_term[0, "Liq", "S_IC"].fix(-0.000482354)
m.fs.R1.liquid_phase.mass_transfer_term[0, "Liq", "H2O"].fix(-0.031733218)

m.fs.R1.liquid_phase.mass_transfer_term[0, "Liq", "S_h2"].unfix()
m.fs.R1.liquid_phase.mass_transfer_term[0, "Liq", "S_ch4"].unfix()
m.fs.R1.liquid_phase.mass_transfer_term[0, "Liq", "S_IC"].unfix()
m.fs.R1.liquid_phase.mass_transfer_term[0, "Liq", "H2O"].unfix()

# TO DO: Fix initialization
m.fs.R1.initialize(outlvl=idaeslog.INFO_HIGH, optarg={"bound_push": 1e-8})
#
# solver = get_solver(options={"bound_push": 1e-8})
# #
# results = solver.solve(m, tee=True)

m.display()
print(large_residuals_set(m))
