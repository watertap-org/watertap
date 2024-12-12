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
__author__ = "Alejandro Garciadiego, Adam Atia"

import pyomo.environ as pyo
from pyomo.environ import (
    units,
)
from idaes.core import FlowsheetBlock
import idaes.core.util.scaling as iscale
from watertap.core.solvers import get_solver
import idaes.logger as idaeslog
from watertap.unit_models.anaerobic_digester import AD
from watertap.property_models.unit_specific.anaerobic_digestion.adm1_properties import (
    ADM1ParameterBlock,
)
from watertap.property_models.unit_specific.anaerobic_digestion.adm1_properties_vapor import (
    ADM1_vaporParameterBlock,
)
from watertap.property_models.unit_specific.anaerobic_digestion.adm1_reactions import (
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
    m.fs.R1.inlet.conc_mass_comp[0, "X_su"].fix(0.0 * pyo.units.mg / pyo.units.liter)
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

    iscale.calculate_scaling_factors(m)

    # TO DO: Fix initialization
    m.fs.R1.initialize(outlvl=idaeslog.INFO_HIGH)

    solver = get_solver()

    results = solver.solve(m, tee=True)

    pyo.assert_optimal_termination(results)

    return m, results
