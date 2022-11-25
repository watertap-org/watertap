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

from idaes.models.unit_models.separator import SplittingType
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


def automate_rescale_variables(m):
    for var, sv in iscale.badly_scaled_var_generator(m):
        if iscale.get_scaling_factor(var) is None:
            continue
        sf = iscale.get_scaling_factor(var)
        iscale.set_scaling_factor(var, sf / sv)
        iscale.calculate_scaling_factors(m)


m = pyo.ConcreteModel()

m.fs = FlowsheetBlock(dynamic=False)

m.fs.props = ADM1ParameterBlock()
m.fs.props_vap = ADM1_vaporParameterBlock()
m.fs.rxn_props = ADM1ReactionParameterBlock(property_package=m.fs.props)

m.fs.R1 = AD(
    liquid_property_package=m.fs.props,
    vapor_property_package=m.fs.props_vap,
    reaction_package=m.fs.rxn_props,
)

# Feed conditions based on manual mass balance of inlet and recycle streams
m.fs.R1.inlet.flow_vol.fix(92230 * pyo.units.m**3 / pyo.units.day)
m.fs.R1.inlet.temperature.fix(298.15 * pyo.units.K)
m.fs.R1.inlet.pressure.fix(1 * pyo.units.atm)
m.fs.R1.inlet.conc_mass_comp[0, "S_su"].fix(10 * pyo.units.mg / pyo.units.liter)
m.fs.R1.inlet.conc_mass_comp[0, "S_aa"].fix(1 * pyo.units.mg / pyo.units.liter)
m.fs.R1.inlet.conc_mass_comp[0, "S_fa"].fix(1 * pyo.units.mg / pyo.units.liter)
m.fs.R1.inlet.conc_mass_comp[0, "S_va"].fix(1 * pyo.units.mg / pyo.units.liter)
m.fs.R1.inlet.conc_mass_comp[0, "S_bu"].fix(1 * pyo.units.mg / pyo.units.liter)
m.fs.R1.inlet.conc_mass_comp[0, "S_pro"].fix(1 * pyo.units.mg / pyo.units.liter)
m.fs.R1.inlet.conc_mass_comp[0, "S_ac"].fix(1 * pyo.units.mg / pyo.units.liter)
m.fs.R1.inlet.conc_mass_comp[0, "S_h2"].fix(1e-5 * pyo.units.mg / pyo.units.liter)
m.fs.R1.inlet.conc_mass_comp[0, "S_ch4"].fix(1e-3 * pyo.units.mg / pyo.units.liter)
m.fs.R1.inlet.conc_mass_comp[0, "S_IC"].fix(40 * pyo.units.mg / pyo.units.liter * 14)
m.fs.R1.inlet.conc_mass_comp[0, "S_IN"].fix(10 * pyo.units.mg / pyo.units.liter * 12)
m.fs.R1.inlet.conc_mass_comp[0, "S_I"].fix(20 * pyo.units.mg / pyo.units.liter)

m.fs.R1.inlet.conc_mass_comp[0, "X_c"].fix(2000 * pyo.units.mg / pyo.units.liter)
m.fs.R1.inlet.conc_mass_comp[0, "X_ch"].fix(5000 * pyo.units.mg / pyo.units.liter)
m.fs.R1.inlet.conc_mass_comp[0, "X_pr"].fix(20000 * pyo.units.mg / pyo.units.liter)
m.fs.R1.inlet.conc_mass_comp[0, "X_li"].fix(5000 * pyo.units.mg / pyo.units.liter)
m.fs.R1.inlet.conc_mass_comp[0, "X_su"].fix(1e-6 * pyo.units.mg / pyo.units.liter)
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

# TO DO: Include in reactor
m.fs.R1.k_p = pyo.Param(
    initialize=5e4,
    mutable=True,
    units=pyo.units.m**3 / pyo.units.day / pyo.units.Pa,
    doc="Component mass concentrations",
)


@m.fs.R1.Constraint(m.fs.time, doc="Pressure")
def outlet_P(self, t):
    return m.fs.R1.vapor_outlet.pressure[t] == (
        m.fs.R1.vapor_phase[0].p_w_sat
        + sum(m.fs.R1.vapor_phase[0].p_sat[j] for j in m.fs.props_vap.solute_set)
        - m.fs.props_vap.pressure_ref
    )


m.fs.R1.KH_co2 = pyo.Param(
    initialize=0.0271 / 1e6,
    units=pyo.units.kmol / pyo.units.m**3 * pyo.units.Pa**-1,
    mutable=True,
    doc="KH_co2",
)
m.fs.R1.KH_ch4 = pyo.Param(
    initialize=0.00116 / 1e6,
    units=pyo.units.kmol / pyo.units.m**3 * pyo.units.Pa**-1,
    mutable=True,
    doc="KH_ch4",
)
m.fs.R1.KH_h2 = pyo.Param(
    initialize=7.38e-4 / 1e6,
    units=pyo.units.kmol / pyo.units.m**3 * pyo.units.Pa**-1,
    mutable=True,
    doc="KH_h2",
)
m.fs.R1.K_La = pyo.Param(
    initialize=200,
    units=pyo.units.day**-1,
    mutable=True,
    doc="K_La",
)


@m.fs.R1.Constraint(m.fs.time, doc="Pressure")
def Sh2_conc(self, t):
    return m.fs.R1.liquid_phase.mass_transfer_term[t, "Liq", "S_h2"] == (
        pyo.units.convert(m.fs.R1.K_La, to_units=1 / pyo.units.s)
        * (
            m.fs.R1.liquid_phase.properties_out[t].conc_mass_comp["S_h2"]
            / (2 * pyo.units.kg / pyo.units.kmole)
            - 16 * m.fs.R1.KH_h2 * m.fs.R1.vapor_phase[t].p_sat["S_h2"]
        )
        * m.fs.R1.volume_liquid[t]
    ) * (2 * pyo.units.kg / pyo.units.kmole)


@m.fs.R1.Constraint(m.fs.time, doc="Pressure")
def Sch4_conc(self, t):
    return m.fs.R1.liquid_phase.mass_transfer_term[t, "Liq", "S_ch4"] == (
        pyo.units.convert(m.fs.R1.K_La, to_units=1 / pyo.units.s)
        * (
            m.fs.R1.liquid_phase.properties_out[t].conc_mass_comp["S_ch4"]
            / (16 * pyo.units.kg / pyo.units.kmole)
            - 64 * m.fs.R1.KH_ch4 * m.fs.R1.vapor_phase[t].p_sat["S_ch4"]
        )
        * m.fs.R1.volume_liquid[t]
    ) * (16 * pyo.units.kg / pyo.units.kmole)


@m.fs.R1.Constraint(m.fs.time, doc="Pressure")
def Sco2_conc(self, t):
    return m.fs.R1.liquid_phase.mass_transfer_term_j[t] == (
        pyo.units.convert(m.fs.R1.K_La, to_units=1 / pyo.units.s)
        * (
            m.fs.R1.liquid_phase.reactions[t].conc_mol_co2
            - m.fs.R1.KH_co2 * m.fs.R1.vapor_phase[t].p_sat["S_co2"]
        )
        * m.fs.R1.volume_liquid[t]
    ) * (44 * pyo.units.kg / pyo.units.kmole)


@m.fs.R1.Constraint(m.fs.time, doc="Pressure")
def flow_vol_vap(self, t):
    return m.fs.R1.vapor_outlet.flow_vol[t] == (
        pyo.units.convert(
            m.fs.R1.k_p * (m.fs.R1.vapor_outlet.pressure[t] - 100000 * pyo.units.Pa),
            to_units=pyo.units.m**3 / pyo.units.s,
        )
    )


# TO DO: Better scaling
def scale_variables(m):
    for var in m.fs.component_data_objects(pyo.Var, descend_into=True):
        if "flow_vol" in var.name:
            iscale.set_scaling_factor(var, 1e2)
        if "temperature" in var.name:
            iscale.set_scaling_factor(var, 1e-1)
        if "pressure" in var.name:
            iscale.set_scaling_factor(var, 1e-3)
        if "enth_mol" in var.name:
            iscale.set_scaling_factor(var, 1e-3)
        if "heat" in var.name:
            iscale.set_scaling_factor(var, 1e-5)
        if "conc_mass_comp" in var.name:
            iscale.set_scaling_factor(var, 1e2)

        # if "p_sat" in var.name:
        #     iscale.set_scaling_factor(var, 1e-3)
        if "p_w_sat" in var.name:
            iscale.set_scaling_factor(var, 1e-3)
        if "conc_mol_" in var.name:
            iscale.set_scaling_factor(var, 1e5)
        if "reaction_rate" in var.name:
            iscale.set_scaling_factor(var, 1e3)
        if "S_H" in var.name:
            iscale.set_scaling_factor(var, 1e7)
        if "S_OH" in var.name:
            iscale.set_scaling_factor(var, 1e7)
        if "rate_reaction_extent" in var.name:
            iscale.set_scaling_factor(var, 1e2)
        if "volume" in var.name:
            iscale.set_scaling_factor(var, 1e-3)


iscale.set_scaling_factor(m.fs.R1.liquid_phase.reactions[0.0].rate_expression, 1e3)
iscale.set_scaling_factor(m.fs.R1.ad_performance_eqn, 1e3)

iscale.set_scaling_factor(
    m.fs.R1.liquid_phase.rate_reaction_stoichiometry_constraint, 1e3
)
iscale.set_scaling_factor(m.fs.R1.liquid_phase.material_balances, 1e3)

# Apply scaling
scale_variables(m)

iscale.calculate_scaling_factors(m)

# TO DO: Fix initialization
m.fs.R1.initialize(outlvl=idaeslog.INFO_HIGH, optarg={"bound_push": 1e-8})

# TO DO: Better scaling
badly_scaled_vars = list(iscale.badly_scaled_var_generator(m.fs.R1))

for var, sv in iscale.badly_scaled_var_generator(m):
    print(var, ",", sv)
if len(badly_scaled_vars) > 0:
    automate_rescale_variables(m.fs.R1)

solver = get_solver(options={"bound_push": 1e-8})
results = solver.solve(m, tee=True)
