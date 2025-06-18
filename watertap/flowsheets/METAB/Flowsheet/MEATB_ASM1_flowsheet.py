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
This flowsheet is WWTP model including METAB system and BSM2
In the BSM2, the ADM1 and ASM 1 is used

"""
__author__ = "Maojian Wang"

import pyomo.environ as pyo
import re
from pyomo.network import Arc, SequentialDecomposition

# from watertap.unit_models.anaerobic_digester import AD
# from watertap.unit_models.thickener import Thickener
# from watertap.unit_models.dewatering import DewateringUnit
# from watertap.unit_models.cstr import CSTR
# from watertap.unit_models.clarifier import Clarifier

# from watertap.unit_models.translators.translator_asm1_adm1 import Translator_ASM1_ADM1
# from watertap.unit_models.translators.translator_adm1_asm1 import Translator_ADM1_ASM1

import idaes.logger as idaeslog
from watertap.core.solvers import get_solver
import idaes.core.util.scaling as iscale
from watertap.property_models.unit_specific.anaerobic_digestion.adm1_properties import (
    ADM1ParameterBlock,
)
from watertap.property_models.unit_specific.anaerobic_digestion.adm1_reactions import (
    ADM1ReactionParameterBlock,
)
from idaes.models.unit_models.mixer import MomentumMixingType
from idaes.models.unit_models.separator import SplittingType
from watertap.property_models.unit_specific.anaerobic_digestion.adm1_properties_vapor import (
    ADM1_vaporParameterBlock,
)

from idaes.core import FlowsheetBlock, UnitModelCostingBlock, UnitModelBlockData
from idaes.models.unit_models import (
    Feed,
    Mixer,
    Separator,
    PressureChanger,
    Product,
)

# from watertap.unit_models.aeration_tank import AerationTank, ElectricityConsumption
from watertap.property_models.unit_specific.activated_sludge.asm1_properties import (
    ASM1ParameterBlock,
)
from watertap.property_models.unit_specific.activated_sludge.asm1_reactions import (
    ASM1ReactionParameterBlock,
)

# from watertap.core.util.initialization import assert_degrees_of_freedom
# from watertap.costing import WaterTAPCosting
# from watertap.costing.unit_models.clarifier import (
#     cost_circular_clarifier,
#     cost_primary_clarifier,
# )
from pyomo.util.check_units import assert_units_consistent
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import propagate_state
from idaes.core.util.tables import (
    create_stream_table_dataframe,
    stream_table_dataframe_to_string,
)
from idaes.core.util import DiagnosticsToolbox
from pyomo.environ import units
import pandas as pd

import ASM1_flowsheet as asm1
import model_connector as metab


def model_checker(model, solver_info):
    dof = degrees_of_freedom(model)
    assert dof == 0, f"Error: Degrees of freedom is {dof}, but it should be 0."
    print("DOF check passed")
    for var in m.fs.component_data_objects(pyo.Var, descend_into=True):
        if "flow_vol" in var.name:
            iscale.set_scaling_factor(var, 1e1)
        if "temperature" in var.name:
            iscale.set_scaling_factor(var, 1e-1)
        if "pressure" in var.name:
            iscale.set_scaling_factor(var, 1e-4)
        if "conc_mass_comp" in var.name:
            iscale.set_scaling_factor(var, 1e3)

    iscale.calculate_scaling_factors(m)
    solver = get_solver()
    results = solver.solve(model, tee=solver_info)
    print(results)
    # pyo.assert_optimal_termination(results)

    return print(" model passed the check")


def build_asm1():
    m, results = asm1.build_flowsheet()

    return m


def add_metab(m):
    m.fs.props_asm1 = ASM1ParameterBlock()
    m.fs.props_adm1 = ADM1ParameterBlock()
    m.fs.adm1_rxn_props = ADM1ReactionParameterBlock(property_package=m.fs.props_adm1)
    m = metab.build_metab(m, mass=True)

    return m


def set_port_values(m):
    # TODO feedwater block no inlet, need to fix
    m.fs.metab_to_asm1 = Arc(
        source=m.fs.metab_effluent.outlet, destination=m.fs.FeedWater.outlet
    )
    propagate_state(m.fs.metab_to_asm1)
    # set_port_values(m.fs.FeedWater.outlet,m.fs.metab_effluent.outlet)
    return m


def reset_influent(m):
    # TODO: Need to update
    m.fs.FeedWater.flow_vol.fix(18446 * pyo.units.m**3 / pyo.units.day)
    # TODO: Need to update
    m.fs.FeedWater.temperature.fix(298.15 * pyo.units.K)
    m.fs.FeedWater.pressure.fix(1 * pyo.units.atm)
    m.fs.FeedWater.conc_mass_comp[0, "S_I"].fix(86.312 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_S"].fix(
        0.039661 * 1000 * pyo.units.g / pyo.units.m**3
    )
    m.fs.FeedWater.conc_mass_comp[0, "X_I"].fix(
        0.11152 * 1000 * pyo.units.g / pyo.units.m**3
    )
    m.fs.FeedWater.conc_mass_comp[0, "X_S"].fix(
        0.74911 * 1000 * pyo.units.g / pyo.units.m**3
    )
    m.fs.FeedWater.conc_mass_comp[0, "X_BH"].fix(0 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_BA"].fix(0 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_P"].fix(0 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_O"].fix(0 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_NO"].fix(0 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_NH"].fix(
        0.20578 * 1000 * pyo.units.g / pyo.units.m**3
    )
    m.fs.FeedWater.conc_mass_comp[0, "S_ND"].fix(
        0.0052525 * 1000 * pyo.units.g / pyo.units.m**3
    )
    m.fs.FeedWater.conc_mass_comp[0, "X_ND"].fix(
        0.058108 * 1000 * pyo.units.g / pyo.units.m**3
    )
    m.fs.FeedWater.alkalinity.fix(7 * pyo.units.mol / pyo.units.m**3)  # 46.6

    return m


def report_st(m):

    stream_table = create_stream_table_dataframe(
        {
            "Feed": m.fs.FeedWater.outlet,
            "R1": m.fs.R1.outlet,
            "R2": m.fs.R2.outlet,
            "R3": m.fs.R3.outlet,
            "R4": m.fs.R4.outlet,
            "R5": m.fs.R5.outlet,
            "from metab": m.fs.metab_effluent.inlet,
            "Influent": m.fs.metab_effluent.outlet,
        },
        time_point=0,
    )
    print(stream_table_dataframe_to_string(stream_table))


if __name__ == "__main__":

    # TODO: NEED to update parameters
    m = build_asm1()
    m = add_metab(m)
    model_checker(m, True)
    report_st(m)
    # m = set_port_values(m)
    m = reset_influent(m)
    model_checker(m, True)
    report_st(m)
