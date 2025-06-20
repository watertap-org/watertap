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
Translate brewery influnt condistions from the ADM1 property to ASM1 porperty

"""
__author__ = "Maojian Wang"

import pyomo.environ as pyo
import re
from pyomo.network import Arc, SequentialDecomposition
from watertap.unit_models.anaerobic_digester import AD
from watertap.unit_models.thickener import Thickener
from watertap.unit_models.dewatering import DewateringUnit
from watertap.unit_models.cstr import CSTR
from watertap.unit_models.clarifier import Clarifier

from watertap.unit_models.translators.translator_asm1_adm1 import Translator_ASM1_ADM1
from watertap.unit_models.translators.translator_adm1_asm1 import Translator_ADM1_ASM1

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

from watertap.unit_models.aeration_tank import AerationTank, ElectricityConsumption
from watertap.property_models.unit_specific.activated_sludge.asm1_properties import (
    ASM1ParameterBlock,
)
from watertap.property_models.unit_specific.activated_sludge.asm1_reactions import (
    ASM1ReactionParameterBlock,
)
from watertap.core.util.initialization import assert_degrees_of_freedom
from watertap.costing import WaterTAPCosting
from watertap.costing.unit_models.clarifier import (
    cost_circular_clarifier,
    cost_primary_clarifier,
)
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


def get_brewwry_influent(m):

    df = pd.read_csv("brewery_ins.csv", index_col=0).dropna(axis=1)
    df.columns = df.columns.str.replace(" ", "")

    print(df)

    return df


def get_outlet(df, idx, compnent):

    outlet = df.loc[idx, compnent]

    return float(outlet)


def build_estimator(m, mass=False):
    df_brewery = get_brewwry_influent(m)
    m.fs.metab_effluent = Translator_ADM1_ASM1(
        inlet_property_package=m.fs.props_adm1,
        outlet_property_package=m.fs.props_asm1,
        reaction_package=m.fs.adm1_rxn_props,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )
    idx = 0
    comp_list = [
        "H2O",
        "S_su",
        "S_aa",
        "S_fa",
        "S_va",
        "S_bu",
        "S_pro",
        "S_ac",
        "S_h2",
        "S_ch4",
        "S_IC",
        "S_IN",
        "S_I",
        "X_c",
        "X_ch",
        "X_pr",
        "X_li",
        # "X_su",
        "X_aa",
        "X_fa",
        "X_c4",
        "X_pro",
        "X_ac",
        "X_h2",
        "X_I",
        "S_cat",
        "S_an",
    ]  #'S_co2' need to fix

    m.fs.metab_effluent.inlet.flow_vol.fix(
        float(get_outlet(df_brewery, idx, "VolumetricFlowrate"))
        * units.m**3
        / units.day
    )
    m.fs.metab_effluent.inlet.temperature.fix(308.15 * units.K)
    m.fs.metab_effluent.inlet.pressure.fix(1 * units.atm)

    for i in comp_list:
        if i == "S_cat":
            m.fs.metab_effluent.inlet.cations[0].fix(
                float(get_outlet(df_brewery, idx, i)) * units.mmol / units.liter
            )

        elif i == "S_an":
            m.fs.metab_effluent.inlet.anions[0].fix(
                float(get_outlet(df_brewery, idx, i)) * units.mmol / units.liter
            )
        elif i != "H2O":
            m.fs.metab_effluent.inlet.conc_mass_comp[0, i].fix(
                get_outlet(df_brewery, idx, i) * units.mg / units.liter
            )

    m.fs.metab_effluent.inlet.conc_mass_comp[0, "X_su"].fix(0 * units.mg / units.liter)

    return m


if __name__ == "__main__":
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.props_asm1 = ASM1ParameterBlock()
    m.fs.props_adm1 = ADM1ParameterBlock()
    m.fs.adm1_rxn_props = ADM1ReactionParameterBlock(property_package=m.fs.props_adm1)
    m = build_estimator(m, mass=True)
    model_checker(m, True)
    output = {
        # "brewery in ADM1": m.fs.metab_effluent.inlet,
        "brewery in ASM1": m.fs.metab_effluent.outlet,
    }
    stream_table = create_stream_table_dataframe(
        output,
        time_point=0,
    )
    print(stream_table_dataframe_to_string(stream_table))

    m.fs.metab_effluent.display()

    stream_table.to_csv("brewery_ins_asm1.csv")

    blk = m.fs.metab_effluent.properties_out[0]

    blk.TSS.display()
    blk.BOD5.display()
    blk.COD.display()
    blk.TKN.display()

    # print("TSS : ", print(blk.TSS[None]))
    # print("BOD5 : ", pyo.value(blk.BOD5[0]) )
    # print("COD : ", pyo.value(blk.COD[0]) )
    # print("TKN : ", pyo.value(blk.TKN[0]) )
