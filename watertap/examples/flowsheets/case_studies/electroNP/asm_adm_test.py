#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
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
from pyomo.environ import (
    value,
    units as pyunits,
    SolverFactory,
    TransformationFactory,
    Suffix,
)
from pyomo.network import Arc, SequentialDecomposition

from idaes.core import (
    FlowsheetBlock,
    UnitModelCostingBlock,
)
from idaes.models.unit_models import (
    CSTR,
    Feed,
    Separator,
    Product,
    Mixer,
    PressureChanger,
)
from idaes.models.unit_models.separator import SplittingType
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import propagate_state
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale
from idaes.core.util.tables import (
    create_stream_table_dataframe,
    stream_table_dataframe_to_string,
)
from watertap.unit_models.cstr_injection import CSTR_Injection
from watertap.property_models.anaerobic_digestion.modified_adm1_properties import (
    ModifiedADM1ParameterBlock,
)
from watertap.property_models.anaerobic_digestion.adm1_properties_vapor import (
    ADM1_vaporParameterBlock,
)
from watertap.property_models.anaerobic_digestion.modified_adm1_reactions import (
    ModifiedADM1ReactionParameterBlock,
)
from watertap.property_models.activated_sludge.modified_asm2d_properties import (
    ModifiedASM2dParameterBlock,
)
from watertap.property_models.activated_sludge.modified_asm2d_reactions import (
    ModifiedASM2dReactionParameterBlock,
)
from watertap.unit_models.translators.translator_adm1_asm2d import (
    Translator_ADM1_ASM2D,
)
from watertap.unit_models.translators.translator_asm2d_adm1 import Translator_ASM2d_ADM1
from watertap.unit_models.anaerobic_digestor import AD
from watertap.unit_models.electroNP_ZO import ElectroNPZO
from watertap.unit_models.dewatering import (
    DewateringUnit,
    ActivatedSludgeModelType as dewater_type,
)
from watertap.unit_models.thickener import (
    Thickener,
    ActivatedSludgeModelType as thickener_type,
)
from watertap.core.util.initialization import check_solve
from watertap.costing import WaterTAPCosting

from watertap.core.util.model_diagnostics.infeasible import *
from idaes.core.util.model_diagnostics import DegeneracyHunter

from idaes.core.initialization import BlockTriangularizationInitializer


# Set up logger
_log = idaeslog.getLogger(__name__)


def main():
    m = build()

    set_operating_conditions(m)
    set_scaling(m)

    badly_scaled_var_list = iscale.badly_scaled_var_generator(m, large=1e2, small=1e-2)
    print("----------------   badly_scaled_var_list   ----------------")
    for x in badly_scaled_var_list:
        print(f"{x[0].name}\t{x[0].value}\tsf: {iscale.get_scaling_factor(x[0])}")

    print("Structural issues after setting operating conditions")
    dt = DiagnosticsToolbox(model=m)
    dt.report_structural_issues()

    initialize_system(m)
    print("Numerical issues after initialization")
    dt.report_numerical_issues()

    results = solve(m)
    print("Numerical issues after solving")
    dt.report_numerical_issues()

    display_results(m)

    return m, results


def build():
    m = pyo.ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.props_ASM2D = ModifiedASM2dParameterBlock()
    m.fs.rxn_props_ASM2D = ModifiedASM2dReactionParameterBlock(
        property_package=m.fs.props_ASM2D
    )
    m.fs.props_ADM1 = ModifiedADM1ParameterBlock()
    m.fs.props_vap_ADM1 = ADM1_vaporParameterBlock()
    m.fs.rxn_props_ADM1 = ModifiedADM1ReactionParameterBlock(
        property_package=m.fs.props_ADM1
    )

    m.fs.costing = WaterTAPCosting()

    # Feed water stream
    m.fs.FeedWater = Feed(property_package=m.fs.props_ASM2D)

    # Translators
    m.fs.translator_asm2d_adm1 = Translator_ASM2d_ADM1(
        inlet_property_package=m.fs.props_ASM2D,
        outlet_property_package=m.fs.props_ADM1,
        inlet_reaction_package=m.fs.rxn_props_ASM2D,
        outlet_reaction_package=m.fs.rxn_props_ADM1,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )

    # Link units
    m.fs.stream1 = Arc(
        source=m.fs.FeedWater.outlet, destination=m.fs.translator_asm2d_adm1.inlet
    )

    pyo.TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def set_scaling(m):
    m.scaling_factor = Suffix(direction=Suffix.EXPORT)

    iscale.set_scaling_factor(m.fs.translator_asm2d_adm1.properties_in[0.0].flow_vol, 1)
    iscale.set_scaling_factor(
        m.fs.translator_asm2d_adm1.properties_out[0.0].flow_vol, 1
    )
    iscale.set_scaling_factor(
        m.fs.translator_asm2d_adm1.properties_in[0.0].temperature, 1e-2
    )
    iscale.set_scaling_factor(
        m.fs.translator_asm2d_adm1.properties_in[0.0].pressure, 1e-5
    )

    # iscale.set_scaling_factor(
    #     m.fs.translator_asm2d_adm1.properties_out[0.0].conc_mass_comp["S_su"], 1
    # )
    # iscale.set_scaling_factor(
    #     m.fs.translator_asm2d_adm1.properties_out[0.0].conc_mass_comp["S_fa"], 1e-1
    # )
    # iscale.set_scaling_factor(
    #     m.fs.translator_asm2d_adm1.properties_out[0.0].conc_mass_comp["S_ch4"], 1e9
    # )
    # iscale.set_scaling_factor(
    #     m.fs.translator_asm2d_adm1.properties_out[0.0].conc_mass_comp["S_IC"], 1
    # )
    # iscale.set_scaling_factor(
    #     m.fs.translator_asm2d_adm1.properties_out[0.0].conc_mass_comp["S_IN"], 1
    # )
    # iscale.set_scaling_factor(
    #     m.fs.translator_asm2d_adm1.properties_out[0.0].conc_mass_comp["S_IP"], 1
    # )
    # iscale.set_scaling_factor(
    #     m.fs.translator_asm2d_adm1.properties_out[0.0].conc_mass_comp["X_su"], 1e9
    # )
    # iscale.set_scaling_factor(
    #     m.fs.translator_asm2d_adm1.properties_out[0.0].conc_mass_comp["X_fa"], 1e9
    # )
    # iscale.set_scaling_factor(
    #     m.fs.translator_asm2d_adm1.properties_out[0.0].conc_mass_comp["X_c4"], 1e9
    # )
    # iscale.set_scaling_factor(
    #     m.fs.translator_asm2d_adm1.properties_out[0.0].conc_mass_comp["X_pro"], 1e9
    # )
    # iscale.set_scaling_factor(
    #     m.fs.translator_asm2d_adm1.properties_out[0.0].conc_mass_comp["X_ac"], 1e9
    # )
    # iscale.set_scaling_factor(
    #     m.fs.translator_asm2d_adm1.properties_out[0.0].conc_mass_comp["X_I"], 1e-1
    # )
    # iscale.set_scaling_factor(
    #     m.fs.translator_asm2d_adm1.properties_out[0.0].conc_mass_comp["X_PHA"], 1
    # )


def set_operating_conditions(m):
    m.fs.FeedWater.flow_vol.fix(0.00018359)
    m.fs.FeedWater.temperature.fix(308.15 * pyo.units.K)
    m.fs.FeedWater.pressure.fix(1 * pyo.units.atm)
    m.fs.FeedWater.conc_mass_comp[0, "S_O2"].fix(0.0079024)
    m.fs.FeedWater.conc_mass_comp[0, "S_F"].fix(0.00030596)
    m.fs.FeedWater.conc_mass_comp[0, "S_A"].fix(4.1385e-05)
    m.fs.FeedWater.conc_mass_comp[0, "S_NH4"].fix(0.0070896)
    m.fs.FeedWater.conc_mass_comp[0, "S_NO3"].fix(1.1750e-07)
    m.fs.FeedWater.conc_mass_comp[0, "S_PO4"].fix(0.0027761)
    m.fs.FeedWater.conc_mass_comp[0, "S_I"].fix(0.03)
    m.fs.FeedWater.conc_mass_comp[0, "S_N2"].fix(5.2387e-07)
    m.fs.FeedWater.conc_mass_comp[0, "X_I"].fix(19.794)
    m.fs.FeedWater.conc_mass_comp[0, "X_S"].fix(1.2980)
    m.fs.FeedWater.conc_mass_comp[0, "X_H"].fix(64.284)
    m.fs.FeedWater.conc_mass_comp[0, "X_PAO"].fix(3.7115e-07)
    m.fs.FeedWater.conc_mass_comp[0, "X_PP"].fix(7.3446e-07)
    m.fs.FeedWater.conc_mass_comp[0, "X_PHA"].fix(1.2442e-06)
    m.fs.FeedWater.conc_mass_comp[0, "X_AUT"].fix(7.6090e-05)
    m.fs.FeedWater.properties[0].conc_mass_comp["S_IC"].fix(0.10470)
    m.fs.FeedWater.properties[0].conc_mass_comp["S_K"].fix(2.4407e-08)
    m.fs.FeedWater.properties[0].conc_mass_comp["S_Mg"].fix(2.4425e-08)


def initialize_system(m):
    seq = SequentialDecomposition()
    seq.options.select_tear_method = "heuristic"
    seq.options.tear_method = "Wegstein"
    seq.options.iterLim = 1

    def function(unit):
        unit.initialize(outlvl=idaeslog.DEBUG, optarg={"bound_push": 1e-2})

    seq.run(m, function)


def solve(m):
    solver = SolverFactory("ipopt")
    results = solver.solve(m, tee=True)

    return results


def display_results(m):
    # m.fs.translator_asm2d_adm1.display()
    stream_table = create_stream_table_dataframe(
        {
            "Feed": m.fs.FeedWater.outlet,
            "ASM-ADM translator inlet": m.fs.translator_asm2d_adm1.inlet,
            "ASM-ADM translator outlet": m.fs.translator_asm2d_adm1.outlet,
        },
        time_point=0,
    )
    print(stream_table_dataframe_to_string(stream_table))


if __name__ == "__main__":
    m, results = main()
