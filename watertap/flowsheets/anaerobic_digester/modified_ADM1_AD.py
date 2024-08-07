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
Flowsheet example full Water Resource Recovery Facility
(WRRF; a.k.a., wastewater treatment plant) with ASM2d and ADM1 with P extension.

The flowsheet follows the same formulation as benchmark simulation model no.2 (BSM2)
but comprises different specifications for default values than BSM2.
"""

# Some more information about this module
__author__ = "Chenyu Wang, Adam Atia, Alejandro Garciadiego, Marcus Holly"

import pyomo.environ as pyo
from pyomo.network import Arc, SequentialDecomposition

from idaes.core import (
    FlowsheetBlock,
    # UnitModelCostingBlock,
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
from watertap.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale
from idaes.core.util.tables import (
    create_stream_table_dataframe,
    stream_table_dataframe_to_string,
)
from watertap.unit_models.cstr_injection import CSTR_Injection
from watertap.unit_models.clarifier import Clarifier
from watertap.property_models.unit_specific.anaerobic_digestion.modified_adm1_properties import (
    ModifiedADM1ParameterBlock,
)
from watertap.property_models.unit_specific.anaerobic_digestion.adm1_properties_vapor import (
    ADM1_vaporParameterBlock,
)
from watertap.property_models.unit_specific.anaerobic_digestion.modified_adm1_reactions import (
    ModifiedADM1ReactionParameterBlock,
)
from watertap.property_models.unit_specific.activated_sludge.modified_asm2d_properties import (
    ModifiedASM2dParameterBlock,
)
from watertap.property_models.unit_specific.activated_sludge.modified_asm2d_reactions import (
    ModifiedASM2dReactionParameterBlock,
)
from watertap.unit_models.translators.translator_adm1_asm2d import (
    Translator_ADM1_ASM2D,
)
from idaes.models.unit_models.mixer import MomentumMixingType
from watertap.unit_models.translators.translator_asm2d_adm1 import Translator_ASM2d_ADM1
from watertap.unit_models.anaerobic_digester import AD
from watertap.unit_models.dewatering import (
    DewateringUnit,
    ActivatedSludgeModelType as dewater_type,
)
from watertap.unit_models.thickener import (
    Thickener,
    ActivatedSludgeModelType as thickener_type,
)

from watertap.core.util.initialization import (
    check_solve,
    # assert_degrees_of_freedom
)

# from watertap.costing import WaterTAPCosting
# from watertap.costing.unit_models.clarifier import (
#     cost_circular_clarifier,
#     cost_primary_clarifier,
# )

from idaes.core.util.initialization import (
    propagate_state as _pro_state,
)
from idaes.core.util import DiagnosticsToolbox

# Set up logger
_log = idaeslog.getLogger(__name__)


def propagate_state(arc):
    _pro_state(arc)
    print(arc.destination.name)
    arc.destination.display()


def main():
    m = build()
    set_operating_conditions(m)

    initialize_system(m)

    print("----------------   Re-scaling V1  ----------------")
    badly_scaled_var_list = iscale.badly_scaled_var_generator(m, large=1e1, small=1e-1)
    for x in badly_scaled_var_list:
        if 1 < x[0].value < 10:
            sf = 1
        else:
            power = round(pyo.log10(abs(x[0].value)))
            sf = 1 / 10**power
        iscale.set_scaling_factor(x[0], sf)

    badly_scaled_var_list = iscale.badly_scaled_var_generator(m, large=1e1, small=1e-1)
    for x in badly_scaled_var_list:
        print(f"{x[0].name}\t{x[0].value}\tsf: {iscale.get_scaling_factor(x[0])}")

    # dt = DiagnosticsToolbox(m)
    # print("---Structural Issues---")
    # dt.report_structural_issues()
    # dt.display_potential_evaluation_errors()

    # print("----------------   Degen Hunter  ----------------")
    # # Use of Degeneracy Hunter for troubleshooting model.
    # m.obj = pyo.Objective(expr=0)
    # solver = get_solver()
    # solver.options["max_iter"] = 10000
    # results = solver.solve(m, tee=True)

    results = solve(m)

    pyo.assert_optimal_termination(results)
    check_solve(
        results,
        checkpoint="re-solve with controls in place",
        logger=_log,
        fail_flag=True,
    )

    print("---Numerical Issues---")
    dt.report_numerical_issues()
    # dt.display_variables_at_or_outside_bounds()
    # dt.display_variables_with_extreme_jacobians()
    # dt.display_constraints_with_extreme_jacobians()

    # add_costing(m)
    # m.fs.costing.initialize()
    #
    # assert_degrees_of_freedom(m, 0)
    #
    # results = solve(m)
    # pyo.assert_optimal_termination(results)
    #
    # display_costing(m)
    # display_performance_metrics(m)

    return m, results


def build():
    m = pyo.ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)

    # Properties
    m.fs.props_ADM1 = ModifiedADM1ParameterBlock()
    m.fs.props_vap_ADM1 = ADM1_vaporParameterBlock()
    m.fs.rxn_props_ADM1 = ModifiedADM1ReactionParameterBlock(
        property_package=m.fs.props_ADM1
    )

    # Anaerobic digester
    m.fs.AD = AD(
        liquid_property_package=m.fs.props_ADM1,
        vapor_property_package=m.fs.props_vap_ADM1,
        reaction_package=m.fs.rxn_props_ADM1,
        has_heat_transfer=True,
        has_pressure_change=False,
    )

    return m


def set_operating_conditions(m):
    # Feed Water Conditions
    print(f"DOF before feed: {degrees_of_freedom(m)}")
    m.fs.AD.inlet.flow_vol.fix(0.003 * pyo.units.m**3 / pyo.units.s)
    m.fs.AD.inlet.temperature.fix(308.15 * pyo.units.K)
    m.fs.AD.inlet.pressure.fix(1 * pyo.units.atm)

    m.fs.AD.inlet.conc_mass_comp[0, "S_I"].fix(0.057450 * pyo.units.kg / pyo.units.m**3)
    m.fs.AD.inlet.conc_mass_comp[0, "S_K"].fix(1.4070 * pyo.units.kg / pyo.units.m**3)
    m.fs.AD.inlet.conc_mass_comp[0, "S_Mg"].fix(1.0553 * pyo.units.kg / pyo.units.m**3)
    m.fs.AD.inlet.conc_mass_comp[0, "S_IC"].fix(1.1566 * pyo.units.kg / pyo.units.m**3)
    m.fs.AD.inlet.conc_mass_comp[0, "X_I"].fix(14.181 * pyo.units.kg / pyo.units.m**3)
    m.fs.AD.inlet.conc_mass_comp[0, "X_PAO"].fix(1e-9 * pyo.units.kg / pyo.units.m**3)
    m.fs.AD.inlet.conc_mass_comp[0, "X_PHA"].fix(1e-9 * pyo.units.kg / pyo.units.m**3)
    m.fs.AD.inlet.conc_mass_comp[0, "X_PP"].fix(1e-9 * pyo.units.kg / pyo.units.m**3)

    m.fs.AD.inlet.conc_mass_comp[0, "S_su"].fix(0.10191 * pyo.units.kg / pyo.units.m**3)
    m.fs.AD.inlet.conc_mass_comp[0, "S_aa"].fix(
        0.044290 * pyo.units.kg / pyo.units.m**3
    )
    m.fs.AD.inlet.conc_mass_comp[0, "S_fa"].fix(1e-9 * pyo.units.kg / pyo.units.m**3)
    m.fs.AD.inlet.conc_mass_comp[0, "S_va"].fix(1e-9 * pyo.units.kg / pyo.units.m**3)
    m.fs.AD.inlet.conc_mass_comp[0, "S_bu"].fix(1e-9 * pyo.units.kg / pyo.units.m**3)
    m.fs.AD.inlet.conc_mass_comp[0, "S_pro"].fix(1e-9 * pyo.units.kg / pyo.units.m**3)
    m.fs.AD.inlet.conc_mass_comp[0, "S_ac"].fix(
        0.091810 * pyo.units.kg / pyo.units.m**3
    )
    m.fs.AD.inlet.conc_mass_comp[0, "S_h2"].fix(1e-9 * pyo.units.kg / pyo.units.m**3)
    m.fs.AD.inlet.conc_mass_comp[0, "S_ch4"].fix(1e-9 * pyo.units.kg / pyo.units.m**3)
    m.fs.AD.inlet.conc_mass_comp[0, "S_IN"].fix(1.7106 * pyo.units.kg / pyo.units.m**3)
    m.fs.AD.inlet.conc_mass_comp[0, "S_IP"].fix(4.5311 * pyo.units.kg / pyo.units.m**3)
    m.fs.AD.inlet.conc_mass_comp[0, "X_ch"].fix(10.296 * pyo.units.kg / pyo.units.m**3)
    m.fs.AD.inlet.conc_mass_comp[0, "X_pr"].fix(10.388 * pyo.units.kg / pyo.units.m**3)
    m.fs.AD.inlet.conc_mass_comp[0, "X_li"].fix(13.346 * pyo.units.kg / pyo.units.m**3)
    m.fs.AD.inlet.conc_mass_comp[0, "X_su"].fix(1e-9 * pyo.units.kg / pyo.units.m**3)
    m.fs.AD.inlet.conc_mass_comp[0, "X_aa"].fix(1e-9 * pyo.units.kg / pyo.units.m**3)
    m.fs.AD.inlet.conc_mass_comp[0, "X_fa"].fix(1e-9 * pyo.units.kg / pyo.units.m**3)
    m.fs.AD.inlet.conc_mass_comp[0, "X_c4"].fix(1e-9 * pyo.units.kg / pyo.units.m**3)
    m.fs.AD.inlet.conc_mass_comp[0, "X_pro"].fix(1e-9 * pyo.units.kg / pyo.units.m**3)
    m.fs.AD.inlet.conc_mass_comp[0, "X_ac"].fix(1e-9 * pyo.units.kg / pyo.units.m**3)
    m.fs.AD.inlet.conc_mass_comp[0, "X_h2"].fix(1e-9 * pyo.units.kg / pyo.units.m**3)

    # AD
    m.fs.AD.volume_liquid.fix(3400)
    m.fs.AD.volume_vapor.fix(300)
    m.fs.AD.liquid_outlet.temperature.fix(308.15)

    def scale_variables(m):
        for var in m.fs.component_data_objects(pyo.Var, descend_into=True):
            if "flow_vol" in var.name:
                iscale.set_scaling_factor(var, 1e3)
            if "temperature" in var.name:
                iscale.set_scaling_factor(var, 1e-2)
            if "pressure" in var.name:
                iscale.set_scaling_factor(var, 1e-4)
            if "conc_mass_comp" in var.name:
                iscale.set_scaling_factor(var, 1e2)
                if var.value >= 1:
                    sf = 1e0
                    iscale.set_scaling_factor(var, sf)
                elif 1e-2 < var.value < 1:
                    sf = 1e1
                    iscale.set_scaling_factor(var, sf)
                else:
                    sf = 1e2
                    iscale.set_scaling_factor(var, sf)

    # Apply scaling
    scale_variables(m)
    iscale.calculate_scaling_factors(m.fs)


def initialize_system(m):
    # # Initialize flowsheet
    # # Apply sequential decomposition - 1 iteration should suffice
    # seq = SequentialDecomposition()
    # seq.options.select_tear_method = "heuristic"
    # seq.options.iterLim = 0
    #
    # G = seq.create_graph(m)
    # # Uncomment this code to see tear set and initialization order
    # order = seq.calculation_order(G)
    # print("Initialization Order")
    # for o in order:
    #     print(o[0].name)
    #
    # def function(unit):
    #     unit.initialize(outlvl=idaeslog.INFO, solver="ipopt-watertap")
    #
    # seq.run(m, function)

    m.fs.AD.initialize(outlvl=idaeslog.INFO, solver="ipopt-watertap")


def solve(m, solver=None):
    if solver is None:
        solver = get_solver()
    results = solver.solve(m, tee=True)
    # pyo.assert_optimal_termination(results)
    return results


if __name__ == "__main__":
    # This method builds and runs a steady state activated sludge flowsheet.
    m, results = main()

    # stream_table = create_stream_table_dataframe(
    #     {
    #         "Feed": m.fs.FeedWater.outlet,
    #         "ASM-ADM translator inlet": m.fs.translator_asm2d_adm1.inlet,
    #         # "R1": m.fs.R1.outlet,
    #         # "R2": m.fs.R2.outlet,
    #         # "R3": m.fs.R3.outlet,
    #         # "R4": m.fs.R4.outlet,
    #         # "R5": m.fs.R5.outlet,
    #         # "R6": m.fs.R6.outlet,
    #         # "R7": m.fs.R7.outlet,
    #         # "thickener outlet": m.fs.thickener.underflow,
    #         # "ADM-ASM translator outlet": m.fs.translator_adm1_asm2d.outlet,
    #         # "dewater outlet": m.fs.dewater.overflow,
    #         # "Treated water": m.fs.Treated.inlet,
    #         # "Sludge": m.fs.Sludge.inlet,
    #     },
    #     time_point=0,
    # )
    # print(stream_table_dataframe_to_string(stream_table))
