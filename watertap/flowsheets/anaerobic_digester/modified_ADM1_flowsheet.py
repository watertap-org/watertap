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
from watertap.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale
from idaes.core.util.tables import (
    create_stream_table_dataframe,
    stream_table_dataframe_to_string,
)
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
from watertap.unit_models.translators.translator_asm2d_adm1 import Translator_ASM2d_ADM1
from watertap.unit_models.anaerobic_digester import AD
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


def main(bio_P=False):
    m = build(bio_P=bio_P)
    set_operating_conditions(m)

    print("----------------   scaling V0  ----------------")
    badly_scaled_var_list = iscale.badly_scaled_var_generator(m, large=1e1, small=1e-1)
    for x in badly_scaled_var_list:
        print(f"{x[0].name}\t{x[0].value}\tsf: {iscale.get_scaling_factor(x[0])}")

    initialize_system(m)

    # print("----------------   Re-scaling V1  ----------------")
    # badly_scaled_var_list = iscale.badly_scaled_var_generator(m, large=1e1, small=1e-1)
    # for x in badly_scaled_var_list:
    #     if 1 < x[0].value < 10:
    #         sf = 1
    #     else:
    #         power = round(pyo.log10(abs(x[0].value)))
    #         sf = 1 / 10**power
    #     iscale.set_scaling_factor(x[0], sf)
    #
    # badly_scaled_var_list = iscale.badly_scaled_var_generator(m, large=1e1, small=1e-1)
    # for x in badly_scaled_var_list:
    #     print(f"{x[0].name}\t{x[0].value}\tsf: {iscale.get_scaling_factor(x[0])}")

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

    # print("---Numerical Issues---")
    # dt.report_numerical_issues()
    # dt.display_variables_at_or_outside_bounds()
    # dt.display_variables_with_extreme_jacobians()
    # dt.display_constraints_with_extreme_jacobians()

    return m, results


def build(bio_P=False):
    m = pyo.ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)

    # Properties
    m.fs.props_ASM2D = ModifiedASM2dParameterBlock()
    m.fs.rxn_props_ASM2D = ModifiedASM2dReactionParameterBlock(
        property_package=m.fs.props_ASM2D
    )
    m.fs.props_ADM1 = ModifiedADM1ParameterBlock()
    m.fs.props_vap_ADM1 = ADM1_vaporParameterBlock()
    m.fs.rxn_props_ADM1 = ModifiedADM1ReactionParameterBlock(
        property_package=m.fs.props_ADM1
    )

    # Feed water stream
    m.fs.FeedWater = Feed(property_package=m.fs.props_ASM2D)

    # ======================================================================
    # Anaerobic digester section
    # ASM2d-ADM1 translator
    m.fs.translator_asm2d_adm1 = Translator_ASM2d_ADM1(
        inlet_property_package=m.fs.props_ASM2D,
        outlet_property_package=m.fs.props_ADM1,
        inlet_reaction_package=m.fs.rxn_props_ASM2D,
        outlet_reaction_package=m.fs.rxn_props_ADM1,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
        bio_P=bio_P,
    )

    # Anaerobic digester
    m.fs.AD = AD(
        liquid_property_package=m.fs.props_ADM1,
        vapor_property_package=m.fs.props_vap_ADM1,
        reaction_package=m.fs.rxn_props_ADM1,
        has_heat_transfer=True,
        has_pressure_change=False,
    )

    # ADM1-ASM2d translator
    m.fs.translator_adm1_asm2d = Translator_ADM1_ASM2D(
        inlet_property_package=m.fs.props_ADM1,
        outlet_property_package=m.fs.props_ASM2D,
        inlet_reaction_package=m.fs.rxn_props_ADM1,
        outlet_reaction_package=m.fs.rxn_props_ASM2D,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )

    # Product Blocks
    m.fs.Treated = Product(property_package=m.fs.props_ASM2D)

    # ======================================================================
    # Link units related to AD section
    m.fs.stream_feed_translator = Arc(
        source=m.fs.FeedWater.outlet, destination=m.fs.translator_asm2d_adm1.inlet
    )
    m.fs.stream_translator_AD = Arc(
        source=m.fs.translator_asm2d_adm1.outlet, destination=m.fs.AD.inlet
    )
    m.fs.stream_AD_translator = Arc(
        source=m.fs.AD.liquid_outlet, destination=m.fs.translator_adm1_asm2d.inlet
    )
    m.fs.stream_translator_product = Arc(
        source=m.fs.translator_adm1_asm2d.outlet, destination=m.fs.Treated.inlet
    )
    pyo.TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def set_operating_conditions(m, bio_P=False):
    # Feed Water Conditions
    print(f"DOF before feed: {degrees_of_freedom(m)}")
    m.fs.FeedWater.flow_vol.fix(0.003 * pyo.units.m**3 / pyo.units.s)
    m.fs.FeedWater.temperature.fix(308.15 * pyo.units.K)
    m.fs.FeedWater.pressure.fix(1 * pyo.units.atm)

    if bio_P is True:
        m.fs.FeedWater.conc_mass_comp[0, "S_A"].fix(
            0.10034 * pyo.units.kg / pyo.units.m**3
        )
        m.fs.FeedWater.conc_mass_comp[0, "S_F"].fix(
            0.16118 * pyo.units.kg / pyo.units.m**3
        )
        m.fs.FeedWater.conc_mass_comp[0, "S_I"].fix(
            0.057450 * pyo.units.kg / pyo.units.m**3
        )
        m.fs.FeedWater.conc_mass_comp[0, "S_N2"].fix(
            0.024884 * pyo.units.kg / pyo.units.m**3
        )
        m.fs.FeedWater.conc_mass_comp[0, "S_NH4"].fix(
            0.040029 * pyo.units.kg / pyo.units.m**3
        )
        m.fs.FeedWater.conc_mass_comp[0, "S_NO3"].fix(
            1e-9 * pyo.units.kg / pyo.units.m**3
        )
        m.fs.FeedWater.conc_mass_comp[0, "S_O2"].fix(
            0.0014165 * pyo.units.kg / pyo.units.m**3
        )
        m.fs.FeedWater.conc_mass_comp[0, "S_PO4"].fix(
            0.026112 * pyo.units.kg / pyo.units.m**3
        )
        m.fs.FeedWater.conc_mass_comp[0, "S_K"].fix(
            0.37917 * pyo.units.kg / pyo.units.m**3
        )
        m.fs.FeedWater.conc_mass_comp[0, "S_Mg"].fix(
            0.027692 * pyo.units.kg / pyo.units.m**3
        )
        m.fs.FeedWater.conc_mass_comp[0, "S_IC"].fix(
            0.075367 * pyo.units.kg / pyo.units.m**3
        )
        m.fs.FeedWater.conc_mass_comp[0, "X_AUT"].fix(
            1e-9 * pyo.units.kg / pyo.units.m**3
        )
        m.fs.FeedWater.conc_mass_comp[0, "X_H"].fix(
            22.038 * pyo.units.kg / pyo.units.m**3
        )
        m.fs.FeedWater.conc_mass_comp[0, "X_I"].fix(
            10.794 * pyo.units.kg / pyo.units.m**3
        )
        m.fs.FeedWater.conc_mass_comp[0, "X_PAO"].fix(
            11.858 * pyo.units.kg / pyo.units.m**3
        )
        m.fs.FeedWater.conc_mass_comp[0, "X_PHA"].fix(
            0.0071834 * pyo.units.kg / pyo.units.m**3
        )
        m.fs.FeedWater.conc_mass_comp[0, "X_PP"].fix(
            3.1705 * pyo.units.kg / pyo.units.m**3
        )
        m.fs.FeedWater.conc_mass_comp[0, "X_S"].fix(
            3.7082 * pyo.units.kg / pyo.units.m**3
        )
    else:
        m.fs.FeedWater.conc_mass_comp[0, "S_A"].fix(
            0.095485 * pyo.units.kg / pyo.units.m**3
        )
        m.fs.FeedWater.conc_mass_comp[0, "S_F"].fix(
            0.14620 * pyo.units.kg / pyo.units.m**3
        )
        m.fs.FeedWater.conc_mass_comp[0, "S_I"].fix(
            0.057450 * pyo.units.kg / pyo.units.m**3
        )
        m.fs.FeedWater.conc_mass_comp[0, "S_N2"].fix(
            0.024893 * pyo.units.kg / pyo.units.m**3
        )
        m.fs.FeedWater.conc_mass_comp[0, "S_NH4"].fix(
            0.041267 * pyo.units.kg / pyo.units.m**3
        )
        m.fs.FeedWater.conc_mass_comp[0, "S_NO3"].fix(
            1e-9 * pyo.units.kg / pyo.units.m**3
        )
        m.fs.FeedWater.conc_mass_comp[0, "S_O2"].fix(
            0.0013782 * pyo.units.kg / pyo.units.m**3
        )
        m.fs.FeedWater.conc_mass_comp[0, "S_PO4"].fix(
            0.85735 * pyo.units.kg / pyo.units.m**3
        )
        m.fs.FeedWater.conc_mass_comp[0, "S_K"].fix(
            0.37612 * pyo.units.kg / pyo.units.m**3
        )
        m.fs.FeedWater.conc_mass_comp[0, "S_Mg"].fix(
            0.024419 * pyo.units.kg / pyo.units.m**3
        )
        m.fs.FeedWater.conc_mass_comp[0, "S_IC"].fix(
            0.074737 * pyo.units.kg / pyo.units.m**3
        )
        m.fs.FeedWater.conc_mass_comp[0, "X_AUT"].fix(
            1e-9 * pyo.units.kg / pyo.units.m**3
        )
        m.fs.FeedWater.conc_mass_comp[0, "X_H"].fix(
            22.310 * pyo.units.kg / pyo.units.m**3
        )
        m.fs.FeedWater.conc_mass_comp[0, "X_I"].fix(
            10.823 * pyo.units.kg / pyo.units.m**3
        )
        m.fs.FeedWater.conc_mass_comp[0, "X_PAO"].fix(
            11.268 * pyo.units.kg / pyo.units.m**3
        )
        m.fs.FeedWater.conc_mass_comp[0, "X_PHA"].fix(
            0.0056535 * pyo.units.kg / pyo.units.m**3
        )
        m.fs.FeedWater.conc_mass_comp[0, "X_PP"].fix(
            3.0925 * pyo.units.kg / pyo.units.m**3
        )
        m.fs.FeedWater.conc_mass_comp[0, "X_S"].fix(
            3.8081 * pyo.units.kg / pyo.units.m**3
        )

    # AD
    m.fs.AD.volume_liquid.fix(3400)
    m.fs.AD.volume_vapor.fix(300)
    m.fs.AD.liquid_outlet.temperature.fix(308.15)

    def scale_variables(m):
        for var in m.fs.component_data_objects(pyo.Var, descend_into=True):
            if "flow_vol" in var.name:
                iscale.set_scaling_factor(var, 1e1)
            if "temperature" in var.name:
                iscale.set_scaling_factor(var, 1e-2)
            if "pressure" in var.name:
                iscale.set_scaling_factor(var, 1e-5)
            if "conc_mass_comp" in var.name:
                if var.value > 1:
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
    # Initialize flowsheet
    m.fs.FeedWater.initialize(outlvl=idaeslog.INFO, solver="ipopt-watertap")
    propagate_state(m.fs.stream_feed_translator)
    m.fs.translator_asm2d_adm1.initialize(outlvl=idaeslog.INFO, solver="ipopt-watertap")
    propagate_state(m.fs.stream_translator_AD)
    m.fs.AD.initialize(outlvl=idaeslog.INFO, solver="ipopt-watertap")
    propagate_state(m.fs.stream_AD_translator)
    m.fs.translator_adm1_asm2d.initialize(outlvl=idaeslog.INFO, solver="ipopt-watertap")
    propagate_state(m.fs.stream_translator_product)
    m.fs.Treated.initialize(outlvl=idaeslog.INFO, solver="ipopt-watertap")


def solve(m, solver=None):
    if solver is None:
        solver = get_solver()
    results = solver.solve(m, tee=True)
    # pyo.assert_optimal_termination(results)
    return results


if __name__ == "__main__":
    # This method builds and runs a steady state activated sludge flowsheet.
    m, results = main(bio_P=True)

    stream_table = create_stream_table_dataframe(
        {
            "Feed": m.fs.FeedWater.outlet,
            # "ASM-ADM translator inlet": m.fs.translator_asm2d_adm1.inlet,
            # "AD inlet": m.fs.AD.inlet,
            "Treated water": m.fs.Treated.inlet,
        },
        time_point=0,
    )
    print(stream_table_dataframe_to_string(stream_table))
