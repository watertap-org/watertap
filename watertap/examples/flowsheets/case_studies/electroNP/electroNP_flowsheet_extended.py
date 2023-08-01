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
__author__ = "Chenyu Wang"

import pyomo.environ as pyo
from pyomo.environ import (
    value,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.network import Arc, SequentialDecomposition
from idaes.core import (
    FlowsheetBlock,
    UnitModelCostingBlock,
)
from idaes.core.solvers import get_solver
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models import (
    Product,
    Feed,
    CSTR,
    Mixer,
    Separator,
)
from idaes.models.unit_models.separator import SplittingType
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
    DecaySwitch,
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
from idaes.core.util.tables import (
    create_stream_table_dataframe,
    stream_table_dataframe_to_string,
)
from idaes.core.util.initialization import propagate_state
from watertap.core.util.initialization import check_solve
from watertap.costing import WaterTAPCosting

# Set up logger
_log = idaeslog.getLogger(__name__)


def main():
    # build, set, and initialize
    m = build_flowsheet()
    set_operating_conditions(m)
    initialize_system(m)

    # solve
    results = solve(m, tee=True)
    assert_optimal_termination(results)

    # display results
    stream_table = create_stream_table_dataframe(
        {
            "Feed inlet": m.fs.feed.outlet,
            "Translator outlet": m.fs.translator_asm2d_adm1.outlet,
            "AD liquid outlet": m.fs.AD.liquid_outlet,
            "AD vapor outlet": m.fs.AD.vapor_outlet,
            "Dewatering Unit outlet": m.fs.dewater.overflow,
            "ElectroNP treated": m.fs.electroNP.treated,
            "ElectroNP byproduct": m.fs.electroNP.byproduct,
        },
        time_point=0,
    )
    print(stream_table_dataframe_to_string(stream_table))
    display_costing(m)

    return m


def build_flowsheet():
    # flowsheet set up
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # properties package
    m.fs.props_ADM1 = ModifiedADM1ParameterBlock()
    m.fs.props_vap_ADM1 = ADM1_vaporParameterBlock()
    m.fs.rxn_props_ADM1 = ModifiedADM1ReactionParameterBlock(
        property_package=m.fs.props_ADM1
    )
    m.fs.props_ASM2D = ModifiedASM2dParameterBlock()
    m.fs.rxn_props_ASM2D = ModifiedASM2dReactionParameterBlock(
        property_package=m.fs.props_ASM2D
    )
    m.fs.costing = WaterTAPCosting()

    # Control volume flow blocks
    m.fs.feed = Feed(property_package=m.fs.props_ASM2D)

    # Unit models

    # Anaerobic digestor
    m.fs.AD = AD(
        liquid_property_package=m.fs.props_ADM1,
        vapor_property_package=m.fs.props_vap_ADM1,
        reaction_package=m.fs.rxn_props_ADM1,
        has_heat_transfer=True,
        has_pressure_change=False,
    )
    m.fs.AD.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    # Translators
    m.fs.translator_asm2d_adm1 = Translator_ASM2d_ADM1(
        inlet_property_package=m.fs.props_ASM2D,
        outlet_property_package=m.fs.props_ADM1,
        inlet_reaction_package=m.fs.rxn_props_ASM2D,
        outlet_reaction_package=m.fs.rxn_props_ADM1,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )

    m.fs.translator_adm1_asm2d = Translator_ADM1_ASM2D(
        inlet_property_package=m.fs.props_ADM1,
        outlet_property_package=m.fs.props_ASM2D,
        reaction_package=m.fs.rxn_props_ADM1,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )

    # ElectroNP
    m.fs.electroNP = ElectroNPZO(property_package=m.fs.props_ASM2D)
    m.fs.electroNP.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    # Dewatering Unit
    m.fs.dewater = DewateringUnit(
        property_package=m.fs.props_ASM2D,
        activated_sludge_model=dewater_type.modified_ASM2D,
    )

    # Thickener
    m.fs.thickener = Thickener(
        property_package=m.fs.props_ASM2D,
        activated_sludge_model=thickener_type.modified_ASM2D,
    )

    # Costing
    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(
        m.fs.electroNP.properties_treated[0].flow_vol
    )
    m.fs.costing.add_LCOW(m.fs.AD.inlet.flow_vol[0])

    # Connections
    m.fs.stream_feed_thickener = Arc(
        source=m.fs.feed.outlet, destination=m.fs.thickener.inlet
    )
    m.fs.stream_thickener_translator = Arc(
        source=m.fs.thickener.underflow, destination=m.fs.translator_asm2d_adm1.inlet
    )
    m.fs.stream_translator_AD = Arc(
        source=m.fs.translator_asm2d_adm1.outlet, destination=m.fs.AD.inlet
    )
    m.fs.stream_AD_translator = Arc(
        source=m.fs.AD.liquid_outlet, destination=m.fs.translator_adm1_asm2d.inlet
    )
    m.fs.stream_translator_dewater = Arc(
        source=m.fs.translator_adm1_asm2d.outlet, destination=m.fs.dewater.inlet
    )
    m.fs.stream_dewater_electroNP = Arc(
        source=m.fs.dewater.overflow, destination=m.fs.electroNP.inlet
    )
    pyo.TransformationFactory("network.expand_arcs").apply_to(m)

    # Scaling
    for var in m.fs.component_data_objects(pyo.Var, descend_into=True):
        if "flow_vol" in var.name:
            iscale.set_scaling_factor(var, 1e3)
        if "temperature" in var.name:
            iscale.set_scaling_factor(var, 1e-1)
        if "pressure" in var.name:
            iscale.set_scaling_factor(var, 1e-3)
        if "alkalinity" in var.name:
            iscale.set_scaling_factor(var, 1e-1)
        if "conc_mass_comp" in var.name:
            iscale.set_scaling_factor(var, 1e1)
        # if "conc_mass_comp[S_IN]" in var.name:
        #     iscale.set_scaling_factor(var, 1e-1)
        # if "conc_mass_comp[S_IP]" in var.name:
        #     iscale.set_scaling_factor(var, 1e-1)
        # if "conc_mass_comp[S_PO4]" in var.name:
        #     iscale.set_scaling_factor(var, 1e-1)
        # if "conc_mass_comp[S_NH4]" in var.name:
        #     iscale.set_scaling_factor(var, 1e1)
        # if "conc_mass_comp[S_F]" in var.name:
        #     iscale.set_scaling_factor(var, 1e-1)
        # if "conc_mass_comp[X_I]" in var.name:
        #     iscale.set_scaling_factor(var, 1e-1)

    iscale.calculate_scaling_factors(m)

    return m


def set_operating_conditions(m):
    # Feed inlet
    m.fs.feed.properties[0].pressure.fix(1 * pyunits.atm)
    m.fs.feed.properties[0].temperature.fix(308.15 * pyunits.K)
    m.fs.feed.properties[0].flow_vol.fix(178.4674 * pyunits.m**3 / pyunits.day)
    eps = 1e-9 * pyunits.kg / pyunits.m**3

    m.fs.feed.properties[0].conc_mass_comp["S_O2"].fix(eps)
    m.fs.feed.properties[0].conc_mass_comp["S_F"].fix(
        0.02644 * pyunits.kg / pyunits.m**3
    )
    m.fs.feed.properties[0].conc_mass_comp["S_A"].fix(
        0.01766 * pyunits.kg / pyunits.m**3
    )
    m.fs.feed.properties[0].conc_mass_comp["S_I"].fix(
        0.02723 * pyunits.kg / pyunits.m**3
    )
    m.fs.feed.properties[0].conc_mass_comp["S_NH4"].fix(
        0.01858 * pyunits.kg / pyunits.m**3
    )
    m.fs.feed.properties[0].conc_mass_comp["S_N2"].fix(
        0.00507 * pyunits.kg / pyunits.m**3
    )
    m.fs.feed.properties[0].conc_mass_comp["S_NO3"].fix(
        0.00002 * pyunits.kg / pyunits.m**3
    )
    m.fs.feed.properties[0].conc_mass_comp["S_PO4"].fix(
        0.00469 * pyunits.kg / pyunits.m**3
    )
    m.fs.feed.properties[0].conc_mass_comp["S_IC"].fix(
        0.07899 * pyunits.kg / pyunits.m**3
    )

    m.fs.feed.properties[0].conc_mass_comp["X_I"].fix(
        10.96441 * pyunits.kg / pyunits.m**3
    )
    m.fs.feed.properties[0].conc_mass_comp["X_S"].fix(
        19.08476 * pyunits.kg / pyunits.m**3
    )
    m.fs.feed.properties[0].conc_mass_comp["X_H"].fix(
        9.47939 * pyunits.kg / pyunits.m**3
    )
    m.fs.feed.properties[0].conc_mass_comp["X_PAO"].fix(
        3.8622 * pyunits.kg / pyunits.m**3
    )
    m.fs.feed.properties[0].conc_mass_comp["X_PP"].fix(
        0.45087 * pyunits.kg / pyunits.m**3
    )
    m.fs.feed.properties[0].conc_mass_comp["X_PHA"].fix(
        0.02464 * pyunits.kg / pyunits.m**3
    )
    m.fs.feed.properties[0].conc_mass_comp["X_AUT"].fix(
        0.33379 * pyunits.kg / pyunits.m**3
    )
    m.fs.feed.properties[0].conc_mass_comp["S_K"].fix(
        0.01979 * pyunits.kg / pyunits.m**3
    )
    m.fs.feed.properties[0].conc_mass_comp["S_Mg"].fix(
        0.18987 * pyunits.kg / pyunits.m**3
    )

    print("DOF before AD:", degrees_of_freedom(m))

    # AD
    m.fs.AD.volume_liquid.fix(3400)
    m.fs.AD.volume_vapor.fix(300)
    m.fs.AD.liquid_outlet.temperature.fix(308.15)

    # ElectroNP
    m.fs.electroNP.energy_electric_flow_mass.fix(0.044 * pyunits.kWh / pyunits.kg)
    m.fs.electroNP.magnesium_chloride_dosage.fix(0.388)

    # Costing
    m.fs.costing.electroNP.phosphorus_recovery_value = 0


def initialize_system(m):
    # Initialize
    seq = SequentialDecomposition()
    seq.options.select_tear_method = "heuristic"
    seq.options.tear_method = "Wegstein"
    seq.options.iterLim = 0

    def function(unit):
        unit.initialize(outlvl=idaeslog.INFO)

    seq.run(m, function)

    # propagate_state(m.fs.stream_feed_thickener)
    # m.fs.thickener.initialize()
    # propagate_state(m.fs.stream_thickener_translator)
    # m.fs.translator_asm2d_adm1.initialize(outlvl=idaeslog.INFO_HIGH)
    # propagate_state(m.fs.stream_translator_AD)
    # m.fs.AD.initialize(outlvl=idaeslog.INFO_HIGH)
    # propagate_state((m.fs.stream_AD_translator))
    # m.fs.translator_adm1_asm2d.initialize(outlvl=idaeslog.INFO_HIGH)
    # propagate_state(m.fs.stream_translator_dewater)
    # m.fs.dewater.initialize(outlvl=idaeslog.INFO_HIGH)
    # propagate_state(m.fs.stream_dewater_electroNP)
    # m.fs.electroNP.initialize(outlvl=idaeslog.INFO_HIGH)
    m.fs.costing.initialize()


def solve(blk, solver=None, checkpoint=None, tee=False, fail_flag=True):
    if solver is None:
        solver = get_solver()
    results = solver.solve(blk, tee=tee)
    check_solve(results, checkpoint=checkpoint, logger=_log, fail_flag=fail_flag)
    return results


def display_costing(m):
    print("\n----------Capital Cost----------")
    total_capital_cost = value(
        pyunits.convert(m.fs.costing.total_capital_cost, to_units=pyunits.USD_2018)
    )
    normalized_capex = total_capital_cost / value(
        pyunits.convert(m.fs.AD.inlet.flow_vol[0], to_units=pyunits.m**3 / pyunits.hr)
    )
    print(f"Total Capital Costs: {total_capital_cost:.3f} $")
    print(f"Normalized Capital Costs: {normalized_capex:.3f} $/m3/hr")
    print("Capital Cost Breakdown")
    for u in m.fs.costing._registered_unit_costing:
        print(
            u.name,
            " : {price:0.3f} $".format(
                price=value(pyunits.convert(u.capital_cost, to_units=pyunits.USD_2018))
            ),
        )
    print("\n----------Operation Cost----------")
    total_operating_cost = value(
        pyunits.convert(
            m.fs.costing.total_operating_cost, to_units=pyunits.USD_2018 / pyunits.year
        )
    )
    print(f"Total Operating Cost: {total_operating_cost:.3f} $/year")

    opex_fraction = value(
        pyunits.convert(
            m.fs.costing.total_operating_cost, to_units=pyunits.USD_2018 / pyunits.year
        )
        / pyunits.convert(
            m.fs.AD.inlet.flow_vol[0], to_units=pyunits.m**3 / pyunits.year
        )
        / m.fs.costing.LCOW
    )
    print(f"Operating cost fraction: {opex_fraction:.3f} $ opex / $ LCOW")

    print("Operating Cost Breakdown")
    for f in m.fs.costing.used_flows:
        print(
            f.title(),
            " :    {price:0.3f} $/m3 feed".format(
                price=value(
                    pyunits.convert(
                        m.fs.costing.aggregate_flow_costs[f],
                        to_units=pyunits.USD_2018 / pyunits.year,
                    )
                    / pyunits.convert(
                        m.fs.AD.inlet.flow_vol[0],
                        to_units=pyunits.m**3 / pyunits.year,
                    )
                )
            ),
        )

    print("\n----------Energy----------")

    electricity_intensity = value(
        pyunits.convert(
            m.fs.costing.aggregate_flow_electricity / m.fs.AD.inlet.flow_vol[0],
            to_units=pyunits.kWh / pyunits.m**3,
        )
    )
    print(f"Electricity Intensity: {electricity_intensity:.4f} kWh/m3")

    print("\n----------Levelized Cost----------")
    LCOW = value(
        pyunits.convert(m.fs.costing.LCOW, to_units=pyunits.USD_2018 / pyunits.m**3)
    )
    print(f"Levelized Cost of Water: {LCOW:.3f} $/m^3")


if __name__ == "__main__":
    m = main()
