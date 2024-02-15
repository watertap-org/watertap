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
    units,
    value,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.network import Arc
from idaes.core import (
    FlowsheetBlock,
    UnitModelCostingBlock,
)
from idaes.core.solvers import get_solver
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale
from watertap.unit_models.anaerobic_digester import AD
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
from watertap.unit_models.translators.translator_adm1_asm2d import (
    Translator_ADM1_ASM2D,
)
from watertap.unit_models.electroNP_ZO import ElectroNPZO
from idaes.core.util.tables import (
    create_stream_table_dataframe,
    stream_table_dataframe_to_string,
)
from idaes.core.util.initialization import propagate_state
from watertap.core.util.initialization import check_solve
from watertap.costing import WaterTAPCosting

# Set up logger
_log = idaeslog.getLogger(__name__)


def build_flowsheet():
    # flowsheet set up
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.props_ADM1 = ModifiedADM1ParameterBlock()
    m.fs.props_vap_ADM1 = ADM1_vaporParameterBlock()
    m.fs.rxn_props_ADM1 = ModifiedADM1ReactionParameterBlock(
        property_package=m.fs.props_ADM1
    )
    m.fs.props_ASM2D = ModifiedASM2dParameterBlock()
    m.fs.costing = WaterTAPCosting()

    # Unit models
    m.fs.AD = AD(
        liquid_property_package=m.fs.props_ADM1,
        vapor_property_package=m.fs.props_vap_ADM1,
        reaction_package=m.fs.rxn_props_ADM1,
        has_heat_transfer=True,
        has_pressure_change=False,
    )
    m.fs.AD.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    m.fs.translator_adm1_asm2d = Translator_ADM1_ASM2D(
        inlet_property_package=m.fs.props_ADM1,
        outlet_property_package=m.fs.props_ASM2D,
        reaction_package=m.fs.rxn_props_ADM1,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )

    m.fs.electroNP = ElectroNPZO(property_package=m.fs.props_ASM2D)
    m.fs.electroNP.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(m.fs.electroNP.treated.flow_vol[0])
    m.fs.costing.add_LCOW(m.fs.AD.inlet.flow_vol[0])

    # connections
    m.fs.stream_adm1_translator = Arc(
        source=m.fs.AD.liquid_outlet, destination=m.fs.translator_adm1_asm2d.inlet
    )
    m.fs.stream_translator_electroNP = Arc(
        source=m.fs.translator_adm1_asm2d.outlet, destination=m.fs.electroNP.inlet
    )
    pyo.TransformationFactory("network.expand_arcs").apply_to(m)

    # Feed conditions based on mass balance in Flores-Alsina, where 0 terms are expressed as 1e-9
    m.fs.AD.inlet.flow_vol[0].fix(
        170 * units.m**3 / units.day
    )  # Double check this value
    m.fs.AD.inlet.temperature.fix(308.15)
    m.fs.AD.inlet.pressure.fix(101325)

    m.fs.AD.inlet.conc_mass_comp[0, "S_su"].fix(0.034597)
    m.fs.AD.inlet.conc_mass_comp[0, "S_aa"].fix(0.015037)
    m.fs.AD.inlet.conc_mass_comp[0, "S_fa"].fix(1e-6)
    m.fs.AD.inlet.conc_mass_comp[0, "S_va"].fix(1e-6)
    m.fs.AD.inlet.conc_mass_comp[0, "S_bu"].fix(1e-6)
    m.fs.AD.inlet.conc_mass_comp[0, "S_pro"].fix(1e-6)
    m.fs.AD.inlet.conc_mass_comp[0, "S_ac"].fix(0.025072)
    m.fs.AD.inlet.conc_mass_comp[0, "S_h2"].fix(1e-6)
    m.fs.AD.inlet.conc_mass_comp[0, "S_ch4"].fix(1e-6)
    m.fs.AD.inlet.conc_mass_comp[0, "S_IC"].fix(0.34628)
    m.fs.AD.inlet.conc_mass_comp[0, "S_IN"].fix(0.60014)
    m.fs.AD.inlet.conc_mass_comp[0, "S_IP"].fix(0.22677)
    m.fs.AD.inlet.conc_mass_comp[0, "S_I"].fix(0.026599)

    m.fs.AD.inlet.conc_mass_comp[0, "X_ch"].fix(7.3687)
    m.fs.AD.inlet.conc_mass_comp[0, "X_pr"].fix(7.7308)
    m.fs.AD.inlet.conc_mass_comp[0, "X_li"].fix(10.3288)
    m.fs.AD.inlet.conc_mass_comp[0, "X_su"].fix(1e-6)
    m.fs.AD.inlet.conc_mass_comp[0, "X_aa"].fix(1e-6)
    m.fs.AD.inlet.conc_mass_comp[0, "X_fa"].fix(1e-6)
    m.fs.AD.inlet.conc_mass_comp[0, "X_c4"].fix(1e-6)
    m.fs.AD.inlet.conc_mass_comp[0, "X_pro"].fix(1e-6)
    m.fs.AD.inlet.conc_mass_comp[0, "X_ac"].fix(1e-6)
    m.fs.AD.inlet.conc_mass_comp[0, "X_h2"].fix(1e-6)
    m.fs.AD.inlet.conc_mass_comp[0, "X_I"].fix(12.7727)
    m.fs.AD.inlet.conc_mass_comp[0, "X_PHA"].fix(0.0022493)
    m.fs.AD.inlet.conc_mass_comp[0, "X_PP"].fix(1.04110)
    m.fs.AD.inlet.conc_mass_comp[0, "X_PAO"].fix(3.4655)
    m.fs.AD.inlet.conc_mass_comp[0, "S_K"].fix(0.02268)
    m.fs.AD.inlet.conc_mass_comp[0, "S_Mg"].fix(0.02893)

    m.fs.AD.inlet.cations[0].fix(0.04)
    m.fs.AD.inlet.anions[0].fix(0.02)

    m.fs.AD.volume_liquid.fix(3400)
    m.fs.AD.volume_vapor.fix(300)
    m.fs.AD.liquid_outlet.temperature.fix(308.15)

    # ElectroNP
    m.fs.electroNP.energy_electric_flow_mass.fix(0.044 * units.kWh / units.kg)
    m.fs.electroNP.magnesium_chloride_dosage.fix(0.388)

    # Costing
    m.fs.costing.electroNP.phosphorus_recovery_value = 0

    # scaling
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
        if "conc_mass_comp[S_IN]" in var.name:
            iscale.set_scaling_factor(var, 1e-1)
        if "conc_mass_comp[S_IP]" in var.name:
            iscale.set_scaling_factor(var, 1e-1)
        if "conc_mass_comp[S_PO4]" in var.name:
            iscale.set_scaling_factor(var, 1e-1)
        if "conc_mass_comp[S_NH4]" in var.name:
            iscale.set_scaling_factor(var, 1e1)
        if "conc_mass_comp[S_F]" in var.name:
            iscale.set_scaling_factor(var, 1e-1)
        if "conc_mass_comp[X_I]" in var.name:
            iscale.set_scaling_factor(var, 1e-1)

    iscale.calculate_scaling_factors(m)

    iscale.set_scaling_factor(m.fs.electroNP.byproduct.flow_vol[0.0], 1e7)
    iscale.set_scaling_factor(m.fs.AD.vapor_phase[0].pressure_sat, 1e-3)

    m.fs.AD.initialize(outlvl=idaeslog.INFO_HIGH)
    propagate_state(m.fs.stream_adm1_translator)
    m.fs.translator_adm1_asm2d.initialize(outlvl=idaeslog.INFO_HIGH)
    propagate_state(m.fs.stream_translator_electroNP)
    m.fs.electroNP.initialize(outlvl=idaeslog.INFO_HIGH)
    m.fs.costing.initialize()

    results = solve(m, tee=True)
    return m, results


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
            to_units=units.kWh / units.m**3,
        )
    )
    print(f"Electricity Intensity: {electricity_intensity:.4f} kWh/m3")

    print("\n----------Levelized Cost----------")
    LCOW = value(
        pyunits.convert(m.fs.costing.LCOW, to_units=pyunits.USD_2018 / pyunits.m**3)
    )
    print(f"Levelized Cost of Water: {LCOW:.3f} $/m^3")


if __name__ == "__main__":
    m, results = build_flowsheet()
    assert_optimal_termination(results)
    stream_table = create_stream_table_dataframe(
        {
            "AD inlet": m.fs.AD.inlet,
            "AD liquid outlet": m.fs.AD.liquid_outlet,
            "AD vapor outlet": m.fs.AD.vapor_outlet,
            "Translator outlet": m.fs.translator_adm1_asm2d.outlet,
            "ElectroNP treated": m.fs.electroNP.treated,
            "ElectroNP byproduct": m.fs.electroNP.byproduct,
        },
        time_point=0,
    )
    print(stream_table_dataframe_to_string(stream_table))
    display_costing(m)
