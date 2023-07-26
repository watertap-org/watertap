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
    # Activated Sludge Process
    # First reactor (anoxic) - standard CSTR
    m.fs.R1 = CSTR(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    # First reactor (anoxic) - standard CSTR
    m.fs.R2 = CSTR(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    # Second reactor (anoxic) - standard CSTR
    m.fs.R3 = CSTR(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    m.fs.R4 = CSTR(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    # Third reactor (aerobic) - CSTR with injection
    m.fs.R5 = CSTR_Injection(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    # Fourth reactor (aerobic) - CSTR with injection
    m.fs.R6 = CSTR_Injection(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    # Fifth reactor (aerobic) - CSTR with injection
    m.fs.R7 = CSTR_Injection(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    m.fs.SP1 = Separator(
        property_package=m.fs.props_ASM2D, outlet_list=["underflow", "overflow"]
    )

    # Clarifier
    # TODO: Replace with more detailed model when available
    m.fs.CL1 = Separator(
        property_package=m.fs.props_ASM2D,
        outlet_list=["underflow", "effluent"],
        split_basis=SplittingType.componentFlow,
    )
    # Mixing sludge recycle and R5 underflow
    m.fs.MX2 = Mixer(
        property_package=m.fs.props_ASM2D, inlet_list=["clarifier", "reactor"]
    )
    # Sludge separator
    m.fs.SP2 = Separator(
        property_package=m.fs.props_ASM2D, outlet_list=["waste", "recycle"]
    )

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
    m.fs.stream1 = Arc(source=m.fs.feed.outlet, destination=m.fs.R1.inlet)
    m.fs.stream3 = Arc(source=m.fs.R1.outlet, destination=m.fs.R2.inlet)
    m.fs.stream4 = Arc(source=m.fs.R2.outlet, destination=m.fs.MX2.reactor)
    m.fs.stream5 = Arc(source=m.fs.MX2.outlet, destination=m.fs.R3.inlet)
    m.fs.stream6 = Arc(source=m.fs.R3.outlet, destination=m.fs.R4.inlet)
    m.fs.stream7 = Arc(source=m.fs.R4.outlet, destination=m.fs.R5.inlet)
    m.fs.stream8 = Arc(source=m.fs.R5.outlet, destination=m.fs.R6.inlet)
    m.fs.stream9 = Arc(source=m.fs.R6.outlet, destination=m.fs.R7.inlet)
    m.fs.stream10 = Arc(source=m.fs.R7.outlet, destination=m.fs.SP1.inlet)
    m.fs.stream11 = Arc(source=m.fs.SP1.overflow, destination=m.fs.CL1.inlet)
    m.fs.stream12 = Arc(source=m.fs.SP1.underflow, destination=m.fs.MX2.clarifier)
    # m.fs.stream13 = Arc(source=m.fs.CL1.effluent, destination=m.fs.Treated.inlet)
    m.fs.stream14 = Arc(source=m.fs.CL1.underflow, destination=m.fs.thickener.inlet)
    # m.fs.stream_feed_thickener = Arc(
    #     source=m.fs.feed.outlet, destination=m.fs.thickener.inlet
    # )
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

    # Oxygen concentration in reactors 3 and 4 is governed by mass transfer
    # Add additional parameter and constraints
    m.fs.R5.KLa = pyo.Var(
        initialize=240,
        units=pyo.units.hour**-1,
        doc="Lumped mass transfer coefficient for oxygen",
    )
    m.fs.R6.KLa = pyo.Var(
        initialize=240,
        units=pyo.units.hour**-1,
        doc="Lumped mass transfer coefficient for oxygen",
    )
    m.fs.R7.KLa = pyo.Var(
        initialize=84,
        units=pyo.units.hour**-1,
        doc="Lumped mass transfer coefficient for oxygen",
    )
    m.fs.S_O_eq = pyo.Param(
        default=8e-3,
        units=pyo.units.kg / pyo.units.m**3,
        mutable=True,
        doc="Dissolved oxygen concentration at equilibrium",
    )

    @m.fs.R5.Constraint(m.fs.time, doc="Mass transfer constraint for R3")
    def mass_transfer_R5(self, t):
        return pyo.units.convert(
            m.fs.R5.injection[t, "Liq", "S_O2"], to_units=pyo.units.kg / pyo.units.hour
        ) == (
            m.fs.R5.KLa
            * m.fs.R5.volume[t]
            * (m.fs.S_O_eq - m.fs.R5.outlet.conc_mass_comp[t, "S_O2"])
        )

    @m.fs.R6.Constraint(m.fs.time, doc="Mass transfer constraint for R4")
    def mass_transfer_R6(self, t):
        return pyo.units.convert(
            m.fs.R6.injection[t, "Liq", "S_O2"], to_units=pyo.units.kg / pyo.units.hour
        ) == (
            m.fs.R6.KLa
            * m.fs.R6.volume[t]
            * (m.fs.S_O_eq - m.fs.R6.outlet.conc_mass_comp[t, "S_O2"])
        )

    @m.fs.R7.Constraint(m.fs.time, doc="Mass transfer constraint for R4")
    def mass_transfer_R7(self, t):
        return pyo.units.convert(
            m.fs.R7.injection[t, "Liq", "S_O2"], to_units=pyo.units.kg / pyo.units.hour
        ) == (
            m.fs.R7.KLa
            * m.fs.R7.volume[t]
            * (m.fs.S_O_eq - m.fs.R7.outlet.conc_mass_comp[t, "S_O2"])
        )

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
    m.fs.feed.properties[0].pressure.fix(1 * units.atm)
    m.fs.feed.properties[0].temperature.fix(308.15 * units.K)
    m.fs.feed.properties[0].flow_vol.fix(178.4674 * units.m**3 / units.day)
    eps = 1e-9 * units.kg / units.m**3

    m.fs.feed.properties[0].conc_mass_comp["S_O2"].fix(eps)
    m.fs.feed.properties[0].conc_mass_comp["S_F"].fix(0.02644 * units.kg / units.m**3)
    m.fs.feed.properties[0].conc_mass_comp["S_A"].fix(0.01766 * units.kg / units.m**3)
    m.fs.feed.properties[0].conc_mass_comp["S_I"].fix(0.02723 * units.kg / units.m**3)
    m.fs.feed.properties[0].conc_mass_comp["S_NH4"].fix(
        0.01858 * units.kg / units.m**3
    )
    m.fs.feed.properties[0].conc_mass_comp["S_N2"].fix(
        0.00507 * units.kg / units.m**3
    )
    m.fs.feed.properties[0].conc_mass_comp["S_NO3"].fix(
        0.00002 * units.kg / units.m**3
    )
    m.fs.feed.properties[0].conc_mass_comp["S_PO4"].fix(
        0.00469 * units.kg / units.m**3
    )
    m.fs.feed.properties[0].conc_mass_comp["S_IC"].fix(
        0.07899 * units.kg / units.m**3
    )

    m.fs.feed.properties[0].conc_mass_comp["X_I"].fix(
        10.96441 * units.kg / units.m**3
    )
    m.fs.feed.properties[0].conc_mass_comp["X_S"].fix(
        19.08476 * units.kg / units.m**3
    )
    m.fs.feed.properties[0].conc_mass_comp["X_H"].fix(9.47939 * units.kg / units.m**3)
    m.fs.feed.properties[0].conc_mass_comp["X_PAO"].fix(
        3.8622 * units.kg / units.m**3
    )
    m.fs.feed.properties[0].conc_mass_comp["X_PP"].fix(
        0.45087 * units.kg / units.m**3
    )
    m.fs.feed.properties[0].conc_mass_comp["X_PHA"].fix(
        0.02464 * units.kg / units.m**3
    )
    m.fs.feed.properties[0].conc_mass_comp["X_AUT"].fix(
        0.33379 * units.kg / units.m**3
    )
    m.fs.feed.properties[0].conc_mass_comp["S_K"].fix(0.01979 * units.kg / units.m**3)
    m.fs.feed.properties[0].conc_mass_comp["S_Mg"].fix(
        0.18987 * units.kg / units.m**3
    )

    print("DOF before AD:", degrees_of_freedom(m))

    # Activated Sludge
    # Reactor sizing
    m.fs.R1.volume.fix(1000 * pyo.units.m**3)
    m.fs.R2.volume.fix(1000 * pyo.units.m**3)
    m.fs.R3.volume.fix(1000 * pyo.units.m**3)
    m.fs.R4.volume.fix(1000 * pyo.units.m**3)
    m.fs.R5.volume.fix(1333 * pyo.units.m**3)
    m.fs.R6.volume.fix(1333 * pyo.units.m**3)
    m.fs.R7.volume.fix(1333 * pyo.units.m**3)

    # Injection rates to Reactions 3, 4 and 5
    for j in m.fs.props_ASM2D.component_list:
        if j != "S_O2":
            # All components except S_O have no injection
            m.fs.R5.injection[:, :, j].fix(0)
            m.fs.R6.injection[:, :, j].fix(0)
            m.fs.R7.injection[:, :, j].fix(0)
    # Then set injections rates for O2
    m.fs.R5.outlet.conc_mass_comp[:, "S_O2"].fix(1.91e-3)
    m.fs.R6.outlet.conc_mass_comp[:, "S_O2"].fix(2.60e-3)
    m.fs.R7.outlet.conc_mass_comp[:, "S_O2"].fix(3.20e-3)

    # Set fraction of outflow from reactor 5 that goes to recycle
    m.fs.SP1.split_fraction[:, "underflow"].fix(0.60)

    # Clarifier
    # TODO: Update once more detailed model available
    m.fs.CL1.split_fraction[0, "effluent", "H2O"].fix(0.49986)
    m.fs.CL1.split_fraction[0, "effluent", "S_A"].fix(0.49986)
    m.fs.CL1.split_fraction[0, "effluent", "S_F"].fix(0.49986)
    m.fs.CL1.split_fraction[0, "effluent", "S_I"].fix(0.49986)
    m.fs.CL1.split_fraction[0, "effluent", "S_N2"].fix(0.49986)
    m.fs.CL1.split_fraction[0, "effluent", "S_NH4"].fix(0.49986)
    m.fs.CL1.split_fraction[0, "effluent", "S_NO3"].fix(0.49986)
    m.fs.CL1.split_fraction[0, "effluent", "S_O2"].fix(0.49986)
    m.fs.CL1.split_fraction[0, "effluent", "S_PO4"].fix(0.49986)
    m.fs.CL1.split_fraction[0, "effluent", "S_IC"].fix(0.49986)
    m.fs.CL1.split_fraction[0, "effluent", "S_K"].fix(0.49986)
    m.fs.CL1.split_fraction[0, "effluent", "S_Mg"].fix(0.49986)
    # m.fs.CL1.split_fraction[0, "effluent", "S_ALK"].fix(0.49986)
    m.fs.CL1.split_fraction[0, "effluent", "X_AUT"].fix(0.022117)
    m.fs.CL1.split_fraction[0, "effluent", "X_H"].fix(0.021922)
    m.fs.CL1.split_fraction[0, "effluent", "X_I"].fix(0.021715)
    # m.fs.CL1.split_fraction[0, "effluent", "X_MeOH"].fix(0.022)
    # m.fs.CL1.split_fraction[0, "effluent", "X_MeP"].fix(0.022)
    m.fs.CL1.split_fraction[0, "effluent", "X_PAO"].fix(0.022)
    m.fs.CL1.split_fraction[0, "effluent", "X_PHA"].fix(0.02147)
    m.fs.CL1.split_fraction[0, "effluent", "X_PP"].fix(0.02144)
    m.fs.CL1.split_fraction[0, "effluent", "X_S"].fix(0.02221)
    # m.fs.CL1.split_fraction[0, "effluent", "X_TSS"].fix(0.02194)

    # Sludge purge separator
    m.fs.SP2.split_fraction[:, "recycle"].fix(0.97955)

    # AD
    m.fs.AD.volume_liquid.fix(3400)
    m.fs.AD.volume_vapor.fix(300)
    m.fs.AD.liquid_outlet.temperature.fix(308.15)

    # ElectroNP
    m.fs.electroNP.energy_electric_flow_mass.fix(0.044 * units.kWh / units.kg)
    m.fs.electroNP.magnesium_chloride_dosage.fix(0.388)

    # Costing
    m.fs.costing.electroNP.phosphorus_recovery_value = 0


def initialize_system(m):
    # Initialize
    m.fs.feed.initialize()
    propagate_state(m.fs.stream1)
    m.fs.R1.initialize()
    propagate_state(m.fs.stream3)
    m.fs.R2.initialize()
    propagate_state(m.fs.stream4)
    propagate_state(m.fs.stream5)
    m.fs.R3.initialize()
    propagate_state(m.fs.stream6)
    m.fs.R4.initialize()
    propagate_state(m.fs.stream7)
    m.fs.R5.initialize()
    propagate_state(m.fs.stream8)
    m.fs.R6.initialize()
    propagate_state(m.fs.stream9)
    m.fs.R7.initialize()
    propagate_state(m.fs.stream10)
    m.fs.SP1.initialize()
    propagate_state(m.fs.stream14)
    # propagate_state(m.fs.stream_feed_thickener)
    m.fs.thickener.initialize()
    propagate_state(m.fs.stream_thickener_translator)
    m.fs.translator_asm2d_adm1.initialize(outlvl=idaeslog.INFO_HIGH)
    propagate_state(m.fs.stream_translator_AD)
    m.fs.AD.initialize(outlvl=idaeslog.INFO_HIGH)
    propagate_state((m.fs.stream_AD_translator))
    m.fs.translator_adm1_asm2d.initialize(outlvl=idaeslog.INFO_HIGH)
    propagate_state(m.fs.stream_translator_dewater)
    m.fs.dewater.initialize(outlvl=idaeslog.INFO_HIGH)
    propagate_state(m.fs.stream_dewater_electroNP)
    m.fs.electroNP.initialize(outlvl=idaeslog.INFO_HIGH)
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
    m = main()
