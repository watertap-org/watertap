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

"""

# Some more information about this module
__author__ = "Alejandro Garciadiego, Andrew Lee"

import pyomo.environ as pyo
from pyomo.environ import (
    value,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.network import Arc, SequentialDecomposition

from idaes.core import FlowsheetBlock
from idaes.models.unit_models import (
    CSTR,
    Feed,
    Mixer,
    Separator,
    PressureChanger,
    Product,
)
from idaes.models.unit_models.separator import SplittingType
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
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
from watertap.core.util.initialization import check_solve

# Set up logger
_log = idaeslog.getLogger(__name__)


def automate_rescale_variables(m):
    for var, sv in iscale.badly_scaled_var_generator(m):
        if iscale.get_scaling_factor(var) is None:
            continue
        sf = iscale.get_scaling_factor(var)
        iscale.set_scaling_factor(var, sf / sv)
        iscale.calculate_scaling_factors(m)


def build_flowsheet():
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

    # Feed water stream
    m.fs.feed = Feed(property_package=m.fs.props_ASM2D)
    # Mixer for feed water and recycled sludge
    # m.fs.MX1 = Mixer(property_package=m.fs.props_ASM2D, inlet_list=["feed_water", "recycle"])
    # First reactor (anoxic) - standard CSTR
    m.fs.R1 = CSTR(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    # Second reactor (anoxic) - standard CSTR
    m.fs.R2 = CSTR(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    # Third reactor (aerobic) - CSTR with injection
    m.fs.R3 = CSTR_Injection(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    # Fourth reactor (aerobic) - CSTR with injection
    m.fs.R4 = CSTR_Injection(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    # Fifth reactor (aerobic) - CSTR with injection
    m.fs.R5 = CSTR_Injection(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    m.fs.SP5 = Separator(
        property_package=m.fs.props_ASM2D, outlet_list=["underflow", "overflow"]
    )
    # Clarifier
    # TODO: Replace with more detailed model when available
    m.fs.CL1 = Separator(
        property_package=m.fs.props_ASM2D,
        outlet_list=["underflow", "effluent"],
        split_basis=SplittingType.componentFlow,
    )
    # Sludge purge splitter
    m.fs.SP6 = Separator(
        property_package=m.fs.props_ASM2D,
        outlet_list=["recycle", "waste"],
        split_basis=SplittingType.totalFlow,
    )
    # Mixing sludge recycle and R5 underflow
    # m.fs.MX6 = Mixer(property_package=m.fs.props_ASM2D, inlet_list=["clarifier", "reactor"])
    # Product Blocks
    m.fs.Treated = Product(property_package=m.fs.props_ASM2D)
    m.fs.Sludge = Product(property_package=m.fs.props_ASM2D)
    # Recycle pressure changer - use a simple isothermal unit for now
    # m.fs.P1 = PressureChanger(property_package=m.fs.props_ASM2D)

    # Thickener
    m.fs.thickener = Thickener(
        property_package=m.fs.props_ASM2D,
        activated_sludge_model=thickener_type.modified_ASM2D,
    )

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

    # Anaerobic digestor
    m.fs.AD = AD(
        liquid_property_package=m.fs.props_ADM1,
        vapor_property_package=m.fs.props_vap_ADM1,
        reaction_package=m.fs.rxn_props_ADM1,
        has_heat_transfer=True,
        has_pressure_change=False,
    )

    # Dewatering Unit
    m.fs.dewater = DewateringUnit(
        property_package=m.fs.props_ASM2D,
        activated_sludge_model=dewater_type.modified_ASM2D,
    )

    # ElectroNP
    m.fs.electroNP = ElectroNPZO(property_package=m.fs.props_ASM2D)
    # m.fs.electroNP.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    # Link units
    m.fs.stream1 = Arc(source=m.fs.feed.outlet, destination=m.fs.R1.inlet)
    m.fs.stream3 = Arc(source=m.fs.R1.outlet, destination=m.fs.R2.inlet)
    m.fs.stream4 = Arc(source=m.fs.R2.outlet, destination=m.fs.R3.inlet)
    m.fs.stream5 = Arc(source=m.fs.R3.outlet, destination=m.fs.R4.inlet)
    m.fs.stream6 = Arc(source=m.fs.R4.outlet, destination=m.fs.R5.inlet)
    m.fs.stream7 = Arc(source=m.fs.R5.outlet, destination=m.fs.SP5.inlet)
    m.fs.stream8 = Arc(source=m.fs.SP5.overflow, destination=m.fs.CL1.inlet)
    m.fs.stream10 = Arc(source=m.fs.CL1.effluent, destination=m.fs.Treated.inlet)
    m.fs.stream11 = Arc(source=m.fs.CL1.underflow, destination=m.fs.SP6.inlet)

    m.fs.stream_SP_thickener = Arc(
        source=m.fs.SP6.waste, destination=m.fs.thickener.inlet
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

    # m.fs.stream1 = Arc(source=m.fs.feed.outlet, destination=m.fs.MX1.feed_water)
    # m.fs.stream2 = Arc(source=m.fs.MX1.outlet, destination=m.fs.R1.inlet)
    # m.fs.stream3 = Arc(source=m.fs.R1.outlet, destination=m.fs.R2.inlet)
    # m.fs.stream4 = Arc(source=m.fs.R2.outlet, destination=m.fs.R3.inlet)
    # m.fs.stream5 = Arc(source=m.fs.R3.outlet, destination=m.fs.R4.inlet)
    # m.fs.stream6 = Arc(source=m.fs.R4.outlet, destination=m.fs.R5.inlet)
    # m.fs.stream7 = Arc(source=m.fs.R5.outlet, destination=m.fs.SP5.inlet)
    # m.fs.stream8 = Arc(source=m.fs.SP5.overflow, destination=m.fs.CL1.inlet)
    # m.fs.stream9 = Arc(source=m.fs.SP5.underflow, destination=m.fs.MX6.reactor)
    # m.fs.stream10 = Arc(source=m.fs.CL1.effluent, destination=m.fs.Treated.inlet)
    # m.fs.stream11 = Arc(source=m.fs.CL1.underflow, destination=m.fs.SP6.inlet)
    # m.fs.stream102 = Arc(source=m.fs.SP6.waste, destination=m.fs.Sludge.inlet)
    # m.fs.stream13 = Arc(source=m.fs.SP6.recycle, destination=m.fs.MX6.clarifier)
    # m.fs.stream14 = Arc(source=m.fs.MX6.outlet, destination=m.fs.P1.inlet)
    # m.fs.stream15 = Arc(source=m.fs.P1.outlet, destination=m.fs.MX1.recycle)
    pyo.TransformationFactory("network.expand_arcs").apply_to(m)

    # Oxygen concentration in reactors 3 and 4 is governed by mass transfer
    # Add additional parameter and constraints
    m.fs.R3.KLa = pyo.Var(
        initialize=3.4,
        units=pyo.units.hour**-1,
        doc="Lumped mass transfer coefficient for oxygen",
    )
    m.fs.R4.KLa = pyo.Var(
        initialize=1.8,
        units=pyo.units.hour**-1,
        doc="Lumped mass transfer coefficient for oxygen",
    )
    m.fs.S_O_eq = pyo.Param(
        default=8e-3,
        units=pyo.units.kg / pyo.units.m**3,
        mutable=True,
        doc="Dissolved oxygen concentration at equilibrium",
    )

    m.fs.R3.Constraint(m.fs.time, doc="Mass transfer constraint for R3")

    def mass_transfer_R3(self, t):
        return pyo.units.convert(
            m.fs.R3.injection[t, "Liq", "S_O2"], to_units=pyo.units.kg / pyo.units.hour
        ) == (
            m.fs.R3.KLa
            * m.fs.R3.volume[t]
            * (m.fs.S_O_eq - m.fs.R3.outlet.conc_mass_comp[t, "S_O2"])
        )

    m.fs.R4.Constraint(m.fs.time, doc="Mass transfer constraint for R4")

    def mass_transfer_R4(self, t):
        return pyo.units.convert(
            m.fs.R4.injection[t, "Liq", "S_O2"], to_units=pyo.units.kg / pyo.units.hour
        ) == (
            m.fs.R4.KLa
            * m.fs.R4.volume[t]
            * (m.fs.S_O_eq - m.fs.R4.outlet.conc_mass_comp[t, "S_O2"])
        )

    # Feed inlet
    m.fs.feed.flow_vol.fix(18446 * pyo.units.m**3 / pyo.units.day)
    m.fs.feed.temperature.fix(298.15 * pyo.units.K)
    m.fs.feed.pressure.fix(1 * pyo.units.atm)
    m.fs.feed.conc_mass_comp[0, "S_O2"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "S_F"].fix(30 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "S_A"].fix(20 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "S_NH4"].fix(16 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "S_NO3"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "S_PO4"].fix(3.6 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "S_I"].fix(30 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "S_N2"].fix(15 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "X_I"].fix(25 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "X_S"].fix(125 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "X_H"].fix(30 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "X_PAO"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "X_PP"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "X_PHA"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "X_AUT"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.properties[0].conc_mass_comp["S_IC"].fix(
        0.07899 * pyunits.kg / pyunits.m**3
    )
    m.fs.feed.properties[0].conc_mass_comp["S_K"].fix(
        1e-6 * pyo.units.g / pyo.units.m**3
    )
    m.fs.feed.properties[0].conc_mass_comp["S_Mg"].fix(
        1e-6 * pyo.units.g / pyo.units.m**3
    )

    # m.fs.feed.properties[0].pressure.fix(1 * pyunits.atm)
    # m.fs.feed.properties[0].temperature.fix(308.15 * pyunits.K)
    # m.fs.feed.properties[0].flow_vol.fix(178.4674 * pyunits.m ** 3 / pyunits.day)

    # eps = 1e-9 * pyunits.kg / pyunits.m ** 3
    # m.fs.feed.properties[0].conc_mass_comp["S_O2"].fix(eps)
    # m.fs.feed.properties[0].conc_mass_comp["S_F"].fix(0.02644 * pyunits.kg / pyunits.m ** 3)
    # m.fs.feed.properties[0].conc_mass_comp["S_A"].fix(0.01766 * pyunits.kg / pyunits.m ** 3)
    # m.fs.feed.properties[0].conc_mass_comp["S_I"].fix(0.02723 * pyunits.kg / pyunits.m ** 3)
    # m.fs.feed.properties[0].conc_mass_comp["S_NH4"].fix(
    #     0.01858 * pyunits.kg / pyunits.m ** 3
    # )
    # m.fs.feed.properties[0].conc_mass_comp["S_N2"].fix(
    #     0.00507 * pyunits.kg / pyunits.m ** 3
    # )
    # m.fs.feed.properties[0].conc_mass_comp["S_NO3"].fix(
    #     0.00002 * pyunits.kg / pyunits.m ** 3
    # )
    # m.fs.feed.properties[0].conc_mass_comp["S_PO4"].fix(
    #     0.00469 * pyunits.kg / pyunits.m ** 3
    # )
    # m.fs.feed.properties[0].conc_mass_comp["S_IC"].fix(
    #     0.07899 * pyunits.kg / pyunits.m ** 3
    # )
    #
    # m.fs.feed.properties[0].conc_mass_comp["X_I"].fix(
    #     10.96441 * pyunits.kg / pyunits.m ** 3
    # )
    # m.fs.feed.properties[0].conc_mass_comp["X_S"].fix(
    #     19.08476 * pyunits.kg / pyunits.m ** 3
    # )
    # m.fs.feed.properties[0].conc_mass_comp["X_H"].fix(9.47939 * pyunits.kg / pyunits.m ** 3)
    # m.fs.feed.properties[0].conc_mass_comp["X_PAO"].fix(
    #     3.8622 * pyunits.kg / pyunits.m ** 3
    # )
    # m.fs.feed.properties[0].conc_mass_comp["X_PP"].fix(
    #     0.45087 * pyunits.kg / pyunits.m ** 3
    # )
    # m.fs.feed.properties[0].conc_mass_comp["X_PHA"].fix(
    #     0.02464 * pyunits.kg / pyunits.m ** 3
    # )
    # m.fs.feed.properties[0].conc_mass_comp["X_AUT"].fix(
    #     0.33379 * pyunits.kg / pyunits.m ** 3
    # )
    # m.fs.feed.properties[0].conc_mass_comp["S_K"].fix(0.01979 * pyunits.kg / pyunits.m ** 3)
    # m.fs.feed.properties[0].conc_mass_comp["S_Mg"].fix(
    #     0.18987 * pyunits.kg / pyunits.m ** 3
    # )

    # Reactor sizing
    m.fs.R1.volume.fix(1000 * pyo.units.m**3)
    m.fs.R2.volume.fix(1000 * pyo.units.m**3)
    m.fs.R3.volume.fix(1333 * pyo.units.m**3)
    m.fs.R4.volume.fix(1333 * pyo.units.m**3)
    m.fs.R5.volume.fix(1333 * pyo.units.m**3)

    # Injection rates to Reactions 3, 4 and 5
    for j in m.fs.props_ASM2D.component_list:
        if j != "S_O2":
            # All components except S_O have no injection
            m.fs.R3.injection[:, :, j].fix(0)
            m.fs.R4.injection[:, :, j].fix(0)
            m.fs.R5.injection[:, :, j].fix(0)
    # Then set injections rates for O2
    m.fs.R3.outlet.conc_mass_comp[:, "S_O2"].fix(2.72e-3)
    m.fs.R4.outlet.conc_mass_comp[:, "S_O2"].fix(2.43e-3)
    m.fs.R5.outlet.conc_mass_comp[:, "S_O2"].fix(4.49e-4)

    # Set fraction of outflow from reactor 5 that goes to recycle
    m.fs.SP5.split_fraction[:, "underflow"].fix(0.1)

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
    m.fs.CL1.split_fraction[0, "effluent", "X_AUT"].fix(0.022117)
    m.fs.CL1.split_fraction[0, "effluent", "X_H"].fix(0.021922)
    m.fs.CL1.split_fraction[0, "effluent", "X_I"].fix(0.021715)
    m.fs.CL1.split_fraction[0, "effluent", "X_PAO"].fix(0.022)
    m.fs.CL1.split_fraction[0, "effluent", "X_PHA"].fix(0.02147)
    m.fs.CL1.split_fraction[0, "effluent", "X_PP"].fix(0.02144)
    m.fs.CL1.split_fraction[0, "effluent", "X_S"].fix(0.02221)

    # Sludge purge separator
    m.fs.SP6.split_fraction[:, "recycle"].fix(0.65)

    # Outlet pressure from recycle pump
    # m.fs.P1.outlet.pressure.fix(101325)

    # AD
    m.fs.AD.volume_liquid.fix(3400)
    m.fs.AD.volume_vapor.fix(300)
    m.fs.AD.liquid_outlet.temperature.fix(298.15)

    # ElectroNP
    m.fs.electroNP.energy_electric_flow_mass.fix(0.044 * pyunits.kWh / pyunits.kg)
    m.fs.electroNP.magnesium_chloride_dosage.fix(0.388)

    # Check degrees of freedom
    print(degrees_of_freedom(m))
    assert degrees_of_freedom(m) == 0

    def scale_variables(m):
        for var in m.fs.component_data_objects(pyo.Var, descend_into=True):
            if "flow_vol" in var.name:
                iscale.set_scaling_factor(var, 1e2)
            if "temperature" in var.name:
                iscale.set_scaling_factor(var, 1e-1)
            if "pressure" in var.name:
                iscale.set_scaling_factor(var, 1e-4)
            if "conc_mass_comp" in var.name:
                iscale.set_scaling_factor(var, 1e1)
            if "enth_mol" in var.name:
                iscale.set_scaling_factor(var, 1e-4)
            if "conc_mass_comp[S_F]" in var.name:
                iscale.set_scaling_factor(var, 1e2)
            if "conc_mass_comp[S_A]" in var.name:
                iscale.set_scaling_factor(var, 1e2)
            if "conc_mass_comp[S_NH4]" in var.name:
                iscale.set_scaling_factor(var, 1e2)
            if "conc_mass_comp[S_PO4]" in var.name:
                iscale.set_scaling_factor(var, 1e3)
            if "conc_mass_comp[S_I]" in var.name:
                iscale.set_scaling_factor(var, 1e2)
            if "conc_mass_comp[S_N2]" in var.name:
                iscale.set_scaling_factor(var, 1e2)
            if "conc_mass_comp[X_I]" in var.name:
                iscale.set_scaling_factor(var, 1e2)
            if "conc_mass_comp[X_S]" in var.name:
                iscale.set_scaling_factor(var, 1e1)
            if "conc_mass_comp[X_H]" in var.name:
                iscale.set_scaling_factor(var, 1e2)
            if "conc_mass_comp[S_IC]" in var.name:
                iscale.set_scaling_factor(var, 1e2)
            if "conc_mass_comp[S_O2]" in var.name:
                iscale.set_scaling_factor(var, 1e5)
            if "conc_mass_comp[S_NO3]" in var.name:
                iscale.set_scaling_factor(var, 1e5)
            if "conc_mass_comp[X_PAO]" in var.name:
                iscale.set_scaling_factor(var, 1e5)
            if "conc_mass_comp[X_PP]" in var.name:
                iscale.set_scaling_factor(var, 1e5)
            if "conc_mass_comp[X_PHA]" in var.name:
                iscale.set_scaling_factor(var, 1e5)
            if "conc_mass_comp[X_AUT]" in var.name:
                iscale.set_scaling_factor(var, 1e5)
            if "conc_mass_comp[S_K]" in var.name:
                iscale.set_scaling_factor(var, 1e5)
            if "conc_mass_comp[S_Mg]" in var.name:
                iscale.set_scaling_factor(var, 1e5)

            # if "conc_mass_comp[S_su]" in var.name:
            #     iscale.set_scaling_factor(var, 1e3)
            # if "conc_mass_comp[S_aa]" in var.name:
            #     iscale.set_scaling_factor(var, 1e3)
            # if "conc_mass_comp[S_fa]" in var.name:
            #     iscale.set_scaling_factor(var, 1e0)
            # if "conc_mass_comp[S_va]" in var.name:
            #     iscale.set_scaling_factor(var, 1e0)
            # if "conc_mass_comp[S_bu]" in var.name:
            #     iscale.set_scaling_factor(var, 1e0)
            # if "conc_mass_comp[S_pro]" in var.name:
            #     iscale.set_scaling_factor(var, 1e0)
            # if "conc_mass_comp[S_ac]" in var.name:
            #     iscale.set_scaling_factor(var, 1e2)
            # if "conc_mass_comp[S_h2]" in var.name:
            #     iscale.set_scaling_factor(var, 1e3)
            # if "conc_mass_comp[S_ch4]" in var.name:
            #     iscale.set_scaling_factor(var, 1e2)
            # if "conc_mass_comp[S_IN]" in var.name:
            #     iscale.set_scaling_factor(var, 1e0)
            # if "conc_mass_comp[S_IP]" in var.name:
            #     iscale.set_scaling_factor(var, 1e1)
            # if "conc_mass_comp[X_ch]" in var.name:
            #     iscale.set_scaling_factor(var, 1e3)
            # if "conc_mass_comp[X_pr]" in var.name:
            #     iscale.set_scaling_factor(var, 1e3)
            # if "conc_mass_comp[X_li]" in var.name:
            #     iscale.set_scaling_factor(var, 1e2)
            # if "conc_mass_comp[X_su]" in var.name:
            #     iscale.set_scaling_factor(var, 1e1)
            # if "conc_mass_comp[X_aa]" in var.name:
            #     iscale.set_scaling_factor(var, 1e1)
            # if "conc_mass_comp[X_fa]" in var.name:
            #     iscale.set_scaling_factor(var, 1e4)
            # if "conc_mass_comp[X_c4]" in var.name:
            #     iscale.set_scaling_factor(var, 1e4)
            # if "conc_mass_comp[X_pro]" in var.name:
            #     iscale.set_scaling_factor(var, 1e4)
            # if "conc_mass_comp[X_ac]" in var.name:
            #     iscale.set_scaling_factor(var, 1e1)
            # if "conc_mass_comp[X_h2]" in var.name:
            #     iscale.set_scaling_factor(var, 1e5)

    scale_variables(m)
    iscale.calculate_scaling_factors(m.fs)

    # Initialize flowsheet
    # Apply sequential decomposition - 1 iteration should suffice
    seq = SequentialDecomposition()
    seq.options.select_tear_method = "heuristic"
    seq.options.tear_method = "Wegstein"
    seq.options.iterLim = 0

    G = seq.create_graph(m)

    # # Uncomment this code to see tear set and initialization order
    # heuristic_tear_set = seq.tear_set_arcs(G, method="heuristic")
    # order = seq.calculation_order(G)
    # for o in heuristic_tear_set:
    #     print(o.name)
    # for o in order:
    #     print(o[0].name)

    # Initial guesses for flow into first reactor
    # tear_guesses = {
    #     "flow_vol": {0: 0.31},
    #     "conc_mass_comp": {
    #         (0, "S_O2"): 0.14,
    #         (0, "S_F"): 0.021,
    #         (0, "S_A"): 0.014,
    #         (0, "S_NH4"): 0.016,
    #         (0, "S_NO3"): 1e-6,
    #         (0, "S_PO4"): 0.0035,
    #         (0, "S_I"): 0.030,
    #         (0, "S_IC"): 0.030,
    #         (0, "S_K"): 0.030,
    #         (0, "S_Mg"): 0.030,
    #         (0, "S_N2"): 0.015,
    #         (0, "X_I"): 0.050,
    #         (0, "X_S"): 0.138,
    #         (0, "X_H"): 0.014,
    #         (0, "X_PAO"): 1e-6,
    #         (0, "X_PP"): 1e-6,
    #         (0, "X_PHA"): 1e-6,
    #         (0, "X_AUT"): 1e-6,
    #     },
    #     "temperature": {0: 298.15},
    #     "pressure": {0: 101325},
    # }
    #
    # # Pass the tear_guess to the SD tool
    # seq.set_guesses_for(m.fs.R1.inlet, tear_guesses)

    def function(unit):
        unit.initialize(outlvl=idaeslog.INFO, optarg={"bound_push": 1e-2})
        badly_scaled_vars = list(iscale.badly_scaled_var_generator(unit))
        if len(badly_scaled_vars) > 0:
            automate_rescale_variables(unit)

    seq.run(m, function)

    results = solve(m, tee=True)

    return m, results


def solve(blk, solver=None, checkpoint=None, tee=False, fail_flag=True):
    if solver is None:
        solver = get_solver()
    results = solver.solve(blk, tee=tee)
    check_solve(results, checkpoint=checkpoint, logger=_log, fail_flag=fail_flag)
    return results


if __name__ == "__main__":
    # This method builds and runs a steady state activated sludge
    # flowsheet.
    m, results = build_flowsheet()

    stream_table = create_stream_table_dataframe(
        {
            "Feed": m.fs.feed.outlet,
            "Mix": m.fs.R1.inlet,
            "R1": m.fs.R1.outlet,
            "R2": m.fs.R2.outlet,
            "R3": m.fs.R3.outlet,
            "R4": m.fs.R4.outlet,
            "R5": m.fs.R5.outlet,
            "thickener": m.fs.thickener.underflow,
            "AD liquid outlet": m.fs.AD.liquid_outlet,
            "AD vapor outlet": m.fs.AD.vapor_outlet,
            "dewater": m.fs.dewater.overflow,
            "ElectroNP treated": m.fs.electroNP.treated,
            "ElectroNP byproduct": m.fs.electroNP.byproduct,
        },
        time_point=0,
    )
    print(stream_table_dataframe_to_string(stream_table))
