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

from watertap.core.util.initialization import check_solve, assert_degrees_of_freedom
from watertap.costing import WaterTAPCosting
from watertap.costing.unit_models.clarifier import (
    cost_circular_clarifier,
    cost_primary_clarifier,
)


# Set up logger
_log = idaeslog.getLogger(__name__)


def main(bio_P=True):
    m = build(bio_P=bio_P)
    set_operating_conditions(m)
    for mx in m.fs.mixers:
        mx.pressure_equality_constraints[0.0, 2].deactivate()
    m.fs.MX3.pressure_equality_constraints[0.0, 2].deactivate()
    m.fs.MX3.pressure_equality_constraints[0.0, 3].deactivate()
    print(f"DOF before initialization: {degrees_of_freedom(m)}")

    initialize_system(m, bio_P=bio_P)
    for mx in m.fs.mixers:
        mx.pressure_equality_constraints[0.0, 2].deactivate()
    m.fs.MX3.pressure_equality_constraints[0.0, 2].deactivate()
    m.fs.MX3.pressure_equality_constraints[0.0, 3].deactivate()
    print(f"DOF after initialization: {degrees_of_freedom(m)}")

    results = solve(m)

    # Switch to fixed KLa in R5, R6, and R7 (S_O concentration is controlled in R5)
    m.fs.R5.KLa.fix(240)
    m.fs.R6.KLa.fix(240)
    m.fs.R7.KLa.fix(84)
    m.fs.R5.outlet.conc_mass_comp[:, "S_O2"].unfix()
    m.fs.R6.outlet.conc_mass_comp[:, "S_O2"].unfix()
    m.fs.R7.outlet.conc_mass_comp[:, "S_O2"].unfix()
    # Resolve with controls in place
    results = solve(m)

    pyo.assert_optimal_termination(results)
    check_solve(
        results,
        checkpoint="re-solve with controls in place",
        logger=_log,
        fail_flag=True,
    )

    add_costing(m)
    m.fs.costing.initialize()

    assert_degrees_of_freedom(m, 0)

    results = solve(m)
    pyo.assert_optimal_termination(results)

    display_costing(m)
    display_performance_metrics(m)

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

    # ====================================================================
    # Primary Clarifier
    m.fs.CL = Clarifier(
        property_package=m.fs.props_ASM2D,
        outlet_list=["underflow", "effluent"],
        split_basis=SplittingType.componentFlow,
    )

    # ======================================================================
    # Activated Sludge Process
    # Mixer for feed water and recycled sludge
    m.fs.MX1 = Mixer(
        property_package=m.fs.props_ASM2D,
        inlet_list=["feed_water", "recycle"],
        momentum_mixing_type=MomentumMixingType.equality,
    )
    # First reactor (anaerobic) - standard CSTR
    m.fs.R1 = CSTR(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    # Second reactor (anaerobic) - standard CSTR
    m.fs.R2 = CSTR(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    # Third reactor (anoxic) - standard CSTR
    m.fs.R3 = CSTR(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    # Fourth reactor (anoxic) - standard CSTR
    m.fs.R4 = CSTR(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    # Fifth reactor (aerobic) - CSTR with injection
    m.fs.R5 = CSTR_Injection(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    # Sixth reactor (aerobic) - CSTR with injection
    m.fs.R6 = CSTR_Injection(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    # Seventh reactor (aerobic) - CSTR with injection
    m.fs.R7 = CSTR_Injection(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    m.fs.SP1 = Separator(
        property_package=m.fs.props_ASM2D, outlet_list=["underflow", "overflow"]
    )
    # Secondary Clarifier
    m.fs.CL2 = Clarifier(
        property_package=m.fs.props_ASM2D,
        outlet_list=["underflow", "effluent"],
        split_basis=SplittingType.componentFlow,
    )
    # Mixing sludge recycle and R5 underflow
    m.fs.MX2 = Mixer(
        property_package=m.fs.props_ASM2D,
        inlet_list=["reactor", "clarifier"],
        momentum_mixing_type=MomentumMixingType.equality,
    )
    # Sludge separator
    m.fs.SP2 = Separator(
        property_package=m.fs.props_ASM2D, outlet_list=["waste", "recycle"]
    )
    # Recycle pressure changer - use a simple isothermal unit for now
    m.fs.P1 = PressureChanger(property_package=m.fs.props_ASM2D)

    # ======================================================================
    # Thickener
    m.fs.thickener = Thickener(
        property_package=m.fs.props_ASM2D,
        activated_sludge_model=thickener_type.modified_ASM2D,
    )
    # Mixing feed and recycle streams from thickener and dewatering unit
    m.fs.MX3 = Mixer(
        property_package=m.fs.props_ASM2D,
        inlet_list=["feed_water", "recycle1", "recycle2"],
        momentum_mixing_type=MomentumMixingType.equality,
    )
    # Mixing sludge from thickener and primary clarifier
    m.fs.MX4 = Mixer(
        property_package=m.fs.props_ASM2D,
        inlet_list=["thickener", "clarifier"],
        momentum_mixing_type=MomentumMixingType.equality,
    )

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
        reaction_package=m.fs.rxn_props_ADM1,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )

    # ======================================================================
    # Dewatering Unit
    m.fs.dewater = DewateringUnit(
        property_package=m.fs.props_ASM2D,
        activated_sludge_model=dewater_type.modified_ASM2D,
    )

    # ======================================================================
    # Product Blocks
    m.fs.Treated = Product(property_package=m.fs.props_ASM2D)
    m.fs.Sludge = Product(property_package=m.fs.props_ASM2D)
    # Mixers
    m.fs.mixers = (m.fs.MX1, m.fs.MX2, m.fs.MX4)

    # ======================================================================
    # Link units related to ASM section
    m.fs.stream2 = Arc(source=m.fs.MX1.outlet, destination=m.fs.R1.inlet)
    m.fs.stream3 = Arc(source=m.fs.R1.outlet, destination=m.fs.R2.inlet)
    m.fs.stream4 = Arc(source=m.fs.R2.outlet, destination=m.fs.MX2.reactor)
    m.fs.stream5 = Arc(source=m.fs.MX2.outlet, destination=m.fs.R3.inlet)
    m.fs.stream6 = Arc(source=m.fs.R3.outlet, destination=m.fs.R4.inlet)
    m.fs.stream7 = Arc(source=m.fs.R4.outlet, destination=m.fs.R5.inlet)
    m.fs.stream8 = Arc(source=m.fs.R5.outlet, destination=m.fs.R6.inlet)
    m.fs.stream9 = Arc(source=m.fs.R6.outlet, destination=m.fs.R7.inlet)
    m.fs.stream10 = Arc(source=m.fs.R7.outlet, destination=m.fs.SP1.inlet)
    m.fs.stream11 = Arc(source=m.fs.SP1.overflow, destination=m.fs.CL2.inlet)
    m.fs.stream12 = Arc(source=m.fs.SP1.underflow, destination=m.fs.MX2.clarifier)
    m.fs.stream13 = Arc(source=m.fs.CL2.effluent, destination=m.fs.Treated.inlet)
    m.fs.stream14 = Arc(source=m.fs.CL2.underflow, destination=m.fs.SP2.inlet)
    m.fs.stream15 = Arc(source=m.fs.SP2.recycle, destination=m.fs.P1.inlet)
    m.fs.stream16 = Arc(source=m.fs.P1.outlet, destination=m.fs.MX1.recycle)

    # Link units related to AD section
    m.fs.stream_AD_translator = Arc(
        source=m.fs.AD.liquid_outlet, destination=m.fs.translator_adm1_asm2d.inlet
    )
    m.fs.stream_SP_thickener = Arc(
        source=m.fs.SP2.waste, destination=m.fs.thickener.inlet
    )
    m.fs.stream3adm = Arc(
        source=m.fs.thickener.underflow, destination=m.fs.MX4.thickener
    )
    m.fs.stream7adm = Arc(source=m.fs.thickener.overflow, destination=m.fs.MX3.recycle2)
    m.fs.stream9adm = Arc(source=m.fs.CL.underflow, destination=m.fs.MX4.clarifier)
    m.fs.stream_translator_dewater = Arc(
        source=m.fs.translator_adm1_asm2d.outlet, destination=m.fs.dewater.inlet
    )
    m.fs.stream1a = Arc(source=m.fs.FeedWater.outlet, destination=m.fs.MX3.feed_water)
    m.fs.stream1b = Arc(source=m.fs.MX3.outlet, destination=m.fs.CL.inlet)
    m.fs.stream1c = Arc(source=m.fs.CL.effluent, destination=m.fs.MX1.feed_water)
    m.fs.stream_dewater_sludge = Arc(
        source=m.fs.dewater.underflow, destination=m.fs.Sludge.inlet
    )
    m.fs.stream_dewater_mixer = Arc(
        source=m.fs.dewater.overflow, destination=m.fs.MX3.recycle1
    )
    m.fs.stream10adm = Arc(
        source=m.fs.MX4.outlet, destination=m.fs.translator_asm2d_adm1.inlet
    )
    m.fs.stream_translator_AD = Arc(
        source=m.fs.translator_asm2d_adm1.outlet, destination=m.fs.AD.inlet
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

    @m.fs.R5.Constraint(m.fs.time, doc="Mass transfer constraint for R5")
    def mass_transfer_R5(self, t):
        return pyo.units.convert(
            m.fs.R5.injection[t, "Liq", "S_O2"], to_units=pyo.units.kg / pyo.units.hour
        ) == (
            m.fs.R5.KLa
            * m.fs.R5.volume[t]
            * (m.fs.S_O_eq - m.fs.R5.outlet.conc_mass_comp[t, "S_O2"])
        )

    @m.fs.R6.Constraint(m.fs.time, doc="Mass transfer constraint for R6")
    def mass_transfer_R6(self, t):
        return pyo.units.convert(
            m.fs.R6.injection[t, "Liq", "S_O2"], to_units=pyo.units.kg / pyo.units.hour
        ) == (
            m.fs.R6.KLa
            * m.fs.R6.volume[t]
            * (m.fs.S_O_eq - m.fs.R6.outlet.conc_mass_comp[t, "S_O2"])
        )

    @m.fs.R7.Constraint(m.fs.time, doc="Mass transfer constraint for R7")
    def mass_transfer_R7(self, t):
        return pyo.units.convert(
            m.fs.R7.injection[t, "Liq", "S_O2"], to_units=pyo.units.kg / pyo.units.hour
        ) == (
            m.fs.R7.KLa
            * m.fs.R7.volume[t]
            * (m.fs.S_O_eq - m.fs.R7.outlet.conc_mass_comp[t, "S_O2"])
        )

    return m


def set_operating_conditions(m):
    # Feed Water Conditions
    print(f"DOF before feed: {degrees_of_freedom(m)}")
    m.fs.FeedWater.flow_vol.fix(20935.15 * pyo.units.m**3 / pyo.units.day)
    m.fs.FeedWater.temperature.fix(308.15 * pyo.units.K)
    m.fs.FeedWater.pressure.fix(1 * pyo.units.atm)
    m.fs.FeedWater.conc_mass_comp[0, "S_O2"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_F"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_A"].fix(70 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_NH4"].fix(26.6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_NO3"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_PO4"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_I"].fix(57.45 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_N2"].fix(25.19 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_I"].fix(84 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_S"].fix(94.1 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_H"].fix(370 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_PAO"].fix(
        51.5262 * pyo.units.g / pyo.units.m**3
    )
    m.fs.FeedWater.conc_mass_comp[0, "X_PP"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_PHA"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_AUT"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_IC"].fix(5.652 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_K"].fix(374.6925 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_Mg"].fix(20 * pyo.units.g / pyo.units.m**3)

    # Primary Clarifier
    # TODO: Update primary clarifier once more detailed model available
    m.fs.CL.split_fraction[0, "effluent", "H2O"].fix(0.993)
    m.fs.CL.split_fraction[0, "effluent", "S_A"].fix(0.993)
    m.fs.CL.split_fraction[0, "effluent", "S_F"].fix(0.993)
    m.fs.CL.split_fraction[0, "effluent", "S_I"].fix(0.993)
    m.fs.CL.split_fraction[0, "effluent", "S_N2"].fix(0.993)
    m.fs.CL.split_fraction[0, "effluent", "S_NH4"].fix(0.993)
    m.fs.CL.split_fraction[0, "effluent", "S_NO3"].fix(0.993)
    m.fs.CL.split_fraction[0, "effluent", "S_O2"].fix(0.993)
    m.fs.CL.split_fraction[0, "effluent", "S_PO4"].fix(0.993)
    m.fs.CL.split_fraction[0, "effluent", "S_IC"].fix(0.993)
    m.fs.CL.split_fraction[0, "effluent", "S_K"].fix(0.993)
    m.fs.CL.split_fraction[0, "effluent", "S_Mg"].fix(0.993)
    m.fs.CL.split_fraction[0, "effluent", "X_AUT"].fix(0.5192)
    m.fs.CL.split_fraction[0, "effluent", "X_H"].fix(0.5192)
    m.fs.CL.split_fraction[0, "effluent", "X_I"].fix(0.5192)
    m.fs.CL.split_fraction[0, "effluent", "X_PAO"].fix(0.5192)
    m.fs.CL.split_fraction[0, "effluent", "X_PHA"].fix(0.5192)
    m.fs.CL.split_fraction[0, "effluent", "X_PP"].fix(0.5192)
    m.fs.CL.split_fraction[0, "effluent", "X_S"].fix(0.5192)

    # Reactor sizing
    m.fs.R1.volume.fix(1000 * pyo.units.m**3)
    m.fs.R2.volume.fix(1000 * pyo.units.m**3)
    m.fs.R3.volume.fix(1500 * pyo.units.m**3)
    m.fs.R4.volume.fix(1500 * pyo.units.m**3)
    m.fs.R5.volume.fix(3000 * pyo.units.m**3)
    m.fs.R6.volume.fix(3000 * pyo.units.m**3)
    m.fs.R7.volume.fix(3000 * pyo.units.m**3)

    # Injection rates to Reactions 5, 6 and 7
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

    # Set fraction of outflow from reactor 7 that goes to recycle
    m.fs.SP1.split_fraction[:, "underflow"].fix(0.60)

    # Secondary Clarifier
    # TODO: Update once more detailed model available
    m.fs.CL2.split_fraction[0, "effluent", "H2O"].fix(0.48956)
    m.fs.CL2.split_fraction[0, "effluent", "S_A"].fix(0.48956)
    m.fs.CL2.split_fraction[0, "effluent", "S_F"].fix(0.48956)
    m.fs.CL2.split_fraction[0, "effluent", "S_I"].fix(0.48956)
    m.fs.CL2.split_fraction[0, "effluent", "S_N2"].fix(0.48956)
    m.fs.CL2.split_fraction[0, "effluent", "S_NH4"].fix(0.48956)
    m.fs.CL2.split_fraction[0, "effluent", "S_NO3"].fix(0.48956)
    m.fs.CL2.split_fraction[0, "effluent", "S_O2"].fix(0.48956)
    m.fs.CL2.split_fraction[0, "effluent", "S_PO4"].fix(0.48956)
    m.fs.CL2.split_fraction[0, "effluent", "S_IC"].fix(0.48956)
    m.fs.CL2.split_fraction[0, "effluent", "S_K"].fix(0.48956)
    m.fs.CL2.split_fraction[0, "effluent", "S_Mg"].fix(0.48956)
    m.fs.CL2.split_fraction[0, "effluent", "X_AUT"].fix(0.00187)
    m.fs.CL2.split_fraction[0, "effluent", "X_H"].fix(0.00187)
    m.fs.CL2.split_fraction[0, "effluent", "X_I"].fix(0.00187)
    m.fs.CL2.split_fraction[0, "effluent", "X_PAO"].fix(0.00187)
    m.fs.CL2.split_fraction[0, "effluent", "X_PHA"].fix(0.00187)
    m.fs.CL2.split_fraction[0, "effluent", "X_PP"].fix(0.00187)
    m.fs.CL2.split_fraction[0, "effluent", "X_S"].fix(0.00187)

    m.fs.CL2.surface_area.fix(1500 * pyo.units.m**2)

    # Sludge purge separator
    m.fs.SP2.split_fraction[:, "recycle"].fix(0.985)

    # Outlet pressure from recycle pump
    m.fs.P1.outlet.pressure.fix(101325)

    # AD
    m.fs.AD.volume_liquid.fix(3400)
    m.fs.AD.volume_vapor.fix(300)
    m.fs.AD.liquid_outlet.temperature.fix(308.15)

    # Dewatering Unit - fix either HRT or volume.
    m.fs.dewater.hydraulic_retention_time.fix(1800 * pyo.units.s)

    # Thickener unit
    m.fs.thickener.hydraulic_retention_time.fix(86400 * pyo.units.s)
    m.fs.thickener.diameter.fix(10 * pyo.units.m)

    def scale_variables(m):
        for var in m.fs.component_data_objects(pyo.Var, descend_into=True):
            if "flow_vol" in var.name:
                iscale.set_scaling_factor(var, 1e1)
            if "temperature" in var.name:
                iscale.set_scaling_factor(var, 1e-2)
            if "pressure" in var.name:
                iscale.set_scaling_factor(var, 1e-4)
            if "conc_mass_comp" in var.name:
                iscale.set_scaling_factor(var, 1e2)

    for unit in ("R1", "R2", "R3", "R4", "R5", "R6", "R7"):
        block = getattr(m.fs, unit)
        iscale.set_scaling_factor(
            block.control_volume.reactions[0.0].rate_expression, 1e3
        )
        iscale.set_scaling_factor(block.cstr_performance_eqn, 1e3)
        iscale.set_scaling_factor(
            block.control_volume.rate_reaction_stoichiometry_constraint, 1e3
        )
        iscale.set_scaling_factor(block.control_volume.material_balances, 1e3)

    # Apply scaling
    scale_variables(m)
    iscale.calculate_scaling_factors(m)


def initialize_system(m, bio_P=False):
    # Initialize flowsheet
    # Apply sequential decomposition - 1 iteration should suffice
    seq = SequentialDecomposition()
    seq.options.tear_method = "Direct"
    seq.options.iterLim = 5
    seq.options.tear_set = [m.fs.stream5, m.fs.stream10adm]

    G = seq.create_graph(m)
    # Uncomment this code to see tear set and initialization order
    order = seq.calculation_order(G)
    print("Initialization Order")
    for o in order:
        print(o[0].name)

    if bio_P:
        # Initial guesses for flow into first reactor
        tear_guesses = {
            "flow_vol": {0: 1.2368},
            "conc_mass_comp": {
                (0, "S_A"): 0.0009,
                (0, "S_F"): 0.00042,
                (0, "S_I"): 0.05745,
                (0, "S_N2"): 0.0534,
                (0, "S_NH4"): 0.0089,
                (0, "S_NO3"): 0.0046,
                (0, "S_O2"): 0.0046,
                (0, "S_PO4"): 0.0118,
                (0, "S_K"): 0.373,
                (0, "S_Mg"): 0.023,
                (0, "S_IC"): 0.1366,
                (0, "X_AUT"): 0.135,
                (0, "X_H"): 3.647,
                (0, "X_I"): 3.314,
                (0, "X_PAO"): 2.984,
                (0, "X_PHA"): 0.083,
                (0, "X_PP"): 0.9907,
                (0, "X_S"): 0.05823,
            },
            "temperature": {0: 308.15},
            "pressure": {0: 101325},
        }

        tear_guesses2 = {
            "flow_vol": {0: 0.003},
            "conc_mass_comp": {
                (0, "S_A"): 0.10,
                (0, "S_F"): 0.16,
                (0, "S_I"): 0.05745,
                (0, "S_N2"): 0.039,
                (0, "S_NH4"): 0.035,
                (0, "S_NO3"): 0.0032,
                (0, "S_O2"): 0.00314,
                (0, "S_PO4"): 0.0238,
                (0, "S_K"): 0.379,
                (0, "S_Mg"): 0.026,
                (0, "S_IC"): 0.078,
                (0, "X_AUT"): 0.342,
                (0, "X_H"): 23.95,
                (0, "X_I"): 11.9,
                (0, "X_PAO"): 9.62,
                (0, "X_PHA"): 0.0033,
                (0, "X_PP"): 2.51,
                (0, "X_S"): 3.92,
            },
            "temperature": {0: 308.15},
            "pressure": {0: 101325},
        }

    else:
        # Initial guesses for flow into first reactor
        tear_guesses = {
            "flow_vol": {0: 1.235},
            "conc_mass_comp": {
                (0, "S_A"): 0.0007,
                (0, "S_F"): 0.0004,
                (0, "S_I"): 0.0575,
                (0, "S_N2"): 0.05,
                (0, "S_NH4"): 0.007,
                (0, "S_NO3"): 0.0035,
                (0, "S_O2"): 0.00192,
                (0, "S_PO4"): 0.02,
                (0, "S_K"): 0.37,
                (0, "S_Mg"): 0.02,
                (0, "S_IC"): 0.11,
                (0, "X_AUT"): 0.12,
                (0, "X_H"): 3.3,
                (0, "X_I"): 3.0,
                (0, "X_PAO"): 2.3,
                (0, "X_PHA"): 0.06,
                (0, "X_PP"): 0.75,
                (0, "X_S"): 0.050,
            },
            "temperature": {0: 308.15},
            "pressure": {0: 101325},
        }

        tear_guesses2 = {
            "flow_vol": {0: 0.0027},
            "conc_mass_comp": {
                (0, "S_A"): 0.044,
                (0, "S_F"): 0.15,
                (0, "S_I"): 0.0575,
                (0, "S_N2"): 0.035,
                (0, "S_NH4"): 0.03,
                (0, "S_NO3"): 0.002,
                (0, "S_O2"): 0.0012,
                (0, "S_PO4"): 0.02,
                (0, "S_K"): 0.38,
                (0, "S_Mg"): 0.023,
                (0, "S_IC"): 0.063,
                (0, "X_AUT"): 0.31,
                (0, "X_H"): 24.8,
                (0, "X_I"): 11.8,
                (0, "X_PAO"): 8.5,
                (0, "X_PHA"): 0.086,
                (0, "X_PP"): 2.1,
                (0, "X_S"): 4.2,
            },
            "temperature": {0: 308.15},
            "pressure": {0: 101325},
        }

    # Pass the tear_guess to the SD tool
    seq.set_guesses_for(m.fs.R3.inlet, tear_guesses)
    seq.set_guesses_for(m.fs.translator_asm2d_adm1.inlet, tear_guesses2)

    def function(unit):
        unit.initialize(outlvl=idaeslog.INFO)

    seq.run(m, function)


def solve(m, solver=None):
    if solver is None:
        solver = get_solver()
    results = solver.solve(m, tee=True)
    check_solve(results, checkpoint="closing recycle", logger=_log, fail_flag=True)
    pyo.assert_optimal_termination(results)
    return results


def add_costing(m):
    m.fs.costing = WaterTAPCosting()
    m.fs.costing.base_currency = pyo.units.USD_2020

    # Costing Blocks
    m.fs.R1.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.R2.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.R3.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.R4.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.R5.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.R6.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.R7.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.CL.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=cost_primary_clarifier,
    )

    m.fs.CL2.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=cost_circular_clarifier,
    )

    m.fs.AD.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.dewater.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.thickener.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    # TODO: Leaving out mixer costs; consider including later

    # process costing and add system level metrics
    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(m.fs.Treated.properties[0].flow_vol)
    m.fs.costing.add_LCOW(m.fs.FeedWater.properties[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(m.fs.FeedWater.properties[0].flow_vol)

    m.fs.objective = pyo.Objective(expr=m.fs.costing.LCOW)
    iscale.set_scaling_factor(m.fs.costing.LCOW, 1e3)
    iscale.set_scaling_factor(m.fs.costing.total_capital_cost, 1e-7)
    iscale.set_scaling_factor(m.fs.costing.total_capital_cost, 1e-5)

    iscale.calculate_scaling_factors(m.fs)


def display_costing(m):
    print("Levelized cost of water: %.2f $/m3" % pyo.value(m.fs.costing.LCOW))

    print(
        "Total operating cost: %.2f $/yr" % pyo.value(m.fs.costing.total_operating_cost)
    )
    print("Total capital cost: %.2f $" % pyo.value(m.fs.costing.total_capital_cost))

    print(
        "Total annualized cost: %.2f $/yr"
        % pyo.value(m.fs.costing.total_annualized_cost)
    )

    print(
        "capital cost R1",
        pyo.value(m.fs.R1.costing.capital_cost),
        pyo.units.get_units(m.fs.R1.costing.capital_cost),
    )
    print(
        "capital cost R2",
        pyo.value(m.fs.R2.costing.capital_cost),
        pyo.units.get_units(m.fs.R2.costing.capital_cost),
    )
    print(
        "capital cost R3",
        pyo.value(m.fs.R3.costing.capital_cost),
        pyo.units.get_units(m.fs.R3.costing.capital_cost),
    )
    print(
        "capital cost R4",
        pyo.value(m.fs.R4.costing.capital_cost),
        pyo.units.get_units(m.fs.R4.costing.capital_cost),
    )
    print(
        "capital cost R5",
        pyo.value(m.fs.R5.costing.capital_cost),
        pyo.units.get_units(m.fs.R5.costing.capital_cost),
    )
    print(
        "capital cost R6",
        pyo.value(m.fs.R6.costing.capital_cost),
        pyo.units.get_units(m.fs.R6.costing.capital_cost),
    )
    print(
        "capital cost R7",
        pyo.value(m.fs.R7.costing.capital_cost),
        pyo.units.get_units(m.fs.R7.costing.capital_cost),
    )
    print(
        "capital cost primary clarifier",
        pyo.value(m.fs.CL.costing.capital_cost),
        pyo.units.get_units(m.fs.CL.costing.capital_cost),
    )
    print(
        "capital cost secondary clarifier",
        pyo.value(m.fs.CL2.costing.capital_cost),
        pyo.units.get_units(m.fs.CL2.costing.capital_cost),
    )
    print(
        "capital cost AD",
        pyo.value(m.fs.AD.costing.capital_cost),
        pyo.units.get_units(m.fs.AD.costing.capital_cost),
    )
    print(
        "capital cost dewatering Unit",
        pyo.value(m.fs.dewater.costing.capital_cost),
        pyo.units.get_units(m.fs.dewater.costing.capital_cost),
    )
    print(
        "capital cost thickener unit",
        pyo.value(m.fs.thickener.costing.capital_cost),
        pyo.units.get_units(m.fs.thickener.costing.capital_cost),
    )


def display_performance_metrics(m):
    print(
        "Specific energy consumption with respect to influent flowrate: %.1f kWh/m3"
        % pyo.value(m.fs.costing.specific_energy_consumption)
    )

    print(
        "electricity consumption R5",
        pyo.value(m.fs.R5.electricity_consumption[0]),
        pyo.units.get_units(m.fs.R5.electricity_consumption[0]),
    )
    print(
        "electricity consumption R6",
        pyo.value(m.fs.R6.electricity_consumption[0]),
        pyo.units.get_units(m.fs.R6.electricity_consumption[0]),
    )
    print(
        "electricity consumption R7",
        pyo.value(m.fs.R7.electricity_consumption[0]),
        pyo.units.get_units(m.fs.R7.electricity_consumption[0]),
    )
    print(
        "electricity consumption primary clarifier",
        pyo.value(m.fs.CL.electricity_consumption[0]),
        pyo.units.get_units(m.fs.CL.electricity_consumption[0]),
    )
    print(
        "electricity consumption secondary clarifier",
        pyo.value(m.fs.CL2.electricity_consumption[0]),
        pyo.units.get_units(m.fs.CL2.electricity_consumption[0]),
    )
    print(
        "electricity consumption AD",
        pyo.value(m.fs.AD.electricity_consumption[0]),
        pyo.units.get_units(m.fs.AD.electricity_consumption[0]),
    )
    print(
        "electricity consumption dewatering Unit",
        pyo.value(m.fs.dewater.electricity_consumption[0]),
        pyo.units.get_units(m.fs.dewater.electricity_consumption[0]),
    )
    print(
        "electricity consumption thickening Unit",
        pyo.value(m.fs.thickener.electricity_consumption[0]),
        pyo.units.get_units(m.fs.thickener.electricity_consumption[0]),
    )
    print(
        "Influent flow",
        pyo.value(m.fs.FeedWater.flow_vol[0]),
        pyo.units.get_units(m.fs.FeedWater.flow_vol[0]),
    )
    print(
        "flow into R3",
        pyo.value(m.fs.R3.control_volume.properties_in[0].flow_vol),
        pyo.units.get_units(m.fs.R3.control_volume.properties_in[0].flow_vol),
    )
    print(
        "flow into RADM",
        pyo.value(m.fs.AD.liquid_phase.properties_in[0].flow_vol),
        pyo.units.get_units(m.fs.AD.liquid_phase.properties_in[0].flow_vol),
    )


if __name__ == "__main__":
    # This method builds and runs a steady state activated sludge flowsheet.
    m, results = main()

    stream_table = create_stream_table_dataframe(
        {
            "Feed": m.fs.FeedWater.outlet,
            "R1": m.fs.R1.outlet,
            "R2": m.fs.R2.outlet,
            "R3": m.fs.R3.outlet,
            "R4": m.fs.R4.outlet,
            "R5": m.fs.R5.outlet,
            "R6": m.fs.R6.outlet,
            "R7": m.fs.R7.outlet,
            "thickener outlet": m.fs.thickener.underflow,
            "ADM-ASM translator outlet": m.fs.translator_adm1_asm2d.outlet,
            "dewater outlet": m.fs.dewater.overflow,
            "Treated water": m.fs.Treated.inlet,
            "Sludge": m.fs.Sludge.inlet,
        },
        time_point=0,
    )
    print(stream_table_dataframe_to_string(stream_table))
