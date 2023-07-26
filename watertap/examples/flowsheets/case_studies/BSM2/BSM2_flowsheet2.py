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
__author__ = "Alejandro Garciadiego"

import pyomo.environ as pyo
from pyomo.environ import (
    units,
)

from pyomo.network import Arc, SequentialDecomposition
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.unit_models.anaerobic_digestor import AD
from watertap.unit_models.thickener import Thickener
from watertap.unit_models.dewatering import DewateringUnit

from watertap.unit_models.translators.translator_asm1_adm1 import Translator_ASM1_ADM1
from watertap.unit_models.translators.translator_adm1_asm1 import Translator_ADM1_ASM1
from idaes.models.unit_models import Separator, Mixer

from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    large_residuals_set,
    activated_inequalities_set,
)
import idaes.logger as idaeslog
from idaes.core.util.initialization import propagate_state
from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale
from idaes.core.util.tables import (
    create_stream_table_dataframe,
    stream_table_dataframe_to_string,
)

from watertap.property_models.anaerobic_digestion.adm1_properties import (
    ADM1ParameterBlock,
)
from watertap.property_models.activated_sludge.asm1_properties import (
    ASM1ParameterBlock,
)
from watertap.property_models.anaerobic_digestion.adm1_reactions import (
    ADM1ReactionParameterBlock,
)
from idaes.models.unit_models.separator import SplittingType
from watertap.property_models.anaerobic_digestion.adm1_properties_vapor import (
    ADM1_vaporParameterBlock,
)

from idaes.core.util.tables import (
    create_stream_table_dataframe,
    stream_table_dataframe_to_string,
)

from watertap.core.util.initialization import check_solve

from idaes.core import FlowsheetBlock
from idaes.models.unit_models import (
    CSTR,
    Feed,
    Mixer,
    Separator,
    PressureChanger,
    Product,
)

from watertap.unit_models.cstr_injection import CSTR_Injection
from watertap.property_models.activated_sludge.asm1_properties import ASM1ParameterBlock
from watertap.property_models.activated_sludge.asm1_reactions import (
    ASM1ReactionParameterBlock,
)

from watertap.core.util.initialization import check_solve

from pyomo.environ import SolverFactory, Objective
from idaes.core.util.model_diagnostics import DegeneracyHunter
from idaes.core.util.scaling import report_scaling_issues, badly_scaled_var_generator,list_badly_scaled_variables


def build_flowsheet():
    m = pyo.ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.props_ASM1 = ASM1ParameterBlock()
    m.fs.props_ADM1 = ADM1ParameterBlock()
    m.fs.props_vap = ADM1_vaporParameterBlock()
    m.fs.ADM1_rxn_props = ADM1ReactionParameterBlock(property_package=m.fs.props_ADM1)
    m.fs.ASM1_rxn_props = ASM1ReactionParameterBlock(property_package=m.fs.props_ASM1)
    # Feed water stream
    m.fs.FeedWater = Feed(property_package=m.fs.props_ASM1)
    # Mixer for feed water and recycled sludge
    m.fs.MX1 = Mixer(
        property_package=m.fs.props_ASM1, inlet_list=["feed_water", "recycle"]
    )
    # First reactor (anoxic) - standard CSTR
    m.fs.R1 = CSTR(
        property_package=m.fs.props_ASM1, reaction_package=m.fs.ASM1_rxn_props
    )
    # Second reactor (anoxic) - standard CSTR
    m.fs.R2 = CSTR(
        property_package=m.fs.props_ASM1, reaction_package=m.fs.ASM1_rxn_props
    )
    # Third reactor (aerobic) - CSTR with injection
    m.fs.R3 = CSTR_Injection(
        property_package=m.fs.props_ASM1, reaction_package=m.fs.ASM1_rxn_props
    )
    # Fourth reactor (aerobic) - CSTR with injection
    m.fs.R4 = CSTR_Injection(
        property_package=m.fs.props_ASM1, reaction_package=m.fs.ASM1_rxn_props
    )
    # Fifth reactor (aerobic) - CSTR with injection
    m.fs.R5 = CSTR_Injection(
        property_package=m.fs.props_ASM1, reaction_package=m.fs.ASM1_rxn_props
    )
    m.fs.SP5 = Separator(
        property_package=m.fs.props_ASM1, outlet_list=["underflow", "overflow"]
    )
    # Clarifier
    # TODO: Replace with more detailed model when available
    m.fs.CL1 = Separator(
        property_package=m.fs.props_ASM1,
        outlet_list=["underflow", "effluent"],
        split_basis=SplittingType.componentFlow,
    )
    # Sludge purge splitter
    m.fs.SP6 = Separator(
        property_package=m.fs.props_ASM1,
        outlet_list=["recycle", "waste"],
        split_basis=SplittingType.totalFlow,
    )
    # Mixing sludge recycle and R5 underflow
    m.fs.MX6 = Mixer(
        property_package=m.fs.props_ASM1, inlet_list=["clarifier", "reactor"]
    )
    # Product Blocks
    # m.fs.Treated = Product(property_package=m.fs.props_ASM1)
    # Recycle pressure changer - use a simple isothermal unit for now
    m.fs.P1 = PressureChanger(property_package=m.fs.props_ASM1)

    # Link units
    m.fs.stream2 = Arc(source=m.fs.MX1.outlet, destination=m.fs.R1.inlet)
    m.fs.stream3 = Arc(source=m.fs.R1.outlet, destination=m.fs.R2.inlet)
    m.fs.stream4 = Arc(source=m.fs.R2.outlet, destination=m.fs.R3.inlet)
    m.fs.stream5 = Arc(source=m.fs.R3.outlet, destination=m.fs.R4.inlet)
    m.fs.stream6 = Arc(source=m.fs.R4.outlet, destination=m.fs.R5.inlet)
    m.fs.stream7 = Arc(source=m.fs.R5.outlet, destination=m.fs.SP5.inlet)
    m.fs.stream8 = Arc(source=m.fs.SP5.overflow, destination=m.fs.CL1.inlet)
    m.fs.stream9 = Arc(source=m.fs.SP5.underflow, destination=m.fs.MX6.reactor)
    # m.fs.stream10 = Arc(source=m.fs.CL1.effluent, destination=m.fs.Treated.inlet)
    m.fs.stream11 = Arc(source=m.fs.CL1.underflow, destination=m.fs.SP6.inlet)
    m.fs.stream13 = Arc(source=m.fs.SP6.recycle, destination=m.fs.MX6.clarifier)
    m.fs.stream14 = Arc(source=m.fs.MX6.outlet, destination=m.fs.P1.inlet)
    m.fs.stream15 = Arc(source=m.fs.P1.outlet, destination=m.fs.MX1.recycle)
    pyo.TransformationFactory("network.expand_arcs").apply_to(m)

    # Oxygen concentration in reactors 3 and 4 is governed by mass transfer
    # Add additional parameter and constraints
    m.fs.R3.KLa = pyo.Var(
        initialize=7.6,
        units=pyo.units.hour**-1,
        doc="Lumped mass transfer coefficient for oxygen",
    )
    m.fs.R4.KLa = pyo.Var(
        initialize=5.7,
        units=pyo.units.hour**-1,
        doc="Lumped mass transfer coefficient for oxygen",
    )
    m.fs.S_O_eq = pyo.Param(
        default=8e-3,
        units=pyo.units.kg / pyo.units.m**3,
        mutable=True,
        doc="Dissolved oxygen concentration at equilibrium",
    )

    @m.fs.R3.Constraint(m.fs.time, doc="Mass transfer constraint for R3")
    def mass_transfer_R3(self, t):
        return pyo.units.convert(
            m.fs.R3.injection[t, "Liq", "S_O"], to_units=pyo.units.kg / pyo.units.hour
        ) == (
            m.fs.R3.KLa
            * m.fs.R3.volume[t]
            * (m.fs.S_O_eq - m.fs.R3.outlet.conc_mass_comp[t, "S_O"])
        )

    @m.fs.R4.Constraint(m.fs.time, doc="Mass transfer constraint for R4")
    def mass_transfer_R4(self, t):
        return pyo.units.convert(
            m.fs.R4.injection[t, "Liq", "S_O"], to_units=pyo.units.kg / pyo.units.hour
        ) == (
            m.fs.R4.KLa
            * m.fs.R4.volume[t]
            * (m.fs.S_O_eq - m.fs.R4.outlet.conc_mass_comp[t, "S_O"])
        )

    # Feed Water Conditions
    m.fs.FeedWater.flow_vol.fix(20648 * pyo.units.m**3 / pyo.units.day)
    m.fs.FeedWater.temperature.fix(308.15 * pyo.units.K)
    m.fs.FeedWater.pressure.fix(1 * pyo.units.atm)
    m.fs.FeedWater.conc_mass_comp[0, "S_I"].fix(27 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_S"].fix(58 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_I"].fix(92 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_S"].fix(363 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_BH"].fix(50 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_BA"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_P"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_O"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_NO"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_NH"].fix(23 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_ND"].fix(5 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_ND"].fix(16 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.alkalinity.fix(7 * pyo.units.mol / pyo.units.m**3)

    # Reactor sizing
    m.fs.R1.volume.fix(1000 * pyo.units.m**3)
    m.fs.R2.volume.fix(1000 * pyo.units.m**3)
    m.fs.R3.volume.fix(1333 * pyo.units.m**3)
    m.fs.R4.volume.fix(1333 * pyo.units.m**3)
    m.fs.R5.volume.fix(1333 * pyo.units.m**3)

    # Injection rates to Reactions 3, 4 and 5
    for j in m.fs.props_ASM1.component_list:
        if j != "S_O":
            # All components except S_O have no injection
            m.fs.R3.injection[:, :, j].fix(0)
            m.fs.R4.injection[:, :, j].fix(0)
            m.fs.R5.injection[:, :, j].fix(0)
    # Then set injections rates for O2
    m.fs.R3.outlet.conc_mass_comp[:, "S_O"].fix(1.72e-3)
    m.fs.R4.outlet.conc_mass_comp[:, "S_O"].fix(2.43e-3)
    m.fs.R5.outlet.conc_mass_comp[:, "S_O"].fix(4.49e-4)

    # Set fraction of outflow from reactor 5 that goes to recycle
    m.fs.SP5.split_fraction[:, "underflow"].fix(0.6)

    # Clarifier
    # TODO: Update once more detailed model available
    m.fs.CL1.split_fraction[0, "effluent", "H2O"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "S_I"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "S_S"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "X_I"].fix(0.00187)
    m.fs.CL1.split_fraction[0, "effluent", "X_S"].fix(0.00187)
    m.fs.CL1.split_fraction[0, "effluent", "X_BH"].fix(0.00187)
    m.fs.CL1.split_fraction[0, "effluent", "X_BA"].fix(0.00187)
    m.fs.CL1.split_fraction[0, "effluent", "X_P"].fix(0.00187)
    m.fs.CL1.split_fraction[0, "effluent", "S_O"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "S_NO"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "S_NH"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "S_ND"].fix(0.48956)
    m.fs.CL1.split_fraction[0, "effluent", "X_ND"].fix(0.00187)
    m.fs.CL1.split_fraction[0, "effluent", "S_ALK"].fix(0.48956)

    # Sludge purge separator
    m.fs.SP6.split_fraction[:, "recycle"].fix(0.985)

    # Outlet pressure from recycle pump
    m.fs.P1.outlet.pressure.fix(101325)

    m.fs.asm_adm = Translator_ASM1_ADM1(
        inlet_property_package=m.fs.props_ASM1,
        outlet_property_package=m.fs.props_ADM1,
        reaction_package=m.fs.ADM1_rxn_props,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )

    m.fs.RADM = AD(
        liquid_property_package=m.fs.props_ADM1,
        vapor_property_package=m.fs.props_vap,
        reaction_package=m.fs.ADM1_rxn_props,
        has_heat_transfer=True,
        has_pressure_change=False,
    )

    m.fs.adm_asm = Translator_ADM1_ASM1(
        inlet_property_package=m.fs.props_ADM1,
        outlet_property_package=m.fs.props_ASM1,
        reaction_package=m.fs.ADM1_rxn_props,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )

    m.fs.CL = Separator(
        property_package=m.fs.props_ASM1,
        outlet_list=["underflow", "effluent"],
        split_basis=SplittingType.componentFlow,
    )

    m.fs.TU = Thickener(property_package=m.fs.props_ASM1)
    m.fs.DU = DewateringUnit(property_package=m.fs.props_ASM1)

    m.fs.MX2 = Mixer(
        property_package=m.fs.props_ASM1, inlet_list=["feed_water1", "recycle1"]
    )
    m.fs.MX3 = Mixer(
        property_package=m.fs.props_ASM1, inlet_list=["feed_water2", "recycle2"]
    )
    m.fs.MX4 = Mixer(
        property_package=m.fs.props_ASM1, inlet_list=["thickener", "clarifier"]
    )

    # Clarifier
    # TODO: Update once more detailed model available
    m.fs.CL.split_fraction[0, "effluent", "H2O"].fix(0.993)
    m.fs.CL.split_fraction[0, "effluent", "S_I"].fix(0.993)
    m.fs.CL.split_fraction[0, "effluent", "S_S"].fix(0.993)
    m.fs.CL.split_fraction[0, "effluent", "X_I"].fix(0.5192)
    m.fs.CL.split_fraction[0, "effluent", "X_S"].fix(0.5192)
    m.fs.CL.split_fraction[0, "effluent", "X_BH"].fix(0.5192)
    m.fs.CL.split_fraction[0, "effluent", "X_BA"].fix(0.5192)
    m.fs.CL.split_fraction[0, "effluent", "X_P"].fix(0.5192)
    m.fs.CL.split_fraction[0, "effluent", "S_O"].fix(0.993)
    m.fs.CL.split_fraction[0, "effluent", "S_NO"].fix(0.993)
    m.fs.CL.split_fraction[0, "effluent", "S_NH"].fix(0.993)
    m.fs.CL.split_fraction[0, "effluent", "S_ND"].fix(0.993)
    m.fs.CL.split_fraction[0, "effluent", "X_ND"].fix(0.5192)
    m.fs.CL.split_fraction[0, "effluent", "S_ALK"].fix(0.993)
    m.fs.RADM.volume_liquid.fix(3400)
    m.fs.RADM.volume_vapor.fix(300)
    m.fs.RADM.liquid_outlet.temperature.fix(308.15)

    # Apply scaling

    m.fs.stream2adm = Arc(
        source=m.fs.RADM.liquid_outlet, destination=m.fs.adm_asm.inlet
    )
    m.fs.stream6adm = Arc(source=m.fs.SP6.waste, destination=m.fs.TU.inlet)
    m.fs.stream3adm = Arc(source=m.fs.TU.underflow, destination=m.fs.MX4.thickener)
    m.fs.stream7adm = Arc(source=m.fs.TU.overflow, destination=m.fs.MX3.recycle2)
    m.fs.stream9adm = Arc(source=m.fs.CL.underflow, destination=m.fs.MX4.clarifier)
    m.fs.stream4adm = Arc(source=m.fs.adm_asm.outlet, destination=m.fs.DU.inlet)
    m.fs.stream5adm = Arc(source=m.fs.DU.overflow, destination=m.fs.MX2.recycle1)
    m.fs.stream01 = Arc(source=m.fs.FeedWater.outlet, destination=m.fs.MX2.feed_water1)
    m.fs.stream02 = Arc(source=m.fs.MX2.outlet, destination=m.fs.MX3.feed_water2)
    m.fs.stream03 = Arc(source=m.fs.MX3.outlet, destination=m.fs.CL.inlet)
    m.fs.stream04 = Arc(source=m.fs.CL.effluent, destination=m.fs.MX1.feed_water)
    m.fs.stream10adm = Arc(source=m.fs.MX4.outlet, destination=m.fs.asm_adm.inlet)
    m.fs.stream1adm = Arc(source=m.fs.asm_adm.outlet, destination=m.fs.RADM.inlet)
    pyo.TransformationFactory("network.expand_arcs").apply_to(m)
    
    m.fs.set_default_scaling("conc_mass_comp", 1)
    iscale.calculate_scaling_factors(m.fs)
    
    print("Degree of Freedom is ", degrees_of_freedom(m))

    # Initialize flowsheet
    # Apply sequential decomposition - 1 iteration should suffice
    seq = SequentialDecomposition()
    seq.options.select_tear_method = "heuristic"
    seq.options.tear_method = "Wegstein"
    seq.options.iterLim = 1
    # seq.options.tear_set = [m.fs.stream2, m.fs.stream10adm]
    # seq.options.tol = 1e-5
    seq.options.solve_tear = False

    G = seq.create_graph(m)
    # # Uncomment this code to see tear set and initialization order
    heuristic_tear_set = seq.tear_set_arcs(G, method="heuristic")
    order = seq.calculation_order(G)
    for o in heuristic_tear_set:
        print(o.name)
    print("order")    
    for o in order:
        print(o[0].name)

    # Initial guesses for flow into first reactor
    tear_guesses = {
        "flow_vol": {0: 103531/24/3600},
        "conc_mass_comp": {
            (0, "S_I"): 0.028,
            (0, "S_S"): 0.012,
            (0, "X_I"): 1.532,
            (0, "X_S"): 0.069,
            (0, "X_BH"): 2.233,
            (0, "X_BA"): 0.167,
            (0, "X_P"): 0.964,
            (0, "S_O"): 0.0011,
            (0, "S_NO"): 0.0073,
            (0, "S_NH"): 0.0072,
            (0, "S_ND"): 0.0016,
            (0, "X_ND"): 0.0040,
        },
        "alkalinity": {0: 0.0052},
        "temperature": {0: 308.15},
        "pressure": {0: 101325},
    }
    
    tear_guesses1 = {
        "flow_vol": {0: 147/24/3600},
        "conc_mass_comp": {
            (0, "S_I"): 0.028,
            (0, "S_S"): 0.059,
            (0, "X_I"): 6.480,
            (0, "X_S"): 24.509,
            (0, "X_BH"): 3.495,
            (0, "X_BA"): 0.0065,
            (0, "X_P"): 0.045,
            (0, "S_O"): 0.000017,
            (0, "S_NO"): 0.000117,
            (0, "S_NH"): 0.035,
            (0, "S_ND"): 0.0055,
            (0, "X_ND"): 1.086,
        },
        "alkalinity": {0: 0.0077},
        "temperature": {0: 308.15},
        "pressure": {0: 101325},
    }
    
    tear_guesses2 = {
        "flow_vol": {0: 178/24/3600},
        "conc_mass_comp": {
            (0, "S_I"): 0.028,
            (0, "S_S"): 0.048,
            (0, "X_I"): 10.362,
            (0, "X_S"): 20.375,
            (0, "X_BH"): 10.210,
            (0, "X_BA"): 0.553,
            (0, "X_P"): 3.204,
            (0, "S_O"): 0.00025,
            (0, "S_NO"): 0.00169,
            (0, "S_NH"): 0.0289,
            (0, "S_ND"): 0.00468,
            (0, "X_ND"): 0.906,
        },
        "alkalinity": {0: 0.00715},
        "temperature": {0: 308.15},
        "pressure": {0: 101325},
    }
    
    tear_guesses3 = {
        "flow_vol": {0: 178/24/3600},
        "conc_mass_comp": {
            (0, "S_su"): 0,
            (0, "S_aa"): 0.04388,
            (0, "S_fa"): 0,
            (0, "S_va"): 0,
            (0, "S_bu"): 0,
            (0, "S_pro"): 0,
            (0, "S_ac"): 0,
            (0, "S_h2"): 0,
            (0, "S_ch4"): 0,
            (0, "S_IC"): 0.0079326*12,
            (0, "S_IN"): 0.0019721*14,
            (0, "S_I"): 0.028067,
            (0, "X_c"): 0,
            (0, "X_ch"): 3.7236,
            (0, "X_pr"): 15.9235,
            (0, "X_li"): 8.047,
            (0, "X_su"): 0,
            (0, "X_aa"): 0,
            (0, "X_fa"): 0,
            (0, "X_c4"): 0,
            (0, "X_pro"): 0,
            (0, "X_ac"): 0,
            (0, "X_h2"): 0,
            (0, "X_I"): 17.0106,
        },
        "cations": {0: 0},
        "anions": {0: 0.005},
        "temperature": {0: 308.15},
        "pressure": {0: 101325},
    }
    
    tear_guesses4 = {
        "flow_vol": {0: 178/24/3600},
        "conc_mass_comp": {
            (0, "S_su"): 0,
            (0, "S_aa"): 0.04388,
            (0, "S_fa"): 0,
            (0, "S_va"): 0,
            (0, "S_bu"): 0,
            (0, "S_pro"): 0,
            (0, "S_ac"): 0,
            (0, "S_h2"): 0,
            (0, "S_ch4"): 0,
            (0, "S_IC"): 0.0079326*12,
            (0, "S_IN"): 0.0019721*14,
            (0, "S_I"): 0.028067,
            (0, "X_c"): 0,
            (0, "X_ch"): 3.7236,
            (0, "X_pr"): 15.9235,
            (0, "X_li"): 8.047,
            (0, "X_su"): 0,
            (0, "X_aa"): 0,
            (0, "X_fa"): 0,
            (0, "X_c4"): 0,
            (0, "X_pro"): 0,
            (0, "X_ac"): 0,
            (0, "X_h2"): 0,
            (0, "X_I"): 17.0106,
        },
        "cations": {0: 0},
        "anions": {0: 0.005},
        "temperature": {0: 308.15},
        "pressure": {0: 101325},
    }

    # Pass the tear_guess to the SD tool
    seq.set_guesses_for(m.fs.R1.inlet, tear_guesses)
    seq.set_guesses_for(m.fs.MX4.clarifier, tear_guesses1)
    # seq.set_guesses_for(m.fs.asm_adm.inlet, tear_guesses2)
    # seq.set_guesses_for(m.fs.asm_adm.outlet, tear_guesses3)
    # seq.set_guesses_for(m.fs.RADM.inlet, tear_guesses4)
    
    def function(unit):
        if unit == m.fs.RADM:
            # i = 0
            unit.initialize(outlvl=idaeslog.DEBUG,optarg={"halt_on_ampl_error": "yes","bound_push": 1e-8})
            print("Degree of Freedom is ", degrees_of_freedom(unit))
            # unit.pprint()
            # results = solver.solve(unit, tee=True)
        elif unit == m.fs.asm_adm:
            # unit.pprint()
            unit.initialize(outlvl=idaeslog.DEBUG)
            print("Degree of Freedom is ", degrees_of_freedom(unit))

            # results = solver.solve(unit, tee=True)
        else:
            unit.initialize(outlvl=idaeslog.DEBUG)
            print("Degree of Freedom is ", degrees_of_freedom(unit))
            # unit.pprint()
            # results = solver.solve(unit, tee=True)

    seq.run(m, function)
    print(degrees_of_freedom(m))
    solver = get_solver(options={"halt_on_ampl_error": "yes","tol": 1e-6, "constr_viol_tol": 1e-6, "bound_push": 1e-8})
    results = solver.solve(m, tee=True)
    m.obj = Objective(expr=0)
    dh = DegeneracyHunter(m, solver=SolverFactory('cbc'))
    dh.check_residuals(tol=1E-14) 
    report_scaling_issues(m)
    
    # pyo.assert_optimal_termination(results)
    # m.display()

    print(large_residuals_set(m))

    return m, results