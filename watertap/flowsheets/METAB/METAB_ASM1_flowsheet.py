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
This flowsheet is a WWTP model where METAB is integrated with the ASM1 flowsheet

"""
__author__ = "Marcus Holly"

import pyomo.environ as pyo
from pyomo.environ import units
from pyomo.network import Arc, SequentialDecomposition

from idaes.core import FlowsheetBlock
from idaes.core.util import DiagnosticsToolbox
from idaes.models.unit_models import (
    CSTR,
    Feed,
    Mixer,
    Separator,
    PressureChanger,
    Product,
)
from idaes.models.unit_models.mixer import MomentumMixingType
from idaes.models.unit_models.separator import SplittingType
from watertap.core.solvers import get_solver
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale
from idaes.core.util.tables import (
    create_stream_table_dataframe,
    stream_table_dataframe_to_string,
)

from watertap.unit_models.cstr_injection import CSTR_Injection
from watertap.property_models.unit_specific.activated_sludge.asm1_properties import (
    ASM1ParameterBlock,
)
from watertap.property_models.unit_specific.activated_sludge.asm1_reactions import (
    ASM1ReactionParameterBlock,
)
from watertap.property_models.unit_specific.anaerobic_digestion.adm1_properties import (
    ADM1ParameterBlock,
)
from watertap.property_models.unit_specific.anaerobic_digestion.adm1_reactions import (
    ADM1ReactionParameterBlock,
)
from watertap.unit_models.translators.translator_adm1_asm1 import Translator_ADM1_ASM1

from watertap.flowsheets.METAB import model_connector as metab
from watertap.core.util.initialization import check_solve

# Set up logger
_log = idaeslog.getLogger(__name__)


def main():
    m = build()
    set_operating_conditions(m)
    scale_system(m)
    initialize_system(m)
    results = solve_system(m)
    report_st(m)

    return m, results


def build():
    m = pyo.ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.props_asm1 = ASM1ParameterBlock()
    m.fs.props_adm1 = ADM1ParameterBlock()
    m.fs.asm1_rxn_props = ASM1ReactionParameterBlock(property_package=m.fs.props_asm1)
    m.fs.adm1_rxn_props = ADM1ReactionParameterBlock(property_package=m.fs.props_adm1)
    # Mixer for feed water and recycled sludge
    m.fs.MX1 = Mixer(
        property_package=m.fs.props_asm1,
        inlet_list=["feed_water", "recycle"],
        momentum_mixing_type=MomentumMixingType.none,
    )
    # First reactor (anoxic) - standard CSTR
    m.fs.R1 = CSTR(
        property_package=m.fs.props_asm1, reaction_package=m.fs.asm1_rxn_props
    )
    # Second reactor (anoxic) - standard CSTR
    m.fs.R2 = CSTR(
        property_package=m.fs.props_asm1, reaction_package=m.fs.asm1_rxn_props
    )
    # Third reactor (aerobic) - CSTR with injection
    m.fs.R3 = CSTR_Injection(
        property_package=m.fs.props_asm1,
        reaction_package=m.fs.asm1_rxn_props,
    )
    # Fourth reactor (aerobic) - CSTR with injection
    m.fs.R4 = CSTR_Injection(
        property_package=m.fs.props_asm1, reaction_package=m.fs.asm1_rxn_props
    )
    # Fifth reactor (aerobic) - CSTR with injection
    m.fs.R5 = CSTR_Injection(
        property_package=m.fs.props_asm1, reaction_package=m.fs.asm1_rxn_props
    )
    m.fs.SP5 = Separator(
        property_package=m.fs.props_asm1, outlet_list=["underflow", "overflow"]
    )
    # Clarifier
    m.fs.CL1 = Separator(
        property_package=m.fs.props_asm1,
        outlet_list=["underflow", "effluent"],
        split_basis=SplittingType.componentFlow,
    )
    # Sludge purge splitter
    m.fs.SP6 = Separator(
        property_package=m.fs.props_asm1,
        outlet_list=["recycle", "waste"],
        split_basis=SplittingType.totalFlow,
    )
    # Mixing sludge recycle and R5 underflow
    m.fs.MX6 = Mixer(
        property_package=m.fs.props_asm1,
        inlet_list=["clarifier", "reactor"],
        momentum_mixing_type=MomentumMixingType.none,
    )
    # Product Blocks
    m.fs.Treated = Product(property_package=m.fs.props_asm1)
    m.fs.Sludge = Product(property_package=m.fs.props_asm1)
    # Recycle pressure changer - use a simple isothermal unit for now
    m.fs.P1 = PressureChanger(property_package=m.fs.props_asm1)

    m.fs.metab_effluent = Translator_ADM1_ASM1(
        inlet_property_package=m.fs.props_adm1,
        outlet_property_package=m.fs.props_asm1,
        reaction_package=m.fs.adm1_rxn_props,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )

    m.fs.mixers = (m.fs.MX1, m.fs.MX6)
    for mixer in m.fs.mixers:
        mixer.outlet.pressure.fix()

    # Link units
    m.fs.stream1 = Arc(
        source=m.fs.metab_effluent.outlet, destination=m.fs.MX1.feed_water
    )
    m.fs.stream2 = Arc(source=m.fs.MX1.outlet, destination=m.fs.R1.inlet)
    m.fs.stream3 = Arc(source=m.fs.R1.outlet, destination=m.fs.R2.inlet)
    m.fs.stream4 = Arc(source=m.fs.R2.outlet, destination=m.fs.R3.inlet)
    m.fs.stream5 = Arc(source=m.fs.R3.outlet, destination=m.fs.R4.inlet)
    m.fs.stream6 = Arc(source=m.fs.R4.outlet, destination=m.fs.R5.inlet)
    m.fs.stream7 = Arc(source=m.fs.R5.outlet, destination=m.fs.SP5.inlet)
    m.fs.stream8 = Arc(source=m.fs.SP5.overflow, destination=m.fs.CL1.inlet)
    m.fs.stream9 = Arc(source=m.fs.SP5.underflow, destination=m.fs.MX6.reactor)
    m.fs.stream10 = Arc(source=m.fs.CL1.effluent, destination=m.fs.Treated.inlet)
    m.fs.stream11 = Arc(source=m.fs.CL1.underflow, destination=m.fs.SP6.inlet)
    m.fs.stream102 = Arc(source=m.fs.SP6.waste, destination=m.fs.Sludge.inlet)
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

    return m


def set_operating_conditions(m):
    # Reactor sizing
    m.fs.R1.volume.fix(1000 * 1e-5 * pyo.units.m**3)
    m.fs.R2.volume.fix(1000 * 1e-5 * pyo.units.m**3)
    m.fs.R3.volume.fix(1333 * 1e-5 * pyo.units.m**3)
    m.fs.R4.volume.fix(1333 * 1e-5 * pyo.units.m**3)
    m.fs.R5.volume.fix(1333 * 1e-5 * pyo.units.m**3)

    # Injection rates to Reactions 3, 4 and 5
    for j in m.fs.props_asm1.component_list:
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
    m.fs.SP6.split_fraction[:, "recycle"].fix(0.97955)

    m.fs.R3.electricity_consumption[0].set_value(20)
    m.fs.R4.electricity_consumption[0].set_value(20)
    m.fs.R5.electricity_consumption[0].set_value(20)

    # Outlet pressure from recycle pump
    m.fs.P1.outlet.pressure.fix(101325)

    idx = 1
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
        "X_su",
        "X_aa",
        "X_fa",
        "X_c4",
        "X_pro",
        "X_ac",
        "X_h2",
        "X_I",
        "S_cat",
        "S_an",
    ]  # 'S_co2' need to fix

    df_metab = metab.get_influent(m)
    m.fs.metab_effluent.inlet.flow_vol.fix(
        float(metab.get_outlet(df_metab, idx, "VolumetricFlowrate"))
        * units.m**3
        / units.day
    )
    m.fs.metab_effluent.inlet.temperature.fix(308.15 * units.K)
    m.fs.metab_effluent.inlet.pressure.fix(1 * units.atm)

    for i in comp_list:
        if i == "S_cat":
            m.fs.metab_effluent.inlet.cations[0].fix(
                float(metab.get_outlet(df_metab, idx, i)) * units.mmol / units.liter
            )

        elif i == "S_an":
            m.fs.metab_effluent.inlet.anions[0].fix(
                float(metab.get_outlet(df_metab, idx, i)) * units.mmol / units.liter
            )
        elif i != "H2O":
            m.fs.metab_effluent.inlet.conc_mass_comp[0, i].fix(
                metab.get_outlet(df_metab, idx, i) * units.mg / units.liter
            )


def scale_system(m):
    # Apply scaling
    for var in m.fs.component_data_objects(pyo.Var, descend_into=True):
        if "flow_vol" in var.name:
            iscale.set_scaling_factor(var, 1e2)
        if "temperature" in var.name:
            iscale.set_scaling_factor(var, 1e-1)
        if "pressure" in var.name:
            iscale.set_scaling_factor(var, 1e-6)
        if "conc_mass_comp" in var.name:
            iscale.set_scaling_factor(var, 1e1)
    iscale.calculate_scaling_factors(m.fs)


def initialize_system(m):
    seq = SequentialDecomposition()
    seq.options.select_tear_method = "heuristic"
    seq.options.tear_method = "Wegstein"
    seq.options.iterLim = 1

    G = seq.create_graph(m)

    # # Uncomment this code to see tear set and initialization order
    heuristic_tear_set = seq.tear_set_arcs(G, method="heuristic")
    order = seq.calculation_order(G)
    for o in heuristic_tear_set:
        print(o.name)
    for o in order:
        print(o[0].name)

    # Initial guesses for flow into first reactor
    tear_guesses = {
        "flow_vol": {0: 5.57e-4},
        "conc_mass_comp": {
            (0, "S_I"): 0.086,
            (0, "S_S"): 0.040,
            (0, "X_I"): 2.50,
            (0, "X_S"): 16.81,
            (0, "X_BH"): 3.14,
            (0, "X_BA"): 2.32e-9,
            (0, "X_P"): 2.26e-9,
            (0, "S_O"): 3.59e-4,
            (0, "S_NO"): 1.14e-10,
            (0, "S_NH"): 0.21,
            (0, "S_ND"): 5.25e-3,
            (0, "X_ND"): 1.30,
        },
        "alkalinity": {0: 0.047},
        "temperature": {0: 308.15},
        "pressure": {0: 101325},
    }

    # Pass the tear_guess to the SD tool
    seq.set_guesses_for(m.fs.R1.inlet, tear_guesses)

    def function(unit):
        unit.initialize(outlvl=idaeslog.INFO)

    seq.run(m, function)


def solve_system(m):
    # Solve overall flowsheet to close recycle loop
    solver = get_solver()
    results = solver.solve(m, tee=True)
    check_solve(results, checkpoint="closing recycle", logger=_log, fail_flag=True)

    # Switch to fixed KLa in R3 and R4 (S_O concentration is controlled in R5)
    m.fs.R3.KLa.fix(10)
    m.fs.R4.KLa.fix(10)
    m.fs.R3.outlet.conc_mass_comp[:, "S_O"].unfix()
    m.fs.R4.outlet.conc_mass_comp[:, "S_O"].unfix()

    # Resolve with controls in place
    results = solver.solve(m, tee=True)
    check_solve(
        results,
        checkpoint="re-solve with controls in place",
        logger=_log,
        fail_flag=True,
    )
    return results


def report_st(m):

    stream_table = create_stream_table_dataframe(
        {
            "MX1": m.fs.MX1.outlet,
            "R1": m.fs.R1.outlet,
            "R2": m.fs.R2.outlet,
            "R3": m.fs.R3.outlet,
            "R4": m.fs.R4.outlet,
            "R5": m.fs.R5.outlet,
            "from metab": m.fs.metab_effluent.inlet,
            "Influent": m.fs.metab_effluent.outlet,
        },
        time_point=0,
    )
    print(stream_table_dataframe_to_string(stream_table))
    stream_table.to_csv("influent_values.csv")


if __name__ == "__main__":
    m, results = main()
    m.fs.Treated.display()
    m.fs.Sludge.display()
