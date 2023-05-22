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

from watertap.examples.flowsheets.case_studies.activated_sludge.ASM1_flowsheet import (
    build_flowsheet as build_ASM,
)

from watertap.core.util.initialization import check_solve


def build_flowsheet():
    # Call ASM model
    m, res = build_ASM()

    m.del_component(m.fs.Sludge)
    m.del_component(m.fs.stream1)
    m.del_component(m.fs.stream102)

    m.fs.props_ASM1 = ASM1ParameterBlock()
    m.fs.props_ADM1 = ADM1ParameterBlock()
    m.fs.props_vap = ADM1_vaporParameterBlock()
    m.fs.ADM1_rxn_props = ADM1ReactionParameterBlock(property_package=m.fs.props_ADM1)

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
        property_package=m.fs.props,
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

    m.fs.FeedWater.temperature.fix(308.15 * pyo.units.K)
    m.fs.SP6.split_fraction[:, "recycle"].fix(0.975)

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

    m.fs.stream2adm = Arc(
        source=m.fs.RADM.liquid_outlet, destination=m.fs.adm_asm.inlet
    )
    pyo.TransformationFactory("network.expand_arcs").apply_to(m)

    m.fs.asm_adm.inlet.flow_vol.fix(170 * pyo.units.m**3 / pyo.units.day)
    m.fs.asm_adm.inlet.temperature.fix(308.15 * units.K)
    m.fs.asm_adm.inlet.pressure.fix(1 * units.atm)

    m.fs.asm_adm.inlet.conc_mass_comp[0, "S_I"].fix(30 * pyo.units.mg / pyo.units.liter)
    m.fs.asm_adm.inlet.conc_mass_comp[0, "S_S"].fix(
        8.89e-1 * pyo.units.mg / pyo.units.liter
    )
    m.fs.asm_adm.inlet.conc_mass_comp[0, "X_I"].fix(
        2247 * pyo.units.mg / pyo.units.liter
    )
    m.fs.asm_adm.inlet.conc_mass_comp[0, "X_S"].fix(
        96.8 * pyo.units.mg / pyo.units.liter
    )
    m.fs.asm_adm.inlet.conc_mass_comp[0, "X_BH"].fix(
        5004 * pyo.units.mg / pyo.units.liter
    )
    m.fs.asm_adm.inlet.conc_mass_comp[0, "X_BA"].fix(
        292 * pyo.units.mg / pyo.units.liter
    )
    m.fs.asm_adm.inlet.conc_mass_comp[0, "X_P"].fix(
        884 * pyo.units.mg / pyo.units.liter
    )
    m.fs.asm_adm.inlet.conc_mass_comp[0, "S_O"].fix(
        4.49e-1 * pyo.units.mg / pyo.units.liter
    )
    m.fs.asm_adm.inlet.conc_mass_comp[0, "S_NO"].fix(
        10.14 * pyo.units.mg / pyo.units.liter
    )
    m.fs.asm_adm.inlet.conc_mass_comp[0, "S_NH"].fix(
        1.86e-1 * pyo.units.mg / pyo.units.liter
    )
    m.fs.asm_adm.inlet.conc_mass_comp[0, "S_ND"].fix(
        6.88e-1 * pyo.units.mg / pyo.units.liter
    )
    m.fs.asm_adm.inlet.conc_mass_comp[0, "X_ND"].fix(
        6.92 * pyo.units.mg / pyo.units.liter
    )
    m.fs.asm_adm.inlet.alkalinity.fix(4.13 * units.mol / units.m**3)

    m.fs.RADM.inlet.flow_vol.fix(170 * pyo.units.m**3 / pyo.units.day)
    m.fs.RADM.inlet.temperature.fix(308.15 * pyo.units.K)
    m.fs.RADM.inlet.pressure.fix(1 * pyo.units.atm)
    m.fs.RADM.inlet.conc_mass_comp[0, "S_su"].fix(1 * pyo.units.mg / pyo.units.liter)
    m.fs.RADM.inlet.conc_mass_comp[0, "S_aa"].fix(1 * pyo.units.mg / pyo.units.liter)
    m.fs.RADM.inlet.conc_mass_comp[0, "S_fa"].fix(1 * pyo.units.mg / pyo.units.liter)
    m.fs.RADM.inlet.conc_mass_comp[0, "S_va"].fix(1 * pyo.units.mg / pyo.units.liter)
    m.fs.RADM.inlet.conc_mass_comp[0, "S_bu"].fix(1 * pyo.units.mg / pyo.units.liter)
    m.fs.RADM.inlet.conc_mass_comp[0, "S_pro"].fix(1 * pyo.units.mg / pyo.units.liter)
    m.fs.RADM.inlet.conc_mass_comp[0, "S_ac"].fix(1 * pyo.units.mg / pyo.units.liter)
    m.fs.RADM.inlet.conc_mass_comp[0, "S_h2"].fix(1e-5 * pyo.units.mg / pyo.units.liter)
    m.fs.RADM.inlet.conc_mass_comp[0, "S_ch4"].fix(1 * pyo.units.mg / pyo.units.liter)
    m.fs.RADM.inlet.conc_mass_comp[0, "S_IC"].fix(
        40 * units.mmol / units.liter * 12 * units.mg / units.mmol
    )
    m.fs.RADM.inlet.conc_mass_comp[0, "S_IN"].fix(
        10 * units.mmol / units.liter * 14 * units.mg / units.mmol
    )
    m.fs.RADM.inlet.conc_mass_comp[0, "S_I"].fix(20 * pyo.units.mg / pyo.units.liter)

    m.fs.RADM.inlet.conc_mass_comp[0, "X_c"].fix(2000 * pyo.units.mg / pyo.units.liter)
    m.fs.RADM.inlet.conc_mass_comp[0, "X_ch"].fix(5000 * pyo.units.mg / pyo.units.liter)
    m.fs.RADM.inlet.conc_mass_comp[0, "X_pr"].fix(
        20000 * pyo.units.mg / pyo.units.liter
    )
    m.fs.RADM.inlet.conc_mass_comp[0, "X_li"].fix(5000 * pyo.units.mg / pyo.units.liter)
    m.fs.RADM.inlet.conc_mass_comp[0, "X_su"].fix(1 * pyo.units.mg / pyo.units.liter)
    m.fs.RADM.inlet.conc_mass_comp[0, "X_aa"].fix(10 * pyo.units.mg / pyo.units.liter)
    m.fs.RADM.inlet.conc_mass_comp[0, "X_fa"].fix(10 * pyo.units.mg / pyo.units.liter)
    m.fs.RADM.inlet.conc_mass_comp[0, "X_c4"].fix(10 * pyo.units.mg / pyo.units.liter)
    m.fs.RADM.inlet.conc_mass_comp[0, "X_pro"].fix(10 * pyo.units.mg / pyo.units.liter)
    m.fs.RADM.inlet.conc_mass_comp[0, "X_ac"].fix(10 * pyo.units.mg / pyo.units.liter)
    m.fs.RADM.inlet.conc_mass_comp[0, "X_h2"].fix(10 * pyo.units.mg / pyo.units.liter)
    m.fs.RADM.inlet.conc_mass_comp[0, "X_I"].fix(25000 * pyo.units.mg / pyo.units.liter)

    m.fs.RADM.inlet.cations[0].fix(40 * pyo.units.mmol / pyo.units.liter)
    m.fs.RADM.inlet.anions[0].fix(20 * pyo.units.mmol / pyo.units.liter)

    m.fs.RADM.volume_liquid.fix(3400 * pyo.units.m**3)
    m.fs.RADM.volume_vapor.fix(300 * pyo.units.m**3)

    m.fs.RADM.liquid_outlet.temperature.fix(308.15 * pyo.units.K)
    print(degrees_of_freedom(m))

    # Apply scaling
    iscale.calculate_scaling_factors(m.fs)

    m.fs.asm_adm.initialize(outlvl=idaeslog.INFO_HIGH, optarg={"bound_push": 1e-8})
    m.fs.RADM.initialize(outlvl=idaeslog.INFO_HIGH, optarg={"bound_push": 1e-8})
    propagate_state(m.fs.stream2adm)
    m.fs.adm_asm.initialize(outlvl=idaeslog.INFO_HIGH, optarg={"bound_push": 1e-8})

    print(degrees_of_freedom(m))
    # assert degrees_of_freedom(m) == 0

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

    pyo.TransformationFactory("network.expand_arcs").apply_to(m)

    propagate_state(m.fs.stream6adm)

    m.fs.TU.initialize(outlvl=idaeslog.INFO_HIGH, optarg={"bound_push": 1e-8})
    propagate_state(m.fs.stream3adm)
    propagate_state(m.fs.stream7adm)

    m.fs.DU.initialize(outlvl=idaeslog.INFO_HIGH, optarg={"bound_push": 1e-8})
    propagate_state(m.fs.stream5adm)

    m.fs.FeedWater.flow_vol.fix(20648 * pyo.units.m**3 / pyo.units.day)
    propagate_state(m.fs.stream01)
    m.fs.MX2.initialize(outlvl=idaeslog.INFO_HIGH, optarg={"bound_push": 1e-8})
    propagate_state(m.fs.stream02)

    m.fs.MX3.initialize(outlvl=idaeslog.INFO_HIGH, optarg={"bound_push": 1e-8})
    propagate_state(m.fs.stream03)

    m.fs.CL.initialize(outlvl=idaeslog.INFO_HIGH, optarg={"bound_push": 1e-8})
    propagate_state(m.fs.stream04)

    propagate_state(m.fs.stream9adm)
    m.fs.MX4.initialize(outlvl=idaeslog.INFO_HIGH, optarg={"bound_push": 1e-8})
    propagate_state(m.fs.stream4adm)

    m.fs.RADM.inlet.flow_vol.unfix()
    m.fs.RADM.inlet.temperature.unfix()
    m.fs.RADM.inlet.pressure.unfix()
    m.fs.RADM.inlet.conc_mass_comp[0, "S_su"].unfix()
    m.fs.RADM.inlet.conc_mass_comp[0, "S_aa"].unfix()
    m.fs.RADM.inlet.conc_mass_comp[0, "S_fa"].unfix()
    m.fs.RADM.inlet.conc_mass_comp[0, "S_va"].unfix()
    m.fs.RADM.inlet.conc_mass_comp[0, "S_bu"].unfix()
    m.fs.RADM.inlet.conc_mass_comp[0, "S_pro"].unfix()
    m.fs.RADM.inlet.conc_mass_comp[0, "S_ac"].unfix()
    m.fs.RADM.inlet.conc_mass_comp[0, "S_h2"].unfix()
    m.fs.RADM.inlet.conc_mass_comp[0, "S_ch4"].unfix()
    m.fs.RADM.inlet.conc_mass_comp[0, "S_IC"].unfix()
    m.fs.RADM.inlet.conc_mass_comp[0, "S_IN"].unfix()
    m.fs.RADM.inlet.conc_mass_comp[0, "S_I"].unfix()
    m.fs.RADM.inlet.conc_mass_comp[0, "X_c"].unfix()
    m.fs.RADM.inlet.conc_mass_comp[0, "X_ch"].unfix()
    m.fs.RADM.inlet.conc_mass_comp[0, "X_pr"].unfix()
    m.fs.RADM.inlet.conc_mass_comp[0, "X_li"].unfix()
    m.fs.RADM.inlet.conc_mass_comp[0, "X_su"].unfix()
    m.fs.RADM.inlet.conc_mass_comp[0, "X_aa"].unfix()
    m.fs.RADM.inlet.conc_mass_comp[0, "X_fa"].unfix()
    m.fs.RADM.inlet.conc_mass_comp[0, "X_c4"].unfix()
    m.fs.RADM.inlet.conc_mass_comp[0, "X_pro"].unfix()
    m.fs.RADM.inlet.conc_mass_comp[0, "X_ac"].unfix()
    m.fs.RADM.inlet.conc_mass_comp[0, "X_h2"].unfix()
    m.fs.RADM.inlet.conc_mass_comp[0, "X_I"].unfix()
    m.fs.RADM.inlet.cations[0].unfix()
    m.fs.RADM.inlet.anions[0].unfix()

    m.fs.stream10adm = Arc(source=m.fs.MX4.outlet, destination=m.fs.asm_adm.inlet)
    pyo.TransformationFactory("network.expand_arcs").apply_to(m)

    m.fs.asm_adm.inlet.flow_vol.unfix()
    m.fs.asm_adm.inlet.temperature.unfix()
    m.fs.asm_adm.inlet.pressure.unfix()
    m.fs.asm_adm.inlet.conc_mass_comp[0, "S_I"].unfix()
    m.fs.asm_adm.inlet.conc_mass_comp[0, "S_S"].unfix()
    m.fs.asm_adm.inlet.conc_mass_comp[0, "X_I"].unfix()
    m.fs.asm_adm.inlet.conc_mass_comp[0, "X_S"].unfix()
    m.fs.asm_adm.inlet.conc_mass_comp[0, "X_BH"].unfix()
    m.fs.asm_adm.inlet.conc_mass_comp[0, "X_BA"].unfix()
    m.fs.asm_adm.inlet.conc_mass_comp[0, "X_P"].unfix()
    m.fs.asm_adm.inlet.conc_mass_comp[0, "S_O"].unfix()
    m.fs.asm_adm.inlet.conc_mass_comp[0, "S_NO"].unfix()
    m.fs.asm_adm.inlet.conc_mass_comp[0, "S_NH"].unfix()
    m.fs.asm_adm.inlet.conc_mass_comp[0, "S_ND"].unfix()
    m.fs.asm_adm.inlet.conc_mass_comp[0, "X_ND"].unfix()
    m.fs.asm_adm.inlet.alkalinity.unfix()

    print("3", degrees_of_freedom(m))
    m.fs.stream1adm = Arc(source=m.fs.asm_adm.outlet, destination=m.fs.RADM.inlet)
    pyo.TransformationFactory("network.expand_arcs").apply_to(m)
    print("3", degrees_of_freedom(m))
    # Initialize flowsheet
    # Apply sequential decomposition - 1 iteration should suffice
    # seq = SequentialDecomposition()
    # seq.options.select_tear_method = "heuristic"
    # seq.options.tear_method = "Wegstein"
    # seq.options.iterLim = 1
    #
    # G = seq.create_graph(m)
    #
    # def function(unit):
    #     unit.initialize(outlvl=idaeslog.INFO_HIGH, optarg={"bound_push": 1e-8})
    #
    # seq.run(m, function)

    # m.fs.pprint()
    # print(large_residuals_set(m.fs))

    # for var, sv in iscale.badly_scaled_var_generator(m):
    #     print(var,sv)

    # solver = get_solver(options={"bound_push": 1e-8,"max_iter":0})
    solver = get_solver(options={"bound_push": 1e-8})
    results = solver.solve(m, tee=True)
    print(degrees_of_freedom(m))
    # pyo.assert_optimal_termination(results)
    # m.display()

    print(large_residuals_set(m))

    return m, results
