###############################################################################
# ProteusLib Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/nawi-hub/proteuslib/"
#
###############################################################################

"""Simple flowsheet with zero order separators for NF and RO examples"""

from pyomo.environ import ConcreteModel, Constraint, TransformationFactory
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models import Separator, Mixer
from idaes.generic_models.unit_models.separator import SplittingType, EnergySplittingType
from idaes.generic_models.unit_models.mixer import MixingType, MomentumMixingType
from idaes.generic_models.unit_models.translator import Translator
from proteuslib.flowsheets.full_treatment_train.example_models import unit_separator
from proteuslib.flowsheets.full_treatment_train.util import solve_with_user_scaling, check_dof
from idaes.core.util.scaling import calculate_scaling_factors


def build_flowsheet_NF_salt_basis_example(m):
    """
    Builds a flowsheet with NF and RO. Both the NF and RO are modeled as separators with
    specified split fractions. The flowsheet includes a translator block to convert between
    the NF's multi-salt property package and the RO's seawater package.
    """
    # build flowsheet
    unit_separator.build_NF_salt_separator_example(m)
    unit_separator.build_RO_separator_example(m)
    m.fs.tb_NF_to_RO = Translator(
        default={"inlet_property_package": m.fs.NF_properties,
                 "outlet_property_package": m.fs.RO_properties})
    m.fs.splitter = Separator(default={
            "property_package": m.fs.NF_properties,
            "outlet_list": ['pretreatment', 'bypass'],
            "split_basis": SplittingType.totalFlow,
            "energy_split_basis": EnergySplittingType.equal_temperature})
    m.fs.mixer = Mixer(default={
            "property_package": m.fs.NF_properties,
            "inlet_list": ['pretreatment', 'bypass']})

    # add constraints for translator block (tb)
    m.fs.tb_NF_to_RO.eq_H2O_balance = Constraint(
        expr=m.fs.tb_NF_to_RO.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O']
        == m.fs.tb_NF_to_RO.outlet.flow_mass_phase_comp[0, 'Liq', 'H2O'])
    m.fs.tb_NF_to_RO.eq_TDS_balance = Constraint(
        expr=sum(m.fs.tb_NF_to_RO.inlet.flow_mass_phase_comp[0, 'Liq', j] for j in ['NaCl', 'CaSO4', 'MgSO4', 'MgCl2'])
        == m.fs.tb_NF_to_RO.outlet.flow_mass_phase_comp[0, 'Liq', 'TDS'])
    m.fs.tb_NF_to_RO.eq_equal_temperature = Constraint(
        expr=m.fs.tb_NF_to_RO.inlet.temperature[0]
        == m.fs.tb_NF_to_RO.outlet.temperature[0])
    m.fs.tb_NF_to_RO.eq_equal_pressure = Constraint(
        expr=m.fs.tb_NF_to_RO.inlet.pressure[0]
        == m.fs.tb_NF_to_RO.outlet.pressure[0])

    # connect models
    m.fs.S1 = Arc(source=m.fs.splitter.bypass, destination=m.fs.mixer.bypass)
    m.fs.S2 = Arc(source=m.fs.splitter.pretreatment, destination=m.fs.NF.inlet)
    m.fs.S3 = Arc(source=m.fs.NF.permeate, destination=m.fs.mixer.pretreatment)
    m.fs.S4 = Arc(source=m.fs.mixer.outlet, destination=m.fs.tb_NF_to_RO.inlet)
    m.fs.S5 = Arc(source=m.fs.tb_NF_to_RO.outlet, destination=m.fs.RO.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # specify (NF and RO models are already specified, mixer has 0 DOF, splitter has 1 DOF)
    # feed
    feed_flow_mass = 1
    feed_mass_frac = {'NaCl': 2.827e-2,
                      'CaSO4': 1.298e-3,
                      'MgSO4': 1.529e-3,
                      'MgCl2': 4.251e-3,
                      'H2O': 0.9647}
    for s in feed_mass_frac:
        m.fs.splitter.inlet.flow_mass_phase_comp[0, 'Liq', s].fix(feed_flow_mass * feed_mass_frac[s])
    m.fs.splitter.inlet.pressure[0].fix(101325)
    m.fs.splitter.inlet.temperature[0].fix(298.15)
    # splitter
    m.fs.splitter.split_fraction[0, 'bypass'].fix(0.25)

    # scaling (NF and RO models are already scaled)
    calculate_scaling_factors(m.fs.splitter)
    calculate_scaling_factors(m.fs.mixer)
    calculate_scaling_factors(m.fs.tb_NF_to_RO)

    # initialize (default values are close enough)

    # unfix NF and RO inlet
    m.fs.NF.inlet.flow_mass_phase_comp.unfix()
    m.fs.NF.inlet.temperature.unfix()
    m.fs.NF.inlet.pressure.unfix()
    m.fs.RO.inlet.flow_mass_phase_comp.unfix()
    m.fs.RO.inlet.temperature.unfix()
    m.fs.RO.inlet.pressure.unfix()
    check_dof(m)

    return m

def run_flowsheet_NF_salt_basis_example():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    build_flowsheet_NF_salt_basis_example(m)
    solve_with_user_scaling(m)

    m.fs.splitter.report()
    m.fs.NF.inlet.display()
    m.fs.NF.retentate.display()
    m.fs.NF.permeate.display()
    m.fs.mixer.report()
    m.fs.RO.inlet.display()
    m.fs.RO.retentate.display()
    m.fs.RO.permeate.display()

    return m


def build_flowsheet_NF_ion_basis_example(m):
    """
    Builds a flowsheet with NF and RO. Both the NF and RO are modeled as separators with
    specified split fractions. The flowsheet includes a translator block to convert between
    the NF's multi-ion property package and the RO's seawater package.
    """
    # build flowsheet
    unit_separator.build_NF_ion_separator_example(m)
    unit_separator.build_RO_separator_example(m)
    m.fs.tb_NF_to_RO = Translator(
        default={"inlet_property_package": m.fs.NF_properties,
                 "outlet_property_package": m.fs.RO_properties})
    m.fs.splitter = Separator(default={
        "property_package": m.fs.NF_properties,
        "outlet_list": ['pretreatment', 'bypass'],
        "split_basis": SplittingType.totalFlow,
        "energy_split_basis": EnergySplittingType.equal_temperature})
    m.fs.mixer = Mixer(default={
        "property_package": m.fs.NF_properties,
        "inlet_list": ['pretreatment', 'bypass']})

    # add constraints for translator block (tb)
    m.fs.tb_NF_to_RO.eq_H2O_balance = Constraint(
        expr=m.fs.tb_NF_to_RO.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O']
        == m.fs.tb_NF_to_RO.outlet.flow_mass_phase_comp[0, 'Liq', 'H2O'])
    m.fs.tb_NF_to_RO.eq_TDS_balance = Constraint(
        expr=sum(m.fs.tb_NF_to_RO.inlet.flow_mass_phase_comp[0, 'Liq', j] for j in ['Na', 'Ca', 'Mg', 'SO4', 'Cl'])
        == m.fs.tb_NF_to_RO.outlet.flow_mass_phase_comp[0, 'Liq', 'TDS'])
    m.fs.tb_NF_to_RO.eq_equal_temperature = Constraint(
        expr=m.fs.tb_NF_to_RO.inlet.temperature[0]
        == m.fs.tb_NF_to_RO.outlet.temperature[0])
    m.fs.tb_NF_to_RO.eq_equal_pressure = Constraint(
        expr=m.fs.tb_NF_to_RO.inlet.pressure[0]
        == m.fs.tb_NF_to_RO.outlet.pressure[0])

    # connect models
    m.fs.S1 = Arc(source=m.fs.splitter.bypass, destination=m.fs.mixer.bypass)
    m.fs.S2 = Arc(source=m.fs.splitter.pretreatment, destination=m.fs.NF.inlet)
    m.fs.S3 = Arc(source=m.fs.NF.permeate, destination=m.fs.mixer.pretreatment)
    m.fs.S4 = Arc(source=m.fs.mixer.outlet, destination=m.fs.tb_NF_to_RO.inlet)
    m.fs.S5 = Arc(source=m.fs.tb_NF_to_RO.outlet, destination=m.fs.RO.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # specify (NF and RO models are already specified, mixer has 0 DOF, splitter has 1 DOF)
    # feed
    feed_flow_mass = 1
    feed_mass_frac = {'Na': 11122e-6,
                      'Ca': 382e-6,
                      'Mg': 1394e-6,
                      'SO4': 2136e-6,
                      'Cl': 20300e-6}
    m.fs.splitter.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O'].fix(
        feed_flow_mass * (1 - sum(x for x in feed_mass_frac.values())))
    for j in feed_mass_frac:
        m.fs.splitter.inlet.flow_mass_phase_comp[0, 'Liq', j].fix(feed_flow_mass * feed_mass_frac[j])
    m.fs.splitter.inlet.pressure[0].fix(101325)
    m.fs.splitter.inlet.temperature[0].fix(298.15)
    # enforce electro-neutrality for the inlet with Cl adjusting
    m.fs.splitter.inlet.flow_mass_phase_comp[0, 'Liq', 'Cl'].unfix()
    charge_dict = {'Na': 1, 'Ca': 2, 'Mg': 2, 'SO4': -2, 'Cl': -1}
    m.fs.splitter.EN_in = Constraint(
        expr=0 ==
             sum(charge_dict[j] * m.fs.splitter.mixed_state[0].flow_mol_phase_comp['Liq', j]
                 for j in charge_dict))
    # splitter
    m.fs.splitter.split_fraction[0, 'bypass'].fix(0.25)

    # scaling (NF and RO models are already scaled)
    calculate_scaling_factors(m.fs.splitter)
    calculate_scaling_factors(m.fs.mixer)
    calculate_scaling_factors(m.fs.tb_NF_to_RO)

    # initialize (default values are close enough)

    # unfix NF and RO inlet
    m.fs.NF.inlet.flow_mass_phase_comp.unfix()
    m.fs.NF.inlet.temperature.unfix()
    m.fs.NF.inlet.pressure.unfix()
    m.fs.NF.EN_in.deactivate()  # electro-neutrality already addressed in feed
    m.fs.RO.inlet.flow_mass_phase_comp.unfix()
    m.fs.RO.inlet.temperature.unfix()
    m.fs.RO.inlet.pressure.unfix()
    check_dof(m)

    return m


def run_flowsheet_NF_ion_basis_example():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    build_flowsheet_NF_ion_basis_example(m)
    solve_with_user_scaling(m)

    m.fs.splitter.report()
    m.fs.NF.inlet.display()
    m.fs.NF.retentate.display()
    m.fs.NF.permeate.display()
    m.fs.mixer.report()
    m.fs.RO.inlet.display()
    m.fs.RO.retentate.display()
    m.fs.RO.permeate.display()

    return m


if __name__ == "__main__":
    run_flowsheet_NF_salt_basis_example()
    # run_flowsheet_NF_ion_basis_example()

