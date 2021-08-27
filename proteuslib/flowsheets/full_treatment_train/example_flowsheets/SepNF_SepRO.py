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

"""Flowsheets with zero order separators for NF and RO examples"""

from pyomo.environ import ConcreteModel, Constraint, TransformationFactory
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models import Separator, Mixer
from idaes.generic_models.unit_models.separator import SplittingType, EnergySplittingType
from proteuslib.flowsheets.full_treatment_train.example_models import unit_separator, feed_specification
from proteuslib.flowsheets.full_treatment_train.example_flowsheets import translator_block
from proteuslib.flowsheets.full_treatment_train.util import solve_with_user_scaling, check_dof
from idaes.core.util.scaling import calculate_scaling_factors


def build_flowsheet_NF_salt_basis_example(m):
    """
    Build a flowsheet with NF and RO. Both the NF and RO are modeled as separators with
    specified split fractions. The flowsheet includes a translator block to convert between
    the NF's multi-salt property package and the RO's seawater package.
    """
    # build flowsheet
    unit_separator.build_NF_salt_example(m)
    unit_separator.build_RO_example(m)
    translator_block.build_tb_salt_to_TDS(m, prop_salt=m.fs.NF_properties, prop_TDS=m.fs.RO_properties)
    m.fs.splitter = Separator(default={
            "property_package": m.fs.NF_properties,
            "outlet_list": ['pretreatment', 'bypass'],
            "split_basis": SplittingType.totalFlow,
            "energy_split_basis": EnergySplittingType.equal_temperature})
    m.fs.mixer = Mixer(default={
            "property_package": m.fs.NF_properties,
            "inlet_list": ['pretreatment', 'bypass']})

    # connect models
    m.fs.S1 = Arc(source=m.fs.splitter.bypass, destination=m.fs.mixer.bypass)
    m.fs.S2 = Arc(source=m.fs.splitter.pretreatment, destination=m.fs.NF.inlet)
    m.fs.S3 = Arc(source=m.fs.NF.permeate, destination=m.fs.mixer.pretreatment)
    m.fs.S4 = Arc(source=m.fs.mixer.outlet, destination=m.fs.tb_salt_to_TDS.inlet)
    m.fs.S5 = Arc(source=m.fs.tb_salt_to_TDS.outlet, destination=m.fs.RO.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # specify (NF and RO models are already specified, mixer has 0 DOF, splitter has 1 DOF)
    # feed
    feed_specification.specify_seawater_salts(m.fs.splitter.mixed_state[0])
    # splitter
    m.fs.splitter.split_fraction[0, 'bypass'].fix(0.25)
    # unfix NF and RO inlet
    m.fs.NF.inlet.flow_mass_phase_comp.unfix()
    m.fs.NF.inlet.temperature.unfix()
    m.fs.NF.inlet.pressure.unfix()
    m.fs.RO.inlet.flow_mass_phase_comp.unfix()
    m.fs.RO.inlet.temperature.unfix()
    m.fs.RO.inlet.pressure.unfix()
    check_dof(m)

    # scaling (NF and RO models and translator are already scaled)
    calculate_scaling_factors(m.fs.splitter)
    calculate_scaling_factors(m.fs.mixer)

    # initialize
    m.fs.splitter.initialize(optarg={'nlp_scaling_method': 'user-scaling'})
    m.fs.mixer.initialize(optarg={'nlp_scaling_method': 'user-scaling'})

    return m


def build_flowsheet_NF_ion_basis_example(m):
    """
    Build a flowsheet with NF and RO. Both the NF and RO are modeled as separators with
    specified split fractions. The flowsheet includes a translator block to convert between
    the NF's multi-ion property package and the RO's seawater package.
    """
    # build flowsheet
    unit_separator.build_NF_ion_example(m)
    unit_separator.build_RO_example(m)
    translator_block.build_tb_ion_to_TDS(m, prop_ion=m.fs.NF_properties, prop_TDS=m.fs.RO_properties)
    m.fs.splitter = Separator(default={
        "property_package": m.fs.NF_properties,
        "outlet_list": ['pretreatment', 'bypass'],
        "split_basis": SplittingType.totalFlow,
        "energy_split_basis": EnergySplittingType.equal_temperature})
    m.fs.mixer = Mixer(default={
        "property_package": m.fs.NF_properties,
        "inlet_list": ['pretreatment', 'bypass']})

    # connect models
    m.fs.S1 = Arc(source=m.fs.splitter.bypass, destination=m.fs.mixer.bypass)
    m.fs.S2 = Arc(source=m.fs.splitter.pretreatment, destination=m.fs.NF.inlet)
    m.fs.S3 = Arc(source=m.fs.NF.permeate, destination=m.fs.mixer.pretreatment)
    m.fs.S4 = Arc(source=m.fs.mixer.outlet, destination=m.fs.tb_ion_to_TDS.inlet)
    m.fs.S5 = Arc(source=m.fs.tb_ion_to_TDS.outlet, destination=m.fs.RO.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # specify (NF and RO models are already specified, mixer has 0 DOF, splitter has 1 DOF)
    # feed
    feed_specification.specify_seawater_ions(m.fs.splitter.mixed_state[0])
    # splitter
    m.fs.splitter.split_fraction[0, 'bypass'].fix(0.25)
    # unfix NF and RO inlet
    m.fs.NF.inlet.flow_mass_phase_comp.unfix()
    m.fs.NF.inlet.temperature.unfix()
    m.fs.NF.inlet.pressure.unfix()
    m.fs.RO.inlet.flow_mass_phase_comp.unfix()
    m.fs.RO.inlet.temperature.unfix()
    m.fs.RO.inlet.pressure.unfix()
    check_dof(m)

    # scaling (NF and RO models and translator are already scaled)
    calculate_scaling_factors(m.fs.splitter)
    calculate_scaling_factors(m.fs.mixer)

    # initialize
    m.fs.splitter.initialize(optarg={'nlp_scaling_method': 'user-scaling'})
    m.fs.mixer.initialize(optarg={'nlp_scaling_method': 'user-scaling'})

    return m


def run_flowsheet_example(case):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    if case == 'salt':
        build_flowsheet_NF_salt_basis_example(m)
    elif case == 'ion':
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
    run_flowsheet_example('salt')
    run_flowsheet_example('ion')
