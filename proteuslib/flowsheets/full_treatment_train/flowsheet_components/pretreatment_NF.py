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

"""Pretreatment flowsheet components"""

from pyomo.environ import ConcreteModel, TransformationFactory
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models import Feed
from proteuslib.flowsheets.full_treatment_train.model_components import Separator, Mixer
from idaes.generic_models.unit_models.separator import SplittingType, EnergySplittingType
from idaes.core.util.scaling import (calculate_scaling_factors,
                                     set_scaling_factor,
                                     get_scaling_factor,
                                     constraint_scaling_transform)
from idaes.core.util.initialization import propagate_state
from proteuslib.unit_models.pump_isothermal import Pump
from proteuslib.flowsheets.full_treatment_train.flowsheet_components import feed_block
from proteuslib.flowsheets.full_treatment_train.model_components import unit_separator, unit_ZONF, property_models
from proteuslib.flowsheets.full_treatment_train.util import solve_with_user_scaling, check_dof


def build_pretreatment_NF(m, has_bypass=True, NF_type='ZO', NF_base='ion'):
    """
    Builds NF pretreatment including specified feed and auxiliary equipment.
    Arguments:
        has_bypass: True or False, default = True
        NF_type: 'Sep' or 'ZO', default = 'ZO'
        NF_base: 'ion' or 'salt', default = 'ion'
    """
    pretrt_port = {}
    prop = property_models.get_prop(m, base=NF_base)

    # build feed
    feed_block.build_feed(m, base=NF_base)

    # build NF
    if NF_type == 'Sep':
        unit_separator.build_SepNF(m, base=NF_base)
    elif NF_type == 'ZO':
        unit_ZONF.build_ZONF(m, base=NF_base)
        m.fs.pump_NF = Pump(default={'property_package': prop})
    else:
        raise ValueError('Unexpected model type {NF_type} provided to build_NF_no_bypass'
                         ''.format(NF_type=NF_type))

    if has_bypass:
        # build auxiliary units
        m.fs.splitter = Separator(default={
            "property_package": prop,
            "outlet_list": ['pretreatment', 'bypass'],
            "split_basis": SplittingType.totalFlow,
            "energy_split_basis": EnergySplittingType.equal_temperature})
        m.fs.mixer = Mixer(default={
            "property_package": prop,
            "inlet_list": ['pretreatment', 'bypass']})

        # connect models
        m.fs.s_pretrt_feed_splitter = Arc(source=m.fs.feed.outlet, destination=m.fs.splitter.inlet)
        m.fs.s_pretrt_splitter_mixer = Arc(source=m.fs.splitter.bypass, destination=m.fs.mixer.bypass)
        if NF_type == 'ZO':
            m.fs.s_pretrt_splitter_pumpNF = Arc(source=m.fs.splitter.pretreatment, destination=m.fs.pump_NF.inlet)
            m.fs.s_pretrt_pumpNF_NF = Arc(source=m.fs.pump_NF.outlet, destination=m.fs.NF.inlet)
        else:
            m.fs.s_pretrt_splitter_NF = Arc(source=m.fs.splitter.pretreatment, destination=m.fs.NF.inlet)
        m.fs.s_pretrt_NF_mixer = Arc(source=m.fs.NF.permeate, destination=m.fs.mixer.pretreatment)

        # specify (NF and feed is already specified, mixer has 0 DOF, splitter has 1 DOF, NF pump has 2 DOF)
        # splitter
        m.fs.splitter.split_fraction[0, 'bypass'].fix(0.1)
        if NF_type == 'ZO':
            m.fs.pump_NF.efficiency_pump.fix(0.80)
            m.fs.pump_NF.control_volume.properties_out[0].pressure.fix(4e5)

        # inlet/outlet ports for pretreatment
        pretrt_port['out'] = m.fs.mixer.outlet
        pretrt_port['waste'] = m.fs.NF.retentate

    else:  # no bypass
        # build auxiliary units (none)

        # connect models
        if NF_type == 'ZO':
            m.fs.s_pretrt_feed_pumpNF = Arc(source=m.fs.feed.outlet, destination=m.fs.pump_NF.inlet)
            m.fs.s_pretrt_pumpNF_NF = Arc(source=m.fs.pump_NF.outlet, destination=m.fs.NF.inlet)
            #TODO: should source be m.fs.pump_NF.outlet? Double-check here and other arcs with pump_NF
        else:
            m.fs.s_pretrt_feed_NF = Arc(source=m.fs.feed.outlet, destination=m.fs.NF.inlet)

        # specify (NF and feed are already specified, NF pump has 2 DOF)
        if NF_type == 'ZO':
            m.fs.pump_NF.efficiency_pump.fix(0.80)
            m.fs.pump_NF.control_volume.properties_out[0].pressure.fix(4e5)

        # inlet/outlet ports for pretreatment
        pretrt_port['out'] = m.fs.NF.permeate
        pretrt_port['waste'] = m.fs.NF.retentate

    return pretrt_port


def scale_pretreatment_NF(m, **kwargs):
    calculate_scaling_factors(m.fs.feed)
    calculate_scaling_factors(m.fs.NF)

    if kwargs['has_bypass']:
        calculate_scaling_factors(m.fs.splitter)
        set_scaling_factor(m.fs.splitter.split_fraction, 1)  # TODO: should have an IDAES default
        constraint_scaling_transform(m.fs.splitter.sum_split_frac[0], 1)  # TODO: should have an IDAES default
        calculate_scaling_factors(m.fs.mixer)
        set_scaling_factor(m.fs.mixer.minimum_pressure,
                           get_scaling_factor(m.fs.mixer.mixed_state[0].pressure)
                           )  # TODO: IDAES should have a default and link to the constraint
        for c in [m.fs.mixer.minimum_pressure_constraint[0, 1],
                  m.fs.mixer.minimum_pressure_constraint[0, 2],
                  m.fs.mixer.mixture_pressure[0.0]]:
            constraint_scaling_transform(c, get_scaling_factor(m.fs.mixer.minimum_pressure))

    if kwargs['NF_type'] == 'ZO':
        set_scaling_factor(m.fs.pump_NF.control_volume.work, 1e-3)
        calculate_scaling_factors(m.fs.pump_NF)
        set_scaling_factor(m.fs.pump_NF.ratioP, 1)  # TODO: IDAES should have a default and link to the constraint

def initialize_pretreatment_NF(m, **kwargs):
    optarg = {'nlp_scaling_method': 'user-scaling'}

    if kwargs['has_bypass']:
        m.fs.feed.initialize(optarg=optarg)
        propagate_state(m.fs.s_pretrt_feed_splitter)
        m.fs.splitter.initialize(optarg=optarg)
        propagate_state(m.fs.s_pretrt_splitter_mixer)
        if kwargs['NF_type'] == 'ZO':
            propagate_state(m.fs.s_pretrt_splitter_pumpNF)
            m.fs.pump_NF.initialize(optarg=optarg)
            propagate_state(m.fs.s_pretrt_pumpNF_NF)
            m.fs.NF.initialize(optarg=optarg)
        else:  # NF_type == 'Sep'
            propagate_state(m.fs.s_pretrt_splitter_NF)
            # m.fs.NF.initialize(optarg=optarg)  # IDAES error when NF is a separator TODO: address in IDAES
        propagate_state(m.fs.s_pretrt_NF_mixer)
        m.fs.mixer.initialize(optarg=optarg)

    else:  # no bypass
        m.fs.feed.initialize(optarg=optarg)
        if kwargs['NF_type'] == 'ZO':
            propagate_state(m.fs.s_pretrt_feed_pumpNF)
            m.fs.pump_NF.initialize(optarg=optarg)
            propagate_state(m.fs.s_pretrt_pumpNF_NF)
            m.fs.NF.initialize(optarg=optarg)
        else:  # NF_type == 'Sep'
            propagate_state(m.fs.s_pretrt_feed_NF)
            # m.fs.NF.initialize(optarg=optarg)  # IDAES error when NF is a separator TODO: address in IDAES


def display_pretreatment_NF(m, **kwargs):

    m.fs.feed.report()

    if kwargs['has_bypass']:
        m.fs.splitter.report()
        if kwargs['NF_type'] == 'ZO':
            m.fs.pump_NF.report()
        m.fs.NF.report()
        m.fs.mixer.report()
    else:  # no bypass
        if kwargs['NF_type'] == 'ZO':
            m.fs.pump_NF.report()
        m.fs.NF.report()


def solve_pretreatment_NF(**kwargs):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    property_models.build_prop(m, base=kwargs['NF_base'])
    build_pretreatment_NF(m, **kwargs)
    TransformationFactory("network.expand_arcs").apply_to(m)

    scale_pretreatment_NF(m, **kwargs)

    initialize_pretreatment_NF(m, **kwargs)

    check_dof(m)
    solve_with_user_scaling(m, tee=True, fail_flag=True)

    display_pretreatment_NF(m, **kwargs)


if __name__ == "__main__":
    solve_pretreatment_NF(has_bypass=True, NF_type='ZO', NF_base='ion')
