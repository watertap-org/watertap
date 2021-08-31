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

"""Desalination flowsheet components"""

from pyomo.environ import ConcreteModel, TransformationFactory
from idaes.core import FlowsheetBlock
from pyomo.network import Arc
from idaes.core.util.scaling import calculate_scaling_factors, set_scaling_factor
from idaes.core.util.initialization import propagate_state
from proteuslib.unit_models.pump_isothermal import Pump
from proteuslib.flowsheets.full_treatment_train.example_flowsheets import feed_block
from proteuslib.flowsheets.full_treatment_train.example_models import unit_separator, unit_0DRO, property_models
from proteuslib.flowsheets.full_treatment_train.util import solve_with_user_scaling, check_dof


def build_desalination_RO(m, has_feed=False, RO_type='0D', RO_base='TDS', RO_level='simple'):
    """
    Builds RO desalination including specified feed and auxiliary equipment.
    The built components are initialized based on default values.
    Arguments:
        has_feed: True or False, default = False
        RO_type: 'Sep' or '0D', default = '0D'
        RO_level: 'simple' or 'detailed', default = 'simple'
        RO_base: 'TDS' only, default = 'ion'
    """
    desal_port = {}
    optarg = {'nlp_scaling_method': 'user-scaling'}
    prop = property_models.get_prop(m, base=RO_base)

    if has_feed:
        # build feed
        feed_block.build_feed(m, base=RO_base)

    # build RO
    if RO_type == 'Sep':
        unit_separator.build_SepRO(m, base=RO_base)
    elif RO_type == '0D':
        unit_0DRO.build_RO(m, base='TDS', level=RO_level)
    else:
        raise ValueError('Unexpected model type {RO_type} provided to build_desalination_RO'
                         ''.format(RO_type=RO_type))

    # auxiliary units
    if RO_type == 'Sep':
        # build auxiliary units (none)

        # connect models
        if has_feed:
            m.fs.s_desal_feed_RO = Arc(source=m.fs.feed.outlet, destination=m.fs.RO.inlet)

        # specify (RO and feed already specified)

        # scaling (RO and feed already scaled)

        # initialize
        if has_feed:
            m.fs.feed.initialize(optarg=optarg)
            propagate_state(m.fs.s_desal_feed_RO)
        # m.fs.RO.mixed_state[0].mass_frac_phase_comp  # touch properties to have a constraint on stateblock
        # m.fs.RO.permeate_state[0].mass_frac_phase_comp
        # m.fs.RO.retentate_state[0].mass_frac_phase_comp
        # m.fs.RO.initialize(optarg=optarg)  # IDAES error on initializing separators, simple enough to not need it

        # inlet/outlet ports for pretreatment
        if not has_feed:
            desal_port['in'] = m.fs.RO.inlet

    elif RO_type == '0D':
        # build auxiliary units
        m.fs.pump_RO = Pump(default={'property_package': prop})

        # connect models
        if has_feed:
            m.fs.s_desal_feed_pumpRO = Arc(source=m.fs.feed.outlet, destination=m.fs.pump_RO.inlet)
        m.fs.s_desal_pumpRO_RO = Arc(source=m.fs.pump_RO.outlet, destination=m.fs.RO.inlet)

        # specify (RO already specified, Pump 2 DOF)
        m.fs.pump_RO.efficiency_pump.fix(0.80)
        m.fs.pump_RO.control_volume.properties_out[0].pressure.fix(50e5)

        # scaling (RO already scaled)
        set_scaling_factor(m.fs.pump_RO.control_volume.work, 1e-3)
        calculate_scaling_factors(m.fs.pump_RO)

        # initialize
        if has_feed:
            m.fs.feed.initialize(optarg=optarg)
            propagate_state(m.fs.s_desal_feed_pumpRO)
        m.fs.pump_RO.initialize(optarg=optarg)
        propagate_state(m.fs.s_desal_pumpRO_RO)
        m.fs.RO.initialize(optarg=optarg)

        # inlet/outlet ports for pretreatment
        if not has_feed:
            desal_port['in'] = m.fs.pump_RO.inlet

    desal_port['out'] = m.fs.RO.permeate
    desal_port['waste'] = m.fs.RO.retentate

    return desal_port


if __name__ == "__main__":
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    property_models.build_prop(m, base='TDS')
    build_desalination_RO(m, has_feed=False, RO_type='Sep', RO_base='TDS', RO_level='simple')
    # property_models.specify_feed(m.fs.pump_RO.control_volume.properties_in[0], base='TDS')
    property_models.specify_feed(m.fs.RO.mixed_state[0], base='TDS')

    TransformationFactory("network.expand_arcs").apply_to(m)
    check_dof(m)
    solve_with_user_scaling(m)

    # m.fs.pump_RO.report()
    m.fs.RO.report()
