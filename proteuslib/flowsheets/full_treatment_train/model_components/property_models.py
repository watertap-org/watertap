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

"""Seawater feed specifications for supported property packages"""

from pyomo.environ import ConcreteModel
from idaes.core import FlowsheetBlock
from idaes.generic_models.properties.core.generic.generic_property import GenericParameterBlock
from idaes.core.util.scaling import calculate_scaling_factors
from proteuslib.property_models import seawater_prop_pack
from proteuslib.flowsheets.full_treatment_train.model_components import seawater_salt_prop_pack, seawater_ion_prop_pack
from proteuslib.flowsheets.full_treatment_train.model_components.eNRTL import entrl_config_FpcTP
from proteuslib.flowsheets.full_treatment_train.util import solve_with_user_scaling


def build_prop(m, base='TDS'):
    """
    Builds a property package for the specified base. Includes default scaling.
    Bases include: 'TDS', 'ion', 'salt'.
    """
    if base == 'TDS':
        m.fs.prop_TDS = seawater_prop_pack.SeawaterParameterBlock()

        m.fs.prop_TDS.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
        m.fs.prop_TDS.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'TDS'))

    elif base == 'ion':
        m.fs.prop_ion = seawater_ion_prop_pack.PropParameterBlock()

        m.fs.prop_ion.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
        m.fs.prop_ion.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'Na'))
        m.fs.prop_ion.set_default_scaling('flow_mass_phase_comp', 1e4, index=('Liq', 'Ca'))
        m.fs.prop_ion.set_default_scaling('flow_mass_phase_comp', 1e3, index=('Liq', 'Mg'))
        m.fs.prop_ion.set_default_scaling('flow_mass_phase_comp', 1e3, index=('Liq', 'SO4'))
        m.fs.prop_ion.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'Cl'))

    elif base == 'salt':
        m.fs.prop_salt = seawater_salt_prop_pack.PropParameterBlock()

        m.fs.prop_salt.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
        m.fs.prop_salt.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'NaCl'))
        m.fs.prop_salt.set_default_scaling('flow_mass_phase_comp', 1e3, index=('Liq', 'CaSO4'))
        m.fs.prop_salt.set_default_scaling('flow_mass_phase_comp', 1e3, index=('Liq', 'MgSO4'))
        m.fs.prop_salt.set_default_scaling('flow_mass_phase_comp', 1e3, index=('Liq', 'MgCl2'))

    elif base == 'eNRTL':
        m.fs.prop_eNRTL = GenericParameterBlock(default=entrl_config_FpcTP.configuration)

        # default scaling in config file

    else:
        raise ValueError('Unexpected property base {base} provided to build_prop'
                         ''.format(base=base))


def get_prop(m, base='TDS'):
    if base == 'TDS':
        prop = m.fs.prop_TDS
    elif base == 'ion':
        prop = m.fs.prop_ion
    elif base == 'salt':
        prop = m.fs.prop_salt
    elif base == 'eNRTL':
        prop = m.fs.prop_eNRTL
    else:
        raise ValueError('Unexpected property base {base} for get_prop'
                         ''.format(base=base))
    return prop



def specify_feed(sb, base='TDS'):
    """
    Fixes the state variables on the stateblock to the base seawater composition for
    the specified base. Bases include: 'TDS', 'ion', 'salt'.
    """
    sb.pressure.fix(101325)
    sb.temperature.fix(298.15)

    feed_flow_mass = 1
    if base == 'TDS':
        feed_mass_frac_TDS = 0.035
        sb.flow_mass_phase_comp['Liq', 'TDS'].fix(feed_flow_mass * feed_mass_frac_TDS)
        sb.flow_mass_phase_comp['Liq', 'H2O'].fix(feed_flow_mass * (1 - feed_mass_frac_TDS))

    elif base == 'ion':
        feed_mass_frac = {'Na': 11122e-6,
                          'Ca': 382e-6,
                          'Mg': 1394e-6,
                          'SO4': 2136e-6,
                          'Cl': 20316.88e-6}
        sb.flow_mass_phase_comp['Liq', 'H2O'].fix(
            feed_flow_mass * (1 - sum(x for x in feed_mass_frac.values())))
        for j in feed_mass_frac:
            sb.flow_mass_phase_comp['Liq', j].fix(feed_flow_mass * feed_mass_frac[j])

    elif base == 'salt':
        feed_mass_frac = {'NaCl': 2.827e-2,
                          'CaSO4': 1.298e-3,
                          'MgSO4': 1.529e-3,
                          'MgCl2': 4.251e-3,
                          'H2O': 0.9647}
        for s in feed_mass_frac:
            sb.flow_mass_phase_comp['Liq', s].fix(feed_flow_mass * feed_mass_frac[s])

    else:
        raise ValueError('Unexpected property base {base} provided to specify_feed'
                         ''.format(base=base))


def solve_specify_feed(base):
    # build state block
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    build_prop(m, base=base)
    if base == 'TDS':
        m.fs.stream = m.fs.prop_TDS.build_state_block([0], default={})
    elif base == 'ion':
        m.fs.stream = m.fs.prop_ion.build_state_block([0], default={})
    elif base == 'salt':
        m.fs.stream = m.fs.prop_salt.build_state_block([0], default={})
    specify_feed(m.fs.stream[0], base=base)

    m.fs.stream[0].mass_frac_phase_comp  # touch a variable to have a model with at least one constraint

    # scale
    calculate_scaling_factors(m.fs)
    # solve
    solve_with_user_scaling(m)
    # display
    m.fs.stream.display()

    return m


if __name__ == "__main__":
    solve_specify_feed('TDS')
    solve_specify_feed('ion')
    solve_specify_feed('salt')
