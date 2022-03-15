###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################

from pyomo.environ import ConcreteModel
from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models import Feed
from idaes.core.util.scaling import calculate_scaling_factors
from watertap.examples.flowsheets.high_pressure_RO.components import ion_prop_pack


def build_prop(m, case='seawater'):
    """
    Builds a property package for the specified case.
    cases include: 'seawater'
    """
    if case == 'seawater':
        ion_list = ['Na_+', 'K_+', 'Mg_2+', 'Ca_2+', 'Cl_-', 'SO4_2-', 'HCO3_-']
        m.fs.prop_feed = ion_prop_pack.PropParameterBlock(
            default={'ion_list': ion_list})
    else:
        raise ValueError('Unexpected feed case {case} provided to build_prop'
                         ''.format(case=case))


def build_specified_feed(m, case='seawater'):
    """
    Build a specified feed block for the given case. The state vars are fixed to the standard condition.
    Feed cases include: 'seawater',
    """

    # build property block
    build_prop(m, case=case)

    # build
    m.fs.feed = Feed(default={'property_package': m.fs.prop_feed})

    # specify
    specify_feed(m.fs.feed.properties[0], case=case)


def specify_feed(sb, case='seawater'):
    """
    Fixes the state variables on the stateblock to the feed case.
    Feed cases include: 'seawater'
    """
    sb.pressure.fix(101325)
    sb.temperature.fix(298.15)
    sb.flow_vol.fix(1e-3)

    if case == 'seawater':
        conc_dict = {'Na_+': 10.556,
                     'K_+': 0.380,
                     'Ca_2+': 0.400,
                     'Mg_2+': 1.262,
                     'Cl_-': 18.980,
                     'SO4_2-': 2.649,
                     'HCO3_-': 0.129}
        for (j, val) in conc_dict.items():
            sb.conc_mass_comp[j].fix(val)
    else:
        raise ValueError('Unexpected feed case {case} provided to specify_feed'
                         ''.format(case=case))
