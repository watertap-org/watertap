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

from pyomo.network import Arc
from idaes.generic_models.unit_models import Feed
from proteuslib.flowsheets.full_treatment_train.example_flowsheets import feed_block
from proteuslib.flowsheets.full_treatment_train.example_models import unit_separator, unit_ZONF, property_models


def build_NF_no_bypass(m, NF_type='ZO', NF_base='ion'):
    """
    Builds an NF model with no bypass (i.e. just the NF model).
    Arguments:
        NF_type: 'Sep' or 'ZO', default = 'ZO'
        NF_base: 'ion' or 'salt', default = 'ion'
    """
    # feed
    feed_block.build_feed(m, base=NF_base)

    # NF
    if NF_type == 'Sep':
        unit_separator.build_SepNF(m, base=NF_base)
    elif NF_type == 'ZO':
        unit_ZONF.build_ZONF(m)
    else:
        raise ValueError('Unexpected model type {NF_type} provided to build_NF_no_bypass'
                         ''.format(NF_type=NF_type))

    # connections
    m.fs.s_feed_NF = Arc(source=m.fs.feed.outlet, destination=m.fs.NF.inlet)


def build_NF_with_bypass():
    pass
