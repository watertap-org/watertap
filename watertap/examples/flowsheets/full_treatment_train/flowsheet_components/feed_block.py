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

"""Feed blocks for supported property packages"""

from pyomo.environ import Constraint
from idaes.models.unit_models import Feed
from watertap.examples.flowsheets.full_treatment_train.model_components import (
    property_models,
)
from idaes.core.util.scaling import calculate_scaling_factors


def build_feed(m, base="TDS"):
    """
    Build a feed block for a specified property base. The state vars are fixed to the standard condition.
    Property bases include: 'TDS', 'ion', 'salt'
    """

    # feed block
    prop = property_models.get_prop(m, base=base)

    # build
    m.fs.feed = Feed(property_package=prop)

    # specify
    property_models.specify_feed(m.fs.feed.properties[0], base=base)
    m.fs.feed.properties[
        0
    ].mass_frac_phase_comp  # touch so the block can be initialized
