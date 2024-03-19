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

from pyomo.environ import Constraint

import idaes.core.util.scaling as iscale
from idaes.models.unit_models import Mixer
from idaes.models.unit_models.mixer import MomentumMixingType, MixingType


__author__ = "Alexander V. Dudchenko"


def build_mixer(m, mixer_name, inlet_streams):
    m.fs.add_component(
        mixer_name,
        Mixer(
            property_package=m.fs.properties,
            inlet_list=inlet_streams,
            energy_mixing_type=MixingType.none,
            momentum_mixing_type=MomentumMixingType.minimize,
        ),
    )
    add_mixer_vars(m.fs.find_component(mixer_name), inlet_streams)


def add_mixer_vars(mixer, inlet_streams):
    mixer.temp_constraint = Constraint(
        expr=mixer.find_component(inlet_streams[0]).temperature[0]
        == mixer.mixed_state[0].temperature
    )
    iscale.constraint_scaling_transform(mixer.temp_constraint, 1e-2)
    inlet_streams.append("mixed")
