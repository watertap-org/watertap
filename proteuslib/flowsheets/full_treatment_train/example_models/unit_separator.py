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

"""Simple zero order separator examples"""

from pyomo.environ import ConcreteModel
from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models import Separator
from idaes.generic_models.unit_models.separator import SplittingType, EnergySplittingType
from idaes.core.util.scaling import calculate_scaling_factors
import proteuslib.property_models.seawater_prop_pack as props
from proteuslib.flowsheets.full_treatment_train.util import solve_with_user_scaling, check_dof


def build_RO_separator_example(m):
    m.fs.RO_properties = props.SeawaterParameterBlock()
    m.fs.RO = Separator(default={
        "property_package": m.fs.RO_properties,
        "outlet_list": ['retentate', 'permeate'],
        "split_basis": SplittingType.componentFlow,
        "energy_split_basis": EnergySplittingType.equal_temperature})

    # specifying
    # feed
    feed_flow_mass = 1
    feed_mass_frac_TDS = 0.035
    m.fs.RO.inlet.flow_mass_phase_comp[0, 'Liq', 'TDS'].fix(feed_flow_mass * feed_mass_frac_TDS)
    m.fs.RO.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O'].fix(feed_flow_mass * (1 - feed_mass_frac_TDS))
    m.fs.RO.inlet.pressure[0].fix(101325)
    m.fs.RO.inlet.temperature[0].fix(298.15)
    # separator
    m.fs.RO.split_fraction[0, 'permeate', 'H2O'].fix(0.5)
    m.fs.RO.split_fraction[0, 'permeate', 'TDS'].fix(0.01)
    check_dof(m)

    # scaling
    m.fs.RO_properties.set_default_scaling('flow_mass_phase_comp', 1/feed_flow_mass, index=('Liq', 'H2O'))
    m.fs.RO_properties.set_default_scaling('flow_mass_phase_comp', 1/feed_flow_mass * 1e2, index=('Liq', 'TDS'))
    calculate_scaling_factors(m.fs.RO)


def run_RO_example():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    build_RO_separator_example(m)
    solve_with_user_scaling(m)

    return m



if __name__ == "__main__":
    run_RO_example()
