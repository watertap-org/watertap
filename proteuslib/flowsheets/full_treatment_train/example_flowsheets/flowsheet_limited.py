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

"""Flowsheet examples that are limited (i.e. do not satisfy minimum viable product requirements)"""

from pyomo.environ import ConcreteModel, TransformationFactory
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from proteuslib.flowsheets.full_treatment_train.example_flowsheets import pretreatment
from proteuslib.flowsheets.full_treatment_train.example_models import unit_separator, property_models
from proteuslib.flowsheets.full_treatment_train.example_flowsheets import feed_block, translator_block
from proteuslib.flowsheets.full_treatment_train.util import solve_with_user_scaling, check_dof


def build_flowsheet_NF_no_bypass(m, NF_type='ZO', NF_base='ion'):
    """
    Build a flowsheet with NF and RO connected without bypass.
    """
    # build flowsheet
    property_models.build_prop(m, base=NF_base)
    pretreatment.build_NF_no_bypass(m, NF_type=NF_type, NF_base=NF_base)

    property_models.build_prop(m, base='TDS')
    unit_separator.build_SepRO(m)

    translator_block.build_tb(m, base_inlet=NF_base, base_outlet='TDS', name_str='tb_pretrt_to_desal')

    # set up Arcs
    m.fs.s_pretrt_tb = Arc(source=m.fs.NF.permeate, destination=m.fs.tb_pretrt_to_desal.inlet)
    m.fs.s_tb_desal = Arc(source=m.fs.tb_pretrt_to_desal.outlet, destination=m.fs.RO.inlet)

    # specify (unit model parameters are already specified)

    # scaling (unit models and translator blocks are already scaled)

    # initialize (default values are close enough for simple separator models)

    return m


def solve_build_flowsheet_NF_no_bypass(NF_type='ZO', NF_base='ion'):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    build_flowsheet_NF_no_bypass(m, NF_type=NF_type, NF_base=NF_base)

    TransformationFactory("network.expand_arcs").apply_to(m)
    check_dof(m)
    solve_with_user_scaling(m)

    m.fs.NF.inlet.display()
    m.fs.NF.retentate.display()
    m.fs.NF.permeate.display()
    m.fs.RO.inlet.display()
    m.fs.RO.retentate.display()
    m.fs.RO.permeate.display()

    return m


if __name__ == "__main__":
    solve_build_flowsheet_NF_no_bypass(NF_type='ZO', NF_base='ion')
    solve_build_flowsheet_NF_no_bypass(NF_type='Sep', NF_base='salt')
    solve_build_flowsheet_NF_no_bypass(NF_type='Sep', NF_base='ion')
