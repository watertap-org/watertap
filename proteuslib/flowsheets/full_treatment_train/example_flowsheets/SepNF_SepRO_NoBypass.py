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

"""Simple flowsheets with zero order separators for NF and RO examples"""

from pyomo.environ import ConcreteModel, Constraint, TransformationFactory
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models.translator import Translator
from proteuslib.flowsheets.full_treatment_train.example_models import unit_separator
from proteuslib.flowsheets.full_treatment_train.example_flowsheets import translator_block
from proteuslib.flowsheets.full_treatment_train.util import solve_with_user_scaling, check_dof
from idaes.core.util.scaling import calculate_scaling_factors


def build_flowsheet_NF_salt_basis_example(m):
    """
    Build a flowsheet with NF and RO connected without bypass. Both the NF and RO are
    modeled as separators with specified split fractions. The flowsheet includes a translator
    block to convert between the NF's multi-salt property package and the RO's seawater package.
    """
    # build flowsheet
    unit_separator.build_NF_salt_example(m)
    unit_separator.build_RO_example(m)
    translator_block.build_tb_salt_to_TDS(m, prop_salt=m.fs.NF_properties, prop_TDS=m.fs.RO_properties)

    # connect models
    m.fs.S1 = Arc(source=m.fs.NF.permeate, destination=m.fs.tb_salt_to_TDS.inlet)
    m.fs.S2 = Arc(source=m.fs.tb_salt_to_TDS.outlet, destination=m.fs.RO.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # specify (unit models are already specified)

    # scaling (unit models and translator blocks are already scaled)

    # initialize (default values are close enough)

    # unfix RO inlet
    m.fs.RO.inlet.flow_mass_phase_comp.unfix()
    m.fs.RO.inlet.temperature.unfix()
    m.fs.RO.inlet.pressure.unfix()
    check_dof(m)

    return m


def build_flowsheet_NF_ion_basis_example(m):
    """
    Build a flowsheet with NF and RO connected without bypass. Both the NF and RO are
    modeled as separators with specified split fractions. The flowsheet includes a translator
    block to convert between the NF's multi-ion property package and the RO's seawater package.
    """
    # build flowsheet
    unit_separator.build_NF_ion_example(m)
    unit_separator.build_RO_example(m)
    translator_block.build_tb_ion_to_TDS(m, prop_ion=m.fs.NF_properties, prop_TDS=m.fs.RO_properties)

    # connect models
    m.fs.S1 = Arc(source=m.fs.NF.permeate, destination=m.fs.tb_ion_to_TDS.inlet)
    m.fs.S2 = Arc(source=m.fs.tb_ion_to_TDS.outlet, destination=m.fs.RO.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # specify (unit models are already specified)

    # scaling (unit models and translator blocks are already scaled)

    # initialize (default values are close enough)

    # unfix RO inlet
    m.fs.RO.inlet.flow_mass_phase_comp.unfix()
    m.fs.RO.inlet.temperature.unfix()
    m.fs.RO.inlet.pressure.unfix()
    check_dof(m)

    return m


def run_flowsheet_example(case):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    if case == 'salt':
        build_flowsheet_NF_salt_basis_example(m)
    elif case == 'ion':
        build_flowsheet_NF_ion_basis_example(m)

    solve_with_user_scaling(m)

    m.fs.NF.inlet.display()
    m.fs.NF.retentate.display()
    m.fs.NF.permeate.display()
    m.fs.RO.inlet.display()
    m.fs.RO.retentate.display()
    m.fs.RO.permeate.display()

    return m


if __name__ == "__main__":
    run_flowsheet_example('salt')
    run_flowsheet_example('ion')

