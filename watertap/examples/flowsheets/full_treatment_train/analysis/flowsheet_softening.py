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

from idaes.core.util.scaling import calculate_scaling_factors
from watertap.examples.flowsheets.full_treatment_train.flowsheet_components import (
    pretreatment_softening,
    costing,
)

# Added import statements for testing.
#       Need the pretreatment_stoich_softening_block functions to setup
#       flowsheet to solve for lime dosage
from watertap.examples.flowsheets.full_treatment_train.flowsheet_components.chemistry.pretreatment_stoich_softening_block import *


def build_components(m):
    # build flowsheet
    pretrt_port = pretreatment_softening.build(m)
    property_models.build_prop(m, base="TDS")
    pretreatment_softening.build_tb(m)

    # Arc to translator block
    m.fs.s_pretrt_tb = Arc(
        source=pretrt_port["out"], destination=m.fs.tb_pretrt_to_desal.inlet
    )


def build(m):
    """
    Build a flowsheet with lime softening as the pretreatment process.
    """
    build_components(m)

    m.fs.treated_flow_vol = Expression(
        expr=m.fs.tb_pretrt_to_desal.properties_out[0].flow_vol
    )
    costing.build_costing(m)

    m.fs.removal_Ca = Expression(
        expr=(
            m.fs.stoich_softening_mixer_unit.inlet_stream_state[0.0].flow_mol_comp[
                "Ca(HCO3)2"
            ]
            - m.fs.stoich_softening_separator_unit.outlet_stream_state[
                0.0
            ].flow_mol_comp["Ca(HCO3)2"]
        )
        / m.fs.stoich_softening_mixer_unit.inlet_stream_state[0.0].flow_mol_comp[
            "Ca(HCO3)2"
        ]
    )
    m.fs.removal_Mg = Expression(
        expr=(
            m.fs.stoich_softening_mixer_unit.inlet_stream_state[0.0].flow_mol_comp[
                "Mg(HCO3)2"
            ]
            - m.fs.stoich_softening_separator_unit.outlet_stream_state[
                0.0
            ].flow_mol_comp["Mg(HCO3)2"]
        )
        / m.fs.stoich_softening_mixer_unit.inlet_stream_state[0.0].flow_mol_comp[
            "Mg(HCO3)2"
        ]
    )

    return m


def scale(m):
    pretreatment_softening.scale(m)
    calculate_scaling_factors(m.fs.tb_pretrt_to_desal)


def initialize(m):
    pretreatment_softening.initialize(m)
    propagate_state(m.fs.s_pretrt_tb)
    m.fs.tb_pretrt_to_desal.initialize()


def report(m):
    pretreatment_softening.display(m)
    m.fs.tb_pretrt_to_desal.report()


def solve_flowsheet():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    build(m)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scale
    scale(m)
    calculate_scaling_factors(m)

    # initialize
    initialize(m)

    check_dof(m)
    solve_block(m, tee=True, fail_flag=True)

    # report
    report(m)

    return m


def simulate(m, check_termination=True):
    return solve_block(m, tee=False, fail_flag=check_termination)


if __name__ == "__main__":
    m = solve_flowsheet()
