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

from idaes.core.util.scaling import (calculate_scaling_factors)
from proteuslib.flowsheets.full_treatment_train.flowsheet_components import (pretreatment_softening)

# Added import statements for testing.
#       Need the pretreatment_stoich_softening_block functions to setup
#       flowsheet to solve for lime dosage
from proteuslib.flowsheets.full_treatment_train.flowsheet_components.chemistry.pretreatment_stoich_softening_block import *

def build_components(m):
    # build flowsheet
    pretrt_port = pretreatment_softening.build(m)
    property_models.build_prop(m, base='TDS')
    pretreatment_softening.build_tb(m)

    # Arc to translator block
    m.fs.s_pretrt_tb = Arc(source=pretrt_port['out'], destination=m.fs.tb_pretrt_to_desal.inlet)


def build(m):
    """
    Build a flowsheet with lime softening as the pretreatment process.
    """
    build_components(m)

    # set up costing
    financials.add_costing_param_block(m.fs)
    # annual water production
    m.fs.annual_water_production = Expression(
        expr=pyunits.convert(m.fs.tb_pretrt_to_desal.properties_out[0].flow_vol, to_units=pyunits.m ** 3 / pyunits.year)
             * m.fs.costing_param.load_factor)
    costing.build_costing(m, module=financials)

    m.fs.removal_Ca = Expression(
        expr=(m.fs.stoich_softening_mixer_unit.inlet_stream_state[0.0].flow_mol_comp['Ca(HCO3)2']
              - m.fs.stoich_softening_separator_unit.outlet_stream_state[0.0].flow_mol_comp['Ca(HCO3)2'])
             / m.fs.stoich_softening_mixer_unit.inlet_stream_state[0.0].flow_mol_comp['Ca(HCO3)2']
    )
    m.fs.removal_Mg = Expression(
        expr=(m.fs.stoich_softening_mixer_unit.inlet_stream_state[0.0].flow_mol_comp['Mg(HCO3)2']
              - m.fs.stoich_softening_separator_unit.outlet_stream_state[0.0].flow_mol_comp['Mg(HCO3)2'])
             / m.fs.stoich_softening_mixer_unit.inlet_stream_state[0.0].flow_mol_comp['Mg(HCO3)2']
    )

    return m


def scale(m):
    pretreatment_softening.scale(m)
    calculate_scaling_factors(m.fs.tb_pretrt_to_desal)


def initialize(m):
    optarg = {'nlp_scaling_method': 'user-scaling', 'halt_on_ampl_error': 'yes'}
    pretreatment_softening.initialize(m)
    propagate_state(m.fs.s_pretrt_tb)
    m.fs.tb_pretrt_to_desal.initialize(optarg=optarg)


def report(m):
    pretreatment_softening.display(m)
    m.fs.tb_pretrt_to_desal.report()


def solve_flowsheet():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    build(m)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scale
    scale(m)
    calculate_scaling_factors(m)

    # initialize
    initialize(m)

    check_dof(m)
    solve_with_user_scaling(m, tee=True, fail_flag=True)

    # report
    report(m)

    return m


def simulate(m):
    solve_with_user_scaling(m, tee=False, fail_flag=True)

if __name__ == "__main__":
    m = solve_flowsheet()

    for x in [0, 0.002, 0.004, 0.006, 0.008, 0.010, 0.012, 0.014, 0.02, 0.03, 0.04, 0.05, 0.1, 0.128]:
        m.fs.stoich_softening_mixer_unit.lime_stream_state[0].flow_mol.fix(x)
        simulate(m)
        print('Lime flow (kg/day): %.3f, Ca removal: %.2f, Mg removal: %.2f, Ca: %.2f, Mg: %.2f, LCOW: %.2f' %
              (value(m.fs.stoich_softening_mixer_unit.lime_stream_state[0].flow_mol)
               * 74.09e-3 * 60 * 60 * 24,
               value(m.fs.removal_Ca),
               value(m.fs.removal_Mg),
               value(m.fs.stoich_softening_separator_unit.outlet_stream_state[0.0].mole_frac_comp['Ca(HCO3)2'] * 1e6),
               value(m.fs.stoich_softening_separator_unit.outlet_stream_state[0.0].mole_frac_comp['Mg(HCO3)2'] * 1e6),
               value(m.fs.costing.LCOW)
               )
              )
