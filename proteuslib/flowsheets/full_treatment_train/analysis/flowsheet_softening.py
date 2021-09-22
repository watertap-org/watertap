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

from pyomo.environ import (ConcreteModel, Objective, Expression, Constraint, Param,
        TransformationFactory, value, units as pyunits)
from pyomo.network import Arc
from pyomo.util import infeasible
from idaes.core import FlowsheetBlock
from idaes.core.util.scaling import (calculate_scaling_factors,
                                     unscaled_constraints_generator,
                                     unscaled_variables_generator,
                                     badly_scaled_var_generator,
                                     constraint_autoscale_large_jac)
from idaes.core.util.initialization import propagate_state
from proteuslib.flowsheets.full_treatment_train.example_flowsheets import (pretreatment_softening,
                                                                           desalination,
                                                                           translator_block,
                                                                           costing,
                                                                           financials)
from proteuslib.flowsheets.full_treatment_train.example_models import property_models
from proteuslib.flowsheets.full_treatment_train.util import solve_with_user_scaling, check_dof

# Added import statements for testing.
#       Need the pretreatment_stoich_softening_block functions to setup
#       flowsheet to solve for lime dosage
from idaes.core.util import scaling as iscale
from proteuslib.flowsheets.full_treatment_train.chemistry_flowsheets.pretreatment_stoich_softening_block import *


def build(m):
    """
    Build a flowsheet with NF pretreatment and RO.
    """

    # build flowsheet
    pretrt_port = pretreatment_softening.build(m)

    property_models.build_prop(m, base='TDS')

    pretreatment_softening.build_tb(m)

    # Arc to translator block
    m.fs.s_pretrt_tb = Arc(source=pretrt_port['out'], destination=m.fs.tb_pretrt_to_desal.inlet)

    # set up costing
    financials.add_costing_param_block(m.fs)
    # annual water production
    m.fs.annual_water_production = Expression(
        expr=pyunits.convert(m.fs.tb_pretrt_to_desal.properties_out[0].flow_vol, to_units=pyunits.m ** 3 / pyunits.year)
             * m.fs.costing_param.load_factor)
    costing.build_costing(m, module=financials)

    # create useful expressions


    return m


def solve():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    build(m)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scale
    pretreatment_softening.scale(m)
    calculate_scaling_factors(m.fs.tb_pretrt_to_desal)
    calculate_scaling_factors(m)

    # initialize
    optarg = {'nlp_scaling_method': 'user-scaling', 'halt_on_ampl_error': 'yes'}
    pretreatment_softening.initialize(m)
    propagate_state(m.fs.s_pretrt_tb)
    m.fs.tb_pretrt_to_desal.initialize(optarg=optarg)

    check_dof(m)
    solve_with_user_scaling(m, tee=True, fail_flag=True, bound_push=1e-10, mu_init=1e-6)

    pretreatment_softening.display(m)
    m.fs.tb_pretrt_to_desal.report()

    return m


def simulate(m):
    solve_with_user_scaling(m, tee=False, fail_flag=True, bound_push=1e-10, mu_init=1e-6)


if __name__ == "__main__":
    m = solve()

    for x in [0, 0.002, 0.004, 0.006, 0.008, 0.010, 0.012, 0.014, 0.02, 0.03, 0.04, 0.05, 0.1]:
        m.fs.stoich_softening_mixer_unit.lime_stream_state[0].flow_mol.fix(x)
        simulate(m)
        print('Lime flow: %.3f, Hardness: %.0f, Ca: %.2f, Mg: %.2f, LCOW: %.2f' %
              (value(m.fs.stoich_softening_mixer_unit.lime_stream_state[0].flow_mol),
               value(m.fs.stoich_softening_separator_unit.hardness),
               value(m.fs.stoich_softening_separator_unit.outlet_stream_state[0.0].mole_frac_comp['Ca(HCO3)2']*1e6),
               value(m.fs.stoich_softening_separator_unit.outlet_stream_state[0.0].mole_frac_comp['Mg(HCO3)2']*1e6),
               value(m.fs.costing.LCOW))
              )
