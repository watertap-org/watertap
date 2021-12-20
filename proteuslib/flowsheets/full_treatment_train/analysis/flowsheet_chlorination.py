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
from pyomo.environ import (ConcreteModel, Objective, Expression, Constraint, Param, Block,
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
from proteuslib.flowsheets.full_treatment_train.flowsheet_components import (pretreatment_NF,
                                                                             desalination,
                                                                             gypsum_saturation_index,
                                                                             translator_block,
                                                                             costing,
                                                                             financials)
from proteuslib.flowsheets.full_treatment_train.model_components import property_models
from proteuslib.flowsheets.full_treatment_train.flowsheet_components.chemistry import (
    posttreatment_ideal_naocl_chlorination_block as chlorination)
from proteuslib.flowsheets.full_treatment_train.util import solve_with_user_scaling, check_dof

def build_components(m):
    # build flowsheet
    # property_models.build_prop(m, base='TDS')
    m.chlor = Block()
    m.chlor.fs = FlowsheetBlock(default={"dynamic": False})
    chlorination.build_ideal_naocl_chlorination_block(m.chlor)
    # chlorination.build_translator_from_RO_to_chlorination_block(m)

    # arcs

    # specify
    chlorination.set_ideal_naocl_mixer_inlets(m.chlor, dosing_rate_of_NaOCl_mg_per_s = 0.4,
                                            inlet_water_density_kg_per_L = 1,
                                            inlet_temperature_K = 298,
                                            inlet_pressure_Pa = 101325,
                                            inlet_flow_mol_per_s = 25)
    chlorination.set_ideal_naocl_chlorination_inlets(m.chlor)
    chlorination.fix_ideal_naocl_mixer_inlets(m.chlor)

    # expressions
    unit = m.chlor.fs.ideal_naocl_chlorination_unit
    m.chlor.fs.OCl_distribution = Expression(
        expr=unit.control_volume.properties_out[0.0].conc_mol_phase_comp["Liq", "OCl_-"] / 1000
        / (unit.free_chlorine/70900)
    )
    m.chlor.fs.HOCl_distribution = Expression(
        expr=unit.control_volume.properties_out[0.0].conc_mol_phase_comp["Liq", "HOCl"] / 1000
        / (unit.free_chlorine/70900)
    )
    m.chlor.fs.NaOCl_dosing = Expression(
        expr=m.chlor.fs.ideal_naocl_mixer_unit.naocl_stream.flow_mol[0] * 74.44e-3 * 1e6
        / (m.chlor.fs.ideal_naocl_chlorination_unit.control_volume.properties_out[0.0].flow_vol_phase['Liq'] * 1e3)
    )


def build(m):
    build_components(m)

    financials.add_costing_param_block(m.chlor.fs)

    # annual water production
    m.chlor.fs.annual_water_production = Expression(
        expr=pyunits.convert(
            m.chlor.fs.ideal_naocl_chlorination_unit.control_volume.properties_out[0.0].flow_vol_phase['Liq'],
            to_units=pyunits.m ** 3 / pyunits.year)
             * m.chlor.fs.costing_param.load_factor)
    costing.build_costing(m.chlor, module=financials)
    return m


def scale_initialize(m):
    # scale and initialize the mixer
    chlorination.scale_ideal_naocl_mixer(m.chlor.fs.ideal_naocl_mixer_unit)
    chlorination.initialize_ideal_naocl_mixer(m.chlor.fs.ideal_naocl_mixer_unit)

    # scale and initialize the chlorination unit
    state_args = chlorination.scale_ideal_naocl_chlorination(
        m.chlor.fs.ideal_naocl_chlorination_unit,
        m.chlor.fs.ideal_naocl_rxn_params,
        m.chlor.fs.ideal_naocl_thermo_params,
        chlorination.ideal_naocl_reaction_config)

    # initialize the unit
    chlorination.initialize_ideal_naocl_chlorination(m.chlor.fs.ideal_naocl_chlorination_unit, state_args, debug_out=False)


def initialize(m):
    pass


def report(m):
    chlorination.display_results_of_ideal_naocl_mixer(m.chlor.fs.ideal_naocl_mixer_unit)
    chlorination.display_results_of_chlorination_unit(m.chlor.fs.ideal_naocl_chlorination_unit)


def solve_flowsheet():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    build(m)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scale
    scale_initialize(m)
    # calculate_scaling_factors(m)
    # constraint_autoscale_large_jac(m)
    constraint_autoscale_large_jac(m.chlor.fs.ideal_naocl_mixer_unit)
    constraint_autoscale_large_jac(m.chlor.fs.ideal_naocl_chlorination_unit)

    # # initialize
    # initialize(m)

    check_dof(m)
    solve_with_user_scaling(m, tee=True, fail_flag=True, bound_push=1e-10, mu_init=1e-6)

    # report
    report(m)

    return m


def simulate(m):
    solve_with_user_scaling(m, tee=False, fail_flag=True, bound_push=1e-10, mu_init=1e-6)


if __name__ == "__main__":
    m = solve_flowsheet()
    m.chlor.fs.ideal_naocl_mixer_unit.naocl_stream.flow_mol[0].display()
    for x in [0.5, 1, 2, 4, 8]:
        m.chlor.fs.ideal_naocl_mixer_unit.naocl_stream.flow_mol[0].fix(x * 1e-5)
        simulate(m)
        print('molar flow (1e-5 mol/s): %.1f, free chlorine (mg/L): %.2f, OCl percent: %.1f, HOCl percent: %.1f '
              'LCOW ($/m3): %.4f, dosing (mg/L): %.1f' % (
            value(m.chlor.fs.ideal_naocl_mixer_unit.naocl_stream.flow_mol[0]*1e5),
            value(m.chlor.fs.ideal_naocl_chlorination_unit.free_chlorine),
            value(m.chlor.fs.OCl_distribution * 100),
            value(m.chlor.fs.HOCl_distribution * 100),
            value(m.chlor.fs.costing.LCOW),
            value(m.chlor.fs.NaOCl_dosing),
        ))
