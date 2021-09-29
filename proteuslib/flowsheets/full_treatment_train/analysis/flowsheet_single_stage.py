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
'''
mutable parameters for optimization:
    m.fs.system_recovery_target
'''
from pyomo.environ import ConcreteModel, Objective, Expression, Constraint, TransformationFactory, value, Param, Var
from pyomo.environ import units as pyunits
from pyomo.network import Arc

from idaes.core import FlowsheetBlock
from idaes.core.util.initialization import propagate_state
from idaes.core.util.scaling import calculate_scaling_factors

from proteuslib.flowsheets.full_treatment_train.flowsheet_components import (desalination,
                                                                             feed_block,
                                                                             gypsum_saturation_index,
                                                                             translator_block,
                                                                             costing,
                                                                             financials, )
from proteuslib.flowsheets.full_treatment_train.model_components import property_models
from proteuslib.flowsheets.full_treatment_train.util import (solve_with_user_scaling,
                                                             check_dof)


desal_kwargs = {'has_desal_feed': False, 'is_twostage': False, 'has_ERD': True,
        'RO_type': '0D', 'RO_base': 'TDS', 'RO_level': 'detailed'}


def build_components(m, pretrt_type='NF', **kwargs):
    kwargs_desalination = {k: kwargs[k] for k in ('has_desal_feed', 'is_twostage', 'has_ERD',
                                                  'RO_type', 'RO_base', 'RO_level',)}

    desal_port = desalination.build_desalination(m, **kwargs_desalination)
    m.fs.s_tb_desal = Arc(source=m.fs.tb_pretrt_to_desal.outlet, destination=desal_port['in'])

    if pretrt_type == 'softening':
        property_models.build_prop(m, base='eNRTL')
    gypsum_saturation_index.build(m, section='desalination', pretrt_type=pretrt_type,  **kwargs)

    m.fs.RO.area.fix(80)
    m.fs.pump_RO.control_volume.properties_out[0].pressure.fix(60e5)

    if kwargs['is_twostage']:
        m.fs.RO2.area.fix(20)
        m.fs.pump_RO2.control_volume.properties_out[0].pressure.fix(90e5)

    # touch some properties used in optimization
    if kwargs['is_twostage']:
            product_water_sb = m.fs.mixer_permeate.mixed_state[0]
    else:
        if kwargs['RO_type'] == '0D':
            product_water_sb = m.fs.RO.permeate_side.properties_mixed[0]
        elif kwargs['RO_type'] == '1D':
            product_water_sb = m.fs.RO.mixed_permeate[0]

    feed_flow_vol = 0.0009769808  # value of feed flowrate using the seawater property package with 1 kg/s mass flowrate
    m.fs.system_recovery = Expression(
        expr=product_water_sb.flow_vol / feed_flow_vol)

    # need load factor from costing_param_block for annual_water_production
    financials.add_costing_param_block(m.fs)

    # RO recovery
    m.fs.RO_recovery = Var(initialize=0.5,
                           bounds=(0.01, 0.99),
                           doc='Total volumetric water recovery for RO')
    m.fs.eq_RO_recovery = Constraint(
        expr=m.fs.RO_recovery
             == product_water_sb.flow_vol / m.fs.tb_pretrt_to_desal.properties_out[0].flow_vol)

    # annual water production
    m.fs.annual_water_production = Expression(
        expr=pyunits.convert(product_water_sb.flow_vol, to_units=pyunits.m ** 3 / pyunits.year)
             * m.fs.costing_param.load_factor)
    costing.build_costing(m, module=financials, **kwargs)

    return desal_port


def build(m, **kwargs):
    """
    build an RO
    """
    assert not kwargs['has_desal_feed']
    property_models.build_prop(m, base='ion')
    feed_block.build_feed(m, base='ion')

    property_models.build_prop(m, base=kwargs['RO_base'])
    translator_block.build_tb(m, base_inlet='ion',
                              base_outlet=kwargs['RO_base'],
                              name_str='tb_pretrt_to_desal')
    m.fs.s_pretrt_tb = Arc(source=m.fs.feed.outlet, destination=m.fs.tb_pretrt_to_desal.inlet)

    property_models.build_prop(m, base='eNRTL')
    desal_port = build_components(m, **kwargs)


def scale(m, **kwargs):
    calculate_scaling_factors(m.fs.s_tb_desal)
    desalination.scale_desalination(m, **kwargs)


def initialize(m, **kwargs):
    propagate_state(m.fs.s_tb_desal)
    desalination.initialize_desalination(m, **kwargs)
    m.fs.desal_saturation.properties.initialize()


def report(m, **kwargs):
    m.fs.tb_pretrt_to_desal.report()
    desalination.display_desalination(m, **kwargs)
    print('desalination solubility index:', value(m.fs.desal_saturation.saturation_index))
    print('water recovery:', value(m.fs.system_recovery))
    costing.display_costing(m)


def set_optimization_components(m, system_recovery, **kwargs):
    '''
    adds max_saturation_index and system_recovery_target
    '''
    m.fs.pump_RO.control_volume.properties_out[0].pressure.unfix()
    m.fs.pump_RO.control_volume.properties_out[0].pressure.setlb(20e5)
    m.fs.pump_RO.control_volume.properties_out[0].pressure.setub(75e5)

    m.fs.RO.area.unfix()
    m.fs.RO.area.setlb(10)
    m.fs.RO.area.setub(300)

    if kwargs['RO_type'] == '0D':
        m.fs.RO.N_Re_io[0, 'in'].unfix()
    elif kwargs['RO_type'] == '1D':
        m.fs.RO.N_Re[0, 0].unfix()

    # Set lower bound for water flux at the RO outlet, based on a minimum net driving pressure, NDPmin
    m.fs.RO.NDPmin = Param(initialize=1e5, mutable=True, units=pyunits.Pa)
    if kwargs['RO_type'] == '0D':
        m.fs.RO.flux_mass_io_phase_comp[0, 'out', 'Liq', 'H2O'].setlb(m.fs.RO.A_comp[0, 'H2O']
                                                                      * m.fs.RO.dens_solvent
                                                                      * m.fs.RO.NDPmin)
    elif kwargs['RO_type'] == '1D':
        m.fs.RO.flux_mass_phase_comp[0, 1, 'Liq', 'H2O'].setlb(m.fs.RO.A_comp[0, 'H2O']
                                                               * m.fs.RO.dens_solvent
                                                               * m.fs.RO.NDPmin)

    # saturation index
    m.fs.max_saturation_index = Param(initialize=1.0, mutable=True)
    m.fs.eq_max_saturation_index_desal = Constraint(
        expr=m.fs.desal_saturation.saturation_index <= m.fs.max_saturation_index)

    m.fs.system_recovery_target = Param(initialize=system_recovery, mutable=True)
    m.fs.system_recovery_tol = Param(initialize=5e-3, mutable=True)
    m.fs.eq_system_recovery = Constraint(
        expr=(m.fs.system_recovery_target,
                m.fs.system_recovery,
                m.fs.system_recovery_target+m.fs.system_recovery_tol))

    # set objective
    m.fs.objective = Objective(expr=m.fs.costing.LCOW)


def set_up_optimization(m, system_recovery=0.50, **kwargs):
    set_optimization_components(m, system_recovery, **kwargs)
    check_dof(m, 3)


def optimize(m):
    solve_with_user_scaling(m, tee=False, fail_flag=True)


def solve_flowsheet(**desal_kwargs):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    build(m, **desal_kwargs)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scale
    calculate_scaling_factors(m.fs.tb_pretrt_to_desal)
    scale(m, **desal_kwargs)
    calculate_scaling_factors(m)

    # initialize
    m.fs.feed.initialize()
    propagate_state(m.fs.s_pretrt_tb)
    optarg = {'nlp_scaling_method': 'user-scaling'}
    m.fs.tb_pretrt_to_desal.initialize(optarg=optarg)
    initialize(m, **desal_kwargs)

    check_dof(m)
    solve_with_user_scaling(m, tee=False, fail_flag=True)

    # report
    report(m, **desal_kwargs)

    return m


def optimize_flowsheet(system_recovery=0.50, **kwargs):
    m = solve_flowsheet(**kwargs)
    set_up_optimization(m, system_recovery=system_recovery, **kwargs)
    optimize(m)

    report(m, **kwargs)

    return m


if __name__ == "__main__":
    import sys
    if len(sys.argv) == 1:
        m = solve_flowsheet(**desal_kwargs)
    else:
        m = optimize_flowsheet(system_recovery=float(sys.argv[1]), **desal_kwargs)
