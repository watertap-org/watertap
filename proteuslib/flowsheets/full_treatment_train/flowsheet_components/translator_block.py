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

"""Translator blocks for supported property packages"""

from pyomo.environ import Constraint
from idaes.generic_models.unit_models.translator import Translator
from idaes.core.util.scaling import calculate_scaling_factors, constraint_scaling_transform, get_scaling_factor
from proteuslib.flowsheets.full_treatment_train.model_components import property_models


def build_tb(m, base_inlet='ion', base_outlet='TDS', name_str=None):
    """
    Build a translator block to convert for the specified base from inlet to outlet.
    """

    if name_str is None:
        name_str = 'tb_' + base_inlet + '_to_' + base_outlet

    if base_inlet not in ['ion', 'salt']:
        raise ValueError('Unexpected property base inlet {base_inlet} for build_tb'
                         ''.format(base_inlet=base_inlet))
    prop_inlet = property_models.get_prop(m, base=base_inlet)

    if base_outlet not in ['TDS']:
        raise ValueError('Unexpected property base outlet {base_outlet} for build_tb'
                         ''.format(base_outlet=base_outlet))
    prop_outlet = property_models.get_prop(m, base=base_outlet)

    # build translator block
    setattr(m.fs, name_str, Translator(default={"inlet_property_package": prop_inlet,
                                                "outlet_property_package": prop_outlet}))
    blk = getattr(m.fs, name_str)

    # scale translator block to get scaling factors
    calculate_scaling_factors(blk)

    # add translator block constraints
    blk.eq_equal_temperature = Constraint(
        expr=blk.inlet.temperature[0]
             == blk.outlet.temperature[0])
    constraint_scaling_transform(
        blk.eq_equal_temperature,
        get_scaling_factor(blk.properties_in[0].temperature))
    blk.eq_equal_pressure = Constraint(
        expr=blk.inlet.pressure[0]
             == blk.outlet.pressure[0])
    constraint_scaling_transform(
        blk.eq_equal_pressure,
        get_scaling_factor(blk.properties_in[0].pressure))

    if base_inlet == 'ion' and base_outlet == 'TDS':
        blk.eq_H2O_balance = Constraint(
            expr=blk.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O']
                 == blk.outlet.flow_mass_phase_comp[0, 'Liq', 'H2O'])
        constraint_scaling_transform(
            blk.eq_H2O_balance,
            get_scaling_factor(blk.properties_out[0].flow_mass_phase_comp['Liq', 'H2O']))

        blk.eq_TDS_balance = Constraint(
            expr=sum(
                blk.inlet.flow_mass_phase_comp[0, 'Liq', j] for j in ['Na', 'Ca', 'Mg', 'SO4', 'Cl'])
                 == blk.outlet.flow_mass_phase_comp[0, 'Liq', 'TDS'])
        constraint_scaling_transform(
            blk.eq_TDS_balance,
            get_scaling_factor(blk.properties_out[0].flow_mass_phase_comp['Liq', 'TDS']))

    elif base_inlet == 'salt' and base_outlet == 'TDS':
        blk.eq_H2O_balance = Constraint(
            expr=blk.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O']
                 == blk.outlet.flow_mass_phase_comp[0, 'Liq', 'H2O'])
        constraint_scaling_transform(
            blk.eq_H2O_balance,
            get_scaling_factor(blk.properties_out[0].flow_mass_phase_comp['Liq', 'H2O']))

        blk.eq_TDS_balance = Constraint(
            expr=sum(
                blk.inlet.flow_mass_phase_comp[0, 'Liq', j] for j in ['NaCl', 'CaSO4', 'MgSO4', 'MgCl2'])
                 == blk.outlet.flow_mass_phase_comp[0, 'Liq', 'TDS'])
        constraint_scaling_transform(
            blk.eq_TDS_balance,
            get_scaling_factor(blk.properties_out[0].flow_mass_phase_comp['Liq', 'TDS']))

    else:
        raise ValueError('Unexpected property base combination for build_tb')

    blk.properties_in[0].mass_frac_phase_comp  # touch for initialization
    blk.properties_out[0].mass_frac_phase_comp