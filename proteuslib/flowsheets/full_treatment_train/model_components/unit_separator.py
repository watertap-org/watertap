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

from pyomo.environ import ConcreteModel, Constraint
from idaes.core import FlowsheetBlock
# from idaes.generic_models.unit_models import Separator  # replaced separator
from idaes.generic_models.unit_models.separator import SplittingType, EnergySplittingType
from idaes.core.util.scaling import calculate_scaling_factors, set_scaling_factor, constraint_scaling_transform
from proteuslib.flowsheets.full_treatment_train.model_components import property_models
from proteuslib.flowsheets.full_treatment_train.util import solve_with_user_scaling, check_dof
from proteuslib.flowsheets.full_treatment_train.model_components import Separator


def build_SepRO(m, base='TDS'):
    """
    Builds RO model based on the IDAES separator.
    Requires prop_TDS property package.
    """
    prop = property_models.get_prop(m, base=base)

    m.fs.RO = Separator(default={
        "property_package": prop,
        "outlet_list": ['retentate', 'permeate'],
        "split_basis": SplittingType.componentFlow,
        "energy_split_basis": EnergySplittingType.equal_temperature})

    # specify
    if base == 'TDS':
        m.fs.RO.split_fraction[0, 'permeate', 'H2O'].fix(0.5)
        m.fs.RO.split_fraction[0, 'permeate', 'TDS'].fix(0.01)
    else:
        raise ValueError('Unexpected property base {base} provided to build_SepRO'
                         ''.format(base=base))

    # scale
    set_scaling_factor(m.fs.RO.split_fraction, 1)  # TODO: IDAES should set these scaling factors by default
    constraint_scaling_transform(m.fs.RO.sum_split_frac[0.0, 'H2O'], 1)
    constraint_scaling_transform(m.fs.RO.sum_split_frac[0.0, 'TDS'], 1)


def build_SepNF(m, base='ion'):
    """
    Builds NF model based on the IDAES separator for a specified property base.
    Requires prop_ion or prop_salt property package.
    """
    prop = property_models.get_prop(m, base=base)

    m.fs.NF = Separator(default={
        "property_package": prop,
        "outlet_list": ['retentate', 'permeate'],
        "split_basis": SplittingType.componentFlow,
        "energy_split_basis": EnergySplittingType.equal_temperature})

    # specify
    if base == 'ion':
        m.fs.NF.split_fraction[0, 'permeate', 'H2O'].fix(0.9)
        m.fs.NF.split_fraction[0, 'permeate', 'Na'].fix(0.9)
        m.fs.NF.split_fraction[0, 'permeate', 'Ca'].fix(0.1)
        m.fs.NF.split_fraction[0, 'permeate', 'Mg'].fix(0.1)
        m.fs.NF.split_fraction[0, 'permeate', 'SO4'].fix(0.1)
        # Cl split fraction determined through electro-neutrality for the retentate
        charge_dict = {'Na': 1, 'Ca': 2, 'Mg': 2, 'SO4': -2, 'Cl': -1}
        m.fs.NF.EN_out = Constraint(
            expr=0 ==
                 sum(charge_dict[j] * m.fs.NF.retentate_state[0].flow_mol_phase_comp['Liq', j]
                     for j in charge_dict))
        constraint_scaling_transform(m.fs.NF.EN_out, 1)
    elif base == 'salt':
        m.fs.NF.split_fraction[0, 'permeate', 'H2O'].fix(0.9)
        m.fs.NF.split_fraction[0, 'permeate', 'NaCl'].fix(0.9)
        m.fs.NF.split_fraction[0, 'permeate', 'CaSO4'].fix(0.1)
        m.fs.NF.split_fraction[0, 'permeate', 'MgSO4'].fix(0.1)
        m.fs.NF.split_fraction[0, 'permeate', 'MgCl2'].fix(0.2)

    # scale
    set_scaling_factor(m.fs.NF.split_fraction, 1)  # TODO: IDAES should set these scaling factors by default
    if base == 'ion':
        constraint_scaling_transform(m.fs.NF.sum_split_frac[0.0, 'H2O'], 1)
        constraint_scaling_transform(m.fs.NF.sum_split_frac[0.0, 'Na'], 1)
        constraint_scaling_transform(m.fs.NF.sum_split_frac[0.0, 'Ca'], 1)
        constraint_scaling_transform(m.fs.NF.sum_split_frac[0.0, 'Mg'], 1)
        constraint_scaling_transform(m.fs.NF.sum_split_frac[0.0, 'SO4'], 1)
        constraint_scaling_transform(m.fs.NF.sum_split_frac[0.0, 'Cl'], 1)
    elif base == 'salt':
        constraint_scaling_transform(m.fs.NF.sum_split_frac[0.0, 'H2O'], 1)
        constraint_scaling_transform(m.fs.NF.sum_split_frac[0.0, 'NaCl'], 1)
        constraint_scaling_transform(m.fs.NF.sum_split_frac[0.0, 'CaSO4'], 1)
        constraint_scaling_transform(m.fs.NF.sum_split_frac[0.0, 'MgSO4'], 1)
        constraint_scaling_transform(m.fs.NF.sum_split_frac[0.0, 'MgCl2'], 1)


def solve_SepRO(base='TDS'):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    property_models.build_prop(m, base=base)
    build_SepRO(m, base=base)
    property_models.specify_feed(m.fs.RO.mixed_state[0], base=base)

    check_dof(m)
    calculate_scaling_factors(m)
    solve_with_user_scaling(m)

    m.fs.RO.inlet.display()
    m.fs.RO.permeate.display()
    m.fs.RO.retentate.display()

    return m


def solve_SepNF(base='ion'):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    property_models.build_prop(m, base=base)
    build_SepNF(m, base=base)
    property_models.specify_feed(m.fs.NF.mixed_state[0], base=base)

    m.fs.NF.mixed_state[0].mass_frac_phase_comp  # touching for tests
    check_dof(m)
    calculate_scaling_factors(m)
    solve_with_user_scaling(m)

    m.fs.NF.inlet.display()
    m.fs.NF.permeate.display()
    m.fs.NF.retentate.display()

    return m


if __name__ == "__main__":
    solve_SepRO(base='TDS')
    solve_SepNF(base='ion')
    solve_SepNF(base='salt')
