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
from idaes.generic_models.unit_models import Separator
from idaes.generic_models.unit_models.separator import SplittingType, EnergySplittingType
from idaes.core.util.scaling import calculate_scaling_factors
from proteuslib.flowsheets.full_treatment_train.example_models import property_models
from proteuslib.flowsheets.full_treatment_train.util import solve_with_user_scaling, check_dof


def build_RO_example(m):
    m.fs.RO = Separator(default={
        "property_package": m.fs.prop_TDS,
        "outlet_list": ['retentate', 'permeate'],
        "split_basis": SplittingType.componentFlow,
        "energy_split_basis": EnergySplittingType.equal_temperature})

    # specify
    m.fs.RO.split_fraction[0, 'permeate', 'H2O'].fix(0.5)
    m.fs.RO.split_fraction[0, 'permeate', 'TDS'].fix(0.01)

    # scaling
    calculate_scaling_factors(m.fs.RO)


def build_NF_salt_example(m):
    m.fs.NF = Separator(default={
        "property_package": m.fs.prop_salt,
        "outlet_list": ['retentate', 'permeate'],
        "split_basis": SplittingType.componentFlow,
        "energy_split_basis": EnergySplittingType.equal_temperature})

    # specify
    m.fs.NF.split_fraction[0, 'permeate', 'H2O'].fix(0.9)
    m.fs.NF.split_fraction[0, 'permeate', 'NaCl'].fix(0.9)
    m.fs.NF.split_fraction[0, 'permeate', 'CaSO4'].fix(0.1)
    m.fs.NF.split_fraction[0, 'permeate', 'MgSO4'].fix(0.1)
    m.fs.NF.split_fraction[0, 'permeate', 'MgCl2'].fix(0.2)

    # scaling
    calculate_scaling_factors(m.fs.NF)


def build_NF_ion_example(m):
    m.fs.NF = Separator(default={
        "property_package": m.fs.prop_ion,
        "outlet_list": ['retentate', 'permeate'],
        "split_basis": SplittingType.componentFlow,
        "energy_split_basis": EnergySplittingType.equal_temperature})

    # specify
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

    # scaling
    calculate_scaling_factors(m.fs.NF)


def run_example(case):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    if case == 'RO':
        property_models.build_prop_TDS(m)
        build_RO_example(m)
        unit = m.fs.RO
        # specify feed
        property_models.specify_feed_TDS(m.fs.RO.mixed_state[0])
    elif case == 'NF_salt':
        property_models.build_prop_salt(m)
        build_NF_salt_example(m)
        unit = m.fs.NF
        # specify feed
        property_models.specify_feed_salt(m.fs.NF.mixed_state[0])
    elif case == 'NF_ion':
        property_models.build_prop_ion(m)
        build_NF_ion_example(m)
        unit = m.fs.NF
        # specify feed
        property_models.specify_feed_ion(m.fs.NF.mixed_state[0])
        m.fs.NF.mixed_state[0].mass_frac_phase_comp  # touching

    check_dof(m)
    solve_with_user_scaling(m)

    unit.inlet.display()
    unit.permeate.display()
    unit.retentate.display()

    return m


if __name__ == "__main__":
    run_example('RO')
    run_example('NF_salt')
    run_example('NF_ion')

