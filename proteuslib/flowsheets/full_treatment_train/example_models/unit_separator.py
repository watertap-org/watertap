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
import proteuslib.property_models.seawater_prop_pack as seawater_prop_pack
import proteuslib.flowsheets.full_treatment_train.example_models.property_seawater_salts as property_seawater_salts
import proteuslib.flowsheets.full_treatment_train.example_models.property_seawater_ions as property_seawater_ions
import proteuslib.flowsheets.full_treatment_train.example_models.feed_specification as feed_specification
from proteuslib.flowsheets.full_treatment_train.util import solve_with_user_scaling, check_dof


def build_RO_example(m):
    m.fs.RO_properties = seawater_prop_pack.SeawaterParameterBlock()
    m.fs.RO = Separator(default={
        "property_package": m.fs.RO_properties,
        "outlet_list": ['retentate', 'permeate'],
        "split_basis": SplittingType.componentFlow,
        "energy_split_basis": EnergySplittingType.equal_temperature})

    # specifying
    # feed
    feed_specification.specify_seawater_TDS(m.fs.RO.mixed_state[0])
    # separator
    m.fs.RO.split_fraction[0, 'permeate', 'H2O'].fix(0.5)
    m.fs.RO.split_fraction[0, 'permeate', 'TDS'].fix(0.01)
    check_dof(m)

    # scaling
    calculate_scaling_factors(m.fs.RO)


def run_RO_example():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    build_RO_separator_example(m)
    solve_with_user_scaling(m)

    m.fs.RO.inlet.display()
    m.fs.RO.permeate.display()
    m.fs.RO.retentate.display()

    return m


def build_NF_salt_example(m):
    m.fs.NF_properties = property_seawater_salts.PropParameterBlock()
    m.fs.NF = Separator(default={
        "property_package": m.fs.NF_properties,
        "outlet_list": ['retentate', 'permeate'],
        "split_basis": SplittingType.componentFlow,
        "energy_split_basis": EnergySplittingType.equal_temperature})

    # specifying
    # feed
    feed_specification.specify_seawater_salts(m.fs.NF.mixed_state[0])
    # separator
    m.fs.NF.split_fraction[0, 'permeate', 'H2O'].fix(0.9)
    m.fs.NF.split_fraction[0, 'permeate', 'NaCl'].fix(0.9)
    m.fs.NF.split_fraction[0, 'permeate', 'CaSO4'].fix(0.1)
    m.fs.NF.split_fraction[0, 'permeate', 'MgSO4'].fix(0.1)
    m.fs.NF.split_fraction[0, 'permeate', 'MgCl2'].fix(0.2)
    check_dof(m)

    # scaling
    calculate_scaling_factors(m.fs.NF)



def build_NF_ion_example(m):
    m.fs.NF_properties = property_seawater_ions.PropParameterBlock()
    m.fs.NF = Separator(default={
        "property_package": m.fs.NF_properties,
        "outlet_list": ['retentate', 'permeate'],
        "split_basis": SplittingType.componentFlow,
        "energy_split_basis": EnergySplittingType.equal_temperature})

    # specifying
    # feed
    feed_specification.specify_seawater_ions(m.fs.NF.mixed_state[0])
    m.fs.NF.mixed_state[0].mass_frac_phase_comp  # touching
    # separator
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
    check_dof(m)

    # scaling
    m.fs.NF_properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
    m.fs.NF_properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'Na'))
    m.fs.NF_properties.set_default_scaling('flow_mass_phase_comp', 1e4, index=('Liq', 'Ca'))
    m.fs.NF_properties.set_default_scaling('flow_mass_phase_comp', 1e3, index=('Liq', 'Mg'))
    m.fs.NF_properties.set_default_scaling('flow_mass_phase_comp', 1e3, index=('Liq', 'SO4'))
    m.fs.NF_properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'Cl'))
    calculate_scaling_factors(m.fs.NF)


def run_example(func, unit_str):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    func(m)
    solve_with_user_scaling(m)

    blk = getattr(m.fs, unit_str)
    blk.inlet.display()
    blk.permeate.display()
    blk.retentate.display()

    return m


if __name__ == "__main__":
    run_example(build_RO_example, 'RO')
    run_example(build_NF_salt_example, 'NF')
    run_example(build_NF_ion_example, 'NF')

