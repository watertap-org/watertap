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
from idaes.core.util.scaling import calculate_scaling_factors


def build_tb_salt_to_TDS(m, name_str=None):
    """
    Build a translator block to convert between the prop_salt and prop_TDS.
    """

    if name_str is None:
        name_str = 'tb_salt_to_TDS'

    setattr(m.fs, name_str, Translator(default={"inlet_property_package": m.fs.prop_salt,
                                                "outlet_property_package": m.fs.prop_TDS}))
    blk = getattr(m.fs, name_str)

    # add constraints for translator block
    blk.eq_H2O_balance = Constraint(
        expr=blk.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O']
             == blk.outlet.flow_mass_phase_comp[0, 'Liq', 'H2O'])
    blk.eq_TDS_balance = Constraint(
        expr=sum(
            blk.inlet.flow_mass_phase_comp[0, 'Liq', j] for j in ['NaCl', 'CaSO4', 'MgSO4', 'MgCl2'])
             == blk.outlet.flow_mass_phase_comp[0, 'Liq', 'TDS'])
    blk.eq_equal_temperature = Constraint(
        expr=blk.inlet.temperature[0]
             == blk.outlet.temperature[0])
    blk.eq_equal_pressure = Constraint(
        expr=blk.inlet.pressure[0]
             == blk.outlet.pressure[0])

    calculate_scaling_factors(blk)


def build_tb_ion_to_TDS(m, name_str=None):
    """
    Build a translator block to convert between the prop_ion and prop_TDS.
    """

    if name_str is None:
        name_str = 'tb_ion_to_TDS'

    setattr(m.fs, name_str, Translator(default={"inlet_property_package": m.fs.prop_ion,
                                                "outlet_property_package": m.fs.prop_TDS}))
    blk = getattr(m.fs, name_str)

    # add constraints for translator block
    blk.eq_H2O_balance = Constraint(
        expr=blk.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O']
             == blk.outlet.flow_mass_phase_comp[0, 'Liq', 'H2O'])
    blk.eq_TDS_balance = Constraint(
        expr=sum(blk.inlet.flow_mass_phase_comp[0, 'Liq', j] for j in ['Na', 'Ca', 'Mg', 'SO4', 'Cl'])
             == blk.outlet.flow_mass_phase_comp[0, 'Liq', 'TDS'])
    blk.eq_equal_temperature = Constraint(
        expr=blk.inlet.temperature[0]
             == blk.outlet.temperature[0])
    blk.eq_equal_pressure = Constraint(
        expr=blk.inlet.pressure[0]
             == blk.outlet.pressure[0])

    calculate_scaling_factors(blk)


