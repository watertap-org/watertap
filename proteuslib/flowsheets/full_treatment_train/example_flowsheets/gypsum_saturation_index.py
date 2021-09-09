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

"""eNRTL property blocks to calculate and constrain gypsum saturation index"""

from pyomo.environ import Block, Expression, Constraint


def build_desalination_saturation(m, pretrt_port=None, desal_por=None):
    m.fs.desal_saturation = Block()
    m.fs.desal_saturation.properties = m.fs.prop_eNRTL.build_state_block([0], default={})

    # populate initial values
    m.fs.desal_saturation.properties[0].temperature = 298.15
    m.fs.desal_saturation.properties[0].pressure = 100 * 1e5
    m.fs.desal_saturation.properties[0].flow_mol = 100
    # 2 times concentration factor
    m.fs.desal_saturation.properties[0].mole_frac_comp["Na_+"] = 0.017327
    m.fs.desal_saturation.properties[0].mole_frac_comp["Ca_2+"] = 0.000341
    m.fs.desal_saturation.properties[0].mole_frac_comp["Mg_2+"] = 0.002054
    m.fs.desal_saturation.properties[0].mole_frac_comp["SO4_2-"] = 0.000796
    m.fs.desal_saturation.properties[0].mole_frac_comp["Cl_-"] = 0.020529
    m.fs.desal_saturation.properties[0].mole_frac_comp["H2O"] = 0.958952

    # constraints
    blk = m.fs.desal_saturation.properties[0]
    sb_dilute = m.fs.tb_pretrt_to_desal.properties_in[0]
    sb_conc = m.fs.RO.feed_side.properties_out[0]

    m.fs.desal_saturation.eq_temperature = Constraint(
        expr=blk.temperature == sb_conc.temperature)
    m.fs.desal_saturation.eq_pressure = Constraint(
        expr=blk.pressure == sb_conc.pressure)

    @m.fs.desal_saturation.Constraint(['Na_+', 'Ca_2+', 'Mg_2+', 'Cl_-', 'SO4_2-', 'H2O'])
    def eq_flow_mol_balance(b, j):
        if j == 'Cl_-':
            return blk.flow_mol_phase_comp['Liq', j] == sb_dilute.flow_mol_phase_comp['Liq', 'Cl'] * 0.99
        elif j == 'Na_+':
            return blk.flow_mol_phase_comp['Liq', j] == sb_dilute.flow_mol_phase_comp['Liq', 'Na'] * 0.99
        elif j == 'Ca_2+':
            return blk.flow_mol_phase_comp['Liq', j] == sb_dilute.flow_mol_phase_comp['Liq', 'Ca']
        elif j == 'Mg_2+':
            return blk.flow_mol_phase_comp['Liq', j] == sb_dilute.flow_mol_phase_comp['Liq', 'Mg']
        elif j == 'SO4_2-':
            return blk.flow_mol_phase_comp['Liq', j] == sb_dilute.flow_mol_phase_comp['Liq', 'SO4']
        elif j == 'H2O':
            return (blk.flow_mol_phase_comp['Liq', j] ==
                    sb_dilute.flow_mol_phase_comp['Liq', 'H2O']
                    * sb_conc.flow_vol / sb_dilute.flow_vol)

    Ksp = {"Gypsum": 3.9e-9}  # Gibbs energy gives 3.9e-8, but this fits expectations better

    blk.saturation_index = Expression(
        expr=blk.act_phase_comp["Liq", "Ca_2+"]
             * blk.act_phase_comp["Liq", "SO4_2-"]
             * blk.act_phase_comp["Liq", "H2O"] ** 2
             / Ksp["Gypsum"])