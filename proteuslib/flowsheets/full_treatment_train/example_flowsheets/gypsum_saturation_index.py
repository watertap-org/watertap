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

from pyomo.environ import Block, Var, Expression, Constraint, value
from pyomo.environ import units as pyunits

def build_desalination_saturation(m, **kwargs):
    m.fs.desal_saturation = Block()
    m.fs.desal_saturation.properties = m.fs.prop_eNRTL.build_state_block([0], default={})
    sb_eNRTL = m.fs.desal_saturation.properties[0]

    # populate initial values
    sb_eNRTL.temperature = 298
    sb_eNRTL.pressure = 101325
    sb_eNRTL.flow_mol = 50
    # 2 times concentration factor
    sb_eNRTL.mole_frac_comp["Na_+"] = 0.017327
    sb_eNRTL.mole_frac_comp["Ca_2+"] = 0.000341
    sb_eNRTL.mole_frac_comp["Mg_2+"] = 0.002054
    sb_eNRTL.mole_frac_comp["SO4_2-"] = 0.000796
    sb_eNRTL.mole_frac_comp["Cl_-"] = 0.020529
    sb_eNRTL.mole_frac_comp["H2O"] = 0.958952

    # constraints
    sb_dilute = m.fs.tb_pretrt_to_desal.properties_in[0]
    if kwargs['is_twostage']:
        sb_perm = m.fs.mixer_permeate.mixed_state[0]
        sb_conc = m.fs.RO2.feed_side.properties_out[0]
    else:
        sb_perm = m.fs.RO.permeate_side.properties_mixed[0]
        sb_conc = m.fs.RO.feed_side.properties_out[0]

    # assumes pretreatment uses the ion property basis
    comp_match_dict = {'Na_+': 'Na',
                       'Ca_2+': 'Ca',
                       'Mg_2+': 'Mg',
                       'SO4_2-': 'SO4',
                       'Cl_-': 'Cl',
                       'H2O': 'H2O'}

    m.fs.desal_saturation.eq_temperature = Constraint(
        expr=sb_eNRTL.temperature == sb_conc.temperature)
    m.fs.desal_saturation.eq_pressure = Constraint(
        expr=sb_eNRTL.pressure == sb_conc.pressure)

    @m.fs.desal_saturation.Constraint(comp_match_dict.keys())
    def eq_flow_mol_balance(b, j):
        if j in ['Cl_-', 'Na_+']:
            return (sb_eNRTL.flow_mol_phase_comp['Liq', j]
                    == sb_dilute.flow_mol_phase_comp['Liq', comp_match_dict[j]]
                    - sb_perm.flow_mass_phase_comp['Liq', 'TDS'] / 58.44e-3)
        elif j in ['Ca_2+', 'Mg_2+', 'SO4_2-']:
            return (sb_eNRTL.flow_mol_phase_comp['Liq', j]
                    == sb_dilute.flow_mol_phase_comp['Liq', comp_match_dict[j]])
        elif j == 'H2O':
            return (sb_eNRTL.flow_mol_phase_comp['Liq', j] ==
                    sb_conc.flow_mol_phase_comp['Liq', 'H2O'])

    # ksp = 3.9e-9  # Gibbs energy gives 3.9e-8, but this fits expectations better
    ksp = 3.2e-9  # This fits expectations even better

    m.fs.desal_saturation.saturation_index = Var(
        initialize=0.5,
        bounds=(1e-8, 10),
        units=pyunits.dimensionless,
        doc="Gypsum saturation index")

    m.fs.desal_saturation.eq_saturation_index = Constraint(
        expr=m.fs.desal_saturation.saturation_index
             == sb_eNRTL.act_phase_comp["Liq", "Ca_2+"]
             * sb_eNRTL.act_phase_comp["Liq", "SO4_2-"]
             * sb_eNRTL.act_phase_comp["Liq", "H2O"] ** 2
             / ksp)
