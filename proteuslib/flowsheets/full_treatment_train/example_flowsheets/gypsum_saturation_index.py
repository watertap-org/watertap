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

from pyomo.environ import Block, Var, Expression, Constraint, value, units as pyunits
from pyomo.environ import units as pyunits

def build(m, section='desalination', pretrt_type='NF', **kwargs):
    if section == 'desalination':
        m.fs.desal_saturation = Block()
        m.fs.desal_saturation.properties = m.fs.prop_eNRTL.build_state_block([0], default={})
        sb_eNRTL = m.fs.desal_saturation.properties[0]
    elif section == 'pretreatment':
        m.fs.pretrt_saturation = Block()
        m.fs.pretrt_saturation.properties = m.fs.prop_eNRTL.build_state_block([0], default={})
        sb_eNRTL = m.fs.pretrt_saturation.properties[0]
    else:
        raise ValueError('{section} is not an expected section for building the saturation index'
                         ''.format(section=section))

    # populate initial values
    populate_eNRTL_state_vars(sb_eNRTL, base='FpcTP')

    # ksp = 3.9e-9  # Gibbs energy gives 3.9e-8, but this fits expectations better
    ksp = 3.2e-9  # This fits expectations even better

    # constraints
    if section == 'desalination':
        sb_dilute = m.fs.tb_pretrt_to_desal.properties_in[0]
        if kwargs['is_twostage']:
            sb_perm = m.fs.mixer_permeate.mixed_state[0]
            if kwargs['RO_type'] == '0D':
                sb_conc = m.fs.RO2.feed_side.properties_out[0]
                sb_conc_inter = m.fs.RO2.feed_side.properties_interface_out[0]
            elif kwargs['RO_type'] == '1D':
                sb_conc = m.fs.RO2.feed_side.properties[0, 1]
                sb_conc_inter = m.fs.RO2.feed_side.properties_interface[0, 1]
        else:
            if kwargs['RO_type'] == '0D':
                sb_perm = m.fs.RO.permeate_side.properties_mixed[0]
                sb_conc = m.fs.RO.feed_side.properties_out[0]
                sb_conc_inter = m.fs.RO.feed_side.properties_interface_out[0]
            elif kwargs['RO_type'] == '1D':
                sb_perm = m.fs.RO.mixed_permeate[0]
                sb_conc = m.fs.RO.feed_side.properties[0, 1]
                sb_conc_inter = m.fs.RO.feed_side.properties_interface[0, 1]

        m.fs.desal_saturation.cp_modulus = Expression(
            expr=sb_conc_inter.conc_mass_phase_comp['Liq', 'TDS'] / sb_conc.conc_mass_phase_comp['Liq', 'TDS'])

        # constraints
        m.fs.desal_saturation.eq_temperature = Constraint(
            expr=sb_eNRTL.temperature == sb_conc.temperature)
        m.fs.desal_saturation.eq_pressure = Constraint(
            expr=sb_eNRTL.pressure == sb_conc.pressure)

        if pretrt_type == 'NF':
            # assumes pretreatment uses the ion property basis
            comp_match_dict = {'Na_+': 'Na',
                               'Ca_2+': 'Ca',
                               'Mg_2+': 'Mg',
                               'SO4_2-': 'SO4',
                               'Cl_-': 'Cl',
                               'H2O': 'H2O'}

            @m.fs.desal_saturation.Constraint(comp_match_dict.keys())
            def eq_flow_mol_balance(b, j):
                if j in ['Cl_-', 'Na_+']:
                    bulk_flow = (sb_dilute.flow_mol_phase_comp['Liq', comp_match_dict[j]]
                                 - sb_perm.flow_mass_phase_comp['Liq', 'TDS'] / (58.44e-3 * pyunits.kg / pyunits.mol))
                    return (sb_eNRTL.flow_mol_phase_comp['Liq', j]
                            == bulk_flow * m.fs.desal_saturation.cp_modulus)
                elif j in ['Ca_2+', 'Mg_2+', 'SO4_2-']:
                    return (sb_eNRTL.flow_mol_phase_comp['Liq', j]
                            == sb_dilute.flow_mol_phase_comp['Liq', comp_match_dict[j]]
                            * m.fs.desal_saturation.cp_modulus)
                elif j == 'H2O':
                    return (sb_eNRTL.flow_mol_phase_comp['Liq', j] ==
                            sb_conc.flow_mol_phase_comp['Liq', 'H2O'])
        elif pretrt_type == 'softening':
            # assumes pretreatment uses the softening property basis
            comp_match_dict = {'Na_+': 'NaCl',
                               'Ca_2+': 'Ca(HCO3)2',
                               'Mg_2+': 'Mg(HCO3)2',
                               'SO4_2-': 'SO4_2-',
                               'Cl_-': 'Cl_-',
                               'H2O': 'H2O'}

            @m.fs.desal_saturation.Constraint(comp_match_dict.keys())
            def eq_flow_mol_balance(b, j):
                if j == 'Na_+':
                    bulk_flow = (sb_dilute.flow_mol_phase_comp['Liq', 'NaCl']
                                 - sb_perm.flow_mass_phase_comp['Liq', 'TDS'] / (58.44e-3 * pyunits.kg / pyunits.mol))
                    return (sb_eNRTL.flow_mol_phase_comp['Liq', j]
                            == bulk_flow * m.fs.desal_saturation.cp_modulus)
                if j in 'Cl_-':
                    bulk_flow = (sb_dilute.flow_mol_phase_comp['Liq', 'Cl_-']
                            + sb_dilute.flow_mol_phase_comp['Liq', 'NaCl']
                            - sb_perm.flow_mass_phase_comp['Liq', 'TDS'] / (58.44e-3 * pyunits.kg / pyunits.mol))
                    return (sb_eNRTL.flow_mol_phase_comp['Liq', j]
                            == bulk_flow * m.fs.desal_saturation.cp_modulus)
                elif j in ['Ca_2+', 'Mg_2+', 'SO4_2-']:
                    return (sb_eNRTL.flow_mol_phase_comp['Liq', j]
                            == sb_dilute.flow_mol_phase_comp['Liq', comp_match_dict[j]]
                            * m.fs.desal_saturation.cp_modulus)
                elif j == 'H2O':
                    return (sb_eNRTL.flow_mol_phase_comp['Liq', j] ==
                            sb_conc.flow_mol_phase_comp['Liq', 'H2O'])

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

    elif section == 'pretreatment':
        comp_match_dict = {'Na_+': 'Na',
                           'Ca_2+': 'Ca',
                           'Mg_2+': 'Mg',
                           'SO4_2-': 'SO4',
                           'Cl_-': 'Cl',
                           'H2O': 'H2O'}

        sb_conc = m.fs.NF.feed_side.properties_out[0]

        m.fs.pretrt_saturation.eq_temperature = Constraint(
            expr=sb_eNRTL.temperature == sb_conc.temperature)
        m.fs.pretrt_saturation.eq_pressure = Constraint(
            expr=sb_eNRTL.pressure == sb_conc.pressure)

        @m.fs.pretrt_saturation.Constraint(comp_match_dict.keys())
        def eq_flow_mol_balance(b, j):
                return (sb_eNRTL.flow_mol_phase_comp['Liq', j] ==
                        sb_conc.flow_mol_phase_comp['Liq', comp_match_dict[j]])

        m.fs.pretrt_saturation.saturation_index = Var(
            initialize=0.5,
            bounds=(1e-8, 1e6),
            units=pyunits.dimensionless,
            doc="Gypsum saturation index")

        m.fs.pretrt_saturation.eq_saturation_index = Constraint(
            expr=m.fs.pretrt_saturation.saturation_index
                 == sb_eNRTL.act_phase_comp["Liq", "Ca_2+"]
                 * sb_eNRTL.act_phase_comp["Liq", "SO4_2-"]
                 * sb_eNRTL.act_phase_comp["Liq", "H2O"] ** 2
                 / ksp)


def populate_eNRTL_state_vars(blk, base='FpcTP'):
    blk.temperature = 298
    blk.pressure = 101325

    if base == 'FpcTP':
        feed_flow_mass = 1  # kg/s
        feed_mass_frac_comp = {'Na_+': 11122e-6,
                               'Ca_2+': 382e-6,
                               'Mg_2+': 1394e-6,
                               'SO4_2-': 2136e-6,
                               'Cl_-': 20316.88e-6}
        feed_mass_frac_comp['H2O'] = 1 - sum(x for x in feed_mass_frac_comp.values())

        mw_comp = {'H2O': 18.015e-3,
                   'Na_+': 22.990e-3,
                   'Ca_2+': 40.078e-3,
                   'Mg_2+': 24.305e-3,
                   'SO4_2-': 96.06e-3,
                   'Cl_-': 35.453e-3}

        for j in feed_mass_frac_comp:
            blk.flow_mol_phase_comp['Liq', j] = feed_flow_mass * feed_mass_frac_comp[j] / mw_comp[j]
            if j == 'H2O':
                blk.flow_mol_phase_comp['Liq', j] /= 2
    elif base == 'FTPx':
        blk.mole_frac_comp["Na_+"] = 0.017327
        blk.mole_frac_comp["Ca_2+"] = 0.000341
        blk.mole_frac_comp["Mg_2+"] = 0.002054
        blk.mole_frac_comp["SO4_2-"] = 0.000796
        blk.mole_frac_comp["Cl_-"] = 0.020529
        blk.mole_frac_comp["H2O"] = 0.958952
