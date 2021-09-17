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

"""Pretreatment flowsheet components"""

from pyomo.environ import ConcreteModel, TransformationFactory
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models import Feed
from proteuslib.flowsheets.full_treatment_train.example_models import Separator, Mixer
from idaes.generic_models.unit_models.separator import SplittingType, EnergySplittingType
from idaes.core.util.scaling import (calculate_scaling_factors,
                                     set_scaling_factor,
                                     get_scaling_factor,
                                     constraint_scaling_transform)
from idaes.core.util.initialization import propagate_state
from proteuslib.unit_models.pump_isothermal import Pump
from proteuslib.flowsheets.full_treatment_train.example_flowsheets import feed_block
from proteuslib.flowsheets.full_treatment_train.example_models import unit_separator, unit_ZONF, property_models
from proteuslib.flowsheets.full_treatment_train.util import solve_with_user_scaling, check_dof

from pyomo.environ import ConcreteModel, TransformationFactory, Constraint
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models import Feed, Translator
from idaes.core.util.scaling import calculate_scaling_factors, constraint_scaling_transform, get_scaling_factor
from idaes.core.util.initialization import propagate_state
from proteuslib.flowsheets.full_treatment_train.chemistry_flowsheets import pretreatment_stoich_softening_block as pssb
from proteuslib.flowsheets.full_treatment_train.util import solve_with_user_scaling, check_dof
from proteuslib.flowsheets.full_treatment_train.example_models import property_models


def build_pretreatment_softening(m):
    """
    Builds softening pretreatment including specified feed and auxiliary equipment.
    """
    pretrt_port = {}

    property_models.build_prop(m, base='TDS')
    pssb.build_stoich_softening_prop(m)

    build_feed_block(m)
    # build_tb(m)
    pssb.build_stoich_softening_block(m)

    # connect feed block
    m.fs.s_prtrt_feed_mixer = Arc(
        source=m.fs.feed.outlet, destination=m.fs.stoich_softening_mixer_unit.inlet_stream)
    # m.fs.s_prtrt_sep_tb = Arc(
    #     source=m.fs.stoich_softening_separator_unit.outlet_stream, destination=m.fs.tb_soft_to_TDS.inlet)

    # set model values
    pssb.set_stoich_softening_mixer_inlets(m, dosing_rate_of_lime_mg_per_s=500)
    pssb.fix_stoich_softening_mixer_lime_stream(m)

    # set ports
    # pretrt_port['out'] = m.fs.mixer.outlet
    # pretrt_port['waste'] = m.fs.NF.retentate




def build_feed_block(m):
    m.fs.feed = Feed(default={'property_package': m.fs.stoich_softening_thermo_params})

    comp_list = ['Na', 'Ca', 'Mg', 'SO4', 'Cl']
    feed_flow_mass = 1  # kg/s
    feed_mass_frac_comp = {'Na': 11122e-6,
                           'Ca': 382e-6,
                           'Mg': 1394e-6,
                           'SO4': 2136e-6,
                           'Cl': 20316.88e-6,
                           'H2O': 1 - sum(x for x in comp_list)}

    mw_comp = {'H2O': 18.015e-3,
               'Na': 22.990e-3,
               'Ca': 40.078e-3,
               'Mg': 24.305e-3,
               'SO4': 96.06e-3,
               'Cl': 35.453e-3}

    feed_flow_mol_comp = {}
    for j in feed_mass_frac_comp:
        feed_flow_mol_comp[j] = feed_flow_mass * feed_mass_frac_comp[j] / mw_comp[j]

    feed_mole_frac_comp = {}
    for j in feed_flow_mol_comp:
        feed_mole_frac_comp[j] = feed_flow_mol_comp[j] / sum(feed_flow_mol_comp[jj] for jj in feed_flow_mol_comp)

    m.fs.feed.properties[0].mole_frac_comp['Ca(HCO3)2'].fix(feed_mole_frac_comp['Ca'])
    m.fs.feed.properties[0].mole_frac_comp['Ca(OH)2'].fix(0)
    m.fs.feed.properties[0].mole_frac_comp['CaCO3'].fix(0)
    m.fs.feed.properties[0].mole_frac_comp['H2O'].fix(feed_mole_frac_comp['H2O'])
    m.fs.feed.properties[0].mole_frac_comp['Mg(HCO3)2'].fix(feed_mole_frac_comp['Mg'])
    m.fs.feed.properties[0].mole_frac_comp['Mg(OH)2'].fix(0)
    m.fs.feed.properties[0].mole_frac_comp['NaCl'].fix(feed_mole_frac_comp['Na'])
    m.fs.feed.properties[0].mole_frac_comp['SO4_2-'].fix(feed_mole_frac_comp['SO4'])
    # this defined feed doesn't have all of the Cl because the apparent species could be broken down as:
    # NaCl, CaSO4, MgSO4, MgCl2, and H2O. So there is more moles of Cl than Na. This can be determined later.

    m.fs.feed.properties[0].flow_mol.fix(55.56)
    m.fs.feed.properties[0].pressure.fix(101325)
    m.fs.feed.properties[0].temperature.fix(298)


def build_tb(m):
    # build translator block
    m.fs.tb_soft_to_TDS = Translator(default={"inlet_property_package": m.fs.stoich_softening_thermo_params,
                                                "outlet_property_package": m.fs.prop_TDS})
    blk = m.fs.tb_soft_to_TDS

    # add translator block constraints
    blk.eq_equal_temperature = Constraint(
        expr=blk.properties_in[0].temperature
             == blk.properties_out[0].temperature)
    blk.eq_equal_pressure = Constraint(
        expr=blk.properties_in[0].pressure
             == blk.properties_out[0].pressure)
    blk.eq_H2O_balance = Constraint(
        expr=blk.properties_in[0].flow_mol * blk.properties_in[0].mole_frac_comp['H2O']
             == blk.properties_out[0].flow_mol_phase_comp['Liq', 'H2O'])

    mw_comp = {'H2O': 18.015e-3,
               'Na': 22.990e-3,
               'Ca': 40.078e-3,
               'Mg': 24.305e-3,
               'SO4': 96.06e-3,
               'Cl': 35.453e-3}

    # TODO: this does not catch the Cl that is not associated with with Na
    blk.eq_TDS_balance = Constraint(
        expr=
        blk.properties_out[0].flow_mass_phase_comp['Liq', 'TDS']
        == blk.properties_in[0].flow_mol * blk.properties_in[0].mole_frac_comp['Ca(HCO3)2'] * mw_comp['Ca']
        + blk.properties_in[0].flow_mol * blk.properties_in[0].mole_frac_comp['Mg(HCO3)2'] * mw_comp['Mg']
        + blk.properties_in[0].flow_mol * blk.properties_in[0].mole_frac_comp['NaCl'] * (mw_comp['Na'] + mw_comp['Cl'])
        + blk.properties_in[0].flow_mol * blk.properties_in[0].mole_frac_comp['SO4_2-'] * mw_comp['SO4']
        )

    # scale translator block to get scaling factors
    calculate_scaling_factors(m)
    constraint_scaling_transform(blk.eq_equal_temperature, get_scaling_factor(blk.properties_out[0].temperature))
    constraint_scaling_transform(blk.eq_equal_pressure, get_scaling_factor(blk.properties_out[0].pressure))
    constraint_scaling_transform(blk.eq_H2O_balance,
        get_scaling_factor(blk.properties_out[0].flow_mol_phase_comp['Liq', 'H2O']))
    constraint_scaling_transform( blk.eq_TDS_balance,
        get_scaling_factor(blk.properties_out[0].flow_mass_phase_comp['Liq', 'TDS']))

    # touch variables
    m.fs.tb_soft_to_TDS.properties_out[0].mass_frac_phase_comp


def scale_pretreatment_softening(m):
    calculate_scaling_factors(m)


def initialize_pretreatment_softening(m):
    # initialize feed
    m.fs.feed.initialize()
    propagate_state(m.fs.s_prtrt_feed_mixer)
    # initializer mixer
    pssb.initialize_stoich_softening_mixer(m.fs.stoich_softening_mixer_unit, debug_out=False)
    # initializer reactor
    propagate_state(m.fs.stoich_softening_arc_mixer_to_reactor)
    pssb.initialize_stoich_softening_reactor(m.fs.stoich_softening_reactor_unit, debug_out=False)
    # initializer separator
    propagate_state(m.fs.stoich_softening_arc_reactor_to_separator)
    pssb.initialize_stoich_softening_separator(m.fs.stoich_softening_separator_unit, debug_out=False)
    # # initialize translator block
    # propagate_state(m.fs.s_prtrt_sep_tb)
    # m.fs.tb_soft_to_TDS.initialize()


def display_pretreatment_softening(m):
    m.fs.feed.report()
    pssb.display_results_of_stoich_softening_mixer(m.fs.stoich_softening_mixer_unit)
    pssb.display_results_of_stoich_softening_reactor(m.fs.stoich_softening_reactor_unit)
    pssb.display_results_of_stoich_softening_separator(m.fs.stoich_softening_separator_unit)


def solve_pretreatment_softening():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    build_pretreatment_softening(m)
    TransformationFactory("network.expand_arcs").apply_to(m)

    scale_pretreatment_softening(m)

    initialize_pretreatment_softening(m)

    check_dof(m)

    # solve
    solve_with_user_scaling(m)

    # display
    display_pretreatment_softening(m)


if __name__ == "__main__":
    solve_pretreatment_softening()
