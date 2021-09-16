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

from pyomo.environ import ConcreteModel, TransformationFactory
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models import Feed
from idaes.core.util import scaling as iscale
from idaes.core.util.initialization import propagate_state
import pretreatment_stoich_softening_block as pssb
from proteuslib.flowsheets.full_treatment_train.util import solve_with_user_scaling, check_dof

def build_feed_block(model):
    model.fs.feed = Feed(default={'property_package': model.fs.stoich_softening_thermo_params})

    comp_list = ['Na', 'Ca', 'Mg', 'SO4', 'Cl']
    feed_flow_mass = 1  # kg/s
    feed_mass_frac_comp = {'Na': 11122e-6,
                           'Ca': 382e-6,
                           'Mg': 1394e-6,
                           'SO4': 2136e-6,
                           'Cl': 20316.88e-6}
    feed_mass_frac_comp['H2O'] = 1 - sum(x for x in feed_mass_frac_comp.values())

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

    model.fs.feed.properties[0].mole_frac_comp['Ca(HCO3)2'].fix(feed_mole_frac_comp['Ca'])
    model.fs.feed.properties[0].mole_frac_comp['Ca(OH)2'].fix(0)
    model.fs.feed.properties[0].mole_frac_comp['CaCO3'].fix(0)
    model.fs.feed.properties[0].mole_frac_comp['H2O'].fix(feed_mole_frac_comp['H2O'])
    model.fs.feed.properties[0].mole_frac_comp['Mg(HCO3)2'].fix(feed_mole_frac_comp['Mg'])
    model.fs.feed.properties[0].mole_frac_comp['Mg(OH)2'].fix(0)
    model.fs.feed.properties[0].mole_frac_comp['NaCl'].fix(feed_mole_frac_comp['Na'])
    model.fs.feed.properties[0].mole_frac_comp['SO4_2-'].fix(feed_mole_frac_comp['SO4'])
    # this defined feed doesn't have all of the Cl because the apparent species could be broken down as:
    # NaCl, CaSO4, MgSO4, MgCl2, and H2O. So there is more moles of Cl than Na. This can be determined later.

    model.fs.feed.properties[0].flow_mol.fix(10)  # TODO: update to 55.56
    model.fs.feed.properties[0].pressure.fix(101325)
    model.fs.feed.properties[0].temperature.fix(298)


def run_feed_mixer():
    model = ConcreteModel()
    model.fs = FlowsheetBlock()
    pssb.build_stoich_softening_prop(model)

    build_feed_block(model)

    # ***build mixer***
    pssb.build_stoich_softening_mixer_unit(model)

    # set model values
    pssb.set_stoich_softening_mixer_inlets(model)
    # pssb.fix_stoich_softening_mixer_inlet_stream(model)  # this will be fixed by the feed block when we expand arcs
    pssb.fix_stoich_softening_mixer_lime_stream(model)

    # connect feed to mixer
    model.fs.s_prtrt_feed_mixer = Arc(
        source=model.fs.feed.outlet, destination=model.fs.stoich_softening_mixer_unit.inlet_stream)

    # transform
    TransformationFactory("network.expand_arcs").apply_to(model)

    # check degrees of freedom
    check_dof(model)

    # scale model
    # pssb.scale_stoich_softening_mixer(model.fs.stoich_softening_mixer_unit)  # not necessary, will still solve without
    # iscale.calculate_scaling_factors(model)  # this scales the feed and the transformed Arcs, also not necessary

    # initialize
    model.fs.feed.initialize()
    propagate_state(model.fs.s_prtrt_feed_mixer)
    pssb.initialize_stoich_softening_mixer(model.fs.stoich_softening_mixer_unit, debug_out=False)

    model.fs.stoich_softening_mixer_unit.outlet.mole_frac_comp.pprint()
    model.fs.stoich_softening_mixer_unit.dosing_rate.pprint()

    pssb.display_results_of_stoich_softening_mixer(model.fs.stoich_softening_mixer_unit)


def run_feed_full_system():
    model = ConcreteModel()
    model.fs = FlowsheetBlock(default={"dynamic": False})
    pssb.build_stoich_softening_prop(model)

    build_feed_block(model)

    # TODO: this build creates multiple property packages and replaces it each time.
    #  Warning: Implicitly replace the Component attribute...
    pssb.build_stoich_softening_block(model)

    # connect feed block
    model.fs.s_prtrt_feed_mixer = Arc(
        source=model.fs.feed.outlet, destination=model.fs.stoich_softening_mixer_unit.inlet_stream)

    # set model values
    pssb.set_stoich_softening_mixer_inlets(model)
    # pssb.fix_stoich_softening_mixer_inlet_stream(model)  # this will be fixed by the feed block when we expand arcs
    pssb.fix_stoich_softening_mixer_lime_stream(model)

    # expand arcs
    TransformationFactory("network.expand_arcs").apply_to(model)

    # scale model
    iscale.calculate_scaling_factors(model)
    # # scale the mixer  # all of these are not necessary
    # pssb.scale_stoich_softening_mixer(model.fs.stoich_softening_mixer_unit)
    # # scale the reactor
    # pssb.scale_stoich_softening_reactor(model.fs.stoich_softening_reactor_unit)
    # # scale the separator
    # pssb.scale_stoich_softening_separator(model.fs.stoich_softening_separator_unit)

    check_dof(model)

    # initialize feed
    model.fs.feed.initialize()
    propagate_state(model.fs.s_prtrt_feed_mixer)
    # initializer mixer
    pssb.initialize_stoich_softening_mixer(model.fs.stoich_softening_mixer_unit, debug_out=False)
    # initializer reactor
    propagate_state(model.fs.stoich_softening_arc_mixer_to_reactor)
    pssb.initialize_stoich_softening_reactor(model.fs.stoich_softening_reactor_unit, debug_out=False)
    # initializer separator
    propagate_state(model.fs.stoich_softening_arc_reactor_to_separator)
    pssb.initialize_stoich_softening_separator(model.fs.stoich_softening_separator_unit, debug_out=False)

    # solve
    solve_with_user_scaling(model)

    # display
    pssb.display_results_of_stoich_softening_mixer(model.fs.stoich_softening_mixer_unit)
    pssb.display_results_of_stoich_softening_reactor(model.fs.stoich_softening_reactor_unit)
    pssb.display_results_of_stoich_softening_separator(model.fs.stoich_softening_separator_unit)


if __name__ == "__main__":
    # run_feed_mixer()
    run_feed_full_system()



