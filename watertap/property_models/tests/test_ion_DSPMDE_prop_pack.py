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

import pytest
from pyomo.environ import ConcreteModel, assert_optimal_termination, value, Set, Param, Var, units as pyunits
from idaes.core import FlowsheetBlock
import idaes.core.util.scaling as iscale
from pyomo.util.check_units import assert_units_consistent
from watertap.property_models.ion_DSPMDE_prop_pack import DSPMDEParameterBlock, DSPMDEStateBlock, ActivityCoefficientModel

from watertap.core.util.initialization import check_dof
from idaes.core.util import get_solver

solver = get_solver()
# -----------------------------------------------------------------------------

@pytest.fixture(scope="module")
def model():
    m = ConcreteModel()

    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = DSPMDEParameterBlock(default={
        "solute_list": ["A", "B", "C", "D"],
        "diffusivity_data": {("Liq", "A"): 1e-9,
                             ("Liq", "B"): 1e-10,
                             ("Liq", "C"): 1e-7,
                             ("Liq", "D"): 1e-11},
        "mw_data": {"H2O": 18e-3,
                    "A": 10e-3,
                    "B": 25e-3,
                    "C": 100e-3,
                    "D": 25e-3},
        "stokes_radius_data": {"A": 1e-9,
                               "B": 1e-9,
                               "C": 1e-9,
                               "D": 1e-10},
        "density_data": {"H2O": 1000,
                         "A": 1200,
                         "B": 1100,
                         "C": 1010,
                         "D": 900},
        "charge": {"A": 1,
                   "B": -2,
                   "C": 2,
                   "D": -1}
    })

    return m

@pytest.mark.unit
def test_parameter_block(model):
    assert isinstance(model.fs.properties.component_list, Set)
    for j in model.fs.properties.component_list:
        assert j in ["H2O", "A", "B", "C", "D"]
    assert isinstance(model.fs.properties.solvent_set, Set)
    for j in model.fs.properties.solvent_set:
        assert j in ["H2O"]
    assert isinstance(model.fs.properties.solute_set, Set)
    for j in model.fs.properties.solute_set:
        assert j in ["A", "B", "C", "D"]

    assert isinstance(model.fs.properties.phase_list, Set)
    for j in model.fs.properties.phase_list:
        assert j in ["Liq"]

    assert model.fs.properties._state_block_class is DSPMDEStateBlock

    assert isinstance(model.fs.properties.dens_mass_comp, Param)
    assert model.fs.properties.dens_mass_comp['A'].value == 1200
    assert model.fs.properties.dens_mass_comp['B'].value == 1100
    assert model.fs.properties.dens_mass_comp['C'].value == 1010
    assert model.fs.properties.dens_mass_comp['H2O'].value == 1000

    assert isinstance(model.fs.properties.mw_comp, Param)
    assert model.fs.properties.mw_comp['A'].value == 10e-3
    assert model.fs.properties.mw_comp['B'].value == 25e-3
    assert model.fs.properties.mw_comp['C'].value == 100e-3
    assert model.fs.properties.mw_comp['H2O'].value == 18e-3

    assert isinstance(model.fs.properties.diffus_phase_comp, Param)
    assert model.fs.properties.diffus_phase_comp['Liq', 'A'].value == 1e-9
    assert model.fs.properties.diffus_phase_comp['Liq', 'B'].value == 1e-10
    assert model.fs.properties.diffus_phase_comp['Liq', 'C'].value == 1e-7

    assert isinstance(model.fs.properties.radius_stokes_comp, Param)
    assert model.fs.properties.radius_stokes_comp['A'].value == 1e-9
    assert model.fs.properties.radius_stokes_comp['B'].value == 1e-9
    assert model.fs.properties.radius_stokes_comp['C'].value == 1e-9

    assert model.fs.properties.config.activity_coefficient_model == ActivityCoefficientModel.ideal


@pytest.mark.component
def test_property_ions(model):
    m = model
    m.fs.stream = m.fs.properties.build_state_block([0], default={'defined_state': True})

    m.fs.stream[0].flow_mol_phase_comp['Liq', 'A'].fix(0.000407)
    m.fs.stream[0].flow_mol_phase_comp['Liq', 'B'].fix(0.010479)
    m.fs.stream[0].flow_mol_phase_comp['Liq', 'C'].fix(0.010479)
    m.fs.stream[0].flow_mol_phase_comp['Liq', 'D'].fix(0.000407)
    m.fs.stream[0].flow_mol_phase_comp['Liq', 'H2O'].fix(0.99046)
    m.fs.stream[0].temperature.fix(298.15)
    m.fs.stream[0].pressure.fix(101325)

    m.fs.stream[0].assert_electroneutrality()

    m.fs.stream[0].mole_frac_phase_comp

    m.fs.stream[0].flow_mass_phase_comp

    m.fs.stream[0].molality_comp
    m.fs.stream[0].pressure_osm
    m.fs.stream[0].dens_mass_phase
    m.fs.stream[0].conc_mol_phase_comp
    m.fs.stream[0].act_coeff_phase_comp

    iscale.calculate_scaling_factors(m.fs)

    assert_units_consistent(m)

    check_dof(m, fail_flag=True)

    m.fs.stream.initialize()

    results = solver.solve(m)
    assert_optimal_termination(results)

    # assert value(m.fs.stream[0].mole_frac_phase_comp['Liq', 'A']) == pytest.approx(4.068e-4,  rel=1e-3)
    # assert value(m.fs.stream[0].mole_frac_phase_comp['Liq', 'B']) == pytest.approx(1.047e-2,  rel=1e-3)
    # assert value(m.fs.stream[0].mole_frac_phase_comp['Liq', 'C']) == pytest.approx(1.047e-2,  rel=1e-3)
    # assert value(m.fs.stream[0].mole_frac_phase_comp['Liq', 'H2O']) == pytest.approx(9.786e-1,  rel=1e-3)
    #
    # assert value(m.fs.stream[0].flow_mass_phase_comp['Liq', 'A']) == pytest.approx(4.07e-6,  rel=1e-3)
    # assert value(m.fs.stream[0].flow_mass_phase_comp['Liq', 'B']) == pytest.approx(2.6198e-4,  rel=1e-3)
    # assert value(m.fs.stream[0].flow_mass_phase_comp['Liq', 'C']) == pytest.approx(1.048e-3,  rel=1e-3)
    # assert value(m.fs.stream[0].flow_mass_phase_comp['Liq', 'H2O']) == pytest.approx(1.762e-2,  rel=1e-3)

    assert value(m.fs.stream[0].conc_mass_phase_comp['Liq','A']) == pytest.approx(2.1288e-1,  rel=1e-3)

    assert value(m.fs.stream[0].conc_mol_phase_comp['Liq','A']) == pytest.approx(21.288,  rel=1e-3)
    assert value(m.fs.stream[0].molality_comp['A']) == pytest.approx(2.2829e-2,  rel=1e-3)

    assert value(m.fs.stream[0].pressure_osm) == pytest.approx(60.546e5,  rel=1e-3)

    assert value(m.fs.stream[0].dens_mass_phase['Liq']) == pytest.approx(1001.76,  rel=1e-3)

    assert value(m.fs.stream[0].act_coeff_phase_comp['Liq', 'A']) == 1

@pytest.fixture(scope="module")
def model2():
    m = ConcreteModel()

    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = DSPMDEParameterBlock(default={
        "solute_list": ["A", "B", "C", "D"]})

    return m

@pytest.mark.component
def test_property_ions(model2):
    m = model2

    stream = m.fs.stream = m.fs.properties.build_state_block([0], default={'defined_state': True})

    stream[0].flow_mol_phase_comp['Liq', 'A'].fix(0.000407)
    stream[0].flow_mol_phase_comp['Liq', 'B'].fix(0.010479)
    stream[0].flow_mol_phase_comp['Liq', 'C'].fix(0.010479)
    stream[0].flow_mol_phase_comp['Liq', 'D'].fix(0.000407)
    stream[0].flow_mol_phase_comp['Liq', 'H2O'].fix(0.99046)
    stream[0].temperature.fix(298.15)
    stream[0].pressure.fix(101325)

    stream[0].diffus_phase_comp['Liq', 'A'] = 1e-9
    stream[0].diffus_phase_comp['Liq', 'B'] = 1e-10
    stream[0].diffus_phase_comp['Liq', 'C'] = 1e-7
    stream[0].diffus_phase_comp['Liq', 'D'] = 1e-11

    stream[0].mw_comp['H2O'] = 18e-3
    stream[0].mw_comp['A'] = 10e-3
    stream[0].mw_comp['B'] = 25e-3
    stream[0].mw_comp['C'] = 100e-3
    stream[0].mw_comp['D'] = 25e-3

    stream[0].radius_stokes_comp['A'] = 1e-9
    stream[0].radius_stokes_comp['B'] = 1e-9
    stream[0].radius_stokes_comp['C'] = 1e-9
    stream[0].radius_stokes_comp['D'] = 1e-10

    stream[0].dens_mass_comp['H2O'] = 1e3
    stream[0].dens_mass_comp['A'] = 1.2e3
    stream[0].dens_mass_comp['B'] = 1.1e3
    stream[0].dens_mass_comp['C'] = 1.01e3
    stream[0].dens_mass_comp['D'] = 0.9e3

    stream[0].charge_comp['A'] = 1
    stream[0].charge_comp['B'] = -2
    stream[0].charge_comp['C'] = 2
    stream[0].charge_comp['D'] = -1

    stream[0].assert_electroneutrality()

    stream[0].mole_frac_phase_comp

    stream[0].flow_mass_phase_comp

    stream[0].molality_comp
    stream[0].pressure_osm
    stream[0].dens_mass_phase
    stream[0].conc_mol_phase_comp
    stream[0].flow_vol
    iscale.calculate_scaling_factors(m.fs)

    assert_units_consistent(m)

    check_dof(m, fail_flag=True)

    stream.initialize()

    results = solver.solve(m)
    assert_optimal_termination(results)

    assert value(stream[0].flow_vol_phase['Liq']) == pytest.approx(1.91524e-5,  rel=1e-3)

@pytest.fixture(scope="module")
def model3():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.properties = DSPMDEParameterBlock(default={
        "solute_list": ["Ca_2+", "SO4_2-", "Na_+", "Cl_-", "Mg_2+"],
        "diffusivity_data": {("Liq", "Ca_2+"): 0.792e-9,
                             ("Liq", "SO4_2-"): 1.06e-9,
                             ("Liq", "Na_+"): 1.33e-9,
                             ("Liq", "Cl_-"): 2.03e-9,
                             ("Liq", "Mg_2+"): 0.706e-9},
        "mw_data": {"H2O": 18e-3,
                    "Na_+": 23e-3,
                    "Ca_2+": 40e-3,
                    "Mg_2+": 24e-3,
                    "Cl_-": 35e-3,
                    "SO4_2-": 96e-3},
        "stokes_radius_data": {"Na_+": 0.184e-9,
                               "Ca_2+": 0.309e-9,
                               "Mg_2+": 0.347e-9,
                               "Cl_-": 0.121e-9,
                               "SO4_2-": 0.230e-9},
        "density_data": {"H2O": 1000,
                         "Na_+": 968,
                         "Ca_2+": 1550,
                         "Mg_2+": 1738,
                         "Cl_-": 3214, #todo: verify valid value; Cl2 specific gravity= 1.41 @ 20 C
                         "SO4_2-": 2553},
        "charge": {"Na_+": 1,
                   "Ca_2+": 2,
                   "Mg_2+": 2,
                   "Cl_-": -1,
                   "SO4_2-": -2},
        })
    return m

@pytest.mark.component
def test_seawater_data(model3):
    m = model3
    stream = m.fs.stream = m.fs.properties.build_state_block([0], default={'defined_state': True})

    mass_flow_in = 1 * pyunits.kg/pyunits.s
    feed_mass_frac = {'Na_+': 11122e-6,
                      'Ca_2+': 382e-6,
                      'Mg_2+': 1394e-6,
                      'SO4_2-': 2136e-6,
                      'Cl_-': 20300e-6}
    for ion, x in feed_mass_frac.items():
        mol_comp_flow = x * pyunits.kg/pyunits.kg * mass_flow_in / stream[0].mw_comp[ion]
        stream[0].flow_mol_phase_comp['Liq', ion].fix(mol_comp_flow)

    H2O_mass_frac = 1 - sum(x for x in feed_mass_frac.values())
    H2O_mol_comp_flow = H2O_mass_frac * pyunits.kg/pyunits.kg * mass_flow_in / stream[0].mw_comp['H2O']

    stream[0].flow_mol_phase_comp['Liq', 'H2O'].fix(H2O_mol_comp_flow)
    stream[0].temperature.fix(298.15)
    stream[0].pressure.fix(101325)

    stream[0].assert_electroneutrality(tol=1e-2)

    temp=0
    for ion in feed_mass_frac.keys():
        temp+= feed_mass_frac[ion]/stream[0].dens_mass_comp[ion]
        print(f" X/rho {ion}:",value(temp))
    temp+=H2O_mass_frac/stream[0].dens_mass_comp['H2O']
    print(f" X/rho H2O:", value(temp))

    print(value(1/temp))
    print(H2O_mass_frac)
    stream[0].mole_frac_phase_comp
    #
    stream[0].flow_mass_phase_comp
    #
    stream[0].molality_comp
    stream[0].pressure_osm
    stream[0].dens_mass_phase
    # stream[0].conc_mol_phase_comp
    # stream[0].flow_vol
    iscale.set_scaling_factor(stream[0].dens_mass_comp['Na_+'], 1e-3)
    iscale.set_scaling_factor(stream[0].dens_mass_comp['Ca_2+'], 1e-3)
    iscale.set_scaling_factor(stream[0].dens_mass_comp['Mg_2+'], 1e-3)
    iscale.set_scaling_factor(stream[0].dens_mass_comp['Cl_-'], 1)
    iscale.set_scaling_factor(stream[0].dens_mass_comp['SO4_2-'], 1e-3)

    iscale.set_scaling_factor(stream[0].mass_frac_phase_comp['Liq', 'Na_+'], 1e6)
    iscale.set_scaling_factor(stream[0].mass_frac_phase_comp['Liq', 'Ca_2+'], 1e6)
    iscale.set_scaling_factor(stream[0].mass_frac_phase_comp['Liq', 'Mg_2+'], 1e6)
    iscale.set_scaling_factor(stream[0].mass_frac_phase_comp['Liq', 'Cl_-'], 1e6)
    iscale.set_scaling_factor(stream[0].mass_frac_phase_comp['Liq', 'SO4_2-'], 1e6)
    iscale.calculate_scaling_factors(m.fs)

    assert_units_consistent(m)

    check_dof(m, fail_flag=True)

    stream.initialize()

    results = solver.solve(m)
    assert_optimal_termination(results)

    # assert value(stream[0].flow_vol_phase['Liq']) == pytest.approx(1.91187e-5,  rel=1e-3)
    assert value(stream[0].flow_mol_phase_comp['Liq', 'H2O']) == pytest.approx(53.59256,  rel=1e-3)
    assert value(stream[0].flow_mol_phase_comp['Liq', 'Na_+']) == pytest.approx(0.4836,  rel=1e-3)
    assert value(stream[0].flow_mol_phase_comp['Liq', 'Ca_2+']) == pytest.approx(0.00955,  rel=1e-3)
    assert value(stream[0].flow_mol_phase_comp['Liq', 'Mg_2+']) == pytest.approx(0.05808,  rel=1e-3)
    assert value(stream[0].flow_mol_phase_comp['Liq', 'Cl_-']) == pytest.approx(0.58,  rel=1e-3)
    assert value(stream[0].flow_mol_phase_comp['Liq', 'SO4_2-']) == pytest.approx(0.02225,  rel=1e-3)
    #assert value(stream[0].dens_mass_phase['Liq']) == pytest.approx(1015.89,  rel=1e-3) #TODO revisit after solution density finalized
