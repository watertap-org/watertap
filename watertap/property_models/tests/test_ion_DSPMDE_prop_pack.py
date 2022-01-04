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
from pyomo.environ import ConcreteModel, assert_optimal_termination, value, Set, Param, Var
from idaes.core import FlowsheetBlock
import idaes.core.util.scaling as iscale
from pyomo.util.check_units import assert_units_consistent
from watertap.property_models.ion_DSPMDE_prop_pack import DSPMDEParameterBlock, DSPMDEStateBlock

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