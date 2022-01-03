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
from pyomo.environ import ConcreteModel, assert_optimal_termination, value
from idaes.core import FlowsheetBlock
import idaes.core.util.scaling as iscale
from pyomo.util.check_units import assert_units_consistent
import watertap.property_models.ion_DSPMDE_prop_pack as props

from watertap.core.util.initialization import check_dof
from idaes.core.util import get_solver

solver = get_solver()
# -----------------------------------------------------------------------------

@pytest.mark.component
def test_property_ions():
    m = ConcreteModel()

    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = props.DSPMDEParameterBlock(default={
        "solute_list": ["A", "B", "C"],
        "diffusivity_data": {("Liq","A"): 1e-9,
                             ("Liq","B"): 1e-10,
                             ("Liq","C"): 1e-7},
        "mw_data": {"H2O": 18e-3,
                    "A": 10e-3,
                    "B": 25e-3,
                    "C": 100e-3},
        "stokes_radius_data": {"A": 1e-9,
                               "B": 1e-9,
                               "C": 1e-9},
        "density_data": {"H2O": 1000,
                         "A": 1200,
                         "B": 1100,
                         "C": 1010}
    })
    m.fs.stream = m.fs.properties.build_state_block([0], default={'defined_state': True})
    m.fs.stream[0].flow_mol_phase_comp['Liq', 'A'].fix(0.000407)
    m.fs.stream[0].flow_mol_phase_comp['Liq', 'B'].fix(0.010479)
    m.fs.stream[0].flow_mol_phase_comp['Liq', 'C'].fix(0.010479)
    m.fs.stream[0].flow_mol_phase_comp['Liq', 'H2O'].fix(0.979046)
    m.fs.stream[0].molality_comp["B"]
    print(m.fs.stream[0].radius_stokes_comp["B"].value)
    assert False
    # # specify
    # m.fs.stream[0].flow_mol_phase_comp['Liq', 'Na_+'].fix(0.008845)
    # m.fs.stream[0].flow_mol_phase_comp['Liq', 'Ca_2+'].fix(0.000174)
    # m.fs.stream[0].flow_mol_phase_comp['Liq', 'Mg_2+'].fix(0.001049)
    # m.fs.stream[0].flow_mol_phase_comp['Liq', 'SO4_2-'].fix(0.000407)
    # m.fs.stream[0].flow_mol_phase_comp['Liq', 'Cl_-'].fix(0.010479)
    # m.fs.stream[0].flow_mol_phase_comp['Liq', 'H2O'].fix(0.979046)
    # m.fs.stream[0].temperature.fix(273.15 + 25)
    # m.fs.stream[0].pressure.fix(101325)
    #
    #
    # # # scaling
    # iscale.calculate_scaling_factors(m.fs)
    #
    # # checking state block
    # assert_units_consistent(m)
    #
    # # check dof = 0
    # check_dof(m, fail_flag=True)
    #
    # # initialize
    # m.fs.stream.initialize()
    #
    # # check solve
    # results = solver.solve(m)
    # assert_optimal_termination(results)
    #
    # # # check values
    # assert value(m.fs.stream[0].mole_frac_phase_comp['Liq', 'Na_+']) == pytest.approx(8.845e-3,  rel=1e-3)
    # assert value(m.fs.stream[0].mole_frac_phase_comp['Liq', 'Ca_2+']) == pytest.approx(1.74e-4,  rel=1e-3)
    # assert value(m.fs.stream[0].mole_frac_phase_comp['Liq', 'Cl_-']) == pytest.approx(1.048e-2,  rel=1e-3)
    # assert value(m.fs.stream[0].mole_frac_phase_comp['Liq', 'H2O']) == pytest.approx(0.9790,  rel=1e-3)
    # assert value(m.fs.stream[0].mole_frac_phase_comp['Liq', 'Mg_2+']) == pytest.approx(1.049e-3,  rel=1e-3)
    # assert value(m.fs.stream[0].mole_frac_phase_comp['Liq', 'SO4_2-']) == pytest.approx(4.07e-4,  rel=1e-3)
    #

