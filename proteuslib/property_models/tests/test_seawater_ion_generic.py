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
from proteuslib.property_models.tests.property_test_harness import \
    (PropertyTestHarness, PropertyRegressionTest)
from pyomo.environ import ConcreteModel, value
from idaes.core import FlowsheetBlock
import idaes.core.util.scaling as iscale
from pyomo.util.check_units import assert_units_consistent
from proteuslib.flowsheets.full_treatment_train.util import solve_with_user_scaling, check_dof
from idaes.generic_models.properties.core.generic.generic_property import GenericParameterBlock
from proteuslib.property_models.seawater_ion_generic import configuration

# -----------------------------------------------------------------------------

@pytest.mark.component
def test_property_seawater_ions():
    m = ConcreteModel()

    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = GenericParameterBlock(default=configuration)
    m.fs.stream = m.fs.properties.build_state_block([0], default={'defined_state': True})

    # specify
    feed_flow_mass = 1
    feed_mass_frac = {'Na': 11122e-6,
                      'Ca': 382e-6,
                      'Mg': 1394e-6,
                      'SO4': 2136e-6,
                      'Cl': 20300e-6}
    m.fs.stream[0].flow_mol_phase_comp['Liq', 'Na_+'].fix(0.008845)
    m.fs.stream[0].flow_mol_phase_comp['Liq', 'Ca_2+'].fix(0.000174)
    m.fs.stream[0].flow_mol_phase_comp['Liq', 'Mg_2+'].fix(0.001049)
    m.fs.stream[0].flow_mol_phase_comp['Liq', 'SO4_2-'].fix(0.000407)
    m.fs.stream[0].flow_mol_phase_comp['Liq', 'Cl_-'].fix(0.010479)
    m.fs.stream[0].flow_mol_phase_comp['Liq', 'H2O'].fix(0.979046)
    m.fs.stream[0].flow_vol
    m.fs.stream[0].temperature.fix(273.15 + 25)
    m.fs.stream[0].pressure.fix(101325)

    m.fs.stream[0].flow_mol_phase_comp
    m.fs.stream[0].mw_comp
    m.fs.stream[0].flow_mass_phase_comp
    # m.fs.stream[0].mass_frac_comp # TODO: mass_frac_comp and mass_frac_phase_comp need to be added to generic property model

    # # scaling
    # m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
    # m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'Na'))
    # m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e4, index=('Liq', 'Ca'))
    # m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e3, index=('Liq', 'Mg'))
    # m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e3, index=('Liq', 'SO4'))
    # m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'Cl'))
    iscale.calculate_scaling_factors(m.fs)
    #
    # # checking state block
    assert_units_consistent(m)
    #TODO: replace check_dof
    check_dof(m)
    #
    # # initialize
    m.fs.stream.initialize()
    #
    # # solve
    #TODO: replace solve_with_user_scaling
    solve_with_user_scaling(m)
    #
    # # check values
    #TODO: mass_frac_phase_comp won't work until generic prop modified

    # assert value(m.fs.stream[0].mass_frac_phase_comp['Liq', 'H2O']) == pytest.approx(0.9647, rel=1e-3)
    #TODO: flow_mol_p_c is being fixed; change to assert flow_mass_phase_comp, flow_vol
    # assert value(m.fs.stream[0].flow_mol_phase_comp['Liq', 'Na_+']) == pytest.approx(0.4838, rel=1e-3)
    # assert value(m.fs.stream[0].flow_mol_phase_comp['Liq', 'Ca_2+']) == pytest.approx(9.531e-3, rel=1e-3)

# class TestPropertySeawaterIons(PropertyTestHarness):
#     def configure(self):
#         self.prop_pack = GenericParameterBlock
#         self.param_args = {'default': configuration}
#         self.scaling_args = {('flow_mass_phase_comp', ('Liq', 'H2O')): 1,
#                              ('flow_mass_phase_comp', ('Liq', 'Na')): 1e2,
#                              ('flow_mass_phase_comp', ('Liq', 'Ca')): 1e4,
#                              ('flow_mass_phase_comp', ('Liq', 'Mg')): 1e3,
#                              ('flow_mass_phase_comp', ('Liq', 'SO4')): 1e3,
#                              ('flow_mass_phase_comp', ('Liq', 'Cl')): 1e2}
#         self.stateblock_statistics = {'number_variables': 27,
#                                       'number_total_constraints': 19,
#                                       'number_unused_variables': 2,  # pressure and temperature are unused
#                                       'default_degrees_of_freedom': 6}  # 8 state vars, pressure and temp are not active
#         self.default_solution = {('mass_frac_phase_comp', ('Liq', 'H2O')): 0.9647,
#                                  ('flow_vol', None): 1e-3,
#                                  ('flow_mol_phase_comp', ('Liq', 'Na')): 0.4838,
#                                  ('flow_mol_phase_comp', ('Liq', 'Ca')): 9.531e-3,
#                                  ('enth_flow', None): 1.05e5}
#
#
# class TestPropertySeawaterSalts(PropertyTestHarness):
#     def configure(self):
#         self.prop_pack = property_seawater_salts.PropParameterBlock
#         self.param_args = {}
#         self.scaling_args = {('flow_mass_phase_comp', ('Liq', 'H2O')): 1,
#                              ('flow_mass_phase_comp', ('Liq', 'NaCl')): 1e2,
#                              ('flow_mass_phase_comp', ('Liq', 'CaSO4')): 1e3,
#                              ('flow_mass_phase_comp', ('Liq', 'MgSO4')): 1e3,
#                              ('flow_mass_phase_comp', ('Liq', 'MgCl2')): 1e3}
#         self.stateblock_statistics = {'number_variables': 18,
#                                       'number_total_constraints': 11,
#                                       'number_unused_variables': 2,  # pressure and temperature are unused
#                                       'default_degrees_of_freedom': 5}  # 7 state vars, pressure and temp are not active
#         self.default_solution = {('mass_frac_phase_comp', ('Liq', 'H2O')): 0.9647,
#                                  ('flow_vol', None): 1e-3,
#                                  ('flow_mol_phase_comp', ('Liq', 'NaCl')): 0.4838,
#                                  ('flow_mol_phase_comp', ('Liq', 'CaSO4')): 9.531e-3,
#                                  ('flow_mol_phase_comp', ('Liq', 'MgSO4')): 1.270e-2,
#                                  ('flow_mol_phase_comp', ('Liq', 'MgCl2')): 4.465e-2,
#                                  ('enth_flow', None): 1.05e5}
