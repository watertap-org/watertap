###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################

import pytest
from pyomo.environ import (ConcreteModel,
                           Constraint,
                           TerminationCondition,
                           SolverStatus,
                           value,
                           Var,
                           units as pyunits,
                           assert_optimal_termination)
from pyomo.network import Port
from idaes.core import (FlowsheetBlock,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        ControlVolume0DBlock)
import watertap.property_models.ion_DSPMDE_prop_pack as props
from watertap.unit_models.nanofiltration_DSPMDE_0D import NanofiltrationDSPMDE0D
from watertap.core.util.initialization import check_dof

from idaes.core.util import get_solver
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_variables,
                                              number_total_constraints,
                                              number_unused_variables,
                                              unused_variables_set)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (calculate_scaling_factors,
                                     unscaled_variables_generator,
                                     unscaled_constraints_generator,
                                     badly_scaled_var_generator)
# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={'dynamic': False})
    m.fs.properties = props.DSPMDEParameterBlock(default={
            "solute_list": ["Ca_2+", "SO4_2-", "Na_+", "Cl_-", "Mg_2+"]})
    m.fs.unit = NanofiltrationDSPMDE0D(default={'property_package': m.fs.properties})


    # check unit config arguments
    assert len(m.fs.unit.config) == 8

    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert m.fs.unit.config.material_balance_type == \
           MaterialBalanceType.useDefault
    assert m.fs.unit.config.momentum_balance_type == \
           MomentumBalanceType.pressureTotal
    assert not m.fs.unit.config.has_pressure_change
    assert m.fs.unit.config.property_package is \
           m.fs.properties

class TestNanoFiltration():
    @pytest.fixture(scope="class")
    def NF_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = props.DSPMDEParameterBlock(default={
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
            "charge": {"Na_+": 1,
                       "Ca_2+": 2,
                       "Mg_2+": 2,
                       "Cl_-": -1,
                       "SO4_2-": -2},
            })

        m.fs.unit = NanofiltrationDSPMDE0D(default={"property_package": m.fs.properties})

        mass_flow_in = 1 * pyunits.kg / pyunits.s
        feed_mass_frac = {'Na_+': 11122e-6,
                          'Ca_2+': 382e-6,
                          'Mg_2+': 1394e-6,
                          'SO4_2-': 2136e-6,
                          'Cl_-': 20300e-6}
        for ion, x in feed_mass_frac.items():
            mol_comp_flow = x * pyunits.kg / pyunits.kg * mass_flow_in / m.fs.properties.mw_comp[ion]
            m.fs.unit.inlet.flow_mol_phase_comp[0, 'Liq', ion].fix(mol_comp_flow)
            # m.fs.unit.permeate.flow_mol_phase_comp[0, 'Liq', ion].fix(0)  # TODO: remove later- temporary to eliminate dof

        H2O_mass_frac = 1 - sum(x for x in feed_mass_frac.values())
        H2O_mol_comp_flow = H2O_mass_frac * pyunits.kg / pyunits.kg * mass_flow_in / \
                            m.fs.unit.feed_side.properties_in[0].mw_comp['H2O']

        m.fs.unit.inlet.flow_mol_phase_comp[0, 'Liq', 'H2O'].fix(H2O_mol_comp_flow)
        m.fs.unit.inlet.temperature[0].fix(298.15)
        m.fs.unit.inlet.pressure[0].fix(5e5)

        m.fs.unit.radius_pore.fix(0.5e-9)
        m.fs.unit.membrane_thickness_effective.fix()
        m.fs.unit.membrane_charge_density.fix()
        m.fs.unit.dielectric_constant_pore.fix(41.3)
        m.fs.unit.mixed_permeate[0].pressure.fix(101325)
        m.fs.unit.area.fix(50)
        m.fs.unit.recovery_vol_phase[0, 'Liq'].fix(0.2)
        m.fs.unit.Kf_comp.fix()

        return m

    @pytest.mark.unit
    def test_build(self, NF_frame):
        m = NF_frame

        # test ports and variables
        port_lst = ['inlet', 'retentate', 'permeate']
        port_vars_lst = ['flow_mol_phase_comp', 'pressure', 'temperature']
        for port_str in port_lst:
            assert hasattr(m.fs.unit, port_str)
            port = getattr(m.fs.unit, port_str)
            assert len(port.vars) == 3
            assert isinstance(port, Port)
            for var_str in port_vars_lst:
                assert hasattr(port, var_str)
                var = getattr(port, var_str)
                assert isinstance(var, Var)

        # test unit objects (including parameters, variables, and constraints)
        unit_objs_lst = []
        for obj_str in unit_objs_lst:
            assert hasattr(m.fs.unit, obj_str)


        # test state block objects
        assert isinstance(m.fs.unit.feed_side, ControlVolume0DBlock)
        cv_stateblock_lst = ['properties_in', 'properties_out', 'properties_interface']
        # feed side
        for sb_str in cv_stateblock_lst:
            sb = getattr(m.fs.unit.feed_side, sb_str)
            assert isinstance(sb, props.DSPMDEStateBlock)
        # test objects added to control volume
        cv_objs_type_dict = {'eq_feed_interface_isothermal': Constraint}
        for (obj_str, obj_type) in cv_objs_type_dict.items():
            obj = getattr(m.fs.unit.feed_side, obj_str)
            assert isinstance(obj, obj_type)
        # permeate side
        assert isinstance(m.fs.unit.permeate_side, props.DSPMDEStateBlock)
        assert isinstance(m.fs.unit.mixed_permeate, props.DSPMDEStateBlock)
        # membrane
        assert isinstance(m.fs.unit.pore_entrance, props.DSPMDEStateBlock)
        assert isinstance(m.fs.unit.pore_exit, props.DSPMDEStateBlock)

        # test statistics
        assert number_variables(m) == 489
        assert number_total_constraints(m) == 450
        assert number_unused_variables(m) == 4  # temperature and pressure unused on particular stateblocks

    @pytest.mark.unit
    def test_dof(self, NF_frame):
        m = NF_frame
        check_dof(m, fail_flag=True)

    @pytest.mark.unit
    def test_calculate_scaling(self, NF_frame):
        m = NF_frame

        # m.fs.properties.set_default_scaling('flow_mol_phase_comp', 1, index=('Liq', 'H2O'))
        # for ion in m.fs.properties.config.solute_list:
        #     m.fs.properties.set_default_scaling('flow_mol_phase_comp', 1e2, index=('Liq', ion))

        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m.fs.unit))
        assert len(unscaled_var_list) == 0

        badly_scaled_var_lst = list(badly_scaled_var_generator(m, include_fixed=True))
        assert len(unscaled_var_list) == 0

        # # check that all constraints have been scaled #TODO: revisit to see if constraints need transformation
        unscaled_constraint_list = list(unscaled_constraints_generator(m.fs.unit))
        [print(i) for i in unscaled_constraint_list]
        assert len(unscaled_constraint_list) == 0


    @pytest.mark.component
    def test_initialize(self, NF_frame):
        m= NF_frame
        # initialization_tester(m)
        m.fs.unit.initialize(fail_on_warning=True)

    @pytest.mark.component
    def test_solve(self, NF_frame):
        m = NF_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

    # @pytest.mark.component
    # def test_conservation(self, NF_frame):
    #     m = NF_frame
    #     b = m.fs.unit
    #     comp_lst = ['NaCl', 'H2O']
    #
    #     flow_mass_inlet = sum(
    #         b.feed_side.properties_in[0].flow_mass_phase_comp['Liq', j] for j in comp_lst)
    #     flow_mass_retentate = sum(
    #         b.feed_side.properties_out[0].flow_mass_phase_comp['Liq', j] for j in comp_lst)
    #     flow_mass_permeate = sum(
    #         b.properties_permeate[0].flow_mass_phase_comp['Liq', j] for j in comp_lst)
    #
    #     assert (abs(value(flow_mass_inlet - flow_mass_retentate - flow_mass_permeate
    #                       )) <= 1e-6)
    #
    #     assert (abs(value(
    #         flow_mass_inlet * b.feed_side.properties_in[0].enth_mass_phase['Liq']
    #         - flow_mass_retentate * b.feed_side.properties_out[0].enth_mass_phase['Liq']
    #         - flow_mass_permeate * b.properties_permeate[0].enth_mass_phase['Liq']
    #     )) <= 1e-6)
    #
    # @pytest.mark.component
    # def test_solution(self, NF_frame):
    #     m = NF_frame
    #     assert (pytest.approx(1.079e-2, rel=1e-3) ==
    #             value(m.fs.unit.flux_mass_phase_comp_avg[0, 'Liq', 'H2O']))
    #     assert (pytest.approx(3.435e-4, rel=1e-3) ==
    #             value(m.fs.unit.flux_mass_phase_comp_avg[0, 'Liq', 'NaCl']))
    #     assert (pytest.approx(0.5396, rel=1e-3) ==
    #             value(m.fs.unit.properties_permeate[0].flow_mass_phase_comp['Liq', 'H2O']))
    #     assert (pytest.approx(1.717e-2, rel=1e-3) ==
    #             value(m.fs.unit.properties_permeate[0].flow_mass_phase_comp['Liq', 'NaCl']))
