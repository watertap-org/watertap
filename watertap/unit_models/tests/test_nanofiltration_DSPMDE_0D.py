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
from pyomo.environ import (
    ConcreteModel,
    Constraint,
    TerminationCondition,
    SolverStatus,
    value,
    Var,
    units as pyunits,
    assert_optimal_termination,
)
from pyomo.util.check_units import assert_units_consistent
from pyomo.network import Port
from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    ControlVolume0DBlock,
)
from watertap.property_models.ion_DSPMDE_prop_pack import (
    DSPMDEParameterBlock,
    ActivityCoefficientModel,
    DensityCalculation,
    DSPMDEStateBlock,
)
from watertap.unit_models.nanofiltration_DSPMDE_0D import (
    NanofiltrationDSPMDE0D,
    MassTransferCoefficient,
    ConcentrationPolarizationType,
)
from watertap.core.util.initialization import check_dof

from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
    unused_variables_set,
)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    unscaled_constraints_generator,
    constraints_with_scale_factor_generator,
    badly_scaled_var_generator,
)
import idaes.logger as idaeslog

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


@pytest.mark.unit
def test_config_with_CP():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = DSPMDEParameterBlock(
        solute_list=["Ca_2+", "SO4_2-", "Na_+", "Cl_-", "Mg_2+"],
        charge={"Ca_2+": 2, "SO4_2-": -2, "Na_+": 1, "Cl_-": -1, "Mg_2+": 2},
    )
    m.fs.unit = NanofiltrationDSPMDE0D(property_package=m.fs.properties)

    # check unit config arguments
    assert len(m.fs.unit.config) == 9

    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
    assert m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
    assert not m.fs.unit.config.has_pressure_change
    assert not hasattr(m.fs.unit.config, "energy_balance_type")
    assert m.fs.unit.config.property_package is m.fs.properties
    assert (
        m.fs.unit.config.property_package.config.activity_coefficient_model
        == ActivityCoefficientModel.ideal
    )
    assert (
        m.fs.unit.config.property_package.config.density_calculation
        == DensityCalculation.constant
    )
    assert (
        m.fs.unit.config.mass_transfer_coefficient
        == MassTransferCoefficient.spiral_wound
    )
    assert (
        m.fs.unit.config.concentration_polarization_type
        == ConcentrationPolarizationType.calculated
    )
    assert hasattr(m.fs.unit.config.property_package, "solute_set")
    assert hasattr(m.fs.unit.config.property_package, "ion_set")
    assert hasattr(m.fs.unit.config.property_package, "cation_set")
    assert hasattr(m.fs.unit.config.property_package, "anion_set")
    assert_units_consistent(m)


@pytest.mark.unit
def test_config_without_CP():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = DSPMDEParameterBlock(
        solute_list=["Ca_2+", "SO4_2-", "Na_+", "Cl_-", "Mg_2+"],
        charge={"Ca_2+": 2, "SO4_2-": -2, "Na_+": 1, "Cl_-": -1, "Mg_2+": 2},
    )
    m.fs.unit = NanofiltrationDSPMDE0D(
        property_package=m.fs.properties,
        mass_transfer_coefficient=MassTransferCoefficient.none,
        concentration_polarization_type=ConcentrationPolarizationType.none,
    )

    # check unit config arguments
    assert len(m.fs.unit.config) == 9

    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
    assert m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
    assert not m.fs.unit.config.has_pressure_change
    assert not hasattr(m.fs.unit.config, "energy_balance_type")
    assert m.fs.unit.config.property_package is m.fs.properties
    assert (
        m.fs.unit.config.property_package.config.activity_coefficient_model
        == ActivityCoefficientModel.ideal
    )
    assert (
        m.fs.unit.config.property_package.config.density_calculation
        == DensityCalculation.constant
    )
    assert m.fs.unit.config.mass_transfer_coefficient == MassTransferCoefficient.none
    assert (
        m.fs.unit.config.concentration_polarization_type
        == ConcentrationPolarizationType.none
    )
    assert hasattr(m.fs.unit.config.property_package, "solute_set")
    assert hasattr(m.fs.unit.config.property_package, "ion_set")
    assert hasattr(m.fs.unit.config.property_package, "cation_set")
    assert hasattr(m.fs.unit.config.property_package, "anion_set")
    assert_units_consistent(m)


class TestNanoFiltration_with_CP_5ions:
    @pytest.fixture(scope="class")
    def NF_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = DSPMDEParameterBlock(
            solute_list=["Ca_2+", "SO4_2-", "Mg_2+", "Na_+", "Cl_-"],
            diffusivity_data={
                ("Liq", "Ca_2+"): 9.2e-10,
                ("Liq", "SO4_2-"): 1.06e-09,
                ("Liq", "Mg_2+"): 7.06e-10,
                ("Liq", "Na_+"): 1.33e-09,
                ("Liq", "Cl_-"): 2.03e-09,
            },
            mw_data={
                "H2O": 0.018,
                "Ca_2+": 0.04,
                "Mg_2+": 0.024,
                "SO4_2-": 0.096,
                "Na_+": 0.023,
                "Cl_-": 0.035,
            },
            stokes_radius_data={
                "Ca_2+": 3.09e-10,
                "Mg_2+": 3.47e-10,
                "SO4_2-": 2.3e-10,
                "Cl_-": 1.21e-10,
                "Na_+": 1.84e-10,
            },
            charge={"Ca_2+": 2, "Mg_2+": 2, "SO4_2-": -2, "Na_+": 1, "Cl_-": -1},
            activity_coefficient_model=ActivityCoefficientModel.davies,
            density_calculation=DensityCalculation.constant,
        )

        m.fs.unit = NanofiltrationDSPMDE0D(property_package=m.fs.properties)
        b = m.fs.unit
        mass_flow_in = 1 * pyunits.kg / pyunits.s
        feed_mass_frac = {
            "Ca_2+": 382e-6,
            "Mg_2+": 1394e-6,
            "SO4_2-": 2136e-6,
            "Cl_-": 20101.6e-6,
            "Na_+": 11122e-6,
        }

        # Fix mole flow rates of each ion and water
        for ion, x in feed_mass_frac.items():
            mol_comp_flow = (
                x
                * pyunits.kg
                / pyunits.kg
                * mass_flow_in
                / m.fs.unit.feed_side.properties_in[0].mw_comp[ion]
            )
            m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", ion].fix(mol_comp_flow)
        H2O_mass_frac = 1 - sum(x for x in feed_mass_frac.values())
        H2O_mol_comp_flow = (
            H2O_mass_frac
            * pyunits.kg
            / pyunits.kg
            * mass_flow_in
            / m.fs.unit.feed_side.properties_in[0].mw_comp["H2O"]
        )
        m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix(H2O_mol_comp_flow)

        # Use assert electroneutrality method from property model to ensure the ion concentrations provided
        # obey electroneutrality condition
        m.fs.unit.feed_side.properties_in[0].assert_electroneutrality(
            defined_state=True,
            adjust_by_ion="Cl_-",
            get_property="mass_frac_phase_comp",
        )

        # Fix other inlet state variables
        m.fs.unit.inlet.temperature[0].fix(298.15)
        m.fs.unit.inlet.pressure[0].fix(4e5)

        # Fix the membrane variables that are usually fixed for the DSPM-DE model
        m.fs.unit.radius_pore.fix(0.5e-9)
        m.fs.unit.membrane_thickness_effective.fix(1.33e-6)
        m.fs.unit.membrane_charge_density.fix(-27)
        m.fs.unit.dielectric_constant_pore.fix(41.3)

        # Fix final permeate pressure to be ~atmospheric
        m.fs.unit.mixed_permeate[0].pressure.fix(101325)

        m.fs.unit.spacer_porosity.fix(0.85)
        m.fs.unit.channel_height.fix(5e-4)
        m.fs.unit.velocity[0, 0].fix(0.25)
        m.fs.unit.area.fix(50)
        # Fix additional variables for calculating mass transfer coefficient with spiral wound correlation
        m.fs.unit.spacer_mixing_efficiency.fix()
        m.fs.unit.spacer_mixing_length.fix()
        return m

    @pytest.mark.unit
    def test_build(self, NF_frame):
        m = NF_frame

        # test ports and variables
        port_lst = ["inlet", "retentate", "permeate"]
        port_vars_lst = ["flow_mol_phase_comp", "pressure", "temperature"]
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
        cv_stateblock_lst = [
            "properties_in",
            "properties_out",
            "properties_interface",
        ]
        # feed side
        for sb_str in cv_stateblock_lst:
            sb = getattr(m.fs.unit.feed_side, sb_str)
            assert isinstance(sb, DSPMDEStateBlock)
        # test objects added to control volume
        cv_objs_type_dict = {"eq_feed_interface_isothermal": Constraint}
        for (obj_str, obj_type) in cv_objs_type_dict.items():
            obj = getattr(m.fs.unit.feed_side, obj_str)
            assert isinstance(obj, obj_type)
        # permeate side
        assert isinstance(m.fs.unit.permeate_side, DSPMDEStateBlock)
        assert isinstance(m.fs.unit.mixed_permeate, DSPMDEStateBlock)
        # membrane
        assert isinstance(m.fs.unit.pore_entrance, DSPMDEStateBlock)
        assert isinstance(m.fs.unit.pore_exit, DSPMDEStateBlock)

        # test statistics
        assert number_variables(m) == 561
        assert number_total_constraints(m) == 526
        assert number_unused_variables(m) == 11

    @pytest.mark.unit
    def test_dof(self, NF_frame):
        m = NF_frame
        check_dof(m, fail_flag=True)

    @pytest.mark.unit
    def test_calculate_scaling(self, NF_frame):
        m = NF_frame

        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e2, index=("Liq", "Ca_2+")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e2, index=("Liq", "SO4_2-")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e2, index=("Liq", "Mg_2+")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e2, index=("Liq", "Cl_-")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e2, index=("Liq", "Na_+")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e0, index=("Liq", "H2O")
        )

        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m.fs.unit))
        assert len(unscaled_var_list) == 0

        badly_scaled_var_lst = list(badly_scaled_var_generator(m.fs.unit))
        assert len(badly_scaled_var_lst) == 0

        # not all constraints have scaling factor so skipping the check for unscaled constraints
        unscaled_con_lst = list(unscaled_constraints_generator(m))
        assert len(unscaled_con_lst) == 0

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_initialize(self, NF_frame):
        m = NF_frame
        # Using the 'initialize' function so that I can view the logs on failure
        # m.fs.unit.initialize(optarg=solver.options, outlvl=idaeslog.DEBUG)
        initialization_tester(m)

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_solve(self, NF_frame):
        m = NF_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_conservation(self, NF_frame):
        m = NF_frame
        b = m.fs.unit
        comp_lst = m.fs.properties.solute_set

        flow_mass_inlet = sum(
            b.feed_side.properties_in[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_retentate = sum(
            b.feed_side.properties_out[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_permeate = sum(
            b.mixed_permeate[0].flow_mass_phase_comp["Liq", j] for j in comp_lst
        )

        assert (
            abs(value(flow_mass_inlet - flow_mass_retentate - flow_mass_permeate))
            <= 1e-6
        )

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_solution(self, NF_frame):
        m = NF_frame

        mole_flux_dict = {
            "Na_+": 0.00029831,
            "Cl_-": 0.00031137,
            "Ca_2+": 1.9289e-06,
            "SO4_2-": 1.95859e-06,
            "Mg_2+": 6.56667e-06,
            "H2O": 0.032671,
        }
        for j, val in mole_flux_dict.items():
            assert pytest.approx(val, rel=1e-3) == value(
                m.fs.unit.flux_mol_phase_comp_avg[0, "Liq", j]
            )

        # TODO: subsequently focus on the segment below during the validation and refinement phase
        intrinsic_rejection_dict = {
            "Na_+": -0.0265173,
            "Cl_-": 0.11538899,
            "Ca_2+": 0.6905496,
            "SO4_2-": 0.8687813,
            "Mg_2+": 0.83309579,
        }
        for j, val in intrinsic_rejection_dict.items():
            assert pytest.approx(val, rel=1e-3) == value(
                m.fs.unit.rejection_intrinsic_phase_comp[0, "Liq", j]
            )

    @pytest.mark.unit
    def test_report(self, NF_frame):
        NF_frame.fs.unit.report()


class TestNanoFiltration_without_CP_5ions:
    @pytest.fixture(scope="class")
    def NF_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = DSPMDEParameterBlock(
            solute_list=["Ca_2+", "SO4_2-", "Mg_2+", "Na_+", "Cl_-"],
            diffusivity_data={
                ("Liq", "Ca_2+"): 9.2e-10,
                ("Liq", "SO4_2-"): 1.06e-09,
                ("Liq", "Mg_2+"): 7.06e-10,
                ("Liq", "Na_+"): 1.33e-09,
                ("Liq", "Cl_-"): 2.03e-09,
            },
            mw_data={
                "H2O": 0.018,
                "Ca_2+": 0.04,
                "Mg_2+": 0.024,
                "SO4_2-": 0.096,
                "Na_+": 0.023,
                "Cl_-": 0.035,
            },
            stokes_radius_data={
                "Ca_2+": 3.09e-10,
                "Mg_2+": 3.47e-10,
                "SO4_2-": 2.3e-10,
                "Cl_-": 1.21e-10,
                "Na_+": 1.84e-10,
            },
            charge={"Ca_2+": 2, "Mg_2+": 2, "SO4_2-": -2, "Na_+": 1, "Cl_-": -1},
            activity_coefficient_model=ActivityCoefficientModel.davies,
            density_calculation=DensityCalculation.constant,
        )

        m.fs.unit = NanofiltrationDSPMDE0D(
            property_package=m.fs.properties,
            mass_transfer_coefficient=MassTransferCoefficient.none,
            concentration_polarization_type=ConcentrationPolarizationType.none,
        )
        b = m.fs.unit
        mass_flow_in = 1 * pyunits.kg / pyunits.s
        feed_mass_frac = {
            "Ca_2+": 382e-6,
            "Mg_2+": 1394e-6,
            "SO4_2-": 2136e-6,
            "Cl_-": 20101.6e-6,
            "Na_+": 11122e-6,
        }

        # Fix mole flow rates of each ion and water
        for ion, x in feed_mass_frac.items():
            mol_comp_flow = (
                x
                * pyunits.kg
                / pyunits.kg
                * mass_flow_in
                / m.fs.unit.feed_side.properties_in[0].mw_comp[ion]
            )
            m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", ion].fix(mol_comp_flow)
        H2O_mass_frac = 1 - sum(x for x in feed_mass_frac.values())
        H2O_mol_comp_flow = (
            H2O_mass_frac
            * pyunits.kg
            / pyunits.kg
            * mass_flow_in
            / m.fs.unit.feed_side.properties_in[0].mw_comp["H2O"]
        )
        m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix(H2O_mol_comp_flow)

        # Use assert electroneutrality method from property model to ensure the ion concentrations provided
        # obey electroneutrality condition
        m.fs.unit.feed_side.properties_in[0].assert_electroneutrality(
            defined_state=True,
            adjust_by_ion="Cl_-",
            get_property="mass_frac_phase_comp",
        )

        # Fix other inlet state variables
        m.fs.unit.inlet.temperature[0].fix(298.15)
        m.fs.unit.inlet.pressure[0].fix(4e5)

        # Fix the membrane variables that are usually fixed for the DSPM-DE model
        m.fs.unit.radius_pore.fix(0.5e-9)
        m.fs.unit.membrane_thickness_effective.fix(1.33e-6)
        m.fs.unit.membrane_charge_density.fix(-27)
        m.fs.unit.dielectric_constant_pore.fix(41.3)

        # Fix final permeate pressure to be ~atmospheric
        m.fs.unit.mixed_permeate[0].pressure.fix(101325)

        m.fs.unit.spacer_porosity.fix(0.85)
        m.fs.unit.channel_height.fix(5e-4)
        m.fs.unit.velocity[0, 0].fix(0.25)
        m.fs.unit.area.fix(50)

        return m

    @pytest.mark.unit
    def test_build(self, NF_frame):
        m = NF_frame

        # test ports and variables
        port_lst = ["inlet", "retentate", "permeate"]
        port_vars_lst = ["flow_mol_phase_comp", "pressure", "temperature"]
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
        cv_stateblock_lst = [
            "properties_in",
            "properties_out",
            "properties_interface",
        ]
        # feed side
        for sb_str in cv_stateblock_lst:
            sb = getattr(m.fs.unit.feed_side, sb_str)
            assert isinstance(sb, DSPMDEStateBlock)
        # test objects added to control volume
        cv_objs_type_dict = {"eq_feed_interface_isothermal": Constraint}
        for (obj_str, obj_type) in cv_objs_type_dict.items():
            obj = getattr(m.fs.unit.feed_side, obj_str)
            assert isinstance(obj, obj_type)
        # permeate side
        assert isinstance(m.fs.unit.permeate_side, DSPMDEStateBlock)
        assert isinstance(m.fs.unit.mixed_permeate, DSPMDEStateBlock)
        # membrane
        assert isinstance(m.fs.unit.pore_entrance, DSPMDEStateBlock)
        assert isinstance(m.fs.unit.pore_exit, DSPMDEStateBlock)

        # test statistics
        assert number_variables(m) == 527
        assert number_total_constraints(m) == 494
        assert number_unused_variables(m) == 11

    @pytest.mark.unit
    def test_dof(self, NF_frame):
        m = NF_frame
        check_dof(m, fail_flag=True)

    @pytest.mark.unit
    def test_calculate_scaling(self, NF_frame):
        m = NF_frame

        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e2, index=("Liq", "Ca_2+")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e2, index=("Liq", "SO4_2-")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e2, index=("Liq", "Mg_2+")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e2, index=("Liq", "Cl_-")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e2, index=("Liq", "Na_+")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e0, index=("Liq", "H2O")
        )

        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m.fs.unit))
        assert len(unscaled_var_list) == 0

        badly_scaled_var_lst = list(badly_scaled_var_generator(m.fs.unit))
        assert len(badly_scaled_var_lst) == 0

        # not all constraints have scaling factor so skipping the check for unscaled constraints
        unscaled_con_lst = list(unscaled_constraints_generator(m))
        assert len(unscaled_con_lst) == 0

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_initialize(self, NF_frame):
        m = NF_frame
        # Using the 'initialize' function so that I can view the logs on failure
        # m.fs.unit.initialize(optarg=solver.options, outlvl=idaeslog.DEBUG)
        initialization_tester(m)

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_solve(self, NF_frame):
        m = NF_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_conservation(self, NF_frame):
        m = NF_frame
        b = m.fs.unit
        comp_lst = m.fs.properties.solute_set

        flow_mass_inlet = sum(
            b.feed_side.properties_in[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_retentate = sum(
            b.feed_side.properties_out[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_permeate = sum(
            b.mixed_permeate[0].flow_mass_phase_comp["Liq", j] for j in comp_lst
        )

        assert (
            abs(value(flow_mass_inlet - flow_mass_retentate - flow_mass_permeate))
            <= 1e-6
        )

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_solution(self, NF_frame):
        m = NF_frame

        mole_flux_dict = {
            "Na_+": 0.0003574,
            "Cl_-": 0.0003694,
            "Ca_2+": 1.860337e-6,
            "SO4_2-": 1.850918e-6,
            "Mg_2+": 5.972139e-06,
            "H2O": 0.039150,
        }
        for j, val in mole_flux_dict.items():
            assert pytest.approx(val, rel=1e-3) == value(
                m.fs.unit.flux_mol_phase_comp_avg[0, "Liq", j]
            )

        # TODO: subsequently focus on the segment below during the validation and refinement phase
        intrinsic_rejection_dict = {
            "Na_+": -0.017782,
            "Cl_-": 0.1143856,
            "Ca_2+": 0.731772,
            "SO4_2-": 0.8854559,
            "Mg_2+": 0.8584226,
        }
        for j, val in intrinsic_rejection_dict.items():
            assert pytest.approx(val, rel=1e-3) == value(
                m.fs.unit.rejection_intrinsic_phase_comp[0, "Liq", j]
            )

    @pytest.mark.unit
    def test_report(self, NF_frame):
        NF_frame.fs.unit.report()


class TestNanoFiltration_with_CP_2ions:
    @pytest.fixture(scope="class")
    def NF_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = DSPMDEParameterBlock(
            solute_list=["Na_+", "Cl_-"],
            diffusivity_data={("Liq", "Na_+"): 1.33e-09, ("Liq", "Cl_-"): 2.03e-09},
            mw_data={"H2O": 0.018, "Na_+": 0.023, "Cl_-": 0.035},
            stokes_radius_data={"Cl_-": 1.21e-10, "Na_+": 1.84e-10},
            charge={"Na_+": 1, "Cl_-": -1},
            activity_coefficient_model=ActivityCoefficientModel.davies,
            density_calculation=DensityCalculation.constant,
        )

        m.fs.unit = NanofiltrationDSPMDE0D(property_package=m.fs.properties)
        b = m.fs.unit
        mass_flow_in = 1 * pyunits.kg / pyunits.s
        feed_mass_frac = {
            "Cl_-": 20101.6e-6,
            "Na_+": 11122e-6,
        }

        # Fix mole flow rates of each ion and water
        for ion, x in feed_mass_frac.items():
            mol_comp_flow = (
                x
                * pyunits.kg
                / pyunits.kg
                * mass_flow_in
                / m.fs.unit.feed_side.properties_in[0].mw_comp[ion]
            )
            m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", ion].fix(mol_comp_flow)
        H2O_mass_frac = 1 - sum(x for x in feed_mass_frac.values())
        H2O_mol_comp_flow = (
            H2O_mass_frac
            * pyunits.kg
            / pyunits.kg
            * mass_flow_in
            / m.fs.unit.feed_side.properties_in[0].mw_comp["H2O"]
        )
        m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix(H2O_mol_comp_flow)

        # Use assert electroneutrality method from property model to ensure the ion concentrations provided
        # obey electroneutrality condition
        m.fs.unit.feed_side.properties_in[0].assert_electroneutrality(
            defined_state=True,
            adjust_by_ion="Cl_-",
            get_property="mass_frac_phase_comp",
        )

        # Fix other inlet state variables
        m.fs.unit.inlet.temperature[0].fix(298.15)
        m.fs.unit.inlet.pressure[0].fix(4e5)

        # Fix the membrane variables that are usually fixed for the DSPM-DE model
        m.fs.unit.radius_pore.fix(0.5e-9)
        m.fs.unit.membrane_thickness_effective.fix(1.33e-6)
        m.fs.unit.membrane_charge_density.fix(-27)
        m.fs.unit.dielectric_constant_pore.fix(41.3)

        # Fix final permeate pressure to be ~atmospheric
        m.fs.unit.mixed_permeate[0].pressure.fix(101325)

        m.fs.unit.spacer_porosity.fix(0.85)
        m.fs.unit.channel_height.fix(5e-4)
        m.fs.unit.velocity[0, 0].fix(0.25)
        m.fs.unit.area.fix(50)
        # Fix additional variables for calculating mass transfer coefficient with spiral wound correlation
        m.fs.unit.spacer_mixing_efficiency.fix()
        m.fs.unit.spacer_mixing_length.fix()
        return m

    @pytest.mark.unit
    def test_build(self, NF_frame):
        m = NF_frame

        # test ports and variables
        port_lst = ["inlet", "retentate", "permeate"]
        port_vars_lst = ["flow_mol_phase_comp", "pressure", "temperature"]
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
        cv_stateblock_lst = [
            "properties_in",
            "properties_out",
            "properties_interface",
        ]
        # feed side
        for sb_str in cv_stateblock_lst:
            sb = getattr(m.fs.unit.feed_side, sb_str)
            assert isinstance(sb, DSPMDEStateBlock)
        # test objects added to control volume
        cv_objs_type_dict = {"eq_feed_interface_isothermal": Constraint}
        for (obj_str, obj_type) in cv_objs_type_dict.items():
            obj = getattr(m.fs.unit.feed_side, obj_str)
            assert isinstance(obj, obj_type)
        # permeate side
        assert isinstance(m.fs.unit.permeate_side, DSPMDEStateBlock)
        assert isinstance(m.fs.unit.mixed_permeate, DSPMDEStateBlock)
        # membrane
        assert isinstance(m.fs.unit.pore_entrance, DSPMDEStateBlock)
        assert isinstance(m.fs.unit.pore_exit, DSPMDEStateBlock)

        # test statistics
        assert number_variables(m) == 318
        assert number_total_constraints(m) == 286
        assert number_unused_variables(m) == 11

    @pytest.mark.unit
    def test_dof(self, NF_frame):
        m = NF_frame
        check_dof(m, fail_flag=True)

    @pytest.mark.unit
    def test_calculate_scaling(self, NF_frame):
        m = NF_frame

        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e2, index=("Liq", "Cl_-")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e2, index=("Liq", "Na_+")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e0, index=("Liq", "H2O")
        )

        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m.fs.unit))
        assert len(unscaled_var_list) == 0

        badly_scaled_var_lst = list(badly_scaled_var_generator(m.fs.unit))
        assert len(badly_scaled_var_lst) == 0

        # not all constraints have scaling factor so skipping the check for unscaled constraints
        unscaled_con_lst = list(unscaled_constraints_generator(m))
        assert len(unscaled_con_lst) == 0

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_initialize(self, NF_frame):
        m = NF_frame
        # Using the 'initialize' function so that I can view the logs on failure
        # m.fs.unit.initialize(optarg=solver.options, outlvl=idaeslog.DEBUG)
        initialization_tester(m)

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_solve(self, NF_frame):
        m = NF_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_conservation(self, NF_frame):
        m = NF_frame
        b = m.fs.unit
        comp_lst = m.fs.properties.solute_set

        flow_mass_inlet = sum(
            b.feed_side.properties_in[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_retentate = sum(
            b.feed_side.properties_out[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_permeate = sum(
            b.mixed_permeate[0].flow_mass_phase_comp["Liq", j] for j in comp_lst
        )

        assert (
            abs(value(flow_mass_inlet - flow_mass_retentate - flow_mass_permeate))
            <= 1e-6
        )

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_solution(self, NF_frame):
        m = NF_frame

        mole_flux_dict = {
            "Na_+": 0.0010739,
            "Cl_-": 0.0010739,
            "H2O": 0.129305,
        }
        for j, val in mole_flux_dict.items():
            assert pytest.approx(val, rel=1e-3) == value(
                m.fs.unit.flux_mol_phase_comp_avg[0, "Liq", j]
            )

        # TODO: subsequently focus on the segment below during the validation and refinement phase
        intrinsic_rejection_dict = {
            "Na_+": 0.0964865,
            "Cl_-": 0.0964865,
        }
        for j, val in intrinsic_rejection_dict.items():
            assert pytest.approx(val, rel=1e-3) == value(
                m.fs.unit.rejection_intrinsic_phase_comp[0, "Liq", j]
            )

        flow_mol_out_dict = {
            "Na_+": 0.429868,
            "Cl_-": 0.429868,
        }
        for j, val in flow_mol_out_dict.items():
            assert pytest.approx(val, rel=1e-3) == value(
                m.fs.unit.retentate.flow_mol_phase_comp[0, "Liq", j]
            )

    @pytest.mark.unit
    def test_report(self, NF_frame):
        NF_frame.fs.unit.report()


class TestNanoFiltration_without_CP_2ions:
    @pytest.fixture(scope="class")
    def NF_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = DSPMDEParameterBlock(
            solute_list=["Na_+", "Cl_-"],
            diffusivity_data={("Liq", "Na_+"): 1.33e-09, ("Liq", "Cl_-"): 2.03e-09},
            mw_data={"H2O": 0.018, "Na_+": 0.023, "Cl_-": 0.035},
            stokes_radius_data={"Cl_-": 1.21e-10, "Na_+": 1.84e-10},
            charge={"Na_+": 1, "Cl_-": -1},
            activity_coefficient_model=ActivityCoefficientModel.davies,
            density_calculation=DensityCalculation.constant,
        )

        m.fs.unit = NanofiltrationDSPMDE0D(
            property_package=m.fs.properties,
            mass_transfer_coefficient=MassTransferCoefficient.none,
            concentration_polarization_type=ConcentrationPolarizationType.none,
        )
        b = m.fs.unit
        mass_flow_in = 1 * pyunits.kg / pyunits.s
        feed_mass_frac = {
            "Cl_-": 20101.6e-6,
            "Na_+": 11122e-6,
        }

        # Fix mole flow rates of each ion and water
        for ion, x in feed_mass_frac.items():
            mol_comp_flow = (
                x
                * pyunits.kg
                / pyunits.kg
                * mass_flow_in
                / m.fs.unit.feed_side.properties_in[0].mw_comp[ion]
            )
            m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", ion].fix(mol_comp_flow)
        H2O_mass_frac = 1 - sum(x for x in feed_mass_frac.values())
        H2O_mol_comp_flow = (
            H2O_mass_frac
            * pyunits.kg
            / pyunits.kg
            * mass_flow_in
            / m.fs.unit.feed_side.properties_in[0].mw_comp["H2O"]
        )
        m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix(H2O_mol_comp_flow)

        # Use assert electroneutrality method from property model to ensure the ion concentrations provided
        # obey electroneutrality condition
        m.fs.unit.feed_side.properties_in[0].assert_electroneutrality(
            defined_state=True,
            adjust_by_ion="Cl_-",
            get_property="mass_frac_phase_comp",
        )

        # Fix other inlet state variables
        m.fs.unit.inlet.temperature[0].fix(298.15)
        m.fs.unit.inlet.pressure[0].fix(4e5)

        # Fix the membrane variables that are usually fixed for the DSPM-DE model
        m.fs.unit.radius_pore.fix(0.5e-9)
        m.fs.unit.membrane_thickness_effective.fix(1.33e-6)
        m.fs.unit.membrane_charge_density.fix(-27)
        m.fs.unit.dielectric_constant_pore.fix(41.3)

        # Fix final permeate pressure to be ~atmospheric
        m.fs.unit.mixed_permeate[0].pressure.fix(101325)

        m.fs.unit.spacer_porosity.fix(0.85)
        m.fs.unit.channel_height.fix(5e-4)
        m.fs.unit.velocity[0, 0].fix(0.25)
        m.fs.unit.area.fix(50)

        return m

    @pytest.mark.unit
    def test_build(self, NF_frame):
        m = NF_frame

        # test ports and variables
        port_lst = ["inlet", "retentate", "permeate"]
        port_vars_lst = ["flow_mol_phase_comp", "pressure", "temperature"]
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
        cv_stateblock_lst = [
            "properties_in",
            "properties_out",
            "properties_interface",
        ]
        # feed side
        for sb_str in cv_stateblock_lst:
            sb = getattr(m.fs.unit.feed_side, sb_str)
            assert isinstance(sb, DSPMDEStateBlock)
        # test objects added to control volume
        cv_objs_type_dict = {"eq_feed_interface_isothermal": Constraint}
        for (obj_str, obj_type) in cv_objs_type_dict.items():
            obj = getattr(m.fs.unit.feed_side, obj_str)
            assert isinstance(obj, obj_type)
        # permeate side
        assert isinstance(m.fs.unit.permeate_side, DSPMDEStateBlock)
        assert isinstance(m.fs.unit.mixed_permeate, DSPMDEStateBlock)
        # membrane
        assert isinstance(m.fs.unit.pore_entrance, DSPMDEStateBlock)
        assert isinstance(m.fs.unit.pore_exit, DSPMDEStateBlock)

        # test statistics
        assert number_variables(m) == 302
        assert number_total_constraints(m) == 272
        assert number_unused_variables(m) == 11

    @pytest.mark.unit
    def test_dof(self, NF_frame):
        m = NF_frame
        check_dof(m, fail_flag=True)

    @pytest.mark.unit
    def test_calculate_scaling(self, NF_frame):
        m = NF_frame

        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e2, index=("Liq", "Cl_-")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e2, index=("Liq", "Na_+")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e0, index=("Liq", "H2O")
        )

        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m.fs.unit))
        assert len(unscaled_var_list) == 0

        badly_scaled_var_lst = list(badly_scaled_var_generator(m.fs.unit))
        assert len(badly_scaled_var_lst) == 0

        # not all constraints have scaling factor so skipping the check for unscaled constraints
        unscaled_con_lst = list(unscaled_constraints_generator(m))
        assert len(unscaled_con_lst) == 0

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_initialize(self, NF_frame):
        m = NF_frame
        # Using the 'initialize' function so that I can view the logs on failure
        # m.fs.unit.initialize(optarg=solver.options, outlvl=idaeslog.DEBUG)
        initialization_tester(m)

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_solve(self, NF_frame):
        m = NF_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_conservation(self, NF_frame):
        m = NF_frame
        b = m.fs.unit
        comp_lst = m.fs.properties.solute_set

        flow_mass_inlet = sum(
            b.feed_side.properties_in[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_retentate = sum(
            b.feed_side.properties_out[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_permeate = sum(
            b.mixed_permeate[0].flow_mass_phase_comp["Liq", j] for j in comp_lst
        )

        assert (
            abs(value(flow_mass_inlet - flow_mass_retentate - flow_mass_permeate))
            <= 1e-6
        )

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_solution(self, NF_frame):
        m = NF_frame

        mole_flux_dict = {
            "Na_+": 0.0010636,
            "Cl_-": 0.0010636,
            "H2O": 0.1313564,
        }
        for j, val in mole_flux_dict.items():
            assert pytest.approx(val, rel=1e-3) == value(
                m.fs.unit.flux_mol_phase_comp_avg[0, "Liq", j]
            )

        # TODO: subsequently focus on the segment below during the validation and refinement phase
        intrinsic_rejection_dict = {
            "Na_+": 0.096269,
            "Cl_-": 0.096269,
        }
        for j, val in intrinsic_rejection_dict.items():
            assert pytest.approx(val, rel=1e-3) == value(
                m.fs.unit.rejection_intrinsic_phase_comp[0, "Liq", j]
            )

    @pytest.mark.unit
    def test_report(self, NF_frame):
        NF_frame.fs.unit.report()


class TestNanoFiltration_with_CP_5ions_double_concentration:
    @pytest.fixture(scope="class")
    def NF_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = DSPMDEParameterBlock(
            solute_list=["Ca_2+", "SO4_2-", "Mg_2+", "Na_+", "Cl_-"],
            diffusivity_data={
                ("Liq", "Ca_2+"): 9.2e-10,
                ("Liq", "SO4_2-"): 1.06e-09,
                ("Liq", "Mg_2+"): 7.06e-10,
                ("Liq", "Na_+"): 1.33e-09,
                ("Liq", "Cl_-"): 2.03e-09,
            },
            mw_data={
                "H2O": 0.018,
                "Ca_2+": 0.04,
                "Mg_2+": 0.024,
                "SO4_2-": 0.096,
                "Na_+": 0.023,
                "Cl_-": 0.035,
            },
            stokes_radius_data={
                "Ca_2+": 3.09e-10,
                "Mg_2+": 3.47e-10,
                "SO4_2-": 2.3e-10,
                "Cl_-": 1.21e-10,
                "Na_+": 1.84e-10,
            },
            charge={"Ca_2+": 2, "Mg_2+": 2, "SO4_2-": -2, "Na_+": 1, "Cl_-": -1},
            activity_coefficient_model=ActivityCoefficientModel.davies,
            density_calculation=DensityCalculation.constant,
        )

        m.fs.unit = NanofiltrationDSPMDE0D(property_package=m.fs.properties)
        b = m.fs.unit
        mass_flow_in = 1 * pyunits.kg / pyunits.s
        feed_mass_frac = {
            "Ca_2+": 2 * 382e-6,
            "Mg_2+": 2 * 1394e-6,
            "SO4_2-": 2 * 2136e-6,
            "Cl_-": 2 * 20101.6e-6,
            "Na_+": 2 * 11122e-6,
        }

        # Fix mole flow rates of each ion and water
        for ion, x in feed_mass_frac.items():
            mol_comp_flow = (
                x
                * pyunits.kg
                / pyunits.kg
                * mass_flow_in
                / m.fs.unit.feed_side.properties_in[0].mw_comp[ion]
            )
            m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", ion].fix(mol_comp_flow)
        H2O_mass_frac = 1 - sum(x for x in feed_mass_frac.values())
        H2O_mol_comp_flow = (
            H2O_mass_frac
            * pyunits.kg
            / pyunits.kg
            * mass_flow_in
            / m.fs.unit.feed_side.properties_in[0].mw_comp["H2O"]
        )
        m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix(H2O_mol_comp_flow)

        # Use assert electroneutrality method from property model to ensure the ion concentrations provided
        # obey electroneutrality condition
        m.fs.unit.feed_side.properties_in[0].assert_electroneutrality(
            defined_state=True,
            adjust_by_ion="Cl_-",
            get_property="mass_frac_phase_comp",
        )

        # Fix other inlet state variables
        m.fs.unit.inlet.temperature[0].fix(298.15)
        m.fs.unit.inlet.pressure[0].fix(4e5)

        # Fix the membrane variables that are usually fixed for the DSPM-DE model
        m.fs.unit.radius_pore.fix(0.5e-9)
        m.fs.unit.membrane_thickness_effective.fix(1.33e-6)
        m.fs.unit.membrane_charge_density.fix(-27)
        m.fs.unit.dielectric_constant_pore.fix(41.3)

        # Fix final permeate pressure to be ~atmospheric
        m.fs.unit.mixed_permeate[0].pressure.fix(101325)

        m.fs.unit.spacer_porosity.fix(0.85)
        m.fs.unit.channel_height.fix(5e-4)
        m.fs.unit.velocity[0, 0].fix(0.25)
        m.fs.unit.area.fix(50)
        # Fix additional variables for calculating mass transfer coefficient with spiral wound correlation
        m.fs.unit.spacer_mixing_efficiency.fix()
        m.fs.unit.spacer_mixing_length.fix()
        return m

    @pytest.mark.unit
    def test_build(self, NF_frame):
        m = NF_frame

        # test ports and variables
        port_lst = ["inlet", "retentate", "permeate"]
        port_vars_lst = ["flow_mol_phase_comp", "pressure", "temperature"]
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
        cv_stateblock_lst = [
            "properties_in",
            "properties_out",
            "properties_interface",
        ]
        # feed side
        for sb_str in cv_stateblock_lst:
            sb = getattr(m.fs.unit.feed_side, sb_str)
            assert isinstance(sb, DSPMDEStateBlock)
        # test objects added to control volume
        cv_objs_type_dict = {"eq_feed_interface_isothermal": Constraint}
        for (obj_str, obj_type) in cv_objs_type_dict.items():
            obj = getattr(m.fs.unit.feed_side, obj_str)
            assert isinstance(obj, obj_type)
        # permeate side
        assert isinstance(m.fs.unit.permeate_side, DSPMDEStateBlock)
        assert isinstance(m.fs.unit.mixed_permeate, DSPMDEStateBlock)
        # membrane
        assert isinstance(m.fs.unit.pore_entrance, DSPMDEStateBlock)
        assert isinstance(m.fs.unit.pore_exit, DSPMDEStateBlock)

        # test statistics
        assert number_variables(m) == 561
        assert number_total_constraints(m) == 526
        assert number_unused_variables(m) == 11

    @pytest.mark.unit
    def test_dof(self, NF_frame):
        m = NF_frame
        check_dof(m, fail_flag=True)

    @pytest.mark.unit
    def test_calculate_scaling(self, NF_frame):
        m = NF_frame

        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e2, index=("Liq", "Ca_2+")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e2, index=("Liq", "SO4_2-")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e2, index=("Liq", "Mg_2+")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e2, index=("Liq", "Cl_-")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e2, index=("Liq", "Na_+")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e0, index=("Liq", "H2O")
        )

        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m.fs.unit))
        assert len(unscaled_var_list) == 0

        badly_scaled_var_lst = list(badly_scaled_var_generator(m.fs.unit))
        assert len(badly_scaled_var_lst) == 0

        # not all constraints have scaling factor so skipping the check for unscaled constraints
        unscaled_con_lst = list(unscaled_constraints_generator(m))
        assert len(unscaled_con_lst) == 0

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_initialize(self, NF_frame):
        m = NF_frame
        # Using the 'initialize' function so that I can view the logs on failure
        # m.fs.unit.initialize(optarg=solver.options, outlvl=idaeslog.DEBUG)
        initialization_tester(m)

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_solve(self, NF_frame):
        m = NF_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_conservation(self, NF_frame):
        m = NF_frame
        b = m.fs.unit
        comp_lst = m.fs.properties.solute_set

        flow_mass_inlet = sum(
            b.feed_side.properties_in[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_retentate = sum(
            b.feed_side.properties_out[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_permeate = sum(
            b.mixed_permeate[0].flow_mass_phase_comp["Liq", j] for j in comp_lst
        )

        assert (
            abs(value(flow_mass_inlet - flow_mass_retentate - flow_mass_permeate))
            <= 1e-6
        )

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_solution(self, NF_frame):
        m = NF_frame

        mole_flux_dict = {
            "Na_+": 0.0001269756,
            "Cl_-": 0.0001352,
            "Ca_2+": 1.4258536e-06,
            "SO4_2-": 2.9210354e-06,
            "Mg_2+": 5.5916186e-06,
            "H2O": 0.0069096,
        }
        for j, val in mole_flux_dict.items():
            assert pytest.approx(val, rel=1e-3) == value(
                m.fs.unit.flux_mol_phase_comp_avg[0, "Liq", j]
            )

        # TODO: subsequently focus on the segment below during the validation and refinement phase
        intrinsic_rejection_dict = {
            "Na_+": 0.00753,
            "Cl_-": 0.11317679,
            "Ca_2+": 0.441427,
            "SO4_2-": 0.510229,
            "Mg_2+": 0.64273,
        }
        for j, val in intrinsic_rejection_dict.items():
            assert pytest.approx(val, rel=1e-3) == value(
                m.fs.unit.rejection_intrinsic_phase_comp[0, "Liq", j]
            )

    @pytest.mark.unit
    def test_report(self, NF_frame):
        NF_frame.fs.unit.report()


@pytest.mark.requires_idaes_solver
@pytest.mark.component
def test_inverse_solve():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = DSPMDEParameterBlock(
        solute_list=["Na_+", "Cl_-"],
        diffusivity_data={("Liq", "Na_+"): 1.33e-09, ("Liq", "Cl_-"): 2.03e-09},
        mw_data={"H2O": 0.018, "Na_+": 0.023, "Cl_-": 0.035},
        stokes_radius_data={"Cl_-": 1.21e-10, "Na_+": 1.84e-10},
        charge={"Na_+": 1, "Cl_-": -1},
        activity_coefficient_model=ActivityCoefficientModel.davies,
        density_calculation=DensityCalculation.constant,
    )

    m.fs.unit = NanofiltrationDSPMDE0D(property_package=m.fs.properties)

    m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(0.429868)
    m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(0.429868)
    m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix(47.356)

    # Fix other inlet state variables
    m.fs.unit.inlet.temperature[0].fix(298.15)
    m.fs.unit.inlet.pressure[0].fix(4e5)

    # Fix the membrane variables that are usually fixed for the DSPM-DE model
    m.fs.unit.radius_pore.fix(0.5e-9)
    m.fs.unit.membrane_thickness_effective.fix(1.33e-6)
    m.fs.unit.membrane_charge_density.fix(-27)
    m.fs.unit.dielectric_constant_pore.fix(41.3)

    # Fix final permeate pressure to be ~atmospheric
    m.fs.unit.mixed_permeate[0].pressure.fix(101325)

    m.fs.unit.spacer_porosity.fix(0.85)
    m.fs.unit.channel_height.fix(5e-4)
    m.fs.unit.velocity[0, 0].fix(0.25)
    m.fs.unit.area.fix(50)
    # Fix additional variables for calculating mass transfer coefficient with spiral wound correlation
    m.fs.unit.spacer_mixing_efficiency.fix()
    m.fs.unit.spacer_mixing_length.fix()

    check_dof(m, fail_flag=True)

    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e2, index=("Liq", "Cl_-")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e2, index=("Liq", "Na_+")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e0, index=("Liq", "H2O")
    )

    calculate_scaling_factors(m)

    # check that all variables have scaling factors
    unscaled_var_list = list(unscaled_variables_generator(m.fs.unit))
    assert len(unscaled_var_list) == 0

    badly_scaled_var_lst = list(badly_scaled_var_generator(m.fs.unit))
    assert len(badly_scaled_var_lst) == 0

    # not all constraints have scaling factor so skipping the check for unscaled constraints
    unscaled_con_lst = list(unscaled_constraints_generator(m))
    assert len(unscaled_con_lst) == 0

    initialization_tester(m)

    m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Na_+"].unfix()
    m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Cl_-"].unfix()
    m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].unfix()
    m.fs.unit.inlet.temperature[0].unfix()
    m.fs.unit.inlet.pressure[0].unfix()

    m.fs.unit.retentate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(0.429868)
    m.fs.unit.retentate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(0.429868)
    m.fs.unit.retentate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(47.356)
    m.fs.unit.retentate.temperature[0].fix(298.15)
    m.fs.unit.retentate.pressure[0].fix(4e5)

    results = solver.solve(m)

    # Check for optimal solution
    assert_optimal_termination(results)

    b = m.fs.unit
    comp_lst = m.fs.properties.solute_set

    flow_mass_inlet = sum(
        b.feed_side.properties_in[0].flow_mass_phase_comp["Liq", j] for j in comp_lst
    )
    flow_mass_retentate = sum(
        b.feed_side.properties_out[0].flow_mass_phase_comp["Liq", j] for j in comp_lst
    )
    flow_mass_permeate = sum(
        b.mixed_permeate[0].flow_mass_phase_comp["Liq", j] for j in comp_lst
    )

    assert (
        abs(value(flow_mass_inlet - flow_mass_retentate - flow_mass_permeate)) <= 1e-6
    )

    mole_flux_dict = {
        "Na_+": 0.0010739,
        "Cl_-": 0.0010739,
        "H2O": 0.129305,
    }
    for j, val in mole_flux_dict.items():
        assert pytest.approx(val, rel=1e-3) == value(
            m.fs.unit.flux_mol_phase_comp_avg[0, "Liq", j]
        )

    intrinsic_rejection_dict = {
        "Na_+": 0.0964865,
        "Cl_-": 0.0964865,
    }
    for j, val in intrinsic_rejection_dict.items():
        assert pytest.approx(val, rel=1e-3) == value(
            m.fs.unit.rejection_intrinsic_phase_comp[0, "Liq", j]
        )
    flow_mol_in_dict = {"Na_+": 0.483564, "Cl_-": 0.483564, "H2O": 53.821}
    for j, val in flow_mol_in_dict.items():
        assert pytest.approx(val, rel=1e-3) == value(
            m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", j]
        )

    m.fs.unit.report()


@pytest.mark.requires_idaes_solver
@pytest.mark.component
def test_mass_transfer_coeff_fixed():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = DSPMDEParameterBlock(
        solute_list=["Na_+", "Cl_-"],
        diffusivity_data={("Liq", "Na_+"): 1.33e-09, ("Liq", "Cl_-"): 2.03e-09},
        mw_data={"H2O": 0.018, "Na_+": 0.023, "Cl_-": 0.035},
        stokes_radius_data={"Cl_-": 1.21e-10, "Na_+": 1.84e-10},
        charge={"Na_+": 1, "Cl_-": -1},
        activity_coefficient_model=ActivityCoefficientModel.davies,
        density_calculation=DensityCalculation.constant,
    )

    m.fs.unit = NanofiltrationDSPMDE0D(
        property_package=m.fs.properties,
        mass_transfer_coefficient=MassTransferCoefficient.fixed,
    )

    m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(0.429868)
    m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(0.429868)
    m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix(47.356)

    # Fix other inlet state variables
    m.fs.unit.inlet.temperature[0].fix(298.15)
    m.fs.unit.inlet.pressure[0].fix(4e5)

    # Fix the membrane variables that are usually fixed for the DSPM-DE model
    m.fs.unit.radius_pore.fix(0.5e-9)
    m.fs.unit.membrane_thickness_effective.fix(1.33e-6)
    m.fs.unit.membrane_charge_density.fix(-27)
    m.fs.unit.dielectric_constant_pore.fix(41.3)

    # Fix final permeate pressure to be ~atmospheric
    m.fs.unit.mixed_permeate[0].pressure.fix(101325)

    m.fs.unit.spacer_porosity.fix(0.85)
    m.fs.unit.channel_height.fix(5e-4)
    m.fs.unit.velocity[0, 0].fix(0.25)
    m.fs.unit.area.fix(50)

    # Fix mass transfer coefficient for each ion
    m.fs.unit.Kf_comp.fix(1e-5)

    check_dof(m, fail_flag=True)

    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e2, index=("Liq", "Cl_-")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e2, index=("Liq", "Na_+")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e0, index=("Liq", "H2O")
    )

    calculate_scaling_factors(m)

    # check that all variables have scaling factors
    unscaled_var_list = list(unscaled_variables_generator(m.fs.unit))
    assert len(unscaled_var_list) == 0

    badly_scaled_var_lst = list(badly_scaled_var_generator(m.fs.unit))
    assert len(badly_scaled_var_lst) == 0

    # not all constraints have scaling factor so skipping the check for unscaled constraints
    unscaled_con_lst = list(unscaled_constraints_generator(m))
    assert len(unscaled_con_lst) == 0

    initialization_tester(m)

    m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Na_+"].unfix()
    m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Cl_-"].unfix()
    m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].unfix()
    m.fs.unit.inlet.temperature[0].unfix()
    m.fs.unit.inlet.pressure[0].unfix()

    m.fs.unit.retentate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(0.429868)
    m.fs.unit.retentate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(0.429868)
    m.fs.unit.retentate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(47.356)
    m.fs.unit.retentate.temperature[0].fix(298.15)
    m.fs.unit.retentate.pressure[0].fix(4e5)

    results = solver.solve(m)

    # Check for optimal solution
    assert_optimal_termination(results)

    b = m.fs.unit
    comp_lst = m.fs.properties.solute_set

    flow_mass_inlet = sum(
        b.feed_side.properties_in[0].flow_mass_phase_comp["Liq", j] for j in comp_lst
    )
    flow_mass_retentate = sum(
        b.feed_side.properties_out[0].flow_mass_phase_comp["Liq", j] for j in comp_lst
    )
    flow_mass_permeate = sum(
        b.mixed_permeate[0].flow_mass_phase_comp["Liq", j] for j in comp_lst
    )

    assert (
        abs(value(flow_mass_inlet - flow_mass_retentate - flow_mass_permeate)) <= 1e-6
    )

    mole_flux_dict = {
        "Na_+": 0.00107097,
        "Cl_-": 0.00107097,
        "H2O": 0.12989686,
    }
    for j, val in mole_flux_dict.items():
        assert pytest.approx(val, rel=1e-3) == value(
            m.fs.unit.flux_mol_phase_comp_avg[0, "Liq", j]
        )

    intrinsic_rejection_dict = {
        "Na_+": 0.09700,
        "Cl_-": 0.09700,
    }
    for j, val in intrinsic_rejection_dict.items():
        assert pytest.approx(val, rel=1e-3) == value(
            m.fs.unit.rejection_intrinsic_phase_comp[0, "Liq", j]
        )
    flow_mol_in_dict = {"Na_+": 0.483564, "Cl_-": 0.483564, "H2O": 53.821}
    for j, val in flow_mol_in_dict.items():
        assert pytest.approx(val, rel=1e-3) == value(
            m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", j]
        )


@pytest.mark.unit
def test_mass_transfer_CP_config_errors():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = DSPMDEParameterBlock(
        solute_list=["Na_+", "Cl_-"],
        diffusivity_data={("Liq", "Na_+"): 1.33e-09, ("Liq", "Cl_-"): 2.03e-09},
        mw_data={"H2O": 0.018, "Na_+": 0.023, "Cl_-": 0.035},
        stokes_radius_data={"Cl_-": 1.21e-10, "Na_+": 1.84e-10},
        charge={"Na_+": 1, "Cl_-": -1},
        activity_coefficient_model=ActivityCoefficientModel.davies,
        density_calculation=DensityCalculation.constant,
    )

    with pytest.raises(
        ConfigurationError,
        match="\nConflict between configuration options:\n"
        "'mass_transfer_coefficient' cannot be set to MassTransferCoefficient.fixed "
        "while 'concentration_polarization_type' is set to ConcentrationPolarizationType.none.\n\n"
        "'mass_transfer_coefficient' must be set to MassTransferCoefficient.none\nor "
        "'concentration_polarization_type' must be set to ConcentrationPolarizationType.calculated",
    ):
        m.fs.unit = NanofiltrationDSPMDE0D(
            property_package=m.fs.properties,
            mass_transfer_coefficient=MassTransferCoefficient.fixed,
            concentration_polarization_type=ConcentrationPolarizationType.none,
        )

    with pytest.raises(
        ConfigurationError,
        match="\nConflict between configuration options:\n"
        "'mass_transfer_coefficient' cannot be set to MassTransferCoefficient.none "
        "while 'concentration_polarization_type' is set to ConcentrationPolarizationType.calculated.\n\n"
        "'mass_transfer_coefficient' must be set to MassTransferCoefficient.none\nor "
        "'concentration_polarization_type' must be set to ConcentrationPolarizationType.calculated",
    ):
        m.fs.unit = NanofiltrationDSPMDE0D(
            property_package=m.fs.properties,
            mass_transfer_coefficient=MassTransferCoefficient.none,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
        )

    with pytest.raises(
        ConfigurationError,
        match="\nConflict between configuration options:\n"
        "'mass_transfer_coefficient' cannot be set to MassTransferCoefficient.spiral_wound "
        "while 'concentration_polarization_type' is set to ConcentrationPolarizationType.none.\n\n"
        "'mass_transfer_coefficient' must be set to MassTransferCoefficient.none\nor "
        "'concentration_polarization_type' must be set to ConcentrationPolarizationType.calculated",
    ):
        m.fs.unit = NanofiltrationDSPMDE0D(
            property_package=m.fs.properties,
            mass_transfer_coefficient=MassTransferCoefficient.spiral_wound,
            concentration_polarization_type=ConcentrationPolarizationType.none,
        )
