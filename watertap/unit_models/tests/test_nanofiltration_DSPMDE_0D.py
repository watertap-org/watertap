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
)
from watertap.core.util.initialization import check_dof

from idaes.core.util import get_solver
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
    unused_variables_set,
)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    unscaled_constraints_generator,
    constraints_with_scale_factor_generator,
    badly_scaled_var_generator,
)

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = DSPMDEParameterBlock(
        default={"solute_list": ["Ca_2+", "SO4_2-", "Na_+", "Cl_-", "Mg_2+"]}
    )
    m.fs.unit = NanofiltrationDSPMDE0D(default={"property_package": m.fs.properties})

    # check unit config arguments
    assert len(m.fs.unit.config) == 8

    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
    assert m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
    assert not m.fs.unit.config.has_pressure_change
    assert m.fs.unit.config.property_package is m.fs.properties
    assert (
        m.fs.unit.config.property_package.config.activity_coefficient_model
        == ActivityCoefficientModel.ideal
    )
    assert (
        m.fs.unit.config.property_package.config.density_calculation
        == DensityCalculation.constant
    )
    assert m.fs.unit.config.mass_transfer_coefficient == MassTransferCoefficient.fixed


class TestNanoFiltration:
    @pytest.fixture(scope="class")
    def NF_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.properties = DSPMDEParameterBlock(
            default={
                "solute_list": [
                    "Ca_2+",
                    "SO4_2-",
                    "Mg_2+",
                    "Na_+",
                    "Cl_-",
                ],
                "diffusivity_data": {
                    ("Liq", "Ca_2+"): 9.2e-10,
                    ("Liq", "SO4_2-"): 1.06e-9,
                    ("Liq", "Mg_2+"): 0.706e-9,
                    ("Liq", "Na_+"): 1.33e-9,
                    ("Liq", "Cl_-"): 2.03e-9,
                },
                "mw_data": {
                    "H2O": 18e-3,
                    "Ca_2+": 40e-3,
                    "Mg_2+": 24e-3,
                    "SO4_2-": 96e-3,
                    "Na_+": 23e-3,
                    "Cl_-": 35e-3,
                },
                "stokes_radius_data": {
                    "Ca_2+": 0.309e-9,
                    "Mg_2+": 0.347e-9,
                    "SO4_2-": 0.230e-9,
                    "Cl_-": 0.121e-9,
                    "Na_+": 0.184e-9,
                },
                "charge": {
                    "Ca_2+": 2,
                    "Mg_2+": 2,
                    "SO4_2-": -2,
                    "Na_+": 1,
                    "Cl_-": -1,
                },
                "activity_coefficient_model": ActivityCoefficientModel.davies,
                "density_calculation": DensityCalculation.constant,
            }
        )

        m.fs.unit = NanofiltrationDSPMDE0D(
            default={
                "property_package": m.fs.properties,
                "mass_transfer_coefficient": MassTransferCoefficient.spiral_wound,
            }
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

        # Fix additional variables for calculating mass transfer coefficient with spiral wound correlation
        m.fs.unit.spacer_porosity.fix(0.85)
        m.fs.unit.spacer_mixing_efficiency.fix()
        m.fs.unit.spacer_mixing_length.fix()
        m.fs.unit.channel_height.fix(5e-4)
        m.fs.unit.recovery_vol_phase[0, "Liq"].fix(0.5)

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
        assert number_variables(m) == 558
        assert number_total_constraints(m) == 524
        assert number_unused_variables(m) == 11

    @pytest.mark.unit
    def test_dof(self, NF_frame):
        m = NF_frame
        check_dof(m, fail_flag=True)

    @pytest.mark.unit
    def test_calculate_scaling(self, NF_frame):
        m = NF_frame

        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e5, index=("Liq", "Ca_2+")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e5, index=("Liq", "SO4_2-")
        )

        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m.fs.unit))
        assert len(unscaled_var_list) == 0

        badly_scaled_var_lst = list(badly_scaled_var_generator(m, include_fixed=True))
        assert len(unscaled_var_list) == 0

        # not all constraints have scaling factor so skipping the check for unscaled constraints

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_initialize(self, NF_frame):
        m = NF_frame
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
            "Na_+": 0.0036645,
            "Cl_-": 0.004354,
            "Ca_2+": 7.3313e-5,
            "SO4_2-": 1.6845e-4,
            "Mg_2+": 0.000440,
            "H2O": 0.41437,
        }
        for j, val in mole_flux_dict.items():
            assert pytest.approx(val, rel=5e-2) == value(
                m.fs.unit.flux_mol_phase_comp_avg[0, "Liq", j]
            )

        # TODO: subsequently focus on the segment below during the validation and refinement phase
        intrinsic_rejection_dict = {
            "Na_+": 0.017432,
            "Cl_-": 0.015837,
            "Ca_2+": -0.023145,
            "SO4_2-": 0.015924,
            "Mg_2+": 0.01546,
        }
        for j, val in intrinsic_rejection_dict.items():
            assert pytest.approx(val, rel=5e-2) == value(
                m.fs.unit.rejection_intrinsic_phase_comp[0, "Liq", j]
            )

    @pytest.mark.unit
    def test_report(self, NF_frame):
        NF_frame.fs.unit.report()
