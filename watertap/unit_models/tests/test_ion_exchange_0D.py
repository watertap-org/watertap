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
    Param,
    Constraint,
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
from watertap.property_models.IX_prop_pack import (
    IXParameterBlock,
    IXStateBlock,
)
from watertap.unit_models.ion_exchange_0D import IonExchange0D
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
import idaes.logger as idaeslog

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

ix_in = {
    "solute_list": [
        "Na_+",
    ],
    "diffusivity_data": {
        ("Liq", "Na_+"): 1.33e-9,
    },
    "mw_data": {
        "H2O": 18e-3,
        "Na_+": 23e-3,
    },
    "charge": {"Na_+": 1},
}

ion = "Na_+"


# @pytest.mark.unit
# def test_config():
#     m = ConcreteModel()
#     m.fs = FlowsheetBlock(default={"dynamic": False})

#     m.fs.properties = IXParameterBlock(
#         default=ix_in
#     )
#     m.fs.unit = IonExchange0D(default={"property_package": m.fs.properties})
#     # check unit config arguments
#     assert len(m.fs.unit.config) == 5

#     assert not m.fs.unit.config.dynamic
#     assert not m.fs.unit.config.has_holdup

#     assert m.fs.unit.config.property_package is m.fs.properties


class TestIonExchange:
    @pytest.fixture(scope="class")
    def IX_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.properties = IXParameterBlock(default=ix_in)

        m.fs.unit = IonExchange0D(default={"property_package": m.fs.properties})
        mass_flow_in = 1 * pyunits.kg / pyunits.s
        feed_mass_frac = {
            "Na_+": 5e-4,
        }
        ix = m.fs.unit
        # Fix mole flow rates of each ion and water
        for ion, x in feed_mass_frac.items():
            mol_comp_flow = (
                x
                * pyunits.kg
                / pyunits.kg
                * mass_flow_in
                / m.fs.unit.config.property_package.mw_comp[ion]
            )
            ix.inlet.flow_mol_phase_comp[0, "Liq", ion].fix(mol_comp_flow)
        H2O_mass_frac = 1 - sum(x for x in feed_mass_frac.values())
        H2O_mol_comp_flow = (
            H2O_mass_frac
            * pyunits.kg
            / pyunits.kg
            * mass_flow_in
            / m.fs.unit.config.property_package.mw_comp["H2O"]
        )
        ix.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix(H2O_mol_comp_flow)
        ix.resin_diam.fix()
        ix.K_eq[ion].fix(1.5)
        ix.resin_max_capacity.fix(5)
        ix.bed_porosity.fix()
        ix.resin_bulk_dens.fix()
        ix.dimensionless_time.fix(1)
        ix.lh.fix(0)
        ix.bed_depth.fix(3)

        return m

    @pytest.mark.unit
    def test_config(self, IX_frame):
        m = IX_frame

        assert len(m.fs.unit.config) == 5

        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup

        assert m.fs.unit.config.property_package is m.fs.properties

    @pytest.mark.unit
    def test_build(self, IX_frame):
        m = IX_frame
        ix = m.fs.unit
        # test ports and variables
        port_lst = ["inlet", "outlet", "waste"]
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
        ix_params = [
            "diff_ion_resin",
            "underdrain_h",
            "distributor_h",
            "holdup_A",
            "holdup_B",
            "holdup_exp",
            "Pe_p_A",
            "Pe_p_exp",
            "Sh_A",
            "Sh_exp",
            "t_waste_param",
            "vel_bed_ratio",
            "bed_depth_to_diam_ratio",
        ]
        for p in ix_params:
            assert hasattr(ix, p)
            param = getattr(ix, p)
            assert isinstance(param, Param)

        ix_vars = [
            "resin_max_capacity",
            "resin_eq_capacity",
            "resin_diam",
            "resin_bulk_dens",
            "resin_particle_dens",
            "K_eq",
            "R_eq",
            "resin_surf_per_vol",
            "bed_vol",
            "bed_diam",
            "bed_depth",
            "bed_area",
            "bed_porosity",
            "col_height",
            "col_vol",
            "partition_ratio",
            "fluid_mass_transfer_coeff",
            "rate_coeff",
            "t_breakthru",
            "t_contact",
            "t_waste",
            "num_transfer_units",
            "HTU",
            "dimensionless_time",
            "lh",
            "mass_in",
            "mass_removed",
            "mass_out",
            "vel_min",
            "vel_bed",
            "vel_inter",
            "sfr",
            "Re",
            "Sc",
            "Sh",
            "Pe_p",
            "Pe_bed",
            "c_norm",
        ]

        for v in ix_vars:
            assert hasattr(ix, v)
            param = getattr(ix, v)
            assert isinstance(param, Var)

        # # test state block objects
        cv_stateblock_lst = [
            "properties_in",
            "properties_out",
            "properties_waste",
        ]
        # feed side
        for sb_str in cv_stateblock_lst:
            sb = getattr(ix, sb_str)
            assert isinstance(sb, IXStateBlock)
        # # test objects added to control volume
        # cv_objs_type_dict = {"eq_feed_interface_isothermal": Constraint}
        # for (obj_str, obj_type) in cv_objs_type_dict.items():
        #     obj = getattr(m.fs.unit.feed_side, obj_str)
        #     assert isinstance(obj, obj_type)
        # # permeate side
        # assert isinstance(m.fs.unit.permeate_side, IXStateBlock)
        # assert isinstance(m.fs.unit.mixed_permeate, IXStateBlock)
        # # membrane
        # assert isinstance(m.fs.unit.pore_entrance, IXStateBlock)
        # assert isinstance(m.fs.unit.pore_exit, IXStateBlock)

        # # test statistics
        # assert number_variables(m) == 558
        # assert number_total_constraints(m) == 524
        # assert number_unused_variables(m) == 11

    # @pytest.mark.unit
    # def test_dof(self, NF_frame):
    #     m = NF_frame
    #     check_dof(m, fail_flag=True)

    # @pytest.mark.unit
    # def test_calculate_scaling(self, NF_frame):
    #     m = NF_frame

    #     m.fs.properties.set_default_scaling(
    #         "flow_mol_phase_comp", 1e3, index=("Liq", "Ca_2+")
    #     )
    #     m.fs.properties.set_default_scaling(
    #         "flow_mol_phase_comp", 1e3, index=("Liq", "SO4_2-")
    #     )

    #     calculate_scaling_factors(m)

    #     # check that all variables have scaling factors
    #     unscaled_var_list = list(unscaled_variables_generator(m.fs.unit))
    #     assert len(unscaled_var_list) == 0

    #     badly_scaled_var_lst = list(badly_scaled_var_generator(m, include_fixed=True))
    #     assert len(unscaled_var_list) == 0

    #     # not all constraints have scaling factor so skipping the check for unscaled constraints

    # @pytest.mark.skip(
    #     reason="Scaling and/or formulation of unit model needs to be revisited"
    # )
    # @pytest.mark.requires_idaes_solver
    # @pytest.mark.component
    # def test_initialize(self, NF_frame):
    #     m = NF_frame
    #     # Using the 'initialize' function so that I can view the logs on failure
    #     m.fs.unit.initialize(optarg=solver.options, outlvl=idaeslog.DEBUG)
    #     # initialization_tester(m)

    # @pytest.mark.skip(
    #     reason="Scaling and/or formulation of unit model needs to be revisited"
    # )
    # @pytest.mark.requires_idaes_solver
    # @pytest.mark.component
    # def test_solve(self, NF_frame):
    #     m = NF_frame
    #     results = solver.solve(m)

    #     # Check for optimal solution
    #     assert_optimal_termination(results)

    # @pytest.mark.skip(
    #     reason="Scaling and/or formulation of unit model needs to be revisited"
    # )
    # @pytest.mark.requires_idaes_solver
    # @pytest.mark.component
    # def test_conservation(self, NF_frame):
    #     m = NF_frame
    #     b = m.fs.unit
    #     comp_lst = m.fs.properties.solute_set

    #     flow_mass_inlet = sum(
    #         b.feed_side.properties_in[0].flow_mass_phase_comp["Liq", j]
    #         for j in comp_lst
    #     )
    #     flow_mass_retentate = sum(
    #         b.feed_side.properties_out[0].flow_mass_phase_comp["Liq", j]
    #         for j in comp_lst
    #     )
    #     flow_mass_permeate = sum(
    #         b.mixed_permeate[0].flow_mass_phase_comp["Liq", j] for j in comp_lst
    #     )

    #     assert (
    #         abs(value(flow_mass_inlet - flow_mass_retentate - flow_mass_permeate))
    #         <= 1e-6
    #     )

    # @pytest.mark.skip(
    #     reason="Scaling and/or formulation of unit model needs to be revisited"
    # )
    # @pytest.mark.requires_idaes_solver
    # @pytest.mark.component
    # def test_solution(self, NF_frame):
    #     m = NF_frame

    #     mole_flux_dict = {
    #         "Na_+": 0.0036645,
    #         "Cl_-": 0.004354,
    #         "Ca_2+": 7.3313e-5,
    #         "SO4_2-": 1.6845e-4,
    #         "Mg_2+": 0.000440,
    #         "H2O": 0.41437,
    #     }
    #     for j, val in mole_flux_dict.items():
    #         assert pytest.approx(val, rel=5e-2) == value(
    #             m.fs.unit.flux_mol_phase_comp_avg[0, "Liq", j]
    #         )

    #     # TODO: subsequently focus on the segment below during the validation and refinement phase
    #     intrinsic_rejection_dict = {
    #         "Na_+": 0.017432,
    #         "Cl_-": 0.014704,
    #         "Ca_2+": -0.034499,
    #         "SO4_2-": 0.01435907,
    #         "Mg_2+": 0.013672,
    #     }
    #     for j, val in intrinsic_rejection_dict.items():
    #         assert pytest.approx(val, rel=5e-2) == value(
    #             m.fs.unit.rejection_intrinsic_phase_comp[0, "Liq", j]
    #         )

    # @pytest.mark.unit
    # def test_report(self, NF_frame):
    #     NF_frame.fs.unit.report()
