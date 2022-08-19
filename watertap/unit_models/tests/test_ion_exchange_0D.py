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
from watertap.property_models.ion_DSPMDE_prop_pack import (
    DSPMDEParameterBlock,
    ActivityCoefficientModel,
    DensityCalculation,
    DSPMDEStateBlock,
)
from watertap.unit_models.ion_exchange_0D import IonExchange0D, IonExchangeType

from watertap.core.util.initialization import check_dof

from idaes.core.solvers.get_solver import get_solver
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
    get_scaling_factor,
    unscaled_variables_generator,
    unscaled_constraints_generator,
    constraints_with_scale_factor_generator,
    badly_scaled_var_generator,
    set_scaling_factor,
)
import idaes.logger as idaeslog

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

ix_in = {
    "solute_list": ["Ca_2+"],
    "diffusivity_data": {("Liq", "Ca_2+"): 9.2e-10},
    "mw_data": {"H2O": 0.018, "Ca_2+": 0.04},
    "charge": {"Ca_2+": 2},
}

target_ion = "Ca_2+"


def ix_scaling(m, sf=1e5, est_recov=0.95, est_removal=0.99):
    ix = m.fs.unit
    ions = ix.config.property_package.ion_set
    target_ion = ix.config.target_ion
    prop_in = ix.properties_in[0]
    prop_out = ix.properties_out[0]
    prop_regen = ix.properties_regen[0]

    sf_flow_in = 1 / prop_in.flow_mol_phase_comp["Liq", "H2O"].value
    # sf_flow_out = 1 / (prop_in.flow_mol_phase_comp['Liq', 'H2O'].value * est_recov)
    sf_flow_out = 1
    sf_flow_regen = 1 / (
        prop_in.flow_mol_phase_comp["Liq", "H2O"].value * (1 - est_recov)
    )

    set_scaling_factor(prop_in.flow_mol_phase_comp["Liq", "H2O"], sf_flow_in)
    set_scaling_factor(prop_out.flow_mol_phase_comp["Liq", "H2O"], sf_flow_out)
    set_scaling_factor(prop_regen.flow_mol_phase_comp["Liq", "H2O"], sf_flow_regen)

    for ion in ions:
        set_scaling_factor(
            prop_in.flow_mol_phase_comp["Liq", ion],
            1 / prop_in.flow_mol_phase_comp["Liq", ion].value,
        )
        if ion == target_ion:
            set_scaling_factor(prop_out.mass_frac_phase_comp["Liq", ion], sf)
            set_scaling_factor(prop_out.flow_mass_phase_comp["Liq", ion], sf)
            set_scaling_factor(
                prop_out.flow_mol_phase_comp["Liq", ion],
                (
                    1
                    / (
                        prop_in.flow_mol_phase_comp["Liq", ion].value
                        * (1 - est_removal)
                    )
                ),
            )
            set_scaling_factor(
                prop_regen.flow_mol_phase_comp["Liq", ion],
                1 / (prop_in.flow_mol_phase_comp["Liq", ion].value * est_removal),
            )
        else:
            set_scaling_factor(
                prop_out.flow_mol_phase_comp["Liq", ion],
                1 / (prop_in.flow_mol_phase_comp["Liq", "H2O"].value * est_recov),
            )
            set_scaling_factor(prop_regen.flow_mol_phase_comp["Liq", ion], sf)

    return m


class TestIonExchange:
    @pytest.fixture(scope="class")
    def IX_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.properties = DSPMDEParameterBlock(default=ix_in)

        ix = m.fs.unit = IonExchange0D(
            default={"property_package": m.fs.properties, "target_ion": target_ion}
        )
        mass_flow_in = 5 * pyunits.kg / pyunits.s
        feed_mass_frac = {
            "Ca_2+": 2e-4,
        }
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

        ix.inlet.pressure[0].fix(101325)
        ix.inlet.temperature[0].fix(298.15)
        ix.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix(H2O_mol_comp_flow)
        ix.selectivity[ion].fix(3)
        ix.resin_max_capacity.fix(5)
        ix.service_flow_rate.fix(10)
        ix.number_columns.fix(4)
        ix.resin_diam.fix()
        ix.resin_bulk_dens.fix()
        ix.bed_porosity.fix()
        ix.dimensionless_time.fix()
        ix.lh.fix()
        ix.regen_dose.fix()
        ix.regen_ww.fix()
        ix.regen_sg.fix()
        ix.t_regen.fix()
        ix.t_bw.fix()
        ix.bw_rate.fix()
        ix.rinse_bv.fix()
        ix.pump_efficiency.fix()
        ix.p_drop_A.fix()
        ix.p_drop_B.fix()
        ix.p_drop_C.fix()
        ix.bed_expansion_frac_A.fix()
        ix.bed_expansion_frac_B.fix()
        ix.bed_expansion_frac_C.fix()
        ix.service_to_regen_flow_ratio.fix()

        return m

    @pytest.mark.unit
    def test_config(self, IX_frame):
        m = IX_frame

        assert len(m.fs.unit.config) == 7

        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert m.fs.unit.config.ion_exchange_type is IonExchangeType.cation
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
        ]

        for p in ix_params:
            assert hasattr(ix, p)
            param = getattr(ix, p)
            assert isinstance(param, Param)

        ix_vars = [
            "bed_depth_to_diam_ratio",
            "bed_expansion_frac_A",
            "bed_expansion_frac_B",
            "bed_expansion_frac_C",
            "p_drop_A",
            "p_drop_B",
            "p_drop_C",
            "resin_max_capacity",
            "resin_eq_capacity",
            "resin_unused_capacity",
            "resin_diam",
            "resin_bulk_dens",
            "resin_particle_dens",
            "selectivity",
            "separation_factor",
            "resin_surf_per_vol",
            "bed_vol",
            "bed_depth",
            "bed_area",
            "bed_porosity",
            "col_height",
            "col_vol_tot",
            "col_vol_per",
            "col_diam",
            "number_columns",
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
            "vel_bed",
            "vel_inter",
            "holdup",
            "service_flow_rate",
            "pressure_drop",
            "Re",
            "Sc",
            "Sh",
            "Pe_p",
            "Pe_bed",
            "c_norm",
            "regen_flow_rate",
            "service_to_regen_flow_ratio",
            "regen_dose",
            "regen_sg",
            "regen_density",
            "regen_ww",
            "regen_conc",
            "regen_vol_per_bv",
            "regen_flow",
            "t_regen",
            "bw_rate",
            "bw_flow",
            "t_bw",
            "bed_expansion_frac",
            "bed_expansion_h",
            "rinse_bv",
            "rinse_flow",
            "t_rinse",
            "main_pump_power",
            "regen_pump_power",
            "bw_pump_power",
            "rinse_pump_power",
            "pump_efficiency",
        ]

        for v in ix_vars:
            assert hasattr(ix, v)
            var = getattr(ix, v)
            assert isinstance(var, Var)

        # # test state block objects
        stateblock_lst = [
            "properties_in",
            "properties_out",
            "properties_regen",
        ]

        for sb_str in stateblock_lst:
            sb = getattr(ix, sb_str)
            assert isinstance(sb, DSPMDEStateBlock)

        # test statistics
        assert number_variables(m) == 128
        assert number_total_constraints(m) == 87
        assert number_unused_variables(m) == 13

    @pytest.mark.unit
    def test_dof(self, IX_frame):
        m = IX_frame
        check_dof(m, fail_flag=True)

    @pytest.mark.unit
    def test_calculate_scaling(self, IX_frame):
        m = IX_frame
        m = ix_scaling(m)

        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m.fs.unit))
        assert len(unscaled_var_list) == 0

    # @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_initialize(self, IX_frame):
        m = IX_frame
        initialization_tester(m)

    @pytest.mark.component
    def test_solve(self, IX_frame):
        m = IX_frame
        calculate_scaling_factors(m)
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_conservation(self, IX_frame):
        m = IX_frame
        ix = m.fs.unit
        comp_lst = m.fs.properties.solute_set

        flow_mass_inlet = sum(
            ix.inlet.properties_in[0].flow_mass_phase_comp["Liq", j] for j in comp_lst
        )
        flow_mass_out = sum(
            ix.outlet.properties_out[0].flow_mass_phase_comp["Liq", j] for j in comp_lst
        )
        flow_mass_waste = sum(ix.waste.flow_mass_phase_comp["Liq", j] for j in comp_lst)

        assert abs(value(flow_mass_inlet - flow_mass_out - flow_mass_waste)) <= 1e-6

    @pytest.mark.component
    def test_solution(self, IX_frame):
        m = IX_frame
        ix = m.fs.unit

        results_dict = {
            "resin_eq_capacity": 3.452724142139907,
            "separation_factor": 0.33333333333333304,
            "bed_vol": 1.8000000000001537,
            "bed_depth": 1.792345170757359,
            "partition_ratio": 241.69068996276124,
            "fluid_mass_transfer_coeff": 3.671987345548717e-05,
            "rate_coeff": 0.00017485654026422467,
            "mass_in": 4359.4324190056495,
            "mass_out": 8.999999908993056,
            "mass_removed": 4350.432419096656,
            "vel_bed": 0.004978736585437118,
            "vel_inter": 0.009957473170876306,
            "service_flow_rate": 10,
            "Re": 4.480862926893026,
            "Sc": 1086.956521739223,
            "Sh": 35.921615336889616,
            "col_diam": 0.5653930831884962,
            "col_vol_tot": 3.6393808454643657,
            "col_vol_per": 0.9098452113660914,
            "number_columns": 4,
            "t_breakthru": 87188.64838478691,
            "t_contact": 179.99999999996228,
            "t_waste": 3299.9999999998113,
            "t_regen": 1800,
            "t_rinse": 899.9999999998113,
        }

        for v, val in results_dict.items():
            var = getattr(ix, v)
            if var.is_indexed():
                assert pytest.approx(val, rel=5e-2) == value(var[target_ion])
            else:
                assert pytest.approx(val, rel=5e-2) == value(var)
