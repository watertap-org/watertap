import pytest
from pyomo.environ import (
    ConcreteModel,
    value,
    assert_optimal_termination,
)

from idaes.core import (
    FlowsheetBlock,
    UnitModelCostingBlock,
)
import idaes.core.util.scaling as iscale
from idaes.core.util.testing import initialization_tester

from watertap.costing import WaterTAPCosting
from watertap.core.util.initialization import check_dof
from watertap.core.solvers import get_solver
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from watertap.unit_models import IonExchangeCP
from watertap.unit_models.tests.unit_test_harness import UnitTestHarness

__author__ = "Kurban Sitterley"

solver = get_solver()
zero = 1e-8
relative_tolerance = 1e-3


def build_constant_pattern():

    target_component = "Ca_2+"
    ion_props = {
        "solute_list": [target_component],
        "diffusivity_data": {("Liq", target_component): 9.2e-10},
        "mw_data": {"H2O": 0.018, target_component: 0.04},
        "charge": {target_component: 2},
    }
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(**ion_props)
    ix_config = {
        "property_package": m.fs.properties,
        "target_component": target_component,
    }
    m.fs.unit = ix = IonExchangeCP(**ix_config)

    ix.process_flow.properties_in.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): 0.5,
            ("conc_mass_phase_comp", ("Liq", target_component)): 0.1,
            ("pressure", None): 101325,
            ("temperature", None): 298,
        },
        hold_state=True,
    )

    ix.process_flow.properties_in[0].flow_mass_phase_comp[...]
    ix.process_flow.properties_out[0].flow_mass_phase_comp[...]
    ix.regeneration_stream[0].flow_mass_phase_comp[...]

    # ix.service_flow_rate.fix(15)
    ix.ebct.fix(240)
    ix.loading_rate.fix(0.00708333)
    ix.langmuir[target_component].fix(0.9)
    ix.resin_max_capacity.fix(3)
    ix.bed_depth.fix(1.7)
    # ix.bed_volume_total.fix(120)
    ix.dimensionless_time.fix()
    ix.number_columns.fix(8)
    ix.resin_diam.fix()
    ix.resin_density.fix()
    ix.bed_porosity.fix(0.5)
    # ix.eq_number_columns_redundant.deactivate()
    ix.number_columns_redundant.fix(2)

    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e-5, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 10, index=("Liq", "Ca_2+")
    )

    iscale.set_scaling_factor(
        m.fs.unit.process_flow.properties_out[0].flow_mass_phase_comp["Liq", "Ca_2+"],
        1e4,
    )

    iscale.calculate_scaling_factors(m)

    return m


class TestIXConstantPattern(UnitTestHarness):
    def configure(self):
        m = build_constant_pattern()
        self.default_zero = zero
        self.default_relative_tolerance = relative_tolerance

        self.unit_solutions[m.fs.unit.c_norm["Ca_2+"]] = 0.4919
        self.unit_solutions[m.fs.unit.column_height] = 3.488
        self.unit_solutions[m.fs.unit.bed_diameter] = 3.351
        self.unit_solutions[m.fs.unit.number_columns] = 8
        self.unit_solutions[m.fs.unit.breakthrough_time] = 52360
        self.unit_solutions[m.fs.unit.loading_rate] = 0.007083
        self.unit_solutions[m.fs.unit.N_Re] = 4.958
        self.unit_solutions[m.fs.unit.N_Sc["Ca_2+"]] = 1086
        self.unit_solutions[m.fs.unit.N_Sh["Ca_2+"]] = 26.296
        self.unit_solutions[m.fs.unit.N_Pe_particle] = 0.10782
        self.unit_solutions[m.fs.unit.N_Pe_bed] = 261.8
        self.unit_solutions[m.fs.unit.resin_eq_capacity] = 1.554
        self.unit_solutions[m.fs.unit.resin_unused_capacity] = 1.445
        self.unit_solutions[m.fs.unit.langmuir["Ca_2+"]] = 0.9
        self.unit_solutions[m.fs.unit.mass_removed["Ca_2+"]] = 65300
        self.unit_solutions[m.fs.unit.num_transfer_units] = 35.5483
        self.unit_solutions[m.fs.unit.dimensionless_time] = 1
        self.unit_solutions[m.fs.unit.partition_ratio] = 217.669
        self.unit_solutions[m.fs.unit.fluid_mass_transfer_coeff["Ca_2+"]] = 3.4561e-05
        self.unit_solutions[m.fs.unit.bw_flow] = 0.09803
        self.unit_solutions[m.fs.unit.bed_expansion_frac] = 0.46395
        self.unit_solutions[m.fs.unit.column_volume_total] = 246.26
        self.unit_solutions[m.fs.unit.lh] = 0.0
        self.unit_solutions[m.fs.unit.separation_factor["Ca_2+"]] = 1.111
        self.unit_solutions[m.fs.unit.rate_coeff["Ca_2+"]] = 0.00021159
        self.unit_solutions[m.fs.unit.HTU["Ca_2+"]] = 0.04782
        self.unit_solutions[
            m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp["Liq", "H2O"]
        ] = 499.95
        self.unit_solutions[
            m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                "Liq", "Ca_2+"
            ]
        ] = 0.05
        self.unit_solutions[
            m.fs.unit.process_flow.properties_out[0.0].flow_mass_phase_comp[
                "Liq", "H2O"
            ]
        ] = 499.95
        self.unit_solutions[
            m.fs.unit.process_flow.properties_out[0.0].flow_mass_phase_comp[
                "Liq", "Ca_2+"
            ]
        ] = 1.1458e-4
        self.unit_solutions[
            m.fs.unit.regeneration_stream[0.0].flow_mass_phase_comp["Liq", "Ca_2+"]
        ] = 0.049885
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "Ca_2+"]
        ] = -1.24713

        return m

    @pytest.mark.component
    def test_costing(self):

        m = build_constant_pattern()
        ix = m.fs.unit

        m.fs.costing = WaterTAPCosting()
        ix.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(ix.process_flow.properties_out[0].flow_vol_phase["Liq"])
        m.fs.costing.add_specific_energy_consumption(
            ix.process_flow.properties_out[0].flow_vol_phase["Liq"], name="SEC"
        )

        check_dof(m, fail_flag=True)
        initialization_tester(m)
        results = solver.solve(m, tee=True)
        assert_optimal_termination(results)

        sys_cost_results = {
            "aggregate_capital_cost": 4359307.7,
            "aggregate_fixed_operating_cost": 40992.5,
            "aggregate_variable_operating_cost": 0.0,
            "aggregate_flow_electricity": 39.4984,
            "aggregate_flow_NaCl": 25376620.7345,
            "aggregate_flow_costs": {"electricity": 24237.0, "NaCl": 2310328.0},
            "total_capital_cost": 4359307.7,
            "total_operating_cost": 2272880.4,
            "LCOW": 0.19076,
            "SEC": 0.02194,
        }

        for v, r in sys_cost_results.items():
            mv = getattr(m.fs.costing, v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)

        ix_cost_results = {
            "capital_cost": 4359307.7327,
            "fixed_operating_cost": 40992.5908,
            "capital_cost_vessel": 101131.903,
            "capital_cost_resin": 81985.1816,
            "capital_cost_regen_tank": 215778.2617,
            "capital_cost_backwash_tank": 132704.7583,
            "flow_mass_regen_soln": 25376620.7345,
            "total_pumping_power": 39.4984,
            "regeneration_tank_vol": 79251.6157,
            "backwash_tank_vol": 174042.7712,
            "direct_capital_cost": 2179653.8663,
        }

        for v, r in ix_cost_results.items():
            mv = getattr(m.fs.unit.costing, v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)
