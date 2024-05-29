import pytest
from pyomo.environ import (
    ConcreteModel,
    value,
    assert_optimal_termination,
    units as pyunits,
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
from watertap.unit_models import IonExchangeClark
from watertap.unit_models.tests.unit_test_harness import UnitTestHarness

__author__ = "Kurban Sitterley"

solver = get_solver()
zero = 1e-8
relative_tolerance = 1e-3


def build_freundlich():

    target_component = "Cl_-"

    ion_props = {
        "solute_list": [target_component],
        "diffusivity_data": {("Liq", target_component): 1e-9},
        "mw_data": {"H2O": 0.018, target_component: 35.45e-3},
        "charge": {target_component: -1},
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(**ion_props)
    ix_config = {
        "property_package": m.fs.properties,
        "target_component": target_component,
        "isotherm": "freundlich",
        "regenerant": "NaOH",
        "hazardous_waste": True,
    }
    m.fs.unit = ix = IonExchangeClark(**ix_config)

    ix.process_flow.properties_in.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): 0.5,
            ("conc_mass_phase_comp", ("Liq", target_component)): 1e-6,
            ("pressure", None): 101325,
            ("temperature", None): 298,
        },
        hold_state=True,
    )

    ix.process_flow.properties_in[0].flow_mass_phase_comp[...]
    ix.process_flow.properties_out[0].flow_mass_phase_comp[...]
    ix.regeneration_stream[0].flow_mass_phase_comp[...]

    ix.freundlich_n.fix(1.2)
    ix.bv_50.fix(20000)
    ix.bv.fix(18000)
    ix.resin_density.fix(0.72)
    ix.bed_porosity.fix(0.5)
    ix.loading_rate.fix(6.15e-3)
    ix.resin_diam.fix(6.75e-4)
    ix.number_columns.fix(16)
    ix.c_norm.fix(0.25)
    ix.number_columns_redundant.fix(4)
    ix.bed_depth.fix(1.476)
    ix.ebct.fix(240)

    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e-4, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e6, index=("Liq", "Cl_-")
    )
    iscale.calculate_scaling_factors(m)

    return m


class TestIXFreundlich(UnitTestHarness):
    def configure(self):
        m = build_freundlich()

        self.default_zero = zero
        self.default_relative_tolerance = relative_tolerance
        self.unit_solutions[m.fs.unit.column_height] = 3.160
        self.unit_solutions[m.fs.unit.bed_diameter] = 2.5435
        self.unit_solutions[m.fs.unit.N_Sh["Cl_-"]] = 24.083
        self.unit_solutions[m.fs.unit.N_Pe_particle] = 0.09901
        self.unit_solutions[m.fs.unit.N_Pe_bed] = 216.51
        self.unit_solutions[m.fs.unit.tb_traps[0]] = 0
        self.unit_solutions[m.fs.unit.tb_traps[1]] = 3344557
        self.unit_solutions[m.fs.unit.tb_traps[2]] = 3825939
        self.unit_solutions[m.fs.unit.tb_traps[3]] = 4034117
        self.unit_solutions[m.fs.unit.tb_traps[4]] = 4188551
        self.unit_solutions[m.fs.unit.tb_traps[5]] = 4320000
        self.unit_solutions[m.fs.unit.traps[1]] = 0.0038710157
        self.unit_solutions[m.fs.unit.traps[2]] = 0.0044572367
        self.unit_solutions[m.fs.unit.traps[3]] = 0.0048189406
        self.unit_solutions[m.fs.unit.traps[4]] = 0.0057197676
        self.unit_solutions[m.fs.unit.traps[5]] = 0.0066941563
        self.unit_solutions[m.fs.unit.c_norm_avg["Cl_-"]] = 0.02556
        self.unit_solutions[m.fs.unit.mass_transfer_coeff] = 0.159346
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "Cl_-"]
        ] = -1.3744e-05
        self.unit_solutions[
            m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp["Liq", "H2O"]
        ] = 500
        self.unit_solutions[
            m.fs.unit.process_flow.properties_out[0.0].flow_mass_phase_comp[
                "Liq", "H2O"
            ]
        ] = 500
        self.unit_solutions[
            m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                "Liq", "Cl_-"
            ]
        ] = 5e-07
        self.unit_solutions[
            m.fs.unit.process_flow.properties_out[0.0].flow_mass_phase_comp[
                "Liq", "Cl_-"
            ]
        ] = 1.278e-08
        self.unit_solutions[
            m.fs.unit.regeneration_stream[0.0].flow_mass_phase_comp["Liq", "Cl_-"]
        ] = 4.872e-07

        return m

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_solution(self):
        m = build_freundlich()
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
            "aggregate_capital_cost": 5895763.7,
            "aggregate_fixed_operating_cost": 379781.6,
            "aggregate_variable_operating_cost": 0,
            "aggregate_flow_electricity": 30.805,
            "aggregate_flow_NaOH": 328451.29,
            "aggregate_flow_costs": {"electricity": 18902.5, "NaOH": 653430.0},
            "total_capital_cost": 5895763.7,
            "total_operating_cost": 1161753.8,
            "aggregate_direct_capital_cost": 2947881.8,
            "maintenance_labor_chemical_operating_cost": 176872.9,
            "total_fixed_operating_cost": 556654.61,
            "total_variable_operating_cost": 605099.2,
            "total_annualized_cost": 1751330.2,
            "LCOW": 0.12332,
            "SEC": 0.01711,
        }

        for v, r in sys_cost_results.items():
            mv = getattr(m.fs.costing, v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)

        ix_cost_results = {
            "capital_cost": 5895763.7,
            "fixed_operating_cost": 379781.6,
            "capital_cost_vessel": 75000.3,
            "capital_cost_resin": 54924.687,
            "capital_cost_regen_tank": 215778.2,
            "capital_cost_backwash_tank": 133603.4,
            "operating_cost_hazardous": 324857.0,
            "flow_mass_regen_soln": 328451.2,
            "total_pumping_power": 30.8049,
            "backwash_tank_vol": 176401.1,
            "regeneration_tank_vol": 79251.6,
            "direct_capital_cost": 2947881.8,
        }

        for v, r in ix_cost_results.items():
            mv = getattr(m.fs.unit.costing, v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)


def build_inert():

    target_component = "Cl_-"
    inert = "Ca_2+"

    ion_props = {
        "solute_list": [target_component, inert],
        "diffusivity_data": {("Liq", target_component): 1e-9, ("Liq", inert): 9.2e-10},
        "mw_data": {"H2O": 0.018, target_component: 35.45e-3, inert: 40e-3},
        "charge": {target_component: -1, inert: 2},
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(**ion_props)
    ix_config = {
        "property_package": m.fs.properties,
        "target_component": target_component,
        "regenerant": "single_use",
        "isotherm": "freundlich",
    }
    m.fs.unit = ix = IonExchangeClark(**ix_config)

    ix.process_flow.properties_in.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): 0.5,
            ("conc_mass_phase_comp", ("Liq", target_component)): 1e-6,
            ("conc_mass_phase_comp", ("Liq", inert)): 1e-7,
            ("pressure", None): 101325,
            ("temperature", None): 298,
        },
        hold_state=True,
    )

    ix.process_flow.properties_in[0].flow_mass_phase_comp[...]
    ix.process_flow.properties_out[0].flow_mass_phase_comp[...]
    ix.regeneration_stream[0].flow_mass_phase_comp[...]

    ix.freundlich_n.fix(1.2)
    ix.bv_50.fix(20000)
    ix.bv.fix(18000)
    ix.resin_density.fix(0.72)
    ix.bed_porosity.fix(0.5)
    ix.loading_rate.fix(6.15e-3)
    ix.resin_diam.fix(6.75e-4)
    ix.number_columns.fix(16)
    ix.c_norm.fix(0.25)
    ix.number_columns_redundant.fix(4)
    ix.bed_depth.fix(1.476)
    ix.ebct.fix(240)

    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e-4, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e5, index=("Liq", "Cl_-")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e5, index=("Liq", "Ca_2+")
    )
    iscale.calculate_scaling_factors(m)

    return m


class TestIXInert(UnitTestHarness):
    def configure(self):
        m = build_inert()

        self.default_zero = zero
        self.default_relative_tolerance = relative_tolerance

        self.unit_solutions[m.fs.unit.number_columns] = 16
        self.unit_solutions[m.fs.unit.column_height] = 3.160
        self.unit_solutions[m.fs.unit.bed_diameter] = 2.5435
        self.unit_solutions[m.fs.unit.N_Sh["Cl_-"]] = 24.083
        self.unit_solutions[m.fs.unit.N_Pe_particle] = 0.09901
        self.unit_solutions[m.fs.unit.N_Pe_bed] = 216.51
        self.unit_solutions[m.fs.unit.tb_traps[1]] = 3344557
        self.unit_solutions[m.fs.unit.tb_traps[2]] = 3825939
        self.unit_solutions[m.fs.unit.tb_traps[3]] = 4034117
        self.unit_solutions[m.fs.unit.tb_traps[4]] = 4188551
        self.unit_solutions[m.fs.unit.tb_traps[5]] = 4320000
        self.unit_solutions[m.fs.unit.c_traps[2]] = 0.07
        self.unit_solutions[m.fs.unit.c_traps[3]] = 0.13
        self.unit_solutions[m.fs.unit.c_traps[4]] = 0.19
        self.unit_solutions[m.fs.unit.c_traps[5]] = 0.25
        self.unit_solutions[m.fs.unit.traps[4]] = 0.0057197676
        self.unit_solutions[m.fs.unit.traps[5]] = 0.0066941563
        self.unit_solutions[m.fs.unit.c_norm_avg["Cl_-"]] = 0.02556
        self.unit_solutions[m.fs.unit.mass_transfer_coeff] = 0.159346
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "Cl_-"]
        ] = -1.37435e-05
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "Ca_2+"]
        ] = 0
        self.unit_solutions[
            m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp["Liq", "H2O"]
        ] = 500
        self.unit_solutions[
            m.fs.unit.process_flow.properties_out[0.0].flow_mass_phase_comp[
                "Liq", "H2O"
            ]
        ] = 500
        self.unit_solutions[
            m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                "Liq", "Cl_-"
            ]
        ] = 5e-07
        self.unit_solutions[
            m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                "Liq", "Ca_2+"
            ]
        ] = 5e-08
        self.unit_solutions[
            m.fs.unit.process_flow.properties_out[0.0].flow_mass_phase_comp[
                "Liq", "Cl_-"
            ]
        ] = 1.2781e-08
        self.unit_solutions[
            m.fs.unit.process_flow.properties_out[0.0].flow_mass_phase_comp[
                "Liq", "Ca_2+"
            ]
        ] = 5e-08
        self.unit_solutions[
            m.fs.unit.regeneration_stream[0.0].flow_mass_phase_comp["Liq", "Cl_-"]
        ] = 4.872e-07

        return m

    @pytest.mark.component
    def test_costing(self):

        m = build_inert()
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
            "aggregate_capital_cost": 5464207.216,
            "aggregate_fixed_operating_cost": 6419597.4,
            "aggregate_variable_operating_cost": 0.0,
            "aggregate_flow_electricity": 30.813,
            "aggregate_flow_costs": {"electricity": 18907.7},
            "total_capital_cost": 5464207.2,
            "total_operating_cost": 6600540.6,
            "aggregate_direct_capital_cost": 2732103.6,
            "maintenance_labor_chemical_operating_cost": 163926.2,
            "total_fixed_operating_cost": 6583523.67,
            "total_variable_operating_cost": 17016.9,
            "total_annualized_cost": 7146961.3,
            "LCOW": 0.50327,
            "SEC": 0.017118,
        }

        for v, r in sys_cost_results.items():
            mv = getattr(m.fs.costing, v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)

        ix_cost_results = {
            "capital_cost": 5464207.216,
            "fixed_operating_cost": 6419597.454,
            "capital_cost_vessel": 75000.32,
            "capital_cost_resin": 54924.687,
            "capital_cost_backwash_tank": 133603.452,
            "total_pumping_power": 30.813,
            "flow_vol_resin": 876.6,
            "single_use_resin_replacement_cost": 6419597.454,
            "backwash_tank_vol": 176401.066,
            "direct_capital_cost": 2732103.608,
        }

        for v, r in ix_cost_results.items():
            mv = getattr(m.fs.unit.costing, v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)
