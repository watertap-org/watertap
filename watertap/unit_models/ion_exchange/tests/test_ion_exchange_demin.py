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
    MaterialFlowBasis,
)
import idaes.core.util.scaling as iscale
from idaes.core.util.testing import initialization_tester

from watertap.costing import WaterTAPCosting
from watertap.core.util.initialization import check_dof
from watertap.core.solvers import get_solver
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from watertap.unit_models import IonExchangeDemin
from watertap.unit_models.ion_exchange import IonExchangeBaseData
from watertap.unit_models.tests.unit_test_harness import UnitTestHarness

__author__ = "Kurban Sitterley"

solver = get_solver()
zero = 1e-8
relative_tolerance = 1e-3


def build_cation_exchange():
    """
    Adapted from 0.74 MGD exmample from EPA-WBS Cation Exchange model
    """
    ion_props = dict(
        solute_list=["Ca_2+", "Mg_2+", "Inert"],
        material_flow_basis=MaterialFlowBasis.mass,
        diffusivity_data={
            ("Liq", "Ca_2+"): 7.92e-10,
            ("Liq", "Mg_2+"): 7.06e-10,
            ("Liq", "Inert"): 7.06e-10,
        },
        mw_data={"H2O": 0.018, "Ca_2+": 0.04, "Mg_2+": 0.024, "Inert": 0.10},
        stokes_radius_data={"Ca_2+": 3.09e-10, "Mg_2+": 3.47e-10, "Inert": 1e-10},
        charge={"Ca_2+": 2, "Mg_2+": 2, "Inert": 0},
    )
    feed_mass_frac = {
        "Ca_2+": 13.5e-6,
        "Mg_2+": 39.94e-6,
        "Inert": 100e-6,
    }
    q = 0.74 * pyunits.Mgallons / pyunits.day
    rho = 1000 * pyunits.kg / pyunits.m**3

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(**ion_props)

    ix_config = {
        "property_package": m.fs.properties,
    }
    m.fs.unit = ix = IonExchangeDemin(**ix_config)
    pf = m.fs.unit.process_flow
    prop_in = pf.properties_in[0]
    prop_in.pressure.fix()
    prop_in.temperature.fix()

    flow_mass_water = pyunits.convert((q * rho), to_units=pyunits.kg / pyunits.s)
    prop_in.flow_mass_phase_comp["Liq", "H2O"].fix(flow_mass_water)
    for ion, x in feed_mass_frac.items():
        mass_comp_flow = x * pyunits.kg / pyunits.kg * flow_mass_water
        prop_in.flow_mass_phase_comp["Liq", ion].fix(mass_comp_flow)

    loading_rate = 9.08755 * pyunits.gallon / pyunits.minute / pyunits.ft**2
    loading_rate = value(pyunits.convert(loading_rate, to_units=pyunits.m / pyunits.s))

    ix.loading_rate.fix(loading_rate)
    # ix.ebct.fix(150)
    ix.bed_depth.fix(0.92569949427)
    ix.bed_diameter.fix(1.82894037)
    ix.resin_capacity_ax.fix(1.2366)
    ix.resin_capacity_cx.fix(1.2366)
    ix.resin_diam.fix()
    ix.resin_density.fix()
    ix.bed_porosity.fix()
    ix.number_columns_redundant.fix(1)

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e6, index=("Liq", "Ca_2+")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e6, index=("Liq", "Mg_2+")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e6, index=("Liq", "Inert")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 0.1, index=("Liq", "H2O")
    )

    iscale.calculate_scaling_factors(m)

    return m


class TestIXDeminCation(UnitTestHarness):
    def configure(self):
        m = build_cation_exchange()

        self.default_zero = zero
        self.default_relative_tolerance = relative_tolerance

        self.unit_solutions[m.fs.unit.resin_capacity_op] = 1.2366
        self.unit_solutions[m.fs.unit.bv] = 308.97
        self.unit_solutions[m.fs.unit.ebct] = 150
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "H2O"]
        ] = 0
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "Ca_2+"]
        ] = -0.000437644
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "Mg_2+"]
        ] = -0.001294779
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "Inert"]
        ] = 0
        self.unit_solutions[m.fs.unit.column_height] = 2.3551
        self.unit_solutions[m.fs.unit.number_columns] = 2.000
        self.unit_solutions[m.fs.unit.breakthrough_time] = 46345.6
        self.unit_solutions[m.fs.unit.resin_capacity_op] = 1.2366
        self.unit_solutions[
            m.fs.unit.process_flow.properties_in[0.0].total_hardness
        ] = 200.30

        return m

    @pytest.mark.component
    def test_costing(self):
        m = build_cation_exchange()
        ix = m.fs.unit

        m.fs.costing = WaterTAPCosting()
        ix.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        ix.costing.regen_dose.set_value(240)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(ix.process_flow.properties_out[0].flow_vol_phase["Liq"])
        m.fs.costing.add_specific_energy_consumption(
            ix.process_flow.properties_out[0].flow_vol_phase["Liq"], name="SEC"
        )

        sys_cost_results = {
            "aggregate_capital_cost": 483892.9,
            "aggregate_fixed_operating_cost": 1993.8,
            "aggregate_variable_operating_cost": 0.0,
            "aggregate_flow_electricity": 1.2152,
            "aggregate_flow_NaCl": 1116422.5,
            "aggregate_flow_costs": {"electricity": 745.72, "NaCl": 101640.8},
            "total_capital_cost": 483892.9,
            "total_operating_cost": 108658.5,
            "LCOW": 0.170524,
        }

        check_dof(m, fail_flag=True)
        initialization_tester(m)
        results = solver.solve(m)
        assert_optimal_termination(results)

        for v, r in sys_cost_results.items():
            mv = getattr(m.fs.costing, v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)

        ix_cost_results = {
            "capital_cost": 483892.941,
            "capital_cost_vessel": 48385.243,
            "capital_cost_resin": 13292.384,
            "capital_cost_regen_tank": 29343.804,
            "capital_cost_backwash_tank": 27569.78,
            "flow_mass_regen_soln": 1116422.518,
            "total_pumping_power": 1.215287,
            "regeneration_tank_vol": 5139.677,
            "backwash_tank_vol": 7581.307,
            "direct_capital_cost": 241946.47,
        }

        for v, r in ix_cost_results.items():
            mv = getattr(m.fs.unit.costing, v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)