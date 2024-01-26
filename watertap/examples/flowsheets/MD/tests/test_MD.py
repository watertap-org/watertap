#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
import pytest
from pyomo.environ import (
    ConcreteModel,
    value,
    Var,
)
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.util.check_units import assert_units_consistent
from idaes.models.unit_models import Heater, Separator, Mixer, Product, Feed

from idaes.models.unit_models.heat_exchanger import (
    HeatExchanger,
)
from watertap.unit_models.pressure_changer import Pump
from watertap.unit_models.MD.membrane_distillation_0D import MembraneDistillation0D
import watertap.property_models.seawater_prop_pack as props_sw
import watertap.property_models.water_prop_pack as props_w
from watertap.examples.flowsheets.MD.MD_single_stage_continuous_recirculation import (
    main,
    build,
    set_operating_conditions,
    optimize_set_up,
    initialize_system,
    solve,
)

solver = get_solver()

# -----------------------------------------------------------------------------
class TestMDContinuousRecirculation:
    @pytest.fixture(scope="class")
    def system_frame(self):
        m = build()

        return m

    @pytest.mark.unit
    def test_build(self, system_frame):
        m = system_frame
        # model set up
        assert isinstance(m, ConcreteModel)
        assert isinstance(m.fs, FlowsheetBlock)
        assert isinstance(m.fs.properties_hot_ch, props_sw.SeawaterParameterBlock)
        assert isinstance(m.fs.properties_cold_ch, props_sw.SeawaterParameterBlock)
        assert isinstance(m.fs.properties_vapor, props_w.WaterParameterBlock)

        # unit models
        assert isinstance(m.fs.feed, Feed)
        assert isinstance(m.fs.pump_feed, Pump)
        assert isinstance(m.fs.pump_permeate, Pump)
        assert isinstance(m.fs.pump_brine, Pump)
        assert isinstance(m.fs.mixer, Mixer)
        assert isinstance(m.fs.heater, Heater)
        assert isinstance(m.fs.chiller, Heater)
        assert isinstance(m.fs.hx, HeatExchanger)
        assert isinstance(m.fs.MD, MembraneDistillation0D)
        assert isinstance(m.fs.separator_concentrate, Separator)
        assert isinstance(m.fs.separator_permeate, Separator)
        assert isinstance(m.fs.permeate, Product)

        # additional variables
        assert isinstance(m.fs.overall_recovery, Var)
        assert isinstance(m.fs.recycle_ratio, Var)

        # arcs
        arc_dict = {
            # brine (MD hot side) loop
            m.fs.s01: (m.fs.feed.outlet, m.fs.pump_feed.inlet),
            m.fs.s02: (m.fs.pump_feed.outlet, m.fs.mixer.feed),
            m.fs.s03: (m.fs.mixer.outlet, m.fs.hx.cold_inlet),
            m.fs.s04: (m.fs.hx.cold_outlet, m.fs.pump_brine.inlet),
            m.fs.s05: (m.fs.pump_brine.outlet, m.fs.heater.inlet),
            m.fs.s06: (m.fs.heater.outlet, m.fs.MD.hot_ch_inlet),
            m.fs.s07: (m.fs.MD.hot_ch_outlet, m.fs.separator_concentrate.inlet),
            m.fs.s08: (m.fs.separator_concentrate.reject, m.fs.reject.inlet),
            m.fs.s09: (m.fs.separator_concentrate.recycle, m.fs.mixer.recycle),
            # pemeate (MD cold side) loop
            m.fs.s10: (m.fs.chiller.outlet, m.fs.MD.cold_ch_inlet),
            m.fs.s11: (m.fs.MD.cold_ch_outlet, m.fs.hx.hot_inlet),
            m.fs.s12: (m.fs.hx.hot_outlet, m.fs.separator_permeate.inlet),
            m.fs.s13: (m.fs.separator_permeate.permeate, m.fs.permeate.inlet),
            m.fs.s14: (
                m.fs.separator_permeate.cold_loop_stream,
                m.fs.pump_permeate.inlet,
            ),
            m.fs.s15: (m.fs.pump_permeate.outlet, m.fs.chiller.inlet),
        }
        for arc, port_tpl in arc_dict.items():
            assert arc.source is port_tpl[0]
            assert arc.destination is port_tpl[1]

    @pytest.mark.component
    def test_units(self, system_frame):
        assert_units_consistent(system_frame)

    @pytest.mark.component
    def test_set_operating_conditions(self, system_frame):
        m = system_frame
        set_operating_conditions(m)
        assert degrees_of_freedom(m) == 0

    @pytest.mark.component
    @pytest.mark.requires_idaes_solver
    def test_initialize_system(self, system_frame):
        m = system_frame
        initialize_system(m, solver=solver)
        assert degrees_of_freedom(m) == 0

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_solution(self, system_frame):
        m = system_frame
        optimize_set_up(m)
        solve(m, solver=solver)
        assert pytest.approx(0, rel=1e-5) == value(
            m.fs.permeate.flow_mass_phase_comp[0, "Liq", "TDS"]
        )
        assert pytest.approx(0.5, rel=1e-5) == value(m.fs.overall_recovery)
        assert pytest.approx(6.169, rel=1e-3) == value(m.fs.recycle_ratio[0])
        assert pytest.approx(15.44, rel=1e-3) == value(m.fs.costing.LCOW)
        assert pytest.approx(182.9, rel=1e-3) == value(
            m.fs.costing.specific_energy_consumption
        )

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_main(self):
        main()
