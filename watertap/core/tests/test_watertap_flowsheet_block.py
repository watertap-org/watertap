#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

__author__ = "Alexander V. Dudchenko"

from idaes.models.unit_models import Feed, Product

from pyomo.environ import (
    TransformationFactory,
    ConcreteModel,
    Var,
    Constraint,
    units as pyunits,
)
from idaes.core import (
    FlowsheetBlock,
)

from watertap.core.solvers import get_solver
from pyomo.common.config import ConfigValue
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.core.watertap_flowsheet_block import WaterTapFlowsheetBlockData
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import (
    declare_process_block_class,
)
from pyomo.environ import (
    assert_optimal_termination,
)
import pytest

from io import StringIO


@declare_process_block_class("FlowsheetFeed")
class FlowsheetFeedData(WaterTapFlowsheetBlockData):
    """Test class for a simple feed unit within a WaterTapFlowsheetBlock"""

    CONFIG = WaterTapFlowsheetBlockData.CONFIG()
    CONFIG.declare(
        "solute_concentration",
        ConfigValue(
            doc="Total Dissolved Solids (TDS) in the feed water (mg/L)",
            default=35 * pyunits.g / pyunits.L,
        ),
    )
    CONFIG.declare(
        "feed_flow_rate",
        ConfigValue(
            doc="Feed flow rate (m^3/h)",
            default=0.001 * pyunits.m**3 / pyunits.s,
        ),
    )

    def build(self):
        super().build()
        self.feed = Feed(property_package=self.config.default_property_package)
        self.solute_type = list(self.config.default_property_package.solute_set)[0]
        self.feed.ph = Var(initialize=7.0, units=pyunits.dimensionless)
        self.register_port("outlet", self.feed.outlet, {"pH": self.feed.ph})

    def set_fixed_operation(self):
        self.feed.ph.fix(8.5)

        # fix feed flow and concentration based on config arguments
        self.feed.properties[0].flow_vol_phase["Liq"].fix(self.config.feed_flow_rate)
        self.feed.properties[0].conc_mass_phase_comp["Liq", self.solute_type].fix(
            self.config.solute_concentration
        )
        self.feed.properties[0].temperature.fix(298.15)
        self.feed.properties[0].pressure.fix(101325)

        # make sure mass flows are not fixed
        self.feed.properties[0].flow_mass_phase_comp["Liq", self.solute_type].unfix()
        self.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
        assert degrees_of_freedom(self) == 0

        # solve the block so we can get mass flows of solute and water
        solver = get_solver()
        results = solver.solve(self.feed, tee=False)
        # ensure we terminate okay.
        assert_optimal_termination(results)

    def initialize_unit(self, **kwargs):
        """custom initialize routine for the feed unit as we
        are fixing feed concentration and flowrate rather than mass flowrates which will
        throw an error if we just call self.feed.initialize() since feed mass flow and tds flow as unfixed

        This is probably redundant if the model was already "fixed" as result will not change, but
        do't want to assume that user has not change flow/concentration since calling set_fixed_operation (or fix_and_scale)
        """
        solver = get_solver()
        result = solver.solve(self.feed, tee=False)

        assert_optimal_termination(result)

    def scale_before_initialization(self):
        """apply standard scaling factors to feed unit model"""
        self.config.default_property_package.set_default_scaling(
            "flow_mass_phase_comp",
            1 / self.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].value,
            index=("Liq", "H2O"),
        )
        self.config.default_property_package.set_default_scaling(
            "flow_mass_phase_comp",
            1
            / self.feed.properties[0]
            .flow_mass_phase_comp["Liq", self.solute_type]
            .value,
            index=("Liq", self.solute_type),
        )

    def get_unit_name(self):
        return "test feed unit"

    def get_model_state_dict(self):
        model_state = {
            "Composition": {},
            "Physical state": {},
        }
        model_state["Composition"]["Mass flow of H2O"] = self.feed.properties[
            0
        ].flow_mass_phase_comp["Liq", "H2O"]
        model_state["Composition"]["Mass flow of TDS"] = self.feed.properties[
            0
        ].flow_mass_phase_comp["Liq", "TDS"]
        for phase, ion in self.feed.properties[0].conc_mass_phase_comp:
            if ion != "H2O":
                model_state["Composition"][ion] = self.feed.properties[
                    0
                ].conc_mass_phase_comp[phase, ion]
        model_state["Physical state"]["Temperature"] = self.feed.properties[
            0
        ].temperature
        model_state["Physical state"]["Pressure"] = self.feed.properties[0].pressure
        model_state["Physical state"]["Volumetric flowrate"] = self.feed.properties[
            0
        ].flow_vol_phase["Liq"]
        return model_state


@pytest.fixture
def feed_flowsheet_fixture():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.seawater_props = SeawaterParameterBlock()
    m.fs.feed = FlowsheetFeed(
        default_property_package=m.fs.seawater_props,
        solute_concentration=10 * pyunits.kg / pyunits.m**3,
        feed_flow_rate=1000 * pyunits.m**3 / pyunits.s,
    )

    m.fs.product = Product(property_package=m.fs.seawater_props)

    m.fs.feed.outlet.connect_to(m.fs.product.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    m.fs.feed.fix_and_scale()
    m.fs.feed.initialize()
    m.fs.product.initialize()
    return m


def test_feed_flowsheet_initialization(feed_flowsheet_fixture):
    m = feed_flowsheet_fixture

    # should be zero DOF after initialization
    assert degrees_of_freedom(m) == 0

    # this is simpel model, so product should have same mass flow as input
    assert (
        pytest.approx(
            m.fs.product.flow_mass_phase_comp[0, "Liq", "H2O"].value, rel=1e-5
        )
        == m.fs.feed.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].value
    )
    assert (
        pytest.approx(
            m.fs.product.flow_mass_phase_comp[0, "Liq", "TDS"].value, rel=1e-5
        )
        == m.fs.feed.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"].value
    )


def test_report(feed_flowsheet_fixture):
    m = feed_flowsheet_fixture
    os = StringIO()
    m.fs.feed.report(ostream=os)
    result = """
------------------------------------------------------------------------------------
    test feed unit state

    Composition:
    Key              : Value      : Units   : Fixed : Bounds
    Mass flow of H2O : 9.9448e+05 :    kg/s : False : (0.0, None)
    Mass flow of TDS :     10000. :    kg/s : False : (0.0, None)
                 TDS :     10.000 : kg/m**3 :  True : (0.0, 1000000.0)

    Physical state:
    Key                 : Value      : Units  : Fixed : Bounds
                Pressure : 1.0132e+05 :     Pa :  True : (1000.0, 50000000.0)
            Temperature :     298.15 :      K :  True : (273.15, 1000)
    Volumetric flowrate :     1000.0 : m**3/s :  True : (0.0, None)

------------------------------------------------------------------------------------
"""
    print(os.getvalue().replace(" ", ""))
    # testing with out spaces, as they are hard to control in the report output
    assert os.getvalue().replace(" ", "") == result.replace(" ", "")
