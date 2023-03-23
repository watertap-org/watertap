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
import pyomo.environ as pyo
from pyomo.environ import (
    ConcreteModel,
    check_optimal_termination,
    value,
)
from pyomo.network import Port
from pyomo.util.check_units import assert_units_consistent

from idaes.core import (
    FlowsheetBlock,
    EnergyBalanceType,
    MaterialBalanceType,
    MomentumBalanceType,
)
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    badly_scaled_var_generator,
)
from idaes.core.util.exceptions import ConfigurationError
from idaes.core import UnitModelCostingBlock

from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
    DiffusivityCalculation,
)
from watertap.unit_models.electrolyzer import Electrolyzer
from watertap.costing import WaterTAPCosting

__author__ = "Hunter Barber"

solver = get_solver()

# inputs for badly_scaled_var_generator used across test frames
sv_large = 1e2
sv_small = 1e-2
sv_zero = 1e-8

# -----------------------------------------------------------------------------
class TestElectrolyzer:
    @pytest.fixture(scope="class")
    def electrolyzer_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = MCASParameterBlock(
            solute_list=[
                "H+",
                "OH-",
                "NA+",
                "CL-",
                "CL2",
                "H2",
            ],
            mw_data={
                "H+": 0.001008,
                "OH-": 0.017007,
                "NA+": 0.022989,
                "CL-": 0.03545,
                "CL2": 0.0709,
                "H2": 0.002016,
            },
        )
        m.fs.unit = Electrolyzer(
            property_package=m.fs.properties,
        )

        # feed specifications
        anode_blk = m.fs.unit.anode
        cathode_blk = m.fs.unit.cathode
        anode_blk.properties_in[0].pressure.fix(101325)
        anode_blk.properties_in[0].temperature.fix(273.15 + 25)
        anode_blk.properties_in[0].flow_mol_phase_comp["Liq", "H2O"].fix()
        anode_blk.properties_in[0].flow_mol_phase_comp["Liq", "H+"].fix()
        anode_blk.properties_in[0].flow_mol_phase_comp["Liq", "OH-"].fix()
        anode_blk.properties_in[0].flow_mol_phase_comp["Liq", "NA+"].fix()
        anode_blk.properties_in[0].flow_mol_phase_comp["Liq", "CL-"].fix()
        anode_blk.properties_in[0].flow_mol_phase_comp["Liq", "CL2"].fix()
        anode_blk.properties_in[0].flow_mol_phase_comp["Liq", "H2"].fix()
        cathode_blk.properties_in[0].pressure.fix(101325)
        cathode_blk.properties_in[0].temperature.fix(273.15 + 25)
        cathode_blk.properties_in[0].flow_mol_phase_comp["Liq", "H2O"].fix()
        cathode_blk.properties_in[0].flow_mol_phase_comp["Liq", "H+"].fix()
        cathode_blk.properties_in[0].flow_mol_phase_comp["Liq", "OH-"].fix()
        cathode_blk.properties_in[0].flow_mol_phase_comp["Liq", "NA+"].fix()
        cathode_blk.properties_in[0].flow_mol_phase_comp["Liq", "CL-"].fix()
        cathode_blk.properties_in[0].flow_mol_phase_comp["Liq", "CL2"].fix()
        cathode_blk.properties_in[0].flow_mol_phase_comp["Liq", "H2"].fix()

        return m

    @pytest.mark.unit
    def test_config(self, electrolyzer_frame):
        m = electrolyzer_frame

    @pytest.mark.unit
    def test_build(self, electrolyzer_frame):
        m = electrolyzer_frame

    @pytest.mark.unit
    def test_dof(self, electrolyzer_frame):
        m = electrolyzer_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_scaling_calc(self, electrolyzer_frame):
        m = electrolyzer_frame

        prop = m.fs.properties
        # prop.set_default_scaling("flow_mol_phase_comp", 1e-4, index=("Liq", "H2O"))
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

    @pytest.mark.component
    def test_initialize(self, electrolyzer_frame):
        initialization_tester(electrolyzer_frame)

    @pytest.mark.component
    def test_scaling_init(self, electrolyzer_frame):
        m = electrolyzer_frame
        badly_scaled_var_lst = list(
            badly_scaled_var_generator(m, large=sv_large, small=sv_small, zero=sv_zero)
        )
        print([(x[0].name, x[1]) for x in badly_scaled_var_lst])
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_solve(self, electrolyzer_frame):
        m = electrolyzer_frame

        results = solver.solve(m)

        # Check for optimal solution
        assert check_optimal_termination(results)
