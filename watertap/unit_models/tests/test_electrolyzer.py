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
from watertap.unit_models.electrolyzer import (
    Electrolyzer,
)
from watertap.costing import WaterTAPCosting

__author__ = "Hunter Barber"

solver = get_solver()

# -----------------------------------------------------------------------------
class TestElectrolyzer:
    @pytest.fixture(scope="class")
    def chlor_alkali_elec(self):

        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = MCASParameterBlock(
            solute_list=[
                "NA+",
                "CL-",
                "CL2-v",
                "H2",
                "OH-",
            ],
            mw_data={
                "H2O": 0.018015,
                "NA+": 0.022989,
                "CL-": 0.3545,
                "CL2-v": 0.0709,
                "H2": 2.016,
                "OH-": 0.017007,
            },
            charge={
                "NA+": 1,
                "CL-": -1,
                "OH-": -1,
            },
        )
        m.fs.unit = Electrolyzer(
            property_package=m.fs.properties,
        )

        # feed specifications

        # fix variables

        return m

    @pytest.mark.unit
    def test_dof(self, chlor_alkali_elec):

        m = chlor_alkali_elec

        m.display()

        assert degrees_of_freedom(m) == 0
