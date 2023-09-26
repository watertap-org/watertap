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
    Var,
    value,
)
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.util.check_units import assert_units_consistent
import watertap.property_models.NaCl_prop_pack as props
from watertap.examples.flowsheets.oaro.oaro_multi import (
    main,
    build,
    set_operating_conditions,
    initialize_system,
    solve,
    ERDtype,
)

solver = get_solver()

# -----------------------------------------------------------------------------
class TestOAROwithTurbine:
    @pytest.fixture(scope="class")
    def system_frame(self):
        m = build(number_of_stages=3, erd_type=ERDtype.pump_as_turbine)

        return m

    @pytest.mark.unit
    def test_build(self, system_frame):
        m = system_frame
        # model set up
        assert isinstance(m, ConcreteModel)
        assert isinstance(m.fs, FlowsheetBlock)
        assert isinstance(m.fs.properties, props.NaClParameterBlock)

        # additional variables
        assert isinstance(m.fs.water_recovery, Var)

    @pytest.mark.component
    def test_units(self, system_frame):
        assert_units_consistent(system_frame)

    @pytest.mark.component
    def test_set_operating_conditions(self, system_frame):
        m = system_frame
        set_operating_conditions(m)
        assert degrees_of_freedom(m) == 0

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_initialize_system(self, system_frame):
        m = system_frame
        initialize_system(
            m,
            number_of_stages=3,
            solvent_multiplier=0.5,
            solute_multiplier=0.7,
            solver=solver,
        )
        assert degrees_of_freedom(m) == 0

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_solution(self, system_frame):
        m = system_frame
        solve(m, solver=solver)
        fs = m.fs
        assert pytest.approx(3.98455e-3, rel=1e-5) == value(
            fs.product.flow_mass_phase_comp[0, "Liq", "NaCl"]
        )
        assert pytest.approx(0.2888, rel=1e-5) == value(fs.water_recovery)

    @pytest.mark.component
    def test_config_error(self, system_frame):
        with pytest.raises(Exception):
            build(number_of_stages=3, erd_type="not_a_configuration")

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_main(self):
        main(number_of_stages=3, system_recovery=0.5)
