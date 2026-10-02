#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
import pytest
from pyomo.environ import value, units as pyunits

from watertap.property_models.NaCl_T_dep_prop_pack import NaClParameterBlock
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.unit_models.reverse_osmosis_0D import ReverseOsmosis0D
from watertap.unit_models.pressure_changer import Pump
import watertap.flowsheets.multistage_RO.multistage_RO as multistage
import watertap.flowsheets.multistage_RO.utils as utils

# Results for salinity range tests for RO-0D
lcow_results_salinity = {
    5: {
        1: {False: 0.3522, True: 0.3328},
        2: {False: 0.349, True: 0.3301},
        3: {False: 0.3482, True: 0.3295},
    },
    100: {
        1: {False: 1.4469, True: 1.4046},
        2: {False: 1.465, True: 1.4281},
        3: {False: 1.4901, True: 1.4614},
    },
}
sec_results_salinity = {
    5: {
        1: {False: 1.1893, True: 0.8037},
        2: {False: 1.177, True: 0.7952},
        3: {False: 1.1742, True: 0.7939},
    },
    100: {
        1: {False: 7.0295, True: 6.2209},
        2: {False: 7.1696, True: 6.3749},
        3: {False: 7.3433, True: 6.6009},
    },
}

# Results for flow range tests for RO-0D
lcow_results_flow = {
    1: {
        1: {False: 0.7308, True: 0.6671},
        2: {False: 0.7368, True: 0.6718},
        3: {False: 0.7434, True: 0.6921},
    },
    1000: {
        1: {False: 0.7815, True: 0.6902},
        2: {False: 0.774, True: 0.6883},
        3: {False: 0.7742, True: 0.689},
    },
}
sec_results_flow = {
    1: {
        1: {False: 2.9621, True: 2.0312},
        2: {False: 2.9927, True: 2.054},
        3: {False: 3.0317, True: 2.2759},
    },
    1000: {
        1: {False: 3.7141, True: 2.418},
        2: {False: 3.5608, True: 2.3535},
        3: {False: 3.562, True: 2.3556},
    },
}

# Results for recovery range tests for RO-0D
lcow_results_recov = {
    0.4: {
        1: {False: 0.8571, True: 0.7457},
        2: {False: 0.8638, True: 0.7483},
        3: {False: 0.8835, True: 0.7513},
    },
    0.5: {
        1: {False: 0.7308, True: 0.6671},
        2: {False: 0.7368, True: 0.6718},
        3: {False: 0.7434, True: 0.6921},
    },
    0.6: {
        1: {False: 0.6429, True: 0.6033},
        2: {False: 0.6482, True: 0.6079},
        3: {False: 0.6546, True: 0.6644},
    },
}
sec_results_recov = {
    0.4: {
        1: {False: 3.5393, True: 2.2241},
        2: {False: 3.5795, True: 2.2414},
        3: {False: 3.9585, True: 2.2641},
    },
    0.5: {
        1: {False: 2.9621, True: 2.0312},
        2: {False: 2.9927, True: 2.054},
        3: {False: 3.0317, True: 2.2759},
    },
    0.6: {
        1: {False: 2.5666, True: 1.9514},
        2: {False: 2.5914, True: 1.9709},
        3: {False: 2.6261, True: 2.1895},
    },
}

salinity = [5, 35, 75]  # g/L
flows = [1, 5, 10, 50, 100, 500, 1000]  # L/s
recovery = [0.4, 0.5, 0.6]
n_stages = [1, 2, 3]
add_erd = [True, False]

# salinity = [35]
# flows = [1]
# recovery = [0.5]
# n_stages = [2]
# add_erd = [True]

default_ro_op_dict = {
    "A_comp": 1.5 * pyunits.liter / pyunits.m**2 / pyunits.hour / pyunits.bar,
    "B_comp": 0.1 * pyunits.liter / pyunits.m**2 / pyunits.hour,
}


class TestMultiStageRO_0D:

    @pytest.mark.unit
    def test_wrong_prop_pack(self):
        with pytest.raises(
            ValueError,
            match="Only NaCl or Seawater property models can be used but MCAS was passed.",
        ):
            _ = multistage.build_n_stage_system(prop_pack="MCAS", RO_1D=False)

    @pytest.mark.unit
    def test_that_ro_param_doesnt_exist(self):
        with pytest.raises(
            ValueError,
            match="Component rainjacket not found in RO unit model",
        ):
            _ = multistage.build_n_stage_system(
                n_stages=3,
                add_costing=False,
                prop_pack="NACL",
                ro_op_dict={"rainjacket": 7983},
            )

    @pytest.mark.unit
    def test_NaCl_build(self):
        m = multistage.build_n_stage_system(
            n_stages=5, add_costing=False, prop_pack="NACL", RO_1D=False
        )
        assert m.fs.properties.solute_set.at(1) == "NaCl"
        assert isinstance(m.fs.properties, NaClParameterBlock)
        assert len(m.fs.stage) == 5
        assert len(m.fs.stages_set) == 5
        assert len(m.fs.product_mixer.config.inlet_list) == len(m.fs.stages_set)

        for n, stage in m.fs.stage.items():
            assert isinstance(stage.RO, ReverseOsmosis0D)
            if not n == m.fs.stages_set.first():
                assert not stage.has_pump
                assert stage.find_component("pump") is None
            else:
                assert stage.has_pump
                assert isinstance(stage.find_component("pump"), Pump)

        assert not hasattr(m.fs, "costing")
        assert hasattr(m.fs, "SEC")
        assert hasattr(m.fs, "SEC_constraint")
        assert hasattr(m.fs, "ERD")

    @pytest.mark.unit
    def test_SW_build(self):
        m = multistage.build_n_stage_system(
            n_stages=3,
            pump_dict={1: True, 3: True},
            prop_pack="seawater",
            add_erd=False,
            RO_1D=False,
        )
        assert m.fs.properties.solute_set.at(1) == "TDS"
        assert isinstance(m.fs.properties, SeawaterParameterBlock)
        assert len(m.fs.stage) == 3
        for n, stage in m.fs.stage.items():
            assert isinstance(stage.RO, ReverseOsmosis0D)
            if n in [1, 3]:
                assert stage.has_pump
                assert isinstance(stage.find_component("pump"), Pump)
            elif n == 2:
                assert not stage.has_pump
                assert stage.find_component("pump") is None

        assert hasattr(m.fs, "costing")
        assert not hasattr(m.fs, "SEC")
        assert not hasattr(m.fs, "SEC_constraint")
        assert not hasattr(m.fs, "ERD")

    @pytest.mark.unit
    def test_reporting(self):
        m = multistage.run_n_stage_system(RO_1D=False)
        utils.report_n_stage_system(m)

    @pytest.mark.parametrize("n", n_stages)
    @pytest.mark.parametrize("erd", add_erd)
    @pytest.mark.parametrize("salt", lcow_results_salinity.keys())
    @pytest.mark.integration
    def test_multistage_ro_salinity_range(self, salt, n, erd):
        """
        Test up to 3-stage with and without ERD, with booster pump for salinity range
        for 1 L/s, 50% recovery
        """
        lcow_results = lcow_results_salinity
        sec_results = sec_results_salinity
        m = multistage.run_n_stage_system(
            n_stages=n,
            salinity=salt,
            add_erd=erd,
            flow_vol=1,
            pump_dict={1: True, 2: True, 3: False},
            ro_op_dict=default_ro_op_dict,
            RO_1D=False,
        )

        m = multistage.set_system_recovery(m, 0.5)

        _ = utils.solve(model=m, tee=False)

        assert pytest.approx(lcow_results[salt][n][erd], rel=1e-3) == value(
            m.fs.costing.LCOW
        )
        assert pytest.approx(sec_results[salt][n][erd], rel=1e-3) == value(
            m.fs.costing.SEC
        )

    @pytest.mark.parametrize("n", n_stages)
    @pytest.mark.parametrize("erd", add_erd)
    @pytest.mark.parametrize("flow", lcow_results_flow.keys())
    @pytest.mark.integration
    def test_multistage_ro_flow_range(self, flow, n, erd):
        """
        Test up to 3-stage with and without ERD, with booster pump for flow range
        for 35 g/L, 50% recovery
        """
        lcow_results = lcow_results_flow
        sec_results = sec_results_flow
        m = multistage.run_n_stage_system(
            n_stages=n,
            salinity=35,
            add_erd=erd,
            flow_vol=flow,
            pump_dict={1: True, 2: True, 3: False},
            ro_op_dict=default_ro_op_dict,
            RO_1D=False,
        )

        m = multistage.set_system_recovery(m, 0.5)

        _ = utils.solve(model=m, tee=False)

        assert pytest.approx(lcow_results[flow][n][erd], rel=1e-3) == value(
            m.fs.costing.LCOW
        )
        assert pytest.approx(sec_results[flow][n][erd], rel=1e-3) == value(
            m.fs.costing.SEC
        )

    @pytest.mark.parametrize("n", n_stages)
    @pytest.mark.parametrize("erd", add_erd)
    @pytest.mark.parametrize("recov", lcow_results_recov.keys())
    @pytest.mark.integration
    def test_multistage_ro_recovery_range(self, recov, n, erd):
        """
        Test up to 3-stage with and without ERD, with booster pump for recovery range
        for 35 g/L, 1 L/s
        """
        lcow_results = lcow_results_recov
        sec_results = sec_results_recov
        m = multistage.run_n_stage_system(
            n_stages=n,
            salinity=35,
            add_erd=erd,
            flow_vol=1,
            pump_dict={1: True, 2: True, 3: False},
            ro_op_dict=default_ro_op_dict,
            RO_1D=False,
        )

        m = multistage.set_system_recovery(m, recov)

        _ = utils.solve(model=m, tee=False)

        assert pytest.approx(lcow_results[recov][n][erd], rel=1e-3) == value(
            m.fs.costing.LCOW
        )
        assert pytest.approx(sec_results[recov][n][erd], rel=1e-3) == value(
            m.fs.costing.SEC
        )
