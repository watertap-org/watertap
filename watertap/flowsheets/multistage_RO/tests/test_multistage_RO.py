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
from watertap.unit_models.reverse_osmosis_1D import ReverseOsmosis1D
from watertap.unit_models.pressure_changer import Pump
import watertap.flowsheets.multistage_RO.multistage_RO as multistage
import watertap.flowsheets.multistage_RO.utils as utils

lcow_results_salinity = {
    5: {
        1: {False: 0.3501, True: 0.3307},
        2: {False: 0.3485, True: 0.3294},
        3: {False: 0.3478, True: 0.329},
    },
    100: {
        1: {False: 2.8075, True: 2.2792},
        2: {False: 2.6304, True: 2.105},
        3: {False: 2.5822, True: 2.0704},
    },
}

sec_results_salinity = {
    5: {
        1: {False: 1.1852, True: 0.7986},
        2: {False: 1.1822, True: 0.7975},
        3: {False: 1.1788, True: 0.7974},
    },
    100: {
        1: {False: 15.647, True: 9.688},
        2: {False: 14.6124, True: 8.6606},
        3: {False: 14.3329, True: 8.5305},
    },
}


lcow_results_flow = {
    1: {
        1: {False: 0.8456, True: 0.7303},
        2: {False: 0.831, True: 0.7145},
        3: {False: 0.826, True: 0.7109},
    },
    1000: {
        1: {False: 0.8496, True: 0.7313},
        2: {False: 0.8333, True: 0.7167},
        3: {False: 0.8284, True: 0.7132},
    },
}
sec_results_flow = {
    1: {
        1: {False: 3.9539, True: 2.5279},
        2: {False: 3.8437, True: 2.3956},
        3: {False: 3.8237, True: 2.3882},
    },
    1000: {
        1: {False: 4.1014, True: 2.5725},
        2: {False: 3.8789, True: 2.4304},
        3: {False: 3.8618, True: 2.4258},
    },
}


lcow_results_recov = {
    0.4: {
        1: {False: 0.9153, True: 0.7719},
        2: {False: 0.9072, True: 0.7654},
        3: {False: 0.904, True: 0.7624},
    },
    0.5: {
        1: {False: 0.8456, True: 0.7303},
        2: {False: 0.831, True: 0.7145},
        3: {False: 0.826, True: 0.7109},
    },
    0.6: {
        1: {False: 0.8388, True: 0.7419},
        2: {False: 0.8003, True: 0.701},
        3: {False: 0.7954, True: 0.6972},
    },
}


sec_results_recov = {
    0.4: {
        1: {False: 4.2266, True: 2.3989},
        2: {False: 4.1896, True: 2.3768},
        3: {False: 4.1773, True: 2.3472},
    },
    0.5: {
        1: {False: 3.9539, True: 2.5279},
        2: {False: 3.8437, True: 2.3956},
        3: {False: 3.8237, True: 2.3882},
    },
    0.6: {
        1: {False: 4.0167, True: 2.8528},
        2: {False: 3.7355, True: 2.5399},
        3: {False: 3.7172, True: 2.5317},
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


class TestMultiStageRO:

    @pytest.mark.unit
    def test_wrong_prop_pack(self):
        with pytest.raises(
            ValueError,
            match="Only NaCl or Seawater property models can be used but MCAS was passed.",
        ):
            _ = multistage.build_n_stage_system(prop_pack="MCAS")

    @pytest.mark.unit
    def test_NaCl_build(self):
        m = multistage.build_n_stage_system(
            n_stages=5, add_costing=False, prop_pack="NACL"
        )
        assert m.fs.properties.solute_set.at(1) == "NaCl"
        assert isinstance(m.fs.properties, NaClParameterBlock)
        assert len(m.fs.stage) == 5
        assert len(m.fs.stages_set) == 5
        assert len(m.fs.product_mixer.config.inlet_list) == len(m.fs.stages_set)

        for n, stage in m.fs.stage.items():
            assert isinstance(stage.RO, ReverseOsmosis1D)
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
        )
        assert m.fs.properties.solute_set.at(1) == "TDS"
        assert isinstance(m.fs.properties, SeawaterParameterBlock)
        assert len(m.fs.stage) == 3
        for n, stage in m.fs.stage.items():
            assert isinstance(stage.RO, ReverseOsmosis1D)
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
        m = multistage.run_n_stage_system()
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
        )

        m = multistage.set_system_recovery(m, recov)

        _ = utils.solve(model=m, tee=False)

        assert pytest.approx(lcow_results[recov][n][erd], rel=1e-3) == value(
            m.fs.costing.LCOW
        )
        assert pytest.approx(sec_results[recov][n][erd], rel=1e-3) == value(
            m.fs.costing.SEC
        )
