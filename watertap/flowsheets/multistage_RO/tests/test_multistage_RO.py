import pytest
from pyomo.environ import value, units as pyunits

from watertap.property_models.NaCl_T_dep_prop_pack import NaClParameterBlock
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.unit_models.reverse_osmosis_1D import ReverseOsmosis1D
from watertap.unit_models.pressure_changer import Pump
import watertap.flowsheets.multistage_RO.multistage_RO as multistage
import watertap.flowsheets.multistage_RO.utils as utils

# lcow_results[salinity][n_stages][add_erd] = LCOW
lcow_results_salinity = {
    5: {
        1: {False: 0.2646417561974558, True: 0.24997536766294023},
        2: {False: 0.263486123447289, True: 0.24904023931106112},
        3: {False: 0.26290582661705675, True: 0.2486799517624164},
    },
    35: {
        1: {False: 0.639173643425559, True: 0.5520735281058546},
        2: {False: 0.6281747023336776, True: 0.5400820338676561},
        3: {False: 0.6244117389129622, True: 0.5373567250248367},
    },
    75: {
        1: {False: 1.3707996721914295, True: 1.1280897471800224},
        2: {False: 1.310031037609782, True: 1.066010790989584},
        3: {False: 1.2957645082194817, True: 1.056715235134523},
    },
}

# sec_results[salinity][n_stages][add_erd] = SEC
sec_results_salinity = {
    5: {
        1: {False: 1.185200975449009, True: 0.7986981264221501},
        2: {False: 1.1822460224194102, True: 0.7975791321931648},
        3: {False: 1.1788834876349958, True: 0.7974824333627738},
    },
    35: {
        1: {False: 3.9539366172570785, True: 2.527912638003297},
        2: {False: 3.843786781180363, True: 2.395669592457694},
        3: {False: 3.823722058124385, True: 2.3882162553839885},
    },
    75: {
        1: {False: 9.857465150091395, True: 6.137098051963456},
        2: {False: 9.363590468062696, True: 5.608104313741416},
        3: {False: 9.250152434952605, True: 5.564988640713179},
    },
}


lcow_results_flow = {
    1: {
        1: {False: 0.639173643425559, True: 0.5520735281058546},
        2: {False: 0.6281747023336776, True: 0.5400820338676561},
        3: {False: 0.6244117389129622, True: 0.5373567250248367},
    },
    5: {
        1: {False: 0.6391736470342906, True: 0.5520735310135787},
        2: {False: 0.6281747057787588, True: 0.5400820365865827},
        3: {False: 0.6244117423389374, True: 0.5373567277301287},
    },
    10: {
        1: {False: 0.6391736470431408, True: 0.5520735310135786},
        2: {False: 0.6281747057787544, True: 0.5400820365865802},
        3: {False: 0.6244117423435839, True: 0.5373567277298809},
    },
}
sec_results_flow = {
    1: {
        1: {False: 3.9539366172570785, True: 2.527912638003297},
        2: {False: 3.843786781180363, True: 2.395669592457694},
        3: {False: 3.823722058124385, True: 2.3882162553839885},
    },
    5: {
        1: {False: 3.95393664516048, True: 2.5279126552463063},
        2: {False: 3.8437868076866253, True: 2.3956696081618216},
        3: {False: 3.823722084869795, True: 2.3882158722428675},
    },
    10: {
        1: {False: 3.9539366463023233, True: 2.5279126557388425},
        2: {False: 3.843786809243369, True: 2.39566960911794},
        3: {False: 3.8237226074695707, True: 2.38821627099504},
    },
}

lcow_results_recov = {
    0.4: {
        1: {False: 0.6919036149545981, True: 0.5835019346358996},
        2: {False: 0.6857187546581422, True: 0.5785921691490907},
        3: {False: 0.6833326235243236, True: 0.5763050336260196},
    },
    0.5: {
        1: {False: 0.639173643425559, True: 0.5520735281058546},
        2: {False: 0.6281747023336776, True: 0.5400820338676561},
        3: {False: 0.6244117389129622, True: 0.5373567250248367},
    },
    0.6: {
        1: {False: 0.634032081033701, True: 0.5608013704587148},
        2: {False: 0.6049856059766691, True: 0.5298662384592617},
        3: {False: 0.6012572538785682, True: 0.5270390053395566},
    },
}


sec_results_recov = {
    0.4: {
        1: {False: 4.226692698047057, True: 2.398919476685519},
        2: {False: 4.18962300828809, True: 2.3768661562930875},
        3: {False: 4.177362916705994, True: 2.347255523779293},
    },
    0.5: {
        1: {False: 3.9539366172570785, True: 2.527912638003297},
        2: {False: 3.843786781180363, True: 2.395669592457694},
        3: {False: 3.823722058124385, True: 2.3882162553839885},
    },
    0.6: {
        1: {False: 4.016746989547713, True: 2.8528281788348218},
        2: {False: 3.7355669156680023, True: 2.539915317385019},
        3: {False: 3.7172824375971394, True: 2.531721996885654},
    },
}

salinity = [5, 35, 75]  # g/L
flows = [1, 5, 10]  # L/s
recovery = [0.4, 0.5, 0.6]
n_stages = [1, 2, 3]
add_erd = [True, False]

salinity = [35]
flows = [1]
recovery = [0.5]
n_stages = [2]
add_erd = [True]

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
    @pytest.mark.parametrize("salt", salinity)
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
    @pytest.mark.parametrize("flow", flows)
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
    @pytest.mark.parametrize("recov", recovery)
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
