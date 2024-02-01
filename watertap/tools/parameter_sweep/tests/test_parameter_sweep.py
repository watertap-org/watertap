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
import os
import numpy as np
import pyomo.environ as pyo

from pyomo.environ import value

from watertap.tools.parameter_sweep.sampling_types import *
from watertap.tools.parameter_sweep import ParameterSweep, parameter_sweep
from watertap.tools.parameter_sweep.parameter_sweep_writer import *

import watertap.tools.MPI as MPI
from watertap.tools.parallel.parallel_manager_factory import has_mpi_peer_processes

# -----------------------------------------------------------------------------


def build_model():
    m = pyo.ConcreteModel()
    m.fs = fs = pyo.Block()

    fs.input = pyo.Var(["a", "b"], within=pyo.UnitInterval, initialize=0.5)
    fs.output = pyo.Var(["c", "d"], within=pyo.UnitInterval, initialize=0.5)

    fs.slack = pyo.Var(["ab_slack", "cd_slack"], bounds=(0, 0), initialize=0.0)
    fs.slack_penalty = pyo.Param(default=1000.0, mutable=True, within=pyo.PositiveReals)

    fs.ab_constr = pyo.Constraint(
        expr=(fs.output["c"] + fs.slack["ab_slack"] == 2 * fs.input["a"])
    )
    fs.cd_constr = pyo.Constraint(
        expr=(fs.output["d"] + fs.slack["cd_slack"] == 3 * fs.input["b"])
    )

    fs.performance = pyo.Expression(expr=pyo.summation(fs.output))

    m.objective = pyo.Objective(
        expr=m.fs.performance - m.fs.slack_penalty * pyo.summation(m.fs.slack),
        sense=pyo.maximize,
    )
    return m


def build_model_for_tps():
    model = build_model()
    model.fs.slack_penalty = 1000.0
    model.fs.slack.setub(0)
    return model


def build_sweep_params_for_tps(model):
    A = model.fs.input["a"]
    B = model.fs.input["b"]
    sweep_params = {A.name: (A, 0.1, 0.9, 3), B.name: (B, 0.0, 0.5, 3)}
    return sweep_params


def build_outputs_for_tps(model, output_keys):
    outputs = {}
    for key, pyo_object in output_keys.items():
        outputs[key] = model.find_component(pyo_object)
    return outputs


def dummy_kernel_logic(solution_succesful):
    init_state = [True]
    solved_state = [False]
    for sf in solution_succesful:
        if sf and init_state[-1]:
            # we solved model from init and/or prior solved state
            init_state.append(True)
            solved_state.append(True)
        elif sf and init_state[-1] == False:
            # we try to solve ,but first init
            init_state.append(True)
            solved_state.append(False)
            # solve is succesful
            init_state.append(True)
            solved_state.append(True)
        else:
            if solved_state[-1]:
                # this means solution failed after
                # on model that was solved with applied sweep params
                # we add failed states
                init_state.append(False)
                solved_state.append(False)
                # kernel reinits model and then tries solving again
                init_state.append(True)
                solved_state.append(False)
                # but it fails again
                init_state.append(False)
                solved_state.append(False)
            else:
                # we are solving from failed solve, so we reinit
                init_state.append(True)
                solved_state.append(False)
                # but fail again
                init_state.append(False)
                solved_state.append(False)
                # kernel solving from inited state, so we failed we move on
    return init_state, solved_state


class TestParameterSweep:
    @pytest.fixture(scope="class")
    def model(self):
        return build_model()

    @pytest.mark.unit
    def test_single_index_unrolled(self):
        indexed_var = pyo.Var(["a"])
        indexed_var.construct()

        ls = LinearSample(indexed_var, None, None, None)

        assert ls.pyomo_object is indexed_var["a"]

    @pytest.mark.unit
    def test_multiple_indices_error(self):
        indexed_var = pyo.Var(["a", "b"])
        indexed_var.construct()

        with pytest.raises(Exception):
            ls = LinearSample(indexed_var, None, None, None)

    @pytest.mark.component
    def test_linear_build_combinations(self):
        ps = ParameterSweep()

        A_param = pyo.Param(initialize=0.0, mutable=True)
        B_param = pyo.Param(initialize=1.0, mutable=True)
        C_param = pyo.Param(initialize=2.0, mutable=True)

        range_A = [0.0, 10.0]
        range_B = [1.0, 20.0]
        range_C = [2.0, 30.0]

        nn_A = 4
        nn_B = 5
        nn_C = 6

        param_dict = dict()
        param_dict["var_A"] = LinearSample(A_param, range_A[0], range_A[1], nn_A)
        param_dict["var_B"] = LinearSample(B_param, range_B[0], range_B[1], nn_B)
        param_dict["var_C"] = LinearSample(C_param, range_C[0], range_C[1], nn_C)

        global_combo_array = ps._build_combinations(
            param_dict, SamplingType.FIXED, None
        )

        assert np.shape(global_combo_array)[0] == nn_A * nn_B * nn_C
        assert np.shape(global_combo_array)[1] == len(param_dict)

        assert global_combo_array[0, 0] == pytest.approx(range_A[0])
        assert global_combo_array[0, 1] == pytest.approx(range_B[0])
        assert global_combo_array[0, 2] == pytest.approx(range_C[0])

        assert global_combo_array[-1, 0] == pytest.approx(range_A[1])
        assert global_combo_array[-1, 1] == pytest.approx(range_B[1])
        assert global_combo_array[-1, 2] == pytest.approx(range_C[1])

    @pytest.mark.component
    def test_mode_type_build_combinations(self):
        ps = ParameterSweep()
        model = build_model_for_tps()
        model.A_param = pyo.Var(initialize=0.0, bounds=(None, None))
        model.D_param = pyo.Var(initialize=0.0, bounds=(None, None))
        model.B_param = pyo.Param(initialize=1.0, mutable=True)
        model.C_param = pyo.Var(list(range(50)), initialize=1.0)
        range_A = [0.0, 10.0]

        range_A_fix = [0, 1]

        nn_A = 4
        sample_values = np.linspace(range_A[0], range_A[1], nn_A)
        # test ub_mode
        param_dict = dict()
        param_dict["var_A"] = LinearSample(model.A_param, range_A[0], range_A[1], nn_A)
        param_dict["var_A"].set_variable_update_mode(SetMode.SET_UB)

        global_combo_array = ps._build_combinations(
            param_dict, SamplingType.FIXED, None
        )
        for k in range(len(sample_values)):
            ps._update_model_values(model, param_dict, global_combo_array[k])
            assert model.A_param.ub == pytest.approx(sample_values[k])
        param_dict = dict()
        param_dict["var_A"] = LinearSample(model.A_param, range_A[0], range_A[1], nn_A)
        param_dict["var_A"].set_variable_update_mode(SetMode.SET_LB)

        global_combo_array = ps._build_combinations(
            param_dict, SamplingType.FIXED, None
        )
        for k in range(len(sample_values)):
            ps._update_model_values(model, param_dict, global_combo_array[k])
            assert model.A_param.lb == pytest.approx(sample_values[k])

        param_dict = dict()
        param_dict["var_A"] = PredeterminedFixedSample(model.A_param, [False, True])
        param_dict["var_A"].set_variable_update_mode(SetMode.SET_FIXED_STATE, 10)

        global_combo_array = ps._build_combinations(
            param_dict, SamplingType.FIXED, None
        )
        for k in range(2):
            ps._update_model_values(model, param_dict, global_combo_array[k])
            if sample_values[k] < 0.5:
                fixed_state = False
            else:
                fixed_state = True
            assert model.A_param.fixed == fixed_state
            assert model.A_param.value == pytest.approx(10)

        # test ablity to fully speciy a var an unfix it
        param_dict = dict()
        param_dict["var_A_ub"] = PredeterminedFixedSample(model.A_param, [range_A[1]])
        param_dict["var_A_ub"].set_variable_update_mode(SetMode.SET_UB)
        param_dict["var_A_lb"] = PredeterminedFixedSample(model.A_param, [range_A[0]])
        param_dict["var_A_lb"].set_variable_update_mode(SetMode.SET_LB)
        param_dict["var_A_fixed_state"] = PredeterminedFixedSample(
            model.A_param, [False]
        )
        param_dict["var_A_fixed_state"].set_variable_update_mode(
            SetMode.SET_FIXED_STATE, 5
        )
        global_combo_array = ps._build_combinations(
            param_dict, SamplingType.FIXED, None
        )

        ps._update_model_values(model, param_dict, global_combo_array[0])
        assert model.A_param.lb == pytest.approx(range_A[0])
        assert model.A_param.ub == pytest.approx(range_A[1])
        assert model.A_param.value == pytest.approx(5)
        assert model.A_param.fixed == False

        # test ablity to fully speciy a var and fix it
        param_dict = dict()
        param_dict["var_A_ub"] = PredeterminedFixedSample(model.A_param, [range_A[1]])
        param_dict["var_A_ub"].set_variable_update_mode(SetMode.SET_UB)
        param_dict["var_A_lb"] = PredeterminedFixedSample(model.A_param, [range_A[0]])
        param_dict["var_A_lb"].set_variable_update_mode(SetMode.SET_LB)
        param_dict["var_A_fixed_state"] = PredeterminedFixedSample(
            model.A_param, [True]
        )
        param_dict["var_A_fixed_state"].set_variable_update_mode(
            SetMode.SET_FIXED_STATE, 5
        )
        global_combo_array = ps._build_combinations(
            param_dict, SamplingType.FIXED, None
        )

        ps._update_model_values(model, param_dict, global_combo_array[0])
        assert model.A_param.lb == pytest.approx(range_A[0])
        assert model.A_param.ub == pytest.approx(range_A[1])
        assert model.A_param.value == pytest.approx(5)
        assert model.A_param.fixed == True

        # test generating large number of changed vars and single sweep var
        param_dict = dict()
        param_dict["var_A"] = LinearSample(model.A_param, range_A[0], range_A[1], nn_A)
        param_dict["var_A"].set_variable_update_mode(SetMode.FIX_VALUE)
        param_dict["var_A_lb"] = PredeterminedFixedSample(model.A_param, [range_A[0]])
        param_dict["var_A_lb"].set_variable_update_mode(SetMode.SET_LB)

        for k in range(50):
            param_dict["var_b_{}".format(k)] = PredeterminedFixedSample(
                model.C_param[k], [k]
            )
        global_combo_array = ps._build_combinations(
            param_dict, SamplingType.FIXED, None
        )
        for i in range(nn_A):
            ps._update_model_values(model, param_dict, global_combo_array[i])
            assert model.A_param.value == pytest.approx(sample_values[i])
            assert model.A_param.ub == pytest.approx(range_A[1])
            assert model.A_param.fixed == True
            for k in range(50):
                assert model.C_param[k].value == pytest.approx(k)
        # test generating large number of changed vars and multi sweep var
        sample_values_d = np.linspace(range_A[0] * 2, range_A[1] * 2, nn_A)
        temp_global_combo_array = np.array(
            np.meshgrid(*[sample_values, sample_values_d], indexing="ij")
        )
        temp_global_combo_array = temp_global_combo_array.reshape(2, -1).T
        param_dict = dict()
        param_dict["var_A"] = LinearSample(model.A_param, range_A[0], range_A[1], nn_A)
        param_dict["var_A"].set_variable_update_mode(SetMode.FIX_VALUE)
        param_dict["var_D"] = LinearSample(
            model.D_param, range_A[0] * 2, range_A[1] * 2, nn_A
        )
        param_dict["var_D"].set_variable_update_mode(SetMode.FIX_VALUE)
        param_dict["var_A_lb"] = PredeterminedFixedSample(model.A_param, [range_A[0]])
        param_dict["var_A_lb"].set_variable_update_mode(SetMode.SET_LB)

        for k in range(50):
            param_dict["var_b_{}".format(k)] = PredeterminedFixedSample(
                model.C_param[k], [k]
            )
        global_combo_array = ps._build_combinations(
            param_dict, SamplingType.FIXED, None
        )
        for i in range(temp_global_combo_array.shape[1]):
            ps._update_model_values(model, param_dict, global_combo_array[i])
            assert model.A_param.value == pytest.approx(temp_global_combo_array[i][0])
            assert model.D_param.value == pytest.approx(temp_global_combo_array[i][1])
            assert model.A_param.ub == pytest.approx(range_A[1])
            assert model.A_param.fixed == True
            for k in range(50):
                assert model.C_param[k].value == pytest.approx(k)

    @pytest.mark.component
    def test_geom_build_combinations(self):
        ps = ParameterSweep()

        A_param = pyo.Param(initialize=0.0, mutable=True)
        B_param = pyo.Param(initialize=1.0, mutable=True)
        C_param = pyo.Param(initialize=2.0, mutable=True)

        range_A = [1.0, 10.0]
        range_B = [2.0, 20.0]
        range_C = [3.0, 30.0]

        nn_A = 2
        nn_B = 3
        nn_C = 4

        param_dict = dict()
        param_dict["var_A"] = GeomSample(A_param, range_A[0], range_A[1], nn_A)
        param_dict["var_B"] = GeomSample(B_param, range_B[0], range_B[1], nn_B)
        param_dict["var_C"] = GeomSample(C_param, range_C[0], range_C[1], nn_C)

        global_combo_array = ps._build_combinations(
            param_dict, SamplingType.FIXED, None
        )

        assert np.shape(global_combo_array)[0] == nn_A * nn_B * nn_C
        assert np.shape(global_combo_array)[1] == len(param_dict)

        assert global_combo_array[0, 0] == pytest.approx(range_A[0])
        assert global_combo_array[0, 1] == pytest.approx(range_B[0])
        assert global_combo_array[0, 2] == pytest.approx(range_C[0])

        assert global_combo_array[-1, 0] == pytest.approx(range_A[1])
        assert global_combo_array[-1, 1] == pytest.approx(range_B[1])
        assert global_combo_array[-1, 2] == pytest.approx(range_C[1])

    @pytest.mark.component
    def test_reverse_geom_build_combinations(self):
        ps = ParameterSweep()

        A_param = pyo.Param(initialize=0.0, mutable=True)
        B_param = pyo.Param(initialize=1.0, mutable=True)
        C_param = pyo.Param(initialize=2.0, mutable=True)

        range_A = [1.0, 10.0]
        range_B = [2.0, 20.0]
        range_C = [3.0, 30.0]

        nn_A = 2
        nn_B = 3
        nn_C = 4

        param_dict = dict()
        param_dict["var_A"] = ReverseGeomSample(A_param, range_A[0], range_A[1], nn_A)
        param_dict["var_B"] = ReverseGeomSample(B_param, range_B[0], range_B[1], nn_B)
        param_dict["var_C"] = ReverseGeomSample(C_param, range_C[0], range_C[1], nn_C)

        global_combo_array = ps._build_combinations(
            param_dict, SamplingType.FIXED, None
        )

        assert np.shape(global_combo_array)[0] == nn_A * nn_B * nn_C
        assert np.shape(global_combo_array)[1] == len(param_dict)

        assert global_combo_array[0, 0] == pytest.approx(range_A[0])
        assert global_combo_array[0, 1] == pytest.approx(range_B[0])
        assert global_combo_array[0, 2] == pytest.approx(range_C[0])

        assert global_combo_array[-1, 0] == pytest.approx(range_A[1])
        assert global_combo_array[-1, 1] == pytest.approx(range_B[1])
        assert global_combo_array[-1, 2] == pytest.approx(range_C[1])

    @pytest.mark.component
    def test_predetermined_fixed_build_combinations(self):
        ps = ParameterSweep()

        A_param = pyo.Param(initialize=0.0, mutable=True)
        B_param = pyo.Param(initialize=1.0, mutable=True)
        C_param = pyo.Param(initialize=2.0, mutable=True)

        nn_A = 2
        nn_B = 3
        nn_C = 4

        A_values = nn_A + np.arange(nn_A)
        B_values = nn_B + np.arange(nn_B)
        C_values = nn_C + np.arange(nn_C)

        param_dict = dict()
        param_dict["var_A"] = PredeterminedFixedSample(A_param, A_values)
        param_dict["var_B"] = PredeterminedFixedSample(B_param, B_values)
        param_dict["var_C"] = PredeterminedFixedSample(C_param, C_values)

        global_combo_array = ps._build_combinations(
            param_dict, SamplingType.FIXED, None
        )

        assert np.shape(global_combo_array)[0] == nn_A * nn_B * nn_C
        assert np.shape(global_combo_array)[1] == len(param_dict)
        assert (global_combo_array[0] == np.array([nn_A, nn_B, nn_C])).all()
        assert (global_combo_array[-1] == 2 * np.array([nn_A, nn_B, nn_C]) - 1).all()

    @pytest.mark.component
    def test_status_publishing(self):
        requests = pytest.importorskip(
            "requests",
            reason="requests (parameter_sweep optional dependency) not available",
        )
        ps = ParameterSweep(
            publish_progress=True, publish_address="http://localhost:8888"
        )

        with pytest.raises(requests.exceptions.ConnectionError):
            r = ps._publish_updates(1, True, 5.0)

    def test_random_build_combinations(self):
        ps = ParameterSweep()

        nn = int(1e5)

        # Uniform random sampling [lower_limit, upper_limit]
        A_param = pyo.Param(initialize=-10.0, mutable=True)
        B_param = pyo.Param(initialize=0.0, mutable=True)
        C_param = pyo.Param(initialize=10.0, mutable=True)

        range_A = [-10.0, 0.0]
        range_B = [0.0, 10.0]
        range_C = [10.0, 20.0]

        param_dict = dict()
        param_dict["var_A"] = UniformSample(A_param, range_A[0], range_A[1], nn)
        param_dict["var_B"] = UniformSample(B_param, range_B[0], range_B[1], nn)
        param_dict["var_C"] = UniformSample(C_param, range_C[0], range_C[1], nn)

        global_combo_array = ps._build_combinations(param_dict, SamplingType.RANDOM, nn)

        assert np.shape(global_combo_array)[0] == nn
        assert np.shape(global_combo_array)[1] == len(param_dict)

        assert np.all(range_A[0] < global_combo_array[:, 0])
        assert np.all(range_B[0] < global_combo_array[:, 1])
        assert np.all(range_C[0] < global_combo_array[:, 2])

        assert np.all(global_combo_array[:, 0] < range_A[1])
        assert np.all(global_combo_array[:, 1] < range_B[1])
        assert np.all(global_combo_array[:, 2] < range_C[1])

        # Normal random sampling [mean, stdev]
        A_param = pyo.Param(initialize=10.0, mutable=True)
        B_param = pyo.Param(initialize=100.0, mutable=True)
        C_param = pyo.Param(initialize=1000.0, mutable=True)

        range_A = [10.0, 5.0]
        range_B = [100.0, 50.0]
        range_C = [1000.0, 0.0]

        param_dict = dict()
        param_dict["var_A"] = NormalSample(A_param, range_A[0], range_A[1], nn)
        param_dict["var_B"] = NormalSample(B_param, range_B[0], range_B[1], nn)
        param_dict["var_C"] = NormalSample(C_param, range_C[0], range_C[1], nn)

        global_combo_array = ps._build_combinations(param_dict, SamplingType.RANDOM, nn)

        assert np.shape(global_combo_array)[0] == nn
        assert np.shape(global_combo_array)[1] == len(param_dict)

        assert np.mean(global_combo_array[:, 0]) < (range_A[0] + range_A[1])
        assert np.mean(global_combo_array[:, 1]) < (range_B[0] + range_B[1])

        assert (range_A[0] - range_A[1]) < np.mean(global_combo_array[:, 0])
        assert (range_B[0] - range_B[1]) < np.mean(global_combo_array[:, 1])

        assert np.all(global_combo_array[:, 2] == range_C[0])

        A_param = pyo.Param(initialize=0.0, mutable=True)
        B_param = pyo.Param(initialize=0.0, mutable=True)
        C_param = pyo.Param(initialize=0.0, mutable=True)

        A_values = np.arange(nn)
        B_values = 10.0 + A_values
        C_values = 20.0 + A_values

        param_dict = dict()
        param_dict["var_A"] = PredeterminedRandomSample(A_param, A_values)
        param_dict["var_B"] = PredeterminedRandomSample(B_param, B_values)
        param_dict["var_C"] = PredeterminedRandomSample(C_param, C_values)

        global_combo_array = ps._build_combinations(param_dict, SamplingType.RANDOM, nn)

        assert np.shape(global_combo_array)[0] == nn
        assert (global_combo_array == np.array([A_values, B_values, C_values]).T).all()

    @pytest.mark.component
    def test_divide_combinations(self):
        ps = ParameterSweep()

        A_param = pyo.Param(initialize=0.0, mutable=True)
        B_param = pyo.Param(initialize=1.0, mutable=True)
        C_param = pyo.Param(initialize=2.0, mutable=True)

        range_A = [0.0, 10.0]
        range_B = [1.0, 20.0]
        range_C = [2.0, 30.0]

        nn_A = 4
        nn_B = 5
        nn_C = 6

        param_dict = dict()
        param_dict["var_A"] = LinearSample(A_param, range_A[0], range_A[1], nn_A)
        param_dict["var_B"] = LinearSample(B_param, range_B[0], range_B[1], nn_B)
        param_dict["var_C"] = LinearSample(C_param, range_C[0], range_C[1], nn_C)

        global_combo_array = ps._build_combinations(
            param_dict, SamplingType.FIXED, None
        )

        num_procs = ps.parallel_manager.number_of_worker_processes()
        rank = ps.parallel_manager.get_rank()
        test = np.array_split(global_combo_array, num_procs, axis=0)[rank]

        local_combo_array = ps._divide_combinations(global_combo_array)

        assert np.shape(local_combo_array)[1] == 3

        assert np.allclose(test[:, 0], local_combo_array[:, 0])
        assert np.allclose(test[:, 1], local_combo_array[:, 1])
        assert np.allclose(test[:, 2], local_combo_array[:, 2])

        if rank == 0:
            assert local_combo_array[0, 0] == pytest.approx(range_A[0])
            assert local_combo_array[0, 1] == pytest.approx(range_B[0])
            assert local_combo_array[0, 2] == pytest.approx(range_C[0])

        if rank == num_procs - 1:
            assert local_combo_array[-1, 0] == pytest.approx(range_A[1])
            assert local_combo_array[-1, 1] == pytest.approx(range_B[1])
            assert local_combo_array[-1, 2] == pytest.approx(range_C[1])

    @pytest.mark.component
    def test_update_model_values(self, model):
        m = model
        ps = ParameterSweep()

        param_dict = dict()
        param_dict["input_a"] = LinearSample(m.fs.input["a"], None, None, None)
        param_dict["input_b"] = LinearSample(m.fs.input["b"], None, None, None)

        original_a = value(m.fs.input["a"])
        original_b = value(m.fs.input["b"])

        new_values = [1.1 * original_a, 1.1 * original_b]

        ps._update_model_values(m, param_dict, new_values)

        assert value(m.fs.input["a"]) == pytest.approx(new_values[0])
        assert value(m.fs.input["b"]) == pytest.approx(new_values[1])

    @pytest.mark.unit
    def test_aggregate_results_arr(self):
        ps = ParameterSweep()

        input_dict = {
            "outputs": {
                "fs.output[c]": {
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 1,
                    "value": np.array([0.2, 0.2, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0]),
                },
                "fs.output[d]": {
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 1,
                    "value": np.array([0.0, 0.75, 0.0, 0.0, 0.75, 0.0, 0.0, 0.0, 0.0]),
                },
                "fs.slack[ab_slack]": {
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 1,
                    "value": np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                },
                "fs.slack[cd_slack]": {
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 1,
                    "value": np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                },
            },
            "solve_successful": [True] * 9,
            "sweep_params": {
                "fs.input[a]": {
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 1,
                    "value": np.array([0.1, 0.1, 0.0, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0]),
                },
                "fs.input[b]": {
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 1,
                    "value": np.array([0.0, 0.25, 0.0, 0.0, 0.25, 0.0, 0.0, 0.0, 0.0]),
                },
            },
        }

        global_num_cases = len(input_dict["sweep_params"]["fs.input[a]"]["value"])
        global_results_arr = ps._aggregate_results_arr(input_dict, global_num_cases)

        reference_results_arr = np.array(
            [
                [0.2, 0.0, 0.0, 0.0],
                [0.2, 0.75, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0, 0.0],
                [1.0, 0.75, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
            ]
        )

        assert np.array_equal(global_results_arr, reference_results_arr)

    @pytest.mark.unit
    def test_create_local_output_skeleton(self, model):
        ps = ParameterSweep()

        m = model
        m.fs.slack_penalty = 1000.0
        m.fs.slack.setub(0)

        sweep_params = {
            "input_a": (m.fs.input["a"], 0.1, 0.9, 3),
            "input_b": (m.fs.input["b"], 0.0, 0.5, 3),
        }
        # outputs = {
        #     "output_c": m.fs.output["c"],
        #     "output_d": m.fs.output["d"],
        #     "performance": m.fs.performance,
        # }

        sweep_params, sampling_type = ps._process_sweep_params(sweep_params)
        values = ps._build_combinations(sweep_params, sampling_type, None)
        num_cases = np.shape(values)[0]
        output_dict = ps._create_local_output_skeleton(
            model, sweep_params, None, num_cases
        )

        truth_dict = {
            "outputs": {
                "input_a": {
                    "full_name": "fs.input[a]",
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 1,
                    "value": np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                },
                "input_b": {
                    "full_name": "fs.input[b]",
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 1,
                    "value": np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                },
                "fs.output[c]": {
                    "full_name": "fs.output[c]",
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 1,
                    "value": np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                },
                "fs.output[d]": {
                    "full_name": "fs.output[d]",
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 1,
                    "value": np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                },
                "fs.performance": {
                    "full_name": "fs.performance",
                    "value": np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                },
                "fs.slack[ab_slack]": {
                    "full_name": "fs.slack[ab_slack]",
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 0,
                    "value": np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                },
                "fs.slack[cd_slack]": {
                    "full_name": "fs.slack[cd_slack]",
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 0,
                    "value": np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                },
                "fs.slack_penalty": {
                    "full_name": "fs.slack_penalty",
                    "units": "None",
                    "value": np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                },
                "objective": {
                    "full_name": "objective",
                    "value": np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                },
            },
            "sweep_params": {
                "input_a": {
                    "full_name": "fs.input[a]",
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 1,
                    "value": np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                },
                "input_b": {
                    "full_name": "fs.input[b]",
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 1,
                    "value": np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                },
            },
        }

        _assert_dictionary_correctness(truth_dict, output_dict)

    @pytest.mark.unit
    def test_create_global_output(self, model):
        ps = ParameterSweep()

        m = model
        m.fs.slack_penalty = 1000.0
        m.fs.slack.setub(0)

        global_num_cases = 2 * ps.parallel_manager.number_of_worker_processes()
        sweep_params = {
            "input_a": NormalSample(m.fs.input["a"], 0.1, 0.9, global_num_cases),
            "input_b": NormalSample(m.fs.input["b"], 0.0, 0.5, global_num_cases),
        }
        outputs = {
            "output_c": m.fs.output["c"],
            "output_d": m.fs.output["d"],
            "performance": m.fs.performance,
        }

        sweep_params, sampling_type = ps._process_sweep_params(sweep_params)

        # Get the globale sweep param values
        global_values = ps._build_combinations(
            sweep_params, sampling_type, global_num_cases
        )
        # divide the workload between processors
        local_values = ps._divide_combinations(global_values)
        local_num_cases = np.shape(local_values)[0]

        local_output_dict = ps._create_local_output_skeleton(
            model, sweep_params, None, local_num_cases
        )

        # Manually update the values in the numpy array
        for key, value in local_output_dict.items():
            for subkey, subvalue in value.items():
                subvalue["value"][:] = ps.parallel_manager.get_rank()

        # Local output dict also contains the solve_successful. The solve status is
        # based on the
        local_output_dict["solve_successful"] = [True] * local_num_cases

        # Get the global output dictionary, This is properly created only on rank 0
        global_output_dict = ps._create_global_output(
            local_output_dict, global_num_cases
        )

        if ps.parallel_manager.get_rank() > 0:
            assert global_output_dict == local_output_dict
        else:
            test_array = np.repeat(
                np.arange(
                    0, ps.parallel_manager.number_of_worker_processes(), dtype=float
                ),
                2,
            )
            test_list = [True] * global_num_cases
            for key, value in global_output_dict.items():
                if key != "solve_successful":
                    for subkey, subvalue in value.items():
                        assert np.allclose(subvalue["value"], test_array)
                elif key == "solve_successful":
                    assert list(value) == test_list

    @pytest.mark.component
    @pytest.mark.parametrize("number_of_subprocesses", [1, 2, 3])
    def test_parameter_sweep(self, model, tmp_path, number_of_subprocesses):
        comm = MPI.COMM_WORLD

        tmp_path = _get_rank0_path(comm, tmp_path)

        results_fname = os.path.join(tmp_path, "global_results")
        csv_results_file_name = str(results_fname) + ".csv"
        h5_results_file_name = str(results_fname) + ".h5"

        ps = ParameterSweep(
            optimize_function=_optimization,
            csv_results_file_name=csv_results_file_name,
            h5_results_file_name=h5_results_file_name,
            debugging_data_dir=tmp_path,
            interpolate_nan_outputs=True,
            number_of_subprocesses=number_of_subprocesses,
        )

        # Call the parameter_sweep function
        _ = ps.parameter_sweep(
            build_model_for_tps,
            build_sweep_params_for_tps,
            build_outputs=build_outputs_for_tps,
            build_outputs_kwargs={
                "output_keys": {
                    "output_c": "fs.output[c]",
                    "output_d": "fs.output[d]",
                    "performance": "fs.performance",
                }
            },
        )

        # NOTE: rank 0 "owns" tmp_path, so it needs to be
        #       responsible for doing any output file checking
        #       tmp_path can be deleted as soon as this method
        #       returns
        if ps.parallel_manager.is_root_process():
            # Check that the global results file is created
            assert os.path.isfile(csv_results_file_name)
            assert os.path.isfile(
                os.path.join(tmp_path, "interpolated_global_results.csv")
            )

            # Check that all local output files have been created
            for k in range(ps.parallel_manager.number_of_worker_processes()):
                assert os.path.isfile(
                    os.path.join(tmp_path, f"local_results_{k:03}.h5")
                )
                assert os.path.isfile(
                    os.path.join(tmp_path, f"local_results_{k:03}.csv")
                )

            # Attempt to read in the data
            data = np.genfromtxt(csv_results_file_name, skip_header=1, delimiter=",")

            # Compare the last row of the imported data to truth
            truth_data = [0.9, 0.5, np.nan, np.nan, np.nan]
            assert np.allclose(data[-1], truth_data, equal_nan=True)

        # Check for the h5 output
        if ps.parallel_manager.is_root_process():
            truth_dict = {
                "outputs": {
                    "output_c": {
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array(
                            [0.2, 0.2, np.nan, 1.0, 1.0, np.nan, np.nan, np.nan, np.nan]
                        ),
                    },
                    "output_d": {
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array(
                            [
                                0.0,
                                0.75,
                                np.nan,
                                0.0,
                                0.75,
                                np.nan,
                                np.nan,
                                np.nan,
                                np.nan,
                            ]
                        ),
                    },
                    "performance": {
                        "value": np.array(
                            [
                                0.2,
                                0.95,
                                np.nan,
                                1.0,
                                1.75,
                                np.nan,
                                np.nan,
                                np.nan,
                                np.nan,
                            ]
                        )
                    },
                },
                "solve_successful": [
                    True,
                    True,
                    False,
                    True,
                    True,
                    False,
                    False,
                    False,
                    False,
                ],
                "sweep_params": {
                    "fs.input[a]": {
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array(
                            [0.1, 0.1, 0.1, 0.5, 0.5, 0.5, 0.9, 0.9, 0.9]
                        ),
                    },
                    "fs.input[b]": {
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array(
                            [0.0, 0.25, 0.5, 0.0, 0.25, 0.5, 0.0, 0.25, 0.5]
                        ),
                    },
                },
            }

            read_dict = _read_output_h5(h5_results_file_name)
            _assert_dictionary_correctness(truth_dict, read_dict)
            _assert_h5_csv_agreement(csv_results_file_name, read_dict)

            # Check if there is a text file created
            import ast

            truth_txt_dict = {
                "outputs": [
                    "output_c",
                    "output_d",
                    "performance",
                ],
                "sweep_params": ["fs.input[a]", "fs.input[b]"],
            }

            txt_fpath = h5_results_file_name + ".txt"
            assert os.path.exists(txt_fpath)
            f = open(txt_fpath)
            f_contents = f.read()
            read_txt_dict = ast.literal_eval(f_contents)
            assert read_txt_dict == truth_txt_dict

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_parameter_sweep_optimize(self, model, tmp_path):
        comm = MPI.COMM_WORLD

        tmp_path = _get_rank0_path(comm, tmp_path)
        results_fname = os.path.join(tmp_path, "global_results")
        csv_results_file_name = str(results_fname) + ".csv"
        h5_results_file_name = str(results_fname) + ".h5"

        ps = ParameterSweep(
            optimize_function=_optimization,
            optimize_kwargs={"relax_feasibility": True},
            probe_function=_good_test_function,
            csv_results_file_name=csv_results_file_name,
            h5_results_file_name=h5_results_file_name,
            debugging_data_dir=tmp_path,
            interpolate_nan_outputs=True,
        )

        m = model
        m.fs.slack_penalty = 1000.0
        m.fs.slack.setub(0)

        A = m.fs.input["a"]
        B = m.fs.input["b"]
        sweep_params = {A.name: (A, 0.1, 0.9, 3), B.name: (B, 0.0, 0.5, 3)}
        outputs = {
            "output_c": m.fs.output["c"],
            "output_d": m.fs.output["d"],
            "performance": m.fs.performance,
            "objective": m.objective,
        }
        results_fname = os.path.join(tmp_path, "global_results")
        csv_results_file_name = str(results_fname) + ".csv"
        h5_results_file_name = str(results_fname) + ".h5"

        # Call the parameter_sweep function
        ps.parameter_sweep(
            m,
            sweep_params,
            build_outputs=outputs,
        )

        # NOTE: rank 0 "owns" tmp_path, so it needs to be
        #       responsible for doing any output file checking
        #       tmp_path can be deleted as soon as this method
        #       returns
        if ps.parallel_manager.is_root_process():
            # Check that the global results file is created
            assert os.path.isfile(csv_results_file_name)

            # Attempt to read in the data
            data = np.genfromtxt(csv_results_file_name, skip_header=1, delimiter=",")
            # Compare the last row of the imported data to truth
            truth_data = [
                0.9,
                0.5,
                1.0,
                1.0,
                2.0,
                2.0 - 1000.0 * ((2.0 * 0.9 - 1.0) + (3.0 * 0.5 - 1.0)),
            ]
            assert np.allclose(data[-1], truth_data, equal_nan=True)

        # Check the h5
        if ps.parallel_manager.is_root_process():
            truth_dict = {
                "outputs": {
                    "output_c": {
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array(
                            [0.2, 0.2, 0.2, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
                        ),
                    },
                    "output_d": {
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array(
                            [
                                9.98580690e-09,
                                0.75,
                                1.0,
                                9.99872731e-09,
                                0.75,
                                1.0,
                                9.99860382e-09,
                                0.75,
                                1.0,
                            ]
                        ),
                    },
                    "performance": {
                        "value": np.array(
                            [0.2, 0.95, 1.2, 1.0, 1.75, 2.0, 1.0, 1.75, 2.0]
                        )
                    },
                    "objective": {
                        "value": np.array(
                            [
                                0.2,
                                9.50000020e-01,
                                -4.98799990e02,
                                1.0,
                                1.75,
                                -4.97999990e02,
                                -7.98999990e02,
                                -7.98249990e02,
                                2.0 - 1000.0 * ((2.0 * 0.9 - 1.0) + (3.0 * 0.5 - 1.0)),
                            ]
                        )
                    },
                },
                "solve_successful": [True] * 9,
                "sweep_params": {
                    "fs.input[a]": {
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array(
                            [0.1, 0.1, 0.1, 0.5, 0.5, 0.5, 0.9, 0.9, 0.9]
                        ),
                    },
                    "fs.input[b]": {
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array(
                            [0.0, 0.25, 0.5, 0.0, 0.25, 0.5, 0.0, 0.25, 0.5]
                        ),
                    },
                },
            }

            read_dict = _read_output_h5(h5_results_file_name)
            _assert_dictionary_correctness(truth_dict, read_dict)
            _assert_h5_csv_agreement(csv_results_file_name, read_dict)

    @pytest.mark.component
    def test_parameter_sweep_optimize_with_added_var(self, model, tmp_path):
        # this will run a solve function that adds a variable but only in some
        # of the solves.
        """THIS TEST IS DESIGNED FOR 2 Parallel workers!!!!!!"""
        comm = MPI.COMM_WORLD

        tmp_path = _get_rank0_path(comm, tmp_path)
        results_fname = os.path.join(tmp_path, "global_results")
        csv_results_file_name = str(results_fname) + ".csv"
        h5_results_file_name = str(results_fname) + ".h5"

        ps = ParameterSweep(
            optimize_function=_optimization,
            initialize_function=_initialize_with_added_var,
            update_sweep_params_before_init=True,
            initialize_before_sweep=True,
            optimize_kwargs={"relax_feasibility": True},
            probe_function=_good_test_function,
            csv_results_file_name=csv_results_file_name,
            h5_results_file_name=h5_results_file_name,
            debugging_data_dir=tmp_path,
            interpolate_nan_outputs=False,
            number_of_subprocesses=2,
            parallel_back_end="MultiProcessing",
        )

        results_fname = os.path.join(tmp_path, "global_results")
        csv_results_file_name = str(results_fname) + ".csv"
        h5_results_file_name = str(results_fname) + ".h5"

        # Call the parameter_sweep function
        ps.parameter_sweep(
            build_model_for_tps,
            build_sweep_params_for_tps,
        )

        # NOTE: rank 0 "owns" tmp_path, so it needs to be
        #       responsible for doing any output file checking
        #       tmp_path can be deleted as soon as this method
        #       returns
        # Check that test var array was created
        if ps.parallel_manager.is_root_process():
            truth_dict = {
                "outputs": {
                    "fs.test_var": {
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 10,
                        "value": np.array(
                            [
                                np.nan,
                                np.nan,
                                np.nan,
                                np.nan,
                                np.nan,
                                5,
                                np.nan,
                                np.nan,
                                np.nan,
                            ]
                        ),
                    },
                    "objective": {
                        "value": np.array(
                            [
                                0.2,
                                9.50000020e-01,
                                -4.98799990e02,
                                1.0,
                                1.75,
                                -4.97999990e02,
                                -7.98999990e02,
                                -7.98249990e02,
                                2.0 - 1000.0 * ((2.0 * 0.9 - 1.0) + (3.0 * 0.5 - 1.0)),
                            ]
                        )
                    },
                },
                "solve_successful": [True] * 9,
                "sweep_params": {
                    "fs.input[a]": {
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array(
                            [0.1, 0.1, 0.1, 0.5, 0.5, 0.5, 0.9, 0.9, 0.9]
                        ),
                    },
                    "fs.input[b]": {
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array(
                            [0.0, 0.25, 0.5, 0.0, 0.25, 0.5, 0.0, 0.25, 0.5]
                        ),
                    },
                },
            }

            read_dict = _read_output_h5(h5_results_file_name)
            _assert_dictionary_correctness(truth_dict, read_dict, rtol=1e-2)

    @pytest.mark.component
    def test_parameter_sweep_bad_initialize_call_2(self, model, tmp_path):
        comm = MPI.COMM_WORLD

        tmp_path = _get_rank0_path(comm, tmp_path)
        results_fname = os.path.join(tmp_path, "global_results")
        csv_results_file_name = str(results_fname) + ".csv"
        h5_results_file_name = str(results_fname) + ".h5"

        ps = ParameterSweep(
            optimize_function=_optimization,
            initialize_function=_reinitialize,
            initialize_kwargs={"slack_penalty": 10.0, "foo": "bar"},
            csv_results_file_name=csv_results_file_name,
            h5_results_file_name=h5_results_file_name,
            debugging_data_dir=tmp_path,
            interpolate_nan_outputs=True,
        )

        m = model
        m.fs.slack_penalty = 1000.0
        m.fs.slack.setub(0)

        A = m.fs.input["a"]
        B = m.fs.input["b"]
        sweep_params = {A.name: (A, 0.1, 0.9, 3), B.name: (B, 0.0, 0.5, 3)}

        with pytest.raises(TypeError):
            # Call the parameter_sweep function
            _ = ps.parameter_sweep(
                m,
                sweep_params,
                build_outputs=None,
            )

    @pytest.mark.component
    def test_parameter_sweep_recover(self, model, tmp_path):
        comm = MPI.COMM_WORLD

        tmp_path = _get_rank0_path(comm, tmp_path)
        results_fname = os.path.join(tmp_path, "global_results")
        csv_results_file_name = str(results_fname) + ".csv"
        h5_results_file_name = str(results_fname) + ".h5"

        ps = ParameterSweep(
            optimize_function=_optimization,
            initialize_function=_reinitialize,
            initialize_kwargs={"slack_penalty": 10.0},
            csv_results_file_name=csv_results_file_name,
            h5_results_file_name=h5_results_file_name,
            debugging_data_dir=tmp_path,
            interpolate_nan_outputs=True,
        )
        # Call the parameter_sweep function
        _ = ps.parameter_sweep(
            build_model_for_tps,
            build_sweep_params_for_tps,
            build_outputs=None,
        )

        # NOTE: rank 0 "owns" tmp_path, so it needs to be
        #       responsible for doing any output file checking
        #       tmp_path can be deleted as soon as this method
        #       returns
        if ps.parallel_manager.is_root_process():
            # Check that the global results file is created
            assert os.path.isfile(csv_results_file_name)

            # Attempt to read in the data
            data = np.genfromtxt(csv_results_file_name, skip_header=1, delimiter=",")
            # Compare the last row of the imported data to truth
            truth_data = [0.9, 0.5, -11.0, 0.9, 0.5, 1.0, 1.0, 0.8, 0.5, 10.0, 2.0]
            assert np.allclose(data[-1], truth_data, equal_nan=True)

            # H5 dictionary test
            truth_dict = {
                "outputs": {
                    "fs.input[a]": {
                        "full_name": "fs.input[a]",
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array(
                            [0.1, 0.1, 0.1, 0.5, 0.5, 0.5, 0.9, 0.9, 0.9]
                        ),
                    },
                    "fs.input[b]": {
                        "full_name": "fs.input[b]",
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array(
                            [0.0, 0.25, 0.5, 0.0, 0.25, 0.5, 0.0, 0.25, 0.5]
                        ),
                    },
                    "fs.output[c]": {
                        "full_name": "fs.output[c]",
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array(
                            [
                                0.2,
                                0.2,
                                0.2,
                                1.0,
                                1.0,
                                1.0,
                                1.0,
                                1.0,
                                1.0,
                            ]
                        ),
                    },
                    "fs.output[d]": {
                        "full_name": "fs.output[d]",
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array(
                            [
                                0.0,
                                7.50000000e-01,
                                1.0,
                                9.77756334e-09,
                                7.50000010e-01,
                                1.0,
                                9.98605188e-09,
                                7.50000010e-01,
                                1.0,
                            ]
                        ),
                    },
                    "fs.performance": {
                        "full_name": "fs.performance",
                        "value": np.array(
                            [
                                0.2,
                                0.95,
                                1.2,
                                1.0,
                                1.75,
                                2.0,
                                1.0,
                                1.75,
                                2.0,
                            ]
                        ),
                    },
                    # removed upper bound tests, as in pracice
                    # user should not be chanaing upper and lower bounds
                    # during reinitilization and current frame work
                    # this bounds get updated based on existing model in each worker
                    "fs.slack[ab_slack]": {
                        "full_name": "fs.slack[ab_slack]",
                        "lower bound": 0,
                        "units": "None",
                        "value": np.array(
                            [
                                0.0,
                                0.0,
                                -9.77208577e-09,
                                -9.54438119e-09,
                                -9.54432513e-09,
                                -9.54432505e-09,
                                0.8,
                                0.8,
                                0.8,
                            ]
                        ),
                    },
                    "fs.slack[cd_slack]": {
                        "full_name": "fs.slack[cd_slack]",
                        "lower bound": 0,
                        "units": "None",
                        "value": np.array(
                            [
                                0.0,
                                0.0,
                                0.5,
                                -9.77756334e-09,
                                -9.77226571e-09,
                                0.5,
                                -9.98605188e-09,
                                -9.77226571e-09,
                                0.5,
                            ]
                        ),
                    },
                    "fs.slack_penalty": {
                        "full_name": "fs.slack_penalty",
                        "units": "None",
                        "value": np.array(
                            [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0]
                        ),
                    },
                    "objective": {
                        "full_name": "objective",
                        "value": np.array(
                            [
                                0.2,
                                0.95,
                                -3.79999979,
                                1.0,
                                1.75,
                                -2.99999979,
                                -6.99999978,
                                -6.24999979,
                                -10.99999979,
                            ]
                        ),
                    },
                },
                "solve_successful": [
                    True,
                    True,
                    True,
                    True,
                    True,
                    True,
                    True,
                    True,
                    True,
                ],
                "sweep_params": {
                    "fs.input[a]": {
                        "full_name": "fs.input[a]",
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array(
                            [0.1, 0.1, 0.1, 0.5, 0.5, 0.5, 0.9, 0.9, 0.9]
                        ),
                    },
                    "fs.input[b]": {
                        "full_name": "fs.input[b]",
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array(
                            [0.0, 0.25, 0.5, 0.0, 0.25, 0.5, 0.0, 0.25, 0.5]
                        ),
                    },
                },
            }

            read_dict = _read_output_h5(h5_results_file_name)
            _assert_dictionary_correctness(truth_dict, read_dict)
            _assert_h5_csv_agreement(csv_results_file_name, read_dict)

    @pytest.mark.component
    def test_parameter_sweep_bad_recover(self, model, tmp_path):
        comm = MPI.COMM_WORLD

        tmp_path = _get_rank0_path(comm, tmp_path)
        results_fname = os.path.join(tmp_path, "global_results")
        csv_results_file_name = str(results_fname) + ".csv"
        h5_results_file_name = str(results_fname) + ".h5"

        ps = ParameterSweep(
            optimize_function=_optimization,
            initialize_function=_bad_reinitialize,
            initialize_kwargs={"slack_penalty": 10.0},
            csv_results_file_name=csv_results_file_name,
            h5_results_file_name=h5_results_file_name,
            debugging_data_dir=tmp_path,
            interpolate_nan_outputs=True,
        )
        ps.config.log_model_states = True
        # Call the parameter_sweep function
        _ = ps.parameter_sweep(
            build_model_for_tps,
            build_sweep_params_for_tps,
            build_outputs=None,
        )

        # NOTE: rank 0 "owns" tmp_path, so it needs to be
        #       responsible for doing any output file checking
        #       tmp_path can be deleted as soon as this method
        #       returns
        if ps.parallel_manager.is_root_process():
            # Check that the global results file is created
            assert os.path.isfile(csv_results_file_name)

            # Attempt to read in the data
            data = np.genfromtxt(csv_results_file_name, skip_header=1, delimiter=",")

            # Compare the last row of the imported data to truth
            truth_data = [
                0.9,
                0.5,
                np.nan,
                0.9,
                0.5,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
            ]
            assert np.allclose(data[-1], truth_data, equal_nan=True)
            # H5 dictionary test
            truth_dict = {
                "outputs": {
                    "fs.input[a]": {
                        "full_name": "fs.input[a]",
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array(
                            [0.1, 0.1, 0.1, 0.5, 0.5, 0.5, 0.9, 0.9, 0.9]
                        ),
                    },
                    "fs.input[b]": {
                        "full_name": "fs.input[b]",
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array(
                            [0.0, 0.25, 0.5, 0.0, 0.25, 0.5, 0.0, 0.25, 0.5]
                        ),
                    },
                    "fs.output[c]": {
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array(
                            [0.2, 0.2, np.nan, 1.0, 1.0, np.nan, np.nan, np.nan, np.nan]
                        ),
                    },
                    "fs.output[d]": {
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array(
                            [
                                0.0,
                                0.75,
                                np.nan,
                                0.0,
                                0.75,
                                np.nan,
                                np.nan,
                                np.nan,
                                np.nan,
                            ]
                        ),
                    },
                    "fs.performance": {
                        "value": np.array(
                            [
                                0.2,
                                0.95,
                                np.nan,
                                1.0,
                                1.75,
                                np.nan,
                                np.nan,
                                np.nan,
                                np.nan,
                            ]
                        )
                    },
                    "fs.slack[ab_slack]": {
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 0,
                        "value": np.array(
                            [0.0, 0.0, np.nan, 0.0, 0.0, np.nan, np.nan, np.nan, np.nan]
                        ),
                    },
                    "fs.slack[cd_slack]": {
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 0,
                        "value": np.array(
                            [0.0, 0.0, np.nan, 0.0, 0.0, np.nan, np.nan, np.nan, np.nan]
                        ),
                    },
                    "objective": {
                        "value": np.array(
                            [
                                0.2,
                                0.95,
                                np.nan,
                                1.0,
                                1.75,
                                np.nan,
                                np.nan,
                                np.nan,
                                np.nan,
                            ]
                        )
                    },
                },
                "solve_successful": [
                    True,
                    True,
                    False,
                    True,
                    True,
                    False,
                    False,
                    False,
                    False,
                ],
                "sweep_params": {
                    "fs.input[a]": {
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array(
                            [0.1, 0.1, 0.1, 0.5, 0.5, 0.5, 0.9, 0.9, 0.9]
                        ),
                    },
                    "fs.input[b]": {
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array(
                            [0.0, 0.25, 0.5, 0.0, 0.25, 0.5, 0.0, 0.25, 0.5]
                        ),
                    },
                },
            }
            true_init_state, true_solved_state = dummy_kernel_logic(
                truth_dict["solve_successful"]
            )
            if has_mpi_peer_processes == False:
                assert true_init_state == ps.model_manager.initialized_states["state"]
                assert true_solved_state == ps.model_manager.solved_states["state"]
            read_dict = _read_output_h5(h5_results_file_name)
            _assert_dictionary_correctness(truth_dict, read_dict)
            _assert_h5_csv_agreement(csv_results_file_name, read_dict)

    @pytest.mark.component
    def test_parameter_sweep_force_initialize(self, model, tmp_path):
        results_fname = os.path.join(tmp_path, "global_results_force_initialize")
        csv_results_file_name = str(results_fname) + ".csv"
        h5_results_file_name = str(results_fname) + ".h5"
        ps = ParameterSweep(
            optimize_function=_optimization,
            initialize_before_sweep=True,
            initialize_function=_reinitialize,
            initialize_kwargs={"slack_penalty": 10.0},
            csv_results_file_name=csv_results_file_name,
            h5_results_file_name=h5_results_file_name,
            debugging_data_dir=tmp_path,
            interpolate_nan_outputs=False,
        )

        # Call the parameter_sweep function
        _ = ps.parameter_sweep(
            build_model_for_tps,
            build_sweep_params_for_tps,
            build_outputs=None,
        )

        # NOTE: rank 0 "owns" tmp_path, so it needs to be
        #       responsible for doing any output file checking
        #       tmp_path can be deleted as soon as this method
        #       returns
        if ps.parallel_manager.is_root_process():
            # Check that the global results file is created
            assert os.path.isfile(csv_results_file_name)

            # Attempt to read in the data
            data = np.genfromtxt(csv_results_file_name, skip_header=1, delimiter=",")

            # Compare the last row of the imported data to truth
            truth_data = [0.9, 0.5, -11.0, 0.9, 0.5, 1.0, 1.0, 0.8, 0.5, 10.0, 2.0]
            assert np.allclose(data[-1], truth_data, equal_nan=True)

            # H5 dictionary test
            truth_dict = {
                "outputs": {
                    "fs.input[a]": {
                        "full_name": "fs.input[a]",
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array(
                            [0.1, 0.1, 0.1, 0.5, 0.5, 0.5, 0.9, 0.9, 0.9]
                        ),
                    },
                    "fs.input[b]": {
                        "full_name": "fs.input[b]",
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array(
                            [0.0, 0.25, 0.5, 0.0, 0.25, 0.5, 0.0, 0.25, 0.5]
                        ),
                    },
                    "fs.output[c]": {
                        "full_name": "fs.output[c]",
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array(
                            [
                                0.2,
                                0.2,
                                0.2,
                                1.0,
                                1.0,
                                1.0,
                                1.0,
                                1.0,
                                1.0,
                            ]
                        ),
                    },
                    "fs.output[d]": {
                        "full_name": "fs.output[d]",
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array(
                            [
                                9.77899282e-09,
                                7.50000010e-01,
                                1.00000001e00,
                                9.77868510e-09,
                                7.50000010e-01,
                                1.00000001e00,
                                9.90100346e-09,
                                7.50000010e-01,
                                1.00000001e00,
                            ]
                        ),
                    },
                    "fs.performance": {
                        "full_name": "fs.performance",
                        "value": np.array(
                            [
                                0.2,
                                0.95,
                                1.2,
                                1.00,
                                1.75,
                                2.0,
                                1.0,
                                1.75,
                                2.0,
                            ]
                        ),
                    },
                    "fs.slack[ab_slack]": {
                        "full_name": "fs.slack[ab_slack]",
                        "lower bound": 0,
                        "units": "None",
                        # "upper bound": 0,
                        "value": np.array(
                            [
                                -9.77219059e-09,
                                -9.77208576e-09,
                                -9.77208577e-09,
                                -9.54438113e-09,
                                -9.54432509e-09,
                                -9.54432508e-09,
                                7.99999990e-01,
                                7.99999990e-01,
                                7.99999990e-01,
                            ]
                        ),
                    },
                    "fs.slack[cd_slack]": {
                        "full_name": "fs.slack[cd_slack]",
                        "lower bound": 0,
                        "units": "None",
                        # "upper bound": 0,
                        "value": np.array(
                            [
                                -9.77899282e-09,
                                -9.77226572e-09,
                                4.99999990e-01,
                                -9.77868510e-09,
                                -9.77226572e-09,
                                4.99999990e-01,
                                -9.90100346e-09,
                                -9.77226572e-09,
                                4.99999990e-01,
                            ]
                        ),
                    },
                    "fs.slack_penalty": {
                        "full_name": "fs.slack_penalty",
                        "units": "None",
                        "value": np.array(
                            [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0]
                        ),
                    },
                    "objective": {
                        "full_name": "objective",
                        "value": np.array(
                            [
                                0.20000022,
                                0.95000021,
                                -3.79999979,
                                1.00000021,
                                1.75000021,
                                -2.99999979,
                                -6.99999978,
                                -6.24999979,
                                -10.99999979,
                            ]
                        ),
                    },
                },
                "solve_successful": [
                    True,
                    True,
                    True,
                    True,
                    True,
                    True,
                    True,
                    True,
                    True,
                ],
                "sweep_params": {
                    "fs.input[a]": {
                        "full_name": "fs.input[a]",
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array(
                            [0.1, 0.1, 0.1, 0.5, 0.5, 0.5, 0.9, 0.9, 0.9]
                        ),
                    },
                    "fs.input[b]": {
                        "full_name": "fs.input[b]",
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array(
                            [0.0, 0.25, 0.5, 0.0, 0.25, 0.5, 0.0, 0.25, 0.5]
                        ),
                    },
                },
            }

            read_dict = _read_output_h5(h5_results_file_name)
            _assert_dictionary_correctness(truth_dict, read_dict)
            _assert_h5_csv_agreement(csv_results_file_name, read_dict)

    @pytest.mark.component
    def test_parameter_sweep_bad_sweep_update_before_initialize(self, model, tmp_path):
        ps = ParameterSweep(
            optimize_function=_optimization,
            initialize_before_sweep=True,
            initialize_function=None,
            initialize_kwargs=None,
            update_sweep_params_before_init=True,
        )
        m = model
        m.fs.slack_penalty = 1000.0
        m.fs.slack.setub(0)

        A = m.fs.input["a"]
        B = m.fs.input["b"]
        sweep_params = {A.name: (A, 0.1, 0.9, 3), B.name: (B, 0.0, 0.5, 3)}

        with pytest.raises(ValueError):
            # Call the parameter_sweep function
            ps.parameter_sweep(
                m,
                sweep_params,
                build_outputs=None,
            )

    @pytest.mark.component
    def test_parameter_sweep_bad_optimization_call(self, model, tmp_path):
        ps = ParameterSweep(
            optimize_function=_optimization,
            optimize_kwargs={"foo": "bar"},
        )

        m = model
        m.fs.slack_penalty = 1000.0
        m.fs.slack.setub(0)

        A = m.fs.input["a"]
        B = m.fs.input["b"]
        sweep_params = {A.name: (A, 0.1, 0.9, 3), B.name: (B, 0.0, 0.5, 3)}

        with pytest.raises(TypeError):
            # Call the parameter_sweep function
            ps.parameter_sweep(
                m,
                sweep_params,
                build_outputs=None,
            )

    @pytest.mark.component
    def test_parameter_sweep_bad_reinitialize_call(self, model, tmp_path):
        def reinit(a=42):
            pass

        ps = ParameterSweep(
            optimize_function=_optimization,
            initialize_before_sweep=True,
            initialize_function=reinit,
            initialize_kwargs={"foo": "bar"},
        )

        m = model
        m.fs.slack_penalty = 1000.0
        m.fs.slack.setub(0)

        A = m.fs.input["a"]
        B = m.fs.input["b"]
        sweep_params = {A.name: (A, 0.1, 0.9, 3), B.name: (B, 0.0, 0.5, 3)}

        with pytest.raises(TypeError):
            # Call the parameter_sweep function
            ps.parameter_sweep(
                m,
                sweep_params,
                build_outputs=None,
            )

    @pytest.mark.component
    def test_parameter_sweep_probe_fail(self, model, tmp_path):
        comm = MPI.COMM_WORLD

        tmp_path = _get_rank0_path(comm, tmp_path)
        results_fname = os.path.join(tmp_path, "global_results")
        csv_results_file_name = str(results_fname) + ".csv"
        h5_results_file_name = str(results_fname) + ".h5"

        ps = ParameterSweep(
            optimize_function=_optimization,
            optimize_kwargs={"relax_feasibility": True},
            probe_function=_bad_test_function,
            csv_results_file_name=csv_results_file_name,
            h5_results_file_name=h5_results_file_name,
            debugging_data_dir=tmp_path,
            interpolate_nan_outputs=True,
        )

        results_fname = os.path.join(tmp_path, "global_results")
        csv_results_file_name = str(results_fname) + ".csv"
        h5_results_file_name = str(results_fname) + ".h5"

        # Call the parameter_sweep function
        ps.parameter_sweep(
            build_model_for_tps,
            build_sweep_params_for_tps,
            build_outputs=build_outputs_for_tps,
            build_outputs_kwargs={
                "output_keys": {
                    "output_c": "fs.output[c]",
                    "output_d": "fs.output[d]",
                    "performance": "fs.performance",
                    "objective": "objective",
                }
            },
        )

        # NOTE: rank 0 "owns" tmp_path, so it needs to be
        #       responsible for doing any output file checking
        #       tmp_path can be deleted as soon as this method
        #       returns
        if ps.parallel_manager.is_root_process():
            # Check that the global results file is created
            assert os.path.isfile(csv_results_file_name)

            # Attempt to read in the data
            data = np.genfromtxt(csv_results_file_name, skip_header=1, delimiter=",")
            # Compare the last row of the imported data to truth
            truth_data = [
                0.9,
                0.5,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
            ]
            assert np.allclose(data[-1], truth_data, equal_nan=True)

            truth_dict = {
                "outputs": {
                    "output_c": {
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array([np.nan] * 9),
                    },
                    "output_d": {
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array(
                            np.array([np.nan] * 9),
                        ),
                    },
                    "performance": {"value": np.array([np.nan] * 9)},
                    "objective": {
                        "value": np.array(
                            np.array([np.nan] * 9),
                        )
                    },
                },
                "solve_successful": [False] * 9,
                "sweep_params": {
                    "fs.input[a]": {
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array(
                            [0.1, 0.1, 0.1, 0.5, 0.5, 0.5, 0.9, 0.9, 0.9]
                        ),
                    },
                    "fs.input[b]": {
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array(
                            [0.0, 0.25, 0.5, 0.0, 0.25, 0.5, 0.0, 0.25, 0.5]
                        ),
                    },
                },
            }

            read_dict = _read_output_h5(h5_results_file_name)
            _assert_dictionary_correctness(truth_dict, read_dict)
            _assert_h5_csv_agreement(csv_results_file_name, read_dict)

    @pytest.mark.component
    def test_parameter_sweep_function(self, model, tmp_path):
        comm = MPI.COMM_WORLD
        rank = comm.rank
        tmp_path = _get_rank0_path(comm, tmp_path)

        m = model
        m.fs.slack_penalty = 1000.0
        m.fs.slack.setub(0)

        A = m.fs.input["a"]
        B = m.fs.input["b"]
        sweep_params = {A.name: (A, 0.1, 0.9, 3), B.name: (B, 0.0, 0.5, 3)}
        outputs = {
            "output_c": m.fs.output["c"],
            "output_d": m.fs.output["d"],
            "performance": m.fs.performance,
        }

        results_fname = os.path.join(tmp_path, "global_results")
        csv_results_file_name = str(results_fname) + ".csv"
        h5_results_file_name = str(results_fname) + ".h5"

        # Call the parameter_sweep function
        parameter_sweep(
            m,
            sweep_params,
            build_outputs=outputs,
            csv_results_file_name=csv_results_file_name,
            h5_results_file_name=h5_results_file_name,
            optimize_function=_optimization,
            debugging_data_dir=tmp_path,
            interpolate_nan_outputs=True,
        )

        # NOTE: rank 0 "owns" tmp_path, so it needs to be
        #       responsible for doing any output file checking
        #       tmp_path can be deleted as soon as this method
        #       returns
        if rank == 0:
            # Check that the global results file is created
            assert os.path.isfile(csv_results_file_name)
            assert os.path.isfile(
                os.path.join(tmp_path, "interpolated_global_results.csv")
            )

            # Check that all local output files have been created
            for k in range(comm.size):
                assert os.path.isfile(
                    os.path.join(tmp_path, f"local_results_{k:03}.h5")
                )
                assert os.path.isfile(
                    os.path.join(tmp_path, f"local_results_{k:03}.csv")
                )

            # Attempt to read in the data
            data = np.genfromtxt(csv_results_file_name, skip_header=1, delimiter=",")

            # Compare the last row of the imported data to truth
            truth_data = [0.9, 0.5, np.nan, np.nan, np.nan]
            assert np.allclose(data[-1], truth_data, equal_nan=True)

        # Check for the h5 output
        if rank == 0:
            truth_dict = {
                "outputs": {
                    "output_c": {
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array(
                            [0.2, 0.2, np.nan, 1.0, 1.0, np.nan, np.nan, np.nan, np.nan]
                        ),
                    },
                    "output_d": {
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array(
                            [
                                0.0,
                                0.75,
                                np.nan,
                                0.0,
                                0.75,
                                np.nan,
                                np.nan,
                                np.nan,
                                np.nan,
                            ]
                        ),
                    },
                    "performance": {
                        "value": np.array(
                            [
                                0.2,
                                0.95,
                                np.nan,
                                1.0,
                                1.75,
                                np.nan,
                                np.nan,
                                np.nan,
                                np.nan,
                            ]
                        )
                    },
                },
                "solve_successful": [
                    True,
                    True,
                    False,
                    True,
                    True,
                    False,
                    False,
                    False,
                    False,
                ],
                "sweep_params": {
                    "fs.input[a]": {
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array(
                            [0.1, 0.1, 0.1, 0.5, 0.5, 0.5, 0.9, 0.9, 0.9]
                        ),
                    },
                    "fs.input[b]": {
                        "lower bound": 0,
                        "units": "None",
                        "upper bound": 1,
                        "value": np.array(
                            [0.0, 0.25, 0.5, 0.0, 0.25, 0.5, 0.0, 0.25, 0.5]
                        ),
                    },
                },
            }

            read_dict = _read_output_h5(h5_results_file_name)
            _assert_dictionary_correctness(truth_dict, read_dict)
            _assert_h5_csv_agreement(csv_results_file_name, read_dict)

            # Check if there is a text file created
            import ast

            truth_txt_dict = {
                "outputs": [
                    "output_c",
                    "output_d",
                    "performance",
                ],
                "sweep_params": ["fs.input[a]", "fs.input[b]"],
            }

            txt_fpath = h5_results_file_name + ".txt"
            assert os.path.exists(txt_fpath)
            f = open(txt_fpath)
            f_contents = f.read()
            read_txt_dict = ast.literal_eval(f_contents)
            assert read_txt_dict == truth_txt_dict

    @pytest.mark.unit
    def test_parameter_sweep_custom_do_param_sweep(self, model, tmp_path):
        def custom_do_param_sweep(model, sweep_params, outputs, local_values, **kwargs):
            return kwargs

        custom_kwargs = {"val1": 2.0}
        ps = ParameterSweep(
            custom_do_param_sweep=custom_do_param_sweep,
            custom_do_param_sweep_kwargs=custom_kwargs,
        )

        m = model
        m.fs.slack_penalty = 1000.0
        m.fs.slack.setub(0)

        A = m.fs.input["a"]
        B = m.fs.input["b"]
        sweep_params = {A.name: (A, 0.1, 0.9, 3), B.name: (B, 0.0, 0.5, 3)}
        outputs = {
            "output_c": m.fs.output["c"],
            "output_d": m.fs.output["d"],
            "performance": m.fs.performance,
            "objective": m.objective,
        }

        assert ps.config.custom_do_param_sweep is not None
        return_dict = ps.config.custom_do_param_sweep(
            model,
            sweep_params,
            outputs,
            0.0,
            val1=2.0,
        )
        assert return_dict == custom_kwargs


def _optimization(m, relax_feasibility=False):
    if relax_feasibility:
        m.fs.slack.setub(None)

    solver = pyo.SolverFactory("ipopt")
    results = solver.solve(m)

    return results


def _initialize_with_added_var(m):
    if (abs(m.fs.input["a"].value - 0.5) < 1e-6) and abs(
        m.fs.input["b"].value - 0.5
    ) < 1e-6:
        m.fs.test_var = pyo.Var(initialize=5, bounds=(0, 10))
        m.fs.test_var.fix()


def _reinitialize(m, slack_penalty=10.0):
    m.fs.slack.setub(None)
    m.fs.slack_penalty.value = slack_penalty


def _bad_reinitialize(m, **kwargs):
    pass


def _get_rank0_path(comm, tmp_path):
    if comm is None:
        return tmp_path
    return comm.bcast(tmp_path, root=0)


def _good_test_function(m):
    return True


def _bad_test_function(m):
    return False


def _assert_dictionary_correctness(truth_dict, test_dict, rtol=1e-05, atol=1e-08):
    assert truth_dict.keys() == test_dict.keys()

    for key, item in truth_dict.items():
        if key in ["sweep_params", "outputs"]:
            for subkey, subitem in item.items():
                for subsubkey, subsubitem in subitem.items():
                    if subsubkey == "value":
                        assert np.allclose(
                            test_dict[key][subkey]["value"],
                            subitem["value"],
                            equal_nan=True,
                            rtol=rtol,
                            atol=atol,
                        )
                    else:
                        if (
                            subsubkey in ["upper bound", "lower bound"]
                            and test_dict[key][subkey][subsubkey] is None
                        ):
                            assert subsubitem in [np.finfo("d").min, np.finfo("d").max]
                        else:
                            assert subsubitem == test_dict[key][subkey][subsubkey]
        elif key == "solve_successful":
            assert item == test_dict[key]
        elif key in ["nominal_idx", "differential_idx"]:
            assert np.allclose(
                test_dict[key], item, equal_nan=True, rtol=rtol, atol=atol
            )


def _assert_h5_csv_agreement(csv_filename, h5_dict):
    csv_header = _build_header_list_from_csv(csv_filename)
    csv_data = np.genfromtxt(csv_filename, skip_header=1, delimiter=",")

    for output_type in ["sweep_params", "outputs"]:
        for key in h5_dict[output_type].keys():
            # Find the matching label in the CSV file
            idx = csv_header.index(key)
            assert np.allclose(
                h5_dict[output_type][key]["value"], csv_data[:, idx], equal_nan=True
            )


def _build_header_list_from_csv(csv_filename):
    with open(csv_filename, "r") as fp:
        csv_header_raw = fp.readline()

    csv_header_raw = csv_header_raw.replace("# ", "").strip()
    csv_header_raw = csv_header_raw.split(",")

    csv_header = []

    for k, item in enumerate(csv_header_raw):
        if k == 0:
            csv_header.append(item)

        else:
            if csv_header[-1].count("[") == csv_header[-1].count("]"):
                csv_header.append(item)
            else:
                csv_header[-1] += f",{item}"

    return csv_header


def _read_output_h5(filevar):
    if isinstance(filevar, str):
        f = h5py.File(filevar, "r")
    elif isinstance(filevar, h5py.Group):
        f = filevar
    else:
        raise RuntimeError(f"Unrecognized type for filepath {type(filevar)}")

    l1_keys = list(f.keys())
    output_dict = {}
    for key in l1_keys:  # Input or Output
        if key in ["sweep_params", "outputs"]:  #  "solve_successful":
            output_dict[key] = {}
            l2_keys = list(f[key].keys())
            for subkey in l2_keys:  # Variable name
                output_dict[key][subkey] = {}
                l3_keys = list(f[key][subkey].keys())
                for subsubkey in l3_keys:  # variable metadata
                    output_dict[key][subkey][subsubkey] = f[key][subkey][subsubkey][()]
                    if subsubkey == "units" or subsubkey == "full_name":
                        # The strings are recovered in bytes. we choose to convert it to utf-8
                        output_dict[key][subkey][subsubkey] = output_dict[key][subkey][
                            subsubkey
                        ].decode("utf-8")
        elif key == "solve_successful":
            output_dict[key] = list(f[key]["solve_successful"][()])
        elif key in ["nominal_idx", "differential_idx"]:
            output_dict[key] = f[key][key][()]

    if isinstance(filevar, str):
        f.close()

    return output_dict
