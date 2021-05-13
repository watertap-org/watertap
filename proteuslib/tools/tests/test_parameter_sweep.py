##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################

import pytest
import os
import numpy as np

from pyomo.environ import (ConcreteModel,
                           TerminationCondition,
                           SolverStatus,
                           value,
                           Objective,
                           SolverFactory)

from idaes.core import FlowsheetBlock
from idaes.core.util.testing import get_default_solver
from idaes.core.util.scaling import calculate_scaling_factors

from proteuslib.unit_models.reverse_osmosis_0D import ReverseOsmosis0D
import proteuslib.property_models.NaCl_prop_pack as props
from proteuslib.tools.parameter_sweep import (_init_mpi,
                                               _build_and_divide_combinations,
                                               _update_model_values,
                                               _aggregate_results,
                                               parameter_sweep)

# -----------------------------------------------------------------------------

class TestParallelManager():
    @pytest.fixture(scope="class")
    def RO_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = ReverseOsmosis0D(default={
            "property_package": m.fs.properties,
            "has_pressure_change": True, })

        # fully specify system
        feed_flow_mass = 1
        feed_mass_frac_NaCl = 0.035
        feed_pressure = 50e5
        feed_temperature = 273.15 + 25
        membrane_pressure_drop = 3e5
        membrane_area = 50
        A = 4.2e-12
        B = 3.5e-8
        pressure_atmospheric = 101325

        feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
        m.fs.unit.inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl'].fix(
            feed_flow_mass * feed_mass_frac_NaCl)
        m.fs.unit.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O'].fix(
            feed_flow_mass * feed_mass_frac_H2O)
        m.fs.unit.inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.deltaP.fix(-membrane_pressure_drop)
        m.fs.unit.area.fix(membrane_area)
        m.fs.unit.A_comp.fix(A)
        m.fs.unit.B_comp.fix(B)
        m.fs.unit.permeate.pressure[0].fix(pressure_atmospheric)
        return m

    @pytest.mark.unit
    def test_init_mpi(self):
        comm, rank, num_procs = _init_mpi()

        assert type(rank) == int
        assert type(num_procs) == int
        assert 0 <= rank < num_procs

    @pytest.mark.component
    def test_build_and_divide_combinations(self):
        comm, rank, num_procs = _init_mpi()

        range_A = [0.0, 10.0]
        range_B = [1.0, 20.0]
        range_C = [2.0, 30.0]

        nn_A = 4
        nn_B = 5
        nn_C = 6

        param_dict = dict()
        param_dict['var_A'] = (None, range_A[0], range_A[1], nn_A)
        param_dict['var_B']  = (None, range_B[0], range_B[1], nn_B)
        param_dict['var_C']  = (None, range_C[0], range_C[1], nn_C)

        local_combo_array, full_combo_array = _build_and_divide_combinations(param_dict, rank, num_procs)

        test = np.array_split(full_combo_array, num_procs, axis=0)[rank]

        assert np.shape(full_combo_array)[0] == nn_A*nn_B*nn_C
        assert np.shape(full_combo_array)[1] == 3
        assert np.shape(local_combo_array)[1] == 3

        assert test[-1, 0] == pytest.approx(local_combo_array[-1, 0])
        assert test[-1, 1] == pytest.approx(local_combo_array[-1, 1])
        assert test[-1, 2] == pytest.approx(local_combo_array[-1, 2])

        if rank == 0:
            assert test[0, 0] == pytest.approx(range_A[0])
            assert test[0, 1] == pytest.approx(range_B[0])
            assert test[0, 2] == pytest.approx(range_C[0])

        if rank == num_procs - 1:
            assert test[-1, 0] == pytest.approx(range_A[1])
            assert test[-1, 1] == pytest.approx(range_B[1])
            assert test[-1, 2] == pytest.approx(range_C[1])

    @pytest.mark.component
    def test_update_model_values(self, RO_frame):
        m = RO_frame

        param_dict = dict()
        param_dict['pressure'] = (m.fs.unit.inlet.pressure[0], None, None, None)
        param_dict['area'] = (m.fs.unit.area, None, None, None)

        assert hasattr(m.fs.unit.inlet, 'pressure')
        assert hasattr(m.fs.unit, 'area')

        original_pressure = value(m.fs.unit.inlet.pressure[0])
        original_area = value(m.fs.unit.area)

        new_values = [1.1*original_pressure, 1.1*original_area]

        _update_model_values(m, param_dict, new_values)

        assert value(m.fs.unit.inlet.pressure[0]) == pytest.approx(new_values[0])
        assert value(m.fs.unit.area) == pytest.approx(new_values[1])

    @pytest.mark.unit
    def test_aggregate_results(self):
        comm, rank, num_procs = _init_mpi()

        # print('Rank %d, num_procs %d' % (rank, num_procs))

        nn = 5
        np.random.seed(1)
        local_results = (rank+1)*np.random.rand(nn, 2)
        global_values = np.random.rand(nn*num_procs, 4)

        global_results = _aggregate_results(local_results, global_values, comm, num_procs)

        assert np.shape(global_results)[1] == np.shape(local_results)[1]
        assert np.shape(global_results)[0] == np.shape(global_values)[0]

        if rank == 0:
            assert global_results[0, 0] == pytest.approx(local_results[0, 0])
            assert global_results[0, 1] == pytest.approx(local_results[0, 1])
            assert global_results[-1, 0] == pytest.approx(num_procs*local_results[-1, 0])
            assert global_results[-1, 1] == pytest.approx(num_procs*local_results[-1, 1])

    @pytest.mark.component
    def test_parameter_sweep(self, RO_frame):
        def optimization(m, objective, **kwargs):
            if objective == 'Area':
                m.fs.unit.area.unfix()
                m.fs.unit.area.setlb(20)
                m.fs.unit.area.setub(60)
                m.fs.objective = Objective(expr=m.fs.unit.area)


            m.fs.unit.B_comp.unfix()
            m.fs.unit.B_comp.setlb(3.5e-8)
            m.fs.unit.B_comp.setub(3.5e-8 * 1e3)
        
            m.fs.unit.A_comp.unfix()
            m.fs.unit.A_comp.setlb(4.2e-12)
            m.fs.unit.A_comp.setub(4.2e-12 * 1e3)

            # solver = SolverFactory('ipopt')
            solver = SolverFactory('ipopt')
            solver.options = {'nlp_scaling_method': 'user-scaling'}
            results = solver.solve(m, tee=True)

            return m

        def display_metrics(m, outputs, **kwargs):
            metrics = dict()

            metrics['A_comp'] = value(m.fs.unit.A_comp[0.0, 'H2O'])
            metrics['B_comp'] = value(m.fs.unit.B_comp[0.0, 'NaCl'])
            metrics['Area'] = value(m.fs.unit.area)

            output_data = []

            for k in outputs:
                output_data.append(metrics[k])

            return output_data

        comm, rank, num_procs = _init_mpi()

        m = RO_frame

        calculate_scaling_factors(m)

        solver = get_default_solver()
        solver.options = {'nlp_scaling_method': 'user-scaling'}
        results = solver.solve(m)
        # print(results)

        assert results.solver.termination_condition == \
               TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

        sweep_params = dict()
        sweep_params['pressure_drop'] = (m.fs.unit.deltaP, 2e5, 4e5, 3)
        outputs = ['A_comp', 'B_comp', 'Area']

        # Call the parameter_sweep function
        parameter_sweep(m, sweep_params, outputs, objective='Area', 
            output_dir='pytest_output', mpi_comm=comm, optimization=optimization,
            display_metrics=display_metrics, save_debugging_data=True)

        # Check that the global results file is created
        assert os.path.isfile('pytest_output/global_results.csv')

        if rank == 0:
            # Check that all local output files have been created
            for k in range(num_procs):
                assert os.path.isfile('pytest_output/local_results_%03d.csv' % (k))

            # Attempt to read in the data
            data = np.genfromtxt('pytest_output/global_results.csv', skip_header=1, delimiter=',')

            truth_data = [4.000000e+05, 1.469657e-11, 2.922988e-05, 2.000000e+01]

            # Compare the last row of the imported data to truth
            for k in range(len(truth_data)):
                assert data[-1, k] == pytest.approx(truth_data[k])
