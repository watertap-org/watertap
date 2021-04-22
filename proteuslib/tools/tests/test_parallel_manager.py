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
from proteuslib.tools.parallel_manager import (init_mpi,
                                               build_and_divide_combinations,
                                               update_model_values,
                                               aggregate_results,
                                               run_param_sweep)

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
        comm, rank, num_procs = init_mpi()

        assert type(rank) == int
        assert type(num_procs) == int
        assert 0 <= rank < num_procs

    @pytest.mark.component
    def test_build_and_divide_combinations(self):
        comm, rank, num_procs = init_mpi()

        nn_A = 4
        nn_B = 5
        nn_C = 6

        param_dict = dict()
        param_dict['var_A'] = (None, 0.0, 1.0, nn_A)
        param_dict['var_B']  = (None, 1.0, 2.0, nn_B)
        param_dict['var_C']  = (None, 2.0, 3.0, nn_C)

        local_combo_array, full_combo_array = build_and_divide_combinations(param_dict, rank, num_procs)

        test = np.array_split(full_combo_array, num_procs, axis=0)[rank]

        assert np.shape(full_combo_array)[0] == nn_A*nn_B*nn_C
        assert np.shape(full_combo_array)[1] == 3
        assert np.shape(local_combo_array)[1] == 3

        assert test[-1, 0] == pytest.approx(local_combo_array[-1, 0])
        assert test[-1, 1] == pytest.approx(local_combo_array[-1, 1])
        assert test[-1, 2] == pytest.approx(local_combo_array[-1, 2])

        if rank == 0:
            assert test[0, 0] == pytest.approx(0.0)
            assert test[0, 1] == pytest.approx(1.0)
            assert test[0, 2] == pytest.approx(2.0)

        if rank == num_procs - 1:
            assert test[-1, 0] == pytest.approx(1.0)
            assert test[-1, 1] == pytest.approx(2.0)
            assert test[-1, 2] == pytest.approx(3.0)

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

        update_model_values(m, param_dict=param_dict, values=new_values)

        assert value(m.fs.unit.inlet.pressure[0]) == pytest.approx(new_values[0])
        assert value(m.fs.unit.area) == pytest.approx(new_values[1])

    @pytest.mark.unit
    def test_aggregate_results(self):
        comm, rank, num_procs = init_mpi()

        nn = 5
        np.random.seed(1)
        local_results = (rank+1)*np.random.rand(nn, 2)
        global_values = np.random.rand(nn*num_procs, 4)

        global_results = aggregate_results(local_results, global_values, comm, num_procs)

        assert np.shape(global_results)[1] == np.shape(local_results)[1]
        assert np.shape(global_results)[0] == np.shape(global_values)[0]

        if rank == 0:
            assert global_results[0, 0] == pytest.approx(local_results[0, 0])
            assert global_results[0, 1] == pytest.approx(local_results[0, 1])
            assert global_results[-1, 0] == pytest.approx(num_procs*local_results[-1, 0])
            assert global_results[-1, 1] == pytest.approx(num_procs*local_results[-1, 1])

    @pytest.mark.component
    def test_run_param_sweep(self, RO_frame):
        # m, sweep_params, outputs, output_dir='output', mpi_comm=None, num_stages=2, optimization=None, display_metrics=None

        m = RO_frame

        calculate_scaling_factors(m)

        solver = get_default_solver()
        solver.options = {'nlp_scaling_method': 'user-scaling'}
        results = solver.solve(m)
        print(results)

        assert results.solver.termination_condition == \
               TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

        m.fs.objective = Objective(expr=m.fs.unit.area)

        m.fs.unit.B_comp.unfix()
        m.fs.unit.B_comp.setlb(3.5e-8)
        m.fs.unit.B_comp.setub(3.5e-8 * 1e3)

        # solver = SolverFactory('ipopt')
        solver = SolverFactory('ipopt')
        solver.options = {'nlp_scaling_method': 'user-scaling'}
        results = solver.solve(m, tee=True)

        assert results.solver.termination_condition == TerminationCondition.optimal
