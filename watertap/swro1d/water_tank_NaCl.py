import pyomo.environ as pyo
from idaes.core import FlowsheetBlock
from idaes.models_extra.power_generation.unit_models.watertank import WaterTank
from watertap.property_models.NaCl_prop_pack import NaClParameterBlock
# import idaes.core.util.scaling as iscale
from idaes.core.solvers import petsc
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.core.util.scaling as iscale

import idaes.logger as idaeslog
from IPython.display import SVG, display
from idaes.core.util import DiagnosticsToolbox
from enum import Enum
import numpy as np
import matplotlib.pyplot as plt

class OperatingScenario(Enum):
    steady_state = 1
    dynamic = 2

def main():
    operating_scenario = OperatingScenario.steady_state

    m = pyo.ConcreteModel()
    t_start = 0 * 60 # s
    t_end = 1 * 3 # s
    dt_set = [t_start, t_end]
    time_set = [sum(dt_set[:j]) for j in range(len(dt_set)+1)]
    # time_set = np.arange(100)*0.0001
    
    if operating_scenario == OperatingScenario.steady_state:
        m.fs = FlowsheetBlock(dynamic=False)
    else:
        m.fs = FlowsheetBlock(
            dynamic=True,
            time_set=time_set,
            time_units=pyo.units.s
        )
    m.fs.properties = NaClParameterBlock()

    m.fs.wtank = WaterTank(
        property_package=m.fs.properties,
        tank_type='simple_tank',
        # has_holdup='True',
        has_pressure_change=True,
        has_heat_transfer=False
        )
    print("Degrees of Freedom before = ", degrees_of_freedom(m))

    # fix inputs for simple tank
    m.fs.wtank.tank_cross_sect_area.fix(1.14)  # tank cross sectional area
    # m.fs.wtank.tank_level[:].fix(0.6)  # tank level
    m.fs.wtank.volume.fix(0.684)
    m.fs.wtank.inlet.flow_mass_phase_comp[:, "Liq", "H2O"].fix(0.965)
    m.fs.wtank.inlet.flow_mass_phase_comp[:, "Liq", "NaCl"].fix(0.035)
    m.fs.wtank.inlet.pressure.fix(101325)
    m.fs.wtank.inlet.temperature.fix(298.15)
    
    print("Degrees of Freedom after specifying = ", degrees_of_freedom(m))

    # m.fs.wtank.display()

    # iscale.calculate_scaling_factors(m)
    m.fs.wtank.initialize(outlvl=3)

    # Solve the simulation using ipopt
    # Note: If the degrees of freedom = 0, we have a square problem
    opt = pyo.SolverFactory("ipopt")
    optarg = {"tol": 1e-7, "linear_solver": "mumps", "max_iter": 40}
    opt.options = optarg
    solve_status = opt.solve(m, tee=True)

    print("Degrees of Freedom after initialization = ", degrees_of_freedom(m))
    m.fs.wtank.display()
    m.fs.wtank.report()
    # print(pyo.value(m.fs.wtank.outlet.enth_mol))
    # m.fs.wtank.outlet.flow_mol[0].fix()
    # m.fs.wtank.outlet.pressure[0].fix()
    # m.fs.wtank.outlet.enth_mol[0].fix()
    # print("Degrees of Freedom after fixing outlet.enth_mol = ", degrees_of_freedom(m))
    # assert False
    
    print('DeltaP = ', pyo.value(m.fs.wtank.deltaP[0]))
    print('Level = ', pyo.value(m.fs.wtank.tank_level[0]))
    print('V = ', pyo.value(m.fs.wtank.volume[0]))
    # print('Outlet H2O, NaCl, pressure: ', pyo.value(m.fs.wtank.outlet[:].flow_mass_phase_comp[0.0, 'Liq', 'H2O'])[0], pyo.value(m.fs.wtank.outlet[:].flow_mass_phase_comp[0.0, 'Liq', 'NaCl'])[0], pyo.value(m.fs.wtank.outlet[:].pressure[0])[0])
    print('Inlet flow enthalpy, Outlet flow enthalpy, volumetric flow rate, pressure: ', pyo.value(m.fs.wtank.control_volume.properties_in[:].enth_flow)[0], pyo.value(m.fs.wtank.control_volume.properties_out[:].enth_flow)[0], pyo.value(m.fs.wtank.control_volume.properties_out[:].flow_vol)[0], pyo.value(m.fs.wtank.outlet[:].pressure[0])[0])

    # assert False

    if operating_scenario == OperatingScenario.dynamic:
        # m.fs.wtank.display()
        # m.fs.wtank.set_initial_condition()
        print(pyo.value(m.fs.wtank.control_volume.material_accumulation[:, :, :]))
        # m.fs.wtank.display()
        time_nfe = len(m.fs.time) - 1
        

        pyo.TransformationFactory("dae.finite_difference").apply_to(
            m.fs, nfe=time_nfe, wrt=m.fs.time, scheme="BACKWARD"
        )
        scaling_log = idaeslog.getLogger("idaes.core.util.scaling")
        scaling_log.setLevel(idaeslog.ERROR)
        iscale.calculate_scaling_factors(m)

        idaeslog.solver_log.tee = True
        results = petsc.petsc_dae_by_time_element(
            m,
            time=m.fs.time,
            keepfiles=True,
            symbolic_solver_labels=True,
            ts_options={
                "--ts_type": "beuler",
                # "-ts_arkimex_type": "1bee",
                "--ts_dt": 0.1,
                "--ts_rtol": 1e-3,
                # "--ts_adapt_clip":"0.001,3600",
                # "--ksp_monitor":"",
                "--ts_adapt_dt_min": 1e-3,
                "--ts_adapt_dt_max": 3600,
                "--snes_type": "newtontr",
                # "--ts_max_reject": 200,
                "--ts_monitor": "",
                "-ts_adapt_monitor": "",
                # "--snes_monitor":"",
                "-snes_converged_reason": "",
                # "-ksp_monitor_true_residual": "",
                # "-ksp_converged_reason": "",
                # "-snes_test_jacobian": "",
                "snes_grid_sequence": "",
                "-pc_type": "lu",
                # "-mat_view": "",
                "--ts_save_trajectory": 1,
                "--ts_trajectory_type": "visualization",
                "--ts_max_snes_failures": 25,
                # "--show_cl":"",
                "-snes_max_it": 50,
                "-snes_rtol": 0,
                "-snes_stol": 0,
                "-snes_atol": 1e-6,
            },
            skip_initial=True,
            initial_solver="ipopt",
            initial_solver_options={
                "constr_viol_tol": 1e-8,
                "nlp_scaling_method": "user-scaling",
                "linear_solver": "mumps",
                "OF_ma57_automatic_scaling": "yes",
                "max_iter": 300,
                "tol": 1e-8,
                "halt_on_ampl_error": "no",
            },
        )
        for result in results.results:
            pyo.assert_optimal_termination(result)
        
        print(pyo.value(m.fs.wtank.control_volume.material_accumulation[:, :, :]))
        traj = results.trajectory
        time_set = m.fs.time.ordered_data()
        tf = time_set[-1]
        results_dict = {
            "time": np.array(traj.time),
            "tank_level": np.array(traj.vecs[str(m.fs.wtank.tank_level[tf])]),
            "outlet.flow_mol": np.array(traj.vecs[str(m.fs.wtank.outlet.flow_mol[tf])]),
            "outlet.enth_mol": np.array(traj.vecs[str(m.fs.wtank.outlet.enth_mol[tf])]),
            "outlet.pressure": np.array(traj.vecs[str(m.fs.wtank.outlet.pressure[tf])])
        }
        for key, value in results_dict.items():
            # Turn n by 1 arrays in into vectors
            results_dict[key] = np.squeeze(value)
        time = results_dict["time"]
        print(time)
        print(np.array(traj.vecs[str(m.fs.wtank.outlet.flow_mol[tf])]))

        fig = plt.figure()
        ax = fig.subplots(4, 1, sharex = True)
        ax[0].plot(time, results_dict["tank_level"])
        ax[1].plot(time, results_dict["outlet.flow_mol"])
        ax[2].plot(time, results_dict["outlet.enth_mol"])
        ax[3].plot(time, results_dict["outlet.pressure"])
        ax[0].set_xlim(time[0], time[-1])
        # ax.set_ylim((0.65, 1.45))
        ax[3].set_xlabel("Time (s)", fontsize=14)
        ax[0].set_ylabel("tank_level (m)", fontsize=14)
        ax[1].set_ylabel("outlet.flow_mol (m)", fontsize=14)
        ax[2].set_ylabel("outlet.enth_mol (m)", fontsize=14)
        ax[3].set_ylabel("outlet.pressure (m)", fontsize=14)
        
        # ax.set_title("SOEC Voltage", fontsize=16)
        plt.show()

if __name__ == "__main__":
    m = main()