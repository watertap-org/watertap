import watertap.examples.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery as ro_erd
from watertap.tools.parameter_sweep.parameter_sweep import (
    ParameterSweep,
    RecursiveParameterSweep,
)
from watertap.tools.parameter_sweep.parameter_sweep_reader import ParameterSweepReader

from watertap.tools.parameter_sweep.parameter_sweep_differential import (
    DifferentialParameterSweep,
)


from idaes.core.solvers import get_solver


def ro_build(**kwargs):
    m = ro_erd.build(**kwargs)
    # ro_erd.set_operating_conditions(m)
    # ro_erd.optimize_set_up(m)
    return m


def ro_init(m, solver=None, **kwargs):
    # make sure these are fixed by init routine instad..
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].unfix()
    m.fs.RO.recovery_mass_phase_comp[0, "Liq", "H2O"].unfix()

    ro_erd.set_operating_conditions(
        m, water_recovery=0.5, over_pressure=0.3, solver=solver
    )
    ro_erd.initialize_system(m, solver=solver)
    ro_erd.optimize_set_up(m)
    print(
        'm.fs.feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"]',
        m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"].value,
    )
    m.fs.feed.properties[0].display()


def ro_solve(m, solver=None, **kwargs):
    result = solver.solve(m)
    return result


if __name__ == "__main__":
    sweep_params = {
        "feed_mass_nacl": {
            "type": "LinearSample",
            "param": "fs.feed.properties[0].flow_mass_phase_comp[Liq,NaCl]",
            "lower_limit": 0.03,
            "upper_limit": 0.04,
            "num_samples": 3,
        },
        "ro_recovery": {
            "type": "LinearSample",
            "param": "fs.RO.recovery_mass_phase_comp[0,Liq,H2O]",
            "lower_limit": 0.3,
            "upper_limit": 0.5,
            "num_samples": 3,
        },
    }
    solver = get_solver()
    opt_kwargs = {"solver": solver}
    ps_kwargs = {}
    ps_kwargs["csv_results_file_name"] = None
    ps_kwargs["h5_results_file_name"] = "test.h5"
    ps_kwargs["h5_parent_group_name"] = "test/1"

    ps_kwargs["optimize_function"] = ro_solve
    ps_kwargs["optimize_kwargs"] = opt_kwargs

    ps_kwargs["reinitialize_function"] = ro_init
    ps_kwargs["reinitialize_kwargs"] = opt_kwargs
    ps_kwargs["reinitialize_before_sweep"] = False

    ps_kwargs["custom_do_param_sweep"] = None
    ps_kwargs["custom_do_param_sweep_kwargs"] = None

    ps_kwargs["probe_function"] = None

    ps_kwargs["number_of_subprocesses"] = 1

    ps = ParameterSweep(**ps_kwargs)
    ps.parameter_sweep(
        ro_build,
        ParameterSweepReader()._dict_to_params,
        build_outputs=None,
        num_samples=3,
        build_model_kwargs={},
        build_sweep_params_kwargs={"input_dict": sweep_params},
    )
