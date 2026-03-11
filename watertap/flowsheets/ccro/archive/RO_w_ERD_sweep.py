import pandas as pd
from pyomo.environ import units as pyunits
from parameter_sweep import (
    LinearSample,
    parameter_sweep,
)
import watertap.flowsheets.ccro.ro_w_erd.RO_w_ERD_stage as ro
from watertap.core.solvers import get_solver

solver = get_solver()


def build_and_solve(n_stages=2, salt_mass_frac=35e-3, **kwargs):

    m = ro.run_n_stage_system(
        n_stages=n_stages, salt_mass_frac=salt_mass_frac, **kwargs
    )

    for n, stage in m.fs.stage.items():
        if stage.add_pump:
            stage.pump.control_volume.properties_out[0].pressure.unfix()
        stage.RO.area.unfix()

    m.fs.system_recovery.unfix()

    return m


def build_sweep_params(m, num_samples=20):

    sweep_params = dict()
    sweep_params["system_recovery"] = LinearSample(
        m.fs.system_recovery, 0.25, 0.75, num_samples
    )

    return sweep_params


def build_outputs(m):

    outputs = dict()
    outputs["feed_conc"] = m.fs.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"]
    outputs["perm_conc"] = m.fs.product.properties[0].conc_mass_phase_comp[
        "Liq", "NaCl"
    ]
    outputs["brine_conc"] = m.fs.disposal.properties[0].conc_mass_phase_comp[
        "Liq", "NaCl"
    ]
    outputs["system_recovery"] = m.fs.system_recovery
    outputs["LCOW"] = m.fs.costing.LCOW
    outputs["SEC"] = m.fs.costing.SEC

    for n, stage in m.fs.stage.items():
        outputs[f"stage_{n}_inlet_flow"] = stage.feed.properties[0].flow_vol_phase[
            "Liq"
        ]
        outputs[f"stage_{n}_inlet_conc"] = stage.feed.properties[
            0
        ].conc_mass_phase_comp["Liq", "NaCl"]
        outputs[f"stage_{n}_permeate_flow"] = stage.product.properties[
            0
        ].flow_vol_phase["Liq"]
        outputs[f"stage_{n}_permeate_conc"] = stage.product.properties[
            0
        ].conc_mass_phase_comp["Liq", "NaCl"]
        outputs[f"stage_{n}_brine_flow"] = stage.disposal.properties[0].flow_vol_phase[
            "Liq"
        ]
        outputs[f"stage_{n}_brine_conc"] = stage.disposal.properties[
            0
        ].conc_mass_phase_comp["Liq", "NaCl"]
        outputs[f"stage_{n}_recovery"] = stage.recovery

        if stage.add_pump:
            outputs[f"stage_{n}_pressure"] = stage.pump.control_volume.properties_out[
                0
            ].pressure
        outputs[f"stage_{n}_area"] = stage.RO.area

    outputs["erd_power"] = m.fs.ERD.work_mechanical[0]

    return outputs


if __name__ == "__main__":

    salt_fracs = [5e-3, 10e-3, 20e-3, 35e-3, 50e-3, 60e-3, 70e-3]
    df = pd.DataFrame()
    for s in salt_fracs:
        x = int(s * 1000)
        save_file = f"tds_{x}.csv"

        results_array, results_dict = parameter_sweep(
            build_model=build_and_solve,
            build_model_kwargs={"n_stages": 2, "salt_mass_frac": s},
            build_sweep_params=build_sweep_params,
            build_outputs=build_outputs,
            csv_results_file_name=save_file,
        )
        df = pd.concat([df, pd.read_csv(save_file)], ignore_index=True)

    df.to_csv("salt_mass_frac_sweep.csv", index=False)
