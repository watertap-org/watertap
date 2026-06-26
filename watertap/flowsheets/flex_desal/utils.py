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

"""
This module contains some utility functions
"""

from pyomo.environ import SolverFactory, value
from pyomo.common.dependencies import attempt_import
import numpy as np
import matplotlib.pyplot as plt

# pylint: disable = import-error
gurobipy, gurobipy_available = attempt_import("gurobipy", defer_import=False)
if gurobipy_available:
    from gurobipy import nlfunc


# NOTE: This function is needed until Pyomo's Gurobi writer can handle
# general nonlinear constraints. Once Pyomo supports this feature, this
# function can be removed from the codebase. Since this function is
# needed temporarily, tests are not added for this function.
def get_gurobi_solver_model(m, mip_gap=0.01, time_limit=3600, tee=True):
    """
    Returns a Pyomo SolverFactory object that is compatible with Gurobi.
    This function is needed only when the RO recovery is a variable.
    """
    if not gurobipy_available:
        pass
    else:
        solver = SolverFactory("gurobi_persistent")
        solver.options["MIPGap"] = mip_gap
        solver.options["TimeLimit"] = time_limit
        solver.options["OutputFlag"] = int(tee)

        if (
            not m.period[1, 1]
            .reverse_osmosis.ro_skid[1]
            .calculate_energy_intensity.active
        ):
            # If the nonlinear constraint is not active, then return the solver
            # object directly
            solver.set_instance(m)
            return solver

        if m.params.surrogate_type == "quadratic_surrogate":
            # Model is quadratic, so Pyomo's writer can handle it.
            solver.set_instance(m)
            return solver

        # Nonlinear constraint is present.
        # Step 1: Deactivate the nonlinear constraint
        for p in m.period:
            for skid in m.period[p].reverse_osmosis.set_ro_skids:
                m.period[p].reverse_osmosis.ro_skid[
                    skid
                ].calculate_energy_intensity.deactivate()

        # pylint: disable = protected-access
        # Step 2: Build the gurobipy model
        solver.set_instance(m)
        gm = solver._solver_model  # Gurobipy model
        pm_to_gm = solver._pyomo_var_to_solver_var_map

        # Step 3: Add the nonlinear constraint
        coeffs = m.period[1, 1].reverse_osmosis.ro_skid[1].coeffs
        for p in m.period:
            for skid in m.period[p].reverse_osmosis.set_ro_skids:
                ro_skid = m.period[p].reverse_osmosis.ro_skid[skid]
                recovery_var = pm_to_gm[ro_skid.recovery]
                # pylint: disable = possibly-used-before-assignment
                gm.addConstr(
                    (
                        pm_to_gm[ro_skid.energy_intensity]
                        == coeffs["a"] * nlfunc.exp(-coeffs["b"] * recovery_var)
                        + coeffs["c"] * recovery_var * recovery_var
                        + coeffs["d"]
                    ),
                    name=f"ro_energy_intensity_{p[0]}_{p[1]}_{skid}",
                )

    return solver


def fix_recovery(m, recovery):
    """Modifies the model for the fixed recovery case"""
    # Compute the energy intensity
    ro_skid = m.period[1, 1].reverse_osmosis.ro_skid[1]
    energy_intensity = m.params.ro.get_energy_intensity(recovery)

    for p in m.period:
        for skid in m.period[p].reverse_osmosis.set_ro_skids:
            ro_skid = m.period[p].reverse_osmosis.ro_skid[skid]
            ro_skid.recovery.fix(recovery)
            ro_skid.energy_intensity.fix(energy_intensity)
            ro_skid.calculate_energy_intensity.deactivate()


def wrd_fix_uf_recovery(m, uf_recovery):
    """Fixes the recovery of the UF pretreatment"""
    for p in m.period:
        for pump in m.period[p].pretreatment.set_uf_pumps:
            uf_pump = m.period[p].pretreatment.uf_pumps[pump]
            uf_pump.recovery.fix(uf_recovery)


def wrd_fix_ro_recovery(m, ro_recovery):
    """Modifies the model for the fixed recovery case"""
    # For WRD model, energy intensity depends on flowrate and recovery.
    # If there is a case where recovery needs to be fixed, use this function
    for p in m.period:
        for skid in m.period[p].reverse_osmosis.set_ro_skids:
            ro_skid = m.period[p].reverse_osmosis.ro_skid[skid]
            ro_skid.recovery.fix(ro_recovery)


def update_recovery_bounds(m, lb, ub):
    """Updates the bounds on the recovery variable"""
    ro_skid = m.period[1, 1].reverse_osmosis.ro_skid[1]
    ei_lb, ei_ub = m.params.ro.get_energy_intensity_bounds(lb, ub)

    for p in m.period:
        for skid in m.period[p].reverse_osmosis.set_ro_skids:
            ro_skid = m.period[p].reverse_osmosis.ro_skid[skid]
            ro_skid.recovery.setlb(lb)
            ro_skid.recovery.setub(ub)
            ro_skid.energy_intensity.setlb(ei_lb)
            ro_skid.energy_intensity.setub(ei_ub)


def get_baseline_model(m):
    """Returns a baseline model from the given model"""
    bm = m.clone()

    # Ensure that the pretreatment unit is always on, and no leakage from it
    bm.fix_operation_var("pretreatment.op_mode", 1)
    bm.fix_operation_var("pretreatment.recovery", 1)

    # Ensure that the first three skids all always on
    for skid in [1, 2, 3]:
        bm.fix_operation_var(f"reverse_osmosis.ro_skid[{skid}].op_mode", 1)

    if hasattr(m.period[1, 1], "battery"):
        # If the battery model exists, ensure that there is no charging
        # and discharging
        bm.fix_operation_var("battery.power_charge", 0)
        bm.fix_operation_var("battery.power_discharge", 0)

    return bm


def wrd_plot_function(m, n_time_points, output_stem=None, peak_hours=None):
    time = np.linspace(0, n_time_points - 1, n_time_points)
    fig = plt.figure(figsize=(12, 12))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 1])
    ax_energy = fig.add_subplot(gs[0])
    ax_trains = fig.add_subplot(gs[1], sharex=ax_energy)
    ax_energy.set_facecolor("#f5f5f5")
    ax_trains.set_facecolor("#f5f5f5")

    if peak_hours is not None:
        peak_legend_added = False
        for i, is_peak in enumerate(peak_hours):
            if is_peak:
                # Shade full hourly intervals where variable demand charges apply.
                span_label = "Peak Hours" if not peak_legend_added else None
                ax_energy.axvspan(
                    i,
                    i + 1,
                    color="grey",
                    alpha=0.2,
                    linewidth=0,
                    zorder=-1,
                    label=span_label,
                )
                ax_trains.axvspan(
                    i,
                    i + 1,
                    color="grey",
                    alpha=0.2,
                    linewidth=0,
                    zorder=-1,
                    label=span_label,
                )
                peak_legend_added = True

    # First subplot: Stacked energy consumption by major equipment
    total_energy = [
        value(v[None]) for v in m.period[:, :].net_power_consumption.extract_values()
    ]

    ro1_energy = [
        v[None]
        for v in m.period[:, :]
        .reverse_osmosis.ro_skid[1]
        .power_consumption.extract_values()
    ]
    ro2_energy = [
        v[None]
        for v in m.period[:, :]
        .reverse_osmosis.ro_skid[2]
        .power_consumption.extract_values()
    ]
    ro3_energy = [
        v[None]
        for v in m.period[:, :]
        .reverse_osmosis.ro_skid[3]
        .power_consumption.extract_values()
    ]
    ro4_energy = [
        v[None]
        for v in m.period[:, :]
        .reverse_osmosis.ro_skid[4]
        .power_consumption.extract_values()
    ]

    uf1_energy = [
        v[None]
        for v in m.period[:, :]
        .pretreatment.uf_pumps[1]
        .power_consumption.extract_values()
    ]
    uf2_energy = [
        v[None]
        for v in m.period[:, :]
        .pretreatment.uf_pumps[2]
        .power_consumption.extract_values()
    ]
    uf3_energy = [
        v[None]
        for v in m.period[:, :]
        .pretreatment.uf_pumps[3]
        .power_consumption.extract_values()
    ]

    other_energy = np.array(total_energy) - (
        np.array(ro1_energy)
        + np.array(ro2_energy)
        + np.array(ro3_energy)
        + np.array(ro4_energy)
        + np.array(uf1_energy)
        + np.array(uf2_energy)
        + np.array(uf3_energy)
    )
    # Clip tiny negatives from solver tolerances so stackplot remains well-defined.
    other_energy = np.maximum(other_energy, 0.0)

    ax_energy.stackplot(
        time + 0.5,
        ro1_energy,
        ro2_energy,
        ro3_energy,
        ro4_energy,
        uf1_energy,
        uf2_energy,
        uf3_energy,
        other_energy,
        labels=[
            "RO Train 1",
            "RO Train 2",
            "RO Train 3",
            "RO Train 4",
            "UF Pump 1",
            "UF Pump 2",
            "UF Pump 3",
            "Post-Treatment",
        ],
        alpha=0.5,
    )

    ax_energy.plot(
        time + 0.5,
        total_energy,
        label="Total Energy Consumption",
        color="black",
        linestyle="--",
        linewidth=2,
    )
    ax_energy.set_ylim(0, 2500)
    ax_energy.set_ylabel("Energy Consumption (kWh)", fontsize=12)
    ax_energy.set_title(
        "Energy Consumption and Electricity Price", fontsize=14, fontweight="bold"
    )
    ax_energy.grid(False)

    ax_price = ax_energy.twinx()
    elec_price = m._config.lmp_data
    ax_price.plot(
        time + 0.5,
        elec_price,
        label="Electricity Cost ($/kWh)",
        color="orange",
        linestyle="--",
        linewidth=2,
    )
    ax_price.set_ylabel("Electricity Cost ($/kWh)", fontsize=12)
    ax_price.set_ylim(0, max(elec_price) + 0.03)

    handle1, label1 = ax_energy.get_legend_handles_labels()
    handle2, label2 = ax_price.get_legend_handles_labels()
    handles = handle2 + handle1
    labels = label2 + label1
    leg1 = ax_price.legend(
        handles, labels, loc="lower left", framealpha=1.0, ncol=2, fontsize=10
    )
    leg1.set_zorder(1000)
    leg1.get_frame().set_facecolor("white")
    ax_energy.xaxis.set_major_locator(plt.MaxNLocator(24))

    # Second subplot: Water production and RO train flow rates
    prod = [
        v[None] for v in m.period[:, :].posttreatment.product_flowrate.extract_values()
    ]
    ax_trains.plot(
        time + 0.5,
        prod,
        label="Water Production (m$^3$/h)",
        color="black",
        linestyle="--",
        linewidth=2,
        alpha=0.75,
    )
    ax_trains.set_ylim(0, 2500)
    ax_trains.axhline(
        y=602 * 4,
        color="blue",
        linestyle="--",
        linewidth=2,
        alpha=0.75,
        label="Nominal Plant Capacity (m$^3$/h)",
        zorder=0,
    )
    ax_trains.set_ylabel("Water Production (m$^3$/h)", fontsize=12)
    ax_trains.set_xlabel("Hours", fontsize=12)
    ax_trains.set_title(
        "Water Production & RO Train Flow Rates", fontsize=14, fontweight="bold"
    )
    ax_trains.xaxis.set_major_locator(plt.MaxNLocator(24))
    ax_trains.grid(False)

    # Extract RO train flow rates (m3/hr) for stacked plotting
    train_1_flows = [
        v[None]
        for v in m.period[:, :]
        .reverse_osmosis.ro_skid[1]
        .product_flowrate.extract_values()
    ]
    train_2_flows = [
        v[None]
        for v in m.period[:, :]
        .reverse_osmosis.ro_skid[2]
        .product_flowrate.extract_values()
    ]
    train_3_flows = [
        v[None]
        for v in m.period[:, :]
        .reverse_osmosis.ro_skid[3]
        .product_flowrate.extract_values()
    ]
    train_4_flows = [
        v[None]
        for v in m.period[:, :]
        .reverse_osmosis.ro_skid[4]
        .product_flowrate.extract_values()
    ]

    ax_trains.stackplot(
        time + 0.5,
        train_1_flows,
        train_2_flows,
        train_3_flows,
        train_4_flows,
        labels=["RO Train 1", "RO Train 2", "RO Train 3", "RO Train 4"],
        alpha=0.5,
    )

    handle_t, label_t = ax_trains.get_legend_handles_labels()
    leg3 = ax_trains.legend(
        handle_t,
        label_t,
        loc="lower left",
        fontsize=11,
        framealpha=1.0,
        ncol=2,
    )
    leg3.get_frame().set_facecolor("white")

    # Set consistent x-axis limits and formatting
    for a in (ax_energy, ax_trains):
        a.set_xlim(0, n_time_points)
        a.xaxis.set_major_locator(plt.MaxNLocator(24))

    # Tick labels for all axes
    for a in (ax_energy, ax_price, ax_trains):
        a.tick_params(axis="both", labelsize=11)

    fig.tight_layout()
    if output_stem is not None:
        fig.savefig(f"{output_stem}.png", dpi=600)
    plt.show()
