import numpy as np
import pandas as pd
from pyomo.environ import value
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap, BoundaryNorm

def plot_contour(
    df,
    x="fs.flow_mgd",
    y="tds",
    z="fs.costing.LCOW",
    min_x=None,
    min_y=None,
    max_x=None,
    max_y=None,
    z_adj=1,
    approach="pivot",
    cmap="turbo",
    interp_method="cubic",
    grid_len=100,
    levels=10,
    levelsf=None,
    set_dict=dict(),
    cb_title="",
    contour_label_fmt="  %#.1f \$/m$^3$  ",
    add_contour_labels=True,
    fig=None,
    ax=None,
    figsize=(6, 4),
):
    """
    Generic function to create contour plot from
    three columns in DataFrame
    """
    if (fig, ax) == (None, None):
        fig, ax = plt.subplots(figsize=figsize)
        fig.set_size_inches(5, 5, forward=True)

    if min_x is not None:
        df = df[df[x] >= min_x].copy()
    if min_y is not None:
        df = df[df[y] >= min_y].copy()

    if max_x is not None:
        df = df[df[x] <= min_x].copy()
    if max_y is not None:
        df = df[df[y] <= min_y].copy()

    if levelsf is None:
        levelsf = levels

    if approach == "pivot":
        pivot = df.pivot(index=y, columns=x, values=z * z_adj)

        x = pivot.columns.values
        y = pivot.index.values
        z = pivot.values

        contourf = ax.contourf(x, y, z, levels=levels, cmap=cmap)
        # cb = plt.colorbar(contourf, label=cb_title)

    if approach == "interpolate":

        x = df[x]
        y = df[y]
        print(z)
        z = df[z] * z_adj

        xi = np.linspace(min(x), max(x), grid_len)
        yi = np.linspace(min(y), max(y), grid_len)
        xi, yi = np.meshgrid(xi, yi)

        zi = griddata((x, y), z, (xi, yi), method=interp_method)

        contourf = ax.contourf(xi, yi, zi, levels=levels, cmap=cmap)
        # cb = plt.colorbar(contourf, label=cb_title)

    ax.set(**set_dict)

    if add_contour_labels:
        contour = ax.contour(
            contourf,
            levelsf,
            colors="k",
            linestyles="dashed",
        )
        ax.clabel(contour, colors="black", fmt=contour_label_fmt, fontsize=9)

    ax.set_xticklabels([])
    ax.set_yticklabels([])
    # ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x*100:.0f}%"))
    # ax.tick_params(axis="both", labelsize=16)

    # plt.tight_layout()

    return fig, ax

def plot_min_lcow(df1, df2, df3, fig, axs):
    ax = axs["combo"]

    Z1_df = (
        df1.pivot(index="steam_cost", columns="electricity_cost", values="LCOW")
        .fillna(method="ffill")
        .fillna(method="bfill")
    )
    Z2_df = (
        df2.pivot(index="steam_cost", columns="electricity_cost", values="LCOW")
        .fillna(method="ffill")
        .fillna(method="bfill")
    )
    Z3_df = (
        df3.pivot(index="steam_cost", columns="electricity_cost", values="LCOW")
        .fillna(method="ffill")
        .fillna(method="bfill")
    )

    y_vals = np.sort(Z1_df.index.values)
    x_vals = np.sort(Z1_df.columns.values)
    Z1 = Z1_df.reindex(index=y_vals, columns=x_vals).to_numpy()
    Z2 = Z2_df.reindex(index=y_vals, columns=x_vals).to_numpy()
    Z3 = Z3_df.reindex(index=y_vals, columns=x_vals).to_numpy()

    stack = np.stack([Z1, Z2, Z3], axis=0)
    Z_min = stack.min(axis=0)
    Dom = stack.argmin(axis=0)

    X, Y = np.meshgrid(x_vals, y_vals)

    cmap = ListedColormap(["tab:red", "lightgreen", "royalblue"])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)
    ax.pcolormesh(X, Y, Dom, cmap=cmap, norm=norm, shading="nearest")
    cf = ax.contour(X, Y, Z_min, levels=10, colors="k", linestyles="dashed")
    ax.clabel(cf, colors="black", fmt="%#.1f \$/m$^3$", fontsize=12)

    colors = ["tab:red", "lightgreen", "royalblue"]
    labels = ["FC", "MVC", "TVC"]
    patches = list(mpatches.Patch(color=c, label=l) for c, l in zip(colors, labels))

    ax.legend(
        handles=patches,
        loc="lower center",
        bbox_to_anchor=(0, 0.9, 1, 1.0),
        frameon=True,
        ncols=3,
        fontsize=16,
    )

    fig.supxlabel("Electricity Cost (¢/kWh)", fontsize=16)

    fig.supylabel("Steam Cost (¢/kg)", fontsize=16)
    # ax.set_title("Crystallizer Technology Screening\nMinimum LCOW", fontsize=16)
    fig.suptitle("Crystallizer Technology Screening\nMinimum LCOW", fontsize=16)
    ax.tick_params(axis="both", labelsize=12)


def kpis_mass_energy(m):
    # Feed (kg/s): Liq H2O + Liq NaCl
    feed_kg_s = value(m.fs.feed.flow_mass_phase_comp[0, "Liq", "H2O"]) + value(
        m.fs.feed.flow_mass_phase_comp[0, "Liq", "NaCl"]
    )

    # Distillate (kg/s)
    if hasattr(m.fs, "distillate"):  # Steam-Driven & MVC
        d = m.fs.distillate.properties[0]
        try:
            dist_kg_s = value(d.flow_mass_phase_comp["Liq", "H2O"])
        except Exception:
            dist_kg_s = value(d.flow_mass_phase_comp["Vap", "H2O"])
    else:  # TVC: sum heater+condenser product

        def _water_mass(sb):
            try:
                return value(sb.flow_mass_phase_comp["Liq", "H2O"])
            except Exception:
                return value(sb.flow_mass_phase_comp["Vap", "H2O"])

        dh = m.fs.distillate_heater.properties[0]
        dc = m.fs.distillate_condenser.properties[0]
        dist_kg_s = _water_mass(dh) + _water_mass(dc)

    # Solids NaCl (kg/s)
    solids_kg_s = value(m.fs.crystallizer.solids.flow_mass_phase_comp[0, "Sol", "NaCl"])

    # Steam (kg/s) – live steam to heater if present (MVC may not use it)
    steam_kg_s = (
        value(m.fs.heater.hot_side_inlet.flow_mass_phase_comp[0, "Vap", "H2O"])
        if hasattr(m.fs.heater, "hot_side_inlet")
        else float("nan")
    )

    # Heat duty (MW) & heater area (m²)
    heat_MW = value(m.fs.heater.hot.heat[0]) / 1e6
    area_m2 = value(m.fs.heater.area) if hasattr(m.fs.heater, "area") else float("nan")

    # LCOW ($/m³)
    lcow = value(m.fs.costing.LCOW)
    capex = value(m.fs.costing.total_capital_cost) * 1e-6
    opex = value(m.fs.costing.total_operating_cost) * 1e-6
    agg_flow_elec = value(m.fs.costing.aggregate_flow_electricity) * 1e-6
    SEC = value(m.fs.costing.specific_energy_consumption)

    return pd.DataFrame(
        [
            {
                "Feed (kg/s)": feed_kg_s,
                "Distillate (kg/s)": dist_kg_s,
                "Solids NaCl (kg/s)": solids_kg_s,
                "Steam (kg/s)": steam_kg_s,
                "Heat duty (MW)": heat_MW,
                "Heater area (m²)": area_m2,
                "LCOW ($/m³)": lcow,
                "SEC (kWh/m³)": SEC,
                "CAPEX ($M)": capex,
                "OPEX ($M/yr)": opex,
            }
        ]
    )

