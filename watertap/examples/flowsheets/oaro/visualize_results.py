import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.colors as colors
from matplotlib import cm
import matplotlib.patheffects as PathEffects

plt.rcParams.update({"font.size": 13})

output_dir = "param_sweep_output"

# filename_list = sorted(glob.glob(f'{output_dir}/*_stage/interpolated_results_OARO.csv'))
filename_list = sorted(glob.glob(f"{output_dir}/*/results_OARO.csv"))

viz_data = {}
viz_keys = [
    "LCOW",
    "SEC",
]

for k, fn in enumerate(filename_list):
    df = pd.read_csv(fn)

    if k == 0:
        xpts = np.unique(df["# Feed Concentration"].values)
        ypts = np.unique(df["Volumetric Recovery Rate"].values)
        nx = int(len(df) ** 0.5)

    for key in viz_keys:
        if key not in viz_data:
            viz_data[key] = []

        data = df[key].values
        data = data.reshape(nx, nx, order="F")
        viz_data[key].append(data)

#         plt.contourf(xpts, ypts, data)
#         plt.title(f'{key} - {fn}')
#         plt.show()

for key in viz_keys:
    viz_data[key] = np.array(viz_data[key])

min_id_data = np.array(viz_data["LCOW"])

min_id = np.zeros((nx, nx))
min_lcow = np.zeros((nx, nx))
min_sec = np.zeros((nx, nx))

for k in range(nx):
    for j in range(nx):
        try:
            a = np.nanargmin(viz_data["LCOW"][:, k, j])
            min_id[k, j] = a + 1
            min_lcow[k, j] = viz_data["LCOW"][a, k, j]
            min_sec[k, j] = viz_data["SEC"][a, k, j]
        except:
            min_id[k, j] = 9
            min_lcow[k, j] = np.nan
            min_sec[k, j] = np.nan

# Smooth out holes in the min_id data

made_replacement = True
iter_ct = 0
min_id_0 = np.copy(min_id)

while made_replacement == True and iter_ct < 10:
    made_replacement = False
    iter_ct += 1
    print(iter_ct)
    for k in range(1, nx - 1):
        for j in range(1, nx - 1):
            neigh = np.array(
                [min_id[k, j - 1], min_id[k, j + 1], min_id[k - 1, j], min_id[k + 1, j]]
            )

            u, counts = np.unique(neigh, return_counts=True)
            idx = np.argsort(counts)
            counts = counts[idx]
            u = u[idx]

            if counts[-1] >= 3 and u[-1] != min_id[k, j]:
                min_id[k, j] = u[-1]
                made_replacement = True

            if counts[-1] >= 2 and u[-1] < min_id[k, j]:
                min_id[k, j] = u[-1]
                made_replacement = True

hatch_color = "k"
plt.rcParams["hatch.color"] = hatch_color
mpl.rcParams["hatch.linewidth"] = 0.5
# plt.contourf(xpts, ypts, min_id)
# plt.colorbar()
# plt.show()
hatch_prog = ["/", "\\", "|", "//", "\\\\", "||", "///"]

fig, ax = plt.subplots(2, 1, figsize=(8, 12), dpi=300)

case_letter = {1: "A", 2: "B", 3: "C"}

cont_levels = [1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.1]
gnbu = cm.get_cmap("GnBu")
ylgnbu = cm.get_cmap("GnBu")
# ylgnbu_shift = ylgnbu(np.linspace(0, 1, 256)**(0.5))
# ylgnbu_shift = colors.ListedColormap(ylgnbu_shift)

for a, key in enumerate(viz_keys):
    if key == "LCOW":
        data = min_lcow
        title = r"Levelized Cost of Water (\$/m$^3$)"
        contf_levels = [0, 0.5, 1, 2, 4, 8, 16, 32]
    elif key == "SEC":
        data = min_sec
        title = r"Specific Energy Consumption (kWh/m$^3$)"
        contf_levels = [0, 2, 5, 10, 20, 40, 60, 90]

    #     data_copy = np.zeros(np.shape(data))
    #     for k in range(1, len(contf_levels)):
    #         lo_lim = contf_levels[k-1]
    #         hi_lim = contf_levels[k]
    #         print(f'range = {lo_lim} to {hi_lim}')
    #         if k < len(contf_levels)-1:
    #             mask = np.logical_and(lo_lim < data, data <= hi_lim)
    #         else:
    #             mask = lo_lim < data
    #         data_copy[mask] = data[mask]/hi_lim+k

    #     data_copy[np.isnan(data)] = np.nan

    #     print(np.nanmin(data_copy), np.nanmax(data_copy))

    cbar = ax[a].contourf(
        xpts,
        ypts,
        data,
        contf_levels,
        colors=ylgnbu(np.linspace(0, 1, len(contf_levels))),
        extend="max",
    )

    cbar = fig.colorbar(cbar, ax=ax[a], fraction=0.044)
    #     cbar.set_ticks(np.arange(1, len(contf_levels)))
    #     cbar.set_ticklabels(contf_levels)

    # ax[a].contour(xpts, ypts, min_id, cont_levels, colors='w', width=1)
    # ax[a].contour(xpts, ypts, min_id_0, cont_levels, colors='r')

    use_hatches = False

    if use_hatches == True:
        ax[a].contourf(
            xpts, ypts, min_id, cont_levels, colors="none", hatches=hatch_prog
        )

    else:
        label_rot = np.linspace(-78, -55, 8)
        label_pos = np.linspace(50, 180, 8)

        for k in range(8):
            row_num = 15
            mask = min_id[row_num, :] == (k + 1)
            if k >= 6:
                loc = np.nanmean(xpts[mask]) - 15
            elif k >= 4:
                loc = np.nanmean(xpts[mask]) - 13
            elif k >= 2:
                loc = np.nanmean(xpts[mask]) - 10
            elif k >= 1:
                loc = np.nanmean(xpts[mask]) - 5
            else:
                loc = np.nanmean(xpts[mask]) - 1

            ax[a].text(
                loc,
                ypts[row_num] + 0.05,
                "%d-Stage" % (k + 1),
                ha="center",
                va="center",
                rotation=label_rot[k],
            )

    ax[a].contour(xpts, ypts, min_id, cont_levels, colors=hatch_color, linewidths=0.5)

    ax[a].set_title(title)
    ax[a].set_xlabel("Feed Concentration (g/L)")
    ax[a].set_xticks(np.linspace(5, 245, 9))
    ax[a].set_xlim(5, 245)
    if a == 0:
        ax[a].set_ylabel("Volumetric Recovery Rate (%)")
    ax[a].set_yticks(np.linspace(0.3, 0.9, 7))
    ax[a].set_yticklabels([f"{perc:.0f}" for perc in 100 * np.linspace(0.3, 0.9, 7)])
    ax[a].set_box_aspect(1.0)

    specific_case_studies = [[35, 0.70], [70, 0.55], [125, 0.35]]

    for k, cs in enumerate(specific_case_studies):
        #         ax[a].plot([cs[0], cs[0]+35], [cs[1], cs[1]], 'w-', linewidth=16, alpha=0.4)
        ax[a].plot(
            cs[0],
            cs[1],
            "o",
            markersize=6,
            markeredgewidth=1,
            markerfacecolor="C3",
            markeredgecolor="w",
        )
        txt = ax[a].text(
            cs[0] + 5,
            cs[1],
            f"Case {case_letter[k + 1]}",
            ha="left",
            va="center",
            fontsize=12,
            color="C3",
            fontweight="bold",
        )
        txt.set_path_effects([PathEffects.withStroke(linewidth=3, foreground="w")])

# Make swatches of hatch pattern for the legend
patch_handles = []
patch_handles.append(
    mpatches.Patch(
        edgecolor=hatch_color, facecolor="none", hatch="", label="%d-Stage" % 1
    )
)

for hatch_num, hatch in enumerate(hatch_prog):
    patch_handles.append(
        mpatches.Patch(
            edgecolor=hatch_color,
            facecolor="none",
            hatch=hatch + hatch,
            label="%d-Stage" % (hatch_num + 2),
        )
    )

patch_handles_copy = patch_handles.copy()
patch_handles[0::2] = patch_handles_copy[0:4]
patch_handles[1::2] = patch_handles_copy[4:]

fig.legend(handles=patch_handles, loc=8, ncol=4, bbox_to_anchor=(0.5, -0.08))

fig.tight_layout()

plt.savefig("lcow_v1.jpg", bbox_inches="tight", facecolor="w")
plt.show()

# print(xpts, ypts)
print(xpts[1] - xpts[0], (ypts[1] - ypts[0]) * 100)
