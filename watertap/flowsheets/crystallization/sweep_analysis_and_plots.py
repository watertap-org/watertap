import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap
from scipy.ndimage import label, center_of_mass
from sweep_electricity import run_electricity_price_sweep as run_electricity_sweep
from sweep_steam import run_steam_cost_sweep as run_steam_hr_sweep
from sweep_steam_no_heat_recovery import run_steam_cost_sweep as run_steam_nhr_sweep
from matplotlib.ticker import MaxNLocator

def main():
    n_points = 500
    run_electricity_sweep(nx=n_points, output_filename="electricity_price_sweep_elec.csv")
    df_elec = pd.read_csv("electricity_price_sweep_elec.csv")
    LCOW_elec = df_elec["LCOW"].values

    run_steam_hr_sweep(nx=n_points, output_filename="steam_price_sweep_hr.csv")
    df_steam_hr = pd.read_csv("steam_price_sweep_hr.csv")
    LCOW_steam_hr = df_steam_hr["LCOW"].values

    run_steam_nhr_sweep(nx=n_points, output_filename="steam_price_sweep_nhr.csv")
    df_steam_nhr = pd.read_csv("steam_price_sweep_nhr.csv")
    LCOW_steam_nhr = df_steam_nhr["LCOW"].values

    elec_vals = np.linspace(0.0, 0.5, n_points)
    steam_vals = np.linspace(0.0, 0.008, n_points)
    elec_vals_scaled = elec_vals * 100
    steam_vals_scaled = steam_vals * 1000

    nx_val = len(steam_vals)
    ny_val = len(elec_vals)
    min_LCOW = np.empty((nx_val, ny_val))
    cheapest_sheet = np.empty((nx_val, ny_val), dtype=object)

    for i in range(nx_val):
        for j in range(ny_val):
            cost_elec = LCOW_elec[j]
            cost_steam_hr = LCOW_steam_hr[i]
            cost_steam_nhr = LCOW_steam_nhr[i]
            costs = [cost_elec, cost_steam_hr, cost_steam_nhr]
            min_cost = min(costs)
            min_LCOW[i, j] = min_cost
            if costs.index(min_cost) == 0:
                cheapest_sheet[i, j] = "FC+MVC"
            elif costs.index(min_cost) == 1:
                cheapest_sheet[i, j] = "FC+TVC"
            else:
                cheapest_sheet[i, j] = "FC"

    x_step = elec_vals_scaled[1] - elec_vals_scaled[0]
    y_step = steam_vals_scaled[1] - steam_vals_scaled[0]
    x_min = elec_vals_scaled[0] - x_step / 2
    x_max = elec_vals_scaled[-1] + x_step / 2
    y_min = steam_vals_scaled[0] - y_step / 2
    y_max = steam_vals_scaled[-1] + y_step / 2

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(min_LCOW, origin="lower", extent=[x_min, x_max, y_min, y_max], aspect="auto", cmap="viridis")
    plt.colorbar(im, label="Minimum LCOW ($/m³)")
    ax.set_xlabel("Electricity Cost (cents/kWh)")
    ax.set_ylabel("Steam Cost ($/metric ton)")
    ax.set_title("Cheapest Flowsheet by LCOW over Steam and Electricity Costs")

    label_to_int = {"FC+MVC": 0, "FC+TVC": 1, "FC": 2}
    int_grid = np.vectorize(label_to_int.get)(cheapest_sheet)
    X, Y = np.meshgrid(elec_vals_scaled, steam_vals_scaled)
    ax.contour(X, Y, int_grid, levels=[0.5, 1.5], colors=['white', 'red'], linewidths=2)

    for sheet, int_val in label_to_int.items():
        mask = int_grid == int_val
        labeled_array, num_features = label(mask)
        for region in range(1, num_features + 1):
            region_mask = (labeled_array == region)
            com = center_of_mass(region_mask)
            center_x = np.interp(com[1], np.arange(ny_val), elec_vals_scaled)
            center_y = np.interp(com[0], np.arange(nx_val), steam_vals_scaled)
            ax.text(center_x, center_y, sheet, ha="center", va="center", color="black", fontsize=12, bbox=dict(facecolor="white", alpha=0.5))

    legend_elements = [
        Line2D([0], [0], color='white', lw=2, label='FC+MVC ↔ FC+TVC'),
        Line2D([0], [0], color='red', lw=2, label='FC+TVC ↔ FC')
    ]
    ax.legend(handles=legend_elements, loc='upper left')
    ax.set_xticks(elec_vals_scaled)
    ax.set_yticks(steam_vals_scaled)
    ax.xaxis.set_major_locator(MaxNLocator(10))
    ax.yaxis.set_major_locator(MaxNLocator(10))
    plt.tight_layout()
        
    plt.savefig("steam-electricity-heat-map.png", dpi=200)

    # --- Build & save LCOW vs Steam (FC vs TVC) from the steam CSVs ---
    
    plt.figure(figsize=(7,4))
    plt.plot(steam_vals_scaled, LCOW_steam_nhr, label="FC (no HR)", lw=2, color="#d62728")
    plt.plot(steam_vals_scaled, LCOW_steam_hr,  label="FC+TVC",     lw=2, color="#1f77b4")
    plt.xlabel("Steam Cost ($/metric ton)")
    plt.ylabel("LCOW ($/m³)")
    plt.title("LCOW vs Steam Cost")
    plt.legend(); plt.tight_layout()
    plt.savefig("steam_cost_vs_lcow.png", dpi=200)

    # --- Build & save Least-Cost Technology Map (categorical) ---
   
    cmap = ListedColormap(["#2ca02c", "#1f77b4", "#d62728"])  # MVC, TVC, FC
    fig2, ax2 = plt.subplots(figsize=(8,6))
    int_grid = np.vectorize({"FC+MVC":0, "FC+TVC":1, "FC":2}.get)(cheapest_sheet)
    ax2.imshow(int_grid, origin="lower",
               extent=[x_min, x_max, y_min, y_max],
               aspect="auto", cmap=cmap)
    ax2.set_xlabel("Electricity Cost (cents/kWh)")
    ax2.set_ylabel("Steam Cost ($/metric ton)")
    ax2.set_title("Least-Cost Flowsheet by Energy Prices")
    
    legend_elems = [Line2D([0],[0], color=c, lw=6, label=lab)
                    for lab,c in zip(["FC+MVC","FC+TVC","FC"], cmap.colors)]
    ax2.legend(handles=legend_elems, loc="upper left")
    plt.tight_layout()
    plt.savefig("least-cost-map.png", dpi=200)


    plt.show()

if __name__ == "__main__":
    main()
