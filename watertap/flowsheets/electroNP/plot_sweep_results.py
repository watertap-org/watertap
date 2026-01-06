#!/usr/bin/env python
"""
Plot parameter sweep results from multi_sweep.py
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys


def plot_2d_sweep(csv_file="sensitivity_2.csv"):
    """Plot results from 2D parameter sweep (two parameters varied) as heatmaps"""
    df = pd.read_csv(csv_file)

    # Identify sweep parameters
    sweep_params = []
    for col in df.columns:
        if df[col].nunique() > 1 and col not in ["solve_successful", "sweep_index"]:
            if (
                "LCOW" not in col
                and "Capital" not in col
                and "Electricity" not in col
                and "S_PO4" not in col
            ):
                sweep_params.append(col)

    param1, param2 = sweep_params

    # Get output columns
    outputs = [
        col
        for col in df.columns
        if col not in sweep_params + ["solve_successful", "sweep_index"]
    ]

    # Create heatmaps
    n_outputs = len(outputs)
    n_cols = min(2, n_outputs)
    n_rows = (n_outputs + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(8 * n_cols, 6 * n_rows))
    if n_outputs == 1:
        axes = [axes]
    else:
        axes = axes.flatten()

    for idx, output in enumerate(outputs):
        # Create heatmap from pivoted data
        pivot_data = df.pivot(index=param2, columns=param1, values=output)
        sns.heatmap(
            pivot_data,
            annot=True,
            fmt=".3g",
            cmap="viridis",
            ax=axes[idx],
            cbar_kws={"label": output},
        )
        axes[idx].set_title(f"{output}", fontsize=14, fontweight="bold")
        axes[idx].set_xlabel(param1, fontsize=12)
        axes[idx].set_ylabel(param2, fontsize=12)

    # Hide unused subplots
    for idx in range(n_outputs, len(axes)):
        axes[idx].set_visible(False)

    plt.tight_layout()
    plt.savefig(csv_file.replace(".csv", "_heatmap.png"), dpi=300, bbox_inches="tight")
    print(f"Saved heatmap to {csv_file.replace('.csv', '_heatmap.png')}")
    plt.show()


if __name__ == "__main__":
    # Check if a specific case number was provided
    case_num = 1
    if len(sys.argv) > 1:
        case_num = int(sys.argv[1])

    csv_file = f"sensitivity_{case_num}.csv"
    interpolated_file = f"interpolated_sensitivity_{case_num}.csv"

    # Try to use interpolated file if it exists
    import os

    if os.path.exists(interpolated_file):
        csv_file = interpolated_file
        print(f"Using interpolated results: {csv_file}")

    # Print summary
    df = pd.read_csv(csv_file)

    # Identify sweep parameters
    sweep_params = []
    for col in df.columns:
        if df[col].nunique() > 1 and col not in ["solve_successful", "sweep_index"]:
            if (
                "LCOW" not in col
                and "Capital" not in col
                and "Electricity" not in col
                and "S_PO4" not in col
            ):
                sweep_params.append(col)

    # Determine number of sweep parameters and plot accordingly
    df = pd.read_csv(csv_file)
    sweep_params = []
    for col in df.columns:
        if df[col].nunique() > 1 and col not in ["solve_successful", "sweep_index"]:
            if (
                "LCOW" not in col
                and "Capital" not in col
                and "Electricity" not in col
                and "S_PO4" not in col
            ):
                sweep_params.append(col)

    if len(sweep_params) == 2:
        print(f"\n2D sweep detected. Creating heatmaps...")
        plot_2d_sweep(csv_file)
    else:
        print(f"sweep parameters: {len(sweep_params)}")
