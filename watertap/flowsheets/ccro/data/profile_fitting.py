from scipy.optimize import rosen, differential_evolution
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math


def get_data(file_path):
    """Load data from a CSV file and return as a pandas DataFrame."""
    data = pd.read_csv(file_path)
    t = data["t"].to_numpy()
    f_t = data["F_t"].to_numpy()
    return t, f_t


def ft_profile(params, t):
    """Calculate the profile function F(t) based on the given parameters."""
    t_m, N = params
    tau = t_m / int(N)
    N = int(N)
    sum_terms = [np.zeros(len(t))]
    for n in range(N - 1):
        term = ((t / tau) ** n) / math.factorial(n)
        sum_terms.append(term)
    # print(np.sum(sum_terms, axis=0))
    if N == 1:
        F_t = 1 - np.exp(-t / tau)
    else:
        F_t = 1 - np.exp(-t / tau) * np.sum(sum_terms, axis=0)
    # print(F_t)
    return F_t


def objective_function(params, t, f_t):
    """Calculate the sum of squared errors between the model and data."""
    f_t_model = ft_profile(params, t)

    print(f_t_model)
    return np.sum((f_t - f_t_model) ** 2)


if __name__ == "__main__":

    fig, ax = plt.subplots()
    fig.set_dpi(300)
    # Set figure size
    fig.set_size_inches(3.25, 3.25, forward=True)
    for data in ["swc_1_new.csv", "swc_1_worn.csv", "BW2540.csv"]:
        t, f_t = get_data(data)
        bounds = [(1e-5, 500), (1, 50)]  # Bounds for tau and N
        result = differential_evolution(
            objective_function, bounds, workers=8, args=(t, f_t), tol=0.0001
        )
        print(result)
        tau_opt, N_opt = result.x
        # tau_opt = 5
        # N_opt = 2
        print(f"Optimized parameters: tau = {tau_opt}, N = {int(N_opt)}")
        f_t_opt = ft_profile((tau_opt, N_opt), t)

        ax.plot(t, f_t, "o", label="Data", markersize=3)
        ax.plot(
            t,
            f_t_opt,
            "-",
            linewidth=1,
            label=f"Data: {data}, $\\t_m$={tau_opt:.2f}, N={int(N_opt)}",
        )
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("F(t)")
    ax.legend()
    # plt.tight_layout()
    plt.show()
