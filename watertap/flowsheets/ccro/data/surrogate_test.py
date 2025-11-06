from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    Expression,
    Block,
    Param,
    value,
    Var,
    NonNegativeReals,
    assert_optimal_termination,
    Objective,
    units as pyunits,
)
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from idaes.core.surrogate.pysmo_surrogate import PysmoRBFTrainer, PysmoSurrogate

from idaes.core.surrogate.surrogate_block import SurrogateBlock

import numpy as np
import matplotlib.pyplot as plt
from watertap.core.solvers import get_solver


solver = get_solver()


def main():
    f = "flushing_surrogate_multiple_tau_n_2.json"
    m = ConcreteModel()
    m.flushing_time = Var()
    m.mean_residence_time = Var()
    m.flushing_efficiency = Var()

    m.surrogate_blk = SurrogateBlock(concrete=True)
    m.surrogate = PysmoSurrogate.load_from_file(f)
    m.surrogate_blk.build_model(
        m.surrogate,
        input_vars=[m.flushing_time, m.mean_residence_time],
        output_vars=[m.flushing_efficiency],
    )
    m.surrogate_blk.pysmo_constraint.display()  # display()
    # assert False
    m.surrogate_blk.pysmo_constraint["F_t"].pprint()
    assert False
    minx, maxx = m.flushing_time.bounds
    miny, maxy = m.mean_residence_time.bounds

    x_vals = np.linspace(1, maxx, 10)
    y_vals = np.linspace(1, maxy, 10)
    X, Y = np.meshgrid(x_vals, y_vals)
    Z = np.zeros_like(X)
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            m.flushing_time.fix(X[i, j])
            m.mean_residence_time.fix(Y[i, j])
            calculate_variable_from_constraint(
                m.flushing_efficiency, m.surrogate_blk.pysmo_constraint["F_t"]
            )
            m.surrogate_blk.display()
            # results = solver.solve(m, tee=True)
            # m.surrogate_blk.display()
            print(m.flushing_time.value, m.mean_residence_time.value)
            # assert False
            print(f"flushing efficiency: {value(m.flushing_efficiency)}")
        #  Z[i, j] = value(m.flushing_efficiency)

    # fig, ax = plt.subplots(figsize=(8, 6))
    # ax.contourf(X, Y, Z, levels=20, cmap="viridis")
    # ax.set_xlabel("Flushing Time")
    # ax.set_ylabel("Mean Residence Time")
    # ax.set_title("Flushing Efficiency Contour Plot")
    # plt.colorbar(
    #     ax.contourf(X, Y, Z, levels=20, cmap="viridis"),
    #     ax=ax,
    #     label="Flushing Efficiency",
    # )
    # plt.show()


if __name__ == "__main__":
    main()
