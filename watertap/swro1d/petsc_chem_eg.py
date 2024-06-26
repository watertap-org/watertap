import numpy as np
import matplotlib.pyplot as plt

import pyomo.dae as pyodae
import pyomo.environ as pyo
import idaes.core.solvers.petsc as petsc  # petsc utilities module
# from idaes.core.solvers.features import dae  # DAE example/test problem

def dae(nfe=1):
    """This provides a DAE model for solver testing.

    The problem and expected result are from the problem given here:
    https://archimede.dm.uniba.it/~testset/report/chemakzo.pdf.

    Args:
        None

    Returns:
        (tuple): Pyomo ConcreteModel, correct solved value for y[1] to y[5] and y6
    """
    model = pyo.ConcreteModel(name="chemakzo")

    # Set problem parameter values
    model.k = pyo.Param([1, 2, 3, 4], initialize={1: 18.7, 2: 0.58, 3: 0.09, 4: 0.42})
    model.Ke = pyo.Param(initialize=34.4)
    model.klA = pyo.Param(initialize=3.3)
    model.Ks = pyo.Param(initialize=115.83)
    model.pCO2 = pyo.Param(initialize=0.9)
    model.H = pyo.Param(initialize=737)

    # Problem variables ydot = dy/dt,
    #    (dy6/dt is not explicitly in the equations, so only 5 ydots i.e.
    #    y6 is an algebraic variable and y1 to y5 are differential variables)
    model.t = pyodae.ContinuousSet(bounds=(0, 180))
    model.y = pyo.Var(model.t, [1, 2, 3, 4, 5], initialize=1.0)  #
    model.y6 = pyo.Var(model.t, initialize=1.0)  #
    model.ydot = pyodae.DerivativeVar(model.y, wrt=model.t)  # dy/dt
    model.r = pyo.Var(model.t, [1, 2, 3, 4, 5], initialize=1.0)
    model.Fin = pyo.Var(model.t, initialize=1.0)

    # Equations
    @model.Constraint(model.t)
    def eq_ydot1(b, t):
        return b.ydot[t, 1] == -2.0 * b.r[t, 1] + b.r[t, 2] - b.r[t, 3] - b.r[t, 4]

    @model.Constraint(model.t)
    def eq_ydot2(b, t):
        return b.ydot[t, 2] == -0.5 * b.r[t, 1] - b.r[t, 4] - 0.5 * b.r[t, 5] + b.Fin[t]

    @model.Constraint(model.t)
    def eq_ydot3(b, t):
        return b.ydot[t, 3] == b.r[t, 1] - b.r[t, 2] + b.r[t, 3]

    @model.Constraint(model.t)
    def eq_ydot4(b, t):
        return b.ydot[t, 4] == -b.r[t, 2] + b.r[t, 3] - 2.0 * b.r[t, 4]

    @model.Constraint(model.t)
    def eq_ydot5(b, t):
        return b.ydot[t, 5] == b.r[t, 2] - b.r[t, 3] + b.r[t, 5]

    @model.Constraint(model.t)
    def eq_y6(b, t):
        return 0 == b.Ks * b.y[t, 1] * b.y[t, 4] - b.y6[t]

    @model.Constraint(model.t)
    def eq_r1(b, t):
        return b.r[t, 1] == b.k[1] * b.y[t, 1] ** 4 * b.y[t, 2] ** 0.5

    @model.Constraint(model.t)
    def eq_r2(b, t):
        return b.r[t, 2] == b.k[2] * b.y[t, 3] * b.y[t, 4]

    @model.Constraint(model.t)
    def eq_r3(b, t):
        return b.r[t, 3] == b.k[2] / b.Ke * b.y[t, 1] * b.y[t, 5]

    @model.Constraint(model.t)
    def eq_r4(b, t):
        return b.r[t, 4] == b.k[3] * b.y[t, 1] * b.y[t, 4] ** 2

    @model.Constraint(model.t)
    def eq_r5(b, t):
        return b.r[t, 5] == b.k[4] * b.y6[t] ** 2 * b.y[t, 2] ** 0.5

    @model.Constraint(model.t)
    def eq_Fin(b, t):
        return b.Fin[t] == b.klA * (b.pCO2 / b.H - b.y[t, 2])

    # Set initial conditions and solve initial from the values of differential
    # variables.
    y0 = {1: 0.444, 2: 0.00123, 3: 0.0, 4: 0.007, 5: 0.0}  # initial differential vars
    for i, v in y0.items():
        model.y[0, i].fix(v)

    discretizer = pyo.TransformationFactory("dae.finite_difference")
    discretizer.apply_to(model, nfe=nfe, scheme="BACKWARD")

    return (
        model,
        0.1150794920661702,
        0.1203831471567715e-2,
        0.1611562887407974,
        0.3656156421249283e-3,
        0.1708010885264404e-1,
        0.4873531310307455e-2,
    )

# Get the model and known solution for y variables at t=180 minutes.
m, y1, y2, y3, y4, y5, y6 = dae(nfe=10)

# See the initial conditions:
print("at t = 0:")
print(f"    y1 = {pyo.value(m.y[0, 1])}")
print(f"    y2 = {pyo.value(m.y[0, 2])}")
print(f"    y3 = {pyo.value(m.y[0, 3])}")
print(f"    y4 = {pyo.value(m.y[0, 4])}")
print(f"    y5 = {pyo.value(m.y[0, 5])}")

# The command below will solve the problem.  In this case, we want to read the saved
# trajectory for each time element in the Pyomo.DAE problem (in this case there is
# only 1) so we will need to provide solver options to save the trajectory to the PETSc
# solver, a file name stub for variable information files, and a file stub for saving
# the trajectory information.  The options shown below will delete the trajectory
# information written by PETSc and resave it as json.  This allows us to cleanly read
# the trajectory data for multiple time elements.

result = petsc.petsc_dae_by_time_element(
    m,
    time=m.t,
    between=[m.t.first(), m.t.last()],
    ts_options={
        "--ts_type": "cn",  # Crank–Nicolson
        "--ts_adapt_type": "basic",
        "--ts_dt": 0.01,
        "--ts_save_trajectory": 1,
    },
)
tj = result.trajectory
res = result.results

print(abs(y1 - pyo.value(m.y[180, 1])) / y1)
print(abs(y2 - pyo.value(m.y[180, 2])) / y2)
print(abs(y3 - pyo.value(m.y[180, 3])) / y3)
print(abs(y4 - pyo.value(m.y[180, 4])) / y4)
print(abs(y5 - pyo.value(m.y[180, 5])) / y5)
print(abs(y6 - pyo.value(m.y6[180])) / y6)

# Verify results
# assert abs(y1 - pyo.value(m.y[180, 1])) / y1 < 1e-3
# assert abs(y2 - pyo.value(m.y[180, 2])) / y2 < 1e-3
# assert abs(y3 - pyo.value(m.y[180, 3])) / y3 < 1e-3
# assert abs(y4 - pyo.value(m.y[180, 4])) / y4 < 1e-3
# assert abs(y5 - pyo.value(m.y[180, 5])) / y5 < 1e-3
# assert abs(y6 - pyo.value(m.y6[180])) / y6 < 1e-3

a = plt.plot(m.t, [pyo.value(m.y[t, 1]) for t in m.t], label="y1")
a = plt.plot(m.t, [pyo.value(m.y[t, 2]) for t in m.t], label="y2")
a = plt.plot(m.t, [pyo.value(m.y[t, 3]) for t in m.t], label="y3")
a = plt.plot(m.t, [pyo.value(m.y[t, 4]) for t in m.t], label="y4")
a = plt.plot(m.t, [pyo.value(m.y[t, 5]) for t in m.t], label="y5")
a = plt.plot(m.t, [pyo.value(m.y6[t]) for t in m.t], label="y6")
a = plt.legend()
a = plt.ylabel("Concentration (mol/l)")
a = plt.xlabel("time (min)")
plt.show()

# First plot all y's on one plot

a = plt.plot(tj.time, tj.get_vec(m.y[180, 1]), label="y1")
a = plt.plot(tj.time, tj.get_vec(m.y[180, 2]), label="y2")
a = plt.plot(tj.time, tj.get_vec(m.y[180, 3]), label="y3")
a = plt.plot(tj.time, tj.get_vec(m.y[180, 4]), label="y4")
a = plt.plot(tj.time, tj.get_vec(m.y[180, 5]), label="y5")
a = plt.plot(tj.time, tj.get_vec(m.y6[180]), label="y6")
a = plt.legend()
a = plt.ylabel("Concentration (mol/l)")
a = plt.xlabel("time (min)")
plt.show()

# 2 and 4 are pretty low concentration, so plot those so we can see better
a = plt.plot(tj.time, tj.get_vec(m.y[180, 2]), label="y2")
a = plt.plot(tj.time, tj.get_vec(m.y[180, 4]), label="y4")
a = plt.legend()
a = plt.ylabel("Concentration (mol/l)")
a = plt.xlabel("time (min)")
plt.show()

# 2 seems to have some fast dynamics so plot a shorter time
a = plt.plot(tj.vecs["_time"], tj.vecs[str(m.y[180, 2])], label="y2")
a = plt.legend()
a = plt.ylabel("Concentration (mol/l)")
a = plt.xlabel("time (min)")
a = plt.xlim(0, 2)

# This creates a new trajectory data set with data every minute.
tji = tj.interpolate(np.linspace(0, 180, 181))

# The plot of this new data should look the same as the original, although some of the
# fast dynamics of component 2 will be obscured.

a = plt.plot(tji.time, tji.get_vec(m.y[180, 1]), label="y1")
a = plt.plot(tji.time, tji.get_vec(m.y[180, 2]), label="y2")
a = plt.plot(tji.time, tji.get_vec(m.y[180, 3]), label="y3")
a = plt.plot(tji.time, tji.get_vec(m.y[180, 4]), label="y4")
a = plt.plot(tji.time, tji.get_vec(m.y[180, 5]), label="y5")
a = plt.plot(tji.time, tji.get_vec(m.y6[180]), label="y6")
a = plt.legend()
a = plt.ylabel("Concentration (mol/l)")
a = plt.xlabel("time (min)")

a = plt.plot(tji.time, tji.get_vec(m.y[180, 2]), label="y2 interpolate dt=1")
a = plt.plot(tj.time, tj.get_vec(m.y[180, 2]), label="y2 original")
a = plt.legend()
a = plt.ylabel("Concentration (mol/l)")
a = plt.xlabel("time (min)")
a = plt.xlim(0, 2)