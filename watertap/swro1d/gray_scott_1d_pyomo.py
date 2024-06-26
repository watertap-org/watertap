#  ___________________________________________________________________________
#
#  Pyomo: Python Optimization Modeling Objects
#  Copyright (c) 2008-2024
#  National Technology and Engineering Solutions of Sandia, LLC
#  Under the terms of Contract DE-NA0003525 with National Technology and
#  Engineering Solutions of Sandia, LLC, the U.S. Government retains certain
#  rights in this software.
#  This software is distributed under the 3-clause BSD License.
#  ___________________________________________________________________________

# Example 1 from http://www.mathworks.com/help/matlab/ref/pdepe.html

from pyomo.environ import *
from pyomo.dae import *

m = ConcreteModel()
m.pi = Param(initialize=3.1416)
m.f = Param(initialize=0.024)
m.k = Param(initialize=0.055)
m.Du = Param(initialize=0.01)
m.Dv = Param(initialize=0.005)
m.dx = Param(initialize=0.025)

m.t = ContinuousSet(bounds=(0, 4))
m.x = ContinuousSet(bounds=(0, 1))
m.u = Var(m.x, m.t)
m.v = Var(m.x, m.t)

m.dudt = DerivativeVar(m.u, wrt=m.t)
m.dvdt = DerivativeVar(m.v, wrt=m.t)
m.dudx = DerivativeVar(m.u, wrt=m.x)
m.dvdx = DerivativeVar(m.v, wrt=m.x)
m.dudx2 = DerivativeVar(m.u, wrt=(m.x, m.x))
m.dvdx2 = DerivativeVar(m.v, wrt=(m.x, m.x))

def G(u, v, f, k):
    return f * (1 - u) - u*v**2

def H(u, v, f, k):
    return -(f + k) * v + u*v**2

def _diffeq1(m, i, j):
    if i == 0 or i == 1 or j == 0:
        return Constraint.Skip
    return m.dudt[i, j] == m.Du * m.dudx2[i, j] + G(m.u[i, j], m.v[i, j], m.f, m.k)

m.diffeq1 = Constraint(m.x, m.t, rule=_diffeq1)

def _diffeq2(m, i, j):
    if i == 0 or i == 1 or j == 0:
        return Constraint.Skip
    return m.dvdt[i, j] == m.Dv * m.dvdx2[i, j] + H(m.u[i, j], m.v[i, j], m.f, m.k)

m.diffeq2 = Constraint(m.x, m.t, rule=_diffeq2)

def _initcon_1(m, i):
    if i == 0 or i == 1:
        return Constraint.Skip
    return m.u[i, 0] == 3

m.initcon_1 = Constraint(m.x, rule=_initcon_1)

def _initcon_2(m, i):
    if i == 0 or i == 1:
        return Constraint.Skip
    return m.v[i, 0] == 3

m.initcon_2 = Constraint(m.x, rule=_initcon_2)

def _lowerbound_1(m, j):
    return m.dudx[0, j] == 0

m.lowerbound_1 = Constraint(m.t, rule=_lowerbound_1)

def _lowerbound_2(m, j):
    return m.dvdx[0, j] == 0

m.lowerbound_2 = Constraint(m.t, rule=_lowerbound_2)

def _upperbound_1(m, j):
    return m.dudx[1, j] == 0

m.upperbound_1 = Constraint(m.t, rule=_upperbound_1)

def _upperbound_2(m, j):
    return m.dvdx[1, j] == 0

m.upperbound_2 = Constraint(m.t, rule=_upperbound_2)

m.obj = Objective(expr=1)

# Discretize using Orthogonal Collocation
# discretizer = TransformationFactory('dae.collocation')
# discretizer.apply_to(m,nfe=10,ncp=3,wrt=m.x)
# discretizer.apply_to(m,nfe=20,ncp=3,wrt=m.t)

# Discretize using Finite Difference and Collocation
# discretizer = TransformationFactory('dae.finite_difference')
# discretizer2 = TransformationFactory('dae.collocation')
# discretizer.apply_to(m, nfe=25, wrt=m.x, scheme='BACKWARD')
# discretizer2.apply_to(m, nfe=20, ncp=3, wrt=m.t)

# Discretize using Finite Difference Method
discretizer = TransformationFactory('dae.finite_difference')
discretizer.apply_to(m,nfe=4,wrt=m.x,scheme='BACKWARD')
discretizer.apply_to(m,nfe=4,wrt=m.t,scheme='BACKWARD')

solver = SolverFactory('ipopt')
results = solver.solve(m, tee=True)

x = []
t = []
u = []
v = []

for i in sorted(m.x):
    tempu = []
    tempv = []
    tempx = []
    for j in sorted(m.t):
        tempx.append(i)
        tempu.append(value(m.u[i, j]))
        tempv.append(value(m.v[i, j]))
    x.append(tempx)
    t.append(sorted(m.t))
    u.append(tempu)
    v.append(tempv)

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D

x, t, u, v = np.array(x), np.array(t), np.array(u), np.array(v)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection='3d')
ax.set_xlabel('Distance x')
ax.set_ylabel('Time t')
ax.plot_wireframe(x, t, u, rstride=1, cstride=1)
plt.show()
print(t,x)
print(u)