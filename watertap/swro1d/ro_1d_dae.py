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
from pyomo.contrib.parmest.utils.ipopt_solver_wrapper import ipopt_solve_with_stats

m = ConcreteModel()
m.t_f = Param(initialize=100) # in seconds
m.tf = Param(initialize=0.8e-3) # 8 mm
m.tp = Param(initialize=0.5e-3) # 5 mm
m.L = Param(initialize=0.934)
m.W = Param(initialize=8.4)
m.rhow = Param(initialize=55.56) # kmol/m3
m.R = Param(initialize=0.082) # atm m3/K kmol
m.Pp = Param(initialize=1.0) # atm
m.b = Param(initialize=8529.45) # atm s/m4
m.Aw = Param(initialize=9.5188e-7) # m/atm s
m.Bs = Param(initialize=8.468e-8) # m/s
m.deb = Param(initialize=2*m.tf)
m.dep = Param(initialize=2*m.tp)

# Initial and left boundary conditions
m.Cb0 = Param(initialize=6.226e-3) # kmol/m3
m.Cp0 = Param(initialize=0.) # kmol/m3
m.Fb0 = Param(initialize=2.166e-4) # m3/s
m.Pb0 = Param(initialize=5.83) # atm
m.Tb0 = Param(initialize=31) # C
m.Tp0 = Param(initialize=31) # C
m.Jw0 = Param(initialize=4.44604040e-06)
m.Js0 = Param(initialize=5.303351e-09)
m.Cw0 = Param(initialize=6.262814e-02)
m.kx0 = Param(initialize=1.92595956e-06)

m.t = ContinuousSet(bounds=(0, m.t_f))
m.x = ContinuousSet(bounds=(0, 1))
m.Cb = Var(m.t, m.x, within=NonNegativeReals, initialize=m.Cb0, bounds=(0.0, 0.5))
m.Cp = Var(m.t, m.x, within=NonNegativeReals, initialize=m.Cp0, bounds=(0.0, 0.5))
m.Fb = Var(m.t, m.x, within=NonNegativeReals, initialize=m.Fb0)
m.Pb = Var(m.t, m.x, within=Reals, initialize=m.Pb0)
m.Tb = Var(m.t, m.x, within=NonNegativeReals, initialize=m.Tb0)
m.Tp = Var(m.t, m.x, within=NonNegativeReals, initialize=m.Tp0)
m.Jw = Var(m.t, m.x, within=NonNegativeReals, initialize=m.Jw0, bounds=(1e-10, 1e-3))
m.Js = Var(m.t, m.x, within=NonNegativeReals, initialize=m.Js0)
m.Cw = Var(m.t, m.x, within=NonNegativeReals, initialize=m.Cw0)
m.Db = Var(m.t, m.x, within=NonNegativeReals, initialize=1.8e-9)
m.Dp = Var(m.t, m.x, within=NonNegativeReals, initialize=1.8e-9)
m.kx = Var(m.t, m.x, within=NonNegativeReals, initialize=m.kx0)
m.Reb = Var(m.t, m.x, within=NonNegativeReals, initialize=120.)
m.Rep = Var(m.t, m.x, within=NonNegativeReals, initialize=0.005607)
m.mf = Var(m.t, m.x, within=NonNegativeReals, initialize=1.)
m.mp = Var(m.t, m.x, within=NonNegativeReals, initialize=1.)
m.rhob = Var(m.t, m.x, within=NonNegativeReals, initialize=995.)
m.rhop = Var(m.t, m.x, within=NonNegativeReals, initialize=995.)
m.mub = Var(m.t, m.x, within=NonNegativeReals, initialize=8e-4)
m.mup = Var(m.t, m.x, within=NonNegativeReals, initialize=8e-4)
# for i in range(9):
    # m.add_component('x%s' % i, Var(m.t, m.x))

m.dCbdt = DerivativeVar(m.Cb, wrt=m.t, initialize=-1.0)
m.dCpdt = DerivativeVar(m.Cp, wrt=m.t, initialize=1.0)
m.dFbdt = DerivativeVar(m.Fb, wrt=m.t, initialize=-1.0)
m.dPbdt = DerivativeVar(m.Pb, wrt=m.t, initialize=-1.0)
m.dTbdt = DerivativeVar(m.Tb, wrt=m.t, initialize=-1.0)
m.dTpdt = DerivativeVar(m.Tp, wrt=m.t, initialize=1.0)
m.dCbdx = DerivativeVar(m.Cb, wrt=m.x, initialize=-1.0)
m.dCpdx = DerivativeVar(m.Cp, wrt=m.x, initialize=1.0)
m.dFbdx = DerivativeVar(m.Fb, wrt=m.x, initialize=-1.0)
m.dPbdx = DerivativeVar(m.Pb, wrt=m.x, initialize=-1.0)
m.dTbdx = DerivativeVar(m.Tb, wrt=m.x, initialize=-1.0)
m.dTpdx = DerivativeVar(m.Tp, wrt=m.x, initialize=-1.0)
m.dCbdx2 = DerivativeVar(m.Cb, wrt=(m.x, m.x), initialize=0)
m.dCpdx2 = DerivativeVar(m.Cp, wrt=(m.x, m.x), initialize=0)
m.dDbdx = DerivativeVar(m.Db, wrt=m.x, initialize=-1.0)
m.dDpdx = DerivativeVar(m.Dp, wrt=m.x, initialize=1.0)

def _diffeq1(m, t, x): # Cb PDE
    if x == 0 or x == 1 or t == 0:
        return Constraint.Skip
    return m.dCbdt[t, x] == -(m.Cb[t, x]*m.dFbdx[t, x] + m.Fb[t, x]*m.dCbdx[t, x])/(m.tf*m.W) + m.Db[t, x]*m.dCbdx2[t, x] + m.dDbdx[t, x]*m.dCbdx[t, x] - m.Jw[t, x]*m.Cp[t, x]/m.tf

m.diffeq1 = Constraint(m.t, m.x, rule=_diffeq1)

def _diffeq2(m, t, x): # Cp PDE
    if x == 0 or x == 1 or t == 0:
        return Constraint.Skip
    return m.dCpdt[t, x] == (m.Cp[t, x]*m.dFbdx[t, x] - (m.Fb0 - m.Fb[t, x])*m.dCpdx[t, x])/(m.tp*m.W) + m.Dp[t, x]*m.dCpdx2[t, x] + m.dDpdx[t, x]*m.dCpdx[t, x] + m.Jw[t, x]*m.Cp[t, x]/m.tf

m.diffeq2 = Constraint(m.t, m.x, rule=_diffeq2)

def _diffeq3(m, t, x): # Fb PDE
    if x == 0 or t == 0:
        return Constraint.Skip
    return m.dFbdt[t, x] == (-m.W*m.Jw[t, x] - m.dFbdx[t, x])*m.Fb[t, x]/(m.tf*m.W)

m.diffeq3 = Constraint(m.t, m.x, rule=_diffeq3)

def _diffeq4(m, t, x): # Pb PDE
    if x == 0 or t == 0:
        return Constraint.Skip
    return m.dPbdt[t, x] == (-m.b*m.Fb[t, x] - m.dPbdx[t, x])*m.Fb[t, x]/(m.tf*m.W)

m.diffeq4 = Constraint(m.t, m.x, rule=_diffeq4)

def _diffeq5(m, t, x): # Tb PDE
    if x == 0 or t == 0:
        return Constraint.Skip
    # return m.dTbdt[t, x] == -(m.Fb[t, x]*m.dTbdx[t, x] + m.dFbdx[t, x]*m.Tb[t, x])/(m.tf*m.W) - m.Jw[t, x]*(m.Tb[t, x]-m.Tp[t, x])/m.tf
    return m.dTbdt[t, x] == -m.Fb[t, x]*m.dTbdx[t, x]/(m.tf*m.W) - m.Jw[t, x]*(m.Tb[t, x]-m.Tp[t, x])/m.tf

m.diffeq5 = Constraint(m.t, m.x, rule=_diffeq5)

def _diffeq6(m, t, x): # Tp PDE
    if x == 0 or t == 0:
        return Constraint.Skip
    return m.dTpdt[t, x] == m.Jw[t, x]*(m.Tb[t, x]-m.Tp[t, x])/m.tf

m.diffeq6 = Constraint(m.t, m.x, rule=_diffeq6)

def _algeq1(m, t, x): # Jw Algebraic Equation
    return m.Jw[t, x] == m.Aw*(m.Pb[t, x] - m.Pp - m.R*m.Tb[t, x]*(m.Cw[t, x] - m.Cp[t, x]))

m.algeq1 = Constraint(m.t, m.x, rule=_algeq1)

def _algeq2(m, t, x): # Js Algebraic Equation
    return m.Js[t, x] == m.Bs*exp(m.Jw[t, x]/m.kx[t, x])*(m.Cb[t, x] - m.Cp[t, x])

m.algeq2 = Constraint(m.t, m.x, rule=_algeq2)

def _algeq3(m, t, x): # Cw Algebraic Equation
    return m.Cw[t, x] == m.Cp[t, x] + exp(m.Jw[t, x]/m.kx[t, x])*(m.Cb[t, x] - m.Cp[t, x])
    # return m.Cw[t, x] == m.Js[t, x]/m.Bs + m.Cp[t, x]
    # return m.Js[t, x] == m.Jw[t, x]*m.Cp[t, x]

m.algeq3 = Constraint(m.t, m.x, rule=_algeq3)

def _algeq4(m, t, x): # Db Algebraic Equation
    return m.Db[t, x] == 6.725e-6*exp(0.1546e-3*m.Cb[t, x]*18.01253-2513/(m.Tb[t, x]+273.15))

m.algeq4 = Constraint(m.t, m.x, rule=_algeq4)

def _algeq5(m, t, x): # Dp Algebraic Equation
    return m.Dp[t, x] == 6.725e-6*exp(0.1546e-3*m.Cp[t, x]*18.01253-2513/(m.Tp[t, x]+273.15))

m.algeq5 = Constraint(m.t, m.x, rule=_algeq5)

def _algeq6(m, t, x): # kx Algebraic Equation
    return m.kx[t, x] == 147.4/m.deb*m.Db[t, x]*m.Reb[t, x]**0.13*m.Rep[t, x]**0.739*(m.Cb[t, x]/m.rhow)**0.135

m.algeq6 = Constraint(m.t, m.x, rule=_algeq6)

def _algeq7(m, t, x): # Reb Algebraic Equation
    return m.Reb[t, x] == m.rhob[t, x]*m.deb*m.Fb[t, x]/(m.tf*m.W*m.mub[t, x])

m.algeq7 = Constraint(m.t, m.x, rule=_algeq7)

def _algeq8(m, t, x): # Rep Algebraic Equation
    return m.Rep[t, x] == m.rhop[t, x]*m.dep*m.Jw[t, x]/m.mup[t, x]

m.algeq8 = Constraint(m.t, m.x, rule=_algeq8)

def _algeq9(m, t, x): # rhob Algebraic Equation
    return m.rhob[t, x] == 498.4*m.mf[t, x] + sqrt(248400*m.mf[t, x]**2 + 752.4*m.mf[t, x]*m.Cb[t, x]*18.0153)

m.algeq9 = Constraint(m.t, m.x, rule=_algeq9)

def _algeq10(m, t, x): # rhop Algebraic Equation
    return m.rhop[t, x] == 498.4*m.mp[t, x] + sqrt(248400*m.mp[t, x]**2 + 752.4*m.mp[t, x]*m.Cp[t, x]*18.0153)

m.algeq10 = Constraint(m.t, m.x, rule=_algeq10)

def _algeq11(m, t, x): # mf Algebraic Equation
    return m.mf[t, x] == 1.0069 - 2.757e-4*m.Tb[t, x]

m.algeq11 = Constraint(m.t, m.x, rule=_algeq11)

def _algeq12(m, t, x): # mp Algebraic Equation
    return m.mp[t, x] == 1.0069 - 2.757e-4*m.Tp[t, x]

m.algeq12 = Constraint(m.t, m.x, rule=_algeq12)

def _algeq13(m, t, x): # mub Algebraic Equation
    return m.mub[t, x] == 1.234e-6*exp(0.0212e-3*m.Cb[t, x]*18.0153 + 1965/(m.Tb[t, x] + 273.15))

m.algeq13 = Constraint(m.t, m.x, rule=_algeq13)

def _algeq14(m, t, x): # mup Algebraic Equation
    return m.mup[t, x] == 1.234e-6*exp(0.0212e-3*m.Cp[t, x]*18.0153 + 1965/(m.Tp[t, x] + 273.15))

m.algeq14 = Constraint(m.t, m.x, rule=_algeq14)

# def _algeq15(m, t, x): # Fb Algebraic Equation
#     return m.dFbdx[t, x] == -m.W*m.Jw[t, x]

# m.algeq15 = Constraint(m.t, m.x, rule=_algeq15)

def _initcon_1(m, x):
    if x == 1:
        return Constraint.Skip
    return m.Cb[0, x] == m.Cb0

m.initcon_1 = Constraint(m.x, rule=_initcon_1)

def _initcon_2(m, x):
    if x == 1:
        return Constraint.Skip
    return m.Cp[0, x] == m.Cp0

m.initcon_2 = Constraint(m.x, rule=_initcon_2)

def _initcon_3(m, x):
    return m.Fb[0, x] == m.Fb0

m.initcon_3 = Constraint(m.x, rule=_initcon_3)

def _initcon_4(m, x):
    return m.Pb[0, x] == m.Pb0

m.initcon_4 = Constraint(m.x, rule=_initcon_4)

def _initcon_5(m, x):
    return m.Tb[0, x] == m.Tb0

m.initcon_5 = Constraint(m.x, rule=_initcon_5)

def _initcon_6(m, x):
    return m.Tp[0, x] == m.Tp0

m.initcon_6 = Constraint(m.x, rule=_initcon_6)

def _upperbound_1(m, t):
    return m.dCbdx[t, 1] == 0

m.upperbound_1 = Constraint(m.t, rule=_upperbound_1)

def _upperbound_2(m, t):
    return m.dCpdx[t, 1] == 0

m.upperbound_2 = Constraint(m.t, rule=_upperbound_2)

m.obj = Objective(expr=1)

# Discretize using Orthogonal Collocation
# discretizer = TransformationFactory('dae.collocation')
# discretizer.apply_to(m,nfe=6,ncp=3,wrt=m.x)
# discretizer.apply_to(m,nfe=6,ncp=3,wrt=m.t)

# Discretize using Finite Difference and Collocation
discretizer = TransformationFactory('dae.collocation')
discretizer2 = TransformationFactory('dae.finite_difference')
discretizer.apply_to(m, nfe=3, ncp=3, wrt=m.t)
discretizer2.apply_to(m, nfe=3, wrt=m.x, scheme='BACKWARD')

# Discretize using Finite Difference Method
# discretizer = TransformationFactory('dae.finite_difference')
# discretizer.apply_to(m,nfe=6,wrt=m.x,scheme='BACKWARD')
# discretizer.apply_to(m,nfe=6,wrt=m.t,scheme='BACKWARD')

# create the scaling factors
m.scaling_factor = Suffix(direction=Suffix.EXPORT)
m.scaling_factor[m.diffeq1] = 1e3 # scale Cb
m.scaling_factor[m.diffeq2] = 1e2 # scale Cp
m.scaling_factor[m.algeq1] = 1e3 # scale Jw eq
m.scaling_factor[m.algeq2] = 1e8  # scale Js eq
m.scaling_factor[m.algeq4] = 1e5  # scale kx eq
m.scaling_factor[m.algeq6] = 1e5  # scale 
m.scaling_factor[m.algeq7] = 1e5  # scale
m.scaling_factor[m.Cb] = 1e4    # scale the Cb variable
m.scaling_factor[m.Cp] = 1e4    # scale the Cp variable
m.scaling_factor[m.Fb] = 1e2    # scale the Fb variable
# m.scaling_factor[m.Tb] = 1e-2    # scale the Tb variable
# m.scaling_factor[m.Tp] = 1e-2    # scale the Tp variable
m.scaling_factor[m.Jw] = 1e9    # scale the Jw variable
m.scaling_factor[m.Js] = 1e8    # scale the Js variable
m.scaling_factor[m.kx] = 1e9    # scale the kx variable
# transform the model
scaled_model = TransformationFactory('core.scale_model').create_using(m)

solver = SolverFactory('ipopt')
# solver.options['max_iter'] = 1000
solver.options['nlp_scaling_method'] = 'user-scaling'
solver.options['OF_ma57_automatic_scaling'] = 'yes'
solver.options['halt_on_ampl_error'] = 'no'
# solver.options['mu_max_fact'] = 1e3
results = solver.solve(scaled_model, tee=True)
# scaled_model.pprint()
# status_obj, solved, iters, time, regu = ipopt_solve_with_stats(scaled_model, solver, max_iter=3000, max_cpu_time=120)

TransformationFactory('core.scale_model').propagate_solution(scaled_model, m)
# print('Status obj: ', status_obj, 'Solved: ', solved, 'iters: ', iters, 'regu: ', regu)


x = []
t = []
Cb = []
Cp = []
Jw = []
Js = []
Cw = []
kx = []
Db = []
rhob = []
mub = []

for i in sorted(m.x):
    tempCb = []
    tempCp = []
    tempJw = []
    tempx  = []
    tempJs = []
    tempCw = []
    tempkx = []
    tempDb = []
    temprhob = []
    tempmub = []
    for j in sorted(m.t):
        tempx.append(i)
        tempCb.append(value(m.Cb[i, j]))
        tempCp.append(value(m.Cp[i, j]))
        tempJw.append(value(m.Jw[i, j]))
        tempJs.append(value(m.Js[i, j]))
        tempCw.append(value(m.Cw[i, j]))
        tempkx.append(value(m.kx[i, j]))
        tempDb.append(value(m.Db[i, j]))
        temprhob.append(value(m.rhob[i, j]))
        tempmub.append(value(m.mub[i, j]))
    x.append(tempx)
    t.append(sorted(m.t))
    Cb.append(tempCb)
    Cp.append(tempCp)
    Jw.append(tempJw)
    Js.append(tempJs)
    Cw.append(tempCw)
    kx.append(tempkx)
    Db.append(tempDb)
    rhob.append(temprhob)
    mub.append(tempmub)

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D

x, t, Cb, Cp, Jw, Js, Cw, kx, Db, rhob, mub = np.array(x), np.array(t), np.array(Cb), np.array(Cp), np.array(Jw), np.array(Js), np.array(Cw), np.array(kx), np.array(Db), np.array(rhob), np.array(mub)
Cp_av = np.sum(Cp, axis=0)
print(Cb.shape)

# fig = plt.figure()
# ax = fig.add_subplot(1, 1, 1, projection='3d')
# ax.set_xlabel('Distance x')
# ax.set_ylabel('Time t')
# ax.plot_wireframe(x, t, Cp, rstride=1, cstride=1)
# plt.show()
# fig = plt.figure()
# ax = fig.add_subplot(1, 1, 1)
# ax.set_4xlabel('Time t')
# ax.set_ylabel('Cp(av)')
# ax.plot(t[0,:], Cp_av)
# plt.show()
print('Jw0 model vs real: ', Jw[0,0], 4.44604040e-06)
print('Js0 model vs real: ', Js[0,0], 5.303351e-09)
print('Cw0 model vs real: ', Cw[0,0], 0.06262814)
print('kx0 model vs real: ', kx[0,0], 1.92595956e-06)
print('Cb[x=0,t=1]: ', Cb[0,1])
print('Cb[x=L,t=End]: ', Cb[-1,-1])
print('Cpav[t=End]: ', Cp_av[-1])