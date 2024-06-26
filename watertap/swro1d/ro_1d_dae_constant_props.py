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
from enum import Enum
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.core.solvers import petsc
from idaes.core.util.model_statistics import degrees_of_freedom

class Solver(Enum):
    ipopt = 1
    petsc = 2

SolverChoice = Solver.petsc

m = ConcreteModel()
m.tf = Param(initialize=0.8e-3) # 8 mm
m.tp = Param(initialize=0.5e-3) # 5 mm
m.L = Param(initialize=0.934)
m.W = Param(initialize=8.4)
m.rhow = Param(initialize=55.56) # kmol/m3
m.R = Param(initialize=0.082) # atm m3/K kmol
m.Pp = Param(initialize=1) # atm
m.b = Param(initialize=8529.45) # atm s/m4
m.Aw = Param(initialize=9.5188e-7) # m/atm s
m.Bs = Param(initialize=8.468e-8) # m/s

# Initial and left boundary conditions
m.Cb0 = Param(initialize=6.226e-3) # kmol/m3
m.Cp0 = Param(initialize=0.0) # kmol/m3
m.Fb0 = Param(initialize=2.166e-4) # m3/s
m.Pb0 = Param(initialize=5.83) # atm
m.Tb0 = Param(initialize=304.65) # C
m.Tp0 = Param(initialize=304.65) # C
m.Db = Param(initialize=1.7657e-9) # m2/s
m.Dp = Param(initialize=1.7354e-9) # m2/s
m.Jw0 = Param(initialize=4.41493388e-06)
m.Js0 = Param(initialize=6.39196682e-09)
m.Cw0 = Param(initialize=7.54837839e-02)
m.kx0 = Param(initialize=1.76938225e-06)

m.t = ContinuousSet(bounds=(0, 1800))
m.x = ContinuousSet(bounds=(0, 1))
m.Cb = Var(m.x, m.t, within=NonNegativeReals, initialize=m.Cb0)
m.Cp = Var(m.x, m.t, within=NonNegativeReals, initialize=m.Cp0)
m.Fb = Var(m.x, m.t, within=NonNegativeReals, initialize=m.Fb0)
m.Pb = Var(m.x, m.t, within=NonNegativeReals, initialize=m.Pb0)
m.Tb = Var(m.x, m.t, within=NonNegativeReals, initialize=m.Tb0)
m.Tp = Var(m.x, m.t, within=NonNegativeReals, initialize=m.Tp0)
m.Jw = Var(m.x, m.t, within=NonNegativeReals, initialize=m.Jw0)
m.Js = Var(m.x, m.t, within=NonNegativeReals, initialize=m.Js0)
m.Cw = Var(m.x, m.t, within=NonNegativeReals, initialize=m.Cw0)
m.kx = Var(m.x, m.t, within=NonNegativeReals, initialize=m.kx0)
# m.Fp = Var(m.x, m.t, within=NonNegativeReals, initialize=0)
# for i in range(9):
    # m.add_component('x%s' % i, Var(m.x, m.t))

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
m.dCbdx2 = DerivativeVar(m.Cb, wrt=(m.x, m.x), initialize=1.0)
m.dCpdx2 = DerivativeVar(m.Cp, wrt=(m.x, m.x), initialize=-1.0)
# m.dFpdx = DerivativeVar(m.Fp, wrt=m.x, initialize=1.0)

def _diffeq1(m, i, j): # Cb PDE
    if i == 0 or i == 1 or j == 0:
        return Constraint.Skip
    return m.dCbdt[i, j] == -(m.Cb[i, j]*m.dFbdx[i, j] + m.Fb[i, j]*m.dCbdx[i, j])/(m.tf*m.W) + m.Db*m.dCbdx2[i, j] - m.Jw[i, j]*m.Cp[i, j]/m.tf

m.diffeq1 = Constraint(m.x, m.t, rule=_diffeq1)

def _diffeq2(m, i, j): # Cp PDE
    if i == 0 or i == 1 or j == 0:
        return Constraint.Skip
    return m.dCpdt[i, j] == (m.Cp[i, j]*m.dFbdx[i, j] - (m.Fb0 - m.Fb[i, j])*m.dCpdx[i, j])/(m.tp*m.W) + m.Dp*m.dCpdx2[i, j] + m.Jw[i, j]*m.Cp[i, j]/m.tf
    # return m.dCpdt[i, j] == -(m.Cp[i, j]*m.dFpdx[i, j] + m.Fp[i, j]*m.dCpdx[i, j])/(m.tp*m.W) + m.Dp[i, j]*m.dCpdx2[i, j] + m.Jw[i, j]*m.Cp[i, j]/m.tf

m.diffeq2 = Constraint(m.x, m.t, rule=_diffeq2)

def _diffeq3(m, i, j): # Fb PDE
    if i == 0 or j == 0:
        return Constraint.Skip
    return m.dFbdt[i, j] == (-m.W*m.Jw[i, j] - m.dFbdx[i, j])*m.Fb[i, j]/(m.tf*m.W)

m.diffeq3 = Constraint(m.x, m.t, rule=_diffeq3)

def _diffeq4(m, i, j): # Pb PDE
    if i == 0 or j == 0:
        return Constraint.Skip
    return m.dPbdt[i, j] == (-m.b*m.Fb[i, j] - m.dPbdx[i, j])*m.Fb[i, j]/(m.tf*m.W)

m.diffeq4 = Constraint(m.x, m.t, rule=_diffeq4)

def _diffeq5(m, i, j): # Tb PDE
    if i == 0 or j == 0:
        return Constraint.Skip
    # return m.dTbdt[i, j] == -(m.Fb[i, j]*m.dTbdx[i, j] + m.dFbdx[i, j]*m.Tb[i, j])/(m.tf*m.W) - m.Jw[i, j]*(m.Tb[i, j]-m.Tp[i, j])/m.tf
    return m.dTbdt[i, j] == -m.Fb[i, j]*m.dTbdx[i, j]/(m.tf*m.W) - m.Jw[i, j]*(m.Tb[i, j]-m.Tp[i, j])/m.tf

m.diffeq5 = Constraint(m.x, m.t, rule=_diffeq5)

def _diffeq6(m, i, j): # Tp PDE
    if i == 0 or j == 0:
        return Constraint.Skip
    return m.dTpdt[i, j] == m.Jw[i, j]*(m.Tb[i, j]-m.Tp[i, j])/m.tf

m.diffeq6 = Constraint(m.x, m.t, rule=_diffeq6)

# def _diffeq7(m, i, j): # Fb PDE for dx
#     if i == 0:
#         return Constraint.Skip
#     return m.dFbdx[i, j] == -m.W*m.Jw[i, j]

# m.diffeq7 = Constraint(m.x, m.t, rule=_diffeq7)

def _algeq1(m, i, j): # Jw Algebraic Equation
    return m.Jw[i, j] == m.Aw*(m.Pb[i, j] - m.Pp - m.R*m.Tb[i, j]*(m.Cw[i, j] - m.Cp[i, j]))

m.algeq1 = Constraint(m.x, m.t, rule=_algeq1)

def _algeq2(m, i, j): # Js Algebraic Equation
    return m.Js[i, j] == m.Bs*exp(m.Jw[i, j]/m.kx[i, j])*(m.Cb[i, j] - m.Cp[i, j])

m.algeq2 = Constraint(m.x, m.t, rule=_algeq2)

def _algeq3(m, i, j): # Cw Algebraic Equation
    return m.Cw[i, j] == m.Cp[i, j] + exp(m.Jw[i, j]/m.kx[i, j])*(m.Cb[i, j] - m.Cp[i, j])

m.algeq3 = Constraint(m.x, m.t, rule=_algeq3)

def _algeq4(m, i, j): # kx Algebraic Equation
    # return m.kx[i, j] == 0.0344832*m.Jw[i, j]**0.739*m.Cb[i, j]**0.135
    return m.kx[i, j] == 0.095360572*m.Fb[i, j]**0.13*m.Jw[i, j]**0.739*m.Cb[i, j]**0.135
    # return m.kx[i, j] == 147.4/(2*m.tf)*m.Db[i, j]*2100*(m.Cb[i, j]/m.rhow)**0.135

m.algeq4 = Constraint(m.x, m.t, rule=_algeq4)

# def _algeq5(m, i, j): # Fp Algebraic Equation
#     return m.Fp[i, j] == m.Fb0 - m.Fb[i, j]

# m.algeq5 = Constraint(m.x, m.t, rule=_algeq5)

# def _initcon_1(m, i):
#     if i == 1:
#         return Constraint.Skip
#     return m.Cb[i, 0] == m.Cb0

# m.initcon_1 = Constraint(m.x, rule=_initcon_1)

# def _initcon_2(m, i):
#     if i == 1:
#         return Constraint.Skip
#     return m.Cp[i, 0] == m.Cp0

# m.initcon_2 = Constraint(m.x, rule=_initcon_2)

# def _initcon_3(m, i):
#     return m.Fb[i, 0] == m.Fb0

# m.initcon_3 = Constraint(m.x, rule=_initcon_3)

# def _initcon_4(m, i):
#     return m.Pb[i, 0] == m.Pb0

# m.initcon_4 = Constraint(m.x, rule=_initcon_4)

# def _initcon_5(m, i):
#     return m.Tb[i, 0] == m.Tb0

# m.initcon_5 = Constraint(m.x, rule=_initcon_5)

# def _initcon_6(m, i):
#     return m.Tp[i, 0] == m.Tp0

# m.initcon_6 = Constraint(m.x, rule=_initcon_6)

def _leftbccon_1(m, j):
    return m.Cb[0, j] == m.Cb0

m.leftbccon_1 = Constraint(m.t, rule=_leftbccon_1)

def _leftbccon_2(m, j):
    return m.Cp[0, j] == m.Cp0

m.leftbccon_2 = Constraint(m.t, rule=_leftbccon_2)

def _leftbccon_3(m, j):
    return m.Fb[0, j] == m.Fb0

m.leftbccon_3 = Constraint(m.t, rule=_leftbccon_3)

def _leftbccon_4(m, j):
    return m.Pb[0, j] == m.Pb0

m.leftbccon_4 = Constraint(m.t, rule=_leftbccon_4)

def _leftbccon_5(m, j):
    return m.Tb[0, j] == m.Tb0

m.leftbccon_5 = Constraint(m.t, rule=_leftbccon_5)

def _leftbccon_6(m, j):
    return m.Tp[0, j] == m.Tp0

m.leftbccon_6 = Constraint(m.t, rule=_leftbccon_6)

def _upperbound_1(m, j):
    return m.dCbdx[1, j] == 0

m.upperbound_1 = Constraint(m.t, rule=_upperbound_1)

def _upperbound_2(m, j):
    return m.dCpdx[1, j] == 0

m.upperbound_2 = Constraint(m.t, rule=_upperbound_2)

# def _leftbound_1(m, j):
#     return m.Jw[0, j] == m.Jw0

# m.leftbound_1 = Constraint(m.t, rule=_leftbound_1)

# def _leftbound_2(m, j):
#     return m.Js[0, j] == m.Js0

# m.leftbound_2 = Constraint(m.t, rule=_leftbound_2)

# def _leftbound_3(m, j):
#     return m.Cw[0, j] == m.Cw0

# m.leftbound_3 = Constraint(m.t, rule=_leftbound_3)

# def _leftbound_4(m, j):
#     return m.kx[0, j] == m.kx0

# m.leftbound_4 = Constraint(m.t, rule=_leftbound_4)

m.obj = Objective(expr=1)

# Discretize using Orthogonal Collocation
# discretizer = TransformationFactory('dae.collocation')
# discretizer.apply_to(m,nfe=5,ncp=3,wrt=m.x)
# discretizer.apply_to(m,nfe=6,ncp=3,wrt=m.t)
# print("Degrees of Freedom before Discretization = ", degrees_of_freedom(m))
# Discretize using Finite Difference and Collocation
discretizer = TransformationFactory('dae.collocation')
discretizer2 = TransformationFactory('dae.finite_difference')
discretizer.apply_to(m, nfe=5, ncp=3, wrt=m.x)
discretizer2.apply_to(m, nfe=30, wrt=m.t, scheme='BACKWARD')

# Discretize using Finite Difference Method
# discretizer = TransformationFactory('dae.finite_difference')
# discretizer.apply_to(m,nfe=6,wrt=m.x,scheme='BACKWARD')
# discretizer.apply_to(m,nfe=120,wrt=m.t,scheme='BACKWARD')

# create the scaling factors
m.scaling_factor = Suffix(direction=Suffix.EXPORT)
m.scaling_factor[m.diffeq1] = 1e4 # scale Cb eq
m.scaling_factor[m.diffeq2] = 1e3 # scale Cp eq
m.scaling_factor[m.diffeq3] = 1e5 # scale Fb eq
m.scaling_factor[m.diffeq3] = 1e-3 # scale Pb eq
m.scaling_factor[m.diffeq5] = 1e1 # scale Tb eq
m.scaling_factor[m.diffeq6] = 1e1 # scale Tp eq
m.scaling_factor[m.algeq1] = 1e7 # scale Jw eq
m.scaling_factor[m.algeq2] = 1e9  # scale Js eq
m.scaling_factor[m.algeq4] = 1e6  # scale kx eq
m.scaling_factor[m.Cb] = 1e4    # scale the Cb variable
m.scaling_factor[m.Cp] = 1e3    # scale the Cp variable
m.scaling_factor[m.Fb] = 1e5    # scale the Fb variable
m.scaling_factor[m.Jw] = 1e7    # scale the Jw variable
m.scaling_factor[m.Js] = 1e9    # scale the Js variable
m.scaling_factor[m.kx] = 1e6    # scale the kx variable
m.scaling_factor[m.dCbdt] = 1e4    # scale the Cb variable
m.scaling_factor[m.dCpdt] = 1e3    # scale the Cp variable
m.scaling_factor[m.dFbdt] = 1e3    # scale the Fb variable
m.scaling_factor[m.dPbdt] = 1e3    # scale the Pb variable
m.scaling_factor[m.dTbdt] = 1e2    # scale the Tb variable
m.scaling_factor[m.dTpdt] = 1e2    # scale the Tp variable
m.scaling_factor[m.dCbdx] = 1e3    # scale the Cb variable
m.scaling_factor[m.dCpdx] = 1e3    # scale the Cp variable
# iscale.calculate_scaling_factors(m)
# scaled_model = TransformationFactory('core.scale_model').create_using(m)
# print("Degrees of Freedom after Discretization = ", degrees_of_freedom(m))

if SolverChoice == Solver.ipopt:

    solver = SolverFactory('ipopt')
    print(solver._version)
    # solver.options['max_iter'] = 1000
    solver.options['linear_solver'] = 'mumps'
    # solver.options['hsllib'] = '/usr/local/fya/watertap/hsl_ma97-2.8.1'
    solver.options['nlp_scaling_method'] = 'user-scaling'
    solver.options['OF_ma57_automatic_scaling'] = 'yes'
    solver.options['halt_on_ampl_error'] = 'no'
    results = solver.solve(scaled_model, tee=True)

else:
    idaeslog.solver_log.tee = True
    results = petsc.petsc_dae_by_time_element(
        m,
        time=m.t,
        keepfiles=True,
        symbolic_solver_labels=True,
        ts_options={
            "--ts_type": "beuler",
            # "-ts_arkimex_type": "1bee",
            "--ts_dt": 0.1,
            "--ts_rtol": 1e-5,
            # "--ts_adapt_clip":"0.001,3600",
            # "--ksp_monitor":"",
            "--ts_adapt_dt_min": 1e-3,
            "--ts_adapt_dt_max": 3600,
            "--snes_type": "newtontr",
            # "--ts_max_reject": 200,
            "--ts_monitor": "",
            "-ts_adapt_monitor": "",
            # "--snes_monitor":"",
            "-snes_converged_reason": "",
            # "-ksp_monitor_true_residual": "",
            # "-ksp_converged_reason": "",
            # "-snes_test_jacobian": "",
            "snes_grid_sequence": "",
            "-pc_type": "svd",
            # "-mat_view": "",
            "--ts_save_trajectory": 1,
            "--ts_trajectory_type": "visualization",
            "--ts_max_snes_failures": 25,
            # "--show_cl":"",
            "-snes_max_it": 50,
            "-snes_rtol": 0,
            "-snes_stol": 0,
            "-snes_atol": 1e-8,
        },
        skip_initial=True,
        initial_solver="ipopt",
        initial_solver_options={
            "constr_viol_tol": 1e-8,
            "nlp_scaling_method": "user-scaling",
            "linear_solver": "mumps",
            "OF_ma57_automatic_scaling": "yes",
            "max_iter": 300,
            "tol": 1e-8,
            "halt_on_ampl_error": "no",
        },
    )
    for result in results.results:
        assert_optimal_termination(result)

# TransformationFactory('core.scale_model').propagate_solution(scaled_model, m)

# solver = SolverFactory('ipopt')
# solver.options['max_iter'] = 1000
# solver.options['nlp_scaling_method'] = 'user-scaling'
# solver.options['OF_ma57_automatic_scaling'] = 'yes'
# solver.options['halt_on_ampl_error'] = 'no'
# results = solver.solve(m, tee=True)

x = []
t = []
Cb = []
Cp = []
Fb = []
Jw = []
Js = []
Cw = []
kx = []

for i in sorted(m.x):
    tempCb = []
    tempCp = []
    tempFb = []
    tempJw = []
    tempx  = []
    tempJs = []
    tempCw = []
    tempkx = []
    for j in sorted(m.t):
        tempx.append(i)
        tempCb.append(value(m.Cb[i, j]))
        tempCp.append(value(m.Cp[i, j]))
        tempFb.append(value(m.Fb[i, j]))
        tempJw.append(value(m.Jw[i, j]))
        tempJs.append(value(m.Js[i, j]))
        tempCw.append(value(m.Cw[i, j]))
        tempkx.append(value(m.kx[i, j]))
    x.append(tempx)
    t.append(sorted(m.t))
    Cb.append(tempCb)
    Cp.append(tempCp)
    Fb.append(tempFb)
    Jw.append(tempJw)
    Js.append(tempJs)
    Cw.append(tempCw)
    kx.append(tempkx)

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D

x, t, Cb, Cp, Fb, Jw, Js, Cw, kx = np.array(x), np.array(t), np.array(Cb), np.array(Cp), np.array(Fb), np.array(Jw), np.array(Js), np.array(Cw), np.array(kx)
Cp_av = np.sum(Cp, axis=0)
print(Cb.shape)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection='3d')
ax.set_xlabel('Distance x')
ax.set_ylabel('Time t')
ax.plot_wireframe(x, t, Cb, rstride=1, cstride=1)
plt.show()
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.set_xlabel('Time t')
ax.set_ylabel('Cb(x=L)')
ax.plot(t[0,:], Cb[-1,:])
# ax.set_xlim([0, 100])
plt.show()
print('Jw0 model vs real: ', Jw[0,0], 4.41493388e-06)
print('Js0 model vs real: ', Js[0,0], 6.39196682e-09)
print('Cw0 model vs real: ', Cw[0,0], 7.54837839e-02)
print('kx0 model vs real: ', kx[0,0], 1.76938225e-06)
print('Min of params error at 0: ', min([abs(Jw[0,0] - 4.41493388e-06), abs(Js[0,0] - 6.39196682e-09), abs(Cw[0,0] - 7.54837839e-02), abs(kx[0,0] - 1.76938225e-06)]))
print('Cb[x=0,t=1]: ', Cb[0,1])
print('Cb[x=L,t=End]: ', Cb[-1,-1], '0.00675')
print('Fb[x=L,t=End]: ', Fb[-1,-1], '0.0001951')
print('Cpav[t=End]: ', Cp_av[-1], '0.00182')