import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
from scipy.integrate import odeint

def diff_1st(lst):
    lst = np.diff(lst)
    _sum = [] 
    for i in range(len(lst)-1): 
        # adding the alternate numbers 
        _sum.append(lst[i] + lst[i + 1])
    _sum = np.array(_sum)
    return _sum

def ro1d(y, t, W, Aw, Pp, R, b, tf, dx):
    """
    Differential equations for the 1-D Gray-Scott equations.

    The ODEs are derived using the method of lines.
    """
    # The vectors u and v are interleaved in y.  We define
    # views of u and v by slicing y.
    Jw = y[::9]
    Js = y[1::9]
    Cw = y[2::9]
    Cb = y[3::9]
    Cp = y[4::9]
    Fb = y[5::9]
    Pb = y[6::9]
    Tb = y[7::9]
    Tp = y[8::9]

    # dydt is the return value of this function.
    dydt = np.empty_like(y)
    dydt[:9] = np.zeros(9) # Dirichlet BC on LHS

    # Just like u and v are views of the interleaved vectors
    # in y, dudt and dvdt are views of the interleaved output
    # vectors in dydt.
    dJwdt = dydt[::9]
    dJsdt = dydt[1::9]
    dCwdt = dydt[2::9]
    dCbdt = dydt[3::9]
    dCpdt = dydt[4::9]
    dFbdt = dydt[5::9]
    dPbdt = dydt[6::9]
    dTbdt = dydt[7::9]
    dTpdt = dydt[8::9]
    if any(t < 0 for t in np.concatenate((Fb,Jw,Cb))):
        print('Negative in Fb: ', Fb[Fb<=0])
        print('Negative in Jw: ', Jw[Jw<=0])
        print('Negative in Cb: ', Cb[Cb<=0])
    kx = 0.095360572*Fb**0.13*Jw**0.739*Cb**0.135
    # Compute du/dt and dv/dt.  The end points and the interior points
    # are handled separately.
    dJwdt[1:]    = (Aw*(Pb[1:] - Pp - R*Tb[1:]*(Cw[1:] - Cp[1:]))-Jw[1:])*Fb[1:]/(tf*W*dx)
    dJsdt[1:]    = ((Bs*np.exp(Jw[1:]/kx[1:])*(Cb[1:] - Cp[1:]))-Js[1:])*Fb[1:]/(tf*W*dx)
    dCwdt[1:]    = (Cp[1:]+np.exp(Jw[1:]/kx[1:])*(Cb[1:] - Cp[1:])-Cw[1:])*Fb[1:]/(tf*W*dx)
    dCbdt[1:-1]  = -(Cb[1:-1]*diff_1st(Fb) + Fb[1:-1]*diff_1st(Cb))/(tf*W*2*dx) \
        + Db * np.diff(Cb,2) / dx**2 - Jw[1:-1]*Cp[1:-1]/tf
    dCbdt[-1]    = -(Cb[-1]*(Fb[-1]-Fb[-2]) + Fb[-1]*(Cb[-1]-Cb[-2]))/(tf*W*2*dx) \
        + Db * (Cb[-1]-2*Cb[-2]+Cb[-3]) / dx**2 - Jw[-1]*Cp[-1]/tf
    dCpdt[1:-1]  = (Cp[1:-1]*diff_1st(Fb) - (Fb[0] - Fb[1:-1])*diff_1st(Cp))/(tp*W*2*dx) \
        + Dp * np.diff(Cp,2) / dx**2 + Jw[1:-1]*Cp[1:-1]/tf
    dCpdt[-1]    = (Cp[-1]*(Fb[-1]-Fb[-2]) - (Fb[0] - Fb[-1])*(Cp[-1]-Cp[-2]))/(tp*W*2*dx) \
        + Dp * (Cp[-1]-2*Cp[-2]+Cp[-3]) / dx**2 + Jw[-1]*Cp[-1]/tf
    dFbdt[1:-1]  = ((W*Jw[1:-1]-diff_1st(Fb))/(2*dx))*Fb[1:-1]/(tf*W)
    dFbdt[-1]    = ((W*Jw[-1]-(Fb[-1]-Fb[-2]))/(2*dx))*Fb[-1]/(tf*W)
    dPbdt[1:-1]  = ((b*Fb[1:-1]-diff_1st(Pb))/(2*dx))*Fb[1:-1]/(tf*W)
    dPbdt[-1]    = ((b*Fb[-1]-(Pb[-1]-Pb[-2]))/(2*dx))*Fb[-1]/(tf*W)
    dTbdt[1:]    = Fb[1:]/(tf*W*dx)*(Tb[0]-Tb[1:]) - Jw[1:]/tf*(Tb[1:]-Tp[1:])
    dTpdt[1:]    = Jw[1:]/tf*(Tb[1:]-Tp[1:])
    return dydt

rng = np.random.default_rng(0)
num_discretize_x = 5 #2500 gives their solution
t = np.linspace(0, 180, 100)

tf = 0.8e-3 # 8 mm
tp = 0.5e-3 # 5 mm
L = 0.934
W = 8.4
rhow = 55.56 # kmol/m3
R = 0.082 # atm m3/K kmol
Pp = 1 # atm
b = 8529.45 # atm s/m4
Aw = 9.5188e-7 # m/atm s
Bs = 8.468e-8 # m/s

# Initial and left boundary conditions
Cb0 = 6.226e-3 # kmol/m3
Cp0 = 0.0 # kmol/m3
Fb0 = 2.33e-4 # m3/s
Pb0 = 5.83 # atm
Tb0 = 31 # C
Tp0 = 31 # C
Db = 1.7657e-9 # m2/s
Dp = 1.7354e-9 # m2/s
Jw0 = 4.41493388e-06
Js0 = 6.39196682e-09
Cw0 = 7.54837839e-02
kx0 = 1.76938225e-06
Db = 1.7657e-9 # m2/s
Dp =1.7354e-9 # m2/s
dx = 0.0025
x = np.linspace(0, 1, num_discretize_x)
t0 = time.time()
y0 = [Jw0, Js0, Cw0, Cb0, Cp0, Fb0, Pb0, Tb0, Tp0]*5
# y0[:9] = Jw0, Js0, Cw0, Cb0, Cp0, Fb0, Pb0, Tb0, Tp0
sola = odeint(ro1d, y0, t, args=(W, Aw, Pp, R, b, tf, dx))
t1 = time.time()
# solb = odeint(ro1d, y0, t, args=(W, Aw, Pp, R, b, tf, dx), ml=2, mu=2)
# print('solb shape: ', solb.shape)
# t2 = time.time()
#print(f'No banding takes {t1-t0} s, while banding takes {t2-t1} s.')

Cb = sola[:,3::9]
Cp = sola[:,4::9]
Pb = sola[:,6::9]

# print(Cb.T)

t_grid, x_grid = np.meshgrid(t, x)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(t_grid, x_grid, Pb.T)
plt.xlabel('$t$')
plt.ylabel('$x$')
plt.show()
plt.close('all')

plt.plot(t, Cb[:, -1], label='$Cb(x=L)$')
plt.plot(t, Cp[:, -1], label='$Cp(x=L)$')
plt.xlabel('$t$')
plt.ylabel('$u$ or $v$')
plt.legend(loc=0)
plt.show()