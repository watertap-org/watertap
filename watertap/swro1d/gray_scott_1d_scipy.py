import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
from scipy.integrate import odeint

def G(u, v, f, k):
    return f * (1 - u) - u*v**2

def H(u, v, f, k):
    return -(f + k) * v + u*v**2

def grayscott1d(y, t, f, k, Du, Dv, dx):
    """
    Differential equations for the 1-D Gray-Scott equations.

    The ODEs are derived using the method of lines.
    """
    # The vectors u and v are interleaved in y.  We define
    # views of u and v by slicing y.
    u = y[::2]
    v = y[1::2]

    # dydt is the return value of this function.
    dydt = np.empty_like(y)

    # Just like u and v are views of the interleaved vectors
    # in y, dudt and dvdt are views of the interleaved output
    # vectors in dydt.
    dudt = dydt[::2]
    dvdt = dydt[1::2]

    # Compute du/dt and dv/dt.  The end points and the interior points
    # are handled separately.
    dudt[0]    = G(u[0],    v[0],    f, k) + Du * (-2.0*u[0] + 2.0*u[1]) / dx**2
    dudt[1:-1] = G(u[1:-1], v[1:-1], f, k) + Du * np.diff(u,2) / dx**2
    dudt[-1]   = G(u[-1],   v[-1],   f, k) + Du * (- 2.0*u[-1] + 2.0*u[-2]) / dx**2
    dvdt[0]    = H(u[0],    v[0],    f, k) + Dv * (-2.0*v[0] + 2.0*v[1]) / dx**2
    dvdt[1:-1] = H(u[1:-1], v[1:-1], f, k) + Dv * np.diff(v,2) / dx**2
    dvdt[-1]   = H(u[-1],   v[-1],   f, k) + Dv * (-2.0*v[-1] + 2.0*v[-2]) / dx**2

    return dydt

rng = np.random.default_rng(0)
num_discretize_x = 100 #2500 gives their solution
y0 = rng.standard_normal(2*num_discretize_x)
y0 = np.ones(2*num_discretize_x)*3
t = np.linspace(0, 4, 100)
f = 0.024
k = 0.055
Du = 0.01
Dv = 0.005
dx = 0.01
x = np.linspace(0, 1, num_discretize_x)
t0 = time.time()
#sola = odeint(grayscott1d, y0, t, args=(f, k, Du, Dv, dx))
#t1 = time.time()
solb = odeint(grayscott1d, y0, t, args=(f, k, Du, Dv, dx), ml=2, mu=2)
print('solb shape: ', solb.shape)
t2 = time.time()
#print(f'No banding takes {t1-t0} s, while banding takes {t2-t1} s.')

u = solb[:,::2]
v = solb[:,1::2]

print(u.T)

t_grid, x_grid = np.meshgrid(t, x)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(t_grid, x_grid, u.T)
plt.xlabel('$t$')
plt.ylabel('$x$')
plt.show()
plt.close('all')


plt.plot(t, u[:,-1], label='$u(x=L)$')
plt.plot(t, v[:,-1], label='$v(x=L)$')
plt.xlabel('$t$')
plt.ylabel('$u$ or $v$')
plt.legend(loc=0)
plt.show()
