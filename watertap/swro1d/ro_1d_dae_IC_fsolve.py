import numpy as np
from scipy.optimize import fsolve
def func(x):
    Aw = 9.5188e-7
    Bs = 8.468e-8
    R = 0.082
    Pp = 1
    Pb0 = 5.83
    Tb0 = 31.5+273.15
    Cb0 = 6.226e-3
    Cp0 = 0.0
    Fb0 = 2.166e-4
    return [x[0] - Aw*(Pb0-Pp-R*Tb0*(x[2]-Cp0)),
            x[1] - Bs*np.exp(x[0]/x[3])*(Cb0-Cp0),
            x[2] - Cp0 - np.exp(x[0]/x[3])*(Cb0-Cp0),
            x[3] - 147.4/0.0016*1.7657e-9*(300332.247*Fb0)**0.13*(1261.395437*x[0])**0.739*(Cb0/55.56)**0.135]
            # x[3] - 0.0173698*(x[0])**0.739]
root = fsolve(func, [1e-6, 1e-9, 1e-2, 1e-6], maxfev=1000, full_output=True)
print('Jw(0), Js(0), Cw(0), kx(0): ', root[0])
print('Function vector: ', root[1]['fvec'])
print('Message: ', root[3])