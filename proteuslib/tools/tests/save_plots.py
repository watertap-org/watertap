import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

data = np.genfromtxt('raw_data.csv', skip_header=1, delimiter=',')

# id_x = param_keys.index(data_x)
# id_y = param_keys.index(data_y)


# ln_x = sweep_params[data_x]
# ln_y = sweep_params[data_y]

ln_x = np.unique(data[:, 0])
ln_y = np.unique(data[:, 1])


x_mesh, y_mesh = np.meshgrid(ln_x, ln_y)

points = data[:, 0:2]
xi = (x_mesh, y_mesh)
results_global = data[:, 2]

c_data = griddata(points, results_global, xi, method='nearest')

print(c_data)

plt.contourf(x_mesh, y_mesh, c_data)
plt.plot(points[:, 0], points[:, 1], '.', color='k')
# plt.title('')
# plt.xlabel(data_x)
# plt.ylabel(data_y)
plt.colorbar()
plt.savefig('contour.png')
plt.savefig('contour.pdf')

