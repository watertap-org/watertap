import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata



data_ec = np.genfromtxt('raw_data_ec.csv', skip_header=1, delimiter=',')
data_lcow = np.genfromtxt('raw_data_lcow.csv', skip_header=1, delimiter=',')

points = data_ec[:, 0:2]

ln_x = np.unique(data_ec[:, 0])
ln_y = np.unique(data_ec[:, 1])

x_mesh, y_mesh = np.meshgrid(ln_x, ln_y)

xi = (x_mesh, y_mesh)

ec_data = griddata(points, data_ec[:, 2], xi, method='nearest')
lcow_data = griddata(points, data_lcow[:, 2], xi, method='nearest')

fig, ax = plt.subplots(1, 2, figsize=(9, 5), dpi=100)
cbar = ax[0].contourf(x_mesh, y_mesh, ec_data)
fig.colorbar(cbar, ax=ax[0])
ax[0].set_xlabel('Feed Concentration')
ax[0].set_ylabel('Water Recovery (%)')

cbar = ax[1].contourf(x_mesh, y_mesh, lcow_data)
fig.colorbar(cbar, ax=ax[1])
ax[1].set_xlabel('Feed Concentration')
ax[1].set_ylabel('Water Recovery (%)')

# plt.plot(points[:, 0], points[:, 1], '.', color='k')
# plt.title('')
# plt.xlabel(data_x)
# plt.ylabel(data_y)
# ax[1].colorbar()
plt.savefig('contour_ec_and_lcow.png')

