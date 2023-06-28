# Import plotting functions
import matplotlib.pyplot as plt
import pandas as pd

import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt

dataIP = pd.read_csv("resultsS_IP.csv")
dataPHA = pd.read_csv("resultsX_PHA.csv")
dataPP = pd.read_csv("resultsX_PP.csv")
dataPAO = pd.read_csv("resultsX_PAO.csv")
dataT = pd.read_csv("resultsT.csv")

n_bins = 20

# fig, ax = plt.subplots(figsize=(15, 10))
# plt.plot(dataPAO.Feed, dataPAO.byS_PO4, "*")
# ax.set_xlabel("Inlet X$_{PAO}$ (kg/m$^3$)", fontsize=25, labelpad=20)
# ax.set_ylabel("Outlet S$_{PO4}$ (kg/m$^3$)", fontsize=25, labelpad=20)
# ax.set_ylim(2.0e6, 3.2e6)
# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)
# plt.title("S$_{PO4}$ vs Inlet X$_{PAO}$", fontsize=34)
# plt.savefig("PAOvspo4.png", bbox_inches="tight", dpi=300)
# plt.show()

# fig, ax = plt.subplots(figsize=(15, 10))
# plt.plot(dataPP.Feed, dataPP.byS_PO4, "*")
# ax.set_xlabel("Inlet X$_{PP}$ (kg/m$^3$)", fontsize=25, labelpad=20)
# ax.set_ylabel("Outlet S$_{PO4}$ (kg/m$^3$)", fontsize=25, labelpad=20)
# ax.set_ylim(2.0e6, 3.2e6)
# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)
# plt.title("S$_{PO4}$ vs Inlet X$_{PP}$", fontsize=34)
# plt.savefig("PPvspo4.png", bbox_inches="tight", dpi=300)
# plt.show()

# fig, ax = plt.subplots(figsize=(15, 10))
# plt.plot(dataIP.Feed, dataIP.byS_PO4, "*")
# ax.set_xlabel("Inlet S$_{IP}$ (kg/m$^3$)", fontsize=25, labelpad=20)
# ax.set_ylabel("Outlet S$_{PO4}$ (kg/m$^3$)", fontsize=25, labelpad=20)
# ax.set_ylim(2.0e6, 3.2e6)
# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)
# plt.title("S$_{PO4}$ vs Inlet S$_{IP}$", fontsize=34)
# plt.savefig("IPvspo4.png", bbox_inches="tight", dpi=300)
# plt.show()

# fig, ax = plt.subplots(figsize=(15, 10))
# plt.plot(dataPHA.Feed, dataPHA.byS_PO4, "*")
# ax.set_xlabel("Inlet X$_{PHA}$ (kg/m$^3$)", fontsize=25, labelpad=20)
# ax.set_ylabel("Outlet S$_{PO4}$ (kg/m$^3$)", fontsize=25, labelpad=20)
# ax.set_ylim(2.0e6, 3.2e6)
# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)
# plt.title("S$_{PO4}$ vs Inlet X$_{PHA}$", fontsize=34)
# plt.savefig("PHAvspo4.png", bbox_inches="tight", dpi=300)
# plt.show()

# fig, ax = plt.subplots(figsize=(15, 10))
# plt.plot(dataPP.Feed, dataPP.lcow, "*")
# ax.set_xlabel("Inlet X_$_{PP}$ (kg/m$^3$)", fontsize=25, labelpad=20)
# ax.set_ylabel("LCOW", fontsize=25, labelpad=20)
# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)
# plt.title("LCOW vs Inlet X$_{PP}$", fontsize=34)
# # plt.savefig("cod.png", bbox_inches="tight", dpi=300)
# plt.show()

levels = [
    3.90,
    3.91,
    3.92,
    3.93,
    3.94,
    3.95,
    3.96,
    3.97,
    3.98,
    3.99,
    4,
    4.05,
    4.1,
    4.15,
    4.2,
]
fig, ax = plt.subplots(figsize=(15, 10))
tcf = ax.tricontourf(dataT.Feed, dataT.Temp, dataT.lcow, levels)
ax.set_xlabel("Inlet X$_{PP}$ (kg/m$^3$)", fontsize=25, labelpad=20)
ax.set_ylabel("Temperature (K)", fontsize=25, labelpad=20)
ax.set_xlim(0.95, 1.15)
cbar = fig.colorbar(tcf)
cbar.set_label("LCOW \$/m$^3$", fontsize=25)
# ax.set_ylim(2.0e6, 3.2e6)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.title("S$_{PO4}$ vs Inlet X$_{PP}$", fontsize=34)
plt.savefig("heatT.png", bbox_inches="tight", dpi=300)
plt.show()

# x = np.random.uniform(dataT.Feed)
# y = np.random.uniform(dataT.T)

# triang = tri.Triangulation(x, y)
# interpolator = tri.LinearTriInterpolator(triang, dataT.lcow)
# Xi, Yi = np.meshgrid(dataT.Feed, dataPP.T)
# zi = interpolator(Xi, Yi)

# fig, ax = plt.subplots(figsize=(15, 10))
# ax1.contour(xi, yi, zi, levels=14, linewidths=0.5, colors='k')
# cntr1 = ax1.contourf(xi, yi, zi, levels=14, cmap="RdBu_r")
# fig.colorbar(cntr1, ax=ax1)
# ax.set_xlabel("Inlet X$_{PP}$ (kg/m$^3$)", fontsize=25, labelpad=20)
# ax.set_ylabel("Outlet S$_{PO4}$ (kg/m$^3$)", fontsize=25, labelpad=20)
# ax.set_ylim(2.0e6, 3.2e6)
# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)
# plt.title("S$_{PO4}$ vs Inlet X$_{PP}$", fontsize=34)
# plt.savefig("heatT.png", bbox_inches="tight", dpi=300)
# plt.show()
