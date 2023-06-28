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

n_bins = 20

fig, ax = plt.subplots(figsize=(15, 10))
plt.plot(dataPAO.Feed, dataPAO.byS_PO4, "*")
ax.set_xlabel("Inlet X$_{PAO}$ (kg/m$^3$)", fontsize=25, labelpad=20)
ax.set_ylabel("Outlet S$_{PO4}$ (kg/m$^3$)", fontsize=25, labelpad=20)
ax.set_ylim(2.0e6, 3.2e6)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.title("S$_{PO4}$ vs Inlet X$_{PAO}$", fontsize=34)
plt.savefig("PAOvspo4.png", bbox_inches="tight", dpi=300)
plt.show()

fig, ax = plt.subplots(figsize=(15, 10))
plt.plot(dataPP.Feed, dataPP.byS_PO4, "*")
ax.set_xlabel("Inlet X$_{PP}$ (kg/m$^3$)", fontsize=25, labelpad=20)
ax.set_ylabel("Outlet S$_{PO4}$ (kg/m$^3$)", fontsize=25, labelpad=20)
ax.set_ylim(2.0e6, 3.2e6)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.title("S$_{PO4}$ vs Inlet X$_{PP}$", fontsize=34)
plt.savefig("PPvspo4.png", bbox_inches="tight", dpi=300)
plt.show()

fig, ax = plt.subplots(figsize=(15, 10))
plt.plot(dataIP.Feed, dataIP.byS_PO4, "*")
ax.set_xlabel("Inlet S$_{IP}$ (kg/m$^3$)", fontsize=25, labelpad=20)
ax.set_ylabel("Outlet S$_{PO4}$ (kg/m$^3$)", fontsize=25, labelpad=20)
ax.set_ylim(2.0e6, 3.2e6)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.title("S$_{PO4}$ vs Inlet S$_{IP}$", fontsize=34)
plt.savefig("IPvspo4.png", bbox_inches="tight", dpi=300)
plt.show()

fig, ax = plt.subplots(figsize=(15, 10))
plt.plot(dataPHA.Feed, dataPHA.byS_PO4, "*")
ax.set_xlabel("Inlet X$_{PHA}$ (kg/m$^3$)", fontsize=25, labelpad=20)
ax.set_ylabel("Outlet S$_{PO4}$ (kg/m$^3$)", fontsize=25, labelpad=20)
ax.set_ylim(2.0e6, 3.2e6)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.title("S$_{PO4}$ vs Inlet X$_{PHA}$", fontsize=34)
plt.savefig("PHAvspo4.png", bbox_inches="tight", dpi=300)
plt.show()

# fig, ax = plt.subplots(figsize=(15, 10))
# plt.plot(dataPP.Feed, dataPP.lcow, "*")
# ax.set_xlabel("Inlet X_$_{PP}$ (kg/m$^3$)", fontsize=25, labelpad=20)
# ax.set_ylabel("LCOW", fontsize=25, labelpad=20)
# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)
# plt.title("LCOW vs Inlet X$_{PP}$", fontsize=34)
# # plt.savefig("cod.png", bbox_inches="tight", dpi=300)
# plt.show()
