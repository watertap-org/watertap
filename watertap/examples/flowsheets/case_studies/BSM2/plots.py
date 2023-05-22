# Import plotting functions
import matplotlib.pyplot as plt
import pandas as pd

import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt

data = pd.read_csv("results.csv")
datat = pd.read_csv("resultst.csv")
datag = pd.read_csv("resultstg.csv")

n_bins = 20

total = datag.ch4 * datag.flow
cod = data.S_I + data.S_S + data.X_I + data.X_S + data.X_BH + data.X_BA + data.X_P


fig, ax = plt.subplots(figsize=(15, 10))
ax.hist(cod, bins=n_bins)
ax.set_xlabel("Treated COD concentration (kg/m$^3$)", fontsize=25, labelpad=20)
ax.set_ylabel("Count", fontsize=25, labelpad=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.title("COD variation", fontsize=34)
plt.savefig("cod_error.png", bbox_inches="tight", dpi=300)
plt.show()

fig, ax = plt.subplots(figsize=(15, 10))
plt.plot(datag.Feedg, total, "*")
ax.set_xlabel("Inlet X$_S$ (kg/m$^3$)", fontsize=25, labelpad=20)
ax.set_ylabel("Gas ADM CH$_4$ (kg/s)", fontsize=25, labelpad=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.title("Gas variation vs inlet X$_S$", fontsize=34)
plt.savefig("gas.png", bbox_inches="tight", dpi=300)
plt.show()

fig, ax = plt.subplots(figsize=(15, 10))
ax.hist(total, bins=n_bins)
ax.set_xlabel("Gas ADM CH$_4$ (kg/s)", fontsize=25, labelpad=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
ax.set_ylabel("Count", fontsize=25, labelpad=20)
plt.title("Gas flow variation", fontsize=34)
plt.savefig("gas_error.png", bbox_inches="tight", dpi=300)
plt.show()

fig, ax = plt.subplots(figsize=(15, 10))
plt.plot(data.Feed, cod, "*")
ax.set_xlabel("Inlet X$_S$ (kg/m$^3$)", fontsize=25, labelpad=20)
ax.set_ylabel("Treated COD concentration (kg/m$^3$)", fontsize=25, labelpad=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.title("COD variation vs inlet X$_S$", fontsize=34)
plt.savefig("cod.png", bbox_inches="tight", dpi=300)
plt.show()
