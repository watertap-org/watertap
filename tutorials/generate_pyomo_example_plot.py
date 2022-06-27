###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################
import matplotlib.pyplot as plt
import numpy as np


def main():
    v = np.linspace(-1, 1, 1000)
    x1, x2 = np.meshgrid(v, v)
    z = x1**2 + x2**2

    def _z(x1, x2):
        if x1 + 2 * x2 >= 1:
            return x1**2 + x2**2
        else:
            return float("nan")

    z = np.vectorize(_z)(x1, x2)
    plt.figure(figsize=(10, 8))
    plt.rc("font", size=18)
    plt.xlabel(r"$x_1$")
    plt.ylabel(r"$x_2$")
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.axline((-1, 1), (1, 0), color="r", linewidth=4.0)
    plt.contourf(x1, x2, z, 25)
    plt.plot(0.2, 0.4, "ro", markersize=12, markeredgecolor="black", markeredgewidth=3)
    ticks = np.linspace(0, 2, 5, endpoint=True)
    cbar_label = r"$x_1^2 + x_2^2$"
    cbar = plt.colorbar(ticks=ticks)
    cbar.set_label(cbar_label, rotation=270, labelpad=40)
    plt.show()


if __name__ == "__main__":
    main()
