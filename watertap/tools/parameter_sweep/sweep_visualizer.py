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
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def line_plot(data_path, xlabel, ylabel, xunit=None, yunit=None, filetype=".pdf"):
    """
    Plots a simple line and returns the figure object
    Args:
        data_path : path to raw data csv
        xlabel : title of the independent variable column
        ylabel : title of the dependent variable column
        filetype : filetype of the save file (pdf, csv, etc.)
    Returns:
        fig : figure object
        ax : axes object
    """
    # read csv data into dataframe
    df = pd.read_csv(data_path)

    # create figure object
    plt.rcParams.update({"font.size": 10})
    fig, ax = plt.subplots(figsize=(5.1667, 5.1667))  # default size for SIAM
    ax.plot(df[xlabel].values, df[ylabel].values)

    # Add axes labels
    if xunit is not None:
        xlabel += " [" + xunit + "]"
    ax.set_xlabel(xlabel)

    if yunit is not None:
        ylabel += " [" + yunit + "]"
    ax.set_ylabel(ylabel)

    # save file
    plt.savefig(data_path.split(".")[0] + filetype, dpi=300)

    # return figure and axes objects
    return fig, ax


def contour_plot(
    data_path,
    xlabel,
    ylabel,
    zlabel,
    xunit=None,
    yunit=None,
    zunit=None,
    levels=10,
    cmap="viridis",
    isolines=None,
    filetype=".pdf",
):
    """
    Plots a filled contour and returns the figure object
    Args:
        data_path : path to raw data csv
        xlabel : title of the independent variable column
        ylabel : title of the dependent variable column
        zlabel : title of the dependent variable column
        levels : number of iso-lines
        cmap : color scheme defined by matplotlib at https://matplotlib.org/3.5.0/tutorials/colors/colormaps.html
        filetype : filetype of the save file (pdf, csv, etc.)
    Returns:
        fig : figure object
        ax : axes object
    """
    # read csv data into dataframe
    df = pd.read_csv(data_path)

    # reshape data
    Z = df.pivot_table(index=xlabel, columns=ylabel, values=zlabel).T.values
    X_unique = np.sort(df[xlabel].unique())
    Y_unique = np.sort(df[ylabel].unique())
    X, Y = np.meshgrid(X_unique, Y_unique)

    # create plot
    plt.rcParams.update({"font.size": 10})
    fig, ax = plt.subplots(figsize=(5.1667, 5.1667))  # default size for SIAM
    cs = ax.contourf(X, Y, Z, levels=levels, cmap=cmap)
    cbar = fig.colorbar(cs)

    # add isolines as needed
    if isolines is not None:
        cs_lvl = ax.contour(cs, levels=isolines, colors="k", linewidths=1)
        ax.clabel(cs_lvl, inline=True, fontsize=8)
        cbar.add_lines(cs_lvl)

    # Add axes labels
    if xunit is not None:
        xlabel += " [" + xunit + "]"
    ax.set_xlabel(xlabel)

    if yunit is not None:
        ylabel += " [" + yunit + "]"
    ax.set_ylabel(ylabel)

    if zunit is not None:
        zlabel += " [" + zunit + "]"
    ax.set_title(zlabel)

    plt.savefig(data_path.split(".")[0] + filetype, dpi=300)

    return fig, ax


# def box_and_whisker_plot()
