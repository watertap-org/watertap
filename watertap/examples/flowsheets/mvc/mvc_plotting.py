import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import seaborn as sns
import numpy as np
import math

from pyomo.environ import (
    units as pyunits,
    check_optimal_termination,
    value,
    Expression,
    Param,
    Objective
)

# List of plotting functions
# make_cost_bar_charts - plots bar chart of cost breakdown


def main():

    cases = {}
    cases['evap_hx_cost'] = [(3000, 2000)]
    cases['elec_cost'] = [0.15]
    cases['cv_temp_max'] = [475]
    cases['comp_cost_factor'] = [1]

    dir = "C:/Users/carso/Documents/MVC/watertap_results/full_parameter_sweeps_distillate_hx_only"
    dir = "C:/Users/carso/Documents/MVC/watertap_results/full_parameter_sweeps"
    dir = "C:/Users/carso/Documents/MVC/watertap_results/full_parameter_sweeps_P_out_unfixed"
    map_dir = dir + "/evap_3000_hx_2000/elec_0.15/cv_temp_max_450/comp_cost_1"
    save_dir = map_dir + '/figures'

    make_maps(map_dir, save_dir)

    assert False
    filename_1 = dir + "/evap_2000_hx_2000/elec_0.15/cv_temp_max_450/comp_cost_1"
    filename_2 = dir + "/evap_2000_hx_2000/elec_0.15/cv_temp_max_450/comp_cost_2"
    filename_3 = dir + "/evap_4000_hx_4000/elec_0.3/cv_temp_max_500/comp_cost_1"
    filename_4 = dir + "/evap_4000_hx_4000/elec_0.3/cv_temp_max_500/comp_cost_2"
    wf = 0.1
    rr = 0.5
    cases_list = [(filename_1, wf, rr), (filename_2, wf, rr)]#, (filename_3, wf, rr), (filename_4, wf, rr)]
    make_cost_bar_charts(cases_list, ['2000\nCompressor Cost x1', '2000\nCompressor Cost x2'])#, '4000\nCompressor Cost x1', '4000\nCompressor Cost x2'])
    # dir_opt = "C:/Users/carso/Documents/MVC/watertap_results/opt_full_linearized/evap_6000_hx_3000/"
    # map_dir = "C:/Users/carso/Documents/MVC/watertap_results/opt_full_linearized/T_b_sensitivity/evap_6000_hx_3000/"
    # map_dir_list = [dir_opt,
    #                 map_dir+'T_b_80',
    #                 map_dir+'T_b_60',
    #                 map_dir+'T_b_40']
    # title_list = [r'$T_{b}$ optimized',
    #               r'$T_{b} = 80 C$',
    #               r'$T_{b} = 60 C$',
    #               r'$T_{b} = 40 C$']
    # save_dir = map_dir+"comparison_figures"
    # make_maps_comparison(map_dir_list, save_dir, title_list)
    return

def plot_2D_heat_map_dual_param(map_dir, save_dir, param, label, vmin, vmax, ticks, fmt, xticklabels, yticklabels, xlabel, ylabel, make_ticks=True, dif=False):
    fig = plt.figure()
    ax = plt.axes()
    results_file = map_dir + '/' + param + '.csv'
    df = pd.read_csv(results_file)
    cmap='YlGnBu'
    if dif == True:
        center_x = xticklabels.index("0")
        center_y = yticklabels.index("0")
        df_0 = df[str(center_x)][center_y]
        df = df.divide(df_0)-1
        df = df*100
        label = param + " change (%)"
        param = param + ' percentage difference'
        cmap = 'coolwarm'

    if make_ticks:
        decimal = int(fmt[1])
        df_min = float(np.nanmin(df.values))
        vmin = df_min - 10 ** -decimal
        df_max = float(np.nanmax(df.values))
        vmax = df_max + 10 ** -decimal
        n = 5  # number of ticks
        ticks = np.round(np.linspace(df_min, df_max, n), decimal)  # round based on formatting decimal places

    ax = sns.heatmap(df, cmap=cmap,
                     vmin=vmin, vmax=vmax, annot=True, annot_kws={"fontsize": 7}, fmt=fmt,
                     cbar_kws={'label': label, "ticks": ticks}, xticklabels=xticklabels,
                     yticklabels=yticklabels
                     )  # create heatmap
    ax.invert_yaxis()
    plt.xticks(fontsize=7.5)
    plt.yticks(fontsize=7.5,rotation=0)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.set_size_inches(3.25, 3.25)
    # plt.show()
    fig.savefig(save_dir + '/' + param + '.png', bbox_inches='tight', dpi=300)

def plot_2D_heat_map(map_dir, save_dir, param, label, vmin, vmax, ticks, fmt, make_ticks=True):
    fig = plt.figure()
    ax = plt.axes()
    results_file = map_dir + '/' + param + '.csv'
    df = pd.read_csv(results_file)
    # print(df)
    df = df.drop(df.columns[6],axis=1) # get rid of 175 g/kg case
    # print(df)
    # mask = df.isnull()
    df_wb = pd.read_csv("C:/Users/carso/Documents/MVC/watertap_results/Brine salinity.csv")
    df_wb = df_wb.drop(df_wb.columns[6],axis=1)
    mask = df_wb > 0.26
    if make_ticks:
        decimal = int(fmt[1])
        df_min = float(np.nanmin(df[~mask].values))
        vmin = df_min - 10 ** -decimal
        df_max = float(np.nanmax(df[~mask].values))
        vmax = df_max + 10 ** -decimal
        n = 5  # number of ticks
        ticks = np.round(np.linspace(df_min, df_max, n), decimal)  # round based on formatting decimal places

    xticklabels = ['25', '50', '75', '100', '125', '150']
    yticklabels = ['40', '45', '50', '55', '60', '65', '70', '75', '80']
    # ax = sns.heatmap(df, cmap='Reds', mask=mask,
    #                  vmin=vmin, vmax=vmax, annot=True, annot_kws={"fontsize": 8}, fmt=fmt,
    #                  cbar_kws={'label': label, "ticks": ticks}, xticklabels=xticklabels,
    #                  yticklabels=yticklabels
    #                  )  # create heatmap
    ax = sns.heatmap(df, cmap='YlGnBu', mask=mask,
                     vmin=vmin, vmax=vmax, annot=True, annot_kws={"fontsize": 8}, fmt=fmt,
                     cbar_kws={'label': label, "ticks": ticks}, xticklabels=xticklabels,
                     yticklabels=yticklabels
                     )  # create heatmap
    ax.invert_yaxis()
    plt.yticks(rotation=0)
    ax.set_xlabel('Feed concentration (g/kg)')
    ax.set_ylabel('Water recovery (%)')
    fig.set_size_inches(3.25, 3.25)
    # plt.show()
    fig.savefig(save_dir + '/' + param + '.png', bbox_inches='tight', dpi=300)

def plot_single_param_sensitivity(map_dir, x_param, y_param, x_label, y_label):
    results = map_dir + "/optimize_sweep.csv"
    df = pd.read_csv(results)
    x = df[x_param]
    y = df[y_param]
    fig = plt.figure()
    ax = plt.axes()
    plt.plot(x,y,'r')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    ax.set_xlabel(x_label,fontsize=12)
    ax.set_ylabel(y_label,fontsize=12)
    fig.set_size_inches(3.25, 3.25)
    fig.savefig(map_dir + '/' + y_param + ' vs ' + x_param + '.png', bbox_inches='tight', dpi=300)

def plot_2D_heat_map_subplots(map_dir, title_list, save_dir, param, param_label, vmin, vmax, ticks, fmt,
                              make_ticks=True, show=False):
    # map dir is a list
    n = len(map_dir)

    if make_ticks:
        decimal = int(fmt[1])
        # first file
        results_file = map_dir[0] + '/' + param + '.csv'
        df = pd.read_csv(results_file)
        df_min = float(np.nanmin(df.values))
        df_max = float(np.nanmax(df.values))

        # search for max and mins
        for i in range(1, n):
            results_file = map_dir[i] + '/' + param + '.csv'
            df = pd.read_csv(results_file)
            min_val = float(np.nanmin(df.values))
            if min_val < df_min:
                df_min = min_val
            max_val = float(np.nanmax(df.values))
            if max_val > df_max:
                df_max = max_val
        # make ticks
        vmin = df_min - 10 ** -decimal
        vmax = df_max + 10 ** -decimal
        n_ticks = 5  # number of ticks
        ticks = np.round(np.linspace(df_min, df_max, n_ticks), decimal)  # round based on formatting decimal places

    widths = [1 for i in range(n)]
    heights = [1 for i in range(n)]
    widths.append(0.08)
    fig, ax = plt.subplots(1, n + 1, gridspec_kw={'width_ratios': widths}, figsize=(n * 3, 3))  # figsize=(12,3))
    # ax[0].get_shared_y_axes().join(ax[1],ax[2],ax[3])
    # ax[1].get_shared_y_axes().join(ax[2],ax[3])

    cbar = [False for i in range(n - 1)]
    cbar.append(True)
    xticklabels = ['25', '50', '75', '100', '125', '150', '175']
    yticklabels = ['40', '45', '50', '55', '60', '65', '70', '75', '80']
    yticks = [[] for i in range(n)]
    yticks[0] = yticklabels
    g = {}
    for i in range(n):
        results_file = map_dir[i] + '/' + param + '.csv'
        df = pd.read_csv(results_file)
        mask = df.isnull()
        ax[i].axis('equal')
        g[i] = sns.heatmap(df, cmap='YlGnBu', mask=mask, square=True,
                           vmin=vmin, vmax=vmax, annot=True, annot_kws={"fontsize": 8}, fmt=fmt,
                           cbar_kws={'label': param_label, "ticks": ticks}, xticklabels=xticklabels,
                           yticklabels=yticklabels, ax=ax[i], cbar=cbar[i], cbar_ax=ax[i + 1]
                           )  # create heatmap

        # g[i].set_xlabel('Feed concentration (g/kg)')
        ax[i].set_yticklabels(yticklabels, rotation=0)
        ax[i].set_xticklabels(xticklabels, rotation=0)
        ax[i].set_title(title_list[i], size=8)
        g[i].invert_yaxis()
        g[i].set_xlabel('Feed concentration (g/kg)')

    g[0].set_ylabel('Water recovery (%)')
    if show:
        plt.show()
    fig.savefig(save_dir + '/' + param + '.png', bbox_inches='tight', dpi=300)

def make_maps(map_dir, save_dir):
    var = 'Distillate hx LMTD'
    label = r'Distillate HX LMTD (K)'
    vmin = 0
    vmax = 10
    ticks = [0, 10]
    fmt = '.0f'
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt, make_ticks=True)

    var = 'Brine hx LMTD'
    label = r'Brine HX LMTD (K)'
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt, make_ticks=True)

    var = 'Distillate hx heat transfer kW'
    label = r'Distillate HX Q [kW]'
    vmin = 0
    vmax = 10
    ticks = [0,10]
    fmt = '.0f'
    plot_2D_heat_map(map_dir, save_dir, var,label, vmin, vmax, ticks,fmt, make_ticks=True)

    var = 'Brine hx heat transfer kW'
    label = r'Brine HX Q [kW]'
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt, make_ticks=True)

    var = 'Normalized distillate hx heat transfer kJ per kg'
    label = r'Distillate HX normalized q [kJ/kg-feed]'
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt, make_ticks=True)

    # var = 'Normalized brine hx heat transfer kJ per kg'
    # label = r'Brine HX normalized q [kJ/kg-feed]'
    # plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt, make_ticks=True)

    var = 'LCOW'
    label = r'LCOW [\$/$\rmm^3$ of product]'
    vmin = 2.8  # minimum cost on bar, $/m3
    vmax = 6.1  # maximum cost on bar, $/m3
    ticks = [3, 4, 5, 6]  # tick marks on bar
    fmt = '.1f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'SEC'
    label = r'SEC [kWh/$\rmm^3$ of product]'
    vmin = 15  # minimum cost on bar, $/m3
    vmax = 63  # maximum cost on bar, $/m3
    ticks = [20, 30, 40, 50, 63]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Brine temperature Celsius'
    label = 'Evaporator temperature [C]'
    vmin = 25  # minimum cost on bar, $/m3
    vmax = 150  # maximum cost on bar, $/m3
    ticks = [25, 50, 100, 150]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Brine pressure kPa'
    label = 'Evaporator pressure [kPa]'
    vmin = 0  # minimum cost on bar, $/m3
    vmax = 403  # maximum cost on bar, $/m3
    ticks = [0, 100, 200, 300, 400]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Preheated feed temperature Celsius'
    label = 'Preheated feed temperature [C]'
    vmin = 25  # minimum cost on bar, $/m3
    vmax = 150  # maximum cost on bar, $/m3
    ticks = [50, 100, 150]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Compressed vapor temperature Celsius'
    label = 'Compressed vapor temperature [C]'
    vmin = 87  # minimum cost on bar, $/m3
    vmax = 227  # maximum cost on bar, $/m3
    ticks = [100, 150, 200]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Compressed vapor pressure kPa'
    label = 'Compressed vapor pressure [kPa]'
    vmin = 3  # minimum cost on bar, $/m3
    vmax = 629  # maximum cost on bar, $/m3
    ticks = [200, 400, 600]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Distillate temperature Celsius'
    label = 'Distillate temperature [C]'
    vmin = 34  # minimum cost on bar, $/m3
    vmax = 161  # maximum cost on bar, $/m3
    ticks = [50, 100, 150]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    # var = 'Preheater split ratio'
    # label = 'Preheater split ratio [C]'
    # vmin = 0.4  # minimum cost on bar, $/m3
    # vmax = 1  # maximum cost on bar, $/m3
    # ticks = [0.4, 0.6, 0.8]  # tick marks on bar
    # fmt = '.2f'  # format of annotation
    # plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Distillate hx area'
    label = r'Distillate preheater area [$\rmm^2$]'
    vmin = 0  # minimum cost on bar, $/m3
    vmax = 1630  # maximum cost on bar, $/m3
    ticks = [500, 1000, 2000, 3000, 4000]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    # var = 'Brine hx area'
    # label = r'Brine preheater area [$\rmm^2$]'
    # vmin = 0  # minimum cost on bar, $/m3
    # vmax = 1001  # maximum cost on bar, $/m3
    # ticks = [250, 500, 750]  # tick marks on bar
    # fmt = '.0f'  # format of annotation
    # plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Evaporator area'
    label = r'Evaporator area [$\rmm^2$]'
    vmin = 742  # minimum cost on bar, $/m3
    vmax = 3740  # maximum cost on bar, $/m3
    ticks = [1000, 2000, 3000]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Evaporator LMTD'
    label = r'Evaporator LMTD [K]'
    vmin = 19  # minimum cost on bar, $/m3
    vmax = 61  # maximum cost on bar, $/m3
    ticks = [20, 40, 60]  # tick marks on bar
    fmt = '.1f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Compressor pressure ratio'
    label = r'Compressor pressure ratio [-]'
    vmin = 1.5  # minimum cost on bar, $/m3
    vmax = 3.6  # maximum cost on bar, $/m3
    ticks = [2, 2.5, 3]  # tick marks on bar
    fmt = '.2f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Evaporator-feed temperature difference'
    label = r'$T_{brine}-T_{feed}$ [C]'
    vmin = -7  # minimum cost on bar, $/m3
    vmax = 15  # maximum cost on bar, $/m3
    ticks = [-5, 0, 5, 10, 15]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Mass flux LMH'
    label = r'Product flux over evaporator [LMH]'
    vmin = 28  # minimum cost on bar, $/m3
    vmax = 82  # maximum cost on bar, $/m3
    ticks = [30, 40, 50, 60, 70, 80]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

def make_maps_dual_param(map_dir, save_dir,xticklabels,yticklabels,xlabel,ylabel):
    var = 'LCOW'
    label = r'LCOW [\$/$\rmm^3$ of product]'
    vmin = 2.8  # minimum cost on bar, $/m3
    vmax = 6.1  # maximum cost on bar, $/m3
    ticks = [3, 4, 5, 6]  # tick marks on bar
    fmt = '.1f'  # format of annotation
    plot_2D_heat_map_dual_param(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt,xticklabels,yticklabels,xlabel,ylabel)

    plot_2D_heat_map_dual_param(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt, xticklabels, yticklabels,xlabel,ylabel, dif=True)

    var = 'SEC'
    label = r'SEC [kWh/$\rmm^3$ of product]'
    vmin = 15  # minimum cost on bar, $/m3
    vmax = 63  # maximum cost on bar, $/m3
    ticks = [20, 30, 40, 50, 63]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Brine temperature Celsius'
    label = 'Evaporator temperature [C]'
    vmin = 25  # minimum cost on bar, $/m3
    vmax = 150  # maximum cost on bar, $/m3
    ticks = [25, 50, 100, 150]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Brine pressure kPa'
    label = 'Evaporator pressure [kPa]'
    vmin = 0  # minimum cost on bar, $/m3
    vmax = 403  # maximum cost on bar, $/m3
    ticks = [0, 100, 200, 300, 400]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Preheated feed temperature Celsius'
    label = 'Preheated feed temperature [C]'
    vmin = 25  # minimum cost on bar, $/m3
    vmax = 150  # maximum cost on bar, $/m3
    ticks = [50, 100, 150]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Compressed vapor temperature Celsius'
    label = 'Compressed vapor temperature [C]'
    vmin = 87  # minimum cost on bar, $/m3
    vmax = 227  # maximum cost on bar, $/m3
    ticks = [100, 150, 200]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Compressed vapor pressure kPa'
    label = 'Compressed vapor pressure [kPa]'
    vmin = 3  # minimum cost on bar, $/m3
    vmax = 629  # maximum cost on bar, $/m3
    ticks = [200, 400, 600]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Distillate temperature Celsius'
    label = 'Distillate temperature [C]'
    vmin = 34  # minimum cost on bar, $/m3
    vmax = 161  # maximum cost on bar, $/m3
    ticks = [50, 100, 150]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Preheater split ratio'
    label = 'Preheater split ratio [C]'
    vmin = 0.4  # minimum cost on bar, $/m3
    vmax = 1  # maximum cost on bar, $/m3
    ticks = [0.4, 0.6, 0.8]  # tick marks on bar
    fmt = '.2f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Distillate hx area'
    label = r'Distillate preheater area [$\rmm^2$]'
    vmin = 0  # minimum cost on bar, $/m3
    vmax = 1630  # maximum cost on bar, $/m3
    ticks = [500, 1000, 2000, 3000, 4000]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Brine hx area'
    label = r'Brine preheater area [$\rmm^2$]'
    vmin = 0  # minimum cost on bar, $/m3
    vmax = 1001  # maximum cost on bar, $/m3
    ticks = [250, 500, 750]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Evaporator area'
    label = r'Evaporator area [$\rmm^2$]'
    vmin = 742  # minimum cost on bar, $/m3
    vmax = 3740  # maximum cost on bar, $/m3
    ticks = [1000, 2000, 3000]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Evaporator LMTD'
    label = r'Evaporator LMTD [K]'
    vmin = 19  # minimum cost on bar, $/m3
    vmax = 61  # maximum cost on bar, $/m3
    ticks = [20, 40, 60]  # tick marks on bar
    fmt = '.1f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Compressor pressure ratio'
    label = r'Compressor pressure ratio [-]'
    vmin = 1.5  # minimum cost on bar, $/m3
    vmax = 3.6  # maximum cost on bar, $/m3
    ticks = [2, 2.5, 3]  # tick marks on bar
    fmt = '.2f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Evaporator-feed temperature difference'
    label = r'$T_{brine}-T_{feed}$ [C]'
    vmin = -7  # minimum cost on bar, $/m3
    vmax = 15  # maximum cost on bar, $/m3
    ticks = [-5, 0, 5, 10, 15]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Mass flux LMH'
    label = r'Product flux over evaporator [LMH]'
    vmin = 28  # minimum cost on bar, $/m3
    vmax = 82  # maximum cost on bar, $/m3
    ticks = [30, 40, 50, 60, 70, 80]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

def make_maps_comparison(map_dir_list, save_dir, title_list):
    var = 'LCOW'
    label = r'LCOW [\$/$\rmm^3$ of product]'
    vmin = 2.8  # minimum cost on bar, $/m3
    vmax = 6.1  # maximum cost on bar, $/m3
    ticks = [3, 4, 5, 6]  # tick marks on bar
    fmt = '.1f'  # format of annotation
    plot_2D_heat_map_subplots(map_dir_list, title_list, save_dir, var, label, vmin, vmax, ticks, fmt,
                              make_ticks=True)

    var = 'SEC'
    label = r'SEC [kWh/$\rmm^3$ of product]'
    vmin = 15  # minimum cost on bar, $/m3
    vmax = 63  # maximum cost on bar, $/m3
    ticks = [20, 30, 40, 50, 63]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map_subplots(map_dir_list, title_list, save_dir, var, label, vmin, vmax, ticks, fmt,
                              make_ticks=True)

    var = 'Brine temperature Celsius'
    label = 'Evaporator temperature [C]'
    vmin = 25  # minimum cost on bar, $/m3
    vmax = 150  # maximum cost on bar, $/m3
    ticks = [25, 50, 100, 150]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map_subplots(map_dir_list, title_list, save_dir, var, label, vmin, vmax, ticks, fmt,
                              make_ticks=True)

    var = 'Brine pressure kPa'
    label = 'Evaporator pressure [kPa]'
    vmin = 0  # minimum cost on bar, $/m3
    vmax = 403  # maximum cost on bar, $/m3
    ticks = [0, 100, 200, 300, 400]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map_subplots(map_dir_list, title_list, save_dir, var, label, vmin, vmax, ticks, fmt,
                              make_ticks=True)

    var = 'Preheated feed temperature Celsius'
    label = 'Preheated feed temperature [C]'
    vmin = 25  # minimum cost on bar, $/m3
    vmax = 150  # maximum cost on bar, $/m3
    ticks = [50, 100, 150]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map_subplots(map_dir_list, title_list, save_dir, var, label, vmin, vmax, ticks, fmt,
                              make_ticks=True)

    var = 'Compressed vapor temperature Celsius'
    label = 'Compressed vapor temperature [C]'
    vmin = 87  # minimum cost on bar, $/m3
    vmax = 227  # maximum cost on bar, $/m3
    ticks = [100, 150, 200]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map_subplots(map_dir_list, title_list, save_dir, var, label, vmin, vmax, ticks, fmt,
                              make_ticks=True)

    var = 'Compressed vapor pressure kPa'
    label = 'Compressed vapor pressure [kPa]'
    vmin = 3  # minimum cost on bar, $/m3
    vmax = 629  # maximum cost on bar, $/m3
    ticks = [200, 400, 600]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map_subplots(map_dir_list, title_list, save_dir, var, label, vmin, vmax, ticks, fmt,
                              make_ticks=True)

    var = 'Distillate temperature Celsius'
    label = 'Distillate temperature [C]'
    vmin = 34  # minimum cost on bar, $/m3
    vmax = 161  # maximum cost on bar, $/m3
    ticks = [50, 100, 150]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map_subplots(map_dir_list, title_list, save_dir, var, label, vmin, vmax, ticks, fmt,
                              make_ticks=True)

    var = 'Preheater split ratio'
    label = 'Preheater split ratio [C]'
    vmin = 0.4  # minimum cost on bar, $/m3
    vmax = 1  # maximum cost on bar, $/m3
    ticks = [0.4, 0.6, 0.8]  # tick marks on bar
    fmt = '.2f'  # format of annotation
    plot_2D_heat_map_subplots(map_dir_list, title_list, save_dir, var, label, vmin, vmax, ticks, fmt,
                              make_ticks=True)

    var = 'Distillate hx area'
    label = r'Distillate preheater area [$\rmm^2$]'
    vmin = 0  # minimum cost on bar, $/m3
    vmax = 1630  # maximum cost on bar, $/m3
    ticks = [500, 1000, 2000, 3000, 4000]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map_subplots(map_dir_list, title_list, save_dir, var, label, vmin, vmax, ticks, fmt,
                              make_ticks=True)

    var = 'Brine hx area'
    label = r'Brine preheater area [$\rmm^2$]'
    vmin = 0  # minimum cost on bar, $/m3
    vmax = 1001  # maximum cost on bar, $/m3
    ticks = [250, 500, 750]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map_subplots(map_dir_list, title_list, save_dir, var, label, vmin, vmax, ticks, fmt,
                              make_ticks=True)

    var = 'Evaporator area'
    label = r'Evaporator area [$\rmm^2$]'
    vmin = 742  # minimum cost on bar, $/m3
    vmax = 3740  # maximum cost on bar, $/m3
    ticks = [1000, 2000, 3000]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map_subplots(map_dir_list, title_list, save_dir, var, label, vmin, vmax, ticks, fmt,
                              make_ticks=True)

    var = 'Evaporator LMTD'
    label = r'Evaporator LMTD [K]'
    vmin = 19  # minimum cost on bar, $/m3
    vmax = 61  # maximum cost on bar, $/m3
    ticks = [20, 40, 60]  # tick marks on bar
    fmt = '.1f'  # format of annotation
    plot_2D_heat_map_subplots(map_dir_list, title_list, save_dir, var, label, vmin, vmax, ticks, fmt,
                              make_ticks=True)

    var = 'Compressor pressure ratio'
    label = r'Compressor pressure ratio [-]'
    vmin = 1.5  # minimum cost on bar, $/m3
    vmax = 3.6  # maximum cost on bar, $/m3
    ticks = [2, 2.5, 3]  # tick marks on bar
    fmt = '.2f'  # format of annotation
    plot_2D_heat_map_subplots(map_dir_list, title_list, save_dir, var, label, vmin, vmax, ticks, fmt,
                              make_ticks=True)

    var = 'Evaporator-feed temperature difference'
    label = r'$T_{brine}-T_{feed}$ [C]'
    vmin = -7  # minimum cost on bar, $/m3
    vmax = 15  # maximum cost on bar, $/m3
    ticks = [-5, 0, 5, 10, 15]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map_subplots(map_dir_list, title_list, save_dir, var, label, vmin, vmax, ticks, fmt,
                              make_ticks=True)

    var = 'Mass flux LMH'
    label = r'Product flux over evaporator [LMH]'
    vmin = 28  # minimum cost on bar, $/m3
    vmax = 82  # maximum cost on bar, $/m3
    ticks = [30, 40, 50, 60, 70, 80]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map_subplots(map_dir_list, title_list, save_dir, var, label, vmin, vmax, ticks, fmt,
                              make_ticks=True)

def make_cost_bar_charts(cases_list, x):
    # cases_list: list [(filename, wf, rr)]
    # x: x axis labels
    n_cases = len(cases_list)
    pump_cc = []
    distillate_hx_cc = []
    brine_hx_cc = []
    mixer_cc = []
    evap_cc = []
    comp_cc = []
    mlc_oc = []
    elec_oc = []
    lcow = []
    for dir, wf, rr in cases_list:
        filename = dir + "/sweep_results.csv"
        df = pd.read_csv(filename, index_col=False)
        df.set_index(["# w_f", "recovery"], inplace=True, drop=True)

        # normalized LCOW costs
        pump_cc.append(df.loc[(wf, rr), 'LCOW normalized feed pump']
                       + df.loc[(wf, rr), 'LCOW normalized distillate pump']
                       + df.loc[(wf, rr), 'LCOW normalized brine pump'])
        distillate_hx_cc.append(df.loc[(wf, rr),'LCOW normalized distillate hx'])
        brine_hx_cc.append(df.loc[(wf, rr),'LCOW normalized brine hx'])
        mixer_cc.append(df.loc[(wf, rr),'LCOW normalized mixer'])
        evap_cc.append(df.loc[(wf, rr),'LCOW normalized evaporator'])
        comp_cc.append(df.loc[(wf, rr),'LCOW normalized compressor'])
        mlc_oc.append(df.loc[(wf, rr),'LCOW normalized MLC'])
        elec_oc.append(df.loc[(wf, rr),'LCOW normalized electricity'])
        lcow.append(df.loc[(wf, rr),'LCOW'])

    pump_cc = tuple(pump_cc)
    distillate_hx_cc = tuple(distillate_hx_cc)
    brine_hx_cc = tuple(brine_hx_cc)
    mixer_cc = tuple(mixer_cc)
    evap_cc = tuple(evap_cc)
    comp_cc = tuple(comp_cc)
    mlc_oc = np.array(mlc_oc)
    elec_oc = np.array(elec_oc)

    # Make plots
    width = 0.8
    ind = np.arange(n_cases)
    fig, ax = plt.subplots(figsize=(5, 4))
    p1 = plt.bar(x, elec_oc, width)
    p2 = plt.bar(x, mlc_oc, width, bottom=elec_oc)
    p3 = plt.bar(x, pump_cc, width, bottom=elec_oc + mlc_oc)
    p4 = plt.bar(x, distillate_hx_cc, width, bottom=elec_oc + mlc_oc)
    p5 = plt.bar(x, brine_hx_cc, width, bottom=elec_oc + mlc_oc + distillate_hx_cc)
    p6 = plt.bar(x, mixer_cc, width, bottom=elec_oc + mlc_oc + distillate_hx_cc + brine_hx_cc)
    p7 = plt.bar(x, comp_cc, width, bottom=elec_oc + mlc_oc + distillate_hx_cc + brine_hx_cc + mixer_cc)
    p8 = plt.bar(x, evap_cc, width, bottom=elec_oc + mlc_oc + distillate_hx_cc + brine_hx_cc + mixer_cc + comp_cc)

    for i, data in enumerate(lcow):
        plt.text(x=i, y=1 + 0.01, s=f"{round(data, 2)}" + r" $\$/m^3$", ha="center", fontsize=8)

    plt.ylabel('Normalized LCOW (-)')
    plt.legend((p8[0], p7[0], p6[0], p5[0], p4[0], p3[0], p2[0], p1[0]),
               ('Evaporator', 'Compressor', 'Mixer', 'Brine HX', 'Distillate HX', 'Pumps', 'MLC', 'Electricity'),
               loc='center left',
               bbox_to_anchor=(1, 0.5),
               prop={'size': 8}
               )
    plt.tight_layout()
    ax.xaxis.label.set_size(6)
    plt.show()

    return

def make_cost_bar_charts_old(map_dir):
    results_file = map_dir + '/costs_cases.csv'
    df = pd.read_csv(results_file, index_col=[0])  # to make the row names the keys
    n_cases = df.shape[1]
    x = []
    pump_cc = []
    distillate_hx_cc = []
    brine_hx_cc = []
    mixer_cc = []
    evap_cc = []
    comp_cc = []
    mlc_oc = []
    elec_oc = []
    lcow = []
    for i in range(n_cases):
        case_i = 'Case ' + str(i + 1)
        # label = case_i + '\n' + str(round(df[case_i]['Feed concentration'],0)) + ' g/kg\n' + str(round(df[case_i]['Recovery'],0)) + '%'
        label = str(round(df[case_i]['Feed concentration'])) + ' g/kg\n' + str(round(df[case_i]['Recovery'])) + '%'
        x.append(label)
        # normalized LCOW costs
        pump_cc.append(
            df[case_i]['LCOW normalized feed pump'] + df[case_i]['LCOW normalized distillate pump'] + df[case_i][
                'LCOW normalized brine pump'])
        distillate_hx_cc.append(df[case_i]['LCOW normalized distillate hx'])
        brine_hx_cc.append(df[case_i]['LCOW normalized brine hx'])
        mixer_cc.append(df[case_i]['LCOW normalized mixer'])
        evap_cc.append(df[case_i]['LCOW normalized evaporator'])
        comp_cc.append(df[case_i]['LCOW normalized compressor'])
        mlc_oc.append(df[case_i]['LCOW normalized MLC'])
        elec_oc.append(df[case_i]['LCOW normalized electricity'])
        lcow.append(df[case_i]['LCOW'])

    pump_cc = tuple(pump_cc)
    distillate_hx_cc = tuple(distillate_hx_cc)
    brine_hx_cc = tuple(brine_hx_cc)
    mixer_cc = tuple(mixer_cc)
    evap_cc = tuple(evap_cc)
    comp_cc = tuple(comp_cc)
    mlc_oc = np.array(mlc_oc)
    elec_oc = np.array(elec_oc)

    # Make plots
    width = 0.8
    ind = np.arange(n_cases)
    fig, ax = plt.subplots(figsize=(5, 4))
    p1 = plt.bar(x, elec_oc, width)
    p2 = plt.bar(x, mlc_oc, width, bottom=elec_oc)
    p3 = plt.bar(x, pump_cc, width, bottom=elec_oc + mlc_oc)
    p4 = plt.bar(x, distillate_hx_cc, width, bottom=elec_oc + mlc_oc)
    p5 = plt.bar(x, brine_hx_cc, width, bottom=elec_oc + mlc_oc + distillate_hx_cc)
    p6 = plt.bar(x, mixer_cc, width, bottom=elec_oc + mlc_oc + distillate_hx_cc + brine_hx_cc)
    p7 = plt.bar(x, comp_cc, width, bottom=elec_oc + mlc_oc + distillate_hx_cc + brine_hx_cc + mixer_cc)
    p8 = plt.bar(x, evap_cc, width, bottom=elec_oc + mlc_oc + distillate_hx_cc + brine_hx_cc + mixer_cc + comp_cc)

    for i, data in enumerate(lcow):
        plt.text(x=i, y=1 + 0.01, s=f"{round(data, 2)}" + r" $\$/m^3$", ha="center", fontsize=8)

    plt.ylabel('Normalized LCOW (-)')
    plt.legend((p8[0], p7[0], p6[0], p5[0], p4[0], p3[0], p2[0], p1[0]),
               ('Evaporator', 'Compressor', 'Mixer', 'Brine HX', 'Distillate HX', 'Pumps', 'MLC', 'Electricity'),
               loc='center left',
               bbox_to_anchor=(1, 0.5),
               prop={'size': 8}
               )
    plt.tight_layout()
    ax.xaxis.label.set_size(6)
    plt.show()

def plot_3D_results(results_file):
    df = pd.read_csv(results_file)
    outputs = list(df.columns)
    param_x_name = outputs[0]
    param_y_name = outputs[1]
    param_z_name = outputs[2]
    param_check = outputs[9]

    # get feasible results
    df_feasible = df[~df[param_check].isnull()]

    # make scatter plot
    plt.figure()
    axes=plt.axes(projection='3d')
    x = np.array(df_feasible[param_x_name])
    y = np.array(df_feasible[param_y_name])
    z = np.array(df_feasible[param_z_name])
    #color = np.array(df_feasible['Brine temperature']-273.15)
    color = np.array(df_feasible['LCOW'])
    color = np.array(df_feasible['Q external'])

    tick_min = min(color)
    tick_max = max(color)
    fig = axes.scatter3D(x,y,z,c=color) #, cbar_kws={'label': 'Brine temperature'})
    axes.set_xlabel('Evaporator area (m2)')
    axes.set_ylabel('Pressure ratio (-)')
    axes.set_zlabel('Vapor flow rate (kg/s)')
    #plt.colorbar(fig,label='LCOW ($/m3)')
    plt.colorbar(fig, label='Q external')
    plt.show()

if __name__ == "__main__":
    main()
