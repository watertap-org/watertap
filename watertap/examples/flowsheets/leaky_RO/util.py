from pyomo.environ import units as pyunits
from pyomo.environ import value, units

from idaes.core.util.scaling import *
from watertap.core.util.infeasible import *

import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors as colors
from matplotlib import ticker
import matplotlib.ticker as mtick

liq = "Liq"
h2o = "H2O"
tds = 'TDS'


def print_results(m):
    print('\n')
    print('=============== SIMULATION RESULTS ===============\n')
    print(
            f'{"Water Permeability Coeff":<30s}{f"{value(m.fs.A):<10,.2e}"}{"m/s*kPa":<10s}'
        )
    print('\n')
    print(
            f'{"Feed Pressure":<30s}{f"{value(units.convert(m.fs.feed[0].pressure() * pyunits.Pa, to_units=pyunits.bar)):<10,.2f}"}{"bar":<10s}'
        )
    print(
            f'{"Feed Mass Flow Rate":<30s}{f"{value(m.fs.feed_mass_flow):<10,.2f}"}{str(pyunits.get_units(m.fs.feed_mass_flow)):<10s}'
        )
    print(
            f'{"Feed Solute Mass Flow Rate":<30s}{f"{value(units.convert(m.fs.feed[0].flow_mass_phase_comp[liq, tds], to_units=pyunits.g / pyunits.second)):<10,.2f}"}{"g/s":<10s}'
        )
    print(
            f'{"Feed Vol Flow Rate":<30s}{f"{value(pyunits.convert(m.fs.feed_vol_flow, to_units=pyunits.L * pyunits.second ** -1)):<10,.2f}"}{"L/s":<10s}'
        )
    print(
            f'{"Feed Osm Pressure":<30s}{f"{value(units.convert(m.fs.feed[0].pressure_osm_phase[liq], to_units=pyunits.bar)):<10,.2f}"}{"bar":<10s}'
        )
    print('\n')
    print(
            f'{"Perm Pressure":<30s}{f"{value(units.convert(m.fs.perm[0].pressure() * pyunits.Pa, to_units=pyunits.bar)):<10,.2f}"}{"bar":<10s}'
        )
    print(
            f'{"Perm Mass Flow Rate":<30s}{f"{value(units.convert(m.fs.perm_mass_flow, to_units=pyunits.kg / pyunits.second)):<10,.2f}"}{"kg/s":<10s}'
        )
    print(
            f'{"Perm Solute Mass Flow Rate":<30s}{f"{value(units.convert(m.fs.perm[0].flow_mass_phase_comp[liq, tds], to_units=pyunits.g / pyunits.second)):<10,.2f}"}{"g/s":<10s}'
        )
    print(
            f'{"Perm Vol Flow Rate":<30s}{f"{value(units.convert(m.fs.perm_vol_flow, to_units=pyunits.L / pyunits.second)):<10,.2f}"}{"L/s":<10s}'
        )
    print(
            f'{"Perm Osm Pressure":<30s}{f"{value(units.convert(m.fs.perm[0].pressure_osm_phase[liq], to_units=pyunits.bar)):<10,.2f}"}{"bar":<10s}'
        )
    print('\n')
    print(
            f'{"Perm Density":<30s}{f"{value(m.fs.perm[0].dens_mass_solvent()):<10,.2f}"}{str(pyunits.get_units(m.fs.perm[0].dens_mass_solvent)):<10s}'
        )
    print(
            f'{"Rejection":<30s}{f"{value(m.fs.rejection):<10,.3f}"}{"%":<10s}'
        )
    print('\n')
    print(
            f'{"Water Flux":<30s}{f"{value(units.convert(m.fs.flux_mass_phase_comp / m.fs.perm[0].dens_mass_solvent, to_units=pyunits.liter * pyunits.m ** -2 * pyunits.hour ** -1)):<10,.2f}"}{"LMH":<10s}'
        )
    print(
            f'{"Salt Flux":<30s}{f"{value(units.convert(m.fs.salt_flux_mass_phase_comp, to_units=pyunits.kg * pyunits.m ** -2 * pyunits.hour ** -1)):<10,.2f}"}{"GMH":<10s}'
        )
    print('\n')
    print(
            f'{"Leaky Water Flux":<30s}{f"{value(units.convert(m.fs.leaky_flux_mass_phase_comp / m.fs.perm[0].dens_mass_solvent, to_units=pyunits.liter * pyunits.m ** -2 * pyunits.hour ** -1)):<10,.2f}"}{"LMH":<10s}'
        )
    print(
            f'{"Leaky Salt Flux":<30s}{f"{value(units.convert(m.fs.leaky_salt_flux_mass_phase_comp, to_units=pyunits.kg * pyunits.m ** -2 * pyunits.hour ** -1)):<10,.2f}"}{"GMH":<10s}'
        )
    
    print('\n')
    print(
            f'{"Salt Perm":<30s}{f"{value(m.fs.B):<10,.2e}"}{str(pyunits.get_units(m.fs.B)):<10s}'
        )


    print('\n')


def create_report(m):
    var_labels = ["Water Permeability Coeff", "Feed Pressure", "Feed Mass Flow Rate", "Feed Solute Mass Flow Rate", "Feed Vol Flow Rate", "Feed Osm Pressure", "Perm Pressure", "Perm Mass Flow Rate", "Perm Solute Mass Flow Rate", "Perm Vol Flow Rate", "Perm Osm Pressure", "Perm Density", "Rejection", "Water Flux", "Salt Flux", "Leaky Water Flux", "Leaky Salt Flux", "Obs Salt Perm", "Reflection Coeff.", "Alpha"]
    vars = [m.fs.A, m.fs.feed[0].pressure(), m.fs.feed_mass_flow, m.fs.feed[0].flow_mass_phase_comp[liq, tds], m.fs.feed_vol_flow, m.fs.feed[0].pressure_osm_phase[liq], m.fs.perm[0].pressure(), m.fs.perm_mass_flow, m.fs.perm[0].flow_mass_phase_comp[liq, tds], m.fs.perm_vol_flow, m.fs.perm[0].pressure_osm_phase[liq], m.fs.perm[0].dens_mass_solvent(), m.fs.rejection, m.fs.flux_mass_phase_comp, m.fs.salt_flux_mass_phase_comp, m.fs.leaky_flux_mass_phase_comp, m.fs.leaky_salt_flux_mass_phase_comp, m.fs.B, m.fs.reflect_coeff, m.fs.alpha]
    var_vals = [value(i) for i in vars]
    df = pd.DataFrame.from_dict(dict.fromkeys(var_labels, []))
    df = pd.concat([df,pd.DataFrame([var_vals], columns=var_labels)], ignore_index=True)

    return df

def SKK_RO_report(blk):
    liq = "Liq"
    h2o = "H2O"
    nacl = "NaCl"
    if blk.config.transport_model == 'SKK':
        row = [
            blk.name,
            value(pyunits.convert(blk.inlet.pressure[0], to_units=pyunits.bar)),
            blk.area.value,
            blk.A_comp[0, h2o].value,# * (3.6e11),
            blk.B_comp[0, nacl].value,# * (1000.0 * 3600.0),
            blk.reflect_coeff.value,
            blk.alpha.value,
            value(units.convert(blk.flux_mass_phase_comp_avg[0, liq, h2o] / (1000 * units.kg / units.m ** 3),
                              to_units=units.mm / units.hr)),
            value(units.convert(blk.flux_mass_phase_comp_avg[0, liq, nacl], 
                                to_units=units.kg / units.m ** 2 / units.hr)),
            100 * value(blk.recovery_vol_phase[0.0, liq]),
            100* (1 - value(blk.mixed_permeate[0.0].flow_mass_phase_comp['Liq', 'NaCl'] / 
                            blk.feed_side.properties_in[0.0].flow_mass_phase_comp['Liq', 'NaCl']
                       )
                       ),
            100*value(blk.rejection_phase_comp[0, liq, nacl])
        ]
    else:
        row = [
            blk.name,
            value(pyunits.convert(blk.inlet.pressure[0], to_units=pyunits.bar)),
            blk.area.value,
            blk.A_comp[0, h2o].value,# * (3.6e11),
            blk.B_comp[0, nacl].value,# * (1000.0 * 3600.0),
            1,
            float("NaN"),
            value(units.convert(blk.flux_mass_phase_comp_avg[0, liq, h2o] / (1000 * units.kg / units.m ** 3),
                              to_units=units.mm / units.hr)),
            value(units.convert(blk.flux_mass_phase_comp_avg[0, liq, nacl], 
                                to_units=units.kg / units.m ** 2 / units.hr)),
            100 * value(blk.recovery_vol_phase[0.0, liq]),
            100* (1 - value(blk.mixed_permeate[0.0].flow_mass_phase_comp['Liq', 'NaCl'] / 
                            blk.feed_side.properties_in[0.0].flow_mass_phase_comp['Liq', 'NaCl'])),
            100* value(blk.rejection_phase_comp[0, liq, nacl])
            
        ]
    df = pd.DataFrame(
                [row],
                columns=[
                    "Unit",
                    "Pressure [bar]",
                    "Membrane Area [m2]",
                    "Water Perm. [LMH/bar]",
                    "Salt Perm. [LMH]",
                    "Reflection Coeff.",
                    "Alpha",
                    "Water Flux [LMH]",
                    "Salt Flux [kg/m2/hr]",
                    "Recovery [%]",
                    "Rejection [%]",
                    "Rejection 2 [%]"
                ],
            )
    return df

def auto_scale(m):
     debug(m, automate_rescale=True, verbose=False)

def debug(m, automate_rescale=False, verbose=False):
    badly_scaled_vars = list(badly_scaled_var_generator(m))
    if verbose == True:
        print(f'\n{"=======> DEBUGGING <=======":^60}\n')
        print(f'\n{"=======> BADLY SCALED VARIABLES <=======":^60}\n')
        print([print(i[0], i[1]) for i in badly_scaled_vars])
        print(f'\n{"=======> INFEASIBLE BOUNDS <=======":^60}\n')
        print_infeasible_bounds(m)
        print(f'\n{"=======> INFEASIBLE CONSTRAINTS <=======":^60}\n')
        print_infeasible_constraints(m)
        print(f'\n{"=======> CONSTRAINTS CLOSE TO BOUNDS <=======":^60}\n')
        print_close_to_bounds(m)

    if automate_rescale:
        if verbose == True:
            print(
                f"\n{len(badly_scaled_vars)} poorly scaled "
                f"variable(s) will be rescaled so that each scaled variable value = 1\n"
            )
        automate_rescale_variables(m)
        badly_scaled_vars = list(badly_scaled_var_generator(m))
    if verbose == True:
        print(
                f"\nNow {len(badly_scaled_vars)} poorly scaled\n"
            )
    
def automate_rescale_variables(self, rescale_factor=1, default=1):
        if rescale_factor is None:
            rescale_factor = 1
        for var, sv in badly_scaled_var_generator(self):
            sf = get_scaling_factor(var)
            if get_scaling_factor(var) is None:
                print(f"{var} is missing a scaling factor")
                sf = default
                set_scaling_factor(var, sf, data_objects=False)

            set_scaling_factor(var, sf / sv * rescale_factor)
            calculate_scaling_factors(self)

def plot_data(x, y1, y2, y3, ax=None, xlabel=None, ylabel=None, ylabel2=None, label=None):
    colors = ["#2f5c8c", "#c07432", "#474747"]
    ax.plot(x,y1,'--',linewidth=2, c=colors[0])
    ax.plot(x,y1,'o',label=ylabel,linewidth=0, mfc=colors[0], mec=colors[0])
    # ax.legend(loc='upper left', fontsize=12, frameon = False)
    ax2 = ax.twinx()
    ax2.plot(x,y2,'--',linewidth=2, c=colors[1])
    ax2.plot(x,y2,'o',label=ylabel2,linewidth=0, mfc=colors[1], mec=colors[1])
    # ax2.legend(loc='upper right', fontsize=12, frameon = False)
    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)
    ax.tick_params(axis="x", labelsize=16)
    ax.tick_params(axis="y", labelsize=16)
    ax.yaxis.set_major_locator(plt.MaxNLocator(6))
    ax.yaxis.set_minor_locator(plt.MaxNLocator(6))
    ax.tick_params(axis='y', which='minor', bottom=False)


    ax2.set_ylabel(ylabel2, fontsize=16)
    ax2.tick_params(axis="y", labelsize=16)
    ax.tick_params(axis="y", colors=colors[0])
    ax2.tick_params(axis="y", colors=colors[1])
    ax.yaxis.label.set_color(colors[0])
    ax2.yaxis.label.set_color(colors[1])
    # Set number of tick marks for y-axis
    ax2.yaxis.set_major_locator(plt.MaxNLocator(5))

    ax3 = ax.twinx()
    ax3.spines.right.set_position(("outward", 85))
    ax3.plot(x,y3,'--',linewidth=2, c='k')
    ax3.plot(x,y3,'o',label='Rejection',linewidth=0, mfc="k", mec="k")
    ax3.set_ylabel('Rejection (%)', fontsize=16)
    ax3.tick_params(axis="y", labelsize=16)
    ax3.yaxis.label.set_color('k')
    ax3.set_ylim(90, 100)

    plt.figlegend(loc='upper left', fontsize=12, frameon = False, bbox_to_anchor=(0.1, 0.98), ncols=2)
    # ax.set_xlim(min(x), max(x))
    # ax2.set_xlim(min(x), max(x))

    ax.set_ylim(0, 60)
    ax2.set_ylim(0, 0.1)
    

    annote_lat = 0.28
    annote_long = 0.03
    annote_long2 = 0.65
    
    ax.annotate('',
            xy=(0.95, 0.03), xycoords='axes fraction',
            xytext=(0.65, 0.03), textcoords='axes fraction',
            size=24,
            arrowprops=dict(arrowstyle="simple",
                            fc="0.6", ec="none",
                            patchB=None,
                            connectionstyle=None))
    
    ax.annotate('',
            xy=(0.05, 0.03), xycoords='axes fraction',
            xytext=(0.35, 0.03), textcoords='axes fraction',
            size=24,
            arrowprops=dict(arrowstyle="simple",
                            fc="0.6", ec="none",
                            patchB=None,
                            connectionstyle=None))
    
    # ax.annotate('Leaky\nMembrane',
    ax.annotate('Pristine\nMembrane',
            xy=(0.80, 0.08), xycoords='axes fraction',
            horizontalalignment='center', verticalalignment='center',
            fontsize=12)
    
    # ax.annotate('Pristine\nMembrane',
    ax.annotate('Leaky\nMembrane',
            xy=(0.20, 0.08), xycoords='axes fraction',
            horizontalalignment='center', verticalalignment='center',
            fontsize=12)
        
    ax.annotate(r'$J_w = A(\Delta P\ - \sigma\Delta\pi)$',
            xy=(annote_long, annote_lat), xycoords='axes fraction',
            horizontalalignment='left', verticalalignment='top',
            bbox=dict(boxstyle="round", fc="0.8"),
            fontsize=14)
    
    ax.annotate(r'$J_s = B\Delta C + (1-\sigma)J_wC_f$',
            xy=(annote_long, annote_lat-0.09), xycoords='axes fraction',
            horizontalalignment='left', verticalalignment='top',
            bbox=dict(boxstyle="round", fc="0.8"),
            fontsize=14)
    
    ax.annotate(r'$A = 4.2\bullet10^{-12}  \frac{m}{s\bullet Pa}$',
            xy=(annote_long2, annote_lat-0.03), xycoords='axes fraction',
            horizontalalignment='left', verticalalignment='top',
            fontsize=12)
    
    ax.annotate(r'$B = 3.5\bullet10^{-7} \frac{m}{s}$',
            xy=(annote_long2, annote_lat-0.1), xycoords='axes fraction',
            horizontalalignment='left', verticalalignment='top',
            fontsize=12)
    

def plot_fixed_alpha(x, y1, y2, y3, y4, ax=None, xlabel=None, ylabel=None, ylabel2=None):
    colors = ["#2f5c8c", "#c07432", "#474747", "#26425A"]
    ax.plot(x,y1,'--',linewidth=2, c=colors[0])
    ax.plot(x,y1,'o',label=ylabel,linewidth=0, mfc=colors[0], mec=colors[0])
    # ax.legend(loc='upper left', fontsize=12, frameon = False)
    ax2 = ax.twinx()
    ax2.plot(x,y2,'--',linewidth=2, c=colors[1])
    ax2.plot(x,y2,'o',label=ylabel2,linewidth=0, mfc=colors[1], mec=colors[1])
    # ax2.legend(loc='upper right', fontsize=12, frameon = False)
    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)
    ax.tick_params(axis="x", labelsize=16)
    ax.tick_params(axis="y", labelsize=16)
    ax.yaxis.set_major_locator(plt.MaxNLocator(6))
    ax.yaxis.set_minor_locator(plt.MaxNLocator(6))
    ax.tick_params(axis='y', which='minor', bottom=False)

    ax.set_ylim(0, 1E-6)

    ax2.set_ylabel(ylabel2, fontsize=16)
    ax2.tick_params(axis="y", labelsize=16)
    ax.tick_params(axis="y", colors=colors[0])
    ax2.tick_params(axis="y", colors=colors[1])
    ax.yaxis.label.set_color(colors[0])
    ax2.yaxis.label.set_color(colors[1])
    # Set number of tick marks for y-axis
    ax2.yaxis.set_major_locator(plt.MaxNLocator(5))
    ax2.set_ylim(0, .1)
    plt.gca().yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0, decimals=0))

    ax3 = ax.twinx()
    ax3.plot(x,y3,'s:',linewidth=2, c='k',label='SD Rejection')
    # ax3.plot(x,y3,'o',label='SD Rejection',linewidth=0, mfc="k", mec="k")
    ax3.plot(x,y4,'v:',linewidth=2, c='k',label='SKK Rejection')
    # ax3.plot(x,y4,'o',label='SKK Rejection',linewidth=0, mfc="k", mec="k")
    ax3.set_ylabel('Rejection (%)', fontsize=16)
    ax3.tick_params(axis="y", labelsize=16)
    ax3.yaxis.label.set_color('k')
    ax3.set_ylim(.90, 1)
    plt.gca().yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0, decimals=0))


    ax2.spines.right.set_position(("outward", 85))

    plt.figlegend(loc='upper left', fontsize=10, frameon = False, bbox_to_anchor=(0.12, 0.95), ncols=1)

def plot_contour(x, y, z, levels=None, x_label='', y_label='', z_label='', low=-1, mid=0, high=1):

    cmap_pallete="RdBu_r"
    divnorm = colors.TwoSlopeNorm(vmin=0, vcenter=5, vmax=18)
    divnorm = colors.SymLogNorm(linthresh=1, linscale=0.1, vmin=0, vmax=15)

    fig, ax = plt.subplots(figsize=(6,4.8))
    if levels != None:
        cs2 = ax.contourf(x, y, z, 100, cmap=cmap_pallete, norm=divnorm)
        cs3 = ax.contour(x, y, z, levels, colors='k', linewidths=2, linestyles='dashed', norm=divnorm)
        ax.clabel(cs3, fmt='%1.1f%%', colors='k', fontsize=12, inline_spacing=40)
        cbar = fig.colorbar(cs2)
    else:
        cs1 = ax.contourf(x, y, z, 100, cmap=cmap_pallete, norm=divnorm)

    plt.gca().set_aspect('equal')
    ax.set_aspect('equal')
    tick_locator = ticker.MaxNLocator(nbins=6)
    cbar.locator = tick_locator
    cbar.update_ticks()

    for l in cbar.ax.yaxis.get_ticklabels():
        l.set_fontsize(16)

    cbar.ax.yaxis.set_major_formatter('{x:1.0f}%')
    cbar.ax.get_yaxis().labelpad = 10
    cbar.ax.set_ylabel(z_label, rotation=90, fontsize=16)

    ax.set_xlabel(x_label, fontsize=16)
    ax.set_ylabel(y_label, fontsize=16)
    ax.tick_params(axis="x", labelsize=16)
    ax.tick_params(axis="y", labelsize=16)
    ax.yaxis.set_major_locator(plt.MaxNLocator(6))
    ax.yaxis.set_minor_locator(plt.MaxNLocator(6))
    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax.tick_params(axis='y', which='minor', bottom=False)

    

    # ax.set_xlim(0, 1.5e-12)
    # ax.set_ylim(0, 1e-7)
    # ax.set_xscale('log')
    ax.set_yscale('log')
    # ax.set_aspect('equal')
    # ax.grid(True, which='both', axis='both', linestyle='--')
    fig.tight_layout()

    return fig, ax