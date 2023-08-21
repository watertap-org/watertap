from pyomo.environ import units as pyunits
from pyomo.environ import value, units

from idaes.core.util.scaling import *
from watertap.core.util.infeasible import *

import matplotlib.pyplot as plt
import pandas as pd

liq = "Liq"
h2o = "H2O"
tds = 'TDS'

def print_results(m):
    print('\n')
    print('=============== SIMULATION RESULTS ===============\n')
    print(
            f'{"Water Permeability Coeff":<30s}{f"{value(m.fs.A_prime):<10,.2e}"}{"m/s*kPa":<10s}'
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
            f'{"Obs Salt Perm":<30s}{f"{value(m.fs.obs_salt_perm):<10,.2f}"}{str(pyunits.get_units(m.fs.obs_salt_perm)):<10s}'
        )


    print('\n')


def create_report(m):
    var_labels = ["Water Permeability Coeff", "Feed Pressure", "Feed Mass Flow Rate", "Feed Solute Mass Flow Rate", "Feed Vol Flow Rate", "Feed Osm Pressure", "Perm Pressure", "Perm Mass Flow Rate", "Perm Solute Mass Flow Rate", "Perm Vol Flow Rate", "Perm Osm Pressure", "Perm Density", "Rejection", "Water Flux", "Salt Flux", "Leaky Water Flux", "Leaky Salt Flux", "Obs Salt Perm", "Reflection Coeff.", "Alpha"]
    vars = [m.fs.A_prime, m.fs.feed[0].pressure(), m.fs.feed_mass_flow, m.fs.feed[0].flow_mass_phase_comp[liq, tds], m.fs.feed_vol_flow, m.fs.feed[0].pressure_osm_phase[liq], m.fs.perm[0].pressure(), m.fs.perm_mass_flow, m.fs.perm[0].flow_mass_phase_comp[liq, tds], m.fs.perm_vol_flow, m.fs.perm[0].pressure_osm_phase[liq], m.fs.perm[0].dens_mass_solvent(), m.fs.rejection, m.fs.flux_mass_phase_comp, m.fs.salt_flux_mass_phase_comp, m.fs.leaky_flux_mass_phase_comp, m.fs.leaky_salt_flux_mass_phase_comp, m.fs.obs_salt_perm, m.fs.reflect_coeff, m.fs.alpha]
    var_vals = [value(i) for i in vars]
    df = pd.DataFrame.from_dict(dict.fromkeys(var_labels, []))
    df = pd.concat([df,pd.DataFrame([var_vals], columns=var_labels)], ignore_index=True)

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

def plot_data(x, y1, y2, y3, ax, xlabel=None, ylabel=None, ylabel2=None, label=None):
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
    ax2.set_ylabel(ylabel2, fontsize=16)
    ax2.tick_params(axis="y", labelsize=16)
    ax.tick_params(axis="y", colors=colors[0])
    ax2.tick_params(axis="y", colors=colors[1])
    ax.yaxis.label.set_color(colors[0])
    ax2.yaxis.label.set_color(colors[1])

    ax3 = ax.twinx()
    ax3.spines.right.set_position(("outward", 80))
    ax3.plot(x,y3,'--',linewidth=2, c='k')
    ax3.plot(x,y3,'o',label='Rejection',linewidth=0, mfc="k", mec="k")
    ax3.set_ylabel('Rejection (%)', fontsize=16)
    ax3.tick_params(axis="y", labelsize=16)
    ax3.yaxis.label.set_color('k')

    plt.figlegend(loc='upper left', fontsize=12, frameon = False, bbox_to_anchor=(0.1, 0.98), ncols=2)
    # ax.set_xlim(min(x), max(x))
    # ax2.set_xlim(min(x), max(x))

    ax.set_ylim(0, 30)
    ax2.set_ylim(0, 1)
    ax3.set_ylim(99, 100)

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
    
    ax.annotate('Leaky\nMembrane',
            xy=(0.80, 0.08), xycoords='axes fraction',
            horizontalalignment='center', verticalalignment='center',
            fontsize=12)
    
    ax.annotate('Pristine\nMembrane',
            xy=(0.20, 0.08), xycoords='axes fraction',
            horizontalalignment='center', verticalalignment='center',
            fontsize=12)
        
    ax.annotate(r'$J_v = L_p(\Delta P\ - \sigma\Delta\pi)$',
            xy=(annote_long, annote_lat), xycoords='axes fraction',
            horizontalalignment='left', verticalalignment='top',
            bbox=dict(boxstyle="round", fc="0.8"),
            fontsize=14)
    
    ax.annotate(r'$J_s = \omega^\prime\Delta\pi + (1-\sigma)J_v\mathcal{C}$',
            xy=(annote_long, annote_lat-0.09), xycoords='axes fraction',
            horizontalalignment='left', verticalalignment='top',
            bbox=dict(boxstyle="round", fc="0.8"),
            fontsize=14)
    
    ax.annotate(r'$L_p = 1\bullet10^{-9}  \frac{m}{s\bullet kPa}$',
            xy=(annote_long2, annote_lat), xycoords='axes fraction',
            horizontalalignment='left', verticalalignment='top',
            fontsize=12)
    
    ax.annotate(r'$\omega^{\prime} = 5\bullet10^{-11} \frac{m}{s\bullet kPa}$',
            xy=(annote_long2, annote_lat-0.07), xycoords='axes fraction',
            horizontalalignment='left', verticalalignment='top',
            fontsize=12)