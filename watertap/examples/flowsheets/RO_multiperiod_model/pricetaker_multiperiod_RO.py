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
import os
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rc("font", size=22)
plt.rc("axes", titlesize=22)
scaling_obj = 1
scaling_factor = 1

from pyomo.environ import (
    Param,
    Expression,
    Objective,
    Constraint,
    units as pyunits,
    Block,
    Expr_if,
    Var,
    value,
    ComponentMap,
)
from idaes.core.solvers import get_solver
from idaes.core.util import model_statistics as stats
from watertap.examples.flowsheets.RO_multiperiod_model.multiperiod_RO import (
    create_multiperiod_swro_model,
)
import watertap.examples.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery as swro
from watertap.unit_models.pressure_changer import VariableEfficiency


def main(
    n_steps=24,
    filename="pricesignals_GOLETA_6_N200_20220601.csv",
    cost_of_carbon = 0,
    daily_production = 50,
    variable_efficiency = VariableEfficiency.none,
):
    file_path = os.path.realpath(__file__)
    base_path = os.path.dirname(file_path)
    data_path = os.path.join(base_path, filename)

    # get data
    price_signal, lmp, co2i = _get_lmp(
        time_steps=n_steps,
        data_path=data_path,
        cost_of_carbon = cost_of_carbon,
    )

    mp_swro = build_flowsheet(n_steps,
                              variable_efficiency=variable_efficiency)

    m, t_blocks = set_objective(mp_swro,
                                price_signal,
                                lmp, 
                                co2i,
                                product_requirement=daily_production,)

    m, results = solve(m)

    df = save_results(m)

    return m, df


def _get_lmp(time_steps, 
             data_path, 
             cost_of_carbon):
    """
    Get price signals from data set
        time_steps: Number of time steps considered in MP analysis
        data_path: Price [$/kWh] and co2i [kg/kWh] on the same interval as time_steps, 
        cost-of-carbon: represents the cost of carbon in $/kg
    :return: 
        weighted price signal
        reshaped electricity price
        reshaped carbon intensity
    """
    # read data into array from file
    data = np.genfromtxt(data_path, delimiter=",")

    lmp = data[:time_steps, 0]
    co2i = data[:time_steps, 1]

    return lmp + co2i*cost_of_carbon, lmp, co2i


def build_flowsheet(n_steps, variable_efficiency=VariableEfficiency.none):
    # create mp model
    mp_swro = create_multiperiod_swro_model(n_time_points=n_steps,
                                            variable_efficiency=variable_efficiency)
    return mp_swro


def set_objective(
    mp_swro,
    price_signal,
    lmp, 
    co2i,
    product_requirement=8,
    oversize_factor = 1.4
):
    """
    mp_swro: pyomo model
    lmp: price signal
    co2i: carbon intensity signal (default = 0)
    carbontax: co2i cost weighting parameter
    permeate_tank_quality_constraint: If true, activates a multi-time-step constraint for permeate quality
    yield_equality_constraint: If true, fixes the daily permeate production in an equality constraint. This turns the
    levelized cost optimization into an electricity cost optimization problem.
    """
    # Retrieve pyomo model and active process blocks (i.e. time blocks)
    m = mp_swro.pyomo_model
    t_blocks = mp_swro.get_active_process_blocks()
    
    m.mp = Block()  # create a block to store parameters at the multiperiod level

    m.mp.product_requirement = Param(
        initialize=product_requirement,
        mutable = True, 
        units = pyunits.m**3,
        doc="Water production requirement over all time steps[m3]",
    )

    m.mp.time_steps = len(t_blocks)

    m.mp.baseline_production = Param(initialize = value(
                                                pyunits.convert(
                                                m.mp.product_requirement,
                                                to_units = pyunits.m**3)) / m.mp.time_steps,
                                        mutable = True,
                                        units = pyunits.m**3,
                                        doc = "Baseline hourly production [m3]")
   
    m.mp.baseline_load = Expression(
        expr = m.blocks[0].process.fs.costing.specific_energy_consumption * m.mp.baseline_production,
        doc = "Baseline hourly load [kWh]")
    

    m.mp.oversize_factor = Param(initialize = oversize_factor,
                                mutable = True,
                                doc = "Ratio between maximum and design water production")

    # index the flowsheet for each timestep
    for count, blk in enumerate(t_blocks):
        # set price and carbon signals as parameters
        blk.fs.dynamic.price_signal = Param(
            default = price_signal[count],
            mutable=True,
            units = blk.fs.costing.base_currency / pyunits.MWh
        )

        blk.fs.dynamic.lmp_signal = Param(
            default=lmp[count],
            mutable=True,
            units=blk.fs.costing.base_currency / pyunits.MWh,
        )

        blk.fs.dynamic.carbon_intensity = Param(
            default=co2i[count], 
            mutable=True, 
            units=pyunits.kg / pyunits.MWh
        )

        blk.fs.dynamic.hourly_water_production = Var(initialize = value(m.mp.baseline_production),
                                                     bounds = (0, 
                                                               value(m.mp.oversize_factor 
                                                                     * m.mp.baseline_production)),
                                                     units = pyunits.m**3 / pyunits.hr,
                                                     doc = "Hourly water production [m3/hr]")
        
        blk.fs.production_constraint = Constraint(
            expr = blk.fs.dynamic.hourly_water_production == 
            pyunits.convert(blk.fs.costing.annual_water_production,
                            to_units = pyunits.m**3 / pyunits.hr),
            doc = "Set the hourly water production on the dynamic block to each time block",
        )

        blk.fs.dynamic.hourly_energy_consumption = Expression(
            expr= pyunits.convert(blk.fs.costing.specific_energy_consumption 
                                  * blk.fs.dynamic.hourly_water_production,
                                  to_units = pyunits.kWh / pyunits.hr),
            doc="Energy consumption [kWh] per time step",
        )

        blk.fs.dynamic.hourly_carbon_emissions = Expression(
            expr=pyunits.convert(blk.fs.dynamic.hourly_energy_consumption, to_units=pyunits.kWh / pyunits.hr)
            * pyunits.convert(blk.fs.dynamic.carbon_intensity, to_units=pyunits.kg / pyunits.kWh),
            doc="Equivalent carbon emissions per timestep ",
        )

        # remove objective at the individual time-step level
        blk.fs.objective.deactivate()

    m.mp.product_requirement_constraint = Constraint(
        expr=sum([blk.fs.dynamic.hourly_water_production for blk in t_blocks]) == m.mp.product_requirement,
        doc="Daily water production requirement",
    )

    # compile time block level expressions into a model-level objective
    m.mp.obj = Objective(
        expr = sum([pyunits.convert(blk.fs.dynamic.hourly_energy_consumption, to_units = pyunits.kWh / pyunits.hr)
                     * pyunits.convert(blk.fs.dynamic.price_signal, to_units = blk.fs.costing.base_currency / pyunits.kWh) 
                     for blk in t_blocks]),
        doc = "Daily electricity and carbon cost to produce water", 
    )

    # fix the initial pressure to default operating pressure at 1 kg/s and 50% recovery
    t_blocks[0].fs.previous_pressure.fix(50e5)

    return m, t_blocks


def solve(m):
    # solve
    opt = get_solver()
    results = opt.solve(m, tee=True)
    return m, results


def save_results(m, savepath=None):
    """
    Description: saves results as a dataframe
    m: pyomo model
    t_blocks: model block for each time step
    data: price signal (comprised of LMP and CO2i
    savepath: path to save results (in .csv format)
    """
    t_blocks = m.get_active_process_blocks()

    recovery = np.array(
        [
            value(blk.fs.RO.recovery_mass_phase_comp[0, "Liq", "H2O"])
            for blk in t_blocks
        ]
    )
    pressure = np.array(
        [
            value(blk.fs.P1.control_volume.properties_out[0].pressure)
            for blk in t_blocks
        ]
    )
    energy_consumption = np.array([value(blk.fs.dynamic.hourly_energy_consumption) 
                      for blk in t_blocks])
    
    permeate = np.array([value(blk.fs.dynamic.hourly_water_production)
                          for blk in t_blocks])
    
    price_signal = np.array([value(blk.fs.dynamic.price_signal)
                              for blk in t_blocks])
    
    lmp_signal = np.array([value(blk.fs.dynamic.lmp_signal)
                              for blk in t_blocks])
    
    carbon_intensity = np.array([value(blk.fs.dynamic.carbon_intensity)
                                for blk in t_blocks])
    
    if savepath is None:
        savepath = os.path.join(os.getcwd(), "results.csv")
    
    output_data = np.vstack([recovery, pressure, energy_consumption, permeate, price_signal, lmp_signal, carbon_intensity])

    df = pd.DataFrame(output_data.T, columns = ["recovery", "pressure", "energy_consumption", "permeate", "price_signal", "lmp_signal", "carbon_intensity"])
    df.to_csv(savepath, index=True, index_label="time_step")
    return df


def visualize_results(
                path, 
                erd_type=swro.ERDtype.pump_as_turbine, 
                co2i=False, 
                title_label=""
            ):
    df = pd.read_csv(path)
    time_step = df["time"].values

    fig, ax = plt.subplots(3, 2)
    fig.suptitle(
        f"SWRO MultiPeriod optimization results: {len(time_step)} hours\n"
        + title_label,
        y=0.95,
    )

    ax[0, 0].plot(time_step, df["lmp"].values, color="black", label="LMP")
    ax[0, 0].set_xlabel("Time [hr]")
    ax[0, 0].set_ylabel("Electricity price [$/kWh]", color="black")

    ax2 = ax[0, 0].twinx()
    ax2.plot(time_step, df["co2i"].values, color="forestgreen", label="CO2i")
    ax2.set_ylabel("Carbon intensity [kgCO2/kWh]", color="forestgreen")

    ax[0, 1].plot(time_step, df["power"].values / max(df["power"].values))
    ax[0, 1].set_xlabel("Time [hr]")
    ax[0, 1].set_ylabel("Power [fraction of max]")

    ax[1, 0].plot(time_step, df["recovery"].values * 100)
    ax[1, 0].set_xlabel("Time [hr]")
    ax[1, 0].set_ylabel("Water Recovery [%]")

    ax[1, 1].plot(
        time_step, df["pressure"].values * 1e-5 / 80
    )  # normalized by burst pressure
    ax[1, 1].set_xlabel("Time [hr]")
    ax[1, 1].set_ylabel("Fraction of 80 bar")
    ax[1, 1].set_title("RO Inlet Pressure")

    ax[2, 0].plot(
        time_step,
        df["permeate"].values / max(df["permeate"].values),
        label="Main RO pump",
    )
    ax[2, 0].set_xlabel("Time [hr]")
    ax[2, 0].set_ylabel("Fraction of max flowrate")
    ax[2, 0].set_title("Permeate flowrate [m3/hr]")
    ax[2, 0].legend()

    ax[2, 1].plot(time_step, df["efficiency1"].values / 0.8, label="Main RO pump")
    ax[2, 1].set_xlabel("Time [hr]")
    ax[2, 1].set_ylabel("Fraction of BEP efficiency")
    ax[2, 1].set_title("Pump efficiency [%]")
    ax[2, 1].legend()

    if erd_type is swro.ERDtype.pressure_exchanger:
        ax[1, 0].plot(
            time_step,
            df["flowrate2"].values / max(df["flowrate2"].values),
            label="Booster pump",
        )
        ax[1, 1].plot(time_step, df["efficiency2"].values / 0.8, label="Booster pump")


if __name__ == "__main__":
    m, df = main(
        n_steps=24,
        cost_of_carbon=0.0,
        daily_production=42,
    )

    # path = os.path.join(os.getcwd(),  "temp.csv")
    # save_results(m, t_blocks, data, path)

    # visualize_results(
    #     path,
    #     erd_type=t_blocks[0].ro_mp.fs.erd_type,
    #     title_label="",
    # )


    # case = 2
    # if case == 1:
    #     price_multiplier = [0.25, 1, 3, 10]
    #     file_save = [
    #         "simulation_data_0_25_ubcon.csv",
    #         # "simulation_data_0_25_unc.csv",
    #         "simulation_data_1_00_ubcon.csv",
    #         # "simulation_data_1_00_unc.csv",
    #         "simulation_data_3_00_ubcon.csv",
    #         # "simulation_data_3_00_unc.csv",
    #         "simulation_data_10_00_ubcon.csv",
    #         # "simulation_data_10_00_unc.csv",
    #     ]
    #     for i, multiplier in enumerate(price_multiplier):

    #         m, t_blocks, data = main(
    #             ndays=1,
    #             # filename="pricesignals_GOLETA_6_N200_20220601.csv",
    #             filename="dagget_CA_LMP_hourly_2015.csv",
    #             price_multiplier=multiplier,
    #         )
    #         path = os.path.join(os.getcwd(), file_save[i])
    #         save_results(m, t_blocks, data, path)
    # elif case == 2:
    #     price_multiplier = 10
    #     m, t_blocks, data = main(
    #         ndays=1,
    #         # filename="pricesignals_GOLETA_6_N200_20220601.csv",
    #         filename="dagget_CA_LMP_hourly_2015.csv",
    #         price_multiplier=price_multiplier,
    #     )
    #     file_save = "temp.csv"
    #     path = os.path.join(os.getcwd(), file_save)
    #     save_results(m, t_blocks, data, path)

    #     title_label = f"LMP scaling = {price_multiplier}"
    #     visualize_results(
    #         path,
    #         erd_type=t_blocks[0].ro_mp.fs.erd_type,
    #         co2i=False,
    #         title_label=title_label,
    #     )
