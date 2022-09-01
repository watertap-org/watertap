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
import sys, os, matplotlib
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
import watertap.core.util.infeasible as infeas


def main(
    ndays=1,
    # filename="pricesignals_GOLETA_6_N200_20220601.csv",
    filename="dagget_CA_LMP_hourly_2015.csv",
    price_multiplier=1,
):
    file_path = os.path.realpath(__file__)
    base_path = os.path.dirname(file_path)
    data_path = os.path.join(base_path, filename)
    # number of time steps assuming 1 hour as the base time
    n_steps = int(ndays * 24)

    if filename == "pricesignals_GOLETA_6_N200_20220601.csv":
        carbon_cost = True
    else:
        carbon_cost = False

    # get data
    lmp, co2i = _get_lmp(n_steps, data_path, carbon_cost, price_multiplier)

    mp_swro = build_flowsheet(n_steps)

    initialize_system(mp_swro)

    m, t_blocks = set_objective(mp_swro, lmp, co2i)

    m, _ = solve(m)

    return m, t_blocks, [lmp, co2i]


def _get_lmp(time_steps, data_path, carbon_cost=False, price_multiplier=1):
    """
    Get price signals from data set
    :param time_steps: Number of time steps considered in MP analysis
    :param data_path: Price [$/kWh] on the same interval as time_steps
    :return: reshaped data
    """
    # read data into array from file
    data = np.genfromtxt(data_path, delimiter=",")

    if carbon_cost == True:
        lmp = data[:time_steps, 0]
        co2i = data[:time_steps, 1]
    else:
        lmp = data[:time_steps]
        co2i = np.zeros_like(lmp)

    # index only the desired number of timesteps
    return lmp * price_multiplier, co2i


def build_flowsheet(n_steps):
    # create mp model
    mp_swro = create_multiperiod_swro_model(n_time_points=n_steps)

    return mp_swro


def initialize_system(m):
    print(f"\n----- Initializing Base Model -----")
    # fully initialize a base model
    # base_model = swro.main(erd_type=swro.ERDtype.pump_as_turbine,
    #     variable_efficiency=swro.VariableEfficiency.flow)
    # # create a copy of the base model
    # reinitialize_values = ComponentMap()
    # for v in base_model.component_data_objects(Var):
    #     reinitialize_values[v] = v.value
    #

    # loop through each time step block and set values
    t_blocks = m.get_active_process_blocks()

    for count, blk in enumerate(t_blocks):
        print(f"\n----- Initializing MultiPeriod Time Step {count} -----")
        blk_swro = blk.ro_mp
        # for v, val in reinitialize_values.items():
        #     blk_swro.find_component(v).set_value(val, skip_validation=True)
        blk_swro.fs.RO.area.fix()
        # unfix the pump flow ratios and fix the bep flowrate as the nominal volumetric flowrate
        blk_swro.fs.P1.bep_flow.fix()
        blk_swro.fs.P1.flow_ratio[0].unfix()
        blk_swro.fs.costing.utilization_factor.fix(1)
        # v1 = blk_swro.fs.P1.control_volume.properties_out[0.0].flow_vol_phase["Liq"]
        #
        # blk_swro.fs.P2.flow_ratio[0].unfix()
        # v2 = blk_swro.fs.P2.control_volume.properties_out[0.0].flow_vol_phase["Liq"]
        # blk_swro.fs.P2.bep_flow.fix(v2)

        # unfix RO control volume operational variables
        blk_swro.fs.P1.control_volume.properties_out[0.0].pressure.unfix()
        blk_swro.fs.RO.recovery_mass_phase_comp[0.0, "Liq", "H2O"].unfix()
        blk_swro.fs.product.properties[0].mass_frac_phase_comp["Liq", "NaCl"].setub(
            0.0005
        )
        # unfix feed flow rate and fix concentration instead
        blk_swro.fs.feed.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].unfix()
        blk_swro.fs.feed.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].unfix()
        blk_swro.fs.feed.properties[0.0].mass_frac_phase_comp["Liq", "NaCl"].fix(0.035)


def set_objective(mp_swro, lmp, co2i, carbontax=0):
    # Retrieve pyomo model and active process blocks (i.e. time blocks)
    m = mp_swro.pyomo_model
    t_blocks = mp_swro.get_active_process_blocks()
    daily_water_production = 38.52 * pyunits.m**3 / pyunits.day
    fixed_hourly_cost = 0.6102
    # base_costing_block = t_blocks[0].ro_mp.fs.costing
    # fixed_hourly_cost = Expression(
    #     expr= pyunits.convert(base_costing_block.total_investment_cost
    #         * base_costing_block.factor_capital_annualization
    #         + base_costing_block.maintenance_labor_chemical_operating_cost
    #         + base_costing_block.aggregate_fixed_operating_cost,
    #         to_units=pyunits.USD_2018 / pyunits.hr),
    #     doc="Base cost of the plant as defined by the initialized model at t=0 [$/m3]"
    # )

    # index the flowsheet for each timestep
    for count, blk in enumerate(t_blocks):
        blk_swro = blk.ro_mp

        # set price and carbon signals as parameters
        blk.lmp_signal = Param(
            default=lmp[count],
            mutable=True,
            units=blk_swro.fs.costing.base_currency / pyunits.kWh,
        )
        blk.carbon_intensity = Param(
            default=co2i[count], mutable=True, units=pyunits.kg / pyunits.MWh
        )
        blk.carbon_tax = Param(
            default=carbontax,
            mutable=True,
            units=blk_swro.fs.costing.base_currency / pyunits.kg,
        )

        # set the electricity_price in each flowsheet
        blk_swro.fs.costing.electricity_cost.fix(blk.lmp_signal)

        # combine/place flowsheet level cost metrics on each time block
        blk.water_prod = Expression(
            expr=pyunits.convert(
                blk_swro.fs.costing.annual_water_production,
                to_units=pyunits.m**3 / pyunits.hour,
            ),
            doc="annual water production",
        )
        blk.weighted_permeate_quality = Expression(
            expr=blk.water_prod
            * blk_swro.fs.product.properties[0].mass_frac_phase_comp["Liq", "NaCl"],
            doc="Permeate flow weighted concentration of salt in the permeate",
        )
        blk.energy_consumption = Expression(
            expr=blk_swro.fs.costing.specific_energy_consumption * blk.water_prod,
            doc="Energy consumption per timestep kWh/hr",
        )
        blk.carbon_emission = Expression(
            expr=blk.energy_consumption
            * pyunits.convert(blk.carbon_intensity, to_units=pyunits.kg / pyunits.kWh),
            doc="Equivalent carbon emissions per timestep ",
        )
        blk.annual_carbon_cost = Expression(
            expr=blk.carbon_tax
            * (blk_swro.fs.costing.base_currency / pyunits.kg)
            * pyunits.convert(blk.carbon_emission, to_units=pyunits.kg / pyunits.hour),
            doc="Hourly cost associated with carbon emissions $/hr",
        )
        blk.weighted_LCOW = Expression(
            expr=fixed_hourly_cost + blk.energy_consumption * blk.lmp_signal,
            doc="hourly cost [$/hr]",
        )
        blk.weighted_LCOW_b = Expression(
            expr=blk.ro_mp.fs.costing.LCOW * blk.water_prod,
            doc="hourly cost [$/hr]",
        )

    m.LCOW = Expression(
        expr=sum([blk.weighted_LCOW_b for blk in t_blocks])
        / sum([blk.water_prod for blk in t_blocks]),
        doc="Daily, flow-normalized cost of water",
    )

    # compile time block level expressions into a model-level objective
    m.obj = Objective(
        expr=sum([blk.weighted_LCOW for blk in t_blocks])
        / sum([blk.water_prod for blk in t_blocks]),
        doc="Daily cost to produce water",
    )

    m.permeate_yield = Constraint(
        expr=(
            daily_water_production - 1e-5,
            sum(blk.water_prod for blk in t_blocks),
            daily_water_production + 1e-5,
        ),
        doc="The water production over the day is fixed within a tolerance",
    )

    # m.permeate_quality_ub = Constraint(
    #     expr=(
    #         sum([blk.weighted_permeate_quality for blk in t_blocks])
    #         <= 0.0005 * sum([blk.water_prod for blk in t_blocks])
    #     ),
    #     doc="Flow-averaged permeate quality must be less than 500 ppm",
    # )
    #
    # m.permeate_quality_lb = Constraint(
    #     expr=(sum([blk.weighted_permeate_quality for blk in t_blocks]) >= 0),
    #     doc="Flow-averaged permeate quality must be greater than 0 ppm",
    # )
    # fix the initial pressure to default operating pressure at 1 kg/s and 50% recovery
    t_blocks[0].ro_mp.previous_pressure.fix(50e5)

    return m, t_blocks


def solve(m):
    # solve
    opt = get_solver()
    results = opt.solve(m, tee=True)
    return m, results


def save_results(m, t_blocks, data, savepath=None):
    time_step = np.array(range(len(t_blocks)))
    LCOW = value(m.obj)
    # qual = 1e6 * value(
    #     sum([blk.weighted_permeate_quality for blk in t_blocks])
    #     / sum([blk.water_prod for blk in t_blocks])
    # )
    print(LCOW)

    recovery = np.array(
        [
            blk.ro_mp.fs.RO.recovery_mass_phase_comp[0, "Liq", "H2O"].value
            for blk in t_blocks
        ]
    )
    pressure = np.array(
        [
            blk.ro_mp.fs.P1.control_volume.properties_out[0].pressure.value
            for blk in t_blocks
        ]
    )
    power = np.array([blk.energy_consumption() for blk in t_blocks])

    flowrate1 = np.array(
        [
            blk.ro_mp.fs.P1.control_volume.properties_out[0].flow_vol_phase["Liq"]()
            for blk in t_blocks
        ]
    )
    efficiency1 = np.array([blk.ro_mp.fs.P1.efficiency_pump[0]() for blk in t_blocks])

    if savepath is None:
        savepath = os.path.join(os.getcwd(), "simulation_data_0_25.csv")
    if t_blocks[0].ro_mp.fs.erd_type is swro.ERDtype.pump_as_turbine:
        output_data = np.hstack(
            (
                time_step.reshape(-1, 1),
                data[0][:].reshape(-1, 1),
                data[1][:].reshape(-1, 1),
                power.reshape(-1, 1),
                recovery.reshape(-1, 1),
                pressure.reshape(-1, 1),
                flowrate1.reshape(-1, 1),
                efficiency1.reshape(-1, 1),
            )
        )
        df = pd.DataFrame(
            output_data,
            columns=[
                "time",
                "lmp",
                "co2i",
                "power",
                "recovery",
                "pressure",
                "flowrate1",
                "efficiency1",
            ],
        )

    elif t_blocks[0].ro_mp.fs.erd_type is swro.ERDtype.pressure_exchanger:
        flowrate2 = np.array(
            [
                blk.ro_mp.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"]()
                for blk in t_blocks
            ]
        )
        efficiency2 = np.array(
            [blk.ro_mp.fs.P2.efficiency_pump[0]() for blk in t_blocks]
        )
        output_data = np.hstack(
            (
                time_step.reshape(-1, 1),
                data[0][:].reshape(-1, 1),
                data[1][:].reshape(-1, 1),
                power.reshape(-1, 1),
                recovery.reshape(-1, 1),
                pressure.reshape(-1, 1),
                flowrate1.reshape(-1, 1),
                efficiency1.reshape(-1, 1),
                flowrate2.reshape(-1, 1),
                efficiency2.reshape(-1, 1),
            )
        )
        df = pd.DataFrame(
            output_data,
            columns=[
                "time",
                "lmp",
                "co2i",
                "power",
                "recovery",
                "pressure",
                "flowrate1",
                "efficiency1",
                "flowrate2",
                "efficiency2",
            ],
        )
    else:
        return

    df.to_csv(savepath)


def visualize_results(
    path, erd_type=swro.ERDtype.pump_as_turbine, co2i=False, title_label=""
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
    if co2i:
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
        df["flowrate1"].values / max(df["flowrate1"].values),
        label="Main RO pump",
    )
    ax[2, 0].set_xlabel("Time [hr]")
    ax[2, 0].set_ylabel("Fraction of max flowrate")
    ax[2, 0].set_title("Pump flowrate [m3/hr]")
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
    price_multiplier = 10
    m, t_blocks, data = main(
        ndays=1,
        # filename="pricesignals_GOLETA_6_N200_20220601.csv",
        filename="dagget_CA_LMP_hourly_2015.csv",
        price_multiplier=price_multiplier,
    )
    path = os.path.join(os.getcwd(), "fixed_simulation_data_10_00.csv")
    save_results(m, t_blocks, data, path)
    title_label = f"LMP scaling = {price_multiplier}"
    visualize_results(
        path,
        erd_type=t_blocks[0].ro_mp.fs.erd_type,
        co2i=False,
        title_label=title_label,
    )
