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
import sys, os
import matplotlib.pyplot as plt

from pyomo.environ import (
    Param,
    Expression,
    Objective,
    SolverFactory,
    units as pyunits,
    Var,
)

from watertap.examples.flowsheets.RO_multiperiod_model.multiperiod_RO import (
    create_multiperiod_swro_model,
)


def main(
    ndays=1,
    filename="pricesignals_GOLETA_6_N200_20220601.csv",
):
    file_path = os.path.realpath(__file__)
    base_path = os.path.dirname(file_path)
    data_path = os.path.join(base_path, filename)
    # number of time steps assuming 1 hour as the base time
    n_steps = int(ndays * 24)

    # get data
    lmp, co2i = _get_lmp(n_steps, data_path)

    mp_swro = build_flowsheet(n_steps)

    m, t_blocks = set_objective(mp_swro, lmp, co2i)

    m, _ = solve(m)

    return m, t_blocks, [lmp, co2i]


def _get_lmp(time_steps, data_path):
    """
    Get price signals from data set
    :param time_steps: Number of time steps considered in MP analysis
    :param data_path: Price [$/kWh] on the same interval as time_steps
    :return: reshaped data
    """
    # read data into array from file
    data = np.genfromtxt(data_path, delimiter=",")
    lmp = data[:time_steps, 0]
    co2i = data[:time_steps, 1]

    # index only the desired number of timesteps
    return lmp, co2i


def build_flowsheet(n_steps):
    # create mp model
    mp_swro = create_multiperiod_swro_model(n_time_points=n_steps)

    return mp_swro


def set_objective(mp_swro, lmp, co2i, carbontax=0):
    # Retrieve pyomo model and active process blocks (i.e. time blocks)
    m = mp_swro.pyomo_model
    t_blocks = mp_swro.get_active_process_blocks()
    avg_yield = 14000 * pyunits.m**3 / pyunits.year

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
        blk.weighted_LCOW = Expression(
            expr=blk_swro.fs.costing.LCOW * blk_swro.fs.costing.annual_water_production,
            doc="annual flow weighted LCOW",
        )
        blk.water_prod = Expression(
            expr=blk_swro.fs.costing.annual_water_production,
            doc="annual water production",
        )
        blk.energy_consumption = Expression(
            expr=blk_swro.fs.costing.specific_energy_consumption
            * pyunits.convert(
                blk_swro.fs.costing.annual_water_production,
                to_units=pyunits.m**3 / pyunits.hour,
            ),
            doc="Energy consumption per timestep ",
        )
        blk.carbon_emission = Expression(
            expr=blk.energy_consumption
            * pyunits.convert(blk.carbon_intensity, to_units=pyunits.kg / pyunits.kWh),
            doc="Equivalent carbon emissions per timestep ",
        )
        blk.annual_carbon_cost = Expression(
            expr=blk.carbon_tax
            * pyunits.convert(blk.carbon_emission, to_units=pyunits.kg / pyunits.year)
            * (blk_swro.fs.costing.base_currency / pyunits.kg),
            doc="Annual cost associated with carbon emissions, to be used in the objective",
        )

        # deactivate each block-level objective function
        blk_swro.fs.objective.deactivate()

    # compile time block level expressions into a model-level objective
    m.obj = Objective(
        expr=(
            sum([blk.weighted_LCOW for blk in t_blocks])
            # + sum([blk.annual_carbon_cost for blk in t_blocks])
        )
        / sum([blk.water_prod for blk in t_blocks]),
        doc="Flow-integrated average cost and carbon tax on an annual basis",
    )

    # @m.Constraint
    # def average_flux(b,doc="fixes average flux across all timesteps"):
    #     return sum([blk.water_prod for blk in t_blocks]) == len(b.blocks) * avg_yield

    # fix the initial pressure to default operating pressure at 1 kg/s and 50% recovery
    t_blocks[0].ro_mp.previous_pressure.fix(55e5)

    return m, t_blocks


def solve(m):
    # solve
    opt = SolverFactory("ipopt")
    results = opt.solve(m, tee=True)

    return m, results


def visualize_results(m, t_blocks, data):
    time_step = np.array(range(len(t_blocks)))
    recovery = np.array(
        [
            blk.ro_mp.fs.RO.recovery_mass_phase_comp[0, "Liq", "H2O"].value
            for blk in t_blocks
        ]
    )

    pump1_flow = np.array(
        [
            blk.ro_mp.fs.P1.control_volume.properties_out[0].flow_vol_phase["Liq"]()
            for blk in t_blocks
        ]
    )
    pump2_flow = np.array(
        [
            blk.ro_mp.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"]()
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
    pump1_efficiency = np.array(
        [blk.ro_mp.fs.P1.efficiency_pump[0]() for blk in t_blocks]
    )
    pump2_efficiency = np.array(
        [blk.ro_mp.fs.P2.efficiency_pump[0]() for blk in t_blocks]
    )

    fig, ax = plt.subplots(3, 2)
    fig.suptitle(
        "SWRO MultiPeriod optimization results: {} hours\nLCOW: {} $/m3".format(
            len(t_blocks), round(m.obj(), 2)
        )
    )
    # ax2 = ax[0, 0].twinx()
    ax[0, 0].plot(time_step, data[0][:], color="black", label="LMP")
    # ax2.plot(time_step, data[1][:], color="forestgreen",label="CO2i")
    ax[0, 0].set_xlabel("Time [hr]")
    ax[0, 0].set_ylabel("Electricity price [$/kWh]", color="black")
    # ax2.set_ylabel("Carbon intensity [kgCO2/kWh]", color="forestgreen")

    ax[1, 0].plot(time_step, pump1_flow * 3600, label="Main RO pump")
    ax[1, 0].plot(time_step, pump2_flow * 3600, label="Booster pump")
    ax[1, 0].set_xlabel("Time [hr]")
    ax[1, 0].set_ylabel("Pump flowrate [m3/hr]")
    ax[1, 0].legend()

    ax[2, 0].plot(time_step, recovery * 100)
    ax[2, 0].set_xlabel("Time [hr]")
    ax[2, 0].set_ylabel("Water Recovery [%]")

    ax[0, 1].plot(time_step, power)
    ax[0, 1].set_xlabel("Time [hr]")
    ax[0, 1].set_ylabel("Net Power [kWh]")

    ax[1, 1].plot(time_step, pump1_efficiency * 100, label="Main RO pump")
    ax[1, 1].plot(time_step, pump2_efficiency * 100, label="Booster pump")
    ax[1, 1].set_xlabel("Time [hr]")
    ax[1, 1].set_ylabel("Pump efficiency [%]")
    ax[1, 1].legend()

    ax[2, 1].plot(time_step, pressure * 1e-5)
    ax[2, 1].set_xlabel("Time [hr]")
    ax[2, 1].set_ylabel("RO Inlet Pressure [bar]")


if __name__ == "__main__":
    m, t_blocks, data = main(*sys.argv[1:])
    visualize_results(m, t_blocks, data)
