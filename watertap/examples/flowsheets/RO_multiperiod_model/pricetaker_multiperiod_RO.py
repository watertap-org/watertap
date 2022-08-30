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
    lmp, co2i = _get_lmp(n_steps, data_path, carbon_cost)

    mp_swro = build_flowsheet(n_steps)

    initialize_system(mp_swro)

    m, t_blocks = set_objective(mp_swro, lmp, co2i)

    m, _ = solve(m)

    return m, t_blocks, [lmp, co2i]


def _get_lmp(time_steps, data_path, carbon_cost=False):
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
    return lmp, co2i


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
        blk_swro.fs.P1.flow_ratio[0].unfix()
        # v1 = blk_swro.fs.P1.control_volume.properties_out[0.0].flow_vol_phase["Liq"]
        blk_swro.fs.P1.bep_flow.fix()
        #
        # blk_swro.fs.P2.flow_ratio[0].unfix()
        # v2 = blk_swro.fs.P2.control_volume.properties_out[0.0].flow_vol_phase["Liq"]
        # blk_swro.fs.P2.bep_flow.fix(v2)

        # unfix RO control volume operational variables
        blk_swro.fs.P1.control_volume.properties_out[0.0].pressure.unfix()
        blk_swro.fs.RO.recovery_mass_phase_comp[0.0, "Liq", "H2O"].unfix()

        # unfix feed flow rate and fix concentration instead
        blk_swro.fs.feed.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].unfix()
        blk_swro.fs.feed.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].unfix()
        blk_swro.fs.feed.properties[0.0].mass_frac_phase_comp["Liq", "NaCl"].fix(0.035)


def set_objective(mp_swro, lmp, co2i, carbontax=0):
    # Retrieve pyomo model and active process blocks (i.e. time blocks)
    m = mp_swro.pyomo_model
    t_blocks = mp_swro.get_active_process_blocks()
    daily_water_production = 38 * pyunits.m**3 / pyunits.day
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

        blk.weighted_LCOW = Expression(
            expr=blk_swro.fs.costing.LCOW * blk.water_prod,
            doc="hourly cost [$/hr]",
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

    m.LCOW = Expression(
        expr=sum([blk.weighted_LCOW for blk in t_blocks])
        / sum([blk.water_prod for blk in t_blocks]),
        doc="Daily, flow-normalized cost of water",
    )
    # compile time block level expressions into a model-level objective
    m.obj = Objective(
        expr=sum([blk.weighted_LCOW for blk in t_blocks]),
        doc="Daily cost to produce water",
    )

    m.permeate_yield = Constraint(
        expr=abs(sum([blk.water_prod for blk in t_blocks]) / daily_water_production - 1)
        <= 1e-3,
        doc="The water production over the day is fixed within a tolerance",
    )

    m.permeate_quality_ub = Constraint(
        expr=(
            sum([blk.weighted_permeate_quality for blk in t_blocks])
            <= 0.0005 * sum([blk.water_prod for blk in t_blocks])
        ),
        doc="Flow-averaged permeate quality must be less than 500 ppm",
    )

    m.permeate_quality_lb = Constraint(
        expr=(sum([blk.weighted_permeate_quality for blk in t_blocks]) >= 0),
        doc="Flow-averaged permeate quality must be greater than 0 ppm",
    )
    # fix the initial pressure to default operating pressure at 1 kg/s and 50% recovery
    t_blocks[0].ro_mp.previous_pressure.fix(50e5)

    return m, t_blocks


def solve(m):
    # solve
    opt = get_solver()
    results = opt.solve(m, tee=True)
    return m, results


def visualize_results(m, t_blocks, data):
    time_step = np.array(range(len(t_blocks)))
    LCOW = value(m.LCOW)
    qual = 1e6 * value(
        sum([blk.weighted_permeate_quality for blk in t_blocks])
        / sum([blk.water_prod for blk in t_blocks])
    )

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
    # pump2_flow = np.array(
    #     [
    #         blk.ro_mp.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"]()
    #         for blk in t_blocks
    #     ]
    # )
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
    # pump2_efficiency = np.array(
    #     [blk.ro_mp.fs.P2.efficiency_pump[0]() for blk in t_blocks]
    # )

    fig, ax = plt.subplots(3, 2)
    fig.suptitle(
        "SWRO MultiPeriod optimization results: {} hours\nLCOW: {} $/m3\nPermeate quality: {}ppm".format(
            len(t_blocks), round(LCOW, 3), round(qual, 3)
        ),
        y=0.95,
    )
    # ax2 = ax[0, 0].twinx()
    ax[0, 0].plot(time_step, data[0][:], color="black", label="LMP")
    # ax2.plot(time_step, data[1][:], color="forestgreen",label="CO2i")
    ax[0, 0].set_xlabel("Time [hr]")
    ax[0, 0].set_ylabel("Electricity price [$/kWh]", color="black")
    # ax2.set_ylabel("Carbon intensity [kgCO2/kWh]", color="forestgreen")

    ax[1, 0].plot(time_step, pump1_flow * 3600, label="Main RO pump")
    # ax[1, 0].plot(time_step, pump2_flow * 3600, label="Booster pump")
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
    # ax[1, 1].plot(time_step, pump2_efficiency * 100, label="Booster pump")
    ax[1, 1].set_xlabel("Time [hr]")
    ax[1, 1].set_ylabel("Pump efficiency [%]")
    ax[1, 1].legend()

    ax[2, 1].plot(time_step, pressure * 1e-5)
    ax[2, 1].set_xlabel("Time [hr]")
    ax[2, 1].set_ylabel("RO Inlet Pressure [bar]")


if __name__ == "__main__":
    m, t_blocks, data = main(*sys.argv[1:])
    visualize_results(m, t_blocks, data)
