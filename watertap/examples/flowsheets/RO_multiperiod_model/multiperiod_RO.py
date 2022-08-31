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

"""
This script uses the IDAES multiperiod class to create a quasi-steady-state
model for a reverse osmosis system with energy recovery (PX) and variable
efficiency pumps. The purpose of this script is to demonstrate the use of
IDAES grid-integration tools with WaterTAP process models.
Code is heavily inspired by the DISPATCHES example of Ultra-supercritical
power plant by N.Susarla and S. Rawlings.
Repository link: https://github.com/gmlc-dispatches/dispatches
Path: models/fossil_case/ultra_supercritical_plant/storage/multiperiod_integrated_storage_usc.py
"""

__author__ = "Akshay Rao, Adam Atia"

from pyomo.environ import (
    NonNegativeReals,
    ConcreteModel,
    Var,
    units as pyunits,
    value,
    Param,
    Constraint,
)

from idaes.apps.grid_integration.multiperiod.multiperiod import MultiPeriodModel

import watertap.examples.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery as swro
from watertap.unit_models.pressure_changer import VariableEfficiency


def create_base_model():

    print("\nCreating RO flowsheet and MP concrete model...")

    # call main to create and initialize model with flow-type variable efficiency
    m = ConcreteModel()
    # m.ro_mp = swro.build(erd_type=swro.ERDtype.pump_as_turbine
    #                      ,variable_efficiency=VariableEfficiency.flow)
    # swro.set_operating_conditions(m.ro_mp)
    # swro.initialize_system(m.ro_mp, skipRO=False)
    m.ro_mp = swro.main(
        erd_type=swro.ERDtype.pump_as_turbine,
        variable_efficiency=VariableEfficiency.flow,
    )
    return m


def create_swro_mp_block():
    print(">>> Creating model for each time period")

    m = create_base_model()
    b1 = m.ro_mp

    ramp_time = 60  # seconds (1 min)
    ramping_rate = 0.7e5  # Pa/s

    # Add coupling variables
    b1.previous_pressure = Var(
        domain=NonNegativeReals,
        units=pyunits.Pa,
        bounds=(10e5, 80e5),
        doc="Applied pressure at the previous time step",
    )

    @b1.Constraint(doc="Pressure ramping down constraint")
    def constraint_ramp_down(b):
        return (
            b.previous_pressure - 40e5
            <= b.fs.P1.control_volume.properties_out[0].pressure
        )

    @b1.Constraint(doc="Pressure ramping up constraint")
    def constraint_ramp_up(b):
        return (
            b.previous_pressure + 40e5
            >= b.fs.P1.control_volume.properties_out[0].pressure
        )

    return m


# The tank level and power output are linked between the contiguous time periods
def get_swro_link_variable_pairs(b1, b2):
    """
    b1: current time block
    b2: next time block
    """
    return [
        (
            b1.ro_mp.fs.P1.control_volume.properties_out[0].pressure,
            b2.ro_mp.previous_pressure,
        )
    ]


def get_swro_periodic_variable_pairs(b1, b2):
    """
    b1: final time block
    b2: first time block
    """
    # return
    return []


def create_multiperiod_swro_model(n_time_points=4):
    """
    n_time_points: Number of time blocks
    """
    multiperiod_swro = MultiPeriodModel(
        n_time_points,
        lambda: create_swro_mp_block(),
        get_swro_link_variable_pairs,
        get_swro_periodic_variable_pairs,
    )

    multiperiod_swro.build_multi_period_model()
    return multiperiod_swro


if __name__ == "__main__":
    m = create_multiperiod_swro_model()
