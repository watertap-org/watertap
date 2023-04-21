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
model for a reverse osmosis system. The purpose of this script is to 
demonstrate the use of IDAES grid-integration tools with WaterTAP process models.
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
from idaes.core.solvers import get_solver

import watertap.examples.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery as swro
from watertap.unit_models.pressure_changer import VariableEfficiency


def create_base_model(m = None, solver = None):

    print("\nCreating RO flowsheet and MP concrete model...")

    # call main to create and initialize model with flow-type variable efficiency
    if m is None:
        m = ConcreteModel()

    if solver is None:
        solver = get_solver()
    
    m = swro.build(erd_type=swro.ERDtype.pump_as_turbine, 
                   variable_efficiency=VariableEfficiency.none,)

    return m


def create_swro_mp_block():
    print(">>> Creating model for each time period")

    m = create_base_model()

    ramp_time = 60  # seconds (1 min)
    ramping_rate = 0.7e5  # Pa/s
    pressure_ramp = ramping_rate * ramp_time

    # Add coupling variables
    m.previous_pressure = Var(
        domain=NonNegativeReals,
        units=pyunits.Pa,
        bounds=(10e5, 80e5),
        doc="Applied pressure at the previous time step",
    )

    @m.Constraint(doc="Pressure ramping down constraint")
    def constraint_ramp_down(b):
        return (
            b.previous_pressure - pressure_ramp
            <= b.fs.P1.control_volume.properties_out[0].pressure
        )

    @m.Constraint(doc="Pressure ramping up constraint")
    def constraint_ramp_up(b):
        return (
            b.previous_pressure + pressure_ramp
            >= b.fs.P1.control_volume.properties_out[0].pressure
        )

    return m

def unfix_dof(b):
    """
    Unfixes the degrees of freedom in the model
    """

    # fix the RO membrane area and utilization factor
    b.fs.RO.area.fix()
    b.fs.costing.utilization_factor.fix(1)

    # unfix feed flow rate and fix concentration instead
    b.fs.feed.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].unfix()
    b.fs.feed.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].unfix()
    b.fs.feed.properties[0.0].mass_frac_phase_comp["Liq", "NaCl"].fix(
        0.035
    )
    b.fs.product.properties[0].mass_frac_phase_comp["Liq", "NaCl"].setub(
            0.0005
        )


# The tank level and power output are linked between the contiguous time periods
def get_swro_link_variable_pairs(b1, b2):
    """
    b1: current time block
    b2: next time block
    """
    return [
        (
            b1.fs.P1.control_volume.properties_out[0].pressure,
            b2.previous_pressure,
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
        n_time_points= n_time_points,
        process_model_func=create_swro_mp_block,
        unfix_dof_func = unfix_dof,
        linking_variable_func= get_swro_link_variable_pairs,
        periodic_variable_func= get_swro_periodic_variable_pairs,
    )

    multiperiod_swro.build_multi_period_model()
    return multiperiod_swro

if __name__ == "__main__":
    m = create_multiperiod_swro_model()
