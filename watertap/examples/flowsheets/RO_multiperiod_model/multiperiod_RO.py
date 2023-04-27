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
    Var,
    units as pyunits,
    Param,
    Constraint,
    Block,
)

from idaes.apps.grid_integration.multiperiod.multiperiod import MultiPeriodModel
from idaes.core.solvers import get_solver

import watertap.examples.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery as swro
from watertap.unit_models.pressure_changer import VariableEfficiency


def create_base_model(m = None):
    # TODO - replace call to main with a build function and defer initialization to a separate routine
    if m is None:
        m = swro.main(erd_type=swro.ERDtype.pump_as_turbine, 
                    variable_efficiency=VariableEfficiency.none,)

    return m


def create_swro_mp_block(m=None):

    m = create_base_model(m)

    # Create a dynamic block to store dynamic operation parameters 
    m.fs.dynamic = Block()
    m.fs.dynamic.ramp_time = Param(initialize=300,
                                   units= pyunits.s, 
                                   mutable=True,
                                   doc="Time associated with ramping up or down in system pressure")
    m.fs.dynamic.ramping_rate = Param(initialize=0.7e5,
                                        units= pyunits.Pa/pyunits.s,
                                        mutable=True,
                                        doc="Slowest rate at which pressure can be ramped up or down")

    # Add coupling variables
    m.fs.previous_pressure = Var(
        domain=NonNegativeReals,
        units=pyunits.Pa,
        bounds=(10e5, 80e5),
        doc="Applied pressure at the previous time step",
    )

    @m.Constraint(doc="Pressure ramping down constraint")
    def constraint_ramp_down(b):
        return (
            b.fs.previous_pressure - b.fs.dynamic.ramping_rate * b.fs.dynamic.ramp_time
            <= b.fs.P1.control_volume.properties_out[0].pressure
        )

    @m.Constraint(doc="Pressure ramping up constraint")
    def constraint_ramp_up(b):
        return (
            b.fs.previous_pressure + b.fs.dynamic.ramping_rate * b.fs.dynamic.ramp_time
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

    #unfix the recovery ratio
    b.fs.RO.recovery_mass_phase_comp[0, "Liq", "H2O"].unfix()

    # unfix feed flow rate and fix concentration instead
    b.fs.feed.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].unfix()
    b.fs.feed.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].unfix()
    b.fs.feed.properties[0.0].mass_frac_phase_comp["Liq", "NaCl"].fix(
        0.035
    )
    b.fs.product.properties[0].mass_frac_phase_comp["Liq", "NaCl"].setub(
            0.0005
        )


def get_swro_link_variable_pairs(b1, b2):
    """
    b1: current time block
    b2: next time block
    """
    return [
        (
            b1.fs.P1.control_volume.properties_out[0].pressure,
            b2.fs.previous_pressure,
        )
    ]


def get_swro_periodic_variable_pairs(b1, b2):
    """
    b1: final time block
    b2: first time block
    """
    return []


def create_multiperiod_swro_model(variable_efficiency,
                                  n_time_points=4,
                                  ):
    """
    n_time_points: Number of time blocks
    """
    # TODO - add option to pass flowsheet build and unfix dof options (variable efficiency)
    multiperiod_swro = MultiPeriodModel(
        n_time_points= n_time_points,
        process_model_func=create_swro_mp_block,
        unfix_dof_func = unfix_dof,
        linking_variable_func= get_swro_link_variable_pairs,
    )

    multiperiod_swro.build_multi_period_model()
    return multiperiod_swro

if __name__ == "__main__":
    m = create_multiperiod_swro_model()
