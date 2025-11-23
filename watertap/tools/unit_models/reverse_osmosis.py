#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
"""
This module contains functions to be used with WaterTAP
ReverseOsmosis0D or ReverseOsmosis1D unit models.
"""

from pyomo.environ import (
    ConcreteModel,
    assert_optimal_termination,
    value,
)
from idaes.core.util.initialization import solve_indexed_blocks

from watertap.property_models.seawater_prop_pack import SeawaterStateBlockData
from watertap.property_models.NaCl_prop_pack import NaClStateBlockData


__all__ = ["calculate_operating_pressure"]


def calculate_operating_pressure(
    feed_state_block=None,
    over_pressure=0.15,
    water_recovery_mass=0.5,
    salt_passage=0.01,
    solver=None,
):
    """
    Estimate operating pressure for RO unit model given the following arguments:

    Arguments:
        feed_state_block:   the state block of the RO feed that has the non-pressure state
                            variables initialized to their values (default=None)
        over_pressure:  the amount of operating pressure above the brine osmotic pressure
                        represented as a fraction (default=0.15)
        water_recovery_mass: the mass-based fraction of inlet H2O that becomes permeate
                        (default=0.5)
        salt_passage:   the mass-based fraction of inlet salt that becomes permeate
                        (default=0.01)
        solver:     solver object to be used (default=None)
    """
    if not any(
        isinstance(feed_state_block, cls)
        for cls in [SeawaterStateBlockData, NaClStateBlockData]
    ):
        raise TypeError(
            "feed_state_block argument must be a SeawaterStateBlockData or NaClStateBlockData object"
        )

    if isinstance(feed_state_block, NaClStateBlockData):
        comp = "NaCl"
    if isinstance(feed_state_block, SeawaterStateBlockData):
        comp = "TDS"
        
    tmp = ConcreteModel()  # create temporary model
    prop = feed_state_block.config.parameters
    tmp.brine = prop.build_state_block([0])
    tmp.brine[0].pressure_osm_phase

    # specify state block
    tmp.brine[0].flow_mass_phase_comp["Liq", "H2O"].fix(
        value(feed_state_block.flow_mass_phase_comp["Liq", "H2O"])
        * (1 - water_recovery_mass)
    )
    tmp.brine[0].flow_mass_phase_comp["Liq", comp].fix(
        value(feed_state_block.flow_mass_phase_comp["Liq", comp]) * (1 - salt_passage)
    )
    tmp.brine[0].pressure.fix(
        101325
    )  # valid when osmotic pressure is independent of hydraulic pressure
    tmp.brine[0].temperature.fix(value(feed_state_block.temperature))
    
    # solve state block
    results = solve_indexed_blocks(solver, [tmp.brine])
    assert_optimal_termination(results)

    op_pressure = value(tmp.brine[0].pressure_osm_phase["Liq"]) * (1 + over_pressure)

    return op_pressure
