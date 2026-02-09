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
    check_optimal_termination,
    value,
)
from idaes.core.util.initialization import solve_indexed_blocks

from watertap.property_models.seawater_prop_pack import SeawaterStateBlockData
from watertap.property_models.NaCl_prop_pack import NaClStateBlockData
from watertap.property_models.NaCl_T_dep_prop_pack import (
    NaClStateBlockData as NaClTDepStateBlockData,
)
from watertap.core import MembraneChannel0DBlock, MembraneChannel1DBlock
from watertap.core.solvers import get_solver


__all__ = ["calculate_operating_pressure"]


def calculate_operating_pressure(
    state_block=None,
    over_pressure_factor=1.15,
    water_recovery_mass=0.5,
    salt_passage=0,
    solver=None,
):
    """
    Estimate operating pressure for RO unit model given the following arguments:

    Arguments:
        state_block: the state block of the RO feed that has the non-pressure state variables set to desired values (default=None)
        over_pressure_factor: the amount of operating pressure above the brine osmotic pressure represented as a fraction (default=1.15)
        water_recovery_mass: the mass-based fraction of inlet H2O that becomes permeate (default=0.5)
        salt_passage: the mass-based fraction of inlet salt that becomes permeate (default=0)
        solver: solver object to be used (default=None)
    """

    if any(
        isinstance(state_block, cls)
        for cls in [MembraneChannel0DBlock, MembraneChannel1DBlock]
    ):
        state_block = state_block.properties[0, 0]

    if not any(
        isinstance(state_block, cls)
        for cls in [SeawaterStateBlockData, NaClStateBlockData, NaClTDepStateBlockData]
    ):
        raise TypeError(
            "state_block must be created with SeawaterParameterBlock, NaClParameterBlock, or NaClTDepParameterBlock"
        )

    if not 0 <= salt_passage < 0.999:
        raise ValueError("salt_passage argument must be between 0 and 0.999")

    if not 1e-3 < water_recovery_mass < 0.999:
        raise ValueError("water_recovery_mass argument must be between 0.001 and 0.999")

    if not over_pressure_factor >= 1.0:
        raise ValueError(
            "over_pressure_factor argument must be greater than or equal to 1.0"
        )

    comp = state_block.params.solute_set.first()

    if comp not in ["NaCl", "TDS"]:
        raise ValueError(
            f"salt_passage calculation only supported for NaCl or TDS components but found {comp}"
        )

    if solver is None:
        solver = get_solver()

    tmp = ConcreteModel()  # create temporary model
    prop = state_block.config.parameters

    tmp.feed = prop.build_state_block([0])
    tmp.feed[0].pressure_osm_phase

    # specify state block
    tmp.feed[0].flow_mass_phase_comp["Liq", "H2O"].fix(
        value(state_block.flow_mass_phase_comp["Liq", "H2O"])
        * (1 - water_recovery_mass)
    )
    tmp.feed[0].flow_mass_phase_comp["Liq", comp].fix(
        value(state_block.flow_mass_phase_comp["Liq", comp]) * (1 - salt_passage)
    )
    tmp.feed[0].temperature.fix(value(state_block.temperature))
    tmp.feed[0].pressure.fix(101325)

    # solve state block
    results = solve_indexed_blocks(solver, [tmp.feed])

    if not check_optimal_termination(results):
        raise RuntimeError(
            "Failed to solve temporary state block for operating pressure"
        )

    op_pressure = value(tmp.feed[0].pressure_osm_phase["Liq"]) * over_pressure_factor

    return op_pressure
