#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
"""
This module contains utilities to be used with WaterTAP unit models.
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
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASStateBlockData
from watertap.core import MembraneChannel0DBlock, MembraneChannel1DBlock
from watertap.core.solvers import get_solver


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
        salt_passage: the mass-based fraction of inlet components that become permeate (default=0)
        solver: solver object to be used (default=None)
    """

    if any(
        isinstance(state_block, cls)
        for cls in [MembraneChannel0DBlock, MembraneChannel1DBlock]
    ):
        state_block = state_block.properties[0, 0]

    if not any(
        isinstance(state_block, cls)
        for cls in [
            SeawaterStateBlockData,
            NaClStateBlockData,
            NaClTDepStateBlockData,
            MCASStateBlockData,
        ]
    ):
        raise TypeError(
            "state_block argument must be a SeawaterParameterBlock, NaClParameterBlock, NaClTDepParameterBlock, or MCASParameterBlock"
        )

    if not isinstance(salt_passage, dict):
        assert isinstance(salt_passage, (int, float))
        if not 0 <= salt_passage < 0.999:
            raise ValueError("salt_passage argument must be between 0 and 0.999")

    comps = state_block.params.solute_set

    if not isinstance(salt_passage, dict):
        # Assume same salt passage for all solutes if a single value is provided
        salt_passage = {comp: salt_passage for comp in comps}
    elif isinstance(salt_passage, dict):
        for comp, sp in salt_passage.items():
            if not 0 <= sp < 0.999:
                raise ValueError(
                    f"salt_passage values must be between 0 and 0.999, but found {sp} for solute {comp}"
                )
    if set(salt_passage.keys()) != set(comps):
        # If it is a dict, keys must match the solute set
        raise ValueError(
            f"salt_passage keys must match keys in {comps} but found {list(salt_passage.keys())}"
        )

    if not 1e-3 < water_recovery_mass < 0.999:
        raise ValueError("water_recovery_mass argument must be between 0.001 and 0.999")

    if not over_pressure_factor >= 1.0:
        raise ValueError(
            "over_pressure_factor argument must be greater than or equal to 1.0"
        )

    if solver is None:
        solver = get_solver()

    tmp = ConcreteModel()  # Create temporary model
    prop = state_block.config.parameters

    tmp.feed = prop.build_state_block([0])
    tmp.feed[0].pressure_osm_phase

    # Specify state block
    tmp.feed[0].flow_mass_phase_comp["Liq", "H2O"].fix(
        value(state_block.flow_mass_phase_comp["Liq", "H2O"])
        * (1 - water_recovery_mass)
    )
    for comp in comps:
        tmp.feed[0].flow_mass_phase_comp["Liq", comp].fix(
            value(state_block.flow_mass_phase_comp["Liq", comp])
            * (1 - salt_passage[comp])
        )
    tmp.feed[0].temperature.fix(value(state_block.temperature))
    tmp.feed[0].pressure.fix(101325)

    # Solve state block
    results = solve_indexed_blocks(solver, [tmp.feed])

    if not check_optimal_termination(results):
        raise RuntimeError(
            f"Failed to calculate operating pressure for {state_block.name}"
        )

    op_pressure = value(tmp.feed[0].pressure_osm_phase["Liq"]) * over_pressure_factor

    return op_pressure
