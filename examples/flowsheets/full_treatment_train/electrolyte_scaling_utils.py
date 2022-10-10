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
    Scaling utilities for flowsheets involving electrolyte chemistry

    NOTE: These are specific for the demonstration flowsheet examples
"""

from idaes.core.util import scaling as iscale
from idaes.core.util.initialization import fix_state_vars, revert_state_vars

# Import specific pyomo objects
from pyomo.environ import value, Suffix

__author__ = "Austin Ladshaw"

# Given a unit model with corresponding set of reaction params and reaction config dictionary
#   find an approximate set of state_args and/or stoich_extents to be used in initialization
#   and scaling for electrolyte systems.
#
#       Assumes using FTPx as the state_vars
def approximate_chemical_state_args(
    unit, rxn_params, reaction_config, contains_stoich_reactions=False
):
    state_args = {}
    stoich_extents = {}

    # Set bulk values to their inlets
    state_args["pressure"] = unit.inlet.pressure[0].value
    state_args["temperature"] = unit.inlet.temperature[0].value
    state_args["flow_mol"] = unit.inlet.flow_mol[0].value

    # Set species based on inlets (and outlets for stoich reaction)
    state_args["mole_frac_comp"] = {}
    min = 1e-08
    for i in unit.control_volume.properties_in[0.0].mole_frac_comp:
        # Set state args to inlets on first pass
        if unit.inlet.mole_frac_comp[0, i].value > min:
            state_args["mole_frac_comp"][i] = unit.inlet.mole_frac_comp[0, i].value
        else:
            state_args["mole_frac_comp"][i] = min

    # Iterate through outlet mole fractions and note the fixed variables
    fixed = {}
    for i, species in unit.outlet.mole_frac_comp:
        if unit.outlet.mole_frac_comp[i, species].is_fixed():
            fixed[species] = True
            state_args["mole_frac_comp"][species] = unit.outlet.mole_frac_comp[
                0, species
            ].value

    return state_args, stoich_extents


# Perform scaling transformations for equilibrium reactions (if they exist)
def calculate_chemical_scaling_factors_for_equilibrium_log_reactions(unit, rxn_params):
    try:
        for i in unit.control_volume.equilibrium_reaction_extent_index:
            if i[1] != "dummy":
                scale = value(unit.control_volume.reactions[0.0].k_eq[i[1]].expr)
                iscale.set_scaling_factor(
                    unit.control_volume.equilibrium_reaction_extent[0.0, i[1]],
                    10 / scale,
                )
                iscale.constraint_scaling_transform(
                    unit.control_volume.reactions[0.0].equilibrium_constraint[i[1]], 0.1
                )
    except:
        pass


# # Perform scaling transformations for mass balances
def calculate_chemical_scaling_factors_for_material_balances(unit):
    # Next, try adding scaling for species
    min = 1e-3
    for i in unit.control_volume.properties_out[0.0].mole_frac_phase_comp:
        # i[0] = phase, i[1] = species
        if unit.inlet.mole_frac_comp[0, i[1]].value > min:
            scale = unit.inlet.mole_frac_comp[0, i[1]].value
        else:
            scale = min
        iscale.set_scaling_factor(
            unit.control_volume.properties_out[0.0].mole_frac_comp[i[1]], 10 / scale
        )
        iscale.set_scaling_factor(
            unit.control_volume.properties_out[0.0].mole_frac_phase_comp[i], 10 / scale
        )
        iscale.set_scaling_factor(
            unit.control_volume.properties_out[0.0].flow_mol_phase_comp[i], 10 / scale
        )
        iscale.constraint_scaling_transform(
            unit.control_volume.properties_out[0.0].component_flow_balances[i[1]],
            10 / scale,
        )
        iscale.constraint_scaling_transform(
            unit.control_volume.material_balances[0.0, i[1]], 10 / scale
        )


# # Perform scaling transformations for energy balances
def calculate_chemical_scaling_factors_for_energy_balances(unit):
    max = 1
    min = 1
    try:
        for phase in unit.control_volume.properties_in[0.0].enth_mol_phase:
            val = abs(
                value(unit.control_volume.properties_in[0.0].enth_mol_phase[phase].expr)
            )
            if val >= max:
                max = val
            if val <= min:
                val = min
            iscale.set_scaling_factor(
                unit.control_volume.properties_in[0.0]._enthalpy_flow_term[phase],
                10 / val,
            )
            iscale.set_scaling_factor(
                unit.control_volume.properties_out[0.0]._enthalpy_flow_term[phase],
                10 / val,
            )

        iscale.constraint_scaling_transform(
            unit.control_volume.enthalpy_balances[0.0], 10 / max
        )
    except:
        pass


# Serially calculate all scaling factors needed
def calculate_chemical_scaling_factors(
    unit, thermo_params, rxn_params, state_args, output_jac=False
):
    calculate_chemical_scaling_factors_for_equilibrium_log_reactions(unit, rxn_params)
    calculate_chemical_scaling_factors_for_energy_balances(unit)
    calculate_chemical_scaling_factors_for_material_balances(unit)

    flags = fix_state_vars(unit.control_volume.properties_out, state_args)
    revert_state_vars(unit.control_volume.properties_out, flags)
