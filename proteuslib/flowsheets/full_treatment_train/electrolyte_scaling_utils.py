###############################################################################
# ProteusLib Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/nawi-hub/proteuslib/"
#
###############################################################################

"""
    Scaling utilities for flowsheets involving electrolyte chemistry

    NOTE: Some of these scaling methods will likely need to be updated or
    changed after the new log form for electrolytes is added to IDAES.

    ALSO NOTE: These scaling methods assume you will be using the FTPx state definition.
"""

from idaes.core.util import scaling as iscale
from idaes.core.util.initialization import fix_state_vars, revert_state_vars

# Import specific pyomo objects
from pyomo.environ import (value, Suffix)

__author__ = "Austin Ladshaw"

# Given a unit model with corresponding set of reaction params and reaction config dictionary
#   find an approximate set of state_args and/or stoich_extents to be used in initialization
#   and scaling for electrolyte systems.
#
#       Assumes using FTPx as the state_vars
def approximate_chemical_state_args(unit, rxn_params, reaction_config, contains_stoich_reactions=False):
    state_args = {}
    stoich_extents = {}

    # Set bulk values to their inlets
    state_args['pressure'] = unit.inlet.pressure[0].value
    state_args['temperature'] = unit.inlet.temperature[0].value
    state_args['flow_mol'] = unit.inlet.flow_mol[0].value

    # Set species based on inlets (and outlets for stoich reaction)
    state_args['mole_frac_comp'] = {}
    min = 1e-6
    for i in unit.control_volume.properties_in[0.0].mole_frac_comp:
        # Set state args to inlets on first pass
        if unit.inlet.mole_frac_comp[0, i].value > min:
            state_args['mole_frac_comp'][i] = unit.inlet.mole_frac_comp[0, i].value
        else:
            state_args['mole_frac_comp'][i] = min

    # Iterate through outlet mole fractions and note the fixed variables
    fixed = {}
    for i, species in unit.outlet.mole_frac_comp:
        if unit.outlet.mole_frac_comp[i, species].is_fixed():
            fixed[species] = True
            state_args['mole_frac_comp'][species] = unit.outlet.mole_frac_comp[0, species].value

    # Checking stoich reactions
    was_OH_changed = False
    was_H_changed = False
    if contains_stoich_reactions == True:
        for rid in rxn_params.rate_reaction_idx:
            #First loop establishes reaction extent
            extent = 0
            for phase, species in reaction_config["rate_reactions"][rid]["stoichiometry"]:
                # If a species here has its outlet fixed, then the difference between
                #   that species outlet and inlet values should serve as the basis for
                #   setting the values of the other species used in that reaction
                if species in fixed:
                    extent = unit.inlet.mole_frac_comp[0, species].value \
                            - unit.outlet.mole_frac_comp[0, species].value
                stoich_extents[rid] = extent

            # Loop again to set values based on extent
            for phase, species in reaction_config["rate_reactions"][rid]["stoichiometry"]:
                state_args['mole_frac_comp'][species] = unit.inlet.mole_frac_comp[0, species].value \
                        + extent*reaction_config["rate_reactions"][rid]["stoichiometry"][phase, species]
                if species == "H_+" and extent != 0.0:
                    was_H_changed = True
                if species == "OH_-" and extent != 0.0:
                    was_OH_changed = True

    # Lastly, we need for correct OH and/or H if they are changed by a stoich reaction
    if was_H_changed == False and was_OH_changed == False:
        state_args['mole_frac_comp']['H_+'] = 10**-7/55.6
        state_args['mole_frac_comp']['OH_-'] = 10**-7/55.6
    elif was_H_changed == False and was_OH_changed == True:
        state_args['mole_frac_comp']['H_+'] = 10**-14/(state_args['mole_frac_comp']['OH_-']*55.6)/55.6
    elif was_H_changed == True and was_OH_changed == False:
        state_args['mole_frac_comp']['OH_-'] = 10**-14/(state_args['mole_frac_comp']['H_+']*55.6)/55.6
    else:
        # Uncertain what to do at this point
        pass

    return state_args, stoich_extents


# Perform scaling transformations for inherent reactions (if they exist)
def calculate_chemical_scaling_factors_for_inherent_log_reactions(unit, thermo_params):
    try:
        # Iterate through the reactions to set appropriate eps values
        factor = 1e-4
        for rid in thermo_params.inherent_reaction_idx:
            scale = value(unit.control_volume.properties_out[0.0].k_eq[rid].expr)
            # Want to set eps in some fashion similar to this
            if scale < 1e-16:
                thermo_params.component("reaction_"+rid).eps.value = scale*factor
            else:
                thermo_params.component("reaction_"+rid).eps.value = 1e-16*factor

        for i in unit.control_volume.inherent_reaction_extent_index:
            scale = value(unit.control_volume.properties_out[0.0].k_eq[i[1]].expr)
            iscale.set_scaling_factor(unit.control_volume.inherent_reaction_extent[0.0,i[1]], 10/scale)
            iscale.constraint_scaling_transform(unit.control_volume.properties_out[0.0].
                    inherent_equilibrium_constraint[i[1]], 0.1)
    except:
        pass

# Perform scaling transformations for equilibrium reactions (if they exist)
def calculate_chemical_scaling_factors_for_equilibrium_log_reactions(unit, rxn_params):
    try:
        # Equilibrium reactions have eps in the 'simple_naocl_rxn_params'
        factor = 1e-4
        for rid in rxn_params.equilibrium_reaction_idx:
            if rid != "dummy":
                scale = value(unit.control_volume.reactions[0.0].k_eq[rid].expr)
                # Want to set eps in some fashion similar to this
                if scale < 1e-16:
                    rxn_params.component("reaction_"+rid).eps.value = scale*factor
                else:
                    rxn_params.component("reaction_"+rid).eps.value = 1e-16*factor

        for i in unit.control_volume.equilibrium_reaction_extent_index:
            if i[1] != "dummy":
                scale = value(unit.control_volume.reactions[0.0].k_eq[i[1]].expr)
                iscale.set_scaling_factor(unit.control_volume.equilibrium_reaction_extent[0.0,i[1]], 10/scale)
                iscale.constraint_scaling_transform(unit.control_volume.reactions[0.0].
                        equilibrium_constraint[i[1]], 0.1)
    except:
        pass

# # Perform scaling transformations for mass balances
def calculate_chemical_scaling_factors_for_material_balances(unit):
    # Next, try adding scaling for species
    min = 1e-6
    for i in unit.control_volume.properties_out[0.0].mole_frac_phase_comp:
        # i[0] = phase, i[1] = species
        if unit.inlet.mole_frac_comp[0, i[1]].value > min:
            scale = unit.inlet.mole_frac_comp[0, i[1]].value
        else:
            scale = min
        iscale.set_scaling_factor(unit.control_volume.properties_out[0.0].mole_frac_comp[i[1]], 10/scale)
        iscale.set_scaling_factor(unit.control_volume.properties_out[0.0].mole_frac_phase_comp[i], 10/scale)
        iscale.set_scaling_factor(unit.control_volume.properties_out[0.0].flow_mol_phase_comp[i], 10/scale)
        iscale.constraint_scaling_transform(
            unit.control_volume.properties_out[0.0].component_flow_balances[i[1]], 10/scale)
        iscale.constraint_scaling_transform(unit.control_volume.material_balances[0.0,i[1]], 10/scale)

# # Perform scaling transformations for energy balances
def calculate_chemical_scaling_factors_for_energy_balances(unit):
    max = 1
    min = 1
    try:
        for phase in unit.control_volume.properties_in[0.0].enth_mol_phase:
            val = abs(value(unit.control_volume.properties_in[0.0].enth_mol_phase[phase].expr))
            if val >= max:
                max = val
            if val <= min:
                val = min
            iscale.set_scaling_factor(unit.control_volume.properties_in[0.0]._enthalpy_flow_term[phase], 10/val)
            iscale.set_scaling_factor(unit.control_volume.properties_out[0.0]._enthalpy_flow_term[phase], 10/val)

        iscale.constraint_scaling_transform(unit.control_volume.enthalpy_balances[0.0], 10/max)
    except:
        pass

# Serially calculate all scaling factors needed
# # TODO: Add more scaling calculations as needed for (i) stoich reactions, (ii) rate reactions, etc.
def calculate_chemical_scaling_factors(unit, thermo_params, rxn_params, state_args, output_jac=False):
    calculate_chemical_scaling_factors_for_inherent_log_reactions(unit, thermo_params)
    calculate_chemical_scaling_factors_for_equilibrium_log_reactions(unit, rxn_params)
    calculate_chemical_scaling_factors_for_energy_balances(unit)
    calculate_chemical_scaling_factors_for_material_balances(unit)

    # If calling multiple times, then this causes errors
    #   Catch those errors here and move on (still need to
    #       have the setting of suffixes and the auto_scale_jac function)
    try:
        iscale.calculate_scaling_factors(unit)
    except:
        pass

    flags = fix_state_vars(unit.control_volume.properties_out, state_args)
    revert_state_vars(unit.control_volume.properties_out, flags)

    iscale.constraint_autoscale_large_jac(unit)

    if output_jac == True:
        jac, nlp = iscale.get_jacobian(unit, scaled=True)
        print("Extreme Jacobian entries:")
        for i in iscale.extreme_jacobian_entries(jac=jac, nlp=nlp, large=100):
            print(f"    {i[0]:.2e}, [{i[1]}, {i[2]}]")
        print("Unscaled constraints:")
        for c in iscale.unscaled_constraints_generator(unit):
            print(f"    {c}")
        print("Scaled constraints by factor:")
        for c, s in iscale.constraints_with_scale_factor_generator(unit):
            print(f"    {c}, {s}")
        print("Badly scaled variables:")
        for v, sv in iscale.badly_scaled_var_generator(unit, large=1e2, small=1e-2, zero=1e-12):
            print(f"    {v} -- {sv} -- {iscale.get_scaling_factor(v)}")
        print(f"Jacobian Condition Number: {iscale.jacobian_cond(jac=jac):.2e}")
