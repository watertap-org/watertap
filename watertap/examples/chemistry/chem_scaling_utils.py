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
    This is a set of scaling utility helper functions used by the chemistry tests.
"""
from idaes.core.util import scaling as iscale
from pyomo.environ import value

__author__ = "Austin Ladshaw"

## Helper function for setting eps values associated with solubility_product functions
# NOTE: Function does nothing if no solubility reactions are present
def _set_eps_vals(rxn_params, rxn_config):
    # For solubility reactions, have to set the eps value
    try:
        for rid in rxn_params.equilibrium_reaction_idx:
            # Grab the 'k_eq_ref' value from the reaction config
            scale = rxn_config["equilibrium_reactions"][rid]["parameter_data"]["k_eq_ref"][0]

            # NOTE: ONLY the solubility_product function has an eps value that we need to set
            try:
                # Want to set eps in some fashion similar to this
                if scale < 1e-16:
                    rxn_params.component("reaction_"+rid).eps.value = scale*1e-2
                else:
                    rxn_params.component("reaction_"+rid).eps.value = 1e-16*1e-2
            except:
                pass
    except:
        pass

## Helper function for setting scaling factors for equilibrium reactions
def _set_equ_rxn_scaling(unit, rxn_config):
    #Add scaling factors for reactions (changes depending on if it is a log form or not)
    min_scale = 1e-3
    for i in unit.control_volume.equilibrium_reaction_extent_index:
        # i[0] = time, i[1] = reaction

        # Grab the 'k_eq_ref' value from the reaction config
        scale = max(min_scale, rxn_config["equilibrium_reactions"][i[1]]["parameter_data"]["k_eq_ref"][0])
        iscale.set_scaling_factor(unit.control_volume.equilibrium_reaction_extent[0.0,i[1]], 10/scale)
        iscale.constraint_scaling_transform(
            unit.control_volume.reactions[0.0].equilibrium_constraint[i[1]], 0.1)

## Helper function for setting scaling factors for inherent reactions
def _set_inherent_rxn_scaling(unit, thermo_config):
    #Add scaling factors for reactions (changes depending on if it is a log form or not)
    min_scale = 1e-3
    for i in unit.control_volume.inherent_reaction_extent_index:
        # i[0] = time, i[1] = reaction

        # Grab the 'k_eq_ref' value from the reaction config
        scale = max(min_scale, thermo_config["inherent_reactions"][i[1]]["parameter_data"]["k_eq_ref"][0])
        iscale.set_scaling_factor(unit.control_volume.inherent_reaction_extent[0.0,i[1]], 10/scale)
        iscale.constraint_scaling_transform(unit.control_volume.properties_out[0.0].
                inherent_equilibrium_constraint[i[1]], 0.1)

## Helper function for setting scaling factors for the mass balance for FpcTP state vars
def _set_mat_bal_scaling_FpcTP(unit):
    # For species
    min = 1e-3
    for i in unit.control_volume.properties_out[0.0].mole_frac_phase_comp:
        # i[0] = phase, i[1] = species
        if unit.inlet.flow_mol_phase_comp[0, i[0], i[1]].value > min:
            scale = unit.inlet.flow_mol_phase_comp[0, i[0], i[1]].value
        else:
            scale = min

        iscale.set_scaling_factor(unit.control_volume.properties_out[0.0].mole_frac_comp[i[1]], 10/scale)
        iscale.set_scaling_factor(unit.control_volume.properties_out[0.0].mole_frac_phase_comp[i], 10/scale)
        iscale.set_scaling_factor(unit.control_volume.properties_out[0.0].flow_mol_phase_comp[i], 10/scale)
        iscale.constraint_scaling_transform(unit.control_volume.material_balances[0.0,i[1]], 10/scale)

## Helper function for setting scaling factors for the mass balance for FTPx state vars
def _set_mat_bal_scaling_FTPx(unit):
    # For species
    min = 1e-3
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

## Helper function for setting energy balance scaling factors
def _set_ene_bal_scaling(unit):
    max = 1
    min = 1
    for phase in unit.control_volume.properties_in[0.0].enth_mol_phase:
        val = abs(value(unit.control_volume.properties_in[0.0].enth_mol_phase[phase].expr))
        if val >= max:
            max = val
        if val <= min:
            val = min
        iscale.set_scaling_factor(unit.control_volume.properties_in[0.0]._enthalpy_flow_term[phase], 10/val)
        iscale.set_scaling_factor(unit.control_volume.properties_out[0.0]._enthalpy_flow_term[phase], 10/val)

    iscale.constraint_scaling_transform(unit.control_volume.enthalpy_balances[0.0], 10/max)
