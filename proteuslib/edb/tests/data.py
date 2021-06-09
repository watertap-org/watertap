###############################################################################
# ProteusLib Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/nawi-hub/proteuslib/"
#
###############################################################################
"""
Shared input data for tests
"""
from pyomo.environ import units as pyunits
from idaes.core.phases import PhaseType as PT
import idaes.generic_models.properties.core.pure.Perrys as Perrys
from proteuslib.edb.equations.van_t_hoff_alt_form import van_t_hoff_aqueous
from proteuslib.edb.equations.equil_log_power_form import log_power_law
from idaes.generic_models.properties.core.reactions.dh_rxn import constant_dh_rxn
from idaes.generic_models.properties.core.generic.generic_reaction import (
    ConcentrationForm,
)


_reaction = {
    "stoichiometry": {"Liq": {"NH4 +": -1, "H +": 1, "NH3": 1}},
    "heat_of_reaction": "constant_dh_rxn",
    "equilibrium_constant": "van_t_hoff_aqueous",
    "equilibrium_form": "log_power_law",
    "concentration_form": "ConcentrationForm.molarity",
    "parameter_data": {
        "dh_rxn_ref": [{"v": 52.21, "u": "kJ/mol", "i": 0}],
        "ds_rxn_ref": [{"v": -2.4, "u": "J/mol/K", "i": 0}],
    },
    "type": "equilibrium",
    "name": "NH4_Ka",
    "components": ["NH4", "Ka"],
    "base_units": {
        "time": "s",
        "length": "m",
        "mass": "kg",
        "amount": "mol",
        "temperature": "K",
    },
}
_component = {
    "valid_phase_types": ["PT.liquidPhase"],
    "dens_mol_liq_comp": "Perrys",
    "enth_mol_liq_comp": "Perrys",
    "cp_mol_liq_comp": "Perrys",
    "entr_mol_liq_comp": "Perrys",
    "parameter_data": {
        "mw": [{"v": 74.09, "u": "g/mol"}],
        "dens_mol_liq_comp_coeff": [
            {"v": 13.5, "u": "kmol*m**-3", "i": 1},
            {"v": 1, "u": "dimensionless", "i": 2},
            {"v": 1, "u": "K", "i": 3},
            {"v": 1, "u": "dimensionless", "i": 4},
        ],
        "enth_mol_form_liq_comp_ref": [{"v": -1003, "u": "kJ/mol"}],
        "cp_mol_liq_comp_coeff": [
            {"v": 276370.0, "u": "J/kmol/K", "i": "1"},
            {"v": -2090.1, "u": "J/kmol/K**2", "i": "2"},
            {"v": 8.125, "u": "J/kmol/K**3", "i": "3"},
            {"v": -0.014116, "u": "J/kmol/K**4", "i": "4"},
            {"v": 9.3701e-06, "u": "J/kmol/K**5", "i": "5"},
        ],
        "entr_mol_form_liq_comp_ref": [{"v": -74.5, "u": "J/K/mol"}],
    },
    "name": "Ca[OH]2",
}

reaction_data = [_reaction, _reaction]
component_data = [_component, _component]

Ca_thermo_config = {
    "components": {
        "Ca 2+": {
            "type": "Component",
            "valid_phase_types": PT.liquidPhase,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Perrys,
            "enth_mol_liq_comp": Perrys,
            "cp_mol_liq_comp": Perrys,
            "entr_mol_liq_comp": Perrys,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (40.078, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": {
                    "1": (13.5, pyunits.kmol * pyunits.m ** -3),
                    "2": (1, pyunits.dimensionless),
                    "3": (1, pyunits.K),
                    "4": (1, pyunits.dimensionless),
                },
                "enth_mol_form_liq_comp_ref": (-542.83, pyunits.J / pyunits.mol),
                "cp_mol_liq_comp_coeff": {
                    "1": (2.7637e5, pyunits.J / pyunits.kmol / pyunits.K),
                    "2": (-2.0901e3, pyunits.J / pyunits.kmol / pyunits.K ** 2),
                    "3": (8.125, pyunits.J / pyunits.kmol / pyunits.K ** 3),
                    "4": (-1.4116e-2, pyunits.J / pyunits.kmol / pyunits.K ** 4),
                    "5": (9.3701e-6, pyunits.J / pyunits.kmol / pyunits.K ** 5),
                },
                "entr_mol_form_liq_comp_ref": (
                    -53,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
        }
    }
}

Ca_thermo_data = {
    "_id": {"$oid": "6099cd12507b1e55f1707555"},
    "valid_phase_types": ["PT.liquidPhase"],
    "dens_mol_liq_comp": "Perrys",
    "enth_mol_liq_comp": "Perrys",
    "cp_mol_liq_comp": "Perrys",
    "entr_mol_liq_comp": "Perrys",
    "parameter_data": {
        "mw": [{"v": 40.078, "u": "g/mol", "i": 0}],
        "dens_mol_liq_comp_coeff": [
            {"v": 13.5, "u": "kmol*m**-3", "i": 1},
            {"v": 1, "u": "dimensionless", "i": 2},
            {"v": 1, "u": "K", "i": 3},
            {"v": 1, "u": "dimensionless", "i": 4},
        ],
        "enth_mol_form_liq_comp_ref": [{"v": -542.83, "u": "J/mol"}],
        "cp_mol_liq_comp_coeff": [
            {"v": 276370, "u": "J/kmol/K", "i": 1},
            {"v": -2090.1, "u": "J/kmol/K**2", "i": 2},
            {"v": 8.125, "u": "J/kmol/K**3", "i": 3},
            {"v": -0.014116, "u": "J/kmol/K**4", "i": 4},
            {"v": 9.3701e-06, "u": "J/kmol/K**5", "i": 5},
        ],
        "entr_mol_form_liq_comp_ref": [{"v": -53, "u": "J/K/mol"}],
    },
    "name": "Ca 2+",
    "elements": ["Ca"],
}

bicarbonate_reaction_config = {
    "equilibrium_reactions": {
        "H2CO3_Ka2": {
            "stoichiometry": {
                ("Liq", "HCO3 -"): -1,
                ("Liq", "H +"): 1,
                ("Liq", "CO3 2-"): 1,
            },
            "heat_of_reaction": constant_dh_rxn,
            "equilibrium_constant": van_t_hoff_aqueous,
            # "equilibrium_constant": gibbs_energy,
            "equilibrium_form": log_power_law,
            "concentration_form": ConcentrationForm.molarity,
            "parameter_data": {
                "dh_rxn_ref": (14.9, pyunits.kJ / pyunits.mol),
                "ds_rxn_ref": (-148.1, pyunits.J / pyunits.mol / pyunits.K),
                # "T_eq_ref": (300, pyunits.K),
                # By default, reaction orders follow stoichiometry
                #    manually set reaction order here to override
                "reaction_order": {
                    ("Liq", "HCO3 -"): -1,
                    ("Liq", "H +"): 1,
                    ("Liq", "CO3 2-"): 1,
                },
            },
        }
    },
}

bicarbonate_reaction_data = {
    "_id": {"$oid": "6099cd12507b1e55f1707546"},
    "stoichiometry": {"Liq": {"HCO3-": -1, "H+": 1, "CO3-2": 1}},
    "heat_of_reaction": "constant_dh_rxn",
    "equilibrium_constant": "van_t_hoff_aqueous",
    "equilibrium_form": "log_power_law",
    "concentration_form": "ConcentrationForm.molarity",
    "parameter_data": {
        "dh_rxn_ref": [{"v": 14.9, "u": "kJ/mol"}],
        "ds_rxn_ref": [{"v": -148.1, "u": "J/mol/K"}],
        "reaction_order": {"Liq": {"HCO3 -": -1, "H +": 1, "CO3 2-": 1}},
    },
    "name": "H2CO3_Ka2",
    "components": ["H2CO3", "Ka2"],
    "reactant_elements": ["H", "O", "C"],
}
