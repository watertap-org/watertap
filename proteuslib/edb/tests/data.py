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
from idaes.generic_models.properties.core.reactions.equilibrium_forms import log_power_law_equil
from idaes.generic_models.properties.core.reactions.equilibrium_constant import (
    van_t_hoff,
)
from idaes.generic_models.properties.core.reactions.dh_rxn import constant_dh_rxn
from idaes.generic_models.properties.core.generic.generic_reaction import (
    ConcentrationForm,
)
from idaes.generic_models.properties.core.reactions.equilibrium_constant import (
    gibbs_energy, van_t_hoff
)
from idaes.generic_models.properties.core.reactions.equilibrium_forms import (
    log_power_law_equil,
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
    "valid_phase_types": ["PT.aqueousPhase"],
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
    "type": "solute",
    "elements": ["Ca", "O", "H"]
}

reaction_data = [_reaction, _reaction]
component_data = [_component, _component]

Ca_thermo_config = {
    "components": {
        "Ca 2+": {
            "charge": 2,
            "type": "Component",
            "valid_phase_types": PT.aqueousPhase,
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
    "valid_phase_types": ["PT.aqueousPhase"],
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
            "equilibrium_constant": van_t_hoff,
            # "equilibrium_constant": gibbs_energy,
            "equilibrium_form": log_power_law_equil,
            "concentration_form": ConcentrationForm.molarity,
            "parameter_data": {
                "dh_rxn_ref": (14.9, pyunits.kJ / pyunits.mol),
                "ds_rxn_ref": (-148.1, pyunits.J / pyunits.mol / pyunits.K),
                # "T_eq_ref": (300, pyunits.K),
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
    },
    "name": "H2CO3_Ka2",
    "components": ["H2CO3", "Ka2"],
    "reactant_elements": ["H", "O", "C"],
    "type": "equilibrium"
}

# ==========================================================

from idaes.generic_models.properties.core.pure.Perrys import Perrys
from idaes.generic_models.properties.core.pure.NIST import NIST
# Import statements to be used in the starter config dict
from idaes.core import VaporPhase, AqueousPhase
from idaes.core.components import Solvent, Solute, Cation, Anion
from idaes.core.phases import PhaseType as PT
from idaes.generic_models.properties.core.phase_equil.forms import fugacity
from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.generic_models.properties.core.eos.ideal import Ideal
from idaes.generic_models.properties.core.phase_equil import SmoothVLE
from idaes.generic_models.properties.core.phase_equil.bubble_dew import IdealBubbleDew

recarbonation_thermo_config = {
    "components": {
        "H2O": {
            "type": Solvent,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Perrys,
            "enth_mol_liq_comp": Perrys,
            "cp_mol_liq_comp": Perrys,
            "entr_mol_liq_comp": Perrys,
            "enth_mol_ig_comp": NIST,
            "pressure_sat_comp": NIST,
            "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (18.0153, pyunits.g / pyunits.mol),
                "pressure_crit": (220.64e5, pyunits.Pa),
                "temperature_crit": (647, pyunits.K),
                # Comes from Perry's Handbook:  p. 2-98
                "dens_mol_liq_comp_coeff": {
                    "1": (5.459, pyunits.kmol * pyunits.m ** -3),
                    "2": (0.30542, pyunits.dimensionless),
                    "3": (647.13, pyunits.K),
                    "4": (0.081, pyunits.dimensionless),
                },
                "enth_mol_form_liq_comp_ref": (-285.830, pyunits.kJ / pyunits.mol),
                "enth_mol_form_vap_comp_ref": (0, pyunits.kJ / pyunits.mol),
                # Comes from Perry's Handbook:  p. 2-174
                "cp_mol_liq_comp_coeff": {
                    "1": (2.7637e5, pyunits.J / pyunits.kmol / pyunits.K),
                    "2": (-2.0901e3, pyunits.J / pyunits.kmol / pyunits.K ** 2),
                    "3": (8.125, pyunits.J / pyunits.kmol / pyunits.K ** 3),
                    "4": (-1.4116e-2, pyunits.J / pyunits.kmol / pyunits.K ** 4),
                    "5": (9.3701e-6, pyunits.J / pyunits.kmol / pyunits.K ** 5),
                },
                "cp_mol_ig_comp_coeff": {
                    "A": (30.09200, pyunits.J / pyunits.mol / pyunits.K),
                    "B": (
                        6.832514,
                        pyunits.J
                        * pyunits.mol ** -1
                        * pyunits.K ** -1
                        * pyunits.kiloK ** -1,
                    ),
                    "C": (
                        6.793435,
                        pyunits.J
                        * pyunits.mol ** -1
                        * pyunits.K ** -1
                        * pyunits.kiloK ** -2,
                    ),
                    "D": (
                        -2.534480,
                        pyunits.J
                        * pyunits.mol ** -1
                        * pyunits.K ** -1
                        * pyunits.kiloK ** -3,
                    ),
                    "E": (
                        0.082139,
                        pyunits.J
                        * pyunits.mol ** -1
                        * pyunits.K ** -1
                        * pyunits.kiloK ** 2,
                    ),
                    "F": (-250.8810, pyunits.kJ / pyunits.mol),
                    "G": (223.3967, pyunits.J / pyunits.mol / pyunits.K),
                    "H": (0, pyunits.kJ / pyunits.mol),
                },
                "entr_mol_form_liq_comp_ref": (
                    69.95,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
                "pressure_sat_comp_coeff": {
                    "A": (4.6543, None),  # [1], temperature range 255.9 K - 373 K
                    "B": (1435.264, pyunits.K),
                    "C": (-64.848, pyunits.K),
                },
            },
        },
        "CO2": {
            "type": Solute,
            "dens_mol_liq_comp": Perrys,
            "enth_mol_liq_comp": Perrys,
            "cp_mol_liq_comp": Perrys,
            "entr_mol_liq_comp": Perrys,
            "enth_mol_ig_comp": NIST,
            "pressure_sat_comp": NIST,
            "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
            "parameter_data": {
                "mw": (44.0095, pyunits.g / pyunits.mol),
                "pressure_crit": (73.825e5, pyunits.Pa),
                "temperature_crit": (304.23, pyunits.K),
                "dens_mol_liq_comp_coeff": {
                    "1": (2.768, pyunits.kmol * pyunits.m ** -3),
                    "2": (0.26212, None),
                    "3": (304.21, pyunits.K),
                    "4": (0.2908, None),
                },
                "cp_mol_ig_comp_coeff": {
                    "A": (24.99735, pyunits.J / pyunits.mol / pyunits.K),
                    "B": (
                        55.18696,
                        pyunits.J
                        * pyunits.mol ** -1
                        * pyunits.K ** -1
                        * pyunits.kiloK ** -1,
                    ),
                    "C": (
                        -33.69137,
                        pyunits.J
                        * pyunits.mol ** -1
                        * pyunits.K ** -1
                        * pyunits.kiloK ** -2,
                    ),
                    "D": (
                        7.948387,
                        pyunits.J
                        * pyunits.mol ** -1
                        * pyunits.K ** -1
                        * pyunits.kiloK ** -3,
                    ),
                    "E": (
                        -0.136638,
                        pyunits.J
                        * pyunits.mol ** -1
                        * pyunits.K ** -1
                        * pyunits.kiloK ** 2,
                    ),
                    "F": (-403.6075, pyunits.kJ / pyunits.mol),
                    "G": (228.2431, pyunits.J / pyunits.mol / pyunits.K),
                    "H": (0, pyunits.kJ / pyunits.mol),
                },
                "cp_mol_liq_comp_coeff": {
                    "1": (-8.3043e6, pyunits.J / pyunits.kmol / pyunits.K),
                    "2": (1.0437e5, pyunits.J / pyunits.kmol / pyunits.K ** 2),
                    "3": (4.333e2, pyunits.J / pyunits.kmol / pyunits.K ** 3),
                    "4": (6.0052e-1, pyunits.J / pyunits.kmol / pyunits.K ** 4),
                    "5": (0, pyunits.J / pyunits.kmol / pyunits.K ** 5),
                },
                "enth_mol_form_liq_comp_ref": (-285.83, pyunits.kJ / pyunits.mol),
                "enth_mol_form_vap_comp_ref": (-393.52, pyunits.kJ / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (0, pyunits.J / pyunits.K / pyunits.mol),
                "entr_mol_form_vap_comp_ref": (213.6, pyunits.J / pyunits.mol),
                "pressure_sat_comp_coeff": {
                    "A": (6.81228, None),
                    "B": (1301.679, pyunits.K),
                    "C": (-3.494, pyunits.K),
                },
            },
        },
        "H_+": {
            "type": Cation,
            "charge": 1,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Perrys,
            "enth_mol_liq_comp": Perrys,
            "cp_mol_liq_comp": Perrys,
            "entr_mol_liq_comp": Perrys,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (1.00784, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": {
                    "1": (5.459, pyunits.kmol * pyunits.m ** -3),
                    "2": (0.30542, pyunits.dimensionless),
                    "3": (647.13, pyunits.K),
                    "4": (0.081, pyunits.dimensionless),
                },
                "enth_mol_form_liq_comp_ref": (-230.000, pyunits.kJ / pyunits.mol),
                "cp_mol_liq_comp_coeff": {
                    "1": (2.7637e5, pyunits.J / pyunits.kmol / pyunits.K),
                    "2": (-2.0901e3, pyunits.J / pyunits.kmol / pyunits.K ** 2),
                    "3": (8.125, pyunits.J / pyunits.kmol / pyunits.K ** 3),
                    "4": (-1.4116e-2, pyunits.J / pyunits.kmol / pyunits.K ** 4),
                    "5": (9.3701e-6, pyunits.J / pyunits.kmol / pyunits.K ** 5),
                },
                "entr_mol_form_liq_comp_ref": (
                    -10.75,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
        },
        "OH_-": {
            "type": Anion,
            "charge": -1,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Perrys,
            "enth_mol_liq_comp": Perrys,
            "cp_mol_liq_comp": Perrys,
            "entr_mol_liq_comp": Perrys,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (17.008, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": {
                    "1": (5.459, pyunits.kmol * pyunits.m ** -3),
                    "2": (0.30542, pyunits.dimensionless),
                    "3": (647.13, pyunits.K),
                    "4": (0.081, pyunits.dimensionless),
                },
                "enth_mol_form_liq_comp_ref": (-230.000, pyunits.kJ / pyunits.mol),
                "cp_mol_liq_comp_coeff": {
                    "1": (2.7637e5, pyunits.J / pyunits.kmol / pyunits.K),
                    "2": (-2.0901e3, pyunits.J / pyunits.kmol / pyunits.K ** 2),
                    "3": (8.125, pyunits.J / pyunits.kmol / pyunits.K ** 3),
                    "4": (-1.4116e-2, pyunits.J / pyunits.kmol / pyunits.K ** 4),
                    "5": (9.3701e-6, pyunits.J / pyunits.kmol / pyunits.K ** 5),
                },
                "entr_mol_form_liq_comp_ref": (
                    -10.75,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
        },
        "H2CO3": {
            "type": Solute,
            "valid_phase_types": PT.aqueousPhase,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Perrys,
            "enth_mol_liq_comp": Perrys,
            "cp_mol_liq_comp": Perrys,
            "entr_mol_liq_comp": Perrys,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (62.03, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": {
                    "1": (5.4495, pyunits.kmol * pyunits.m ** -3),
                    "2": (0.427, pyunits.dimensionless),
                    "3": (429.69, pyunits.K),
                    "4": (0.259, pyunits.dimensionless),
                },
                "enth_mol_form_liq_comp_ref": (-699.7, pyunits.kJ / pyunits.mol),
                "cp_mol_liq_comp_coeff": {
                    "1": (135749.9, pyunits.J / pyunits.kmol / pyunits.K),
                    "2": (0, pyunits.J / pyunits.kmol / pyunits.K ** 2),
                    "3": (0, pyunits.J / pyunits.kmol / pyunits.K ** 3),
                    "4": (0, pyunits.J / pyunits.kmol / pyunits.K ** 4),
                    "5": (0, pyunits.J / pyunits.kmol / pyunits.K ** 5),
                },
                "entr_mol_form_liq_comp_ref": (
                    187,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
        },
        "HCO3_-": {
            "type": Anion,
            "charge": -1,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Perrys,
            "enth_mol_liq_comp": Perrys,
            "cp_mol_liq_comp": Perrys,
            "entr_mol_liq_comp": Perrys,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (61.0168, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": {
                    "1": (5.4495, pyunits.kmol * pyunits.m ** -3),
                    "2": (0.427, pyunits.dimensionless),
                    "3": (429.69, pyunits.K),
                    "4": (0.259, pyunits.dimensionless),
                },
                "enth_mol_form_liq_comp_ref": (-692, pyunits.kJ / pyunits.mol),
                "cp_mol_liq_comp_coeff": {
                    "1": (135749.9, pyunits.J / pyunits.kmol / pyunits.K),
                    "2": (0, pyunits.J / pyunits.kmol / pyunits.K ** 2),
                    "3": (0, pyunits.J / pyunits.kmol / pyunits.K ** 3),
                    "4": (0, pyunits.J / pyunits.kmol / pyunits.K ** 4),
                    "5": (0, pyunits.J / pyunits.kmol / pyunits.K ** 5),
                },
                "entr_mol_form_liq_comp_ref": (
                    91.2,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
        },
        "CO3_2-": {
            "type": Anion,
            "charge": -2,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Perrys,
            "enth_mol_liq_comp": Perrys,
            "cp_mol_liq_comp": Perrys,
            "entr_mol_liq_comp": Perrys,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (60.01, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": {
                    "1": (5.4495, pyunits.kmol * pyunits.m ** -3),
                    "2": (0.427, pyunits.dimensionless),
                    "3": (429.69, pyunits.K),
                    "4": (0.259, pyunits.dimensionless),
                },
                "enth_mol_form_liq_comp_ref": (-677.1, pyunits.J / pyunits.mol),
                "cp_mol_liq_comp_coeff": {
                    "1": (135749.9, pyunits.J / pyunits.kmol / pyunits.K),
                    "2": (0, pyunits.J / pyunits.kmol / pyunits.K ** 2),
                    "3": (0, pyunits.J / pyunits.kmol / pyunits.K ** 3),
                    "4": (0, pyunits.J / pyunits.kmol / pyunits.K ** 4),
                    "5": (0, pyunits.J / pyunits.kmol / pyunits.K ** 5),
                },
                "entr_mol_form_liq_comp_ref": (
                    -56.9,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
        },
    },
    "phases": {
        "Liq": {"type": AqueousPhase, "equation_of_state": Ideal},
        "Vap": {"type": VaporPhase, "equation_of_state": Ideal},
    },
    # Set base units of measurement
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },
    # Specifying state definition
    "state_definition": FTPx,
    "state_bounds": {
        "flow_mol": (0, 100, 20, pyunits.mol / pyunits.s),
        "temperature": (273.15, 300, 500, pyunits.K),
        "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
    },
    "pressure_ref": (101325, pyunits.Pa),
    "temperature_ref": (300, pyunits.K),
    # Defining phase equilibria
    "phases_in_equilibrium": [("Vap", "Liq")],
    "phase_equilibrium_state": {("Vap", "Liq"): SmoothVLE},
    "bubble_dew_method": IdealBubbleDew,
}

recarbonation_reaction_config = {
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },
    "equilibrium_reactions": {
        "H2O_Kw": {
            "stoichiometry": {
                ("Liq", "H2O"): -1,
                ("Liq", "H_+"): 1,
                ("Liq", "OH_-"): 1,
            },
            "heat_of_reaction": constant_dh_rxn,
            "equilibrium_constant": gibbs_energy,
            "equilibrium_form": log_power_law_equil,
            "concentration_form": ConcentrationForm.molarity,
            "parameter_data": {
                "dh_rxn_ref": (55.830, pyunits.kJ / pyunits.mol),
                "ds_rxn_ref": (-80.7, pyunits.J / pyunits.mol / pyunits.K),
                "T_eq_ref": (300, pyunits.K),
                # By default, reaction orders follow stoichiometry
                #    manually set reaction order here to override
                "reaction_order": {
                    ("Liq", "H2O"): 0,
                    ("Liq", "H_+"): 1,
                    ("Liq", "OH_-"): 1,
                },
            }
            # End parameter_data
        },
        # End R1
        "CO2_to_H2CO3": {
            "stoichiometry": {
                ("Liq", "H2O"): -1,
                ("Liq", "CO2"): -1,
                ("Liq", "H2CO3"): 1,
            },
            "heat_of_reaction": constant_dh_rxn,
            "equilibrium_constant": van_t_hoff,
            "equilibrium_form": log_power_law_equil,
            "concentration_form": ConcentrationForm.molarity,
            "parameter_data": {
                "dh_rxn_ref": (0, pyunits.kJ / pyunits.mol),
                "k_eq_ref": (1.7 * 10 ** -3, None),
                "T_eq_ref": (300, pyunits.K),
                # By default, reaction orders follow stoichiometry
                #    manually set reaction order here to override
                "reaction_order": {
                    ("Liq", "H2CO3"): 1,
                    ("Liq", "CO2"): -1,
                    ("Liq", "H2O"): 0,
                },
            }
            # End parameter_data
        },
        # End R2
        "H2CO3_Ka1": {
            "stoichiometry": {
                ("Liq", "H2CO3"): -1,
                ("Liq", "H_+"): 1,
                ("Liq", "HCO3_-"): 1,
            },
            "heat_of_reaction": constant_dh_rxn,
            "equilibrium_constant": gibbs_energy,
            "equilibrium_form": log_power_law_equil,
            "concentration_form": ConcentrationForm.molarity,
            "parameter_data": {
                "dh_rxn_ref": (7.7, pyunits.kJ / pyunits.mol),
                "ds_rxn_ref": (-95.8, pyunits.J / pyunits.mol / pyunits.K),
                "T_eq_ref": (300, pyunits.K),
                # By default, reaction orders follow stoichiometry
                #    manually set reaction order here to override
                "reaction_order": {
                    ("Liq", "H2CO3"): -1,
                    ("Liq", "H_+"): 1,
                    ("Liq", "HCO3_-"): 1,
                },
            }
            # End parameter_data
        },
        # End R3
        "H2CO3_Ka2": {
            "stoichiometry": {
                ("Liq", "HCO3_-"): -1,
                ("Liq", "H_+"): 1,
                ("Liq", "CO3_2-"): 1,
            },
            "heat_of_reaction": constant_dh_rxn,
            "equilibrium_constant": gibbs_energy,
            "equilibrium_form": log_power_law_equil,
            "concentration_form": ConcentrationForm.molarity,
            "parameter_data": {
                "dh_rxn_ref": (14.9, pyunits.kJ / pyunits.mol),
                "ds_rxn_ref": (-148.1, pyunits.J / pyunits.mol / pyunits.K),
                "T_eq_ref": (300, pyunits.K),
                # By default, reaction orders follow stoichiometry
                #    manually set reaction order here to override
                "reaction_order": {
                    ("Liq", "HCO3_-"): -1,
                    ("Liq", "H_+"): 1,
                    ("Liq", "CO3_2-"): 1,
                },
            }
            # End parameter_data
        }
        # End R4
    }
    # End equilibrium_reactions
}