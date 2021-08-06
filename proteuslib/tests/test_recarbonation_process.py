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
    This file is to run the example of running a phase change for the dissolution
    of CO2 into water.
"""
# =================== Import Statements ==============================
import pytest

# Import specific pyomo objects
from pyomo.environ import (
    ConcreteModel,
    SolverStatus,
    TerminationCondition,
    value,
    Suffix,
)

# Import pyomo methods to check the system units
from pyomo.util.check_units import assert_units_consistent
from pyomo.environ import units as pyunits

from idaes.core.util import scaling as iscale

from idaes.core.util import get_solver

# Import idaes methods to check the model during construction
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
)

# Import the idaes objects for Generic Properties and Reactions
from idaes.generic_models.properties.core.generic.generic_property import (
    GenericParameterBlock,
)
from idaes.generic_models.properties.core.generic.generic_reaction import (
    GenericReactionParameterBlock,
    ConcentrationForm,
)


from idaes.generic_models.properties.core.reactions.dh_rxn import constant_dh_rxn
from idaes.generic_models.properties.core.reactions.equilibrium_constant import (
    gibbs_energy,
    van_t_hoff,
)

# Import safe log power law equation
from idaes.generic_models.properties.core.reactions.equilibrium_forms import (
    log_power_law_equil,
)

# Import the idaes object for the EquilibriumReactor unit model
from idaes.generic_models.unit_models.equilibrium_reactor import EquilibriumReactor
from idaes.generic_models.properties.core.pure.Perrys import Perrys
from idaes.generic_models.properties.core.pure.NIST import NIST

# Import the core idaes objects for Flowsheets and types of balances
from idaes.core import FlowsheetBlock

# Import statements to be used in the starter config dict
from idaes.core import VaporPhase, AqueousPhase
from idaes.core.components import Solvent, Solute, Cation, Anion
from idaes.core.phases import PhaseType as PT
from idaes.generic_models.properties.core.phase_equil.forms import fugacity
from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.generic_models.properties.core.eos.ideal import Ideal
from idaes.generic_models.properties.core.phase_equil import SmoothVLE
from idaes.generic_models.properties.core.phase_equil.bubble_dew import IdealBubbleDew

# Import log10 function from pyomo
from pyomo.environ import log10

import idaes.logger as idaeslog

__author__ = "Srikanth Allu"

# Create a thermo_config dictionary
thermo_config = {
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
                    "1": (0.000789, pyunits.kmol * pyunits.m ** -3),
                    "2": (0.000956, pyunits.dimensionless),
                    "3": (500.78, pyunits.K),
                    "4": (0.94599, pyunits.dimensionless),
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

reaction_config = {
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
            "equilibrium_constant": van_t_hoff,
            "equilibrium_form": log_power_law_equil,
            "concentration_form": ConcentrationForm.moleFraction,
            "parameter_data": {
                "dh_rxn_ref": (55.830, pyunits.J/pyunits.mol),
                "k_eq_ref": (10**-14/55.2/55.2, pyunits.dimensionless),
                "T_eq_ref": (298, pyunits.K),
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
            "concentration_form": ConcentrationForm.moleFraction,
            "parameter_data": {
                "dh_rxn_ref": (0, pyunits.kJ / pyunits.mol),
                "k_eq_ref": (1.7*10**-3, pyunits.dimensionless),
                "T_eq_ref": (298, pyunits.K),
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
            "equilibrium_constant": van_t_hoff,
            "equilibrium_form": log_power_law_equil,
            "concentration_form": ConcentrationForm.moleFraction,
            "parameter_data": {
                "dh_rxn_ref": (7.7, pyunits.kJ/pyunits.mol),
                "k_eq_ref": (10**-6.35/55.2, pyunits.dimensionless),
                "T_eq_ref": (298, pyunits.K),
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
            "equilibrium_constant": van_t_hoff,
            "equilibrium_form": log_power_law_equil,
            "concentration_form": ConcentrationForm.moleFraction,
            "parameter_data": {
                "dh_rxn_ref": (14.9, pyunits.kJ/pyunits.mol),
                "k_eq_ref": (10**-10.33/55.2, pyunits.dimensionless),
                "T_eq_ref": (298, pyunits.K),
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

# Get default solver for testing
solver = get_solver()


class TestCarbonationProcess:
    @pytest.fixture(scope="class")
    def equilibrium_config(self):
        # Create a pyomo model object
        model = ConcreteModel()
        model.fs = FlowsheetBlock(default={"dynamic": False})
        model.fs.thermo_params = GenericParameterBlock(default=thermo_config)

        model.fs.rxn_params = GenericReactionParameterBlock(
            default={"property_package": model.fs.thermo_params, **reaction_config}
        )

        model.fs.unit = EquilibriumReactor(
            default={
                "property_package": model.fs.thermo_params,
                "reaction_package": model.fs.rxn_params,
                "has_rate_reactions": False,
                "has_equilibrium_reactions": True,
                "has_heat_transfer": False,
                "has_heat_of_reaction": False,
                "has_pressure_change": False,
            }
        )

        model.fs.unit.inlet.flow_mol.fix(10)
        model.fs.unit.inlet.pressure.fix(101325.0)
        model.fs.unit.inlet.temperature.fix(298.0)
        model.fs.unit.inlet.mole_frac_comp[0, "H_+"].fix(0.0)
        model.fs.unit.inlet.mole_frac_comp[0, "OH_-"].fix(0.0)
        model.fs.unit.inlet.mole_frac_comp[0, "H2CO3"].fix(0.0)
        model.fs.unit.inlet.mole_frac_comp[0, "HCO3_-"].fix(0.0)
        model.fs.unit.inlet.mole_frac_comp[0, "CO3_2-"].fix(0.0)
        model.fs.unit.inlet.mole_frac_comp[0, "CO2"].fix(0.0005)
        model.fs.unit.inlet.mole_frac_comp[0, "H2O"].fix(1 - 0.0005)

        print(model.fs.thermo_params.component_list)
        return model

    @pytest.mark.unit
    def test_build_model_equilibrium(self, equilibrium_config):
        model = equilibrium_config

        assert hasattr(model.fs.thermo_params, "component_list")
        assert len(model.fs.thermo_params.component_list) == 7
        assert "H2O" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.H2O, Solvent)
        assert "H_+" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("H_+"), Cation)
        assert "OH_-" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("OH_-"), Anion)

        assert hasattr(model.fs.thermo_params, "phase_list")
        assert len(model.fs.thermo_params.phase_list) == 2
        assert isinstance(model.fs.thermo_params.Liq, AqueousPhase)
        assert isinstance(model.fs.thermo_params.Vap, VaporPhase)

    @pytest.mark.unit
    def test_units_equilibrium(self, equilibrium_config):
        model = equilibrium_config
        assert_units_consistent(model)

    @pytest.mark.unit
    def test_dof_equilibrium(self, equilibrium_config):
        model = equilibrium_config
        assert degrees_of_freedom(model) == 0

    @pytest.mark.component
    def test_scaling_equilibrium(self, equilibrium_config):
        model = equilibrium_config

        # Equilibrium reactions have eps in the 'rxn_params'
        factor = 1e-4
        for rid in model.fs.rxn_params.equilibrium_reaction_idx:
            scale = value(model.fs.unit.control_volume.reactions[0.0].k_eq[rid].expr)
            # Want to set eps in some fashion similar to this
            if scale < 1e-16:
                model.fs.rxn_params.component("reaction_"+rid).eps.value = scale*factor
            else:
                model.fs.rxn_params.component("reaction_"+rid).eps.value = 1e-16*factor

        for i in model.fs.unit.control_volume.equilibrium_reaction_extent_index:
            scale = value(model.fs.unit.control_volume.reactions[0.0].k_eq[i[1]].expr)
            iscale.set_scaling_factor(model.fs.unit.control_volume.equilibrium_reaction_extent[0.0,i[1]], 10/scale)
            iscale.constraint_scaling_transform(model.fs.unit.control_volume.reactions[0.0].
                    equilibrium_constraint[i[1]], 0.1)

        # Next, try adding scaling for species
        min = 1e-6
        for i in model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp:
            # i[0] = phase, i[1] = species
            if model.fs.unit.inlet.mole_frac_comp[0, i[1]].value > min:
                scale = model.fs.unit.inlet.mole_frac_comp[0, i[1]].value
            else:
                scale = min
            iscale.set_scaling_factor(model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp[i[1]], 10/scale)
            iscale.set_scaling_factor(model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[i], 10/scale)
            iscale.set_scaling_factor(model.fs.unit.control_volume.properties_out[0.0].flow_mol_phase_comp[i], 10/scale)
            iscale.constraint_scaling_transform(
                model.fs.unit.control_volume.properties_out[0.0].component_flow_balances[i[1]], 10/scale)
            iscale.constraint_scaling_transform(model.fs.unit.control_volume.material_balances[0.0,i[1]], 10/scale)

        iscale.calculate_scaling_factors(model.fs.unit)

        assert isinstance(model.fs.unit.control_volume.scaling_factor, Suffix)

        assert isinstance(
            model.fs.unit.control_volume.properties_out[0.0].scaling_factor, Suffix
        )

        assert isinstance(
            model.fs.unit.control_volume.properties_in[0.0].scaling_factor, Suffix
        )

        # When using equilibrium reactions, there are another set of scaling factors calculated
        assert isinstance(
            model.fs.unit.control_volume.reactions[0.0].scaling_factor, Suffix
        )

    @pytest.mark.component
    def test_initialize_solver(self, equilibrium_config):
        model = equilibrium_config
        solver.options["bound_push"] = 1e-10
        solver.options["mu_init"] = 1e-6
        model.fs.unit.initialize(optarg=solver.options, outlvl=idaeslog.DEBUG)
        assert degrees_of_freedom(model) == 0

    @pytest.mark.component
    def test_solve_equilibrium(self, equilibrium_config):
        model = equilibrium_config
        solver.options["max_iter"] = 100
        results = solver.solve(model)
        print(results.solver.termination_condition)
        assert results.solver.termination_condition == TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.component
    def test_solution_equilibrium(self, equilibrium_config):
        model = equilibrium_config

        assert pytest.approx(298, rel=1e-5) == value(
            model.fs.unit.outlet.temperature[0]
        )
        assert pytest.approx(10, rel=1e-5) == value(model.fs.unit.outlet.flow_mol[0])
        assert pytest.approx(101325, rel=1e-5) == value(
            model.fs.unit.outlet.pressure[0]
        )

        total_molar_density = (
            value(
                model.fs.unit.control_volume.properties_out[0.0].dens_mol_phase["Liq"]
            )
            / 1000
        )
        assert pytest.approx(55.165246, rel=1e-5) == total_molar_density
        pH = -value(
            log10(model.fs.unit.outlet.mole_frac_comp[0, "H_+"] * total_molar_density)
        )
        pOH = -value(
            log10(model.fs.unit.outlet.mole_frac_comp[0, "OH_-"] * total_molar_density)
        )
        assert pytest.approx(5.33989, rel=1e-4) == pH
        assert pytest.approx(8.660655, rel=1e-4) == pOH

        CO2_sorbed = value(
            model.fs.unit.control_volume.properties_out[0.0].conc_mol_phase_comp[
                ("Liq", "CO2")
            ]
        )
        assert pytest.approx(27.5312738, rel=1e-5) == CO2_sorbed
