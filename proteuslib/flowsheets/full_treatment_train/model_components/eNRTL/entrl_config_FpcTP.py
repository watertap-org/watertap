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
eNRTL property configuration dicts for synthetic hard water

References:

[1] Islam, R.I., et al., Molecular thermodynamics for scaling prediction: Case
of membrane distillation, Separation and Purification Technology, 2021,
Vol. 276.

Author: Andrew Lee
"""

from pyomo.environ import Param, units as pyunits

from idaes.core import (AqueousPhase,
                        Solvent,
                        Apparent,
                        Anion,
                        Cation)
from idaes.generic_models.properties.core.eos.enrtl import ENRTL
from idaes.generic_models.properties.core.eos.enrtl_reference_states import \
    Symmetric, Unsymmetric
from idaes.generic_models.properties.core.generic.generic_property import (
        StateIndex)
from idaes.generic_models.properties.core.state_definitions import FpcTP
from idaes.generic_models.properties.core.pure.electrolyte import \
    relative_permittivity_constant


class ConstantVolMol():
    def build_parameters(b):
        b.vol_mol_pure = Param(initialize=18e-6,
                               units=pyunits.m**3/pyunits.mol,
                               mutable=True)

    def return_expression(b, cobj, T):
        return cobj.vol_mol_pure


configuration = {
    "components": {
        "H2O": {"type": Solvent,
                "vol_mol_liq_comp": ConstantVolMol,
                "relative_permittivity_liq_comp":
                    relative_permittivity_constant,
                "parameter_data": {
                    "mw": (18E-3, pyunits.kg/pyunits.mol),
                    "relative_permittivity_liq_comp": 78.54}},
        "NaCl": {"type": Apparent,
                 "dissociation_species": {"Na_+": 1, "Cl_-": 1}},
        "Na2SO4": {"type": Apparent,
                   "dissociation_species": {"Na_+": 2, "SO4_2-": 1}},
        "CaCl2": {"type": Apparent,
                  "dissociation_species": {"Ca_2+": 1, "Cl_-": 2}},
        "CaSO4": {"type": Apparent,
                  "dissociation_species": {"Ca_2+": 1, "SO4_2-": 1}},
        "MgCl2": {"type": Apparent,
                  "dissociation_species": {"Mg_2+": 1, "Cl_-": 2}},
        "MgSO4": {"type": Apparent,
                  "dissociation_species": {"Mg_2+": 1, "SO4_2-": 1}},
        "Na_+": {"type": Cation,
                 "charge": +1},
        "Ca_2+": {"type": Cation,
                  "charge": +2},
        "Mg_2+": {"type": Cation,
                  "charge": +2},
        "Cl_-": {"type": Anion,
                 "charge": -1},
        "SO4_2-": {"type": Anion,
                   "charge": -2}},
    "phases": {
        "Liq": {"type": AqueousPhase,
                "equation_of_state": ENRTL,
                "equation_of_state_options": {
                    "reference_state": Symmetric}}},
    "base_units": {"time": pyunits.s,
                   "length": pyunits.m,
                   "mass": pyunits.kg,
                   "amount": pyunits.mol,
                   "temperature": pyunits.K},
    "state_definition": FpcTP,
    "state_components": StateIndex.true,
    "pressure_ref": 1e5,
    "temperature_ref": 298.15,
    "parameter_data": {
        "Liq_tau": {  # Table 1 [1]
            ('H2O', 'Na_+, Cl_-'): 8.866,
            ('Na_+, Cl_-', 'H2O'): -4.451,
            ('H2O', 'Ca_2+, Cl_-'): 10.478,
            ('Ca_2+, Cl_-', 'H2O'): -5.231,
            ('H2O', 'Mg_2+, Cl_-'): 10.854,
            ('Mg_2+, Cl_-', 'H2O'): -5.409,
            ('H2O', 'Na_+, SO4_2-'): 8.012,
            ('Na_+, SO4_2-', 'H2O'): -3.903,
            ('H2O', 'Ca_2+, SO4_2-'): 6.932,
            ('Ca_2+, SO4_2-', 'H2O'): -3.466,
            ('H2O', 'Mg_2+, SO4_2-'): 8.808,
            ('Mg_2+, SO4_2-', 'H2O'): -4.383,
            ('Na_+, Cl_-', 'Ca_2+, Cl_-'): -0.468,
            ('Ca_2+, Cl_-', 'Na_+, Cl_-'): 0.41,
            ('Na_+, Cl_-', 'Mg_2+, Cl_-'): -0.328,
            ('Mg_2+, Cl_-', 'Na_+, Cl_-'): -0.981,
            ('Ca_2+, Cl_-', 'Mg_2+, Cl_-'): 0.22,
            ('Mg_2+, Cl_-', 'Ca_2+, Cl_-'): 0.322,
            ('Na_+, SO4_2-', 'Ca_2+, SO4_2-'): -0.761,
            ('Ca_2+, SO4_2-', 'Na_+, SO4_2-'): 0.368,
            ('Na_+, SO4_2-', 'Mg_2+, SO4_2-'): -0.327,
            ('Mg_2+, SO4_2-', 'Na_+, SO4_2-'): 0.799,
            ('Ca_2+, SO4_2-', 'Mg_2+, SO4_2-'): 0,
            ('Mg_2+, SO4_2-', 'Ca_2+, SO4_2-'): 0.383,
            ('Na_+, Cl_-', 'Na_+, SO4_2-'): 0.803,
            ('Na_+, SO4_2-', 'Na_+, Cl_-'): -0.634,
            ('Ca_2+, Cl_-', 'Ca_2+, SO4_2-'): 0,
            ('Ca_2+, SO4_2-', 'Ca_2+, Cl_-'): -0.264,
            ('Mg_2+, Cl_-', 'Mg_2+, SO4_2-'): -0.707,
            ('Mg_2+, SO4_2-', 'Mg_2+, Cl_-'): -0.841
            }},
    "default_scaling_factors": {("flow_mol_phase_comp", ("Liq", "Na_+")): 1e1,
                                ("flow_mol_phase_comp", ("Liq", "Ca_2+")): 1e3,
                                ("flow_mol_phase_comp", ("Liq", "Mg_2+")): 1e2,
                                ("flow_mol_phase_comp", ("Liq", "SO4_2-")): 1e2,
                                ("flow_mol_phase_comp", ("Liq", "Cl_-")): 1e1,
                                ("flow_mol_phase_comp", ("Liq", "H2O")): 1e-1,
                                ("mole_frac_comp", "Na_+"): 1e2,
                                ("mole_frac_comp", "Ca_2+"): 1e4,
                                ("mole_frac_comp", "Mg_2+"): 1e3,
                                ("mole_frac_comp", "SO4_2-"): 1e3,
                                ("mole_frac_comp", "Cl_-"): 1e2,
                                ("mole_frac_comp", "H2O"): 1,
                                ("mole_frac_phase_comp", ("Liq", "Na_+")): 1e2,
                                ("mole_frac_phase_comp", ("Liq", "Ca_2+")): 1e4,
                                ("mole_frac_phase_comp", ("Liq", "Mg_2+")): 1e3,
                                ("mole_frac_phase_comp", ("Liq", "SO4_2-")): 1e3,
                                ("mole_frac_phase_comp", ("Liq", "Cl_-")): 1e2,
                                ("mole_frac_phase_comp", ("Liq", "H2O")): 1,
                                ("flow_mol_phase_comp_apparent", ("Liq", "NaCl")): 1e1,
                                ("flow_mol_phase_comp_apparent", ("Liq", "Na2SO4")): 1e2,
                                ("flow_mol_phase_comp_apparent", ("Liq", "CaCl2")): 1e2,
                                ("flow_mol_phase_comp_apparent", ("Liq", "CaSO4")): 1e3,
                                ("flow_mol_phase_comp_apparent", ("Liq", "MgCl2")): 1e2,
                                ("flow_mol_phase_comp_apparent", ("Liq", "MgSO4")): 1e3,
                                ("flow_mol_phase_comp_apparent", ("Liq", "H2O")): 1e-1,
                                ("mole_frac_phase_comp_apparent", ("Liq", "NaCl")): 1e3,  # TODO: these seem to be 1 orders of magnitude too low
                                ("mole_frac_phase_comp_apparent", ("Liq", "Na2SO4")): 1e4,
                                ("mole_frac_phase_comp_apparent", ("Liq", "CaCl2")): 1e4,
                                ("mole_frac_phase_comp_apparent", ("Liq", "CaSO4")): 1e5,
                                ("mole_frac_phase_comp_apparent", ("Liq", "MgCl2")): 1e4,
                                ("mole_frac_phase_comp_apparent", ("Liq", "MgSO4")): 1e5,
                                ("mole_frac_phase_comp_apparent", ("Liq", "H2O")): 1,
                                }
}
