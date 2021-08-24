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

[1] Local Composition Model for Excess Gibbs Energy of Electrolyte Systems, Pt 1.
Chen, C.-C., Britt, H.I., Boston, J.F., Evans, L.B.,
AIChE Journal, 1982, Vol. 28(4), pgs. 588-596

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
    Unsymmetric
from idaes.generic_models.properties.core.generic.generic_property import (
        StateIndex)
from idaes.generic_models.properties.core.state_definitions import FTPx
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
                    "reference_state": Unsymmetric}}},
    "base_units": {"time": pyunits.s,
                   "length": pyunits.m,
                   "mass": pyunits.kg,
                   "amount": pyunits.mol,
                   "temperature": pyunits.K},
    "state_definition": FTPx,
    "state_components": StateIndex.true,
    "pressure_ref": 1e5,
    "temperature_ref": 298.15,
    "parameter_data": {
        "Liq_tau": {  # Table 1 [1], no parameters for CaSO4
            ("H2O", "Na_+, Cl_-"): 8.885,
            ("Na_+, Cl_-", "H2O"): -4.549,
            ("H2O", "Ca_2+, Cl_-"): 11.396,
            ("Ca_2+, Cl_-", "H2O"): -6.218,
            ("H2O", "Mg_2+, Cl_-"): 11.579,
            ("Mg_2+, Cl_-", "H2O"): -6.338,
            ("H2O", "Na_+, SO4_2-"): 8.389,
            ("Na_+, SO4_2-", "H2O"): -4.539,
            ("H2O", "Mg_2+, SO4_2-"): 11.346,
            ("Mg_2+, SO4_2-", "H2O"): -6.862}}}
