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
    Ideal NaOCl Chlorination posttreatment process

    This will build an ideal NaOCl pretreatment block as a combination of a
    Mixer (where NaOCl is added) and an EquilibriumReactor (where pH and free
    chlorine is calculated)::

                    NaOCl stream
                        |
                        V
       inlet stream ---> [Mixer] --- outlet stream ---> [EquilibriumReactor] ---> exit stream (to distribution)
"""

# Importing the object for units from pyomo
from pyomo.environ import units as pyunits, assert_optimal_termination, NonNegativeReals

# Imports from idaes core
from idaes.core import AqueousPhase
from idaes.core.base.components import Solvent, Solute, Cation, Anion
from idaes.core.base.phases import PhaseType as PT
from idaes.core.solvers import get_solver

# Imports from idaes generic models
import idaes.models.properties.modular_properties.pure.Perrys as Perrys
from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.eos.ideal import Ideal

# Importing the enum for concentration unit basis used in the 'get_concentration_term' function
from idaes.models.properties.modular_properties.base.generic_reaction import (
    ConcentrationForm,
)

# Import the object/function for heat of reaction
from idaes.models.properties.modular_properties.reactions.dh_rxn import constant_dh_rxn

# Import safe log power law equation
from idaes.models.properties.modular_properties.reactions.equilibrium_forms import (
    log_power_law_equil,
)

# Import built-in van't Hoff function
from idaes.models.properties.modular_properties.reactions.equilibrium_constant import (
    van_t_hoff,
)

# Import specific pyomo objects
from pyomo.environ import (
    ConcreteModel,
    Var,
    Constraint,
    SolverStatus,
    TerminationCondition,
    TransformationFactory,
    value,
    Suffix,
    Expression,
)

from pyomo.network import Arc

from idaes.core.util import scaling as iscale
from idaes.core.util.initialization import fix_state_vars, revert_state_vars

# Import pyomo methods to check the system units
from pyomo.util.check_units import assert_units_consistent


from watertap.examples.flowsheets.full_treatment_train.util import (
    solve_block,
    check_dof,
)
from watertap.examples.flowsheets.full_treatment_train.model_components import (
    property_models,
)
from idaes.core.solvers import get_solver

# Import the idaes objects for Generic Properties and Reactions
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models.properties.modular_properties.base.generic_reaction import (
    GenericReactionParameterBlock,
)

# Import the idaes object for the EquilibriumReactor unit model
from idaes.models.unit_models.equilibrium_reactor import EquilibriumReactor

# Import the Mixer unit model
from idaes.models.unit_models import Mixer
from watertap.examples.flowsheets.full_treatment_train.flowsheet_components import (
    costing,
)

from idaes.models.unit_models.translator import Translator

# Import the core idaes objects for Flowsheets and types of balances
from idaes.core import FlowsheetBlock

# Import log10 function from pyomo
from pyomo.environ import log10

import idaes.logger as idaeslog

# Grab the scaling utilities
from watertap.examples.flowsheets.full_treatment_train.electrolyte_scaling_utils import (
    approximate_chemical_state_args,
    calculate_chemical_scaling_factors,
    calculate_chemical_scaling_factors_for_energy_balances,
)

from watertap.examples.flowsheets.full_treatment_train.chemical_flowsheet_util import (
    set_H2O_molefraction,
    zero_out_non_H2O_molefractions,
    fix_all_molefractions,
    unfix_all_molefractions,
    seq_decomp_initializer,
)

from watertap.examples.flowsheets.full_treatment_train.flowsheet_components import (
    desalination,
)
from idaes.core.util.initialization import propagate_state

__author__ = "Austin Ladshaw"

# Configuration dictionary
ideal_naocl_thermo_config = {
    "components": {
        "H2O": {
            "type": Solvent,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Perrys,
            "enth_mol_liq_comp": Perrys,
            "cp_mol_liq_comp": Perrys,
            "entr_mol_liq_comp": Perrys,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (18.0153, pyunits.g / pyunits.mol),
                "pressure_crit": (220.64e5, pyunits.Pa),
                "temperature_crit": (647, pyunits.K),
                # Comes from Perry's Handbook:  p. 2-98
                "dens_mol_liq_comp_coeff": {
                    "1": (5.459, pyunits.kmol * pyunits.m**-3),
                    "2": (0.30542, pyunits.dimensionless),
                    "3": (647.13, pyunits.K),
                    "4": (0.081, pyunits.dimensionless),
                },
                "enth_mol_form_liq_comp_ref": (-285.830, pyunits.kJ / pyunits.mol),
                "enth_mol_form_vap_comp_ref": (0, pyunits.kJ / pyunits.mol),
                # Comes from Perry's Handbook:  p. 2-174
                "cp_mol_liq_comp_coeff": {
                    "1": (2.7637e5, pyunits.J / pyunits.kmol / pyunits.K),
                    "2": (-2.0901e3, pyunits.J / pyunits.kmol / pyunits.K**2),
                    "3": (8.125, pyunits.J / pyunits.kmol / pyunits.K**3),
                    "4": (-1.4116e-2, pyunits.J / pyunits.kmol / pyunits.K**4),
                    "5": (9.3701e-6, pyunits.J / pyunits.kmol / pyunits.K**5),
                },
                "cp_mol_ig_comp_coeff": {
                    "A": (30.09200, pyunits.J / pyunits.mol / pyunits.K),
                    "B": (
                        6.832514,
                        pyunits.J
                        * pyunits.mol**-1
                        * pyunits.K**-1
                        * pyunits.kiloK**-1,
                    ),
                    "C": (
                        6.793435,
                        pyunits.J
                        * pyunits.mol**-1
                        * pyunits.K**-1
                        * pyunits.kiloK**-2,
                    ),
                    "D": (
                        -2.534480,
                        pyunits.J
                        * pyunits.mol**-1
                        * pyunits.K**-1
                        * pyunits.kiloK**-3,
                    ),
                    "E": (
                        0.082139,
                        pyunits.J
                        * pyunits.mol**-1
                        * pyunits.K**-1
                        * pyunits.kiloK**2,
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
            # End parameter_data
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
                    "1": (5.459, pyunits.kmol * pyunits.m**-3),
                    "2": (0.30542, pyunits.dimensionless),
                    "3": (647.13, pyunits.K),
                    "4": (0.081, pyunits.dimensionless),
                },
                "enth_mol_form_liq_comp_ref": (-230.000, pyunits.kJ / pyunits.mol),
                "cp_mol_liq_comp_coeff": {
                    "1": (2.7637e5, pyunits.J / pyunits.kmol / pyunits.K),
                    "2": (-2.0901e3, pyunits.J / pyunits.kmol / pyunits.K**2),
                    "3": (8.125, pyunits.J / pyunits.kmol / pyunits.K**3),
                    "4": (-1.4116e-2, pyunits.J / pyunits.kmol / pyunits.K**4),
                    "5": (9.3701e-6, pyunits.J / pyunits.kmol / pyunits.K**5),
                },
                "entr_mol_form_liq_comp_ref": (
                    -10.75,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
            # End parameter_data
        },
        "Na_+": {
            "type": Cation,
            "charge": 1,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Perrys,
            "enth_mol_liq_comp": Perrys,
            "cp_mol_liq_comp": Perrys,
            "entr_mol_liq_comp": Perrys,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (22.989769, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": {
                    "1": (5.252, pyunits.kmol * pyunits.m**-3),
                    "2": (0.347, pyunits.dimensionless),
                    "3": (1595.8, pyunits.K),
                    "4": (0.6598, pyunits.dimensionless),
                },
                "enth_mol_form_liq_comp_ref": (-240.1, pyunits.J / pyunits.mol),
                "cp_mol_liq_comp_coeff": {
                    "1": (167039, pyunits.J / pyunits.kmol / pyunits.K),
                    "2": (0, pyunits.J / pyunits.kmol / pyunits.K**2),
                    "3": (0, pyunits.J / pyunits.kmol / pyunits.K**3),
                    "4": (0, pyunits.J / pyunits.kmol / pyunits.K**4),
                    "5": (0, pyunits.J / pyunits.kmol / pyunits.K**5),
                },
                "entr_mol_form_liq_comp_ref": (59, pyunits.J / pyunits.K / pyunits.mol),
            },
            # End parameter_data
        },
        "Cl_-": {
            "type": Anion,
            "charge": -1,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Perrys,
            "enth_mol_liq_comp": Perrys,
            "cp_mol_liq_comp": Perrys,
            "entr_mol_liq_comp": Perrys,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (35.453, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": {
                    "1": (4.985, pyunits.kmol * pyunits.m**-3),
                    "2": (0.36, pyunits.dimensionless),
                    "3": (1464.06, pyunits.K),
                    "4": (0.739, pyunits.dimensionless),
                },
                "enth_mol_form_liq_comp_ref": (-167.2, pyunits.kJ / pyunits.mol),
                "cp_mol_liq_comp_coeff": {
                    "1": (83993.8, pyunits.J / pyunits.kmol / pyunits.K),
                    "2": (0, pyunits.J / pyunits.kmol / pyunits.K**2),
                    "3": (0, pyunits.J / pyunits.kmol / pyunits.K**3),
                    "4": (0, pyunits.J / pyunits.kmol / pyunits.K**4),
                    "5": (0, pyunits.J / pyunits.kmol / pyunits.K**5),
                },
                "entr_mol_form_liq_comp_ref": (
                    56.5,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
            # End parameter_data
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
                    "1": (5.459, pyunits.kmol * pyunits.m**-3),
                    "2": (0.30542, pyunits.dimensionless),
                    "3": (647.13, pyunits.K),
                    "4": (0.081, pyunits.dimensionless),
                },
                "enth_mol_form_liq_comp_ref": (-230.000, pyunits.kJ / pyunits.mol),
                "cp_mol_liq_comp_coeff": {
                    "1": (2.7637e5, pyunits.J / pyunits.kmol / pyunits.K),
                    "2": (-2.0901e3, pyunits.J / pyunits.kmol / pyunits.K**2),
                    "3": (8.125, pyunits.J / pyunits.kmol / pyunits.K**3),
                    "4": (-1.4116e-2, pyunits.J / pyunits.kmol / pyunits.K**4),
                    "5": (9.3701e-6, pyunits.J / pyunits.kmol / pyunits.K**5),
                },
                "entr_mol_form_liq_comp_ref": (
                    -10.75,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
            # End parameter_data
        },
        "HOCl": {
            "type": Solute,
            "valid_phase_types": PT.aqueousPhase,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Perrys,
            "enth_mol_liq_comp": Perrys,
            "cp_mol_liq_comp": Perrys,
            "entr_mol_liq_comp": Perrys,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (52.46, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": {
                    "1": (4.985, pyunits.kmol * pyunits.m**-3),
                    "2": (0.36, pyunits.dimensionless),
                    "3": (1464.06, pyunits.K),
                    "4": (0.739, pyunits.dimensionless),
                },
                "enth_mol_form_liq_comp_ref": (-120.9, pyunits.kJ / pyunits.mol),
                "cp_mol_liq_comp_coeff": {
                    "1": (83993.8, pyunits.J / pyunits.kmol / pyunits.K),
                    "2": (0, pyunits.J / pyunits.kmol / pyunits.K**2),
                    "3": (0, pyunits.J / pyunits.kmol / pyunits.K**3),
                    "4": (0, pyunits.J / pyunits.kmol / pyunits.K**4),
                    "5": (0, pyunits.J / pyunits.kmol / pyunits.K**5),
                },
                "entr_mol_form_liq_comp_ref": (
                    142,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
            # End parameter_data
        },
        "OCl_-": {
            "type": Anion,
            "charge": -1,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Perrys,
            "enth_mol_liq_comp": Perrys,
            "cp_mol_liq_comp": Perrys,
            "entr_mol_liq_comp": Perrys,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (51.46, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": {
                    "1": (4.985, pyunits.kmol * pyunits.m**-3),
                    "2": (0.36, pyunits.dimensionless),
                    "3": (1464.06, pyunits.K),
                    "4": (0.739, pyunits.dimensionless),
                },
                "enth_mol_form_liq_comp_ref": (-107.1, pyunits.kJ / pyunits.mol),
                "cp_mol_liq_comp_coeff": {
                    "1": (83993.8, pyunits.J / pyunits.kmol / pyunits.K),
                    "2": (0, pyunits.J / pyunits.kmol / pyunits.K**2),
                    "3": (0, pyunits.J / pyunits.kmol / pyunits.K**3),
                    "4": (0, pyunits.J / pyunits.kmol / pyunits.K**4),
                    "5": (0, pyunits.J / pyunits.kmol / pyunits.K**5),
                },
                "entr_mol_form_liq_comp_ref": (42, pyunits.J / pyunits.K / pyunits.mol),
            },
            # End parameter_data
        },
    },
    # End Component list
    "phases": {
        "Liq": {"type": AqueousPhase, "equation_of_state": Ideal},
    },
    "state_definition": FTPx,
    "state_bounds": {
        "flow_mol": (0, 50, 100),
        "temperature": (273.15, 300, 650),
        "pressure": (5e4, 1e5, 1e6),
    },
    "pressure_ref": 1e5,
    "temperature_ref": 300,
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },
}
# End simple_naocl_thermo_config definition

# This config is REQUIRED to use EquilibriumReactor even if we have no equilibrium reactions
ideal_naocl_reaction_config = {
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
                "dh_rxn_ref": (55.830, pyunits.J / pyunits.mol),
                "k_eq_ref": (10**-14 / 55.2 / 55.2, pyunits.dimensionless),
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
        "HOCl_Ka": {
            "stoichiometry": {
                ("Liq", "HOCl"): -1,
                ("Liq", "H_+"): 1,
                ("Liq", "OCl_-"): 1,
            },
            "heat_of_reaction": constant_dh_rxn,
            "equilibrium_constant": van_t_hoff,
            "equilibrium_form": log_power_law_equil,
            "concentration_form": ConcentrationForm.moleFraction,
            "parameter_data": {
                "dh_rxn_ref": (13.8, pyunits.J / pyunits.mol),
                "k_eq_ref": (10**-7.6422 / 55.2, pyunits.dimensionless),
                "T_eq_ref": (298, pyunits.K),
                # By default, reaction orders follow stoichiometry
                #    manually set reaction order here to override
                "reaction_order": {
                    ("Liq", "HOCl"): -1,
                    ("Liq", "H_+"): 1,
                    ("Liq", "OCl_-"): 1,
                },
            }
            # End parameter_data
        }
        # End R2
    }
    # End equilibrium_reactions
}
# End reaction_config definition

# Get default solver for testing
solver = get_solver(options={"tol": 1e-11})


def build_ideal_naocl_prop(model):
    model.fs.ideal_naocl_thermo_params = GenericParameterBlock(
        **ideal_naocl_thermo_config
    )
    model.fs.ideal_naocl_rxn_params = GenericReactionParameterBlock(
        property_package=model.fs.ideal_naocl_thermo_params,
        **ideal_naocl_reaction_config
    )


def build_ideal_naocl_mixer_unit(model):
    model.fs.ideal_naocl_mixer_unit = Mixer(
        property_package=model.fs.ideal_naocl_thermo_params,
        inlet_list=["inlet_stream", "naocl_stream"],
    )

    # add new constraint for dosing rate (deactivate constraint for OCl_-)
    dr = value(
        model.fs.ideal_naocl_mixer_unit.naocl_stream.flow_mol[0]
        * model.fs.ideal_naocl_mixer_unit.naocl_stream.mole_frac_comp[0, "OCl_-"]
        * model.fs.ideal_naocl_thermo_params.get_component("OCl_-").mw
    )
    model.fs.ideal_naocl_mixer_unit.dosing_rate = Var(
        initialize=dr,
        domain=NonNegativeReals,
        units=pyunits.kg / pyunits.s,
    )

    def _dosing_rate_cons(blk):
        return blk.dosing_rate == (
            blk.naocl_stream.flow_mol[0]
            * blk.naocl_stream.mole_frac_comp[0, "OCl_-"]
            * blk.naocl_stream_state[0].params.get_component("OCl_-").mw
        ) + (
            blk.naocl_stream.flow_mol[0]
            * blk.naocl_stream.mole_frac_comp[0, "Na_+"]
            * blk.naocl_stream_state[0].params.get_component("Na_+").mw
        )

    model.fs.ideal_naocl_mixer_unit.dosing_cons = Constraint(rule=_dosing_rate_cons)


def build_ideal_naocl_chlorination_unit(model):
    model.fs.ideal_naocl_chlorination_unit = EquilibriumReactor(
        property_package=model.fs.ideal_naocl_thermo_params,
        reaction_package=model.fs.ideal_naocl_rxn_params,
        has_rate_reactions=False,
        has_equilibrium_reactions=True,
        has_heat_transfer=False,
        has_heat_of_reaction=False,
        has_pressure_change=False,
    )

    # new var includes an initial calculation (will be overwritten later)
    fc = (
        model.fs.ideal_naocl_chlorination_unit.inlet.mole_frac_comp[0, "HOCl"].value
        * 55.6
    )
    fc += (
        model.fs.ideal_naocl_chlorination_unit.inlet.mole_frac_comp[0, "OCl_-"].value
        * 55.6
    )
    fc = fc * 70900

    model.fs.ideal_naocl_chlorination_unit.free_chlorine = Var(initialize=fc)

    def _free_chlorine_cons(blk):
        return (
            blk.free_chlorine
            == (
                (
                    blk.control_volume.properties_out[0.0].conc_mol_phase_comp[
                        "Liq", "HOCl"
                    ]
                    / 1000
                )
                + (
                    blk.control_volume.properties_out[0.0].conc_mol_phase_comp[
                        "Liq", "OCl_-"
                    ]
                    / 1000
                )
            )
            * 70900
        )

    model.fs.ideal_naocl_chlorination_unit.chlorine_cons = Constraint(
        rule=_free_chlorine_cons
    )


def set_ideal_naocl_mixer_inlets(
    model,
    dosing_rate_of_NaOCl_mg_per_s=0.4,
    inlet_water_density_kg_per_L=1,
    inlet_temperature_K=298,
    inlet_pressure_Pa=101325,
    inlet_flow_mol_per_s=10,
):

    # inlet stream
    model.fs.ideal_naocl_mixer_unit.inlet_stream.flow_mol[0].set_value(
        inlet_flow_mol_per_s
    )
    model.fs.ideal_naocl_mixer_unit.inlet_stream.pressure[0].set_value(
        inlet_pressure_Pa
    )
    model.fs.ideal_naocl_mixer_unit.inlet_stream.temperature[0].set_value(
        inlet_temperature_K
    )

    zero_out_non_H2O_molefractions(model.fs.ideal_naocl_mixer_unit.inlet_stream)
    set_H2O_molefraction(model.fs.ideal_naocl_mixer_unit.inlet_stream)

    # naocl stream
    model.fs.ideal_naocl_mixer_unit.naocl_stream.pressure[0].set_value(
        inlet_pressure_Pa
    )
    model.fs.ideal_naocl_mixer_unit.naocl_stream.temperature[0].set_value(
        inlet_temperature_K
    )
    # Use given dosing rate value to estimate OCl_- molefraction and flow rate for naocl stream
    zero_out_non_H2O_molefractions(model.fs.ideal_naocl_mixer_unit.naocl_stream)
    model.fs.ideal_naocl_mixer_unit.naocl_stream.mole_frac_comp[0, "OCl_-"].set_value(
        0.5
    )
    model.fs.ideal_naocl_mixer_unit.naocl_stream.mole_frac_comp[0, "Na_+"].set_value(
        0.5
    )
    set_H2O_molefraction(model.fs.ideal_naocl_mixer_unit.naocl_stream)
    flow_of_naocl = (
        dosing_rate_of_NaOCl_mg_per_s
        / model.fs.ideal_naocl_mixer_unit.naocl_stream.mole_frac_comp[0, "OCl_-"].value
        / 74.44
        / 1000
    )
    model.fs.ideal_naocl_mixer_unit.naocl_stream.flow_mol[0].set_value(flow_of_naocl)

    model.fs.ideal_naocl_mixer_unit.dosing_rate.set_value(dosing_rate_of_NaOCl_mg_per_s)


def set_ideal_naocl_chlorination_inlets(
    model,
    mg_per_L_NaOCl_added=2,
    inlet_water_density_kg_per_L=1,
    inlet_temperature_K=298,
    inlet_pressure_Pa=101325,
    inlet_flow_mol_per_s=10,
):

    zero_out_non_H2O_molefractions(model.fs.ideal_naocl_chlorination_unit.inlet)

    total_molar_density = inlet_water_density_kg_per_L / 18 * 1000  # mol/L

    # Free Chlorine (mg-Cl2/L) = total_chlorine_inlet (mol/L) * 70,900
    #       Assumes chlorine is added as NaOCl
    free_chlorine_added = mg_per_L_NaOCl_added / 74.44 / 1000 * 70900  # mg/L as NaOCl
    total_chlorine_inlet = free_chlorine_added / 70900  # mol/L
    total_molar_density += total_chlorine_inlet

    model.fs.ideal_naocl_chlorination_unit.inlet.mole_frac_comp[0, "OCl_-"].set_value(
        total_chlorine_inlet / total_molar_density
    )
    model.fs.ideal_naocl_chlorination_unit.inlet.mole_frac_comp[0, "Na_+"].set_value(
        total_chlorine_inlet / total_molar_density
    )

    set_H2O_molefraction(model.fs.ideal_naocl_chlorination_unit.inlet)

    model.fs.ideal_naocl_chlorination_unit.inlet.pressure[0].set_value(
        inlet_pressure_Pa
    )
    model.fs.ideal_naocl_chlorination_unit.inlet.temperature[0].set_value(
        inlet_temperature_K
    )
    model.fs.ideal_naocl_chlorination_unit.inlet.flow_mol[0].set_value(
        inlet_flow_mol_per_s
    )

    model.fs.ideal_naocl_chlorination_unit.free_chlorine.set_value(mg_per_L_NaOCl_added)


def fix_ideal_naocl_mixer_inlets(model):
    model.fs.ideal_naocl_mixer_unit.inlet_stream.flow_mol[0].fix()
    model.fs.ideal_naocl_mixer_unit.inlet_stream.pressure[0].fix()
    model.fs.ideal_naocl_mixer_unit.inlet_stream.temperature[0].fix()
    fix_all_molefractions(model.fs.ideal_naocl_mixer_unit.inlet_stream)

    model.fs.ideal_naocl_mixer_unit.naocl_stream.flow_mol[0].fix()
    model.fs.ideal_naocl_mixer_unit.naocl_stream.pressure[0].fix()
    model.fs.ideal_naocl_mixer_unit.naocl_stream.temperature[0].fix()
    fix_all_molefractions(model.fs.ideal_naocl_mixer_unit.naocl_stream)


def unfix_ideal_naocl_mixer_inlet_stream(model):
    model.fs.ideal_naocl_mixer_unit.inlet_stream.flow_mol[0].unfix()
    model.fs.ideal_naocl_mixer_unit.inlet_stream.pressure[0].unfix()
    model.fs.ideal_naocl_mixer_unit.inlet_stream.temperature[0].unfix()
    unfix_all_molefractions(model.fs.ideal_naocl_mixer_unit.inlet_stream)


def unfix_ideal_naocl_mixer_naocl_stream(model):
    model.fs.ideal_naocl_mixer_unit.naocl_stream.flow_mol[0].unfix()
    model.fs.ideal_naocl_mixer_unit.naocl_stream.pressure[0].unfix()
    model.fs.ideal_naocl_mixer_unit.naocl_stream.temperature[0].unfix()
    unfix_all_molefractions(model.fs.ideal_naocl_mixer_unit.naocl_stream)


def fix_ideal_naocl_chlorination_inlets(model):
    model.fs.ideal_naocl_chlorination_unit.inlet.flow_mol[0].fix()
    model.fs.ideal_naocl_chlorination_unit.inlet.pressure[0].fix()
    model.fs.ideal_naocl_chlorination_unit.inlet.temperature[0].fix()
    fix_all_molefractions(model.fs.ideal_naocl_chlorination_unit.inlet)


def unfix_ideal_naocl_chlorination_inlets(model):
    model.fs.ideal_naocl_chlorination_unit.inlet.flow_mol[0].unfix()
    model.fs.ideal_naocl_chlorination_unit.inlet.pressure[0].unfix()
    model.fs.ideal_naocl_chlorination_unit.inlet.temperature[0].unfix()
    unfix_all_molefractions(model.fs.ideal_naocl_chlorination_unit.inlet)


def scale_ideal_naocl_chlorination(unit, rxn_params, thermo_params, rxn_config):
    state_args, stoich_extents = approximate_chemical_state_args(
        unit, rxn_params, rxn_config
    )
    calculate_chemical_scaling_factors(unit, thermo_params, rxn_params, state_args)
    return state_args


def initialize_ideal_naocl_mixer(unit, debug_out=False):
    was_fixed = False
    if not unit.naocl_stream.flow_mol[0].is_fixed():
        unit.naocl_stream.flow_mol[0].fix()
        was_fixed = True
    unit.initialize(
        optarg=solver.options, outlvl=idaeslog.DEBUG if debug_out else idaeslog.NOTSET
    )
    if was_fixed:
        unit.naocl_stream.flow_mol[0].unfix()


def initialize_ideal_naocl_chlorination(unit, state_args, debug_out=False):
    unit.initialize(
        state_args=state_args,
        optarg=solver.options,
        outlvl=idaeslog.DEBUG if debug_out else idaeslog.NOTSET,
    )


def setup_block_to_solve_naocl_dosing_rate(model, free_chlorine_mg_per_L=2):
    model.fs.ideal_naocl_chlorination_unit.free_chlorine.fix(free_chlorine_mg_per_L)
    model.fs.ideal_naocl_mixer_unit.naocl_stream.flow_mol[0].unfix()


def display_results_of_ideal_naocl_mixer(unit):
    print()
    print("=========== Ideal NaOCl Mixer Results ============")
    print("Outlet Temperature:       \t" + str(unit.outlet.temperature[0].value))
    print("Outlet Pressure:          \t" + str(unit.outlet.pressure[0].value))
    print("Outlet FlowMole:          \t" + str(unit.outlet.flow_mol[0].value))
    print()
    total_molar_density = 55.6
    total_salt = value(unit.outlet.mole_frac_comp[0, "Na_+"]) * total_molar_density * 23
    total_salt += (
        value(unit.outlet.mole_frac_comp[0, "Cl_-"]) * total_molar_density * 35.4
    )
    psu = total_salt / (total_molar_density * 18) * 1000
    print("STP Salinity (PSU):           \t" + str(psu))
    print("NaOCl Dosing Rate (mg/s): \t" + str(unit.dosing_rate.value))
    print("-------------------------------------------------")
    print()


def display_results_of_chlorination_unit(unit):
    print()
    print("=========== Chlorination Results ============")
    print("Outlet Temperature:       \t" + str(unit.outlet.temperature[0].value))
    print("Outlet Pressure:          \t" + str(unit.outlet.pressure[0].value))
    print("Outlet FlowMole:          \t" + str(unit.outlet.flow_mol[0].value))
    print()
    total_molar_density = (
        value(unit.control_volume.properties_out[0.0].dens_mol_phase["Liq"]) / 1000
    )
    pH = -value(log10(unit.outlet.mole_frac_comp[0, "H_+"] * total_molar_density))
    print("pH at Outlet:             \t" + str(pH))
    total_salt = value(unit.outlet.mole_frac_comp[0, "Na_+"]) * total_molar_density * 23
    total_salt += (
        value(unit.outlet.mole_frac_comp[0, "Cl_-"]) * total_molar_density * 35.4
    )
    psu = total_salt / (total_molar_density * 18) * 1000
    print("Salinity (PSU):           \t" + str(psu))
    print("Free Chlorine (mg/L):     \t" + str(unit.free_chlorine.value))
    print("\tDistribution:")
    hocl = (
        value(
            unit.control_volume.properties_out[0.0].conc_mol_phase_comp["Liq", "HOCl"]
        )
        / 1000
    ) / (unit.free_chlorine.value / 70900)
    print("\t % HOCl: \t" + str(hocl * 100))
    ocl = (
        value(
            unit.control_volume.properties_out[0.0].conc_mol_phase_comp["Liq", "OCl_-"]
        )
        / 1000
    ) / (unit.free_chlorine.value / 70900)
    print("\t % OCl-: \t" + str(ocl * 100))
    print("-------------------------------------------")
    print()


def build_ideal_naocl_chlorination_block(model, expand_arcs=False):

    # Add properties to model
    build_ideal_naocl_prop(model)

    # Add mixer to the model
    build_ideal_naocl_mixer_unit(model)

    # Add reactor to the model
    build_ideal_naocl_chlorination_unit(model)

    # Connect the mixer to the chlorination unit with arcs
    model.fs.ideal_nacol_arc_mixer_to_chlor = Arc(
        source=model.fs.ideal_naocl_mixer_unit.outlet,
        destination=model.fs.ideal_naocl_chlorination_unit.inlet,
    )
    if expand_arcs == True:
        TransformationFactory("network.expand_arcs").apply_to(model)


# This method assumes that the flowsheet has a properties object named prop_TDS
def build_translator_from_RO_to_chlorination_block(model):
    # Translator inlet from RO and outlet goes to chlorination
    model.fs.RO_to_Chlor = Translator(
        inlet_property_package=model.fs.prop_TDS,
        outlet_property_package=model.fs.ideal_naocl_thermo_params,
    )

    # Add constraints to define how the translator will function
    model.fs.RO_to_Chlor.eq_equal_temperature = Constraint(
        expr=model.fs.RO_to_Chlor.inlet.temperature[0]
        == model.fs.RO_to_Chlor.outlet.temperature[0]
    )
    model.fs.RO_to_Chlor.eq_equal_pressure = Constraint(
        expr=model.fs.RO_to_Chlor.inlet.pressure[0]
        == model.fs.RO_to_Chlor.outlet.pressure[0]
    )

    model.fs.RO_to_Chlor.total_flow_cons = Constraint(
        expr=model.fs.RO_to_Chlor.outlet.flow_mol[0]
        == (model.fs.RO_to_Chlor.inlet.flow_mass_phase_comp[0, "Liq", "H2O"] / 18e-3)
        + (model.fs.RO_to_Chlor.inlet.flow_mass_phase_comp[0, "Liq", "TDS"] / 58.4e-3)
    )

    model.fs.RO_to_Chlor.H_con = Constraint(
        expr=model.fs.RO_to_Chlor.outlet.mole_frac_comp[0, "H_+"] == 0
    )
    model.fs.RO_to_Chlor.OH_con = Constraint(
        expr=model.fs.RO_to_Chlor.outlet.mole_frac_comp[0, "OH_-"] == 0
    )
    model.fs.RO_to_Chlor.HOCl_con = Constraint(
        expr=model.fs.RO_to_Chlor.outlet.mole_frac_comp[0, "HOCl"] == 0
    )
    model.fs.RO_to_Chlor.OCl_con = Constraint(
        expr=model.fs.RO_to_Chlor.outlet.mole_frac_comp[0, "OCl_-"] == 0
    )

    model.fs.RO_to_Chlor.Cl_con = Constraint(
        expr=model.fs.RO_to_Chlor.outlet.mole_frac_comp[0, "Cl_-"]
        == (model.fs.RO_to_Chlor.inlet.flow_mass_phase_comp[0, "Liq", "TDS"] / 58.4e-3)
        / model.fs.RO_to_Chlor.outlet.flow_mol[0]
    )

    model.fs.RO_to_Chlor.Na_con = Constraint(
        expr=model.fs.RO_to_Chlor.outlet.mole_frac_comp[0, "Na_+"]
        == (model.fs.RO_to_Chlor.inlet.flow_mass_phase_comp[0, "Liq", "TDS"] / 58.4e-3)
        / model.fs.RO_to_Chlor.outlet.flow_mol[0]
        + model.fs.RO_to_Chlor.outlet.mole_frac_comp[0, "OCl_-"]
    )

    model.fs.RO_to_Chlor.H2O_con = Constraint(
        expr=model.fs.RO_to_Chlor.outlet.mole_frac_comp[0, "H2O"]
        == 1
        - sum(
            model.fs.RO_to_Chlor.outlet.mole_frac_comp[0, j]
            for j in ["H_+", "OH_-", "HOCl", "OCl_-", "Cl_-", "Na_+"]
        )
    )

    iscale.calculate_scaling_factors(model.fs.RO_to_Chlor)


def run_ideal_naocl_mixer_example(fixed_dosage=False):
    model = ConcreteModel()
    model.fs = FlowsheetBlock(dynamic=False)

    # Add properties to model
    build_ideal_naocl_prop(model)

    # Add mixer to the model
    build_ideal_naocl_mixer_unit(model)

    # Set some inlet values
    set_ideal_naocl_mixer_inlets(model)

    # Fix the inlets for a solve
    fix_ideal_naocl_mixer_inlets(model)

    # unfix the flow_mol for the naocl_stream and fix dosing rate (alt form)
    if fixed_dosage == True:
        model.fs.ideal_naocl_mixer_unit.naocl_stream.flow_mol[0].unfix()
        model.fs.ideal_naocl_mixer_unit.dosing_rate.fix()

    check_dof(model)

    # initializer mixer
    initialize_ideal_naocl_mixer(model.fs.ideal_naocl_mixer_unit)

    solve_block(model, tee=True)

    display_results_of_ideal_naocl_mixer(model.fs.ideal_naocl_mixer_unit)

    return model


def run_ideal_naocl_chlorination_example():
    model = ConcreteModel()
    model.fs = FlowsheetBlock(dynamic=False)

    # add properties to model
    build_ideal_naocl_prop(model)

    # add equilibrium unit
    build_ideal_naocl_chlorination_unit(model)

    # Set some inlet values
    set_ideal_naocl_chlorination_inlets(model)

    # fix the inlets
    fix_ideal_naocl_chlorination_inlets(model)

    check_dof(model)

    # scale the chlorination unit
    state_args = scale_ideal_naocl_chlorination(
        model.fs.ideal_naocl_chlorination_unit,
        model.fs.ideal_naocl_rxn_params,
        model.fs.ideal_naocl_thermo_params,
        ideal_naocl_reaction_config,
    )

    # initialize the unit
    initialize_ideal_naocl_chlorination(
        model.fs.ideal_naocl_chlorination_unit, state_args, debug_out=False
    )

    results = solver.solve(model, tee=True)
    assert_optimal_termination(results)

    display_results_of_chlorination_unit(model.fs.ideal_naocl_chlorination_unit)

    return model


def run_chlorination_block_example(fix_free_chlorine=False):
    model = ConcreteModel()
    model.fs = FlowsheetBlock(dynamic=False)

    # Build the partial flowsheet of a mixer and chlorination unit
    build_ideal_naocl_chlorination_block(model, expand_arcs=True)

    # test the naocl_chlorination_costing
    model.fs.treated_flow_vol = Expression(expr=0.85 * pyunits.m**3 / pyunits.s)

    # set some values (using defaults for testing)
    set_ideal_naocl_mixer_inlets(
        model,
        dosing_rate_of_NaOCl_mg_per_s=0.4,
        inlet_water_density_kg_per_L=1,
        inlet_temperature_K=298,
        inlet_pressure_Pa=101325,
        inlet_flow_mol_per_s=25,
    )
    set_ideal_naocl_chlorination_inlets(model)

    # fix only the inlets for the mixer
    fix_ideal_naocl_mixer_inlets(model)
    initialize_ideal_naocl_mixer(model.fs.ideal_naocl_mixer_unit)

    # scale and initialize the chlorination unit
    state_args = scale_ideal_naocl_chlorination(
        model.fs.ideal_naocl_chlorination_unit,
        model.fs.ideal_naocl_rxn_params,
        model.fs.ideal_naocl_thermo_params,
        ideal_naocl_reaction_config,
    )

    # initialize the unit
    initialize_ideal_naocl_chlorination(
        model.fs.ideal_naocl_chlorination_unit, state_args, debug_out=False
    )

    check_dof(model)

    # Scale the full model and call the seq_decomp_initializer
    seq_decomp_initializer(model)

    if fix_free_chlorine:
        setup_block_to_solve_naocl_dosing_rate(model)

    results = solver.solve(model)
    costing.build_costing(model)
    model.fs.costing.initialize()

    results = solver.solve(model, tee=True)
    assert_optimal_termination(results)

    display_results_of_ideal_naocl_mixer(model.fs.ideal_naocl_mixer_unit)
    display_results_of_chlorination_unit(model.fs.ideal_naocl_chlorination_unit)

    return model
