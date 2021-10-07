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
    This test is of the core IDAES components that allow for the declaration
    and usage of solids phases in conjunction with aqueous phases.

    In these tests, we assess the convergence behavior of the solubility problem
    under various pH conditions for a solubility reaction that is inherently pH
    dependent. Several cases will be tested under different pH conditions and
    with different levels of complexity.

    Case 1 (remineralization, post-RO):
        Aqueous Rxns:   H2O <--> H + OH             logK = -14
        Solubility:     Ca(OH)2 <--> Ca + 2 OH      logK = -5.26


    Case 2 (softening, without explicit lime):
        Aqueous Rxns:   H2O <--> H + OH             logK = -14
                        H2CO3 <--> H + HCO3         logK = -6.3
                        HCO3 <--> H + CO3           logK = -10.2
        Solubility:     CaCO3 <--> Ca + CO3         logK = -12


    Case 3: (softening, with lime)
        Aqueous Rxns:   H2O <--> H + OH             logK = -14
                        H2CO3 <--> H + HCO3         logK = -6.3
                        HCO3 <--> H + CO3           logK = -10.2
        Solubility:     CaCO3 <--> Ca + CO3         logK = -12
                        Ca(OH)2 <--> Ca + 2 OH      logK = -5.26


    Case 4: (phosphorus removal, most realistic test)
        Aqueous Rxns:
                        H2O <---> H + OH
                        H2CO3 <---> H + HCO3
                        HCO3 <---> H + CO3
                        H2PO4 <---> H + HPO4
                        HPO4 <---> H + PO4
                        FeOH <---> Fe + OH
                        Fe(OH)2 <---> FeOH + OH
                        Fe(OH)3 <---> Fe(OH)2 + OH
                        Fe(OH)4 <---> Fe(OH)3 + OH
        Solubility:
                        FePO4 <---> Fe + PO4
"""

# Importing testing libraries
import pytest

# Importing the object for units from pyomo
from pyomo.environ import units as pyunits
from pyomo.environ import Var

# Imports from idaes core
from idaes.core import AqueousPhase, SolidPhase, FlowsheetBlock, EnergyBalanceType
from idaes.core.components import Solvent, Solute, Cation, Anion, Component
from idaes.core.phases import PhaseType as PT

# Imports from idaes generic models
import idaes.generic_models.properties.core.pure.Perrys as Perrys
from idaes.generic_models.properties.core.pure.ConstantProperties import Constant
from idaes.generic_models.properties.core.state_definitions import FTPx, FpcTP
from idaes.generic_models.properties.core.eos.ideal import Ideal

# Importing the enum for concentration unit basis used in the 'get_concentration_term' function
from idaes.generic_models.properties.core.generic.generic_reaction import ConcentrationForm

# Import the object/function for heat of reaction
from idaes.generic_models.properties.core.reactions.dh_rxn import constant_dh_rxn

# Import safe log power law equation
from idaes.generic_models.properties.core.reactions.equilibrium_forms import log_power_law_equil, power_law_equil

# Import built-in van't Hoff function
from idaes.generic_models.properties.core.reactions.equilibrium_constant import van_t_hoff

from idaes.generic_models.properties.core.reactions.equilibrium_forms import \
    solubility_product, log_solubility_product, log_power_law_equil
from idaes.generic_models.properties.core.reactions.equilibrium_constant import \
    ConstantKeq

# Import specific pyomo objects
from pyomo.environ import (ConcreteModel,
                           SolverStatus,
                           TerminationCondition,
                           value,
                           Suffix)

from idaes.core.util import scaling as iscale
from idaes.core.util.initialization import fix_state_vars, revert_state_vars

import idaes.logger as idaeslog

# Import pyomo methods to check the system units
from pyomo.util.check_units import assert_units_consistent

# Import idaes methods to check the model during construction
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              fixed_variables_set,
                                              activated_constraints_set,
                                              number_variables,
                                              number_total_constraints,
                                              number_unused_variables)

# Import the idaes objects for Generic Properties and Reactions
from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterBlock)
from idaes.generic_models.properties.core.generic.generic_reaction import (
        GenericReactionParameterBlock)

# Import the idaes object for the EquilibriumReactor unit model
from idaes.generic_models.unit_models.equilibrium_reactor import EquilibriumReactor

# Import log10 function from pyomo
from pyomo.environ import log10

__author__ = "Austin Ladshaw"

# Case 1 Config
case1_thermo_config = {
    "components": {
        'H2O': {"type": Solvent, "valid_phase_types": PT.aqueousPhase,
              "dens_mol_liq_comp": Constant,
              "enth_mol_liq_comp": Constant,
              "cp_mol_liq_comp": Constant,
              "entr_mol_liq_comp": Constant,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (18.0153, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol*pyunits.m**-3),
                    "cp_mol_liq_comp_coeff": (75.312, pyunits.J/pyunits.mol/pyunits.K),
                    "enth_mol_form_liq_comp_ref": (-285, pyunits.kJ/pyunits.mol),
                    "entr_mol_form_liq_comp_ref": (69, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'H_+': {"type": Cation, "charge": 1,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Constant,
              "enth_mol_liq_comp": Constant,
              "cp_mol_liq_comp": Constant,
              "entr_mol_liq_comp": Constant,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (1.00784, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol*pyunits.m**-3),
                    "cp_mol_liq_comp_coeff": (75.312, pyunits.J/pyunits.mol/pyunits.K),
                    "enth_mol_form_liq_comp_ref": (0, pyunits.kJ/pyunits.mol),
                    "entr_mol_form_liq_comp_ref": (0, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'OH_-': {"type": Anion, "charge": -1,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Constant,
              "enth_mol_liq_comp": Constant,
              "cp_mol_liq_comp": Constant,
              "entr_mol_liq_comp": Constant,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (17.008, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol*pyunits.m**-3),
                    "cp_mol_liq_comp_coeff": (75.312, pyunits.J/pyunits.mol/pyunits.K),
                    "enth_mol_form_liq_comp_ref": (-230, pyunits.kJ/pyunits.mol),
                    "entr_mol_form_liq_comp_ref": (-10, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'Ca_2+': {"type": Cation, "charge": 2,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Constant,
              "enth_mol_liq_comp": Constant,
              "cp_mol_liq_comp": Constant,
              "entr_mol_liq_comp": Constant,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (40.078, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": (55, pyunits.kmol*pyunits.m**-3),
                    "cp_mol_liq_comp_coeff": (167039, pyunits.J/pyunits.kmol/pyunits.K),
                    "enth_mol_form_liq_comp_ref": (-542.83, pyunits.J/pyunits.mol),
                    "entr_mol_form_liq_comp_ref": (-53, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'Ca(OH)2': {"type": Component, "valid_phase_types": PT.solidPhase,
              "dens_mol_sol_comp": Constant,
              "enth_mol_sol_comp": Constant,
              "cp_mol_sol_comp": Constant,
              "entr_mol_sol_comp": Constant,
              "parameter_data": {
                    "mw": (74.093, pyunits.g/pyunits.mol),
                    "dens_mol_sol_comp_coeff": (55, pyunits.kmol*pyunits.m**-3),
                    "cp_mol_sol_comp_coeff": (167039, pyunits.J/pyunits.kmol/pyunits.K),
                    "enth_mol_form_sol_comp_ref": (-986, pyunits.kJ/pyunits.mol),
                    "entr_mol_form_sol_comp_ref": (83, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    },

              },
              # End Component list
        "phases":  {'Liq': {"type": AqueousPhase,
                            "equation_of_state": Ideal},
                    'Sol': {"type": SolidPhase,
                            "equation_of_state": Ideal}
                    },

        "state_definition": FpcTP,
        "state_bounds": {
                         "temperature": (273.15, 300, 650),
                         "pressure": (5e4, 1e5, 1e6)
                     },

        "pressure_ref": 1e5,
        "temperature_ref": 300,
        "base_units": {"time": pyunits.s,
                       "length": pyunits.m,
                       "mass": pyunits.kg,
                       "amount": pyunits.mol,
                       "temperature": pyunits.K},

    }
    # End thermo_config definition

# Case 1 rxn config (without log_solubility_product)
case1_rxn_config = {
    "base_units": {"time": pyunits.s,
                   "length": pyunits.m,
                   "mass": pyunits.kg,
                   "amount": pyunits.mol,
                   "temperature": pyunits.K},
    "equilibrium_reactions": {
        "CaOH2_Ksp": {
                    "stoichiometry": {  ("Sol", "Ca(OH)2"): -1,
                                        ("Liq", "Ca_2+"): 1,
                                        ("Liq", "OH_-"): 2},
                    "heat_of_reaction": constant_dh_rxn,
                    "equilibrium_constant": van_t_hoff,
                    "equilibrium_form": solubility_product,
                    "concentration_form": ConcentrationForm.moleFraction,
                    "parameter_data": {
                        "dh_rxn_ref": (0.0, pyunits.J/pyunits.mol),
                        "k_eq_ref": (10**-5.26/55.2/55.2/55.2, pyunits.dimensionless),
                        "T_eq_ref": (298.0, pyunits.K),
                        "reaction_order": { ("Sol", "Ca(OH)2"): 0,
                                            ("Liq", "Ca_2+"): 1,
                                            ("Liq", "OH_-"): 2}
                        }
                        # End parameter_data
                },
        "H2O_Kw": {
                "stoichiometry": {("Liq", "H2O"): -1,
                                 ("Liq", "H_+"): 1,
                                 ("Liq", "OH_-"): 1},
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
                   "reaction_order": {("Liq", "H2O"): 0,
                                    ("Liq", "H_+"): 1,
                                    ("Liq", "OH_-"): 1}
                    }
                    # End parameter_data
               }
                # End Reactions
         }
         # End equilibrium_reactions
    }
    # End reaction_config definition

# Case 1 rxn config (with log_solubility_product)
case1_log_rxn_config = {
    "base_units": {"time": pyunits.s,
                   "length": pyunits.m,
                   "mass": pyunits.kg,
                   "amount": pyunits.mol,
                   "temperature": pyunits.K},
    "equilibrium_reactions": {
        "CaOH2_Ksp": {
                    "stoichiometry": {  ("Sol", "Ca(OH)2"): -1,
                                        ("Liq", "Ca_2+"): 1,
                                        ("Liq", "OH_-"): 2},
                    "heat_of_reaction": constant_dh_rxn,
                    "equilibrium_constant": van_t_hoff,
                    "equilibrium_form": log_solubility_product,
                    "concentration_form": ConcentrationForm.moleFraction,
                    "parameter_data": {
                        "dh_rxn_ref": (0.0, pyunits.J/pyunits.mol),
                        "k_eq_ref": (10**-5.26/55.2/55.2/55.2, pyunits.dimensionless),
                        "T_eq_ref": (298.0, pyunits.K),
                        "reaction_order": { ("Sol", "Ca(OH)2"): 0,
                                            ("Liq", "Ca_2+"): 1,
                                            ("Liq", "OH_-"): 2}
                        }
                        # End parameter_data
                },
        "H2O_Kw": {
                "stoichiometry": {("Liq", "H2O"): -1,
                                 ("Liq", "H_+"): 1,
                                 ("Liq", "OH_-"): 1},
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
                   "reaction_order": {("Liq", "H2O"): 0,
                                    ("Liq", "H_+"): 1,
                                    ("Liq", "OH_-"): 1}
                    }
                    # End parameter_data
               }
                # End Reactions
         }
         # End equilibrium_reactions
    }
    # End reaction_config definition

# Get default solver for testing
solver = get_solver()

## Set of helper functions for scaling
def _set_eps_vals(rxn_params, rxn_config):
    # For solubility reactions, have to set the eps value
    try:
        for rid in rxn_params.equilibrium_reaction_idx:
            scale = rxn_config["equilibrium_reactions"][rid]["parameter_data"]["k_eq_ref"][0]

            # NOTE: The solubility_product function STILL has an eps value that we need to set
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

def _set_equ_rxn_scaling(unit, rxn_config):
    #Add scaling factors for reactions (changes depending on if it is a log form or not)
    for i in unit.control_volume.equilibrium_reaction_extent_index:
        # i[0] = time, i[1] = reaction

        scale = rxn_config["equilibrium_reactions"][i[1]]["parameter_data"]["k_eq_ref"][0]

        # this may also need to be different for solubility_product
        if rxn_config["equilibrium_reactions"][i[1]]["equilibrium_form"] == solubility_product:
            iscale.set_scaling_factor(unit.control_volume.equilibrium_reaction_extent[0.0,i[1]], 1e-1/scale)
            iscale.constraint_scaling_transform(
                unit.control_volume.reactions[0.0].equilibrium_constraint[i[1]], 1e-1/scale)
        else:
            iscale.set_scaling_factor(unit.control_volume.equilibrium_reaction_extent[0.0,i[1]], 100/scale)
            iscale.constraint_scaling_transform(
                unit.control_volume.reactions[0.0].equilibrium_constraint[i[1]], 0.1)

def _set_mat_bal_scaling_FpcTP(unit):
    # For species
    min = 1e-6
    for i in unit.control_volume.properties_out[0.0].mole_frac_phase_comp:
        # i[0] = phase, i[1] = species
        if unit.inlet.flow_mol_phase_comp[0, i[0], i[1]].value > min:
            scale = unit.inlet.flow_mol_phase_comp[0, i[0], i[1]].value
        else:
            scale = min

        ## For now, these are both set the same, but may change in testing
        # Scaling factors for liquid
        if i[0] == 'Liq':
            iscale.set_scaling_factor(unit.control_volume.properties_out[0.0].mole_frac_comp[i[1]], 10/scale)
            iscale.set_scaling_factor(unit.control_volume.properties_out[0.0].mole_frac_phase_comp[i], 10/scale)
            iscale.set_scaling_factor(unit.control_volume.properties_out[0.0].flow_mol_phase_comp[i], 10/scale)
            iscale.constraint_scaling_transform(unit.control_volume.material_balances[0.0,i[1]], 10/scale)
        #Scaling factors for solids
        else:
            iscale.set_scaling_factor(unit.control_volume.properties_out[0.0].mole_frac_comp[i[1]], 10/scale)
            iscale.set_scaling_factor(unit.control_volume.properties_out[0.0].mole_frac_phase_comp[i], 10/scale)
            iscale.set_scaling_factor(unit.control_volume.properties_out[0.0].flow_mol_phase_comp[i], 10/scale)
            iscale.constraint_scaling_transform(unit.control_volume.material_balances[0.0,i[1]], 10/scale)

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

# Defaults to pH of 7 with no lime added
def run_case1(xOH=1e-7/55.2, xH=1e-7/55.2, xCaOH2=1e-20, xCa=1e-20,
            thermo_config=None, rxn_config=None, has_energy_balance=True):
    print("==========================================================================")
    print("Case 1: Remineralization via lime dissolution")
    print("xOH = "+str(xOH))
    print("xH = "+str(xH))
    print("xCaOH2 = "+str(xCaOH2))
    print("Initial pH = " +str(-log10(xH*55.2)))
    print()

    model = ConcreteModel()
    model.fs = FlowsheetBlock(default={"dynamic": False})
    model.fs.thermo_params = GenericParameterBlock(default=thermo_config)

    model.fs.rxn_params = GenericReactionParameterBlock(
            default={"property_package": model.fs.thermo_params,
                    **rxn_config
                    })

    args = {"property_package": model.fs.thermo_params,
            "reaction_package": model.fs.rxn_params,
            "has_rate_reactions": False,
            "has_equilibrium_reactions": True,
            "has_heat_transfer": False,
            "has_heat_of_reaction": False,
            "has_pressure_change": False}
    if has_energy_balance == False:
        args["energy_balance_type"] = EnergyBalanceType.none

    model.fs.unit = EquilibriumReactor(default=args)

    total_flow_mol = 10

    # Set flow_mol_phase_comp
    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Ca_2+"].fix( xCa*total_flow_mol )

    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H_+"].fix( xH*total_flow_mol )
    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "OH_-"].fix( xOH*total_flow_mol )
    model.fs.unit.inlet.flow_mol_phase_comp[0, "Sol", "Ca(OH)2"].fix( xCaOH2*total_flow_mol )
    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix( (1-xH-xOH-xCaOH2-xCa)*total_flow_mol )

    model.fs.unit.inlet.pressure.fix(101325.0)
    model.fs.unit.inlet.temperature.fix(298.)
    if has_energy_balance == False:
        model.fs.unit.outlet.temperature.fix(298.)

    assert (degrees_of_freedom(model) == 0)

    ## ==================== Start Scaling for this problem ===========================
    _set_eps_vals(model.fs.rxn_params, rxn_config)
    _set_equ_rxn_scaling(model.fs.unit, rxn_config)
    _set_mat_bal_scaling_FpcTP(model.fs.unit)
    if has_energy_balance == True:
        _set_ene_bal_scaling(model.fs.unit)

    iscale.calculate_scaling_factors(model.fs.unit)
    assert isinstance(model.fs.unit.control_volume.scaling_factor, Suffix)
    assert isinstance(model.fs.unit.control_volume.properties_out[0.0].scaling_factor, Suffix)
    assert isinstance(model.fs.unit.control_volume.properties_in[0.0].scaling_factor, Suffix)

    ## ==================== END Scaling for this problem ===========================

    # First pass, low bound_push, low mu_init
    #solver.options['bound_push'] = 1e-10
    #solver.options['mu_init'] = 1e-6
    #model.fs.unit.initialize(optarg=solver.options, outlvl=idaeslog.DEBUG)

    # Second pass, loosen bound_push to get more acceptable solve
    solver.options['bound_push'] = 1e-5
    solver.options['mu_init'] = 1e-3
    model.fs.unit.initialize(optarg=solver.options, outlvl=idaeslog.DEBUG)

    assert degrees_of_freedom(model) == 0

    # Solve full model
    solver.options['bound_push'] = 1e-5
    solver.options['mu_init'] = 1e-3
    results = solver.solve(model, tee=True)

    assert results.solver.termination_condition == TerminationCondition.optimal
    assert results.solver.status == SolverStatus.ok

    print("comp\toutlet.tot_molfrac")
    for i in model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp:
        print(str(i)+"\t"+str(value( model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp[i] )))
    print()

    # NOTE: Changed all to mole fraction
    for i in model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp:
        print(str(i)+"\t"+str(value(model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[i])))
    print()

    Ca = value(model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp["Liq","Ca_2+"])
    OH = value(model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp["Liq","OH_-"])
    H = value(model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp["Liq","H_+"])
    Ksp = value(model.fs.unit.control_volume.reactions[0.0].k_eq["CaOH2_Ksp"].expr)

    print("Final pH = " +str(-log10(H*55.2)))
    print("Expected max/min pH = " +str(14+log10(xCaOH2*55.2*2)))

    print()
    if Ksp*1.01 >= Ca*OH*OH:
        print("Constraint is satisfied!")
    else:
        print("Constraint is VIOLATED!")
        print("\tRelative error: "+str(Ksp/Ca/OH/OH)+">=1")
        assert False
    print("Ksp =\t"+str(Ksp))
    print("Ca*OH**2 =\t"+str(Ca*OH*OH))

    print("==========================================================================")

    return model

## ================================= Case 1 Tests ===============================
@pytest.mark.component
def test_case_1_no_dissolution():
    model = run_case1(xOH=1e-7/55.2, xH=1e-7/55.2, xCaOH2=1e-20,
                        thermo_config=case1_thermo_config, rxn_config=case1_log_rxn_config,
                        has_energy_balance=True)

@pytest.mark.component
def test_case_1_high_dissolution():
    model = run_case1(xOH=1e-7/55.2, xH=1e-7/55.2, xCaOH2=1e-5,
                        thermo_config=case1_thermo_config, rxn_config=case1_log_rxn_config,
                        has_energy_balance=True)

@pytest.mark.component
def test_case_1_mid_dissolution():
    model = run_case1(xOH=1e-7/55.2, xH=1e-7/55.2, xCaOH2=1e-7,
                        thermo_config=case1_thermo_config, rxn_config=case1_log_rxn_config,
                        has_energy_balance=True)

@pytest.mark.component
def test_case_1_low_dissolution():
    model = run_case1(xOH=1e-7/55.2, xH=1e-7/55.2, xCaOH2=1e-9,
                        thermo_config=case1_thermo_config, rxn_config=case1_log_rxn_config,
                        has_energy_balance=True)

@pytest.mark.component
def test_case_1_high_precipitation():
    model = run_case1(xOH=1e-1/55.2, xH=1e-13/55.2, xCaOH2=1e-20, xCa=1e-1,
                        thermo_config=case1_thermo_config, rxn_config=case1_log_rxn_config,
                        has_energy_balance=True)

@pytest.mark.component
def test_case_1_low_precipitation():
    model = run_case1(xOH=1e-3/55.2, xH=1e-11/55.2, xCaOH2=1e-20, xCa=1e-1,
                        thermo_config=case1_thermo_config, rxn_config=case1_log_rxn_config,
                        has_energy_balance=True)

# This is for additional testing
if __name__ == "__main__":
    ### All tests with NO dissolution passed ###
    '''
    model = run_case1(xOH=1e-7/55.2, xH=1e-7/55.2, xCaOH2=1e-20,
                        thermo_config=case1_thermo_config, rxn_config=case1_log_rxn_config,
                        has_energy_balance=True)
    '''

    ### This test of lime dissolution is working and giving accurate pH
    '''
    model = run_case1(xOH=1e-7/55.2, xH=1e-7/55.2, xCaOH2=1e-5,
                        thermo_config=case1_thermo_config, rxn_config=case1_log_rxn_config,
                        has_energy_balance=True)
    '''
    '''
    model = run_case1(xOH=1e-7/55.2, xH=1e-7/55.2, xCaOH2=1e-7,
                        thermo_config=case1_thermo_config, rxn_config=case1_log_rxn_config,
                        has_energy_balance=True)
    '''
    '''
    model = run_case1(xOH=1e-7/55.2, xH=1e-7/55.2, xCaOH2=1e-9,
                        thermo_config=case1_thermo_config, rxn_config=case1_log_rxn_config,
                        has_energy_balance=True)
    '''

    ### These set of tests demonstrate precipitation of lime
    '''
    model = run_case1(xOH=1e-1/55.2, xH=1e-13/55.2, xCaOH2=1e-20, xCa=1e-1,
                        thermo_config=case1_thermo_config, rxn_config=case1_log_rxn_config,
                        has_energy_balance=True)
    '''
    '''
    model = run_case1(xOH=1e-2/55.2, xH=1e-12/55.2, xCaOH2=1e-20, xCa=1e-1,
                        thermo_config=case1_thermo_config, rxn_config=case1_log_rxn_config,
                        has_energy_balance=True)
    '''
    '''
    model = run_case1(xOH=1e-3/55.2, xH=1e-11/55.2, xCaOH2=1e-20, xCa=1e-1,
                        thermo_config=case1_thermo_config, rxn_config=case1_log_rxn_config,
                        has_energy_balance=True)
    '''
