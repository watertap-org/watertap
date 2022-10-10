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


    [NOTE: IPOPT/IDAES does not handle simultaneous dissolution and
            precipitation very well. To implement this in a flowsheet,
            there should be separate units for dissolution and precipitation.
            For instance, have 1 unit for addition and dissolution of lime,
            then a separate downstream unit for precipitation of CaCO3.]
    Case 3: (softening, with lime: precipitation + dissolution)
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
from idaes.core.base.components import Solvent, Solute, Cation, Anion, Component
from idaes.core.base.phases import PhaseType as PT

# Imports from idaes generic models
import idaes
import idaes.models.properties.modular_properties.pure.Perrys as Perrys
from idaes.models.properties.modular_properties.pure.ConstantProperties import Constant
from idaes.models.properties.modular_properties.state_definitions import FTPx, FpcTP
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
    power_law_equil,
)

# Import built-in van't Hoff function
from idaes.models.properties.modular_properties.reactions.equilibrium_constant import (
    van_t_hoff,
)

from idaes.models.properties.modular_properties.reactions.equilibrium_forms import (
    solubility_product,
    log_solubility_product,
    log_power_law_equil,
)
from idaes.models.properties.modular_properties.reactions.equilibrium_constant import (
    ConstantKeq,
)

# Import specific pyomo objects
from pyomo.environ import (
    ConcreteModel,
    SolverStatus,
    TerminationCondition,
    value,
    Suffix,
)

from idaes.core.util import scaling as iscale
from idaes.core.util.initialization import fix_state_vars, revert_state_vars

import idaes.logger as idaeslog

# Import pyomo methods to check the system units
from pyomo.util.check_units import assert_units_consistent

# Import idaes methods to check the model during construction
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    fixed_variables_set,
    activated_constraints_set,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)

# Import the idaes objects for Generic Properties and Reactions
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models.properties.modular_properties.base.generic_reaction import (
    GenericReactionParameterBlock,
)

# Import the idaes object for the EquilibriumReactor unit model
from idaes.models.unit_models.equilibrium_reactor import EquilibriumReactor

# Import log10 function from pyomo
from pyomo.environ import log10

# Import scaling helper functions
from watertap.examples.chemistry.chem_scaling_utils import (
    _set_eps_vals,
    _set_equ_rxn_scaling,
    _set_mat_bal_scaling_FpcTP,
    _set_mat_bal_scaling_FTPx,
    _set_ene_bal_scaling,
)

__author__ = "Austin Ladshaw"

# Case 1 Config
case1_thermo_config = {
    "components": {
        "H2O": {
            "type": Solvent,
            "valid_phase_types": PT.aqueousPhase,
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (18.0153, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                "cp_mol_liq_comp_coeff": (75.312, pyunits.J / pyunits.mol / pyunits.K),
                "enth_mol_form_liq_comp_ref": (-285, pyunits.kJ / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (69, pyunits.J / pyunits.K / pyunits.mol),
            },
            # End parameter_data
        },
        "H_+": {
            "type": Cation,
            "charge": 1,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (1.00784, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                "cp_mol_liq_comp_coeff": (75.312, pyunits.J / pyunits.mol / pyunits.K),
                "enth_mol_form_liq_comp_ref": (0, pyunits.kJ / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (0, pyunits.J / pyunits.K / pyunits.mol),
            },
            # End parameter_data
        },
        "OH_-": {
            "type": Anion,
            "charge": -1,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (17.008, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                "cp_mol_liq_comp_coeff": (75.312, pyunits.J / pyunits.mol / pyunits.K),
                "enth_mol_form_liq_comp_ref": (-230, pyunits.kJ / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (
                    -10,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
            # End parameter_data
        },
        "Ca_2+": {
            "type": Cation,
            "charge": 2,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (40.078, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55, pyunits.kmol * pyunits.m**-3),
                "cp_mol_liq_comp_coeff": (167039, pyunits.J / pyunits.kmol / pyunits.K),
                "enth_mol_form_liq_comp_ref": (-542.83, pyunits.J / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (
                    -53,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
            # End parameter_data
        },
        "Ca(OH)2": {
            "type": Component,
            "valid_phase_types": PT.solidPhase,
            "dens_mol_sol_comp": Constant,
            "enth_mol_sol_comp": Constant,
            "cp_mol_sol_comp": Constant,
            "entr_mol_sol_comp": Constant,
            "parameter_data": {
                "mw": (74.093, pyunits.g / pyunits.mol),
                "dens_mol_sol_comp_coeff": (55, pyunits.kmol * pyunits.m**-3),
                "cp_mol_sol_comp_coeff": (167039, pyunits.J / pyunits.kmol / pyunits.K),
                "enth_mol_form_sol_comp_ref": (-986, pyunits.kJ / pyunits.mol),
                "entr_mol_form_sol_comp_ref": (83, pyunits.J / pyunits.K / pyunits.mol),
            },
        },
    },
    # End Component list
    "phases": {
        "Liq": {"type": AqueousPhase, "equation_of_state": Ideal},
        "Sol": {"type": SolidPhase, "equation_of_state": Ideal},
    },
    "state_definition": FpcTP,
    "state_bounds": {"temperature": (273.15, 300, 650), "pressure": (5e4, 1e5, 1e6)},
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
# End thermo_config definition

# Case 1 rxn config (with log_solubility_product)
case1_log_rxn_config = {
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },
    "equilibrium_reactions": {
        "CaOH2_Ksp": {
            "stoichiometry": {
                ("Sol", "Ca(OH)2"): -1,
                ("Liq", "Ca_2+"): 1,
                ("Liq", "OH_-"): 2,
            },
            "heat_of_reaction": constant_dh_rxn,
            "equilibrium_constant": van_t_hoff,
            "equilibrium_form": log_solubility_product,
            "concentration_form": ConcentrationForm.moleFraction,
            "parameter_data": {
                "dh_rxn_ref": (0.0, pyunits.J / pyunits.mol),
                "k_eq_ref": (10**-5.26 / 55.2 / 55.2 / 55.2, pyunits.dimensionless),
                "T_eq_ref": (298.0, pyunits.K),
                "reaction_order": {
                    ("Sol", "Ca(OH)2"): 0,
                    ("Liq", "Ca_2+"): 1,
                    ("Liq", "OH_-"): 2,
                },
            }
            # End parameter_data
        },
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
        }
        # End Reactions
    }
    # End equilibrium_reactions
}
# End reaction_config definition

# Get default solver for testing
solver = get_solver()

# Defaults to pH of 7 with no lime added
def run_case1(
    xOH=1e-7 / 55.2,
    xH=1e-7 / 55.2,
    xCaOH2=1e-20,
    xCa=1e-20,
    thermo_config=None,
    rxn_config=None,
    has_energy_balance=True,
):
    print("==========================================================================")
    print("Case 1: Remineralization via lime dissolution")
    print("xOH = " + str(xOH))
    print("xH = " + str(xH))
    print("xCaOH2 = " + str(xCaOH2))
    print("Initial pH = " + str(-log10(xH * 55.2)))
    print()

    model = ConcreteModel()
    model.fs = FlowsheetBlock(dynamic=False)
    model.fs.thermo_params = GenericParameterBlock(**thermo_config)

    model.fs.rxn_params = GenericReactionParameterBlock(
        property_package=model.fs.thermo_params, **rxn_config
    )

    args = {
        "property_package": model.fs.thermo_params,
        "reaction_package": model.fs.rxn_params,
        "has_rate_reactions": False,
        "has_equilibrium_reactions": True,
        "has_heat_transfer": False,
        "has_heat_of_reaction": False,
        "has_pressure_change": False,
    }
    if has_energy_balance == False:
        args["energy_balance_type"] = EnergyBalanceType.none

    model.fs.unit = EquilibriumReactor(**args)

    total_flow_mol = 10

    # Set flow_mol_phase_comp
    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Ca_2+"].fix(xCa * total_flow_mol)

    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H_+"].fix(xH * total_flow_mol)
    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "OH_-"].fix(xOH * total_flow_mol)
    model.fs.unit.inlet.flow_mol_phase_comp[0, "Sol", "Ca(OH)2"].fix(
        xCaOH2 * total_flow_mol
    )
    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix(
        (1 - xH - xOH - xCaOH2 - xCa) * total_flow_mol
    )

    model.fs.unit.inlet.pressure.fix(101325.0)
    model.fs.unit.inlet.temperature.fix(298.0)
    if has_energy_balance == False:
        model.fs.unit.outlet.temperature.fix(298.0)

    assert degrees_of_freedom(model) == 0

    ## ==================== Start Scaling for this problem ===========================
    _set_eps_vals(model.fs.rxn_params, rxn_config)
    _set_equ_rxn_scaling(model.fs.unit, rxn_config)
    _set_mat_bal_scaling_FpcTP(model.fs.unit)
    if has_energy_balance == True:
        _set_ene_bal_scaling(model.fs.unit)

    iscale.calculate_scaling_factors(model.fs.unit)
    assert isinstance(model.fs.unit.control_volume.scaling_factor, Suffix)
    assert isinstance(
        model.fs.unit.control_volume.properties_out[0.0].scaling_factor, Suffix
    )
    assert isinstance(
        model.fs.unit.control_volume.properties_in[0.0].scaling_factor, Suffix
    )

    ## ==================== END Scaling for this problem ===========================

    model.fs.unit.initialize(optarg=solver.options, outlvl=idaeslog.DEBUG)

    assert degrees_of_freedom(model) == 0

    results = solver.solve(model, tee=True)

    assert results.solver.termination_condition == TerminationCondition.optimal
    assert results.solver.status == SolverStatus.ok

    print("comp\toutlet.tot_molfrac")
    for i in model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp:
        print(
            str(i)
            + "\t"
            + str(
                value(
                    model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp[i]
                )
            )
        )
    print()

    # NOTE: Changed all to mole fraction
    for i in model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp:
        print(
            str(i)
            + "\t"
            + str(
                value(
                    model.fs.unit.control_volume.properties_out[
                        0.0
                    ].mole_frac_phase_comp[i]
                )
            )
        )
    print()

    Ca = value(
        model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
            "Liq", "Ca_2+"
        ]
    )
    OH = value(
        model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
            "Liq", "OH_-"
        ]
    )
    H = value(
        model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
            "Liq", "H_+"
        ]
    )
    Ksp = value(model.fs.unit.control_volume.reactions[0.0].k_eq["CaOH2_Ksp"].expr)

    print("Final pH = " + str(-log10(H * 55.2)))
    print("Expected max/min pH = " + str(14 + log10(xCaOH2 * 55.2 * 2)))

    print()
    if Ksp * 1.01 >= Ca * OH * OH:
        print("Constraint is satisfied!")
    else:
        print("Constraint is VIOLATED!")
        print("\tRelative error: " + str(Ksp / Ca / OH / OH) + ">=1")
        assert False
    print("Ksp =\t" + str(Ksp))
    print("Ca*OH**2 =\t" + str(Ca * OH * OH))

    print("==========================================================================")

    return model


## ================================= Case 1 Tests ===============================
@pytest.mark.component
def test_case_1_no_dissolution():
    model = run_case1(
        xOH=1e-7 / 55.2,
        xH=1e-7 / 55.2,
        xCaOH2=1e-20,
        thermo_config=case1_thermo_config,
        rxn_config=case1_log_rxn_config,
        has_energy_balance=True,
    )


@pytest.mark.component
def test_case_1_high_dissolution():
    model = run_case1(
        xOH=1e-7 / 55.2,
        xH=1e-7 / 55.2,
        xCaOH2=1e-5,
        thermo_config=case1_thermo_config,
        rxn_config=case1_log_rxn_config,
        has_energy_balance=True,
    )


@pytest.mark.component
def test_case_1_mid_dissolution():
    model = run_case1(
        xOH=1e-7 / 55.2,
        xH=1e-7 / 55.2,
        xCaOH2=1e-7,
        thermo_config=case1_thermo_config,
        rxn_config=case1_log_rxn_config,
        has_energy_balance=True,
    )


@pytest.mark.component
def test_case_1_low_dissolution():
    model = run_case1(
        xOH=1e-7 / 55.2,
        xH=1e-7 / 55.2,
        xCaOH2=1e-9,
        thermo_config=case1_thermo_config,
        rxn_config=case1_log_rxn_config,
        has_energy_balance=True,
    )


@pytest.mark.component
def test_case_1_high_precipitation():
    model = run_case1(
        xOH=1e-1 / 55.2,
        xH=1e-13 / 55.2,
        xCaOH2=1e-20,
        xCa=1e-1,
        thermo_config=case1_thermo_config,
        rxn_config=case1_log_rxn_config,
        has_energy_balance=True,
    )


@pytest.mark.component
def test_case_1_low_precipitation():
    model = run_case1(
        xOH=1e-3 / 55.2,
        xH=1e-11 / 55.2,
        xCaOH2=1e-20,
        xCa=1e-1,
        thermo_config=case1_thermo_config,
        rxn_config=case1_log_rxn_config,
        has_energy_balance=True,
    )


# Case 2 Config
"""
    Case 2 (softening, without explicit lime):
        Aqueous Rxns:   H2O <--> H + OH             logK = -14
                        H2CO3 <--> H + HCO3         logK = -6.3
                        HCO3 <--> H + CO3           logK = -10.2
        Solubility:     CaCO3 <--> Ca + CO3         logK = -12
"""
case2_thermo_config = {
    "components": {
        "H2O": {
            "type": Solvent,
            "valid_phase_types": PT.aqueousPhase,
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (18.0153, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                "cp_mol_liq_comp_coeff": (75.312, pyunits.J / pyunits.mol / pyunits.K),
                "enth_mol_form_liq_comp_ref": (-285, pyunits.kJ / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (69, pyunits.J / pyunits.K / pyunits.mol),
            },
            # End parameter_data
        },
        "H_+": {
            "type": Cation,
            "charge": 1,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (1.00784, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                "cp_mol_liq_comp_coeff": (75.312, pyunits.J / pyunits.mol / pyunits.K),
                "enth_mol_form_liq_comp_ref": (0, pyunits.kJ / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (0, pyunits.J / pyunits.K / pyunits.mol),
            },
            # End parameter_data
        },
        "OH_-": {
            "type": Anion,
            "charge": -1,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (17.008, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                "cp_mol_liq_comp_coeff": (75.312, pyunits.J / pyunits.mol / pyunits.K),
                "enth_mol_form_liq_comp_ref": (-230, pyunits.kJ / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (
                    -10,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
            # End parameter_data
        },
        "Ca_2+": {
            "type": Cation,
            "charge": 2,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (40.078, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55, pyunits.kmol * pyunits.m**-3),
                "cp_mol_liq_comp_coeff": (167039, pyunits.J / pyunits.kmol / pyunits.K),
                "enth_mol_form_liq_comp_ref": (-542.83, pyunits.J / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (
                    -53,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
            # End parameter_data
        },
        "H2CO3": {
            "type": Solute,
            "valid_phase_types": PT.aqueousPhase,
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (62.03, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                "cp_mol_liq_comp_coeff": (75.312, pyunits.J / pyunits.mol / pyunits.K),
                "enth_mol_form_liq_comp_ref": (-699.7, pyunits.kJ / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (
                    187,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
            # End parameter_data
        },
        "HCO3_-": {
            "type": Anion,
            "charge": -1,
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (62.03, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                "cp_mol_liq_comp_coeff": (75.312, pyunits.J / pyunits.mol / pyunits.K),
                "enth_mol_form_liq_comp_ref": (-692, pyunits.kJ / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (
                    91.2,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
            # End parameter_data
        },
        "CO3_2-": {
            "type": Anion,
            "charge": -2,
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (62.03, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                "cp_mol_liq_comp_coeff": (75.312, pyunits.J / pyunits.mol / pyunits.K),
                "enth_mol_form_liq_comp_ref": (-677.1, pyunits.kJ / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (
                    -56.9,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
            # End parameter_data
        },
        "CaCO3": {
            "type": Component,
            "valid_phase_types": PT.solidPhase,
            "dens_mol_sol_comp": Constant,
            "enth_mol_sol_comp": Constant,
            "cp_mol_sol_comp": Constant,
            "entr_mol_sol_comp": Constant,
            "parameter_data": {
                "mw": (100.0869, pyunits.g / pyunits.mol),
                "dens_mol_sol_comp_coeff": (55, pyunits.kmol * pyunits.m**-3),
                "cp_mol_sol_comp_coeff": (167039, pyunits.J / pyunits.kmol / pyunits.K),
                "enth_mol_form_sol_comp_ref": (-1207, pyunits.kJ / pyunits.mol),
                "entr_mol_form_sol_comp_ref": (
                    91.7,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
        },
    },
    # End Component list
    "phases": {
        "Liq": {"type": AqueousPhase, "equation_of_state": Ideal},
        "Sol": {"type": SolidPhase, "equation_of_state": Ideal},
    },
    "state_definition": FpcTP,
    "state_bounds": {"temperature": (273.15, 300, 650), "pressure": (5e4, 1e5, 1e6)},
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
# End thermo_config definition

# Case 2 rxn config (with log_solubility_product)
case2_log_rxn_config = {
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },
    "equilibrium_reactions": {
        "CaCO3_Ksp": {
            "stoichiometry": {
                ("Sol", "CaCO3"): -1,
                ("Liq", "Ca_2+"): 1,
                ("Liq", "CO3_2-"): 1,
            },
            "heat_of_reaction": constant_dh_rxn,
            "equilibrium_constant": van_t_hoff,
            "equilibrium_form": log_solubility_product,
            "concentration_form": ConcentrationForm.moleFraction,
            "parameter_data": {
                "dh_rxn_ref": (0.0, pyunits.J / pyunits.mol),
                "k_eq_ref": (10**-12 / 55.2 / 55.2, pyunits.dimensionless),
                "T_eq_ref": (298.0, pyunits.K),
                "reaction_order": {
                    ("Sol", "CaCO3"): 0,
                    ("Liq", "Ca_2+"): 1,
                    ("Liq", "CO3_2-"): 1,
                },
            }
            # End parameter_data
        },
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
                "dh_rxn_ref": (7.7, pyunits.kJ / pyunits.mol),
                "k_eq_ref": (10**-6.35 / 55.2, pyunits.dimensionless),
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
                "dh_rxn_ref": (14.9, pyunits.kJ / pyunits.mol),
                "k_eq_ref": (10**-10.33 / 55.2, pyunits.dimensionless),
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
        # End last reaction
    }
    # End equilibrium_reactions
}
# End reaction_config definition

# Defaults to pH of 7 with no lime added
def run_case2(
    xOH=1e-7 / 55.2,
    xH=1e-7 / 55.2,
    xCaCO3=1e-20,
    xCa=1e-20,
    xH2CO3=1e-20,
    xHCO3=1e-20,
    xCO3=1e-20,
    thermo_config=None,
    rxn_config=None,
    has_energy_balance=True,
):
    print("==========================================================================")
    print("Case 2: Water softening via pH changes")
    print("xOH = " + str(xOH))
    print("xH = " + str(xH))
    print("xCaCO3 = " + str(xCaCO3))
    print("xCa = " + str(xCa))
    print("xH2CO3 = " + str(xH2CO3))
    print("xHCO3 = " + str(xHCO3))
    print("xCO3 = " + str(xCO3))
    print()

    model = ConcreteModel()
    model.fs = FlowsheetBlock(dynamic=False)
    model.fs.thermo_params = GenericParameterBlock(**thermo_config)

    model.fs.rxn_params = GenericReactionParameterBlock(
        property_package=model.fs.thermo_params, **rxn_config
    )

    args = {
        "property_package": model.fs.thermo_params,
        "reaction_package": model.fs.rxn_params,
        "has_rate_reactions": False,
        "has_equilibrium_reactions": True,
        "has_heat_transfer": False,
        "has_heat_of_reaction": False,
        "has_pressure_change": False,
    }
    if has_energy_balance == False:
        args["energy_balance_type"] = EnergyBalanceType.none

    model.fs.unit = EquilibriumReactor(**args)

    total_flow_mol = 10

    # Set flow_mol_phase_comp
    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Ca_2+"].fix(xCa * total_flow_mol)

    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H_+"].fix(xH * total_flow_mol)
    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "OH_-"].fix(xOH * total_flow_mol)

    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H2CO3"].fix(
        xH2CO3 * total_flow_mol
    )
    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "HCO3_-"].fix(
        xHCO3 * total_flow_mol
    )
    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "CO3_2-"].fix(
        xCO3 * total_flow_mol
    )

    model.fs.unit.inlet.flow_mol_phase_comp[0, "Sol", "CaCO3"].fix(
        xCaCO3 * total_flow_mol
    )
    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix(
        (1 - xH - xOH - xCaCO3 - xCa - xH2CO3 - xHCO3 - xCO3) * total_flow_mol
    )

    model.fs.unit.inlet.pressure.fix(101325.0)
    model.fs.unit.inlet.temperature.fix(298.0)
    if has_energy_balance == False:
        model.fs.unit.outlet.temperature.fix(298.0)

    assert degrees_of_freedom(model) == 0

    ## ==================== Start Scaling for this problem ===========================
    _set_eps_vals(model.fs.rxn_params, rxn_config)
    _set_equ_rxn_scaling(model.fs.unit, rxn_config)
    _set_mat_bal_scaling_FpcTP(model.fs.unit)
    if has_energy_balance == True:
        _set_ene_bal_scaling(model.fs.unit)

    iscale.calculate_scaling_factors(model.fs.unit)
    assert isinstance(model.fs.unit.control_volume.scaling_factor, Suffix)
    assert isinstance(
        model.fs.unit.control_volume.properties_out[0.0].scaling_factor, Suffix
    )
    assert isinstance(
        model.fs.unit.control_volume.properties_in[0.0].scaling_factor, Suffix
    )

    ## ==================== END Scaling for this problem ===========================

    # for macOS
    init_options = {**solver.options}
    init_options["tol"] = 1.0e-06
    init_options["constr_viol_tol"] = 1.0e-06
    model.fs.unit.initialize(optarg=init_options, outlvl=idaeslog.DEBUG)

    assert degrees_of_freedom(model) == 0

    results = solver.solve(model, tee=True)

    assert results.solver.termination_condition == TerminationCondition.optimal
    assert results.solver.status == SolverStatus.ok

    ## ==================== Check the results ================================
    print("comp\toutlet.tot_molfrac")
    for i in model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp:
        print(
            str(i)
            + "\t"
            + str(
                value(
                    model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp[i]
                )
            )
        )
    print()

    # NOTE: Changed all to mole fraction
    for i in model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp:
        print(
            str(i)
            + "\t"
            + str(
                value(
                    model.fs.unit.control_volume.properties_out[
                        0.0
                    ].mole_frac_phase_comp[i]
                )
            )
        )
    print()

    Ca = value(
        model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
            "Liq", "Ca_2+"
        ]
    )
    OH = value(
        model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
            "Liq", "OH_-"
        ]
    )
    H = value(
        model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
            "Liq", "H_+"
        ]
    )
    CO3 = value(
        model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
            "Liq", "CO3_2-"
        ]
    )
    Ksp = value(model.fs.unit.control_volume.reactions[0.0].k_eq["CaCO3_Ksp"].expr)

    print("Final pH = " + str(-log10(H * 55.2)))

    print()
    print("Ksp =\t" + str(Ksp))
    print("Ca*CO3 =\t" + str(Ca * CO3))
    print()
    if Ksp * 1.01 >= Ca * CO3:
        print("Constraint is satisfied!")
    else:
        print("Constraint is VIOLATED!")
        print("\tRelative error: " + str(Ksp / Ca / CO3) + ">=1")
        assert False

    print("==========================================================================")

    return model


## ================================= Case 1 Tests ===============================
@pytest.mark.component
def test_case_2_do_nothing():
    model = run_case2(
        xOH=1e-7 / 55.2,
        xH=1e-7 / 55.2,
        xCaCO3=1e-20,
        xCa=1e-20,
        xH2CO3=1e-20,
        xHCO3=1e-20,
        xCO3=1e-20,
        thermo_config=case2_thermo_config,
        rxn_config=case2_log_rxn_config,
        has_energy_balance=True,
    )


@pytest.mark.component
def test_case_2_seawater_no_ca():
    model = run_case2(
        xOH=1e-7 / 55.2,
        xH=1e-7 / 55.2,
        xCaCO3=1e-20,
        xCa=1e-20,
        xH2CO3=1e-20,
        xHCO3=0.00206 / 55.2,
        xCO3=1e-20,
        thermo_config=case2_thermo_config,
        rxn_config=case2_log_rxn_config,
        has_energy_balance=True,
    )


@pytest.mark.component
def test_case_2_seawater_with_ca():
    model = run_case2(
        xOH=1e-7 / 55.2,
        xH=1e-7 / 55.2,
        xCaCO3=1e-20,
        xCa=200 / 50000 / 2 / 55.2,
        xH2CO3=1e-20,
        xHCO3=0.00206 / 55.2,
        xCO3=1e-20,
        thermo_config=case2_thermo_config,
        rxn_config=case2_log_rxn_config,
        has_energy_balance=True,
    )


@pytest.mark.component
def test_case_2_seawater_added_carbonates():
    model = run_case2(
        xOH=1e-7 / 55.2,
        xH=1e-7 / 55.2,
        xCaCO3=1e-20,
        xCa=200 / 50000 / 2 / 55.2,
        xH2CO3=1e-20,
        xHCO3=0.00206 / 55.2,
        xCO3=0.1 / 55.2,
        thermo_config=case2_thermo_config,
        rxn_config=case2_log_rxn_config,
        has_energy_balance=True,
    )


@pytest.mark.component
def test_case_2_low_pH_no_precip():
    model = run_case2(
        xOH=1e-7 / 55.2,
        xH=1e-7 / 55.2,
        xCaCO3=1e-20,
        xCa=200 / 50000 / 2 / 55.2,
        xH2CO3=0.00206 / 55.2,
        xHCO3=1e-20,
        xCO3=1e-20,
        thermo_config=case2_thermo_config,
        rxn_config=case2_log_rxn_config,
        has_energy_balance=True,
    )


@pytest.mark.component
def test_case_2_ultra_high_ca():
    model = run_case2(
        xOH=1e-7 / 55.2,
        xH=1e-7 / 55.2,
        xCaCO3=1e-20,
        xCa=6000 / 50000 / 2 / 55.2,
        xH2CO3=1e-20,
        xHCO3=0.00206 / 55.2,
        xCO3=0.1 / 55.2,
        thermo_config=case2_thermo_config,
        rxn_config=case2_log_rxn_config,
        has_energy_balance=True,
    )


# Case 3 Config
"""
    Case 3: (softening, with lime: precipitation + dissolution)
        Aqueous Rxns:   H2O <--> H + OH             logK = -14
                        H2CO3 <--> H + HCO3         logK = -6.3
                        HCO3 <--> H + CO3           logK = -10.2
        Solubility:     CaCO3 <--> Ca + CO3         logK = -12
                        Ca(OH)2 <--> Ca + 2 OH      logK = -5.26
"""
case3_thermo_config = {
    "components": {
        "H2O": {
            "type": Solvent,
            "valid_phase_types": PT.aqueousPhase,
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (18.0153, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                "cp_mol_liq_comp_coeff": (75.312, pyunits.J / pyunits.mol / pyunits.K),
                "enth_mol_form_liq_comp_ref": (-285, pyunits.kJ / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (69, pyunits.J / pyunits.K / pyunits.mol),
            },
            # End parameter_data
        },
        "H_+": {
            "type": Cation,
            "charge": 1,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (1.00784, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                "cp_mol_liq_comp_coeff": (75.312, pyunits.J / pyunits.mol / pyunits.K),
                "enth_mol_form_liq_comp_ref": (0, pyunits.kJ / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (0, pyunits.J / pyunits.K / pyunits.mol),
            },
            # End parameter_data
        },
        "OH_-": {
            "type": Anion,
            "charge": -1,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (17.008, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                "cp_mol_liq_comp_coeff": (75.312, pyunits.J / pyunits.mol / pyunits.K),
                "enth_mol_form_liq_comp_ref": (-230, pyunits.kJ / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (
                    -10,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
            # End parameter_data
        },
        "Ca_2+": {
            "type": Cation,
            "charge": 2,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (40.078, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55, pyunits.kmol * pyunits.m**-3),
                "cp_mol_liq_comp_coeff": (167039, pyunits.J / pyunits.kmol / pyunits.K),
                "enth_mol_form_liq_comp_ref": (-542.83, pyunits.J / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (
                    -53,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
            # End parameter_data
        },
        "H2CO3": {
            "type": Solute,
            "valid_phase_types": PT.aqueousPhase,
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (62.03, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                "cp_mol_liq_comp_coeff": (75.312, pyunits.J / pyunits.mol / pyunits.K),
                "enth_mol_form_liq_comp_ref": (-699.7, pyunits.kJ / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (
                    187,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
            # End parameter_data
        },
        "HCO3_-": {
            "type": Anion,
            "charge": -1,
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (62.03, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                "cp_mol_liq_comp_coeff": (75.312, pyunits.J / pyunits.mol / pyunits.K),
                "enth_mol_form_liq_comp_ref": (-692, pyunits.kJ / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (
                    91.2,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
            # End parameter_data
        },
        "CO3_2-": {
            "type": Anion,
            "charge": -2,
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (62.03, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                "cp_mol_liq_comp_coeff": (75.312, pyunits.J / pyunits.mol / pyunits.K),
                "enth_mol_form_liq_comp_ref": (-677.1, pyunits.kJ / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (
                    -56.9,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
            # End parameter_data
        },
        "CaCO3": {
            "type": Component,
            "valid_phase_types": PT.solidPhase,
            "dens_mol_sol_comp": Constant,
            "enth_mol_sol_comp": Constant,
            "cp_mol_sol_comp": Constant,
            "entr_mol_sol_comp": Constant,
            "parameter_data": {
                "mw": (100.0869, pyunits.g / pyunits.mol),
                "dens_mol_sol_comp_coeff": (55, pyunits.kmol * pyunits.m**-3),
                "cp_mol_sol_comp_coeff": (167039, pyunits.J / pyunits.kmol / pyunits.K),
                "enth_mol_form_sol_comp_ref": (-1207, pyunits.kJ / pyunits.mol),
                "entr_mol_form_sol_comp_ref": (
                    91.7,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
        },
        "Ca(OH)2": {
            "type": Component,
            "valid_phase_types": PT.solidPhase,
            "dens_mol_sol_comp": Constant,
            "enth_mol_sol_comp": Constant,
            "cp_mol_sol_comp": Constant,
            "entr_mol_sol_comp": Constant,
            "parameter_data": {
                "mw": (74.093, pyunits.g / pyunits.mol),
                "dens_mol_sol_comp_coeff": (55, pyunits.kmol * pyunits.m**-3),
                "cp_mol_sol_comp_coeff": (167039, pyunits.J / pyunits.kmol / pyunits.K),
                "enth_mol_form_sol_comp_ref": (-986, pyunits.kJ / pyunits.mol),
                "entr_mol_form_sol_comp_ref": (83, pyunits.J / pyunits.K / pyunits.mol),
            },
        },
    },
    # End Component list
    "phases": {
        "Liq": {"type": AqueousPhase, "equation_of_state": Ideal},
        "Sol": {"type": SolidPhase, "equation_of_state": Ideal},
    },
    "state_definition": FpcTP,
    "state_bounds": {"temperature": (273.15, 300, 650), "pressure": (5e4, 1e5, 1e6)},
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
# End thermo_config definition

# Case 3 rxn config (with log_solubility_product)
case3_log_rxn_config = {
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },
    "equilibrium_reactions": {
        "CaCO3_Ksp": {
            "stoichiometry": {
                ("Sol", "CaCO3"): -1,
                ("Liq", "Ca_2+"): 1,
                ("Liq", "CO3_2-"): 1,
            },
            "heat_of_reaction": constant_dh_rxn,
            "equilibrium_constant": van_t_hoff,
            "equilibrium_form": log_solubility_product,
            "concentration_form": ConcentrationForm.moleFraction,
            "parameter_data": {
                "dh_rxn_ref": (0.0, pyunits.J / pyunits.mol),
                "k_eq_ref": (10**-12 / 55.2 / 55.2, pyunits.dimensionless),
                "T_eq_ref": (298.0, pyunits.K),
                "reaction_order": {
                    ("Sol", "CaCO3"): 0,
                    ("Liq", "Ca_2+"): 1,
                    ("Liq", "CO3_2-"): 1,
                },
            }
            # End parameter_data
        },
        "CaOH2_Ksp": {
            "stoichiometry": {
                ("Sol", "Ca(OH)2"): -1,
                ("Liq", "Ca_2+"): 1,
                ("Liq", "OH_-"): 2,
            },
            "heat_of_reaction": constant_dh_rxn,
            "equilibrium_constant": van_t_hoff,
            "equilibrium_form": log_solubility_product,
            "concentration_form": ConcentrationForm.moleFraction,
            "parameter_data": {
                "dh_rxn_ref": (0.0, pyunits.J / pyunits.mol),
                "k_eq_ref": (10**-5.26 / 55.2 / 55.2 / 55.2, pyunits.dimensionless),
                "T_eq_ref": (298.0, pyunits.K),
                "reaction_order": {
                    ("Sol", "Ca(OH)2"): 0,
                    ("Liq", "Ca_2+"): 1,
                    ("Liq", "OH_-"): 2,
                },
            }
            # End parameter_data
        },
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
                "dh_rxn_ref": (7.7, pyunits.kJ / pyunits.mol),
                "k_eq_ref": (10**-6.35 / 55.2, pyunits.dimensionless),
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
                "dh_rxn_ref": (14.9, pyunits.kJ / pyunits.mol),
                "k_eq_ref": (10**-10.33 / 55.2, pyunits.dimensionless),
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
        # End last reaction
    }
    # End equilibrium_reactions
}
# End reaction_config definition


# Defaults to pH of 7 with no lime added
def run_case3(
    xOH=1e-7 / 55.2,
    xH=1e-7 / 55.2,
    xCaCO3=1e-20,
    xCaOH2=1e-20,
    xCa=1e-20,
    xH2CO3=1e-20,
    xHCO3=1e-20,
    xCO3=1e-20,
    thermo_config=None,
    rxn_config=None,
    has_energy_balance=True,
):
    print("==========================================================================")
    print("Case 3: Water softening through explicit addition of lime")
    print("xOH = " + str(xOH))
    print("xH = " + str(xH))
    print("xCaCO3 = " + str(xCaCO3))
    print("xCaOH2 = " + str(xCaOH2))
    print("xCa = " + str(xCa))
    print("xH2CO3 = " + str(xH2CO3))
    print("xHCO3 = " + str(xHCO3))
    print("xCO3 = " + str(xCO3))
    print()

    model = ConcreteModel()
    model.fs = FlowsheetBlock(dynamic=False)
    model.fs.thermo_params = GenericParameterBlock(**thermo_config)

    model.fs.rxn_params = GenericReactionParameterBlock(
        property_package=model.fs.thermo_params, **rxn_config
    )

    args = {
        "property_package": model.fs.thermo_params,
        "reaction_package": model.fs.rxn_params,
        "has_rate_reactions": False,
        "has_equilibrium_reactions": True,
        "has_heat_transfer": False,
        "has_heat_of_reaction": False,
        "has_pressure_change": False,
    }
    if has_energy_balance == False:
        args["energy_balance_type"] = EnergyBalanceType.none

    model.fs.unit = EquilibriumReactor(**args)

    total_flow_mol = 10

    # Set flow_mol_phase_comp
    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Ca_2+"].fix(xCa * total_flow_mol)

    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H_+"].fix(xH * total_flow_mol)
    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "OH_-"].fix(xOH * total_flow_mol)

    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H2CO3"].fix(
        xH2CO3 * total_flow_mol
    )
    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "HCO3_-"].fix(
        xHCO3 * total_flow_mol
    )
    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "CO3_2-"].fix(
        xCO3 * total_flow_mol
    )

    model.fs.unit.inlet.flow_mol_phase_comp[0, "Sol", "CaCO3"].fix(
        xCaCO3 * total_flow_mol
    )
    model.fs.unit.inlet.flow_mol_phase_comp[0, "Sol", "Ca(OH)2"].fix(
        xCaOH2 * total_flow_mol
    )
    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix(
        (1 - xH - xOH - xCaCO3 - xCaOH2 - xCa - xH2CO3 - xHCO3 - xCO3) * total_flow_mol
    )

    model.fs.unit.inlet.pressure.fix(101325.0)
    model.fs.unit.inlet.temperature.fix(298.0)
    if has_energy_balance == False:
        model.fs.unit.outlet.temperature.fix(298.0)

    assert degrees_of_freedom(model) == 0

    ## ==================== Start Scaling for this problem ===========================
    #   Modify some of the default scaling factors with function args
    _set_eps_vals(model.fs.rxn_params, rxn_config, factor=1e-2, max_k_eq_ref=1e-16)
    _set_equ_rxn_scaling(model.fs.unit, rxn_config, min_k_eq_ref=1e-3)
    _set_mat_bal_scaling_FpcTP(model.fs.unit, min_flow_mol_phase_comp=1e-2)
    if has_energy_balance == True:
        _set_ene_bal_scaling(model.fs.unit)

    iscale.calculate_scaling_factors(model.fs.unit)
    assert isinstance(model.fs.unit.control_volume.scaling_factor, Suffix)
    assert isinstance(
        model.fs.unit.control_volume.properties_out[0.0].scaling_factor, Suffix
    )
    assert isinstance(
        model.fs.unit.control_volume.properties_in[0.0].scaling_factor, Suffix
    )

    ## ==================== END Scaling for this problem ===========================
    #   Loosen the tolerances for initialization stage by passing a temporary
    #       config dict replacing default values (this does so without any
    #       changes to the global solver options in WaterTAP)
    temp_config = {
        "constr_viol_tol": 1e-5,
        "tol": 1e-5,
    }
    model.fs.unit.initialize(optarg=temp_config, outlvl=idaeslog.DEBUG)

    assert degrees_of_freedom(model) == 0

    results = solver.solve(model, tee=True)

    assert results.solver.termination_condition == TerminationCondition.optimal
    assert results.solver.status == SolverStatus.ok

    ## ==================== Check the results ================================
    print("comp\toutlet.tot_molfrac")
    for i in model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp:
        print(
            str(i)
            + "\t"
            + str(
                value(
                    model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp[i]
                )
            )
        )
    print()

    # NOTE: Changed all to mole fraction
    for i in model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp:
        print(
            str(i)
            + "\t"
            + str(
                value(
                    model.fs.unit.control_volume.properties_out[
                        0.0
                    ].mole_frac_phase_comp[i]
                )
            )
        )
    print()

    Ca = value(
        model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
            "Liq", "Ca_2+"
        ]
    )
    OH = value(
        model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
            "Liq", "OH_-"
        ]
    )
    H = value(
        model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
            "Liq", "H_+"
        ]
    )
    CO3 = value(
        model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
            "Liq", "CO3_2-"
        ]
    )
    Ksp_CaCO3 = value(
        model.fs.unit.control_volume.reactions[0.0].k_eq["CaCO3_Ksp"].expr
    )
    Ksp_CaOH2 = value(
        model.fs.unit.control_volume.reactions[0.0].k_eq["CaOH2_Ksp"].expr
    )

    print("Final pH = " + str(-log10(H * 55.2)))

    print()
    print("Ksp_CaCO3 =\t" + str(Ksp_CaCO3))
    print("Ca*CO3 =\t" + str(Ca * CO3))
    print()
    if Ksp_CaCO3 * 1.01 >= Ca * CO3:
        print("Constraint is satisfied!")
    else:
        print("Constraint is VIOLATED!")
        print("\tRelative error: " + str(Ksp_CaCO3 / Ca / CO3) + ">=1")
        assert False

    print()
    print("Ksp_CaOH2 =\t" + str(Ksp_CaOH2))
    print("Ca*OH**2 =\t" + str(Ca * OH**2))
    print()
    if Ksp_CaOH2 * 1.01 >= Ca * OH**2:
        print("Constraint is satisfied!")
    else:
        print("Constraint is VIOLATED!")
        print("\tRelative error: " + str(Ksp_CaOH2 / Ca / OH**2) + ">=1")
        assert False

    print("==========================================================================")

    return model


@pytest.mark.component
def test_case_3_seawater_with_ca():
    model = run_case3(
        xOH=1e-7 / 55.2,
        xH=1e-7 / 55.2,
        xCaCO3=1e-20,
        xCaOH2=1e-20,
        xCa=200 / 50000 / 2 / 55.2,
        xH2CO3=1e-20,
        xHCO3=0.00206 / 55.2,
        xCO3=1e-20,
        thermo_config=case3_thermo_config,
        rxn_config=case3_log_rxn_config,
        has_energy_balance=True,
    )


@pytest.mark.component
def test_case_3_ultra_high_ca_forced_lime_precip():
    model = run_case3(
        xOH=1e-1 / 55.2,
        xH=1e-13 / 55.2,
        xCaCO3=1e-20,
        xCaOH2=1e-20,
        xCa=2000 / 50000 / 2 / 55.2,
        xH2CO3=1e-20,
        xHCO3=0.00206 / 55.2,
        xCO3=1e-20,
        thermo_config=case3_thermo_config,
        rxn_config=case3_log_rxn_config,
        has_energy_balance=True,
    )


@pytest.mark.component
def test_case_3_ultra_faux_added_lime_and_ash():
    added_lime_x = 100 / 50000 / 2 / 55.2
    added_ash_x = 2000 / 50000 / 2 / 55.2
    model = run_case3(
        xOH=(1e-7 / 55.2) + 2 * added_lime_x,
        xH=1e-7 / 55.2,
        xCaCO3=1e-20,
        xCaOH2=1e-20,
        xCa=(2000 / 50000 / 2 / 55.2) + added_lime_x,
        xH2CO3=1e-20,
        xHCO3=0.00206 / 55.2,
        xCO3=1e-20 + added_ash_x,
        thermo_config=case3_thermo_config,
        rxn_config=case3_log_rxn_config,
        has_energy_balance=True,
    )


# Case 4 Config
"""
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
case4_thermo_config = {
    "components": {
        "H2O": {
            "type": Solvent,
            "valid_phase_types": PT.aqueousPhase,
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (18.0153, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                "cp_mol_liq_comp_coeff": (75.312, pyunits.J / pyunits.mol / pyunits.K),
                "enth_mol_form_liq_comp_ref": (-285, pyunits.kJ / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (69, pyunits.J / pyunits.K / pyunits.mol),
            },
            # End parameter_data
        },
        "H_+": {
            "type": Cation,
            "charge": 1,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (1.00784, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                "cp_mol_liq_comp_coeff": (75.312, pyunits.J / pyunits.mol / pyunits.K),
                "enth_mol_form_liq_comp_ref": (0, pyunits.kJ / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (0, pyunits.J / pyunits.K / pyunits.mol),
            },
            # End parameter_data
        },
        "OH_-": {
            "type": Anion,
            "charge": -1,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (17.008, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                "cp_mol_liq_comp_coeff": (75.312, pyunits.J / pyunits.mol / pyunits.K),
                "enth_mol_form_liq_comp_ref": (-230, pyunits.kJ / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (
                    -10,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
            # End parameter_data
        },
        "H2CO3": {
            "type": Solute,
            "valid_phase_types": PT.aqueousPhase,
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (62.03, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                "cp_mol_liq_comp_coeff": (75.312, pyunits.J / pyunits.mol / pyunits.K),
                "enth_mol_form_liq_comp_ref": (-699.7, pyunits.kJ / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (
                    187,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
            # End parameter_data
        },
        "HCO3_-": {
            "type": Anion,
            "charge": -1,
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (62.03, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                "cp_mol_liq_comp_coeff": (75.312, pyunits.J / pyunits.mol / pyunits.K),
                "enth_mol_form_liq_comp_ref": (-692, pyunits.kJ / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (
                    91.2,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
            # End parameter_data
        },
        "CO3_2-": {
            "type": Anion,
            "charge": -2,
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (62.03, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                "cp_mol_liq_comp_coeff": (75.312, pyunits.J / pyunits.mol / pyunits.K),
                "enth_mol_form_liq_comp_ref": (-677.1, pyunits.kJ / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (
                    -56.9,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
            # End parameter_data
        },
        "H2PO4_-": {
            "type": Anion,
            "charge": -1,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (96.98724, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                "enth_mol_form_liq_comp_ref": (-1296.3, pyunits.kJ / pyunits.mol),
                "cp_mol_liq_comp_coeff": (75.312, pyunits.J / pyunits.mol / pyunits.K),
                "entr_mol_form_liq_comp_ref": (
                    90.4,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
            # End parameter_data
        },
        "HPO4_2-": {
            "type": Anion,
            "charge": -2,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (95.9793, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                "enth_mol_form_liq_comp_ref": (-1292.1, pyunits.kJ / pyunits.mol),
                "cp_mol_liq_comp_coeff": (75.312, pyunits.J / pyunits.mol / pyunits.K),
                "entr_mol_form_liq_comp_ref": (
                    -33.4,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
            # End parameter_data
        },
        "PO4_3-": {
            "type": Anion,
            "charge": -3,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (94.97136, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                "enth_mol_form_liq_comp_ref": (-1277.4, pyunits.kJ / pyunits.mol),
                "cp_mol_liq_comp_coeff": (75.312, pyunits.J / pyunits.mol / pyunits.K),
                "entr_mol_form_liq_comp_ref": (
                    -222,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
            # End parameter_data
        },
        "Fe_3+": {
            "type": Cation,
            "charge": 3,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (55.845, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                "cp_mol_liq_comp_coeff": (233467, pyunits.J / pyunits.kmol / pyunits.K),
                "enth_mol_form_liq_comp_ref": (-48.5, pyunits.kJ / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (
                    -315.9,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
            # End parameter_data
        },
        "FeOH_2+": {
            "type": Cation,
            "charge": 2,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (72.8, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                "cp_mol_liq_comp_coeff": (305000, pyunits.J / pyunits.kmol / pyunits.K),
                # NOTE: these parameters below are not well known
                "enth_mol_form_liq_comp_ref": (-229.4, pyunits.kJ / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (0, pyunits.J / pyunits.K / pyunits.mol),
            },
            # End parameter_data
        },
        "Fe(OH)2_+": {
            "type": Cation,
            "charge": 1,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (89.8, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                "cp_mol_liq_comp_coeff": (375000, pyunits.J / pyunits.kmol / pyunits.K),
                # NOTE: these parameters below are not well known
                "enth_mol_form_liq_comp_ref": (-446.7, pyunits.kJ / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (0, pyunits.J / pyunits.K / pyunits.mol),
            },
            # End parameter_data
        },
        "Fe(OH)3": {
            "type": Solute,
            "valid_phase_types": PT.aqueousPhase,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (106.8, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                "cp_mol_liq_comp_coeff": (446000, pyunits.J / pyunits.kmol / pyunits.K),
                # NOTE: these parameters below are not well known
                "enth_mol_form_liq_comp_ref": (-638.5, pyunits.kJ / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (0, pyunits.J / pyunits.K / pyunits.mol),
            },
            # End parameter_data
        },
        "Fe(OH)4_-": {
            "type": Anion,
            "charge": -1,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (123.8, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                "cp_mol_liq_comp_coeff": (518000, pyunits.J / pyunits.kmol / pyunits.K),
                # NOTE: these parameters below are not well known
                "enth_mol_form_liq_comp_ref": (-830.0, pyunits.kJ / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (0, pyunits.J / pyunits.K / pyunits.mol),
            },
            # End parameter_data
        },
        "FePO4": {
            "type": Component,
            "valid_phase_types": PT.solidPhase,
            "dens_mol_sol_comp": Constant,
            "enth_mol_sol_comp": Constant,
            "cp_mol_sol_comp": Constant,
            "entr_mol_sol_comp": Constant,
            "parameter_data": {
                "mw": (150.8, pyunits.g / pyunits.mol),
                "dens_mol_sol_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                "cp_mol_sol_comp_coeff": (635000, pyunits.J / pyunits.kmol / pyunits.K),
                "enth_mol_form_sol_comp_ref": (0, pyunits.kJ / pyunits.mol),
                "entr_mol_form_sol_comp_ref": (0, pyunits.J / pyunits.K / pyunits.mol),
            },
        },
    },
    # End Component list
    "phases": {
        "Liq": {"type": AqueousPhase, "equation_of_state": Ideal},
        "Sol": {"type": SolidPhase, "equation_of_state": Ideal},
    },
    "state_definition": FpcTP,
    "state_bounds": {"temperature": (273.15, 300, 650), "pressure": (5e4, 1e5, 1e6)},
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
# End thermo_config definition

# Case 4 rxn config (with log_solubility_product)
case4_log_rxn_config = {
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },
    "equilibrium_reactions": {
        "FePO4_Ksp": {
            "stoichiometry": {
                ("Sol", "FePO4"): -1,
                ("Liq", "Fe_3+"): 1,
                ("Liq", "PO4_3-"): 1,
            },
            "heat_of_reaction": constant_dh_rxn,
            "equilibrium_constant": van_t_hoff,
            "equilibrium_form": log_solubility_product,
            "concentration_form": ConcentrationForm.moleFraction,
            "parameter_data": {
                "dh_rxn_ref": (0.0, pyunits.J / pyunits.mol),
                "k_eq_ref": (10**-23 / 55.2 / 55.2, pyunits.dimensionless),
                "T_eq_ref": (298.0, pyunits.K),
                "reaction_order": {
                    ("Sol", "FePO4"): 0,
                    ("Liq", "Fe_3+"): 1,
                    ("Liq", "PO4_3-"): 1,
                },
            }
            # End parameter_data
        },
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
                "dh_rxn_ref": (7.7, pyunits.kJ / pyunits.mol),
                "k_eq_ref": (10**-6.35 / 55.2, pyunits.dimensionless),
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
                "dh_rxn_ref": (14.9, pyunits.kJ / pyunits.mol),
                "k_eq_ref": (10**-10.33 / 55.2, pyunits.dimensionless),
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
        },
        "H3PO4_Ka2": {
            "stoichiometry": {
                ("Liq", "H2PO4_-"): -1,
                ("Liq", "H_+"): 1,
                ("Liq", "HPO4_2-"): 1,
            },
            "heat_of_reaction": constant_dh_rxn,
            "equilibrium_constant": van_t_hoff,
            "equilibrium_form": log_power_law_equil,
            "concentration_form": ConcentrationForm.moleFraction,
            "parameter_data": {
                "dh_rxn_ref": (4.2, pyunits.kJ / pyunits.mol),
                "k_eq_ref": (10**-5.73 / 55.2, pyunits.dimensionless),
                "T_eq_ref": (298, pyunits.K),
                # By default, reaction orders follow stoichiometry
                #    manually set reaction order here to override
                "reaction_order": {
                    ("Liq", "H2PO4_-"): -1,
                    ("Liq", "H_+"): 1,
                    ("Liq", "HPO4_2-"): 1,
                },
            }
            # End parameter_data
        },
        "H3PO4_Ka3": {
            "stoichiometry": {
                ("Liq", "HPO4_2-"): -1,
                ("Liq", "H_+"): 1,
                ("Liq", "PO4_3-"): 1,
            },
            "heat_of_reaction": constant_dh_rxn,
            "equilibrium_constant": van_t_hoff,
            "equilibrium_form": log_power_law_equil,
            "concentration_form": ConcentrationForm.moleFraction,
            "parameter_data": {
                "dh_rxn_ref": (14.7, pyunits.kJ / pyunits.mol),
                "k_eq_ref": (10**-7.275 / 55.2, pyunits.dimensionless),
                "T_eq_ref": (298, pyunits.K),
                # By default, reaction orders follow stoichiometry
                #    manually set reaction order here to override
                "reaction_order": {
                    ("Liq", "HPO4_2-"): -1,
                    ("Liq", "H_+"): 1,
                    ("Liq", "PO4_3-"): 1,
                },
            }
            # End parameter_data
        },
        "FeOH_K": {
            "stoichiometry": {
                ("Liq", "FeOH_2+"): -1,
                ("Liq", "Fe_3+"): 1,
                ("Liq", "OH_-"): 1,
            },
            "heat_of_reaction": constant_dh_rxn,
            "equilibrium_constant": van_t_hoff,
            "equilibrium_form": log_power_law_equil,
            "concentration_form": ConcentrationForm.moleFraction,
            "parameter_data": {
                "dh_rxn_ref": (0, pyunits.kJ / pyunits.mol),
                "k_eq_ref": (1.768e-12 / 55.2, pyunits.dimensionless),
                "T_eq_ref": (298, pyunits.K),
                # By default, reaction orders follow stoichiometry
                #    manually set reaction order here to override
                "reaction_order": {
                    ("Liq", "FeOH_2+"): -1,
                    ("Liq", "Fe_3+"): 1,
                    ("Liq", "OH_-"): 1,
                },
            }
            # End parameter_data
        },
        "FeOH2_K": {
            "stoichiometry": {
                ("Liq", "Fe(OH)2_+"): -1,
                ("Liq", "FeOH_2+"): 1,
                ("Liq", "OH_-"): 1,
            },
            "heat_of_reaction": constant_dh_rxn,
            "equilibrium_constant": van_t_hoff,
            "equilibrium_form": log_power_law_equil,
            "concentration_form": ConcentrationForm.moleFraction,
            "parameter_data": {
                "dh_rxn_ref": (0, pyunits.kJ / pyunits.mol),
                "k_eq_ref": (3.757e-11 / 55.2, pyunits.dimensionless),
                "T_eq_ref": (298, pyunits.K),
                # By default, reaction orders follow stoichiometry
                #    manually set reaction order here to override
                "reaction_order": {
                    ("Liq", "Fe(OH)2_+"): -1,
                    ("Liq", "FeOH_2+"): 1,
                    ("Liq", "OH_-"): 1,
                },
            }
            # End parameter_data
        },
        "FeOH3_K": {
            "stoichiometry": {
                ("Liq", "Fe(OH)3"): -1,
                ("Liq", "Fe(OH)2_+"): 1,
                ("Liq", "OH_-"): 1,
            },
            "heat_of_reaction": constant_dh_rxn,
            "equilibrium_constant": van_t_hoff,
            "equilibrium_form": log_power_law_equil,
            "concentration_form": ConcentrationForm.moleFraction,
            "parameter_data": {
                "dh_rxn_ref": (0, pyunits.kJ / pyunits.mol),
                "k_eq_ref": (9.765e-7 / 55.2, pyunits.dimensionless),
                "T_eq_ref": (298, pyunits.K),
                # By default, reaction orders follow stoichiometry
                #    manually set reaction order here to override
                "reaction_order": {
                    ("Liq", "Fe(OH)3"): -1,
                    ("Liq", "Fe(OH)2_+"): 1,
                    ("Liq", "OH_-"): 1,
                },
            }
            # End parameter_data
        },
        "FeOH4_K": {
            "stoichiometry": {
                ("Liq", "Fe(OH)4_-"): -1,
                ("Liq", "Fe(OH)3"): 1,
                ("Liq", "OH_-"): 1,
            },
            "heat_of_reaction": constant_dh_rxn,
            "equilibrium_constant": van_t_hoff,
            "equilibrium_form": log_power_law_equil,
            "concentration_form": ConcentrationForm.moleFraction,
            "parameter_data": {
                "dh_rxn_ref": (0, pyunits.kJ / pyunits.mol),
                "k_eq_ref": (1.097e-6 / 55.2, pyunits.dimensionless),
                "T_eq_ref": (298, pyunits.K),
                # By default, reaction orders follow stoichiometry
                #    manually set reaction order here to override
                "reaction_order": {
                    ("Liq", "Fe(OH)4_-"): -1,
                    ("Liq", "Fe(OH)3"): 1,
                    ("Liq", "OH_-"): 1,
                },
            }
            # End parameter_data
        }
        # End last reaction
    }
    # End equilibrium_reactions
}
# End reaction_config definition

# Defaults to pH of 7
def run_case4(
    xOH=1e-7 / 55.2,
    xH=1e-7 / 55.2,
    xH2CO3=1e-20,
    xHCO3=1e-20,
    xCO3=1e-20,
    xH2PO4=1e-20,
    xHPO4=1e-20,
    xPO4=1e-20,
    xFe=1e-20,
    xFeOH=1e-20,
    xFeOH2=1e-20,
    xFeOH3=1e-20,
    xFeOH4=1e-20,
    xFePO4=1e-20,
    thermo_config=None,
    rxn_config=None,
    has_energy_balance=True,
    remove_precip_rxn=False,
):
    print("==========================================================================")
    print("Case 4: Phophorus removal through iron precipitation")
    print("xOH = " + str(xOH))
    print("xH = " + str(xH))
    print("xH2CO3 = " + str(xH2CO3))
    print("xHCO3 = " + str(xHCO3))
    print("xCO3 = " + str(xCO3))
    print("xH2PO4 = " + str(xH2PO4))
    print("xHPO4 = " + str(xHPO4))
    print("xPO4 = " + str(xPO4))

    print("xFe = " + str(xFe))
    print("xFeOH = " + str(xFeOH))
    print("xFeOH2 = " + str(xFeOH2))
    print("xFeOH3 = " + str(xFeOH3))
    print("xFeOH4 = " + str(xFeOH4))

    print("xFePO4 = " + str(xFePO4))
    print()

    Ksp_hold = rxn_config["equilibrium_reactions"]["FePO4_Ksp"]["parameter_data"][
        "k_eq_ref"
    ][0]
    if remove_precip_rxn == True:
        del rxn_config["equilibrium_reactions"]["FePO4_Ksp"]

    model = ConcreteModel()
    model.fs = FlowsheetBlock(dynamic=False)
    model.fs.thermo_params = GenericParameterBlock(**thermo_config)

    model.fs.rxn_params = GenericReactionParameterBlock(
        property_package=model.fs.thermo_params, **rxn_config
    )

    args = {
        "property_package": model.fs.thermo_params,
        "reaction_package": model.fs.rxn_params,
        "has_rate_reactions": False,
        "has_equilibrium_reactions": True,
        "has_heat_transfer": False,
        "has_heat_of_reaction": False,
        "has_pressure_change": False,
    }
    if has_energy_balance == False:
        args["energy_balance_type"] = EnergyBalanceType.none

    model.fs.unit = EquilibriumReactor(**args)

    total_flow_mol = 10

    # Set flow_mol_phase_comp
    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H_+"].fix(xH * total_flow_mol)
    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "OH_-"].fix(xOH * total_flow_mol)

    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H2CO3"].fix(
        xH2CO3 * total_flow_mol
    )
    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "HCO3_-"].fix(
        xHCO3 * total_flow_mol
    )
    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "CO3_2-"].fix(
        xCO3 * total_flow_mol
    )

    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H2PO4_-"].fix(
        xH2PO4 * total_flow_mol
    )
    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "HPO4_2-"].fix(
        xHPO4 * total_flow_mol
    )
    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "PO4_3-"].fix(
        xPO4 * total_flow_mol
    )

    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Fe_3+"].fix(xFe * total_flow_mol)
    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "FeOH_2+"].fix(
        xFeOH * total_flow_mol
    )
    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Fe(OH)2_+"].fix(
        xFeOH2 * total_flow_mol
    )
    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Fe(OH)3"].fix(
        xFeOH3 * total_flow_mol
    )
    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Fe(OH)4_-"].fix(
        xFeOH4 * total_flow_mol
    )

    model.fs.unit.inlet.flow_mol_phase_comp[0, "Sol", "FePO4"].fix(
        xFePO4 * total_flow_mol
    )
    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix(
        (
            1
            - xH
            - xOH
            - xFePO4
            - xH2CO3
            - xHCO3
            - xCO3
            - xH2PO4
            - xHPO4
            - xPO4
            - xFe
            - xFeOH
            - xFeOH2
            - xFeOH3
            - xFeOH4
        )
        * total_flow_mol
    )

    model.fs.unit.inlet.pressure.fix(101325.0)
    model.fs.unit.inlet.temperature.fix(298.0)
    if has_energy_balance == False:
        model.fs.unit.outlet.temperature.fix(298.0)

    assert degrees_of_freedom(model) == 0

    ## ==================== Start Scaling for this problem ===========================
    # # TODO: NOTE: The _set_eps_vals function may need to change. Some cases run
    #               better without it included.
    _set_eps_vals(model.fs.rxn_params, rxn_config)
    _set_equ_rxn_scaling(model.fs.unit, rxn_config)
    _set_mat_bal_scaling_FpcTP(model.fs.unit)
    if has_energy_balance == True:
        _set_ene_bal_scaling(model.fs.unit)

    iscale.calculate_scaling_factors(model.fs.unit)
    assert isinstance(model.fs.unit.control_volume.scaling_factor, Suffix)
    assert isinstance(
        model.fs.unit.control_volume.properties_out[0.0].scaling_factor, Suffix
    )
    assert isinstance(
        model.fs.unit.control_volume.properties_in[0.0].scaling_factor, Suffix
    )

    ## ==================== END Scaling for this problem ===========================
    model.fs.unit.initialize(optarg=solver.options, outlvl=idaeslog.DEBUG)

    assert degrees_of_freedom(model) == 0

    solver.options["tol"] = 1.0e-12
    results = solver.solve(model, tee=True)
    del solver.options["tol"]

    assert results.solver.termination_condition == TerminationCondition.optimal
    assert results.solver.status == SolverStatus.ok

    ## ==================== Check the results ================================
    print("comp\toutlet.tot_molfrac")
    for i in model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp:
        print(
            str(i)
            + "\t"
            + str(
                value(
                    model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp[i]
                )
            )
        )
    print()

    # NOTE: Changed all to mole fraction
    for i in model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp:
        print(
            str(i)
            + "\t"
            + str(
                value(
                    model.fs.unit.control_volume.properties_out[
                        0.0
                    ].mole_frac_phase_comp[i]
                )
            )
        )
    print()

    Fe = value(
        model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
            "Liq", "Fe_3+"
        ]
    )
    OH = value(
        model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
            "Liq", "OH_-"
        ]
    )
    H = value(
        model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
            "Liq", "H_+"
        ]
    )
    PO4 = value(
        model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
            "Liq", "PO4_3-"
        ]
    )
    if remove_precip_rxn == False:
        Ksp = value(model.fs.unit.control_volume.reactions[0.0].k_eq["FePO4_Ksp"].expr)
    else:
        Ksp = Ksp_hold

    print("Final pH = " + str(-log10(H * 55.2)))

    print()
    print("Ksp =\t" + str(Ksp))
    print("Fe*PO4 =\t" + str(Fe * PO4))
    print()
    if Ksp * 1.01 >= Fe * PO4:
        print("Constraint is satisfied!")
    else:
        print("Constraint is VIOLATED!")
        print("\tRelative error: " + str(Ksp / Fe / PO4) + ">=1")
        if remove_precip_rxn == False:
            assert False

    print("==========================================================================")

    return model


# seawater with standard amount of iron (5.38e-8 M) and phosphorus (3.22e-6 M)
#   Boost Fe, PO4, and OH to induce precipitation
#       Add 1e-4 M of each...
@pytest.mark.component
def test_case_4_seawater_added_Fe_OH_PO4():
    model = run_case4(
        xOH=(1e-7 / 55.2) + (1e-4 / 55.2),
        xH=1e-7 / 55.2,
        xH2CO3=1e-20,
        xHCO3=0.00206 / 55.2,
        xCO3=1e-20,
        xH2PO4=1e-20,
        xHPO4=(3.22e-6 / 55.2) + (1e-4 / 55.2),
        xPO4=1e-20,
        xFe=(5.38e-8 / 55.2) + (1e-4 / 55.2),
        xFeOH=1e-20,
        xFeOH2=1e-20,
        xFeOH3=1e-20,
        xFeOH4=1e-20,
        xFePO4=1e-20,
        thermo_config=case4_thermo_config,
        rxn_config=case4_log_rxn_config,
        has_energy_balance=True,
        remove_precip_rxn=False,
    )


# seawater with standard amount of iron (5.38e-8 M) and phosphorus (3.22e-6 M)
#   Boost Fe, PO4, and OH to induce precipitation
#       Add 1e-4 M of Fe
#       Add 1e-4 M of PO4
#       Add 0 M of OH
@pytest.mark.component
def test_case_4_seawater_added_Fe_PO4():
    model = run_case4(
        xOH=(1e-7 / 55.2) + (0 / 55.2),
        xH=1e-7 / 55.2,
        xH2CO3=1e-20,
        xHCO3=0.00206 / 55.2,
        xCO3=1e-20,
        xH2PO4=1e-20,
        xHPO4=(3.22e-6 / 55.2) + (1e-4 / 55.2),
        xPO4=1e-20,
        xFe=(5.38e-8 / 55.2) + (1e-4 / 55.2),
        xFeOH=1e-20,
        xFeOH2=1e-20,
        xFeOH3=1e-20,
        xFeOH4=1e-20,
        xFePO4=1e-20,
        thermo_config=case4_thermo_config,
        rxn_config=case4_log_rxn_config,
        has_energy_balance=True,
        remove_precip_rxn=False,
    )


# seawater with standard amount of iron (5.38e-8 M) and phosphorus (3.22e-6 M)
#   Boost Fe, PO4, and OH to induce precipitation
#       Add 1e-4 M of Fe
#       Add 0 M of PO4
#       Add 0 M of OH
@pytest.mark.component
def test_case_4_seawater_added_Fe_only():
    model = run_case4(
        xOH=(1e-7 / 55.2) + (0 / 55.2),
        xH=1e-7 / 55.2,
        xH2CO3=1e-20,
        xHCO3=0.00206 / 55.2,
        xCO3=1e-20,
        xH2PO4=1e-20,
        xHPO4=(3.22e-6 / 55.2) + (0 / 55.2),
        xPO4=1e-20,
        xFe=(5.38e-8 / 55.2) + (1e-4 / 55.2),
        xFeOH=1e-20,
        xFeOH2=1e-20,
        xFeOH3=1e-20,
        xFeOH4=1e-20,
        xFePO4=1e-20,
        thermo_config=case4_thermo_config,
        rxn_config=case4_log_rxn_config,
        has_energy_balance=True,
        remove_precip_rxn=False,
    )


# This is for additional testing
if __name__ == "__main__":
    # seawater with standard amount of iron (5.38e-8 M) and phosphorus (3.22e-6 M)
    #
    #       This works, but convergence of the initialization stage is VERY, VERY poor!!!
    #
    #   # TODO:  Figure out how to better initialize this case.
    # # TODO: NOTE: The _set_eps_vals function may need to change. Some cases run
    #               better without it included.

    model = run_case4(
        xOH=1e-7 / 55.2,
        xH=1e-7 / 55.2,
        xH2CO3=1e-20,
        xHCO3=0.00206 / 55.2,
        xCO3=1e-20,
        xH2PO4=1e-20,
        xHPO4=3.22e-6 / 55.2,
        xPO4=1e-20,
        xFe=5.38e-8 / 55.2,
        xFeOH=1e-20,
        xFeOH2=1e-20,
        xFeOH3=1e-20,
        xFeOH4=1e-20,
        xFePO4=1e-20,
        thermo_config=case4_thermo_config,
        rxn_config=case4_log_rxn_config,
        has_energy_balance=True,
        remove_precip_rxn=False,
    )

    exit()
