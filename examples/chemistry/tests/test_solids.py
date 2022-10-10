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
    and usage of solids phases in conjunction with aqueous phases. This test
    is primarily being used to probe for issues that may exist in how IDAES
    handles this new system and properties, since it is brand new to the
    framework. Several tests of the same basic case are observed.

    Case 1: [Combine] A and B are aqueous, declare AB a solid,
                    add precipitation reaction (record IDAES debug)

    Several studies of this case are repeated for numerous inlet values to
    assess IDAES's ability to capture both precipitation and dissolution

    Solubility Reaction:
        AB <---> A + B

    Case 2: Repeat Case 1, but try using new 'log_solubility_product'
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
                "enth_mol_form_liq_comp_ref": (0, pyunits.kJ / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (0, pyunits.J / pyunits.K / pyunits.mol),
            },
            # End parameter_data
        },
        "A": {
            "type": Solute,
            "valid_phase_types": PT.aqueousPhase,
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            "parameter_data": {
                "mw": (18, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                "cp_mol_liq_comp_coeff": (75.312, pyunits.J / pyunits.mol / pyunits.K),
                "enth_mol_form_liq_comp_ref": (0, pyunits.kJ / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (0, pyunits.J / pyunits.K / pyunits.mol),
            },
        },
        "B": {
            "type": Solute,
            "valid_phase_types": PT.aqueousPhase,
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            "parameter_data": {
                "mw": (18, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                "cp_mol_liq_comp_coeff": (75.312, pyunits.J / pyunits.mol / pyunits.K),
                "enth_mol_form_liq_comp_ref": (0, pyunits.kJ / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (0, pyunits.J / pyunits.K / pyunits.mol),
            },
        },
        "AB": {
            "type": Component,
            "valid_phase_types": PT.solidPhase,
            "dens_mol_sol_comp": Constant,
            "enth_mol_sol_comp": Constant,
            "cp_mol_sol_comp": Constant,
            "entr_mol_sol_comp": Constant,
            "parameter_data": {
                "mw": (18, pyunits.g / pyunits.mol),
                "dens_mol_sol_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                "cp_mol_sol_comp_coeff": (75.312, pyunits.J / pyunits.mol / pyunits.K),
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
    # Default for testing = FpcTP
    "state_definition": FpcTP,
    # "state_definition": FTPx,
    "state_bounds": {
        # "flow_mol": (0, 50, 100),
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
# End thermo_config definition

# Actual solubility product
reaction_solubility = {
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },
    "equilibrium_reactions": {
        "AB_Ksp": {
            "stoichiometry": {("Sol", "AB"): -1, ("Liq", "A"): 1, ("Liq", "B"): 1},
            "heat_of_reaction": constant_dh_rxn,
            "equilibrium_constant": van_t_hoff,
            "equilibrium_form": solubility_product,
            "concentration_form": ConcentrationForm.moleFraction,
            "parameter_data": {
                "dh_rxn_ref": (0.0, pyunits.J / pyunits.mol),
                "k_eq_ref": (10**-10, pyunits.dimensionless),
                "T_eq_ref": (300.0, pyunits.K),
                "reaction_order": {("Sol", "AB"): 0, ("Liq", "A"): 1, ("Liq", "B"): 1},
            }
            # End parameter_data
        }
        # End Reaction
    }
    # End equilibrium_reactions
}
# End reaction_config definition

# Actual solubility product
reaction_log_solubility = {
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },
    "equilibrium_reactions": {
        "AB_Ksp": {
            "stoichiometry": {("Sol", "AB"): -1, ("Liq", "A"): 1, ("Liq", "B"): 1},
            "heat_of_reaction": constant_dh_rxn,
            "equilibrium_constant": van_t_hoff,
            "equilibrium_form": log_solubility_product,
            "concentration_form": ConcentrationForm.moleFraction,
            "parameter_data": {
                "dh_rxn_ref": (0.0, pyunits.J / pyunits.mol),
                "k_eq_ref": (10**-10, pyunits.dimensionless),
                "T_eq_ref": (300.0, pyunits.K),
                "reaction_order": {("Sol", "AB"): 0, ("Liq", "A"): 1, ("Liq", "B"): 1},
            }
            # End parameter_data
        }
        # End Reaction
    }
    # End equilibrium_reactions
}
# End reaction_config definition

# Get default solver for testing
solver = get_solver()


def run_case1(xA, xB, xAB=1e-25, scaling=True, rxn_config=None):
    print("==========================================================================")
    print("Case 1: A and B are aqueous, AB is solid that forms from reaction")
    print("xA = " + str(xA))
    print("xB = " + str(xB))
    print("xAB = " + str(xAB))
    print("scaling = " + str(scaling))
    print("including water = " + str(True))
    print()
    model = ConcreteModel()
    model.fs = FlowsheetBlock(dynamic=False)
    model.fs.thermo_params = GenericParameterBlock(**case1_thermo_config)

    model.fs.rxn_params = GenericReactionParameterBlock(
        property_package=model.fs.thermo_params, **rxn_config
    )

    model.fs.unit = EquilibriumReactor(
        property_package=model.fs.thermo_params,
        reaction_package=model.fs.rxn_params,
        has_rate_reactions=False,
        has_equilibrium_reactions=True,
        has_heat_transfer=False,
        has_heat_of_reaction=False,
        has_pressure_change=False,
        energy_balance_type=EnergyBalanceType.none,
    )

    total_flow_mol = 10

    # Set flow_mol_phase_comp
    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "A"].fix(xA * total_flow_mol)
    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "B"].fix(xB * total_flow_mol)
    model.fs.unit.inlet.flow_mol_phase_comp[0, "Sol", "AB"].fix(xAB * total_flow_mol)
    model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix(
        (1 - xA - xB - xAB) * total_flow_mol
    )

    model.fs.unit.inlet.pressure.fix(101325.0)
    model.fs.unit.inlet.temperature.fix(298.0)
    model.fs.unit.outlet.temperature.fix(298.0)

    assert degrees_of_freedom(model) == 0

    assert_units_consistent(model)

    # Scaling
    _set_eps_vals(model.fs.rxn_params, rxn_config)
    # NOTE: We skip reaction scaling because we are NOT using the log_solubility_product form in this test
    # _set_equ_rxn_scaling(model.fs.unit, rxn_config)
    _set_mat_bal_scaling_FpcTP(model.fs.unit)

    iscale.calculate_scaling_factors(model.fs.unit)
    assert isinstance(model.fs.unit.control_volume.scaling_factor, Suffix)
    assert isinstance(
        model.fs.unit.control_volume.properties_out[0.0].scaling_factor, Suffix
    )
    assert isinstance(
        model.fs.unit.control_volume.properties_in[0.0].scaling_factor, Suffix
    )
    # End scaling if statement

    solver.options["max_iter"] = 200
    init_options = {**solver.options}
    init_options["bound_relax_factor"] = 1.0e-02
    model.fs.unit.initialize(optarg=init_options, outlvl=idaeslog.DEBUG)

    assert degrees_of_freedom(model) == 0

    solver.options["bound_relax_factor"] = 1.0e-02
    results = solver.solve(model, tee=True)
    del solver.options["bound_relax_factor"]

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

    A = value(
        model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
            "Liq", "A"
        ]
    )
    B = value(
        model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
            "Liq", "B"
        ]
    )
    Ksp = value(model.fs.unit.control_volume.reactions[0.0].k_eq["AB_Ksp"].expr)

    print()
    if Ksp * 1.01 >= A * B:
        print("Constraint is satisfied!")
    else:
        print("Constraint is VIOLATED!")
        print("\tRelative error: " + str(Ksp / A / B) + ">=1")
        assert False
    print("Ksp =\t" + str(Ksp))
    print("A*B =\t" + str(A * B))

    print("==========================================================================")

    return model


def run_case2(xA, xB, xAB=1e-25, scaling=True, rxn_config=None, state="FpcTP"):
    print("==========================================================================")
    print(
        "Case 2 (log form): A and B are aqueous, AB is solid that forms from reaction"
    )
    print("xA = " + str(xA))
    print("xB = " + str(xB))
    print("xAB = " + str(xAB))
    print("scaling = " + str(scaling))
    print("including water = " + str(True))
    print()
    model = ConcreteModel()
    model.fs = FlowsheetBlock(dynamic=False)

    if state == "FpcTP":
        case1_thermo_config["state_definition"] = FpcTP
    elif state == "FTPx":
        case1_thermo_config["state_definition"] = FTPx
        case1_thermo_config["state_bounds"]["flow_mol"] = (0, 50, 100)
    else:
        print("Error! Undefined state...")
        assert False

    model.fs.thermo_params = GenericParameterBlock(**case1_thermo_config)

    model.fs.rxn_params = GenericReactionParameterBlock(
        property_package=model.fs.thermo_params, **rxn_config
    )

    model.fs.unit = EquilibriumReactor(
        property_package=model.fs.thermo_params,
        reaction_package=model.fs.rxn_params,
        has_rate_reactions=False,
        has_equilibrium_reactions=True,
        has_heat_transfer=False,
        has_heat_of_reaction=False,
        has_pressure_change=False,
        energy_balance_type=EnergyBalanceType.none,
    )

    total_flow_mol = 10

    # Set flow_mol_phase_comp
    if case1_thermo_config["state_definition"] == FpcTP:

        model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "A"].fix(xA * total_flow_mol)
        model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "B"].fix(xB * total_flow_mol)
        model.fs.unit.inlet.flow_mol_phase_comp[0, "Sol", "AB"].fix(
            xAB * total_flow_mol
        )
        model.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix(
            (1 - xA - xB - xAB) * total_flow_mol
        )

    if case1_thermo_config["state_definition"] == FTPx:
        model.fs.unit.inlet.mole_frac_comp[0, "A"].fix(xA)
        model.fs.unit.inlet.mole_frac_comp[0, "B"].fix(xB)
        model.fs.unit.inlet.mole_frac_comp[0, "AB"].fix(xAB)
        model.fs.unit.inlet.mole_frac_comp[0, "H2O"].fix((1 - xA - xB - xAB))
        model.fs.unit.inlet.flow_mol.fix(total_flow_mol)

    model.fs.unit.inlet.pressure.fix(101325.0)
    model.fs.unit.inlet.temperature.fix(298.0)
    model.fs.unit.outlet.temperature.fix(298.0)

    assert degrees_of_freedom(model) == 0

    assert_units_consistent(model)

    # Scaling
    _set_eps_vals(model.fs.rxn_params, rxn_config)
    _set_equ_rxn_scaling(model.fs.unit, rxn_config)
    if case1_thermo_config["state_definition"] == FpcTP:
        _set_mat_bal_scaling_FpcTP(model.fs.unit)
    if case1_thermo_config["state_definition"] == FTPx:
        _set_mat_bal_scaling_FTPx(model.fs.unit)

    iscale.calculate_scaling_factors(model.fs.unit)
    assert isinstance(model.fs.unit.control_volume.scaling_factor, Suffix)

    # Initialize model

    if case1_thermo_config["state_definition"] == FpcTP:
        state_args = {
            "flow_mol_phase_comp": {
                ("Liq", "H2O"): model.fs.unit.inlet.flow_mol_phase_comp[
                    0, "Liq", "H2O"
                ].value,
                ("Liq", "A"): model.fs.unit.inlet.flow_mol_phase_comp[
                    0, "Liq", "A"
                ].value,
                ("Liq", "B"): model.fs.unit.inlet.flow_mol_phase_comp[
                    0, "Liq", "B"
                ].value,
                ("Sol", "AB"): model.fs.unit.inlet.flow_mol_phase_comp[
                    0, "Sol", "AB"
                ].value,
            },
            "pressure": 101325,
            "temperature": 298,
            "flow_mol": 10,
        }

    if case1_thermo_config["state_definition"] == FTPx:
        state_args = {
            "mole_frac_comp": {
                "H2O": model.fs.unit.inlet.mole_frac_comp[0, "H2O"].value,
                "A": model.fs.unit.inlet.mole_frac_comp[0, "A"].value,
                "B": model.fs.unit.inlet.mole_frac_comp[0, "B"].value,
                "AB": model.fs.unit.inlet.mole_frac_comp[0, "AB"].value,
            },
            "pressure": 101325,
            "temperature": 298,
            "flow_mol": 10,
        }

    flags = fix_state_vars(model.fs.unit.control_volume.properties_out, state_args)
    revert_state_vars(model.fs.unit.control_volume.properties_out, flags)

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

    A = value(
        model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
            "Liq", "A"
        ]
    )
    B = value(
        model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
            "Liq", "B"
        ]
    )
    Ksp = value(model.fs.unit.control_volume.reactions[0.0].k_eq["AB_Ksp"].expr)

    print()
    if Ksp * 1.01 >= A * B:
        print("Constraint is satisfied!")
    else:
        print("Constraint is VIOLATED!")
        print("\tRelative error: " + str(Ksp / A / B) + ">=1")
        assert False
    print("Ksp =\t" + str(Ksp))
    print("A*B =\t" + str(A * B))

    print("==========================================================================")

    return model


## ================================= Case 1 Tests ===============================
@pytest.mark.component
def test_case1_low_conc_no_precipitation():
    model = run_case1(
        xA=1e-9, xB=1e-9, xAB=1e-25, scaling=True, rxn_config=reaction_solubility
    )


@pytest.mark.component
def test_case1_mid_conc_no_precipitation():
    model = run_case1(
        xA=1e-9, xB=1e-2, xAB=1e-25, scaling=True, rxn_config=reaction_solubility
    )


@pytest.mark.component
def test_case1_high_conc_with_precipitation():
    model = run_case1(
        xA=1e-2, xB=1e-2, xAB=1e-25, scaling=True, rxn_config=reaction_solubility
    )


@pytest.mark.component
def test_case1_low_conc_with_dissolution():
    model = run_case1(
        xA=1e-9, xB=1e-9, xAB=1e-2, scaling=True, rxn_config=reaction_solubility
    )


@pytest.mark.component
def test_case1_mid_conc_with_dissolution():
    model = run_case1(
        xA=1e-9, xB=1e-2, xAB=1e-2, scaling=True, rxn_config=reaction_solubility
    )


@pytest.mark.component
def test_case1_high_conc_for_all():
    model = run_case1(
        xA=1e-2, xB=1e-2, xAB=1e-2, scaling=True, rxn_config=reaction_solubility
    )


## ================================= Case 2 Tests ===============================
@pytest.mark.component
def test_case2_low_conc_no_precipitation():
    model = model = run_case2(
        xA=1e-9, xB=1e-9, xAB=1e-25, scaling=True, rxn_config=reaction_log_solubility
    )


@pytest.mark.component
def test_case2_mid_conc_no_precipitation():
    model = model = run_case2(
        xA=1e-9, xB=1e-2, xAB=1e-25, scaling=True, rxn_config=reaction_log_solubility
    )


@pytest.mark.component
def test_case2_high_conc_with_precipitation():
    model = model = run_case2(
        xA=1e-2, xB=1e-2, xAB=1e-25, scaling=True, rxn_config=reaction_log_solubility
    )


@pytest.mark.component
def test_case2_low_conc_with_dissolution():
    model = run_case2(
        xA=1e-9, xB=1e-9, xAB=1e-2, scaling=True, rxn_config=reaction_log_solubility
    )


@pytest.mark.component
def test_case2a_mid_conc_with_dissolution():
    model = model = run_case2(
        xA=1e-9, xB=1e-2, xAB=1e-2, scaling=True, rxn_config=reaction_log_solubility
    )


@pytest.mark.component
def test_case2b_mid_conc_with_dissolution():
    model = model = run_case2(
        xA=1e-2, xB=1e-9, xAB=1e-2, scaling=True, rxn_config=reaction_log_solubility
    )


@pytest.mark.component
def test_case2_high_conc_for_all():
    model = model = run_case2(
        xA=1e-2, xB=1e-2, xAB=1e-2, scaling=True, rxn_config=reaction_log_solubility
    )
