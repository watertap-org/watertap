# Importing the object for units from pyomo
from pyomo.environ import units as pyunits

# Imports from idaes core
from idaes.core import AqueousPhase, EnergyBalanceType
from idaes.core.components import Solvent, Solute, Cation, Anion
from idaes.core.phases import PhaseType as PT

# Imports from idaes generic models
import idaes.generic_models.properties.core.pure.ConstantProperties as Constant
from idaes.generic_models.properties.core.state_definitions import FpcTP
from idaes.generic_models.properties.core.eos.ideal import Ideal

# Import the idaes objects for Generic Properties and Reactions
from idaes.generic_models.properties.core.generic.generic_property import (
    GenericParameterBlock,
)
from idaes.generic_models.properties.core.generic.generic_reaction import (
    GenericReactionParameterBlock,
)

from watertap.property_models.ion_DSPMDE_prop_pack import DSPMDEParameterBlock

from watertap.unit_models.boron_removal import BoronRemoval
from pyomo.environ import (
    ConcreteModel,
    assert_optimal_termination,
    value,
    Set,
    Param,
    Var,
    units as pyunits,
    Suffix,
    Constraint,
    SolverFactory,
    SolverStatus,
    TerminationCondition,
)
from idaes.core import (
    FlowsheetBlock,
    MaterialFlowBasis,
    MaterialBalanceType,
    EnergyBalanceType,
)
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.util.check_units import assert_units_consistent
import idaes.core.util.scaling as iscale
from idaes.core.util.scaling import badly_scaled_var_generator
from idaes.core.util.testing import initialization_tester
from idaes.core.util import get_solver
import re

__author__ = "Austin Ladshaw"

solver = get_solver()


def build_ion_model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    # create dict to define ions (the prop pack of Adam requires this)
    ion_dict = {
        "solute_list": ["H_+", "OH_-", "B[OH]3", "B[OH]4_-"],
        "mw_data": {"H2O": 18e-3, "H_+": 1e-3, "OH_-": 17e-3, "B[OH]3": 61.83e-3, "B[OH]4_-": 78.83e-3},
        "charge": {"H_+": 1, "OH_-": -1, "B[OH]3": 0, "B[OH]4_-": -1,},
    }

    # attach prop pack to flowsheet
    m.fs.properties = DSPMDEParameterBlock(default=ion_dict)

    map = {'boron_name': 'B[OH]3', #[is required]
            'borate_name': 'B[OH]4_-', #[is required]
            'proton_name': 'H_+',  #[is optional]
            'hydroxide_name': 'OH_-', #[is optional]
            'caustic_additive':
                {
                    'mw_additive': (23, pyunits.g/pyunits.mol), #[is required]
                    'charge_additive': 1, #[is required]
                },
    }
    m.fs.unit = BoronRemoval(
        default={
            "property_package": m.fs.properties,
            "chemical_mapping_data": map,
        }
    )

    return m

def build_ion_subset_model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    # create dict to define ions (the prop pack of Adam requires this)
    ion_dict = {
        "solute_list": ["B[OH]3", "B[OH]4_-"],
        "mw_data": {"H2O": 18e-3, "B[OH]3": 61.83e-3, "B[OH]4_-": 78.83e-3},
        "charge": {"B[OH]3": 0, "B[OH]4_-": -1,},
    }

    # attach prop pack to flowsheet
    m.fs.properties = DSPMDEParameterBlock(default=ion_dict)

    map = {'boron_name': 'B[OH]3', #[is required]
            'borate_name': 'B[OH]4_-', #[is required]
            'caustic_additive':
                {
                    'mw_additive': (23, pyunits.g/pyunits.mol), #[is required]
                    'charge_additive': 1, #[is required]
                },
    }
    m.fs.unit = BoronRemoval(
        default={
            "property_package": m.fs.properties,
            "chemical_mapping_data": map,
        }
    )

    return m

def build_ion_subset_with_Na_model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    # create dict to define ions (the prop pack of Adam requires this)
    ion_dict = {
        "solute_list": ["B[OH]3", "B[OH]4_-", "Na_+"],
        "mw_data": {"H2O": 18e-3, "B[OH]3": 61.83e-3, "B[OH]4_-": 78.83e-3, "Na_+": 23e-3},
        "charge": {"B[OH]3": 0, "B[OH]4_-": -1, "Na_+": 1,},
    }

    # attach prop pack to flowsheet
    m.fs.properties = DSPMDEParameterBlock(default=ion_dict)

    map = {'boron_name': 'B[OH]3', #[is required]
            'borate_name': 'B[OH]4_-', #[is required]
            'caustic_additive':
                {
                    'cation_name': 'Na_+', #[is optional]
                    'mw_additive': (23, pyunits.g/pyunits.mol), #[is required]
                    'charge_additive': 1, #[is required]
                },
    }
    m.fs.unit = BoronRemoval(
        default={
            "property_package": m.fs.properties,
            "chemical_mapping_data": map,
        }
    )

    return m

def build_generic_model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    # Configuration dictionary for generic
    thermo_config = {
        "components": {
            "H2O": {
                "type": Solvent,
                # Define the methods used to calculate the following properties
                "dens_mol_liq_comp": Constant,
                "enth_mol_liq_comp": Constant,
                "cp_mol_liq_comp": Constant,
                "entr_mol_liq_comp": Constant,
                # Parameter data is always associated with the methods defined above
                "parameter_data": {
                    "mw": (18.0153, pyunits.g / pyunits.mol),
                    "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                    "cp_mol_liq_comp_coeff": (
                        75.312,
                        pyunits.J / pyunits.mol / pyunits.K,
                    ),
                    "enth_mol_form_liq_comp_ref": (0, pyunits.kJ / pyunits.mol),
                    "entr_mol_form_liq_comp_ref": (
                        0,
                        pyunits.J / pyunits.K / pyunits.mol,
                    ),
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
                    "mw": (1, pyunits.g / pyunits.mol),
                    "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                    "cp_mol_liq_comp_coeff": (
                        75.312,
                        pyunits.J / pyunits.mol / pyunits.K,
                    ),
                    "enth_mol_form_liq_comp_ref": (0, pyunits.kJ / pyunits.mol),
                    "entr_mol_form_liq_comp_ref": (
                        0,
                        pyunits.J / pyunits.K / pyunits.mol,
                    ),
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
                    "mw": (17, pyunits.g / pyunits.mol),
                    "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                    "cp_mol_liq_comp_coeff": (
                        75.312,
                        pyunits.J / pyunits.mol / pyunits.K,
                    ),
                    "enth_mol_form_liq_comp_ref": (0, pyunits.kJ / pyunits.mol),
                    "entr_mol_form_liq_comp_ref": (
                        0,
                        pyunits.J / pyunits.K / pyunits.mol,
                    ),
                },
                # End parameter_data
            },
            "B[OH]4_-": {
                "type": Anion,
                "charge": -1,
                # Define the methods used to calculate the following properties
                "dens_mol_liq_comp": Constant,
                "enth_mol_liq_comp": Constant,
                "cp_mol_liq_comp": Constant,
                "entr_mol_liq_comp": Constant,
                # Parameter data is always associated with the methods defined above
                "parameter_data": {
                    "mw": (78.83, pyunits.g / pyunits.mol),
                    "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                    "cp_mol_liq_comp_coeff": (
                        75.312,
                        pyunits.J / pyunits.mol / pyunits.K,
                    ),
                    "enth_mol_form_liq_comp_ref": (0, pyunits.kJ / pyunits.mol),
                    "entr_mol_form_liq_comp_ref": (
                        0,
                        pyunits.J / pyunits.K / pyunits.mol,
                    ),
                },
                # End parameter_data
            },
            "B[OH]3": {
                "type": Solute,
                "valid_phase_types": PT.aqueousPhase,
                # Define the methods used to calculate the following properties
                "dens_mol_liq_comp": Constant,
                "enth_mol_liq_comp": Constant,
                "cp_mol_liq_comp": Constant,
                "entr_mol_liq_comp": Constant,
                # Parameter data is always associated with the methods defined above
                "parameter_data": {
                    "mw": (61.83, pyunits.g / pyunits.mol),
                    "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                    "cp_mol_liq_comp_coeff": (
                        75.312,
                        pyunits.J / pyunits.mol / pyunits.K,
                    ),
                    "enth_mol_form_liq_comp_ref": (0, pyunits.kJ / pyunits.mol),
                    "entr_mol_form_liq_comp_ref": (
                        0,
                        pyunits.J / pyunits.K / pyunits.mol,
                    ),
                },
                # End parameter_data
            },
        },
        # End Component list
        "phases": {
            "Liq": {"type": AqueousPhase, "equation_of_state": Ideal},
        },
        "state_definition": FpcTP,
        "state_bounds": {
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

    # attach prop pack to flowsheet
    m.fs.properties = GenericParameterBlock(default=thermo_config)

    map = {'boron_name': 'B[OH]3', #[is required]
            'borate_name': 'B[OH]4_-', #[is required]
            'proton_name': 'H_+',  #[is optional]
            'hydroxide_name': 'OH_-', #[is optional]
            'caustic_additive':
                {
                    'mw_additive': (23, pyunits.g/pyunits.mol), #[is required]
                    'charge_additive': 1, #[is required]
                },
    }
    m.fs.unit = BoronRemoval(
        default={
            "property_package": m.fs.properties,
            "chemical_mapping_data": map,
        }
    )

    return m

if __name__ == "__main__":
    #m = build_generic_model()
    #m = build_ion_model()
    #m = build_ion_subset_model()
    m = build_ion_subset_with_Na_model()
