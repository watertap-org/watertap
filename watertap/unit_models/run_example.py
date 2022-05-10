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
    log10,
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
import idaes.logger as idaeslog
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


def build_ion_subset_with_alk_model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    # create dict to define ions (the prop pack of Adam requires this)
    ion_dict = {
        "solute_list": ["B[OH]3", "B[OH]4_-", "Na_+", "HCO3_-"],
        "mw_data": {"H2O": 18e-3, "B[OH]3": 61.83e-3, "B[OH]4_-": 78.83e-3, "Na_+": 23e-3, "HCO3_-": 61e-3},
        "charge": {"B[OH]3": 0, "B[OH]4_-": -1, "Na_+": 1, "HCO3_-": -1},
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

def model_setup(m, state={"H2O": 100, "H_+": 1e-7, "OH_-": 1e-7,
                          "B[OH]3": 2e-4, "B[OH]4_-": 1e-6,
                          "Na_+": 1e-3, "HCO3_-": 1e-4}):
    assert_units_consistent(m)
    #print(degrees_of_freedom(m))

    m.fs.unit.inlet.pressure.fix(101325)
    m.fs.unit.inlet.temperature.fix(298.15)
    for j in state:
        idx = (0, "Liq", j)
        if idx in m.fs.unit.inlet.flow_mol_phase_comp:
            m.fs.unit.inlet.flow_mol_phase_comp[idx].fix(state[j])
    m.fs.unit.caustic_dose.fix(5)

    if degrees_of_freedom(m) != 0:
        print(degrees_of_freedom(m))
        m.fs.unit.eq_electroneutrality.pprint()
        m.fs.unit.eq_total_boron.pprint()
        m.fs.unit.eq_water_dissociation.pprint()
        m.fs.unit.eq_boron_dissociation.pprint()
    assert degrees_of_freedom(m) == 0

def scaling_setup(m, state={"H2O": 100, "H_+": 1e-7, "OH_-": 1e-7,
                          "B[OH]3": 2e-4, "B[OH]4_-": 1e-6,
                          "Na_+": 1e-3, "HCO3_-": 1e-4}):
    # Set some scaling factors and look for 'bad' scaling
    for j in state:
        idx = (0, "Liq", j)
        if idx in m.fs.unit.inlet.flow_mol_phase_comp:
            m.fs.properties.set_default_scaling(
                "flow_mol_phase_comp", 1/state[j], index=("Liq", j)
            )

    iscale.calculate_scaling_factors(m.fs)

    # check that all constraints have scaling factors
    unscaled_con_list = list(iscale.unscaled_constraints_generator(m))
    if len(unscaled_con_list) > 0:
        for j in unscaled_con_list:
            print(j)

    # check that all variables have scaling factors
    unscaled_var_list = list(iscale.unscaled_variables_generator(m))
    if len(unscaled_var_list) > 0:
        for j in unscaled_var_list:
            print(j)
    #assert len(unscaled_var_list) == 0

    # check if any variables are badly scaled
    badly_scaled_var_values = {
        var.name: val
        for (var, val) in iscale.badly_scaled_var_generator(
            m, large=1e3, small=1e-3
        )
    }
    if len(badly_scaled_var_values) > 0:
        for j in badly_scaled_var_values:
            print(str(j) + "\t" + str(badly_scaled_var_values[j]))
    #assert not badly_scaled_var_values

def intialization_setup(m):
    m.fs.unit.initialize(optarg=solver.options, outlvl=idaeslog.DEBUG)

def solve_model(m):
    results = solver.solve(m, tee=True)

    print()

    m.fs.unit.eq_electroneutrality.pprint()
    m.fs.unit.eq_total_boron.pprint()
    m.fs.unit.eq_water_dissociation.pprint()
    m.fs.unit.eq_boron_dissociation.pprint()

    print()

    m.fs.unit.mol_H.pprint()
    m.fs.unit.mol_OH.pprint()
    m.fs.unit.mol_Boron.pprint()
    m.fs.unit.mol_Borate.pprint()

    print("conc of cation M")
    conc = m.fs.unit.caustic_cation_charge*m.fs.unit.caustic_dose[0]/m.fs.unit.caustic_mw
    print(value(conc))
    print()
    print("pH = " + str(-log10(value(m.fs.unit.mol_H[0])/1000)))

    print()

    # NOTE: May not get solved for (depends on whether or not outlet concentration is used)
    for ind in m.fs.unit.control_volume.properties_out[0.0].conc_mol_phase_comp:
        print(ind)
        print(value(m.fs.unit.control_volume.properties_out[0.0].conc_mol_phase_comp[ind]))

    print()

    m.fs.unit.eq_mass_transfer_term.pprint()

    print()

    m.fs.unit.report()

    assert_optimal_termination(results)

def display_unit_vars(m):
    ''' Used to display unit var solution '''
    return

if __name__ == "__main__":
    m = build_generic_model() # not building correctly...

    #m = build_ion_model()
    #m = build_ion_subset_model()
    #m = build_ion_subset_with_Na_model()
    #m = build_ion_subset_with_alk_model()

    model_setup(m)
    scaling_setup(m)
    intialization_setup(m)
    solve_model(m)

    # Displays standard unit result at ports
    m.fs.unit.report()

    display_unit_vars(m)

    #m.fs.unit.eq_electroneutrality.pprint()
    #m.fs.unit.eq_total_boron.pprint()
    #m.fs.unit.eq_water_dissociation.pprint()
    #m.fs.unit.eq_boron_dissociation.pprint()
