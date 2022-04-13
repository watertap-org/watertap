from pyomo.environ import ConcreteModel, SolverFactory, TerminationCondition, \
    value, Constraint, Var, Objective, Expression
from pyomo.environ import units as pyunits
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent, check_units_equivalent
import pyomo.util.infeasible as infeas
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom, number_variables
import idaes.core.util.model_statistics as stats
from idaes.core.util.constants import Constants
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.property_models.ion_DSPMDE_prop_pack import DSPMDEParameterBlock
from electrodialysis_1D import Electrodialysis1D

from idaes.core.util import get_solver

solver = get_solver()


def build_model():
    # create model, flowsheet
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    # create dict to define ions (the prop pack of Adam requires this)
    ion_dict = {
        "solute_list": ["Na_+", "Cl_-", "NaCl"],
        "diffusivity_data": {("Liq", "Na_+"): 1.33e-9,
                             ("Liq", "Cl_-"): 2.03e-9,
                             ("Liq", "NaCl"): 1.70e-9},
        "mw_data": {"H2O": 18e-3,
                    "Na_+": 23e-3,
                    "Cl_-": 35e-3,
                    "NaCl": 58e-3},
        "stokes_radius_data": {"Na_+": 0.184e-9,
                               "Cl_-": 0.121e-9,
                               "NaCl": 0.305e-9},
        "charge": {"Na_+": 1,
                   "Cl_-": -1,
                   "NaCl": 0},
    }
    # attach prop pack to flowsheet
    m.fs.properties = DSPMDEParameterBlock(default=ion_dict)

    # build the unit model, pass prop pack to the model
    m.fs.unit = Electrodialysis1D(default={
        "property_package": m.fs.properties })


    print('----------------------------------------------')
    print('DOF before specifying:', degrees_of_freedom(m.fs))

    assert_units_consistent(m)

    return m

def fix_inlets_and_vars(m):
    # specify the feed for each inlet stream
    m.fs.unit.inlet_dilute.pressure.fix(101325)
    m.fs.unit.inlet_dilute.temperature.fix(298.15)
    m.fs.unit.inlet_dilute.flow_mol_phase_comp[0, 'Liq', 'H2O'].fix(1)
    m.fs.unit.inlet_dilute.flow_mol_phase_comp[0, 'Liq', 'Na_+'].fix(0.01)
    m.fs.unit.inlet_dilute.flow_mol_phase_comp[0, 'Liq', 'Cl_-'].fix(0.01)
    m.fs.unit.inlet_dilute.flow_mol_phase_comp[0, 'Liq', 'NaCl'].fix(0.0001)

    m.fs.unit.inlet_concentrate.pressure.fix(101325)
    m.fs.unit.inlet_concentrate.temperature.fix(298.15)
    m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, 'Liq', 'H2O'].fix(1)
    m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, 'Liq', 'Na_+'].fix(0.01)
    m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, 'Liq', 'Cl_-'].fix(0.01)
    m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, 'Liq', 'NaCl'].fix(0.0001)

    m.fs.unit.cell_length.fix(1)
    m.fs.unit.cell_width.fix(0.1)
    m.fs.unit.membrane_thickness['cem'].fix()
    m.fs.unit.membrane_thickness['aem'].fix()

    print('----------------------------------------------')
    print('DOF after specifying:', degrees_of_freedom(m.fs))

    if degrees_of_freedom(m.fs) != 0:
        print("ERROR! Problem not square")
        exit(-1)

def scale_model(m):
    # set scaling factors for state vars and call the 'calculate_scaling_factors' function
    m.fs.properties.set_default_scaling('flow_mol_phase_comp', 1, index=('Liq', 'H2O'))
    m.fs.properties.set_default_scaling('flow_mol_phase_comp', 1e2, index=('Liq', 'Na_+'))
    m.fs.properties.set_default_scaling('flow_mol_phase_comp', 1e2, index=('Liq', 'Cl_-'))
    m.fs.properties.set_default_scaling('flow_mol_phase_comp', 1e4, index=('Liq', 'NaCl'))

    # NOTE: We have to skip this step for now due to an error in Adams' Prop Pack
    iscale.calculate_scaling_factors(m.fs)

def initialize_model(m):
    # Intialize the model
    m.fs.unit.initialize(optarg=solver.options, outlvl=idaeslog.DEBUG)

def run_model(m):
    # Solve the model
    results = solver.solve(m, tee=True)

def display_results(m):
    # Display results 'cleanly'
    m.fs.unit.report()

    # Display full set of model info on the unit
    #m.fs.unit.pprint()

def view_model_constraints(m):
    # Display the material balance constraints
    m.fs.unit.dilute_side.material_balances.pprint()
    m.fs.unit.dilute_side.material_flow_dx_disc_eq.pprint()
    m.fs.unit.dilute_side.mass_transfer_term.pprint()

    print()

    # Display the mass transfer terms
    m.fs.unit.concentrate_side.material_balances.pprint()
    m.fs.unit.concentrate_side.material_flow_dx_disc_eq.pprint()
    m.fs.unit.concentrate_side.mass_transfer_term.pprint()

def view_model_control_volumes(m):
    # Display the full control volume equation set
    m.fs.unit.dilute_side.pprint()
    print()
    m.fs.unit.concentrate_side.pprint()

def view_model_properties(m):
    # Display the full control volume equation set
    m.fs.unit.dilute_side.properties.pprint()
    print()
    m.fs.unit.concentrate_side.properties.pprint()

if __name__ == "__main__":
   m = build_model()
   fix_inlets_and_vars(m)
   scale_model(m)

   # ONLY use these line to solve the model (requires proper solvers)
   initialize_model(m)
   run_model(m)
   display_results(m)


   #view_model_constraints(m)
   #view_model_control_volumes(m)
   #view_model_properties(m)
