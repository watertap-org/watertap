"""
Carbonic acid dissociation in water
"""
from pyomo.environ import (Block,
                           SolverFactory,
                           ConcreteModel,
                           Set,
                           SolverStatus,
                           TerminationCondition,
                           value)
import idaes.logger as idaeslog
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              fixed_variables_set,
                                              activated_constraints_set)

# Import the idaes objects for Generic Properties and Reactions
from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterBlock)
from idaes.generic_models.properties.core.generic.generic_reaction import (
        GenericReactionParameterBlock)

# Import the idaes object for the EquilibriumReactor unit model
from idaes.generic_models.unit_models.equilibrium_reactor import EquilibriumReactor
#from idaes.generic_models.properties.core.pure.Perrys import *

from pyomo.environ import units as pyunits
from idaes.core import LiquidPhase, VaporPhase, Component
from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.generic_models.properties.core.eos.ideal import Ideal

from proteuslib.edb.db_api import ElectrolyteDB


def get_configs(component_names):
    db = ElectrolyteDB(db="edb")
    thermo_base = next(db.get_base("thermo"))
    result = db.get_components(component_names)
    for comp in result:
        thermo_base.add(comp)
    water_reaction_base = next(db.get_base("water_reaction"))
    for react in db.get_reactions(component_names):
        water_reaction_base.add(react)

    return {"thermo_config": thermo_base.config, "reaction_config": water_reaction_base.config}

    # # Pass a list of the thermo config dictionaries to the stitcher
    # # First argument is a starter dict, second arg is a list of component dicts
    # thermo_config = stitch_thermo_configs(starter_thermo_config,
    #                                         [H2O_thermo_config,
    #                                         H_thermo_config,
    #                                         OH_thermo_config,
    #                                         H2CO3_thermo_config,
    #                                         HCO3_thermo_config,
    #                                         CO3_thermo_config])
    #
    # # Pass a list of the reaction config dictionaries to the stitcher
    # reaction_config = stitch_reaction_configs([water_reaction_config,
    #                                             carbonic_acid_reaction_config,
    #                                             bicarbonate_reaction_config
    #                                             ])
    #return thermo_config, reaction


def create_model(thermo_config=None, reaction_config=None):
    # Create a pyomo model object
    model = ConcreteModel()

    # Add an idaes flowsheet object to that model
    model.fs = FlowsheetBlock(default={"dynamic": False})

    # Add the generic properties object to that flowsheet
    #       and pass the custom configuation dictionary to it
    model.fs.thermo_params = GenericParameterBlock(default=thermo_config)

    # Add the generic reactions object to that flowsheet
    #   and pass the corresponding properties/params object to it
    #   along with the reaction_config dictionary
    model.fs.rxn_params = GenericReactionParameterBlock(
                default={"property_package": model.fs.thermo_params, **reaction_config})

    # Add the unit model for an Equilibrium Reactor to the flowsheet
    model.fs.unit = EquilibriumReactor(default={
                "property_package": model.fs.thermo_params,
                "reaction_package": model.fs.rxn_params,
                "has_rate_reactions": False,
                "has_equilibrium_reactions": True,
                "has_heat_transfer": False,
                "has_heat_of_reaction": False,
                "has_pressure_change": False})


    ## Use this function to ensure the units and unit conversions are all correct
    ##   This will throw an error if things are wrong
    #assert_units_consistent(model)

    print("Degrees of freedom = " + str(degrees_of_freedom(model) ) )

    # After building this model, there will be 9 degrees of freedom that MUST be
    #   fixed before calling the solver. Those degrees of freedom are based on
    #   the 'state_definition' defined in the 'thermo_config' dictionary and
    #   the number of components declared.
    #
    # For our problem:
    # ---------------
    #   'state_definition' == 'FTPx'
    #           Corresponds to:
    #           ---------------
    #                   model.fs.unit.inlet.flow_mol    -   mol/s
    #                   model.fs.unit.inlet.pressure    -   Pa
    #                   model.fs.unit.inlet.temperature -   K
    #                   model.fs.unit.inlet.mole_frac_comp[0, "name_of_component"]
    #                       (This one is repeated for each named component in 'thermo_config')
    #                       (The first value is the time stamp associated with the inlet,
    #                           which in our case is 0 because we are not a dynamic flowsheet)
    #
    #   'components' == ['H2O', 'H +', 'OH -', 'H2CO3', 'HCO3 -', 'CO3 2-']
    #           Corresponds to:
    #           ---------------
    #                   model.fs.unit.inlet.mole_frac_comp[0, "H2O"]    -   unitless
    #                   model.fs.unit.inlet.mole_frac_comp[0, "CO2"]    -   unitless
    #                       (Sum of all mole_frac_comp should equal 1.)
    model.fs.unit.inlet.flow_mol.fix(10)
    model.fs.unit.inlet.pressure.fix(101325.0)
    model.fs.unit.inlet.temperature.fix(300.0)

    # Set molefractions for ions... (this may change later)
    model.fs.unit.inlet.mole_frac_comp[0, "H +"].fix( 0. )
    model.fs.unit.inlet.mole_frac_comp[0, "OH -"].fix( 0. )
    model.fs.unit.inlet.mole_frac_comp[0, "HCO3 -"].fix( 0. )
    model.fs.unit.inlet.mole_frac_comp[0, "CO3 2-"].fix( 0. )

    # To come up with molefractions, we can start be pulling the water density from
    #   our thermo_config as follows (NOTE: it will have same units as before)
    #       [This approximation is valid for all aqueous chemistry]
    #
    #       Units:  kmol/m**3 == mol/L
    total_carbonate = 0.00206
    total_molar_density = calc_water_density(300)

    model.fs.unit.inlet.mole_frac_comp[0, "H2CO3"].fix( total_carbonate/total_molar_density )

    # Perform a summation of all non-H2O molefractions to find the H2O molefraction
    sum = 0
    for i in model.fs.unit.inlet.mole_frac_comp:
        # NOTE: i will be a tuple with format (time, component)
        if i[1] != "H2O":
            sum += value(model.fs.unit.inlet.mole_frac_comp[i[0], i[1]])

    # Fix the inlet H2O molefraction based on remaining fraction
    model.fs.unit.inlet.mole_frac_comp[0, "H2O"].fix(1-sum)

    #Setup scaling factors (optional)

    model.fs.thermo_params.default_scaling_factor[("mole_frac_comp", "H +")] = 1e7
    model.fs.thermo_params.default_scaling_factor[("mole_frac_comp", "OH -")] = 1e7
    model.fs.thermo_params.default_scaling_factor[("mole_frac_comp", "HCO3 -")] = 1e7
    model.fs.thermo_params.default_scaling_factor[("mole_frac_comp", "CO3 2-")] = 1e7
    model.fs.thermo_params.default_scaling_factor[("mole_frac_comp", "H2CO3")] = 1e7
    model.fs.thermo_params.default_scaling_factor[("mole_frac_comp", "H2O")] = 1

    # # TODO: ERROR! Cannot use this simple method to calculate scaling factors
    #           Model does not like it when I do this. I believe I have to give
    #           scaling factors for everything in the model for it to work properly.
    #iscale.calculate_scaling_factors(model.fs.unit)

    print("Degrees of freedom = " + str(degrees_of_freedom(model) ) )

    return model


def solve_model(model):
    # At this point, all degrees of freedom are fixed, and we can solve the model
    solver = SolverFactory('ipopt')
    #model.fs.unit.initialize(optarg=solver.options, outlvl=idaeslog.DEBUG)
    model.fs.unit.initialize(optarg=solver.options)

    solver.options['bound_push'] = 1e-10
    results = solver.solve(model, tee=True)

    # Extract solution information
    model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp.display()
    model.fs.unit.control_volume.properties_out[0.0].dens_mol_phase.display()
    model.fs.unit.control_volume.properties_out[0.0].conc_mol_phase_comp.display()
    total_molar_density = value(model.fs.unit.control_volume.properties_out[0.0].dens_mol_phase['Liq'])/1000

    print("Mole fraction H + =\t"+ str(value(model.fs.unit.outlet.mole_frac_comp[0, "H +"])))
    print("Conc H + =\t"+ str(value(model.fs.unit.outlet.mole_frac_comp[0, "H +"]*total_molar_density*1000)))
    pHo = -log10(value(model.fs.unit.outlet.mole_frac_comp[0, "H +"])*total_molar_density)
    pOHo = -log10(value(model.fs.unit.outlet.mole_frac_comp[0, "OH -"])*total_molar_density)
    print("Outlet pH =\t"+ str(pHo) )
    print("Outlet pOH =\t"+ str(pOHo) )

    print("inlet.temperature\toutlet.temperature")
    print(str(value(model.fs.unit.inlet.temperature[0]))+"\t"+str(value(model.fs.unit.outlet.temperature[0]))+"\n")


def main():
    from pprint import pprint
    component_names = ["H +", "H2CO3", "HCO3 -"]
    configs = get_configs(component_names)
    for key, value in configs.items():
        print(f"## {key}")
        pprint(value)
    #model = create_model(**configs)
    #solve_model(model)


if __name__ == "__main__":
    main()
