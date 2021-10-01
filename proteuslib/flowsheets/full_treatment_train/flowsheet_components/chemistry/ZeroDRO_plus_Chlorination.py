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
    Simple example of a flowsheet containing an RO separator unit model and
    a simple NaOCl chlorination post-treatment unit model.

    inlet ---> [  0D RO  ] ---> permeate ---> (Translator) ---> [Chlorination] ---> outlet
                    |
                    |
                    v
                retentate (i.e., waste)


    NOTE: The 2 unit models use a different set of state_vars. Thus, this will need to be
    resolved with some clever constraint formulation.

    Both inlet and outlet streams use K for temperature and Pa for pressure (no change needed)

    The flow from RO Separator uses kg/s for individual "species" [H2O and TDS]

    The inlet for Chlorination uses a total molar flow rate in mol/s and mole fractions
    of individual species. To make the appropriate conversions, we will have to start
    by making some assumptions about the molecular weight of TDS.

    MW H2O = 18e-3 kg/mol   MW TDS = 58.4e-3 kg/mol
                            (just assume all as NaCl? : MW Na = 23 g/mol MW Cl = 35.4 g/mol)

    Total Molar Flow = [ m.fs.RO.permeate.flow_mass_phase_comp[0, 'Liq', 'H2O']/(MW H2O) +
                            m.fs.RO.permeate.flow_mass_phase_comp[0, 'Liq', 'TDS']/(MW TDS) ]

    Molefraction of Na --> Based on TDS
                =  [m.fs.RO.permeate.flow_mass_phase_comp[0, 'Liq', 'TDS']/(MW TDS)] / (Total Molar Flow)
    Molefraction of Cl = Molefraction of Na (1:1 ratio in the salt)
    Molefraction of H2O --> Whatever is remaining

    ---------- NOTE: This is only an example ---------
"""

from proteuslib.flowsheets.full_treatment_train.model_components import property_models
from proteuslib.flowsheets.full_treatment_train.flowsheet_components import desalination
from proteuslib.flowsheets.full_treatment_train.flowsheet_components.chemistry.PostTreatment_SimpleNaOCl_Chlorination import (
    build_simple_naocl_chlorination_unit,
    initialize_chlorination_example,
    display_results_of_chlorination,
    simple_naocl_reaction_config)

# Import specific pyomo objects
from pyomo.environ import (ConcreteModel,
                           Constraint,
                           Expression,
                           Objective,
                           TransformationFactory,
                           value)

# Import the core idaes objects for Flowsheets and types of balances
from idaes.core import FlowsheetBlock

from idaes.generic_models.unit_models.translator import Translator
from pyomo.network import Arc

from proteuslib.flowsheets.full_treatment_train.util import solve_with_user_scaling, check_dof

from proteuslib.flowsheets.full_treatment_train.electrolyte_scaling_utils import (
    approximate_chemical_state_args,
    calculate_chemical_scaling_factors)

from proteuslib.flowsheets.full_treatment_train.chemical_flowsheet_util import seq_decomp_initializer

from idaes.core.util import scaling as iscale

from idaes.core.util import get_solver

from idaes.core.util.initialization import propagate_state

__author__ = "Austin Ladshaw"

# Get default solver for testing
solver = get_solver()

def build_0DRO_Chlorination_flowsheet(model, mg_per_L_NaOCl_added=0, RO_level='detailed'):
    property_models.build_prop(model, base='TDS')

    # Here, we set 'has_feed' to True because RO is our first block in the flowsheet
    kwargs_desal = {'has_desal_feed': True, 'is_twostage': False, 'has_ERD': False,
                    'RO_type': '0D', 'RO_base': 'TDS', 'RO_level': RO_level}
    desal_port = desalination.build_desalination(model, **kwargs_desal)
    desalination.scale_desalination(model, **kwargs_desal)
    desalination.initialize_desalination(model, **kwargs_desal)

    # You can change some RO unit defaults here
    #model.fs.RO.area.set_value(100)

    #unit_0DRO.build_RO(model, base='TDS', level='simple')
    #property_models.specify_feed(model.fs.RO.feed_side.properties_in[0], base='TDS')

    total_molar_density = 1/18*1000 #mol/L
    free_chlorine_added = mg_per_L_NaOCl_added/74.44/1000*70900 #mg/L as NaOCl
    total_chlorine_inlet = free_chlorine_added/70900 # mol/L
    total_molar_density+=total_chlorine_inlet

    # May need to change this build interface
    build_simple_naocl_chlorination_unit(model, mg_per_L_NaOCl_added = mg_per_L_NaOCl_added)

    # Translator inlet from RO and outlet goes to chlorination
    # NOTE: May need to come up with a way to set state_args for Translator for
    #       better convergence behavior. This block seems to be the trouble maker
    #       for the full solve.
    model.fs.RO_to_Chlor = Translator(
        default={"inlet_property_package": model.fs.prop_TDS,
                 "outlet_property_package": model.fs.simple_naocl_thermo_params})

    # Add constraints to define how the translator will function
    model.fs.RO_to_Chlor.eq_equal_temperature = Constraint(
        expr=model.fs.RO_to_Chlor.inlet.temperature[0]
        == model.fs.RO_to_Chlor.outlet.temperature[0])
    model.fs.RO_to_Chlor.eq_equal_pressure = Constraint(
        expr=model.fs.RO_to_Chlor.inlet.pressure[0]
        == model.fs.RO_to_Chlor.outlet.pressure[0])

    model.fs.RO_to_Chlor.total_flow_cons = Constraint(
        expr=model.fs.RO_to_Chlor.outlet.flow_mol[0] ==
            (model.fs.RO_to_Chlor.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O']/18e-3) +
            (model.fs.RO_to_Chlor.inlet.flow_mass_phase_comp[0, 'Liq', 'TDS']/58.4e-3) )

    model.fs.RO_to_Chlor.H_con = Constraint( expr=model.fs.RO_to_Chlor.outlet.mole_frac_comp[0, "H_+"] == 0 )
    model.fs.RO_to_Chlor.OH_con = Constraint( expr=model.fs.RO_to_Chlor.outlet.mole_frac_comp[0, "OH_-"] == 0 )
    model.fs.RO_to_Chlor.HOCl_con = Constraint( expr=model.fs.RO_to_Chlor.outlet.mole_frac_comp[0, "HOCl"] == 0 )


    model.fs.RO_to_Chlor.OCl_con = Constraint(
        expr=model.fs.RO_to_Chlor.outlet.mole_frac_comp[0, "OCl_-"] == total_chlorine_inlet/total_molar_density )

    model.fs.RO_to_Chlor.Cl_con = Constraint(
        expr=model.fs.RO_to_Chlor.outlet.mole_frac_comp[0, "Cl_-"] ==
            (model.fs.RO_to_Chlor.inlet.flow_mass_phase_comp[0, 'Liq', 'TDS']/58.4e-3) /
             model.fs.RO_to_Chlor.outlet.flow_mol[0] )

    model.fs.RO_to_Chlor.Na_con = Constraint(
        expr=model.fs.RO_to_Chlor.outlet.mole_frac_comp[0, "Na_+"] ==
            (model.fs.RO_to_Chlor.inlet.flow_mass_phase_comp[0, 'Liq', 'TDS']/58.4e-3) /
             model.fs.RO_to_Chlor.outlet.flow_mol[0] + model.fs.RO_to_Chlor.outlet.mole_frac_comp[0, "OCl_-"])

    model.fs.RO_to_Chlor.H2O_con = Constraint(
        expr=model.fs.RO_to_Chlor.outlet.mole_frac_comp[0, "H2O"] == 1 -
            sum(model.fs.RO_to_Chlor.outlet.mole_frac_comp[0, j] for j in ["H_+", "OH_-",
                "HOCl", "OCl_-", "Cl_-", "Na_+"]) )

    # Add the connecting arcs
    model.fs.S1 = Arc(source=desal_port['out'], destination=model.fs.RO_to_Chlor.inlet)
    #model.fs.S1 = Arc(source=model.fs.RO.permeate, destination=model.fs.RO_to_Chlor.inlet)
    model.fs.S2 = Arc(source=model.fs.RO_to_Chlor.outlet, destination=model.fs.simple_naocl_unit.inlet)
    TransformationFactory("network.expand_arcs").apply_to(model)

    # Inlet conditions for RO unit already set from the build function

    # Calculate scaling factors and setup each block
    model.fs.RO_to_Chlor.properties_in[0].mass_frac_phase_comp
    iscale.calculate_scaling_factors(model.fs.RO_to_Chlor)
    model.fs.RO_to_Chlor.initialize(optarg={'nlp_scaling_method': 'user-scaling'})

    iscale.calculate_scaling_factors(model.fs.RO_to_Chlor)
    iscale.constraint_autoscale_large_jac(model.fs.RO_to_Chlor)

    state_args, stoich_extents = approximate_chemical_state_args(model.fs.simple_naocl_unit,
                                model.fs.simple_naocl_rxn_params, simple_naocl_reaction_config)
    calculate_chemical_scaling_factors(model.fs.simple_naocl_unit,
                                model.fs.simple_naocl_thermo_params,
                                model.fs.simple_naocl_rxn_params, state_args)

    # initialize each chemical block (REQUIRED)
    initialize_chlorination_example(model.fs.simple_naocl_unit, state_args)

    # unfix inlet conditions for chlorination
    model.fs.simple_naocl_unit.inlet.pressure.unfix()
    model.fs.simple_naocl_unit.inlet.temperature.unfix()
    model.fs.simple_naocl_unit.inlet.flow_mol.unfix()

    model.fs.simple_naocl_unit.inlet.mole_frac_comp[0, "H_+"].unfix()
    model.fs.simple_naocl_unit.inlet.mole_frac_comp[0, "OH_-"].unfix()
    model.fs.simple_naocl_unit.inlet.mole_frac_comp[0, "HOCl"].unfix()
    model.fs.simple_naocl_unit.inlet.mole_frac_comp[0, "Cl_-"].unfix()

    #   Here is where I would generally add NaOCl (not sure how best to handle it here)
    #       May have to remove translator constraints? or put added chlorine in those
    #       constraints? Or add a mixer block before the chlorination?
    model.fs.simple_naocl_unit.inlet.mole_frac_comp[0, "OCl_-"].unfix()
    model.fs.simple_naocl_unit.dosing_rate.unfix()
    model.fs.simple_naocl_unit.inlet.mole_frac_comp[0, "Na_+"].unfix()
    model.fs.simple_naocl_unit.inlet.mole_frac_comp[0, "H2O"].unfix()

    check_dof(model)

def run_0DRO_Chlorination_flowsheet_example(with_seq_decomp=True, RO_level='detailed'):
    model = ConcreteModel()
    model.fs = FlowsheetBlock(default={"dynamic": False})

    # build the flow sheet
    build_0DRO_Chlorination_flowsheet(model, mg_per_L_NaOCl_added=2, RO_level=RO_level)

    # Call the sequential decomposition initializer tool
    if with_seq_decomp == True:
        seq_decomp_initializer(model)
    # Manually initialize the blocks sequentially
    else:
        propagate_state(model.fs.S1)
        model.fs.RO_to_Chlor.initialize(optarg={'nlp_scaling_method': 'user-scaling'})
        propagate_state(model.fs.S2)

        # Use propogated state (which updated chlorination unit inlet port) to initialize
        model.fs.simple_naocl_unit.initialize(optarg={'nlp_scaling_method': 'user-scaling',
                                     'bound_push': 1e-10,
                                     'mu_init': 1e-6})

        # Recall the auto scaling after initialization
        ##  NOTE: There is a bug in the autoscaling or RO unit which does not allow
        #           us to use autoscaling on that unit if using the 'detailed' version
        #iscale.constraint_autoscale_large_jac(model.fs.RO)
        iscale.constraint_autoscale_large_jac(model.fs.simple_naocl_unit)
        iscale.constraint_autoscale_large_jac(model.fs.RO_to_Chlor)

    #End manual seq decomp

    model.fs.simple_naocl_unit.free_chlorine.fix(2)
    model.fs.RO_to_Chlor.OCl_con.deactivate()

    solve_with_user_scaling(model, tee=True)

    model.fs.RO.inlet.display()
    model.fs.RO.permeate.display()
    model.fs.RO.retentate.display()

    model.fs.RO_to_Chlor.inlet.display()
    model.fs.RO_to_Chlor.outlet.display()

    display_results_of_chlorination(model.fs.simple_naocl_unit)

    return model

def run_0DRO_Chlorination_flowsheet_optimization_example(with_seq_decomp=True, RO_flux=20):

    # First step is to build and solve the flowsheet under current conditions
    model = run_0DRO_Chlorination_flowsheet_example(with_seq_decomp=with_seq_decomp, RO_level='simple')

    model.fs.RO.area.display()
    model.fs.RO.permeate_side.properties_mixed[0].flow_vol.display()

    # Unfix RO area and replace with constraint on RO_flux
    model.fs.RO.area.unfix()
    model.fs.RO.area.setlb(10)
    model.fs.RO.area.setub(300)

    # fixed RO water flux
    model.fs.RO_flux = Expression(
        expr=model.fs.RO.permeate_side.properties_mixed[0].flow_vol
             / model.fs.RO.area)
    model.fs.eq_RO_flux = Constraint(
        expr=model.fs.RO_flux*1000*3600 >= RO_flux)

    # Unfix free chlorine and replace with constraint on minimum chlorine
    model.fs.simple_naocl_unit.free_chlorine.unfix()
    model.fs.simple_naocl_unit.free_chlorine.setlb(0.1)
    model.fs.simple_naocl_unit.free_chlorine.setub(5)

    model.fs.exit_chlorine = Constraint(expr=model.fs.simple_naocl_unit.free_chlorine >= 2)

    # Constraint on flow from permeate to be above certain level
    model.fs.flow_mol_cons = Constraint(expr=model.fs.RO.permeate_side.properties_mixed[0].flow_vol >= 0.0002)

    # Add an objective function (something to minimize)
    model.fs.objective = Objective(expr=model.fs.simple_naocl_unit.dosing_rate*10 + model.fs.RO.area*1)

    iscale.constraint_autoscale_large_jac(model.fs.simple_naocl_unit)
    iscale.constraint_autoscale_large_jac(model.fs.RO_to_Chlor)

    #   Can't use this tool because of issues in detailed RO
    try:
        iscale.constraint_autoscale_large_jac(model)
    except:
        pass

    solve_with_user_scaling(model, tee=True, bound_push=1e-10, mu_init=1e-6)

    model.fs.RO.area.display()
    model.fs.RO.permeate_side.properties_mixed[0].flow_vol.display()
    print(value(model.fs.RO_flux.expr)*1000*3600)

    model.fs.RO.inlet.display()
    model.fs.RO.permeate.display()
    model.fs.RO.retentate.display()

    model.fs.RO_to_Chlor.inlet.display()
    model.fs.RO_to_Chlor.outlet.display()

    display_results_of_chlorination(model.fs.simple_naocl_unit)


    return model

if __name__ == "__main__":
    #model = run_0DRO_Chlorination_flowsheet_example(False)
    #model = run_0DRO_Chlorination_flowsheet_example(True)
    model = run_0DRO_Chlorination_flowsheet_optimization_example()
