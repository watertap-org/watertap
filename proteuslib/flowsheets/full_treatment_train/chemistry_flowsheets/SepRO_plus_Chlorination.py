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

    inlet ---> [RO Separator] ---> permeate ---> (Translator) ---> [Chlorination] ---> outlet
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

    ---------- WORK IN PROGRESS ---------
"""

from proteuslib.flowsheets.full_treatment_train.example_models.unit_separator import (
    build_RO_separator_example, run_RO_example)
from proteuslib.flowsheets.full_treatment_train.chemistry_flowsheets.PostTreatment_SimpleNaOCl_Chlorination import (
    build_simple_naocl_chlorination_unit,
    initialize_chlorination_example,
    display_results_of_chlorination, run_chlorination_example)

# Import specific pyomo objects
from pyomo.environ import (ConcreteModel,
                           SolverStatus,
                           TerminationCondition,
                           Constraint,
                           TransformationFactory,
                           value,
                           Suffix)

# Import the core idaes objects for Flowsheets and types of balances
from idaes.core import FlowsheetBlock

from idaes.generic_models.unit_models.translator import Translator
from pyomo.network import Arc

from proteuslib.flowsheets.full_treatment_train.util import solve_with_user_scaling, check_dof

def build_SepRO_Chlorination_flowsheet(model):
    build_RO_separator_example(model)
    build_simple_naocl_chlorination_unit(model)

    # Translator inlet from RO and outlet goes to chlorination
    model.fs.RO_to_Chlor = Translator(
        default={"inlet_property_package": model.fs.RO_properties,
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
    model.fs.RO_to_Chlor.OCl_con = Constraint( expr=model.fs.RO_to_Chlor.outlet.mole_frac_comp[0, "OCl_-"] == 0 )

    model.fs.RO_to_Chlor.Cl_con = Constraint(
        expr=model.fs.RO_to_Chlor.outlet.mole_frac_comp[0, "Cl_-"] ==
            (model.fs.RO_to_Chlor.inlet.flow_mass_phase_comp[0, 'Liq', 'TDS']/58.4e-3) /
             model.fs.RO_to_Chlor.outlet.flow_mol[0] )

    model.fs.RO_to_Chlor.Na_con = Constraint(
        expr=model.fs.RO_to_Chlor.outlet.mole_frac_comp[0, "Na_+"] ==
            (model.fs.RO_to_Chlor.inlet.flow_mass_phase_comp[0, 'Liq', 'TDS']/58.4e-3) /
             model.fs.RO_to_Chlor.outlet.flow_mol[0] )

    model.fs.RO_to_Chlor.H2O_con = Constraint(
        expr=model.fs.RO_to_Chlor.outlet.mole_frac_comp[0, "H2O"] == 1 -
            sum(model.fs.RO_to_Chlor.outlet.mole_frac_comp[0, j] for j in ["H_+", "OH_-",
                "HOCl", "OCl_-", "Cl_-", "Na_+"]) )


    #check_dof(model.fs.RO_to_Chlor)
    model.fs.RO_to_Chlor.pprint()



def run_SepRO_Chlorination_flowsheet_example():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    build_SepRO_Chlorination_flowsheet(m)

    return m

if __name__ == "__main__":
    model = run_SepRO_Chlorination_flowsheet_example()
