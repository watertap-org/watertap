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

from pyomo.environ import (
    ConcreteModel,
    Var,
    value,
    Constraint,
    Objective,
    TransformationFactory,
    assert_optimal_termination,
    units as pyunits,
    NonNegativeReals,
)
from pyomo.network import Arc

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.solvers import get_solver
from idaes.core.util.initialization import propagate_state
#from idaes.core.util.model_statistics import degrees_of_freedom, report_statistics, variables_set
import idaes.core.util.model_statistics as mstat
from pyomo.core.expr.current import identify_variables
from idaes.models.unit_models import Feed, Product, Separator, Mixer
from idaes.models.unit_models.mixer import MixingType, MomentumMixingType
from pandas import DataFrame
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslogger
from pytest import approx
from watertap.core.util.initialization import check_dof
from watertap.unit_models.electrodialysis_1D import ElectricalOperationMode
from watertap.unit_models.electrodialysis_1D import Electrodialysis1D
from watertap.costing.watertap_costing_package import (
    MixerType,
    WaterTAPCosting,
)
from watertap.property_models.ion_DSPMDE_prop_pack import DSPMDEParameterBlock

__author__ = "Xiangyu Bi"


def main():
    # set up solver
    solver = get_solver()
    # Simulate a fully defined operation
    m = build()
    set_operating_conditions(m)
    initialize_system(m, solver=solver)
    solve(m, solver=solver)
    '''
    # Temperature and Pressure var diagnositcs 
    count1=0
    count2=0
    count3=0
    ct1=0
    ct2=0
    ct3=0
    
    for k in mstat.variables_set(m):
        if "temperature" in k.name:
            print("ALL T:",type(k), k.name)
            ct1+=1
    print("NUM of ALL T", ct1)

    for k in mstat.fixed_variables_generator(m):
        if "temperature" in k.name:
            print("FIXED T:", k)
            ct3+=1
    print("NUM of FIXED T", ct3)       
    for k in mstat.activated_equalities_generator(m):
        for l in identify_variables(k.body):
            if "temperature" in l.name and "pressure" not in k.name and"dens" not in k.name:
                print("CONSTR of T:", k)
                ct2+=1
                break
    print("NUM of Constr T", ct2)
            #for j in mstat.activated_constraints_generator(m):
             ##      if k.name == l:
               #         print("Constraints of T:", j.body)

    

    for k in mstat.variables_set(m):
        if "pressure" in k.name and "pressure_" not in k.name and "_pressure" not in k.name:
            print("ALL P:", k.name)
            count1+=1
    print("NUM of ALL P", count1)

    for k in mstat.fixed_variables_generator(m):
        if "pressure" in k.name and "pressure_" not in k.name and "_pressure" not in k.name:
            print("FIXED P:", k.name)
            count2+=1
    print("NUM of fixed P", count2)

    for k in mstat.activated_equalities_generator(m):
        for l in identify_variables(k.body):
            if "pressure" in l.name and "pressure_" not in l.name and "_pressure" not in l.name:
                print("CONSTR of P:", k)
                count3+=1
                break
    print("NUM of Constr P", count3)
    '''
    print("\n***---Simulation results---***")
    display_model_metrics(m)

    # Perform an optimization over selected variables
    #initialize_system(m, solver=solver)
    #optimize_system(m, solver=solver)
    #print("\n***---Optimization results---***")
    #display_model_metrics(m)


def build():
    # ---building model---
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    ion_dict = {
        "solute_list": ["Na_+", "Cl_-"],
        "mw_data": {"H2O": 18e-3, "Na_+": 23e-3, "Cl_-": 35.5e-3},
        "elec_mobility_data": {("Liq", "Na_+"): 5.19e-8, ("Liq", "Cl_-"): 7.92e-8},
        "charge": {"Na_+": 1, "Cl_-": -1},
    }
    m.fs.properties = DSPMDEParameterBlock(**ion_dict)
    m.fs.costing = WaterTAPCosting()
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.sepa0 = Separator(
        property_package=m.fs.properties,
        outlet_list=["to_dil", "to_conc_in0"],
    ) 
    m.fs.mix0 = Mixer(property_package=m.fs.properties, energy_mixing_type = MixingType.none, momentum_mixing_type = MomentumMixingType.none,
        inlet_list=["from_feed", "from_conc_out"])

    # Add electrodialysis (ED) stacks
    m.fs.EDstack = Electrodialysis1D(
        property_package=m.fs.properties,
        operation_mode=ElectricalOperationMode.Constant_Voltage,
        finite_elements=20,
    )
    m.fs.sepa1 = Separator(
        property_package=m.fs.properties,
        outlet_list=["to_disp", "to_conc_in1"],
    ) 
    m.fs.prod = Product(property_package=m.fs.properties)
    m.fs.disp = Product(property_package=m.fs.properties)

    # Touching needed variables for initialization and displaying results
    m.fs.feed.properties[0].conc_mass_phase_comp[...]
    m.fs.sepa0.to_dil_state[0].conc_mass_phase_comp[...]
    m.fs.sepa0.to_conc_in0_state[0].conc_mass_phase_comp[...]
    m.fs.sepa1.to_disp_state[0].conc_mass_phase_comp[...]
    m.fs.sepa1.to_conc_in1_state[0].conc_mass_phase_comp[...]
    m.fs.mix0.from_feed_state[0].conc_mass_phase_comp[...]
    m.fs.mix0.from_conc_out_state[0].conc_mass_phase_comp[...]
    m.fs.prod.properties[0].conc_mass_phase_comp[...]
    m.fs.disp.properties[0].conc_mass_phase_comp[...]

    m.fs.feed.properties[0].flow_vol_phase[...]
    m.fs.sepa0.to_dil_state[0].flow_vol_phase[...]
    m.fs.sepa0.to_conc_in0_state[0].flow_vol_phase[...]
    m.fs.sepa1.to_disp_state[0].flow_vol_phase[...]
    m.fs.sepa1.to_conc_in1_state[0].flow_vol_phase[...]
    m.fs.mix0.from_feed_state[0].flow_vol_phase[...]
    m.fs.mix0.from_conc_out_state[0].flow_vol_phase[...]
    m.fs.prod.properties[0].flow_vol_phase[...]
    m.fs.disp.properties[0].flow_vol_phase[...]
    m.fs.EDstack.diluate.properties[...].flow_vol_phase[...]
    m.fs.EDstack.concentrate.properties[...].flow_vol_phase[...]


    #m.fs.EDstack.inlet_diluate.flow_vol_phase[...]
    #m.fs.EDstack.inlet_concentrate.flow_vol_phase[...]

    # costing
    m.fs.EDstack.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(
        m.fs.prod.properties[0].flow_vol_phase["Liq"]
    )
    m.fs.costing.add_LCOW(m.fs.prod.properties[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(
        m.fs.prod.properties[0].flow_vol_phase["Liq"]
    )

    # add extra variables and constraints
    m.fs.recovery_volume_H2O = Var(
            initialize=0.5,
            bounds=(0, 1),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="flowsheet level water recovery calculated by volumeric flow rate",
        )
    m.fs.mem_area = Var(
        initialize=1,
        bounds=(0, 1e3),
        units=pyunits.meter**2,
        doc="Total membrane area for cem (or aem) in one stack",
    )
    m.fs.prod_salinity = Var(
        initialize=1, bounds=(0, 1000), units=pyunits.kg * pyunits.meter**-3
    )
    m.fs.disp_salinity = Var(
        initialize=1, bounds=(0, 1e6), units=pyunits.kg * pyunits.meter**-3
    )
    m.fs.eq_recovery_volume_H2O = Constraint(
        expr=m.fs.recovery_volume_H2O == 
        m.fs.prod.properties[0].flow_vol_phase['Liq']
        * m.fs.feed.properties[0].flow_vol_phase['Liq']**-1
    )
    m.fs.eq_electrodialysis_equal_flow =  Constraint(
        expr=m.fs.EDstack.diluate.properties[0,0].flow_vol_phase["Liq"]==m.fs.EDstack.concentrate.properties[0,0].flow_vol_phase["Liq"]
    )
    
    m.fs.eq_product_salinity = Constraint(
        expr=m.fs.prod_salinity
        == sum(
            m.fs.prod.properties[0].conc_mass_phase_comp["Liq", j]
            for j in m.fs.properties.ion_set | m.fs.properties.solute_set
        )
    )
    m.fs.eq_disposal_salinity = Constraint(
        expr=m.fs.disp_salinity
        == sum(
            m.fs.disp.properties[0].conc_mass_phase_comp["Liq", j]
            for j in m.fs.properties.ion_set | m.fs.properties.solute_set
        )
    )

    m.fs.eq_mem_area = Constraint(
        expr=m.fs.mem_area
        == m.fs.EDstack.cell_width
        * m.fs.EDstack.cell_length
        * m.fs.EDstack.cell_pair_num
    )

    # Add Arcs
    m.fs.arc0 = Arc(source=m.fs.feed.outlet, destination=m.fs.sepa0.inlet)
    m.fs.arc1 = Arc(
        source=m.fs.sepa0.to_dil, destination=m.fs.EDstack.inlet_diluate
    )
    m.fs.arc2 = Arc(
        source=m.fs.sepa0.to_conc_in0,
        destination=m.fs.mix0.from_feed,
    )
    m.fs.arc3 = Arc(source=m.fs.mix0.outlet, destination=m.fs.EDstack.inlet_concentrate)
    m.fs.arc4 = Arc(
        source=m.fs.EDstack.outlet_diluate, destination=m.fs.prod.inlet
    )
    m.fs.arc5 = Arc(
        source=m.fs.EDstack.outlet_concentrate, destination=m.fs.sepa1.inlet
    )
    m.fs.arc6 = Arc(
        source=m.fs.sepa1.to_disp, destination=m.fs.disp.inlet
    )
    m.fs.arc7 = Arc(
        source=m.fs.sepa1.to_conc_in1, destination=m.fs.mix0.from_conc_out
    )
    TransformationFactory("network.expand_arcs").apply_to(m)

    # Scaling
    m.fs.properties.set_default_scaling("flow_mol_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e1, index=("Liq", "Na_+")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e1, index=("Liq", "Cl_-")
    )
    iscale.set_scaling_factor(m.fs.EDstack.cell_width, 10)
    iscale.set_scaling_factor(m.fs.EDstack.cell_length, 10)
    iscale.calculate_scaling_factors(m)
    return m


def set_operating_conditions(m):

    solver = get_solver()

    # ---specifications---
    # Fix state variables at the origin
    m.fs.feed.properties[0].pressure.fix(101325)  # feed pressure [Pa]
    m.fs.feed.properties[0].temperature.fix(298.15)  # feed temperature [K]
    m.fs.EDstack.concentrate.properties[0,0].temperature.fix(298.15)
    m.fs.EDstack.concentrate.properties[0,0].pressure.fix(101325)
    m.fs.feed.properties[0].flow_mol_phase_comp["Liq", "H2O"].fix(4.8)
    m.fs.feed.properties[0].flow_mol_phase_comp["Liq", "Na_+"].fix(1.476e-2)
    m.fs.feed.properties[0].flow_mol_phase_comp["Liq", "Cl_-"].fix(1.476e-2)
    m.fs.recovery_volume_H2O.fix(0.5)
    # Fix separator's split_fraction to 0.5, i.e., equal flows into the diluate and concentrate channels
    #m.fs.sepa0.split_fraction[0, "inlet_diluate"].fix(0.5)


    # Fix ED unit vars
    m.fs.EDstack.water_trans_number_membrane["cem"].fix(5.8)
    m.fs.EDstack.water_trans_number_membrane["aem"].fix(4.3)
    m.fs.EDstack.water_permeability_membrane["cem"].fix(2.16e-14)
    m.fs.EDstack.water_permeability_membrane["aem"].fix(1.75e-14)
    m.fs.EDstack.voltage_applied.fix(5)
    m.fs.EDstack.electrodes_resistance.fix(0)
    m.fs.EDstack.cell_pair_num.fix(100)
    m.fs.EDstack.current_utilization.fix(1)
    m.fs.EDstack.spacer_thickness.fix(2.7e-4)
    m.fs.EDstack.membrane_areal_resistance["cem"].fix(1.89e-4)
    m.fs.EDstack.membrane_areal_resistance["aem"].fix(1.77e-4)
    m.fs.EDstack.cell_width.fix(0.1)
    m.fs.EDstack.cell_length.fix(0.79)
    m.fs.EDstack.membrane_thickness["aem"].fix(1.3e-4)
    m.fs.EDstack.membrane_thickness["cem"].fix(1.3e-4)
    m.fs.EDstack.solute_diffusivity_membrane["cem", "Na_+"].fix(1.8e-10)
    m.fs.EDstack.solute_diffusivity_membrane["aem", "Na_+"].fix(1.25e-10)
    m.fs.EDstack.solute_diffusivity_membrane["cem", "Cl_-"].fix(1.8e-10)
    m.fs.EDstack.solute_diffusivity_membrane["aem", "Cl_-"].fix(1.25e-10)
    m.fs.EDstack.ion_trans_number_membrane["cem", "Na_+"].fix(1)
    m.fs.EDstack.ion_trans_number_membrane["aem", "Na_+"].fix(0)
    m.fs.EDstack.ion_trans_number_membrane["cem", "Cl_-"].fix(0)
    m.fs.EDstack.ion_trans_number_membrane["aem", "Cl_-"].fix(1)

    # check zero degrees of freedom
    mstat.report_statistics(m)
    check_dof(m)
    



def solve(blk, solver=None, tee=False, check_termination=True):
    if solver is None:
        solver = get_solver()
    results = solver.solve(blk, tee=True)
    if check_termination:
        assert_optimal_termination(results)
    return results


def initialize_system(m, solver=None):

    # set up solver
    if solver is None:
        solver = get_solver()
    optarg = solver.options

    # populate intitial properties throughout the system
    m.fs.feed.initialize(optarg=optarg)
    propagate_state(m.fs.arc0)
    m.fs.sepa0.initialize(optarg=optarg)
    propagate_state(m.fs.arc1)
    propagate_state(m.fs.arc2)
    propagate_state(destination=m.fs.EDstack.inlet_concentrate,source=m.fs.EDstack.inlet_diluate)
    m.fs.EDstack.initialize(optarg=optarg)
    propagate_state(m.fs.arc3,direction="backward")
    propagate_state(m.fs.arc4)
    propagate_state(m.fs.arc5)
    m.fs.prod.initialize()
    m.fs.sepa1.split_fraction[0,"to_conc_in1"].set_value(max(2-value(m.fs.recovery_volume_H2O)**-1, 1e-8))#2-value(m.fs.recovery_volume_H2O)**-1
    m.fs.sepa1.initialize(optarg=optarg)
    propagate_state(m.fs.arc6)
    propagate_state(m.fs.arc7)
    m.fs.mix0.initialize(optarg=optarg)
    m.fs.disp.initialize()
    m.fs.costing.initialize()

'''
def optimize_system(m, solver=None):

    # Below is an example of optimizing the operational voltage and cell pair number (which translates to membrane use)
    # Define a system with zero dof
    set_operating_conditions(m)

    # Set an objective function
    m.fs.objective = Objective(expr=m.fs.costing.LCOW)

    # Choose and unfix variables to be optimized
    m.fs.EDstack.voltage_applied[0].unfix()
    m.fs.EDstack.cell_pair_num.unfix()
    m.fs.EDstack.cell_pair_num.set_value(10)
    # Give narrower bounds to optimizing variables if available
    m.fs.EDstack.voltage_applied[0].setlb(0.01)
    m.fs.EDstack.voltage_applied[0].setub(20)
    m.fs.EDstack.cell_pair_num.setlb(1)
    m.fs.EDstack.cell_pair_num.setub(500)

    # Set a treatment goal
    # Example here is to reach a final product water containing NaCl = 1 g/L (from a 10 g/L feed)
    m.fs.prod.properties[0].conc_mass_phase_comp["Liq", "Na_+"].fix(0.393)

    print("---report model statistics---\n ", report_statistics(m.fs))
    if solver is None:
        solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)

'''
def display_model_metrics(m):

    print("---Flow properties in feed, product and disposal---")

    fp = {
        "Feed": [
            value(m.fs.feed.properties[0].flow_vol_phase["Liq"]),
            value(
                sum(
                    m.fs.feed.properties[0].conc_mass_phase_comp["Liq", j]
                    for j in m.fs.properties.ion_set
                )
            ),
        ],
        "Product": [
            value(m.fs.prod.properties[0].flow_vol_phase["Liq"]),
            value(m.fs.prod_salinity),
        ],
        "Disposal": [
            value(m.fs.disp.properties[0].flow_vol_phase["Liq"]),
            value(m.fs.disp_salinity),
        ],
    }
    fp_table = DataFrame(
        data=fp,
        index=["Volume Flow Rate (m3/s)", "Total Ion Mass Concentration (kg/m3)"],
    )
    print(fp_table)

    print("---Performance Metrics---")

    pm_table = DataFrame(
        data=[
            value(m.fs.EDstack.recovery_mass_H2O[0]),
            value(m.fs.mem_area),
            value(m.fs.EDstack.voltage_applied[0]),
            value(m.fs.costing.specific_energy_consumption),
            value(m.fs.costing.LCOW),
        ],
        columns=["value"],
        index=[
            "Water recovery by mass",
            "Total membrane area (aem or cem), m2",
            "Operation Voltage, V",
            "Specific energy consumption, kWh/m3",
            "Levelized cost of water, $/m3",
        ],
    )
    print(pm_table)


if __name__ == "__main__":
    main()
