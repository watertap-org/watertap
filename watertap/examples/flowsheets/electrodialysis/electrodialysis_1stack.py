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

from xml.etree.ElementInclude import default_loader
from pyomo.environ import (
    ConcreteModel,
    value,
    Param,
    Var,
    Constraint,
    Expression,
    Objective,
    TransformationFactory,
    Block,
    NonNegativeReals,
    RangeSet,
    Set,
    check_optimal_termination,
    assert_optimal_termination,
    units as pyunits,
    Reference
)
from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.solvers import get_solver
from idaes.core.util.initialization import propagate_state
from idaes.core.util.misc import StrEnum
from idaes.core.util.model_statistics import degrees_of_freedom
from sympy import Range
from idaes.models.unit_models import Feed, Product, Separator, Mixer
from idaes.models.unit_models.mixer import MomentumMixingType
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslogger

from watertap.unit_models.electrodialysis_1D import (
    Electrodialysis1D,
    
)
from watertap.unit_models.pressure_changer import Pump, EnergyRecoveryDevice
from watertap.core.util.initialization import (
    assert_no_degrees_of_freedom,
    assert_degrees_of_freedom,
)
from watertap.costing.watertap_costing_package import (
    WaterTAPCosting,
    make_capital_cost_var,
)
from watertap.property_models.ion_DSPMDE_prop_pack import DSPMDEParameterBlock

def main():
    # set up solver
    solver = get_solver()

    # build, set, and initialize
    m = build()
    set_operating_conditions(m)
    initialize_system(m, solver=solver)

    # simulate and display
    solve(m, solver=solver)
    print("\n***---Simulation results---***")
    display_system(m)
    display_design(m)
    display_state(m)

    # optimize and display
    #optimize_set_up(m)
    #optimize(m, solver=solver)
    print("\n***---Optimization results---***")
    display_system(m)
    display_design(m)
    display_state(m)


def build():
    # ---building model---
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    ion_dict = {
            "solute_list": ["Na_+", "Cl_-"],
            "mw_data": {"H2O": 18e-3, "Na_+": 23e-3, "Cl_-": 35.5e-3},
            "electrical_mobility_data": {"Na_+": 5.19e-8, "Cl_-": 7.92e-8},
            "charge": {"Na_+": 1, "Cl_-": -1},
        }
    m.fs.properties = DSPMDEParameterBlock(default=ion_dict)
    m.fs.costing = WaterTAPCosting()
    m.fs.feed = Feed(default={"property_package": m.fs.properties})
    m.fs.Separator_in = Separator(
        default={"property_package": m.fs.properties, "outlet_list": ["F_D", "F_C"]}
    )
    # Add the pumps
    m.fs.Pump_D = Pump(default={"property_package": m.fs.properties})
    m.fs.Pump_C = Pump(default={"property_package": m.fs.properties})
    
    # Add electrodialysis stacks
    m.fs.EDstack = Electrodialysis1D(
        default={
            "property_package": m.fs.properties,
            "operation_mode": "Constant_Voltage",
            "finite_elements": 20,
        },
    )
    m.fs.product_1cp = Product(default={"property_package": m.fs.properties})
    m.fs.disposal_1cp = Product(default={"property_package": m.fs.properties})
    # costing
    m.fs.Pump_D.costing = UnitModelCostingBlock(
        default={"flowsheet_costing_block": m.fs.costing}
    )
    m.fs.Pump_C.costing = UnitModelCostingBlock(
        default={"flowsheet_costing_block": m.fs.costing}
    )
    m.fs.EDstack.costing = UnitModelCostingBlock(
        default={"flowsheet_costing_block": m.fs.costing}
    )
    
    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(m.fs.product_1cp.properties[0].flow_vol * m.fs.EDstack.cell_pair_num)
    m.fs.costing.add_LCOW(m.fs.product_1cp.properties[0].flow_vol * m.fs.EDstack.cell_pair_num)
    m.fs.costing.add_specific_energy_consumption(m.fs.product_1cp.properties[0].flow_vol * m.fs.EDstack.cell_pair_num)


    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.Separator_in.inlet)
    m.fs.s02 = Arc(source=m.fs.Separator_in.F_D, destination=m.fs.Pump_D.inlet)
    m.fs.s03 = Arc(source=m.fs.Separator_in.F_C, destination=m.fs.Pump_C.inlet)
    m.fs.s04 = Arc(source=m.fs.Pump_D.outlet, destination=m.fs.EDstack.inlet_diluate)
    m.fs.s05 = Arc(source=m.fs.Pump_C.outlet, destination=m.fs.EDstack.inlet_concentrate)
    m.fs.s06 = Arc(source=m.fs.EDstack.outlet_diluate, destination=m.fs.product_1cp.inlet)
    m.fs.s07 = Arc(source=m.fs.EDstack.outlet_concentrate, destination=m.fs.disposal_1cp.inlet)
    
    TransformationFactory("network.expand_arcs").apply_to(m)

    #Scaling  
    m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e2, index=("Liq", "H2O")
        )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e4, index=("Liq", "Na_+")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e4, index=("Liq", "Cl_-")
        )
    iscale.set_scaling_factor(m.fs.EDstack.cell_width, 10)
    iscale.set_scaling_factor(m.fs.EDstack.cell_length, 10)
    iscale.calculate_scaling_factors(m)

    return m

def set_operating_conditions(m):

    solver = get_solver()

    # ---specifications---
    # feed
    # state variables
    m.fs.feed.properties[0].pressure.fix(101325)  # feed pressure [Pa]
    m.fs.feed.properties[0].temperature.fix(273.15 + 25)  # feed temperature [K]
    m.fs.Separator_in.F_D.temperature.fix(273.15 + 25)
    '''
    m.fs.Separator_in.F_D.temperature.fix(273.15 + 25)
    m.fs.Separator_in.F_C.temperature.fix(273.15 + 25)
    m.fs.Pump_D.control_volume.properties_out[0].temperature.fix(273.15 + 25)
    m.fs.Pump_C.control_volume.properties_out[0].temperature.fix(273.15 + 25)
    m.fs.product_1cp.properties[0].temperature.fix(273.15 + 25)
    m.fs.disposal_1cp.properties[0].temperature.fix(273.15 + 25)
    '''
    # properties (cannot be fixed for initialization routines, must calculate the state variables)
    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_mol_phase_comp", ("Liq","H2O")): 2.4, 
            ("flow_mol_phase_comp", ("Liq", "Na_+")): 7.38e-3,
            ("flow_mol_phase_comp", ("Liq", "Cl_-")): 7.38e-3,
        },  # feed molar concentration of Na and Cl ions
        hold_state=True,  
    )
    #m.fs.Separator_in.deltaP.fix(0)
    # separator, no degrees of freedom (i.e. equal flow rates in PXR determines split fraction)
    
    

   
    m.fs.Pump_D.efficiency_pump.fix(0.80) 
    m.fs.Pump_C.efficiency_pump.fix(0.80) 
    #m.fs.Pump_D.deltaP.fix(2e4)
    #m.fs.Pump_C.deltaP.fix(2e4)
    m.fs.Pump_D.control_volume.properties_out[0].pressure.fix(121325) #Assumed 0.2 bar pressure drop in the stack for both channels. 
    m.fs.Pump_C.control_volume.properties_out[0].pressure.fix(121325)

    # ED unit
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

    # initialize ED

    m.fs.EDstack.inlet_diluate.flow_mol_phase_comp[0, 'Liq', 'H2O'] = value(
        m.fs.feed.properties[0].flow_mol_phase_comp['Liq', 'H2O']
     / (m.fs.EDstack.cell_pair_num * 2))
    m.fs.EDstack.inlet_diluate.flow_mol_phase_comp[0, 'Liq', 'Na_+'] = value(
        m.fs.feed.properties[0].flow_mol_phase_comp['Liq', 'Na_+']
      / (m.fs.EDstack.cell_pair_num * 2))
    m.fs.EDstack.inlet_diluate.flow_mol_phase_comp[0, 'Liq', 'Cl_-'] = value(
        m.fs.feed.properties[0].flow_mol_phase_comp['Liq', 'Cl_-']
     / (m.fs.EDstack.cell_pair_num * 2))

    m.fs.EDstack.inlet_concentrate.flow_mol_phase_comp[0, 'Liq', 'H2O'] = value(
        m.fs.feed.properties[0].flow_mol_phase_comp['Liq', 'H2O']
      / (m.fs.EDstack.cell_pair_num * 2))
    m.fs.EDstack.inlet_concentrate.flow_mol_phase_comp[0, 'Liq', 'Na_+'] = value(
        m.fs.feed.properties[0].flow_mol_phase_comp['Liq', 'Na_+']
      / (m.fs.EDstack.cell_pair_num * 2))
    m.fs.EDstack.inlet_concentrate.flow_mol_phase_comp[0, 'Liq', 'Cl_-'] = value(
        m.fs.feed.properties[0].flow_mol_phase_comp['Liq', 'Cl_-']
      / (m.fs.EDstack.cell_pair_num * 2))

    m.fs.EDstack.inlet_diluate.temperature = value(
        m.fs.feed.properties[0].temperature
    )
    m.fs.EDstack.inlet_concentrate.temperature = value(
        m.fs.feed.properties[0].temperature
    )
    
    m.fs.EDstack.initialize(optarg=solver.options)

    # check degrees of freedom
    
    if degrees_of_freedom(m) != 0:
        raise RuntimeError(
            "The set_operating_conditions function resulted in {} "
            "degrees of freedom rather than 0. This error suggests "
            "that too many or not enough variables are fixed for a "
            "simulation.".format(degrees_of_freedom(m))
        )

def solve(blk, solver=None, tee=False, check_termination=True):
    if solver is None:
        solver = get_solver()
    results = solver.solve(blk, tee=tee)
    if check_termination:
        assert_optimal_termination(results)
    return results


def initialize_system(m, solver=None):
    if solver is None:
        solver = get_solver()
    optarg = solver.options

    # ---initialize RO---
 

    # ---initialize feed block---
    m.fs.feed.initialize(optarg=optarg)
    propagate_state(m.fs.s01)
    solve(m.fs.Separator_in)
    propagate_state(m.fs.s02)
    solve(m.fs.Pump_D)
    propagate_state(m.fs.s03)
    solve(m.fs.Pump_C)
    propagate_state(m.fs.s04)
    propagate_state(m.fs.s05)
    propagate_state(m.fs.s06)
    propagate_state(m.fs.s07)
    m.fs.EDstack.initialize(optarg=optarg)
    m.fs.costing.initialize()

'''
def optimize_set_up(m):
    # objective
    m.fs.objective = Objective(expr=m.fs.costing.LCOW)

    # unfix decision variables and add bounds
    # pump 1 and pump 2
    m.fs.P1.control_volume.properties_out[0].pressure.unfix()
    m.fs.P1.control_volume.properties_out[0].pressure.setlb(10e5)
    m.fs.P1.control_volume.properties_out[0].pressure.setub(80e5)
    m.fs.P1.deltaP.setlb(0)
    m.fs.P2.control_volume.properties_out[0].pressure.setlb(10e5)
    m.fs.P2.control_volume.properties_out[0].pressure.setub(80e5)
    m.fs.P2.deltaP.setlb(0)

    # RO
    m.fs.RO.area.setlb(1)
    m.fs.RO.area.setub(150)

    # additional specifications
    m.fs.product_salinity = Param(
        initialize=500e-6, mutable=True
    )  # product NaCl mass fraction [-]
    m.fs.minimum_water_flux = Param(
        initialize=1.0 / 3600.0, mutable=True
    )  # minimum water flux [kg/m2-s]

    # additional constraints
    m.fs.eq_product_quality = Constraint(
        expr=m.fs.product.properties[0].mass_frac_phase_comp["Liq", "NaCl"]
        <= m.fs.product_salinity
    )
    iscale.constraint_scaling_transform(
        m.fs.eq_product_quality, 1e3
    )  # scaling constraint
    m.fs.eq_minimum_water_flux = Constraint(
        expr=m.fs.RO.flux_mass_phase_comp[0, 1, "Liq", "H2O"] >= m.fs.minimum_water_flux
    )

    # ---checking model---
    assert_degrees_of_freedom(m, 1)


def optimize(m, solver=None, check_termination=True):
    # --solve---
    return solve(m, solver=solver, check_termination=check_termination)

'''

def display_system(m):
    print("---system metrics---")
    feed_ion_mass_flow = (
        sum(
            m.fs.feed.properties[0].flow_mass_phase_comp["Liq", j] for j in m.fs.properties.ion_set 
        )
    ) # in [kg/s]
        
    feed_ion_mass_conc = (
        sum(
            m.fs.feed.properties[0].conc_mass_phase_comp["Liq", j] for j in m.fs.properties.ion_set 
        )
    ) # in [kg/m3]
    print(f'Feed: total ion mass flow = {value(feed_ion_mass_flow)} kg/s, total ion mass concentration = {value(feed_ion_mass_conc)} kg/m3')

    product_ion_mass_flow = (
        sum(
            m.fs.EDstack.cell_pair_num * m.fs.product_1cp.properties[0].flow_mass_phase_comp["Liq", j] for j in m.fs.properties.ion_set 
        )
    ) # in [kg/s]
        
    product_ion_mass_conc = (
        sum(
            m.fs.product_1cp.properties[0].conc_mass_phase_comp["Liq", j] for j in m.fs.properties.ion_set 
        )
    ) # in [kg/m3]
    print(f'Produced water: total ion mass flow = {value(product_ion_mass_flow)} kg/s, total ion mass concentration = {value(product_ion_mass_conc)} kg/m3')
    print(f'Water recovery by mass: {value(m.fs.EDstack.water_recovery_mass[0])}')
    print(f'Specific energy consumption: {value(m.fs.costing.specific_energy_consumption)} kWh/m3')
    print(f'Levelized cost of water: {value(m.fs.costing.LCOW)} $/m3')

def display_design(m):
    print("---decision variables---")
    print(f'Electrical Voltage on an elecctrodialysis stack: {value(m.fs.EDstack.voltage_applied[0])} volt')
    print(f'Total membrane area: {value(m.fs.EDstack.cell_width * m.fs.EDstack.cell_length * m.fs.EDstack.cell_pair_num)} m2')

    print("---design variables---")
    print("Separator")
    print(f"Diluate channel feed ratio: {value(m.fs.Separator_in.split_fraction[0,'F_D'])}")

    print(
        "Pump Diluate \noutlet pressure: %.1f bar\npower %.2f kW"
        % (
            m.fs.Pump_D.outlet.pressure[0].value / 1e5,
            m.fs.Pump_D.work_mechanical[0].value / 1e3,
        )
    )
    print(
        "Pump Concentrate\noutlet pressure: %.1f bar\npower %.2f kW"
        % (
            m.fs.Pump_C.outlet.pressure[0].value / 1e5,
            m.fs.Pump_C.work_mechanical[0].value / 1e3,
        )
    )


def display_state(m):
    print("---state---")

    def print_state(s, b):
        if b.parent_block() in [m.fs.EDstack.inlet_diluate, m.fs.EDstack.inlet_concentrate, m.fs.EDstack.outlet_diluate, m.fs.EDstack.outlet_concentrate]:
            Ion_mass_flow = sum(
            b.flow_mass_phase_comp[0, "Liq", j] for j in m.fs.properties.ion_set 
            ) * m.fs.EDstack.cell_pair_num
        else:
            Ion_mass_flow = sum(
            b.flow_mass_phase_comp[0, "Liq", j] for j in m.fs.properties.ion_set 
            )
        
        Ion_mass_conc = sum(
            b.conc_mass_phase_comp[ 0, "Liq", j] for j in m.fs.properties.ion_set 
            )
        
        print(f'{s}: ion mass flow rate = {value(Ion_mass_flow)} kg/s; ion mass concentration = {value(Ion_mass_conc)} kg/m3')
            

    print_state("Feed      ", m.fs.feed.outlet)
    print_state("Separator to Diluate ", m.fs.Separator_in.F_D)
    print_state("Separator to Concentrate ", m.fs.Separator_in.F_C)
    print_state("Pump to Diluate  ", m.fs.Pump_D.outlet)
    print_state("Pump to Concentrate  ", m.fs.Pump_C.outlet)
    print_state("ED produced water   ", m.fs.EDstack.outlet_diluate)
    print_state("ED disposal water  ", m.fs.EDstack.outlet_concentrate)
    


if __name__ == "__main__":
    main()