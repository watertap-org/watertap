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

import itertools

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


def run_electrodialysis(
    number_of_stacks,
    #water_recovery=None,
    Conc_in=None,
    #Cbrine=None,
    #A_case=ACase.fixed,
    #B_case=BCase.optimize,
    #AB_tradeoff=ABTradeoff.none,
    #A_value=None,
    #has_NaCl_solubility_limit=None,
    #has_calculated_concentration_polarization=None,
    #has_calculated_ro_pressure_drop=None,
    #permeate_quality_limit=None,
    #AB_gamma_factor=None,
    #B_max=None,
    finite_elements=20,
):
    m = build(
        number_of_stacks,
        #has_NaCl_solubility_limit,
        #has_calculated_concentration_polarization,
        #has_calculated_ro_pressure_drop,
        finite_elements,
        #B_max,
    )
    set_operating_conditions(m, Conc_in)

    initialize(m)
    solve(m)
    print("\n***---Simulation results---***")
    display_system(m)
    display_design(m)
    display_state(m)

    optimize_set_up(
        m,
        water_recovery,
        Cbrine,
        A_case,
        B_case,
        AB_tradeoff,
        A_value,
        permeate_quality_limit,
        AB_gamma_factor,
        B_max,
    )
    res = solve(m, raise_on_failure=False, tee=False)
    print("\n***---Optimization results---***")
    if check_optimal_termination(res):
        display_system(m)
        display_design(m)
        display_state(m)
        display_RO_reports(m)

    return m, res


def build():
    # ---building model---
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = DSPMDEParameterBlock()
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
    m.fs.product = Product(default={"property_package": m.fs.properties})
    m.fs.disposal = Product(default={"property_package": m.fs.properties})
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
    m.fs.costing.add_annual_water_production(m.fs.product.properties[0].flow_vol)
    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(m.fs.product.properties[0].flow_vol)


    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.Separator_in.inlet)
    m.fs.s02 = Arc(source=m.fs.Separator_in.F_D, destination=m.fs.Pump_D.inlet)
    m.fs.s03 = Arc(source=m.fs.Separator_in.F_C, destination=m.fs.Pump_C.inlet)
    m.fs.s04 = Arc(source=m.fs.Pump_D.outlet, destination=m.fs.EDstack.inlet_diluate)
    m.fs.s05 = Arc(source=m.fs.Pump_C.outlet, destination=m.fs.EDstack.inlet_concentrate)
    
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

def set_operating_conditions(m, water_recovery=0.5):
    if solver is None:
        solver = get_solver()

    # ---specifications---
    # feed
    # state variables
    m.fs.feed.properties[0].pressure.fix(101325)  # feed pressure [Pa]
    m.fs.feed.properties[0].temperature.fix(273.15 + 25)  # feed temperature [K]
    # properties (cannot be fixed for initialization routines, must calculate the state variables)
    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): 8.64e-5,  # feed volumetric flow rate [m3/s] assuming 100 cp for now
            ("conc_mol_phase_comp", ("Liq", "Na_+")): 171,
            ("conc_mol_phase_comp", ("Liq", "Cl_-")): 171,
        },  # feed molar concentration of Na and Cl ions
        hold_state=True,  
    )

    # separator, no degrees of freedom (i.e. equal flow rates in PXR determines split fraction)

   
    m.fs.Pump_D.efficiency_pump.fix(0.80)  # pump efficiency [-]
    m.fs.Pump_C.efficiency_pump.fix(0.80)

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
    ) / (m.fs.EDstack.cell_pair_num * 2)
    m.fs.EDstack.inlet_diluate.flow_mol_phase_comp[0, 'Liq', 'Na_+'] = value(
        m.fs.feed.properties[0].flow_mol_phase_comp['Liq', 'Na_+']
    )  / (m.fs.EDstack.cell_pair_num * 2)
    m.fs.EDstack.inlet_diluate.flow_mol_phase_comp[0, 'Liq', 'Cl_-'] = value(
        m.fs.feed.properties[0].flow_mol_phase_comp['Liq', 'Cl_-']
    )  / (m.fs.EDstack.cell_pair_num * 2)

    m.fs.EDstack.inlet_concentrate.flow_mol_phase_comp[0, 'Liq', 'H2O'] = value(
        m.fs.feed.properties[0].flow_mol_phase_comp['Liq', 'H2O']
    )  / (m.fs.EDstack.cell_pair_num * 2)
    m.fs.EDstack.inlet_concentrate.flow_mol_phase_comp[0, 'Liq', 'Na_+'] = value(
        m.fs.feed.properties[0].flow_mol_phase_comp['Liq', 'Na_+']
    )  / (m.fs.EDstack.cell_pair_num * 2)
    m.fs.EDstack.inlet_concentrate.flow_mol_phase_comp[0, 'Liq', 'Cl_-'] = value(
        m.fs.feed.properties[0].flow_mol_phase_comp['Liq', 'Cl_-']
    )  / (m.fs.EDstack.cell_pair_num * 2)

    m.fs.EDstack.inlet_diluate.temperature = value(
        m.fs.feed.properties[0].temperature
    )
    m.fs.EDstack.inlet_concentrate.temperature = value(
        m.fs.feed.properties[0].temperature
    )
    m.fs.EDstack.inlet_diluate.pressure = value(
        m.fs.feed.properties[0].pressure
    )
    m.fs.EDstack.inlet_concentrate.pressure = value(
        m.fs.feed.properties[0].pressure
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
            m.fs.feed.flow_mass_phase_comp[0, "Liq", "j"] for j in m.fs.properties.ion_set 
        )
    ) # in [kg/s]
        
    feed_ion_mass_conc = (
        sum(
            m.fs.feed.conc_mass_phase_comp[0, "Liq", "j"] for j in m.fs.properties.ion_set 
        )
    ) # in [kg/m3]
    print(f'Feed: flow velocity = {feed_ion_mass_flow} kg/s, total ion mass concentration = {feed_ion_mass_conc} kg/m3')

    product_ion_mass_flow = (
        sum(
            m.fs.product.flow_mass_phase_comp[0, "Liq", "j"] for j in m.fs.properties.ion_set 
        )
    ) # in [kg/s]
        
    product_ion_mass_conc = (
        sum(
            m.fs.product.conc_mass_phase_comp[0, "Liq", "j"] for j in m.fs.properties.ion_set 
        )
    ) # in [kg/m3]
    print(f'Produced water: flow velocity = {feed_ion_mass_flow} kg/s, total ion mass concentration = {feed_ion_mass_conc} kg/m3')

    

    print(
        "Volumetric recovery: %.1f%%"
        % (value(m.fs.RO.recovery_vol_phase[0, "Liq"]) * 100)
    )
    print(
        "Water recovery: %.1f%%"
        % (value(m.fs.RO.recovery_mass_phase_comp[0, "Liq", "H2O"]) * 100)
    )
    print(
        "Energy Consumption: %.1f kWh/m3"
        % value(m.fs.costing.specific_energy_consumption)
    )
    print("Levelized cost of water: %.2f $/m3" % value(m.fs.costing.LCOW))


def display_design(m):
    print("---decision variables---")
    print("Operating pressure %.1f bar" % (m.fs.RO.inlet.pressure[0].value / 1e5))
    print("Membrane area %.1f m2" % (m.fs.RO.area.value))

    print("---design variables---")
    print("Separator")
    print("Split fraction %.2f" % (m.fs.S1.split_fraction[0, "PXR"].value * 100))
    print(
        "Pump 1\noutlet pressure: %.1f bar\npower %.2f kW"
        % (
            m.fs.P1.outlet.pressure[0].value / 1e5,
            m.fs.P1.work_mechanical[0].value / 1e3,
        )
    )
    print(
        "Pump 2\noutlet pressure: %.1f bar\npower %.2f kW"
        % (
            m.fs.P2.outlet.pressure[0].value / 1e5,
            m.fs.P2.work_mechanical[0].value / 1e3,
        )
    )


def display_state(m):
    print("---state---")

    def print_state(s, b):
        flow_mass = sum(
            b.flow_mass_phase_comp[0, "Liq", j].value for j in ["H2O", "NaCl"]
        )
        mass_frac_ppm = b.flow_mass_phase_comp[0, "Liq", "NaCl"].value / flow_mass * 1e6
        pressure_bar = b.pressure[0].value / 1e5
        print(
            s
            + ": %.3f kg/s, %.0f ppm, %.1f bar"
            % (flow_mass, mass_frac_ppm, pressure_bar)
        )

    print_state("Feed      ", m.fs.feed.outlet)
    print_state("Split 1   ", m.fs.S1.P1)
    print_state("P1 out    ", m.fs.P1.outlet)
    print_state("Split 2   ", m.fs.S1.PXR)
    print_state("PXR LP out", m.fs.PXR.low_pressure_outlet)
    print_state("P2 out    ", m.fs.P2.outlet)
    print_state("Mix out   ", m.fs.M1.outlet)
    print_state("RO perm   ", m.fs.RO.permeate)
    print_state("RO reten  ", m.fs.RO.retentate)
    print_state("PXR HP out", m.fs.PXR.high_pressure_outlet)


if __name__ == "__main__":
    main()