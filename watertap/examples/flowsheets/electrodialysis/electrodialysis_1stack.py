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
)
from pyomo.network import Arc

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.solvers import get_solver
from idaes.core.util.initialization import propagate_state
from idaes.core.util.model_statistics import degrees_of_freedom, report_statistics
from idaes.models.unit_models import Feed, Product, Separator
from pandas import DataFrame
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslogger
from watertap.core.util.initialization import check_dof
from watertap.unit_models.electrodialysis_1D import Electrodialysis1D
from watertap.costing.watertap_costing_package import (
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

    print("\n***---Simulation results---***")
    display_model_metrics(m)

    # Perform an optimization over selected variables
    initialize_system(m, solver=solver)
    optimize_system(m, solver=solver)
    print("\n***---Optimization results---***")
    display_model_metrics(m)


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
    m.fs.separator = Separator(
        property_package=m.fs.properties,
        outlet_list=["inlet_diluate", "inlet_concentrate"],
    )  # "inlet_diluate" and "inlet_concentrate" are two separator's outlet ports that are connected to the two inlets of the ED stack.

    # Add electrodialysis (ED) stacks
    m.fs.EDstack = Electrodialysis1D(
        property_package=m.fs.properties,
        operation_mode="Constant_Voltage",
        finite_elements=20,
    )
    m.fs.product = Product(property_package=m.fs.properties)
    m.fs.disposal = Product(property_package=m.fs.properties)

    # Touching needed variables for initialization and displaying results
    m.fs.feed.properties[0].conc_mass_phase_comp[...]
    m.fs.separator.inlet_diluate_state[0].conc_mass_phase_comp[...]
    m.fs.separator.inlet_concentrate_state[0].conc_mass_phase_comp[...]
    m.fs.product.properties[0].conc_mass_phase_comp[...]
    m.fs.disposal.properties[0].conc_mass_phase_comp[...]

    m.fs.feed.properties[0].flow_vol_phase[...]
    m.fs.separator.inlet_diluate_state[0].flow_vol_phase[...]
    m.fs.separator.inlet_concentrate_state[0].flow_vol_phase[...]
    m.fs.product.properties[0].flow_vol_phase[...]
    m.fs.disposal.properties[0].flow_vol_phase[...]

    # costing
    m.fs.EDstack.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(
        m.fs.product.properties[0].flow_vol_phase["Liq"]
    )
    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(
        m.fs.product.properties[0].flow_vol_phase["Liq"]
    )

    # Add two variable for reporting
    m.fs.mem_area = Var(
        initialize=1,
        bounds=(0, 1e3),
        units=pyunits.meter**2,
        doc="Total membrane area for cem (or aem) in one stack",
    )
    m.fs.product_salinity = Var(
        initialize=1, bounds=(0, 1000), units=pyunits.kg * pyunits.meter**-3
    )
    m.fs.disposal_salinity = Var(
        initialize=1, bounds=(0, 1e6), units=pyunits.kg * pyunits.meter**-3
    )
    m.fs.eq_product_salinity = Constraint(
        expr=m.fs.product_salinity
        == sum(
            m.fs.product.properties[0].conc_mass_phase_comp["Liq", j]
            for j in m.fs.properties.ion_set | m.fs.properties.solute_set
        )
    )
    m.fs.eq_disposal_salinity = Constraint(
        expr=m.fs.disposal_salinity
        == sum(
            m.fs.disposal.properties[0].conc_mass_phase_comp["Liq", j]
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
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.separator.inlet)
    m.fs.s02 = Arc(
        source=m.fs.separator.inlet_diluate, destination=m.fs.EDstack.inlet_diluate
    )
    m.fs.s03 = Arc(
        source=m.fs.separator.inlet_concentrate,
        destination=m.fs.EDstack.inlet_concentrate,
    )
    m.fs.s04 = Arc(source=m.fs.EDstack.outlet_diluate, destination=m.fs.product.inlet)
    m.fs.s05 = Arc(
        source=m.fs.EDstack.outlet_concentrate, destination=m.fs.disposal.inlet
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
    m.fs.feed.properties[0].temperature.fix(273.15 + 25)  # feed temperature [K]
    m.fs.feed.properties[0].flow_mol_phase_comp["Liq", "H2O"].fix(4.8)
    m.fs.feed.properties[0].flow_mol_phase_comp["Liq", "Na_+"].fix(1.476e-2)
    m.fs.feed.properties[0].flow_mol_phase_comp["Liq", "Cl_-"].fix(1.476e-2)

    # Fix separator's split_fraction to 0.5, i.e., equal flows into the diluate and concentrate channels
    m.fs.separator.split_fraction[0, "inlet_diluate"].fix(0.5)

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
    propagate_state(m.fs.s01)
    m.fs.separator.initialize(optarg=optarg)
    propagate_state(m.fs.s02)
    propagate_state(m.fs.s03)
    m.fs.EDstack.initialize(optarg=optarg)
    propagate_state(m.fs.s04)
    m.fs.product.initialize()
    propagate_state(m.fs.s05)
    m.fs.disposal.initialize()
    m.fs.costing.initialize()


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
    m.fs.product.properties[0].conc_mass_phase_comp["Liq", "Na_+"].fix(0.393)

    print("---report model statistics---\n ", report_statistics(m.fs))
    if solver is None:
        solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)


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
            value(m.fs.product.properties[0].flow_vol_phase["Liq"]),
            value(m.fs.product_salinity),
        ],
        "Disposal": [
            value(m.fs.disposal.properties[0].flow_vol_phase["Liq"]),
            value(m.fs.disposal_salinity),
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
