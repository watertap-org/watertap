#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
from pyomo.environ import (
    ConcreteModel,
    Objective,
    assert_optimal_termination,
    TransformationFactory,
    units as pyunits,
)
from pyomo.network import Arc

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.solvers.get_solver import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.scaling import calculate_scaling_factors
from idaes.core.util.initialization import propagate_state
from idaes.models.unit_models import Product, Feed

from watertap.core.util.initialization import check_dof
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from watertap.unit_models.ion_exchange_0D import IonExchange0D
from watertap.costing import WaterTAPCosting

import math


solver = get_solver()


def main():
    # The IX model currently only supports one "target" ion (i.e., the component in the water source that can be removed by IX)
    # All other ions are inert. This demo does not contain inert ions, but an example can be found in the IX test file:
    # watertap/watertap/unit_models/tests/test_ion_exchange_0D.py
    target_ion = "Ca_2+"
    ions = [target_ion]

    # See ix_build for details on building the model for this demo.
    m = ix_build(ions)
    # See set_operating_conditions for details on operating conditions for this demo.
    set_operating_conditions(m)
    # See initialize_system for details on initializing the models for this demo.
    initialize_system(m)
    # Check the degrees of freedom of the model to ensure it is zero.
    check_dof(m)
    # Solve the model. Store the results in a local variable.
    results = solver.solve(m)
    # Ensure the solve resulted in an optimal termination status.
    assert_optimal_termination(results)
    # Display the degrees of freedom, termination status, and performance metrics of the model.
    print(f"\nDOF = {degrees_of_freedom(m)}")
    print(f"Model solve {results.solver.termination_condition.swapcase()}")
    display_results(m)

    # See optimize_system for details on optimizing this model for a specific condition.
    optimize_system(m)
    ix = m.fs.ion_exchange

    # With our model optimized to new conditions in optimize_system,
    # we can get the new number_columns and bed_depth and fix them in our model.
    num_col = math.ceil(
        ix.number_columns()
    )  # To eliminate fractional number of columns
    bed_depth = ix.bed_depth()
    ix.bed_depth.fix(bed_depth)
    ix.number_columns.fix(num_col)
    check_dof(m)
    results = solver.solve(m)
    assert_optimal_termination(results)
    print(f"\nDOF = {degrees_of_freedom(m)}")
    print(f"Model solve {results.solver.termination_condition.swapcase()}")
    display_results(m)

    return m


def ix_build(ions, target_ion=None, hazardous_waste=False, regenerant="NaCl"):

    if not target_ion:
        target_ion = ions[0]

    # Create the model and flowsheet
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # ion_config is the dictionary needed to configure the property package.
    # For this demo, we only have properties related to the target ion (Ca_2+) and water (H2O)
    ion_props = get_ion_config(ions)

    # The water property package used for the ion exchange model is the multi-component aqueous solution (MCAS) property package
    m.fs.properties = MCASParameterBlock(**ion_props)

    # Add the flowsheet level costing package
    m.fs.costing = WaterTAPCosting()

    # Add feed and product blocks to the flowsheet
    # These are the unit models on the flowsheet that the source water "flows" from/to
    # The must use the same property package as the ion exchange model
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)

    # Configuration dictionary used to instantiate the ion exchange model:
    #   "property_package" indicates which property package to use for the ion exchange model
    #   "target_ion" indicates which ion in the property package will be the reactive ion for the ion exchange model
    #   "hazardous_waste" indicates if the regeneration and spent resin is considered hazardous. If so, it adds costs
    #   "regenerant" indicates the chemical used to regenerate the ion exchange process
    ix_config = {
        "property_package": m.fs.properties,
        "target_ion": target_ion,
        "hazardous_waste": hazardous_waste,
        "regenerant": regenerant,
    }

    # Add the ion exchange model to the flowsheet
    m.fs.ion_exchange = ix = IonExchange0D(**ix_config)

    # Touch concentration properties so they are available for reporting.
    ix.process_flow.properties_in[0].conc_mass_phase_comp[...]
    ix.process_flow.properties_out[0].conc_mass_phase_comp[...]
    m.fs.feed.properties[0].conc_mass_phase_comp[...]
    m.fs.product.properties[0].conc_mass_phase_comp[...]

    # Add costing blocks to the flowsheet
    # Here, the ion exchange model has its own unit-level costing Block
    m.fs.ion_exchange.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing  # Indicating which flowsheet costing block to use to aggregate unit-level costs to the system-level costs
    )
    # Call cost_process() method to create system-wide global parameters and add aggregating constraints to costing model
    m.fs.costing.cost_process()
    # Designate the volumetric flow on the Product block to be the stream used as the annual water production
    m.fs.costing.add_annual_water_production(
        m.fs.product.properties[0].flow_vol_phase["Liq"]
    )
    # Add LCOW variable to costing block
    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol_phase["Liq"])
    # Add specific energy consumption variable to costing block
    m.fs.costing.add_specific_energy_consumption(
        m.fs.product.properties[0].flow_vol_phase["Liq"]
    )

    # Arcs are used to "connect" Ports on unit process models to Ports on other unit process models
    # For example, in this next line the outlet Port on the Feed model is connected to the inlet Port on the ion exchange model
    m.fs.feed_to_ix = Arc(source=m.fs.feed.outlet, destination=ix.inlet)
    m.fs.ix_to_product = Arc(source=ix.outlet, destination=m.fs.product.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m)

    # Scaling variables in the model
    # Here, the molar flow for water ("flow_mol_phase_comp[Liq, H2O]") on the Feed block is scaled by 1e-4.
    # This is because the molar flow rate of water in this demo is ~2777 mol/s
    # and scaling factors are chosen such that the value of the variable multiplied by the scaling factor is ~1
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e-4, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 10, index=("Liq", target_ion)
    )
    # Call calculate_scaling_factors to apply scaling factors for each variable that we haven't set scaling factors for above.
    calculate_scaling_factors(m)

    return m


def set_operating_conditions(m, flow_in=0.05, conc_mass_in=0.1, solver=None):
    if solver is None:
        solver = get_solver()
    ix = m.fs.ion_exchange
    target_ion = ix.config.target_ion

    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): flow_in,  # m3/s
            ("conc_mass_phase_comp", ("Liq", target_ion)): conc_mass_in,  # kg/m3
            ("pressure", None): 101325,
            ("temperature", None): 298,
        },
        hold_state=True,
    )
    # Fix key decision variables in ion exchange model.
    ix.langmuir[target_ion].fix(0.7)
    ix.resin_max_capacity.fix(3)
    ix.service_flow_rate.fix(15)
    # Note: number_columns and bed_depth are typically not known a priori for this model.
    # They can be determined by first running the model without these variables fixed.
    ix.number_columns.fix(4)
    ix.bed_depth.fix(1.7)

    # Fix remaining variables.
    # Using the .fix() method on a Var fixes the variable to its initialized value.
    ix.resin_diam.fix()
    ix.resin_bulk_dens.fix()
    ix.bed_porosity.fix()
    ix.dimensionless_time.fix()


def initialize_system(m):
    # First we initialize the Feed block using values set in set_operating_conditions
    m.fs.feed.initialize()

    # We then propagate the state of the Feed block to the ion exchange model...
    propagate_state(m.fs.feed_to_ix)
    # ... and then initialize the ion exchange model.
    m.fs.ion_exchange.initialize()
    # With the ion exchange model initialized, we have initial guesses for the Product block
    # and can propagate the state of the IX effluent stream.
    propagate_state(m.fs.ix_to_product)
    # Finally, we initialize the product and costing blocks.
    m.fs.product.initialize()
    m.fs.costing.initialize()


def optimize_system(m):
    # Example of optimizing number of IX columns based on desired effluent equivalent concentration

    # Adding an objective to model.
    # In this case, we want to optimze the model to minimize the LCOW.
    m.fs.obj = Objective(expr=m.fs.costing.LCOW)
    ix = m.fs.ion_exchange
    target_ion = m.fs.ion_exchange.config.target_ion

    # For this demo, we are optimizing the model to have an effluent concentration of 25 mg/L.
    # Our initial model resulted in an effluent concentration of 0.21 mg/L.
    # By increasing the effluent concentration, we will have a longer breakthrough time, which will lead to less regeneration solution used,
    # and (hopefully) a lower cost.
    ix.process_flow.properties_out[0].conc_mass_phase_comp["Liq", target_ion].fix(0.025)

    # With the new effluent conditions for our ion exchange model, this will have implications for our downstream models (the Product block)
    # Thus, we must re-propagate the new effluent state to these models...
    propagate_state(m.fs.ix_to_product)
    # ...and re-initialize them to our new conditions.
    m.fs.product.initialize()

    # To adjust solution to fixed-pattern to achieve desired effluent, must unfix dimensionless_time.
    ix.dimensionless_time.unfix()
    # Can optimize around different design variables, e.g., bed_depth, service_flow_rate (or combinations of these)
    # Here demonstrates optimization around column design
    ix.number_columns.unfix()
    ix.bed_depth.unfix()
    solver.solve(m)


def get_ion_config(ions):

    if not isinstance(ions, (list, tuple)):
        ions = [ions]
    diff_data = {
        "Na_+": 1.33e-9,
        "Ca_2+": 9.2e-10,
        "Cl_-": 2.03e-9,
        "Mg_2+": 0.706e-9,
        "SO4_2-": 1.06e-9,
    }
    mw_data = {
        "Na_+": 23e-3,
        "Ca_2+": 40e-3,
        "Cl_-": 35e-3,
        "Mg_2+": 24e-3,
        "SO4_2-": 96e-3,
    }
    charge_data = {"Na_+": 1, "Ca_2+": 2, "Cl_-": -1, "Mg_2+": 2, "SO4_2-": -2}
    ion_config = {
        "solute_list": [],
        "diffusivity_data": {},
        "mw_data": {"H2O": 18e-3},
        "charge": {},
    }
    for ion in ions:
        ion_config["solute_list"].append(ion)
        ion_config["diffusivity_data"][("Liq", ion)] = diff_data[ion]
        ion_config["mw_data"][ion] = mw_data[ion]
        ion_config["charge"][ion] = charge_data[ion]
    return ion_config


def display_results(m):

    ix = m.fs.ion_exchange
    liq = "Liq"
    header = f'{"PARAM":<40s}{"VALUE":<40s}{"UNITS":<40s}\n'

    prop_in = ix.process_flow.properties_in[0]
    prop_out = ix.process_flow.properties_out[0]

    recovery = prop_out.flow_vol_phase["Liq"]() / prop_in.flow_vol_phase["Liq"]()
    target_ion = ix.config.target_ion
    ion_set = ix.config.property_package.ion_set
    bv_to_regen = (ix.vel_bed() * ix.t_breakthru()) / ix.bed_depth()

    title = f'\n{"=======> SUMMARY <=======":^80}\n'
    print(title)
    print(header)
    print(f'{"LCOW":<40s}{f"{m.fs.costing.LCOW():<40.4f}"}{"$/m3":<40s}')
    print(
        f'{"TOTAL Capital Cost":<40s}{f"${ix.costing.capital_cost():<39,.2f}"}{"$":<40s}'
    )
    print(
        f'{"Specific Energy Consumption":<40s}{f"{m.fs.costing.specific_energy_consumption():<39,.5f}"}{"kWh/m3":<40s}'
    )
    print(
        f'{f"Annual Regenerant cost ({ix.config.regenerant})":<40s}{f"${m.fs.costing.aggregate_flow_costs[ix.config.regenerant]():<39,.2f}"}{"$/yr":<40s}'
    )
    print(f'{"BV Until Regen":<40s}{bv_to_regen:<40.3f}{"Bed Volumes":<40s}')
    print(
        f'{f"Breakthrough/Initial Conc. [{target_ion}]":<40s}{ix.c_norm[target_ion]():<40.3%}'
    )
    print(
        f'{"Vol. Flow In [m3/s]":<40s}{prop_in.flow_vol_phase[liq]():<40.5f}{"m3/s":<40s}'
    )
    print(
        f'{"Vol. Flow Out [m3/s]":<40s}{prop_out.flow_vol_phase[liq]():<40.5f}{"m3/s":<40s}'
    )
    print(f'{"Water Vol. Recovery":<40s}{recovery:<40.2%}{"%":<40s}')
    print(f'{"Breakthrough Time [hr]":<40s}{ix.t_breakthru() / 3600:<40.3f}{"hr":<40s}')
    print(f'{"Number Columns":<40s}{ix.number_columns():<40.2f}{"---":<40s}')
    print(f'{"Column Vol.":<40s}{ix.col_vol_per():<40.2f}{"m3":<40s}')
    print(f'{"Bed Depth":<40s}{ix.bed_depth():<40.2f}{"m":<40s}')
    for ion in ion_set:
        print(
            f'{f"Removal [{ion}]":<40s}{1 - prop_out.conc_mass_phase_comp[liq, ion]() / prop_in.conc_mass_phase_comp[liq, ion]():<40.4%}{"%":<40s}'
        )
        print(
            f'{f"Conc. In [{ion}, mg/L]":<40s}{pyunits.convert(prop_in.conc_mass_phase_comp[liq, ion], to_units=pyunits.mg/pyunits.L)():<40.3e}{"mg/L":<40s}'
        )
        print(
            f'{f"Conc. Out [{ion}, mg/L]":<40s}{pyunits.convert(prop_out.conc_mass_phase_comp[liq, ion], to_units=pyunits.mg/pyunits.L)():<40.3e}{"mg/L":<40s}'
        )


if __name__ == "__main__":
    m = main()
