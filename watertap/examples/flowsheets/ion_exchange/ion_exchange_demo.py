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
from idaes.core.util.scaling import set_scaling_factor, calculate_scaling_factors
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

    # The mass fraction is total fraction of influent mass of the target ion
    mass_frac = 1e-4
    feed_mass_frac = {target_ion: mass_frac}

    # See ix_build for details on building the model for this demo.
    m = ix_build(ions)
    # See set_operating_conditions for details on operating conditions for this demo.
    set_operating_conditions(m, feed_mass_frac=feed_mass_frac)
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
    ion_config = get_ion_config(ions)

    # The water property package used for the ion exchange model is the multi-component aqueous solution (MCAS) property package
    m.fs.properties = MCASParameterBlock(**ion_config)

    # Add the flowsheet level costing package
    m.fs.costing = WaterTAPCosting()

    # Add feed, product, and regeneration blocks to the flowsheet
    # These are the unit models on the flowsheet that the source water "flows" from/to
    # The must use the same property package as the ion exchange model
    m.fs.feed = feed = Feed(property_package=m.fs.properties)
    m.fs.product = prod = Product(property_package=m.fs.properties)
    m.fs.regen = regen = Product(property_package=m.fs.properties)

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

    # Touch properties so they are available for initialization and arc propagation...
    # These must be available on each of these models (feed, product, regeneration) because these properties are used in the ion exchange model
    # and are not state variables, so aren't created automatically upon instantiation of the model.
    feed.properties[0].flow_vol_phase[...]
    feed.properties[0].conc_equiv_phase_comp[...]
    prod.properties[0].flow_vol_phase[...]
    prod.properties[0].conc_equiv_phase_comp[...]
    regen.properties[0].flow_vol_phase[...]
    regen.properties[0].conc_equiv_phase_comp[...]

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
    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol)
    # Add specific energy consumption variable to costing block
    m.fs.costing.add_specific_energy_consumption(
        m.fs.product.properties[0].flow_vol_phase["Liq"]
    )

    # Arcs are used to "connect" Ports on unit process models to Ports on other unit process models
    # For example, in this next line the outlet Port on the Feed model is connected to the inlet Port on the ion exchange model
    m.fs.feed_to_ix = Arc(source=m.fs.feed.outlet, destination=ix.inlet)
    m.fs.ix_to_product = Arc(source=ix.outlet, destination=m.fs.product.inlet)
    m.fs.ix_to_regen = Arc(source=ix.regen, destination=m.fs.regen.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m)

    # Scaling variables in the model
    # Here, the molar flow for water ("flow_mol_phase_comp[Liq, H2O]") on the Feed block is scaled by 1e-3.
    # This is because the molar flow rate of water in this demo is ~2777 mol/s
    # and scaling factors are chosen such that the value of the variable multiplied by the scaling factor is ~1
    set_scaling_factor(feed.properties[0].flow_mol_phase_comp["Liq", "H2O"], 1e-3)
    # Similarly, the molar flow rate of our target ion (Ca_2+) is ~0.125 mol/s, so we scale by a factor of 10.
    set_scaling_factor(feed.properties[0].flow_mol_phase_comp["Liq", target_ion], 10)

    # We use the same scaling factor for the molar flow of water and Ca_2+ (target_ion) for the ion exchange model
    # because the Feed block flows directly to the ion exchange model so the molar flow rate is equivalent
    set_scaling_factor(ix.properties_in[0].flow_mol_phase_comp["Liq", "H2O"], 1e-3)
    set_scaling_factor(ix.properties_in[0].flow_mol_phase_comp["Liq", target_ion], 10)

    # We expect the product molar flow rate of water to be approximately the same as the influent, so we use the same scaling factor.
    set_scaling_factor(prod.properties[0].flow_mol_phase_comp["Liq", "H2O"], 1e-3)
    # However, we expect the product molar flow rate of Ca_2+ to be much smaller than the influent, so we use a larger scaling factor.
    set_scaling_factor(prod.properties[0].flow_mol_phase_comp["Liq", target_ion], 1e6)

    set_scaling_factor(ix.properties_out[0].flow_mol_phase_comp["Liq", "H2O"], 1e-3)
    set_scaling_factor(ix.properties_out[0].flow_mol_phase_comp["Liq", target_ion], 1e6)

    # We expect the regeneration molar flow rate of water to be much smaller than the influent, so we use a larger scaling factor.
    set_scaling_factor(regen.properties[0].flow_mol_phase_comp["Liq", "H2O"], 1)
    # However, we expect the regeneration molar flow rate of Ca_2+ to be much higher than the influent, so we use a smaller scaling factor.
    # Similar logic is applied to other property variables for each stream.
    set_scaling_factor(regen.properties[0].flow_mol_phase_comp["Liq", target_ion], 0.01)
    set_scaling_factor(regen.properties[0].conc_mol_phase_comp["Liq", target_ion], 0.01)
    set_scaling_factor(
        regen.properties[0].conc_mass_phase_comp["Liq", target_ion], 0.01
    )
    set_scaling_factor(regen.properties[0].mass_frac_phase_comp["Liq", target_ion], 100)

    set_scaling_factor(ix.properties_regen[0].flow_mol_phase_comp["Liq", "H2O"], 1)
    set_scaling_factor(
        ix.properties_regen[0].flow_mol_phase_comp["Liq", target_ion], 0.01
    )
    set_scaling_factor(
        ix.properties_regen[0].conc_mol_phase_comp["Liq", target_ion], 0.01
    )
    set_scaling_factor(
        ix.properties_regen[0].conc_mass_phase_comp["Liq", target_ion], 1
    )
    set_scaling_factor(
        ix.properties_regen[0].mass_frac_phase_comp["Liq", target_ion], 100
    )

    # Call calculate_scaling_factors to apply scaling factors for each variable that we haven't set scaling factors for above.
    calculate_scaling_factors(m)

    return m


def set_operating_conditions(m, feed_mass_frac={}, mass_flow_in=50, solver=None):
    if solver is None:
        solver = get_solver()
    ix = m.fs.ion_exchange
    feed = m.fs.feed
    target_ion = ix.config.target_ion
    # Adding pyunits to mass_flow_in
    mass_flow_in = mass_flow_in * (pyunits.kg / pyunits.s)

    # For each ion in the feed_mass_frac dictionary (in this demo, it is only Ca_2+)...
    for ion, mass_frac in feed_mass_frac.items():
        # ... calculate the molar flow rate for that ion using the mass fraction...
        mol_flow = (
            mass_frac
            * (pyunits.kg / pyunits.kg)
            * mass_flow_in
            / feed.config.property_package.mw_comp[ion]
        )
        # ... and fix the molar flow rate for that ion in the Feed block.
        feed.properties[0].flow_mol_phase_comp["Liq", ion].fix(mol_flow)

    # Calculate the mass fraction and molar flow rate of water and fix those variables in Feed block.
    h2o_mass_frac = 1 - sum(x for x in feed_mass_frac.values())
    h2o_mol_flow = (
        h2o_mass_frac
        * (pyunits.kg / pyunits.kg)
        * mass_flow_in
        / feed.config.property_package.mw_comp["H2O"]
    )
    feed.properties[0].flow_mol_phase_comp["Liq", "H2O"].fix(h2o_mol_flow)
    # Fix remaining state variables in Feed block.
    # Note: these correspond to the property package, not the unit model.
    # Temperature and pressure variables are not used in the IonExchange0D model,
    # but because they are part of the property model, we must fix them or we end up with degrees of freedom.
    feed.properties[0].pressure.fix(101325)
    feed.properties[0].temperature.fix(298)

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
    ix.regen_dose.fix()
    ix.regen_recycle.fix()
    ix.t_regen.fix()
    ix.t_bw.fix()
    ix.bw_rate.fix()
    ix.rinse_bv.fix()
    ix.pump_efficiency.fix()
    ix.p_drop_A.fix()
    ix.p_drop_B.fix()
    ix.p_drop_C.fix()
    ix.bed_expansion_frac_A.fix()
    ix.bed_expansion_frac_B.fix()
    ix.bed_expansion_frac_C.fix()
    ix.service_to_regen_flow_ratio.fix()
    ix.number_columns_redund.fix()


def initialize_system(m):
    target_ion = m.fs.ion_exchange.config.target_ion
    # First we initialize the Feed block using values set in set_operating_conditions
    m.fs.feed.initialize()
    m.fs.feed.properties[0].flow_vol_phase["Liq"].fix()
    m.fs.feed.properties[0].conc_equiv_phase_comp["Liq", target_ion].fix()

    # We then propagate the state of the Feed block to the ion exchange model...
    propagate_state(m.fs.feed_to_ix)
    # ... and then initialize the ion exchange model.
    m.fs.ion_exchange.initialize()
    # With the ion exchange model initialized, we have initial guesses for the Product and Regen blocks
    # and can propagate the state of the IX effluent and regen streams to those models.
    propagate_state(m.fs.ix_to_regen)
    propagate_state(m.fs.ix_to_product)
    # Finally, we initialize the regen, product, and costing blocks.
    m.fs.regen.initialize()
    m.fs.product.initialize()
    m.fs.costing.initialize()


def optimize_system(m):
    # Example of optimizing number of IX columns based on desired effluent equivalent concentration

    # Adding an objective to model.
    # In this case, we want to optimze the model to minimize the LCOW.
    m.fs.obj = Objective(expr=m.fs.costing.LCOW)
    ix = m.fs.ion_exchange
    target_ion = m.fs.ion_exchange.config.target_ion

    # For this demo, we are optimizing the model to have an effluent concentration of 0.5 equivalents/m3.
    # Our initial model resulted in an effluent concentration of 0.011 equivalent/m3.
    # By increasing the effluent concentration, we will have a longer breakthrough time, which will lead to less regeneration solution used,
    # and (hopefully) a lower cost.
    ix.properties_out[0].conc_equiv_phase_comp["Liq", target_ion].fix(0.5)

    # With the new effluent conditions for our ion exchange model, this will have implications for our downstream models (the Product and Regeneration blocks)
    # Thus, we must re-propagate the new effluent state to these models...
    propagate_state(m.fs.ix_to_product)
    propagate_state(m.fs.ix_to_regen)
    # ...and re-initialize them to our new conditions.
    m.fs.product.initialize()
    m.fs.regen.initialize()

    # To adjust solution to fixed-pattern to achieve desired effluent, must unfix dimensionless_time.
    ix.dimensionless_time.unfix()
    # Can optimize around different design variables, e.g., bed_depth, service_flow_rate (or combinations of these)
    # Here demonstrates optimization around column design
    ix.number_columns.unfix()
    ix.bed_depth.unfix()
    solver.solve(ix)


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

    prop_in = ix.properties_in[0]
    prop_out = ix.properties_out[0]
    prop_regen = ix.properties_regen[0]

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
        f'{"Specific Energy Consumption":<40s}{f"${m.fs.costing.specific_energy_consumption():<39,.2f}"}{"kWh/m3":<40s}'
    )
    print(
        f'{f"Annual Regenerant cost ({ix.regen_chem})":<40s}{f"${m.fs.costing.aggregate_flow_costs[ix.regen_chem]():<39,.2f}"}{"$/yr":<40s}'
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
    print(
        f'{"Vol. Flow Regen [m3/s]":<40s}{prop_regen.flow_vol_phase[liq]():<40.5f}{"m3/s":<40s}'
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
        print(
            f'{f"Conc. Regen [{ion}, mg/L]":<40s}{pyunits.convert(prop_regen.conc_mass_phase_comp[liq, ion], to_units=pyunits.mg/pyunits.L)():<40.3e}{"mg/L":<40s}'
        )


if __name__ == "__main__":
    m = main()
