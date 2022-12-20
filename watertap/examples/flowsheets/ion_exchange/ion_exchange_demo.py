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

import numpy as np


solver = get_solver()


def main():
    target_ion = "Ca_2+"
    ions = [target_ion]
    mass_frac = 1e-4
    feed_mass_frac = {target_ion: mass_frac}
    m = ix_build(ions)
    set_operating_conditions(m, feed_mass_frac=feed_mass_frac)
    initialize_system(m)
    check_dof(m)
    results = solver.solve(m)
    assert_optimal_termination(results)
    print(f"\nDOF = {degrees_of_freedom(m)}")
    print(f"Model solve {results.solver.termination_condition.swapcase()}")
    display_results(m)

    optimize_system(m)
    ix = m.fs.ion_exchange
    num_col = np.ceil(ix.number_columns())  # To eliminate fractional number of columns
    bed_depth = ix.bed_depth()
    ix.bed_depth.fix(bed_depth)
    ix.number_columns.fix(num_col)
    check_dof(m)
    results = solver.solve(m)
    assert_optimal_termination(results)
    print(f"\nDOF = {degrees_of_freedom(m)}")
    print(f"Model solve {results.solver.termination_condition.swapcase()}")
    display_results(m)


def ix_build(ions, target_ion=None, hazardous_waste=False, regenerant="NaCl"):

    if not target_ion:
        target_ion = ions[0]

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    ion_config = get_ion_config(ions)
    m.fs.properties = MCASParameterBlock(**ion_config)
    m.fs.costing = WaterTAPCosting()
    m.fs.feed = feed = Feed(property_package=m.fs.properties)
    m.fs.product = prod = Product(property_package=m.fs.properties)
    m.fs.regen = regen = Product(property_package=m.fs.properties)

    ix_config = {
        "property_package": m.fs.properties,
        "target_ion": target_ion,
        "hazardous_waste": hazardous_waste,
        "regenerant": regenerant,
    }

    m.fs.ion_exchange = ix = IonExchange0D(**ix_config)

    ## Touch properties so they are available for initialization and arc propagation...
    feed.properties[0].flow_vol_phase[...]
    feed.properties[0].conc_equiv_phase_comp[...]
    prod.properties[0].flow_vol_phase[...]
    prod.properties[0].conc_equiv_phase_comp[...]
    regen.properties[0].flow_vol_phase[...]
    regen.properties[0].conc_equiv_phase_comp[...]

    # costing
    m.fs.ion_exchange.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing
    )
    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(
        m.fs.product.properties[0].flow_vol_phase["Liq"]
    )
    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(
        m.fs.product.properties[0].flow_vol_phase["Liq"]
    )

    m.fs.feed_to_ix = Arc(source=m.fs.feed.outlet, destination=ix.inlet)
    m.fs.ix_to_product = Arc(source=ix.outlet, destination=m.fs.product.inlet)
    m.fs.ix_to_regen = Arc(source=ix.regen, destination=m.fs.regen.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m)

    set_scaling_factor(feed.properties[0].flow_mol_phase_comp["Liq", "H2O"], 1e-3)
    set_scaling_factor(feed.properties[0].flow_mol_phase_comp["Liq", target_ion], 10)

    set_scaling_factor(ix.properties_in[0].flow_mol_phase_comp["Liq", "H2O"], 1e-3)
    set_scaling_factor(ix.properties_in[0].flow_mol_phase_comp["Liq", target_ion], 10)

    set_scaling_factor(prod.properties[0].flow_mol_phase_comp["Liq", "H2O"], 1e-3)
    set_scaling_factor(prod.properties[0].flow_mol_phase_comp["Liq", target_ion], 1e6)

    set_scaling_factor(ix.properties_out[0].flow_mol_phase_comp["Liq", "H2O"], 1e-3)
    set_scaling_factor(ix.properties_out[0].flow_mol_phase_comp["Liq", target_ion], 1e6)

    set_scaling_factor(regen.properties[0].flow_mol_phase_comp["Liq", "H2O"], 1)
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

    calculate_scaling_factors(m)

    return m


def set_operating_conditions(m, feed_mass_frac={}, mass_flow_in=50, solver=None):
    if solver is None:
        solver = get_solver()
    ix = m.fs.ion_exchange
    feed = m.fs.feed
    target_ion = ix.config.target_ion
    mass_flow_in = mass_flow_in * (pyunits.kg / pyunits.s)

    for ion, mass_frac in feed_mass_frac.items():
        mol_flow = (
            mass_frac
            * (pyunits.kg / pyunits.kg)
            * mass_flow_in
            / feed.config.property_package.mw_comp[ion]
        )
        feed.properties[0].flow_mol_phase_comp["Liq", ion].fix(mol_flow)

    h2o_mass_frac = 1 - sum(x for x in feed_mass_frac.values())
    h2o_mol_flow = (
        h2o_mass_frac
        * (pyunits.kg / pyunits.kg)
        * mass_flow_in
        / feed.config.property_package.mw_comp["H2O"]
    )
    feed.properties[0].flow_mol_phase_comp["Liq", "H2O"].fix(h2o_mol_flow)
    feed.properties[0].pressure.fix(101325)
    feed.properties[0].temperature.fix(298)

    ix.langmuir[target_ion].fix(0.7)
    ix.resin_max_capacity.fix(3)
    ix.service_flow_rate.fix(15)
    ix.number_columns.fix(4)
    ix.bed_depth.fix(1.7)
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
    m.fs.feed.initialize()
    m.fs.feed.properties[0].flow_vol_phase["Liq"].fix()
    m.fs.feed.properties[0].conc_equiv_phase_comp["Liq", target_ion].fix()
    propagate_state(m.fs.feed_to_ix)
    m.fs.ion_exchange.initialize()
    propagate_state(m.fs.ix_to_regen)
    propagate_state(m.fs.ix_to_product)
    m.fs.regen.initialize()
    m.fs.product.initialize()
    m.fs.costing.initialize()


def optimize_system(m):
    # Example of optimizing number of IX columns based on desired effluent equivalent concentration
    m.fs.obj = Objective(expr=m.fs.costing.LCOW)
    ix = m.fs.ion_exchange
    target_ion = m.fs.ion_exchange.config.target_ion
    ix.properties_out[0].conc_equiv_phase_comp["Liq", target_ion].fix(0.5)
    propagate_state(m.fs.ix_to_product)
    propagate_state(m.fs.ix_to_regen)
    m.fs.product.initialize()
    m.fs.regen.initialize()

    # To adjust solution to fixed-pattern to achieve desired effluent, must unfix dimensionless_time
    ix.dimensionless_time.unfix()
    # Can optimize around different design variables, e.g., bed_depth, service_flow_rate (or combinations of these)
    # Here demonstrates optimization around column design
    ix.number_columns.unfix()
    ix.bed_depth.unfix()
    solver.solve(ix)
    # num_col = np.ceil(ix.number_columns())  # To eliminate fractional number of columns
    # bed_depth = ix.bed_depth()
    # ix.bed_depth.fix(bed_depth)
    # ix.number_columns.fix(num_col)


def get_ion_config(ions):

    if not isinstance(ions, list):
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
    main()
