#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

from pyomo.environ import (
    ConcreteModel,
    value,
    assert_optimal_termination,
    check_optimal_termination,
    units as pyunits,
)

from idaes.core.util.initialization import solve_indexed_blocks
import idaes.core.util.scaling as iscale
from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.core.solvers import get_solver
from watertap.core.util.model_diagnostics.infeasible import (
    print_close_to_bounds,
    print_infeasible_bounds,
    print_infeasible_constraints,
)


def calculate_operating_pressure(
    feed_state_block=None,
    over_pressure=0.15,
    water_recovery=0.5,
    salt_passage=0.01,
    solver=None,
):
    t = ConcreteModel()
    prop = feed_state_block.config.parameters
    comp = prop.solute_set.at(1)
    t.brine = prop.build_state_block([0])

    t.brine[0].flow_mass_phase_comp["Liq", "H2O"].fix(
        value(feed_state_block.flow_mass_phase_comp["Liq", "H2O"])
        * (1 - water_recovery)
    )
    t.brine[0].flow_mass_phase_comp["Liq", comp].fix(
        value(feed_state_block.flow_mass_phase_comp["Liq", comp]) * (1 - salt_passage)
    )
    t.brine[0].pressure.fix(
        101325
    )  # valid when osmotic pressure is independent of hydraulic pressure
    t.brine[0].temperature.fix(value(feed_state_block.temperature))

    t.brine[0].pressure_osm_phase
    results = solve_indexed_blocks(solver, [t.brine])
    assert_optimal_termination(results)

    return value(t.brine[0].pressure_osm_phase["Liq"]) * (1 + over_pressure)


def report_pump(blk, w=30):
    title = "Pump Report"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")

    pin = blk.pump.control_volume.properties_in[0].pressure
    deltaP = blk.pump.deltaP[0]
    pout = blk.pump.control_volume.properties_out[0].pressure
    work = blk.pump.work_mechanical[0]

    print(
        f'{f"Inlet Pressure (bar)":<{w}s}{value(pyunits.convert(pin, to_units=pyunits.bar)):<{w}.3f}{"bar"}'
    )
    print(
        f'{f"Pressure Change (bar)":<{w}s}{value(pyunits.convert(deltaP, to_units=pyunits.bar)):<{w}.3f}{"bar"}'
    )
    print(
        f'{f"Outlet Pressure (bar)":<{w}s}{value(pyunits.convert(pout, to_units=pyunits.bar)):<{w}.3f}{"bar"}'
    )
    print(
        f'{f"Power (kW)":<{w}s}{value(pyunits.convert(work, to_units=pyunits.kW)):<{w}.3f}{"kW"}'
    )
    print(f'{f"Efficiency (-)":<{w}s}{value(blk.pump.efficiency_pump[0]):<{w}.3f}{"-"}')


def report_ro(blk, w=30):
    comp = blk.RO.config.property_package.solute_set.at(1)
    title = "RO Report"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")

    print(
        f'{f"Inlet Pressure":<{w}s}{value(pyunits.convert(blk.RO.feed_side.properties[0, 0].pressure, to_units=pyunits.bar)):<{w}.3f}{f"bar"}'
    )
    print(
        f'{f"Inlet Flow":<{w}s}{value(blk.feed.properties[0].flow_vol_phase["Liq"])*1000:<{w}.3f}{"L/s"}'
    )
    print(
        f'{f"Brine Flow":<{w}s}{value(blk.disposal.properties[0].flow_vol_phase["Liq"])*1000:<{w}.3f}{"L/s"}'
    )
    print(
        f'{f"Product Flow":<{w}s}{value(blk.product.properties[0].flow_vol_phase["Liq"])*1000:<{w}.3f}{"L/s"}'
    )
    print(f'{f"Recovery":<{w}s}{value(blk.recovery)*100:<{w}.3f}{"%"}')
    print(
        f'{f"Inlet Conc.":<{w}s}{value(blk.feed.properties[0].conc_mass_phase_comp["Liq", comp]):<{w}.3f}{"g/L"}'
    )
    print(
        f'{f"Perm Conc.":<{w}s}{value(blk.product.properties[0].conc_mass_phase_comp["Liq", comp]):<{w}.3f}{"g/L"}'
    )
    print(
        f'{f"Brine Conc.":<{w}s}{value(blk.disposal.properties[0].conc_mass_phase_comp["Liq", comp]):<{w}.3f}{"g/L"}'
    )
    print(
        f'{f"Rejection":<{w}s}{value(blk.RO.rejection_phase_comp[0, "Liq", comp])*100:<{w}.3f}{"%"}'
    )
    print(
        f'{f"deltaP":<{w}s}{value(pyunits.convert(blk.RO.deltaP[0], to_units=pyunits.bar)):<{w}.3f}{f"bar"}'
    )
    acomp = blk.RO.A_comp[0, "H2O"]
    print(
        f'{f"Water Perm (A)":<{w}s}{value(pyunits.convert(acomp, to_units=pyunits.liter / pyunits.m**2 / pyunits.hour / pyunits.bar)):<{w}.3f}{f"LMH/bar"}'
    )
    print(f'{f"Water Perm (A)":<{w}s}{value(acomp):<{w}.3e}{f"m/s/Pa"}')
    bcomp = blk.RO.B_comp[0, comp]
    print(
        f'{f"Salt Perm (B)":<{w}s}{value(pyunits.convert(bcomp, to_units=pyunits.liter / pyunits.m**2 / pyunits.hour)):<{w}.3f}{f"LMH"}'
    )
    print(f'{f"Salt Perm (B)":<{w}s}{value(bcomp):<{w}.3e}{f"m/s"}')
    print(f'{f"Average Flux (LMH)":<{w}s}{value(blk.average_flux_LMH):<{w}.3f}{f"LMH"}')

    print(
        f'{f"Membrane Area":<{w}s}{value(blk.RO.area):<{w}.3f}{f"{pyunits.get_units(blk.RO.area)}"}'
    )
    print(
        f'{f"Membrane Width":<{w}s}{value(blk.RO.width):<{w}.3f}{f"{pyunits.get_units(blk.RO.width)}"}'
    )
    print(
        f'{f"Membrane Length":<{w}s}{value(blk.RO.length):<{w}.3f}{f"{pyunits.get_units(blk.RO.length)}"}'
    )
    print(
        f'{f"Inlet Crossflow Vel.":<{w}s}{value(blk.RO.feed_side.velocity[0, 0])*100:<{w}.3f}{f"cm/s"}'
    )
    print(
        f'{f"Outlet Crossflow Vel.":<{w}s}{value(blk.RO.feed_side.velocity[0, 1])*100:<{w}.3f}{f"cm/s"}'
    )
    print(
        f'{f"Perm Backpressure":<{w}s}{value(pyunits.convert(blk.RO.mixed_permeate[0].pressure, to_units=pyunits.bar)):<{w}.3f}{f"bar"}'
    )


def report_costing(blk, w=30):
    title = "Costing Report"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")
    print(f'{f"LCOW":<{w}s}{value(blk.LCOW):<{w}.3f}{f"{pyunits.get_units(blk.LCOW)}"}')
    print(f'{f"SEC":<{w}s}{value(blk.SEC):<{w}.3f}{f"{pyunits.get_units(blk.SEC)}"}')
    print(
        f'{f"CAPEX":<{w}s}{value(blk.total_capital_cost):<{w}.3f}{f"{pyunits.get_units(blk.total_capital_cost)}"}'
    )
    print(
        f'{f"OPEX":<{w}s}{value(blk.total_operating_cost):<{w}.3f}{f"{pyunits.get_units(blk.total_operating_cost)}"}'
    )


def report_n_stage_system(m, w=30):
    comp = m.fs.properties.solute_set.at(1)
    title = f"N-Stage System Report (N={len(m.fs.stage)})"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    perm_conc = value(m.fs.product.properties[0].conc_mass_phase_comp["Liq", comp])
    feed_conc = value(m.fs.feed.properties[0].conc_mass_phase_comp["Liq", comp])
    brine_conc = value(m.fs.disposal.properties[0].conc_mass_phase_comp["Liq", comp])
    print(f'{f"System Recovery":<{w}s}{value(m.fs.system_recovery)*100:<{w}.3f}{"%"}')
    print(f'{f"Feed {comp} Conc.":<{w}s}{feed_conc:<{w}.3f}{"g/L"}')
    print(f'{f"Perm {comp} Conc.":<{w}s}{perm_conc:<{w}.3f}{"g/L"}')
    print(f'{f"Brine {comp} Conc.":<{w}s}{brine_conc:<{w}.3f}{"g/L"}')

    for n, stage in m.fs.stage.items():

        title = f"STAGE {n}"
        side = int(((3 * w) - len(title)) / 2) - 1
        header = "(" * side + f" {title} " + ")" * side
        print(f"\n{header}\n")
        if stage.has_pump:
            report_pump(stage, w=w)
        report_ro(stage, w=w)

    if hasattr(m.fs, "costing"):
        report_costing(m.fs.costing, w=w)


def relax_bounds_for_low_salinity_waters(blk):

    comp = blk.config.property_package.solute_set.at(1)
    for x in list(blk.length_domain):
        blk.feed_side.velocity[0, x].setub(0.25)
        blk.feed_side.velocity[0, x].setlb(0.01)

    for i, v in blk.feed_side.K.items():
        v.setub(0.01)

    for i, v in blk.feed_side.friction_factor_darcy.items():
        v.setub(200)

    for i, v in blk.feed_side.cp_modulus.items():
        v.setlb(1e-5)
        v.setub(15)

    for (x, p, c), v in blk.recovery_mass_phase_comp.items():
        if c == comp:
            v.setlb(1e-9)
            v.setub(1e-1)

    for (t, x, p, c), v in blk.flux_mass_phase_comp.items():
        if c == comp:
            v.setlb(1e-9)
            v.setub(1e-1)
        if c == "H2O":
            v.setlb(1e-5)
            v.setub(0.999)

    for (x, p, c), v in blk.recovery_mass_phase_comp.items():
        if c == "H2O":
            v.setlb(1e-4)
            v.setub(0.999)


def overscale_ro(ro, props, full_scaling=True):
    comp = props.solute_set.at(1)
    h2o_scale = props._default_scaling_factors["flow_mass_phase_comp", ("Liq", "H2O")]
    comp_scale = props._default_scaling_factors["flow_mass_phase_comp", ("Liq", comp)]
    scales = {"H2O": h2o_scale, comp: comp_scale}

    iscale.set_scaling_factor(ro.area, 1e-2)
    iscale.constraint_scaling_transform(ro.eq_area, 1 / 100)
    iscale.set_scaling_factor(ro.width, 1)
    iscale.set_scaling_factor(ro.length, 1)

    for e in ro.feed_side.velocity:
        iscale.set_scaling_factor(ro.feed_side.velocity[e], 1)

    iscale.constraint_scaling_transform(ro.eq_area, 1 / 10)

    for temp_stream in [
        ro.eq_permeate_isothermal,
        ro.feed_side.eq_equal_temp_interface,
        ro.feed_side.eq_feed_isothermal,
        ro.eq_permeate_outlet_isothermal,
    ]:
        for e in temp_stream:
            iscale.constraint_scaling_transform(temp_stream[e], 1e-2)

    if full_scaling:
        for pressure_stream in [
            ro.eq_permeate_outlet_isobaric,
            ro.feed_side.eq_equal_pressure_interface,
        ]:
            for e in pressure_stream:
                iscale.constraint_scaling_transform(pressure_stream[e], 1e-5)
        for e in ro.eq_pressure_drop:
            iscale.constraint_scaling_transform(ro.eq_pressure_drop[e], 1e-4)
        for e in ro.feed_side.eq_K:
            iscale.constraint_scaling_transform(ro.feed_side.eq_K[e], 1e4)

        for e in ro.feed_side.eq_N_Sh_comp:
            iscale.constraint_scaling_transform(ro.feed_side.eq_N_Sh_comp[e], 1e-2)
        for e in ro.feed_side.eq_N_Re:
            iscale.constraint_scaling_transform(ro.feed_side.eq_N_Re[e], 1e2)

        for e in ro.feed_side.eq_friction_factor:
            iscale.constraint_scaling_transform(
                ro.feed_side.eq_friction_factor[e], 1e-2
            )
        for e in ro.feed_side.eq_dP_dx:
            iscale.constraint_scaling_transform(ro.feed_side.eq_dP_dx[e], 1e-3)

        for e in ro.feed_side.eq_equal_flow_vol_interface:
            iscale.constraint_scaling_transform(
                ro.feed_side.eq_equal_flow_vol_interface[e], 1e1
            )

        for e in ro.eq_mass_transfer_term:
            sf = scales[e[-1]]
            iscale.constraint_scaling_transform(ro.eq_mass_transfer_term[e], sf * 10)
        for e in ro.feed_side.mass_transfer_term:
            sf = scales[e[-1]]
            iscale.set_scaling_factor(ro.feed_side.mass_transfer_term[e], sf * 10)
        for e in ro.eq_mass_flux_equal_mass_transfer:
            sf = scales[e[-1]]
            iscale.constraint_scaling_transform(
                ro.eq_mass_flux_equal_mass_transfer[e], sf
            )
        for e in ro.eq_connect_mass_transfer:
            if e[-1] == "H2O":
                sf = sf * 10
            if e[-1] == comp:
                sf = sf * 100
            sf = scales[e[-1]]
            iscale.constraint_scaling_transform(ro.eq_connect_mass_transfer[e], sf)
        for e in ro.eq_recovery_mass_phase_comp:
            sf = scales[e[-1]]
            if e[-1] == "H2O":
                sf = sf * 10
            if e[-1] == comp:
                sf = sf * 100
            iscale.set_scaling_factor(ro.eq_recovery_mass_phase_comp[e], sf)
        for e in ro.eq_permeate_production:
            sf = scales[e[-1]]
            if e[-1] == "H2O":
                sf = sf * 10
            if e[-1] == comp:
                sf = sf * 100
            iscale.constraint_scaling_transform(ro.eq_permeate_production[e], sf)
        for e in ro.eq_flux_mass:
            sf = scales[e[-1]]
            if e[-1] == "H2O":
                sf = sf * 10
            if e[-1] == comp:
                sf = sf * 100
            iscale.constraint_scaling_transform(ro.eq_flux_mass[e], sf)


def solve(
    model=None,
    solver=None,
    tee=True,
    raise_on_failure=True,
    max_iter=1000,
):

    if solver is None:
        solver = get_solver()

    solver.options["max_iter"] = max_iter

    print("\n--------- SOLVING ---------\n")
    print(f"SOLVING with {degrees_of_freedom(model)} DOF")
    results = solver.solve(model, tee=tee)

    if check_optimal_termination(results):
        print("\n--------- OPTIMAL SOLVE!!! ---------\n")
        return results

    if raise_on_failure:
        print("\n--------- INFEASIBLE SOLVE!!! ---------\n")

        print("\n--------- CLOSE TO BOUNDS ---------\n")
        print_close_to_bounds(model)

        print("\n--------- INFEASIBLE BOUNDS ---------\n")
        print_infeasible_bounds(model)

        print("\n--------- INFEASIBLE CONSTRAINTS ---------\n")
        print_infeasible_constraints(model)
        msg = "The current configuration is infeasible. Please adjust the decision variables."

        raise RuntimeError(msg)
    return results
