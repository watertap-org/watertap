from pyomo.environ import (
    ConcreteModel,
    value,
    assert_optimal_termination,
    check_optimal_termination,
    units as pyunits,
)
from idaes.core.util.initialization import Constraint, solve_indexed_blocks
import idaes.core.util.scaling as iscale
from watertap.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.core.util.model_diagnostics.infeasible import (
    print_close_to_bounds,
    print_infeasible_bounds,
    print_infeasible_constraints,
)


def calculate_operating_pressure(
    feed_state_block=None,
    over_pressure=0.15,
    water_recovery=0.5,
    NaCl_passage=0.01,
    solver=None,
):
    t = ConcreteModel()
    prop = feed_state_block.config.parameters
    t.brine = prop.build_state_block([0])

    t.brine[0].flow_mass_phase_comp["Liq", "H2O"].fix(
        value(feed_state_block.flow_mass_phase_comp["Liq", "H2O"])
        * (1 - water_recovery)
    )
    t.brine[0].flow_mass_phase_comp["Liq", "NaCl"].fix(
        value(feed_state_block.flow_mass_phase_comp["Liq", "NaCl"]) * (1 - NaCl_passage)
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
        f'{f"Inlet Pressure (bar)":<{w}s}{value(pyunits.convert(pin,to_units=pyunits.bar)):<{w}.3f}{"bar"}'
    )
    print(
        f'{f"Pressure Change (bar)":<{w}s}{value(pyunits.convert(deltaP,to_units=pyunits.bar)):<{w}.3f}{"bar"}'
    )
    print(
        f'{f"Outlet Pressure (bar)":<{w}s}{value(pyunits.convert(pout,to_units=pyunits.bar)):<{w}.3f}{"bar"}'
    )
    print(
        f'{f"Power (kW)":<{w}s}{value(pyunits.convert(work, to_units=pyunits.kW)):<{w}.3f}{"kW"}'
    )
    print(f'{f"Efficiency (-)":<{w}s}{value(blk.pump.efficiency_pump[0]):<{w}.3f}{"-"}')


def report_ro(blk, w=30):
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
        f'{f"Inlet Conc.":<{w}s}{value(blk.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"]):<{w}.3f}{"g/L"}'
    )
    print(
        f'{f"Perm Conc.":<{w}s}{value(blk.product.properties[0].conc_mass_phase_comp["Liq", "NaCl"]):<{w}.3f}{"g/L"}'
    )
    print(
        f'{f"Brine Conc.":<{w}s}{value(blk.disposal.properties[0].conc_mass_phase_comp["Liq", "NaCl"]):<{w}.3f}{"g/L"}'
    )
    print(
        f'{f"Rejection":<{w}s}{value(blk.RO.rejection_phase_comp[0, "Liq", "NaCl"])*100:<{w}.3f}{"%"}'
    )
    print(
        f'{f"deltaP":<{w}s}{value(pyunits.convert(blk.RO.deltaP[0], to_units=pyunits.bar)):<{w}.3f}{f"bar"}'
    )
    acomp = blk.RO.A_comp[0, "H2O"]
    print(
        f'{f"Water Perm (A)":<{w}s}{value(pyunits.convert(acomp, to_units=pyunits.liter / pyunits.m**2 / pyunits.hour / pyunits.bar)):<{w}.3f}{f"LMH/bar"}'
    )
    print(f'{f"Water Perm (A)":<{w}s}{value(acomp):<{w}.3e}{f"m/s/Pa"}')
    bcomp = blk.RO.B_comp[0, "NaCl"]
    print(
        f'{f"Salt Perm (B)":<{w}s}{value(pyunits.convert(bcomp, to_units=pyunits.liter / pyunits.m**2 / pyunits.hour)):<{w}.3f}{f"LMH"}'
    )
    print(f'{f"Salt Perm (B)":<{w}s}{value(bcomp):<{w}.3e}{f"m/s"}')
    print(f'{f"Total Flux":<{w}s}{value(blk.flux):<{w}.3f}{f"LMH"}')

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
    # for (_, i, _), x in blk.RO.feed_side.cp_modulus.items():
    #     # print(i, x)
    #     print(
    #         f'{f"CP Modulus @ {i}":<{w}s}{value(x):<{w}.3f}{"-"}'
    #     )
    # for (_, i, _, c), x in blk.RO.flux_mass_phase_comp.items():
    #     # print(i, x)
    #     if c != "H2O":
    #         continue
    #     print(
    #         f'{f"Flux {c} @ {i}":<{w}s}{value(x):<{w}.3f}{"-"}'
    #     )
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
    print(
        f'{f"Levelized Cost of Water":<{w}s}{value(blk.LCOW):<{w}.3f}{f"${pyunits.get_units(blk.LCOW)}"}'
    )
    print(
        f'{f"Specific Energy Consumption":<{w}s}{value(blk.SEC):<{w}.3f}{f"{pyunits.get_units(blk.SEC)}"}'
    )


def report_n_stage_system(m, w=30):

    title = f"N-Stage System Report (N={len(m.fs.stage)})"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    perm_conc = value(m.fs.product.properties[0].conc_mass_phase_comp["Liq", "NaCl"])
    feed_conc = value(m.fs.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"])
    brine_conc = value(m.fs.disposal.properties[0].conc_mass_phase_comp["Liq", "NaCl"])
    print(f'{f"System Recovery":<{w}s}{value(m.fs.system_recovery)*100:<{w}.3f}{"%"}')
    print(f'{f"Feed Conc":<{w}s}{feed_conc:<{w}.3f}{"g/L"}')
    print(f'{f"Perm Conc":<{w}s}{perm_conc:<{w}.3f}{"g/L"}')
    print(f'{f"Brine Conc":<{w}s}{brine_conc:<{w}.3f}{"g/L"}')

    for n, stage in m.fs.stage.items():

        title = f"STAGE {n}"
        side = int(((3 * w) - len(title)) / 2) - 1
        header = "(" * side + f" {title} " + ")" * side
        print(f"\n{header}\n")
        if stage.add_pump:
            report_pump(stage, w=w)
        report_ro(stage, w=w)

    if hasattr(m.fs, "costing"):
        report_costing(m.fs.costing, w=w)



def relax_bounds_for_low_salinity_waters(blk):

    # blk.feed_side.velocity[0, 0].setub(0.3)
    # blk.feed_side.velocity[0, 0].setlb(0.05)

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
        if c == "NaCl":
            v.setlb(1e-9)
            v.setub(1e-1)

    for (t, x, p, c), v in blk.flux_mass_phase_comp.items():
        if c == "NaCl":
            v.setlb(1e-9)
            v.setub(1e-1)
        if c == "H2O":
            v.setlb(1e-5)
            v.setub(0.999)

    for (x, p, c), v in blk.recovery_mass_phase_comp.items():
        if c == "H2O":
            v.setlb(1e-4)
            v.setub(0.999)


def overscale_ro(blk, props, full_scaling=True):

    h2o_scale = props._default_scaling_factors["flow_mass_phase_comp", ("Liq", "H2O")]
    NaCl_scale = props._default_scaling_factors["flow_mass_phase_comp", ("Liq", "NaCl")]
    scales = {"H2O": h2o_scale, "NaCl": NaCl_scale}

    iscale.set_scaling_factor(blk.area, 1e-2)
    iscale.constraint_scaling_transform(blk.eq_area, 1 / 100)
    iscale.set_scaling_factor(blk.width, 1)
    iscale.set_scaling_factor(blk.length, 1)
    stage = ro
    for e in blk.feed_side.velocity:
        iscale.set_scaling_factor(blk.feed_side.velocity[e], 1)

    iscale.constraint_scaling_transform(blk.eq_area, 1 / 10)

    for temp_stream in [
        blk.eq_permeate_isothermal,
        blk.feed_side.eq_equal_temp_interface,
        blk.feed_side.eq_feed_isothermal,
        blk.eq_permeate_outlet_isothermal,
    ]:
        for e in temp_stream:
            iscale.constraint_scaling_transform(temp_stream[e], 1e-2)
    if full_scaling:
        for pressure_stream in [
            blk.eq_permeate_outlet_isobaric,
            blk.feed_side.eq_equal_pressure_interface,
        ]:
            for e in pressure_stream:
                iscale.constraint_scaling_transform(pressure_stream[e], 1e-5)
        for e in blk.eq_pressure_drop:
            iscale.constraint_scaling_transform(blk.eq_pressure_drop[e], 1e-4)
        for e in blk.feed_side.eq_K:
            iscale.constraint_scaling_transform(blk.feed_side.eq_K[e], 1e4)

        for e in blk.feed_side.eq_N_Sh_comp:
            iscale.constraint_scaling_transform(blk.feed_side.eq_N_Sh_comp[e], 1e-2)
        for e in blk.feed_side.eq_N_Re:
            iscale.constraint_scaling_transform(blk.feed_side.eq_N_Re[e], 1e2)

        for e in blk.feed_side.eq_friction_factor:
            iscale.constraint_scaling_transform(
                blk.feed_side.eq_friction_factor[e], 1e-2
            )
        for e in blk.feed_side.eq_dP_dx:
            iscale.constraint_scaling_transform(blk.feed_side.eq_dP_dx[e], 1e-3)

        for e in blk.feed_side.eq_equal_flow_vol_interface:
            iscale.constraint_scaling_transform(
                blk.feed_side.eq_equal_flow_vol_interface[e], 1e1
            )

        for e in blk.eq_mass_transfer_term:
            sf = scales[e[-1]]
            iscale.constraint_scaling_transform(blk.eq_mass_transfer_term[e], sf * 10)
        for e in blk.feed_side.mass_transfer_term:
            sf = scales[e[-1]]
            iscale.set_scaling_factor(blk.feed_side.mass_transfer_term[e], sf * 10)
        for e in blk.eq_mass_flux_equal_mass_transfer:
            sf = scales[e[-1]]
            iscale.constraint_scaling_transform(
                blk.eq_mass_flux_equal_mass_transfer[e], sf
            )
        for e in blk.eq_connect_mass_transfer:
            if e[-1] == "H2O":
                sf = sf * 10
            if e[-1] == "NaCl":
                sf = sf * 100
            sf = scales[e[-1]]
            iscale.constraint_scaling_transform(blk.eq_connect_mass_transfer[e], sf)
        for e in blk.eq_recovery_mass_phase_comp:
            sf = scales[e[-1]]
            if e[-1] == "H2O":
                sf = sf * 10
            if e[-1] == "NaCl":
                sf = sf * 100
            iscale.set_scaling_factor(blk.eq_recovery_mass_phase_comp[e], sf)
        for e in blk.eq_permeate_production:
            sf = scales[e[-1]]
            if e[-1] == "H2O":
                sf = sf * 10
            if e[-1] == "NaCl":
                sf = sf * 100
            iscale.constraint_scaling_transform(blk.eq_permeate_production[e], sf)
        for e in blk.eq_flux_mass:
            sf = scales[e[-1]]
            if e[-1] == "H2O":
                sf = sf * 10
            if e[-1] == "NaCl":
                sf = sf * 100
            iscale.constraint_scaling_transform(blk.eq_flux_mass[e], sf)


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
    # print_close_to_bounds(model)
    if raise_on_failure:
        print("\n--------- INFEASIBLE SOLVE!!! ---------\n")

        print("\n--------- CLOSE TO BOUNDS ---------\n")
        print_close_to_bounds(model)

        print("\n--------- INFEASIBLE BOUNDS ---------\n")
        print_infeasible_bounds(model)

        print("\n--------- INFEASIBLE CONSTRAINTS ---------\n")
        print_infeasible_constraints(model)
        msg = (
            "The current configuration is infeasible. Please adjust the decision variables."
        )

        raise RuntimeError(msg)
    return results
