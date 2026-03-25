from pyomo.environ import (
    ConcreteModel,
    value,
    Var,
    assert_optimal_termination,
    units as pyunits,
)
from watertap.flowsheets.ccro.utils.model_state_tool import ModelState
from idaes.core import UnitModelCostingBlock
from idaes.core.util.initialization import Constraint, solve_indexed_blocks

import idaes.core.util.scaling as iscale


def print_table(data, title=None, index_header="Index", precision=6, min_col_width=10):
    """
    Print a formatted table from a dictionary.

    Supports three input shapes:

    1. **Multi-row** ``{block_index: {key1: val1, key2: val2}, ...}``
       → horizontal table with one column per unique key and one row per block.

    2. **Single block_index** ``{"label": {key1: val1, key2: val2}}``
       → vertical two-column table (key | value), one row per key.

    3. **Flat dict** ``{key1: val1, key2: val2}``
       → vertical two-column table (key | value), one row per key.

    Args:
        data: dictionary (see above).
        title: optional string printed as a centered title above the table.
        index_header: header label for the leftmost column (default "Index").
              Used as the index column header in multi-row mode, or as the
              left-column header ("Metric" / "Parameter" etc.) in vertical mode.
        precision: number of decimal places for float values (default 6).
        min_col_width: minimum column width (default 10).

    Example::

        # Multi-row
        data = {
            1: {"Flow (L/min)": 3.14, "Pressure (Pa)": 101325},
            2: {"Flow (L/min)": 2.72, "Pressure (Pa)": 202650},
        }
        print_table(data, title="Results")

        # Flat / single-block  →  vertical layout
        print_table({"Recovery": 0.85, "LCOW": 1.23}, title="Summary")
    """
    if not data:
        print("(no data)")
        return

    # Format cell values
    def _fmt(val):
        if val is None:
            return ""
        if isinstance(val, float):
            return f"{val:.{precision}f}"
        return str(val)

    # ── Detect whether we should use vertical (key-value) layout ────
    # Case 1: flat dict – none of the values are dicts
    is_flat = all(not isinstance(v, dict) for v in data.values())
    # Case 2: single-entry nested dict
    is_single_block = (
        not is_flat and len(data) == 1 and isinstance(next(iter(data.values())), dict)
    )

    if is_flat or is_single_block:
        # Normalise to a single {key: value} mapping
        if is_single_block:
            kv = next(iter(data.values()))
            value_header = str(next(iter(data.keys())))
        else:
            kv = data
            value_header = "Value"

        key_width = max(min_col_width, len(index_header), *(len(str(k)) for k in kv))
        val_width = max(
            min_col_width, len(value_header), *(len(_fmt(v)) for v in kv.values())
        )
        total_width = key_width + val_width + 7  # "| " + " | " + " |"
        sep = "-" * total_width

        if title:
            padding = max(0, (total_width - len(title)) // 2 - 1)
            print(f"\n{'=' * padding} {title} {'=' * padding}")

        print(sep)
        print(f"| {index_header:<{key_width}s} | {value_header:<{val_width}s} |")
        print(sep)
        for key, val in kv.items():
            print(f"| {str(key):<{key_width}s} | {_fmt(val):<{val_width}s} |")
        print(sep)
        return

    # ── Multi-row horizontal table ──────────────────────────────────
    # Collect all unique column keys preserving first-seen order
    columns = []
    seen = set()
    for row_vals in data.values():
        if isinstance(row_vals, dict):
            for key in row_vals:
                if key not in seen:
                    columns.append(key)
                    seen.add(key)

    # Determine column widths (max of header length and all cell lengths)
    idx_width = max(min_col_width, len(index_header), *(len(str(k)) for k in data))
    col_widths = {}
    for col in columns:
        max_val_len = max(
            (len(_fmt(data[row].get(col))) for row in data),
            default=0,
        )
        col_widths[col] = max(min_col_width, len(col), max_val_len)

    # Build separator line
    total_width = idx_width + sum(col_widths[c] for c in columns) + len(columns) * 3 + 4
    sep = "-" * total_width

    # Print title
    if title:
        padding = max(0, (total_width - len(title)) // 2 - 1)
        print(f"\n{'=' * padding} {title} {'=' * padding}")

    # Print header
    print(sep)
    header = f"| {index_header:<{idx_width}s} |"
    for col in columns:
        header += f" {col:<{col_widths[col]}s} |"
    print(header)
    print(sep)

    # Print rows
    for idx, row_vals in data.items():
        row = f"| {str(idx):<{idx_width}s} |"
        for col in columns:
            val = _fmt(row_vals.get(col))
            row += f" {val:<{col_widths[col]}s} |"
        print(row)

    print(sep)


def copy_time_period_links(m_old, m_new, vars_to_copy=None):
    """
    Copy linking variables between time periods.
    Args:
        m_old: model of the old time period
        m_new: model of the new time period
        vars_to_copy: list of dictionaries with old and new model variable names to copy for example:
            {"old_model_var": 'fs.dead_volume.dead_volume.properties_out[0].mass_frac_phase_comp["Liq", "NaCl"]',
             "new_model_var": 'fs.dead_volume.delta_state.mass_frac_phase_comp[0, "Liq", "NaCl"]'}
    """

    if vars_to_copy is not None:
        for var in vars_to_copy:
            assert (
                m_new.find_component(var["new_model_var"]) is not None
            ), f"{var} not found in new model"
            assert (
                m_new.find_component(var["old_model_var"]) is not None
            ), f"{var} not found in new model"
            m_new.find_component(var["new_model_var"]).fix(
                m_old.find_component(var["old_model_var"]).value
            )


def copy_state(old_model, new_model):
    model_state = ModelState()
    model_state.get_model_state(old_model)
    model_state.set_model_state(new_model)


def register_costed_unit(
    mp,
    unit,
    costing_method_arguments={},
    register_electricity_cost=False,
    register_capital_cost=True,
    utilization_factor=1.0,
    power_expression=None,
):
    if register_electricity_cost and power_expression is None:
        print("Registered electricity cost for ", unit.name)
        lb = unit.work_mechanical[0.0].lb
        # set lower bound to 0 to avoid negative defined flow warning when lb is not >= 0
        unit.work_mechanical.setlb(0)
        utilization_factor = pyunits.convert(
            utilization_factor, to_units=pyunits.dimensionless
        )

        mp.costing.cost_flow(
            utilization_factor
            * pyunits.convert(unit.work_mechanical[0.0], to_units=pyunits.kW),
            "electricity",
        )
        # set lower bound back to its original value that was assigned to lb
        unit.work_mechanical.setlb(lb)
    if register_capital_cost:
        print("Registered capital cost for ", unit.name, costing_method_arguments)
        unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=mp.costing,
            costing_method_arguments=costing_method_arguments,
        )
    if register_electricity_cost and power_expression is not None:
        print(f"Registered power expression for {unit.name}")
        unit.total_power = Var(initialize=0, bounds=(0, None), units=pyunits.kW)
        iscale.set_scaling_factor(unit.total_power, 1e-3)
        unit.total_power_eq = Constraint(
            expr=unit.total_power
            == pyunits.convert(power_expression, to_units=pyunits.kW)
        )

        iscale.constraint_scaling_transform(unit.total_power_eq, 1e-3)
        # assert False        mp.costing.cost_flow(
        mp.costing.cost_flow(
            pyunits.convert(unit.total_power, to_units=pyunits.kW),
            "electricity",
        )
        # set lower bound back to its original value that was assigned to lb


def set_pump_operating_pressure(unit, osmotic_pressure, configuration_options):
    """
    Set pump operating pressure based on configuration option.
    Args:
        unit: pump unit model
        osm_pressure: osmotic pressure at the pump inlet
        configuration_options: CCROConfiguration object
    """
    if configuration_options["p1_pressure_start"] == "osmotic_pressure":
        set_pressure = (
            value(osmotic_pressure)
            * value(configuration_options["osmotic_overpressure"])
            + 2e5
        )  # Pa
        print(
            f"Setting {unit.name} outlet pressure based on osmotic pressure of: {value(osmotic_pressure)/1e5} bar, pressure is set to {(set_pressure)/1e5} bar"
        )
        unit.control_volume.properties_out[0].pressure.fix(set_pressure)
    else:
        unit.control_volume.properties_out[0].pressure.fix(
            configuration_options["p1_pressure_start"]
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


def check_jac(m, print_extreme_jacobian_values=True):
    jac, jac_scaled, nlp = iscale.constraint_autoscale_large_jac(m, min_scale=1e-8)
    try:
        cond_number = iscale.jacobian_cond(m, jac=jac_scaled) / 1e10
        print("--------------------------")
        print("COND NUMBER:", cond_number)
    except:
        print("Cond number failed")
        cond_number = None
    if print_extreme_jacobian_values:
        print("--------------------------")
        print("Extreme Jacobian entries:")
        extreme_entries = iscale.extreme_jacobian_entries(
            m, jac=jac_scaled, nlp=nlp, zero=1e-20, large=100
        )
        for val, var, con in extreme_entries:
            print(val, var.name, con.name)
        print("--------------------------")
        print("Extreme Jacobian columns:")
        extreme_cols = iscale.extreme_jacobian_columns(
            m, jac=jac_scaled, nlp=nlp, small=1e-3
        )
        for val, var in extreme_cols:
            print(val, var.name)
        print("------------------------")
        print("Extreme Jacobian rows:")
        extreme_rows = iscale.extreme_jacobian_rows(
            m, jac=jac_scaled, nlp=nlp, small=1e-3
        )
        for val, con in extreme_rows:
            print(val, con.name)

    for i in iscale.list_unscaled_variables(m):
        print("Var with no scale:", i)
    for i in iscale.list_unscaled_constraints(m):
        print("Constraint with no scale:", i)
    for var, scale in iscale.badly_scaled_var_generator(m):
        print(
            "Badly scaled variable:",
            var.name,
            var.value,
            iscale.get_scaling_factor(var),
        )
    return cond_number
