from pyomo.environ import (
    ConcreteModel,
    value,
    assert_optimal_termination,
    units as pyunits,
)
from watertap.flowsheets.ccro.multiperiod.model_state_tool import ModelState
from idaes.core import UnitModelCostingBlock
from idaes.core.util.initialization import solve_indexed_blocks


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
    register_electricity_cost=True,
    register_capital_cost=True,
    utilization_factor=1.0,
):
    if register_electricity_cost:
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
        unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=mp.costing,
            costing_method_arguments=costing_method_arguments,
        )

        # assert False


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

    deltaP = blk.pump.deltaP[0]
    work = blk.pump.work_mechanical[0]

    print(
        f'{f"Pressure Change (bar)":<{w}s}{value(pyunits.convert(deltaP,to_units=pyunits.bar)):<{w}.3f}{"bar"}'
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
        f'{f"Inlet Flow":<{w}s}{value(blk.feed.properties[0].flow_vol_phase["Liq"]):<{w}.4f}{"m3/s"}'
    )
    print(
        f'{f"Inlet Conc.":<{w}s}{value(pyunits.convert(blk.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"], to_units=pyunits.mg / pyunits.liter)):<{w}.3f}{"mg/L"}'
    )
    print(
        f'{f"Brine Flow":<{w}s}{value(blk.disposal.properties[0].flow_vol_phase["Liq"]):<{w}.4f}{"m3/s"}'
    )
    print(
        f'{f"Product Flow":<{w}s}{value(blk.product.properties[0].flow_vol_phase["Liq"]):<{w}.6f}{"m3/s"}'
    )
    print(
        f'{f"Recovery":<{w}s}{value(blk.RO.recovery_vol_phase[0, "Liq"])*100:<{w}.3f}{"%"}'
    )
    print(
        f'{f"Rejection":<{w}s}{value(blk.RO.rejection_phase_comp[0, "Liq", "NaCl"])*100:<{w}.3f}{"%"}'
    )
    print(
        f'{f"Perm Conc":<{w}s}{value(pyunits.convert(blk.RO.mixed_permeate[0].conc_mass_phase_comp["Liq", "NaCl"], to_units=pyunits.mg / pyunits.liter)):<{w}.3f}{f"mg/L"}'
    )
    print(
        f'{f"Brine Conc":<{w}s}{value(pyunits.convert(blk.RO.feed_side.properties[0, 1].conc_mass_phase_comp["Liq", "NaCl"], to_units=pyunits.mg / pyunits.liter)):<{w}.3f}{f"mg/L"}'
    )
    print(
        f'{f"∆P":<{w}s}{value(pyunits.convert(blk.RO.deltaP[0], to_units=pyunits.bar)):<{w}.3f}{f"bar"}'
    )
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
    blk.feed_side.cp_modulus.setub(5)
    for e in blk.feed_side.K:
        blk.feed_side.K[e].setub(0.01)
        # blk.feed_side.K[e].setlb(1e-7)

    for e in blk.feed_side.cp_modulus:
        blk.feed_side.cp_modulus[e].setlb(1e-5)

    for e in blk.recovery_mass_phase_comp:
        if e[-1] == "NaCl":
            blk.recovery_mass_phase_comp[e].setlb(1e-9)
            blk.recovery_mass_phase_comp[e].setub(1e-1)

    for e in blk.flux_mass_phase_comp:
        if e[-1] == "NaCl":
            blk.flux_mass_phase_comp[e].setlb(1e-9)
            blk.flux_mass_phase_comp[e].setub(1e-1)

    for e in blk.recovery_mass_phase_comp:
        if e[-1] == "H2O":
            blk.recovery_mass_phase_comp[e].setlb(1e-4)
            blk.recovery_mass_phase_comp[e].setub(0.999)

    for e in blk.flux_mass_phase_comp:
        if e[-1] == "H2O":
            blk.flux_mass_phase_comp[e].setlb(1e-5)
            blk.flux_mass_phase_comp[e].setub(0.999)
