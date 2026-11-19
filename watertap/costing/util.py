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

import functools
from math import isclose, ceil
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors
import pyomo.environ as pyo
from watertap.custom_exceptions import FrozenPipes


def register_costing_parameter_block(build_rule, parameter_block_name):
    """
    Decorator to create a global-level parameter block on the costing package.

    This decorator is meant to be utilized on a costing method or function in order
    to register parameters (including flows specific to a unit model) which should
    have global-scope within a costing package. In this way, all instances of a unit
    model can share the same global parameters, such as baseline capital costs.

    The decorator takes a `build_rule`, which is a rule for constructing a Pyomo
    `Block`, and a `parameter_block_name`, which determines the name of the constructed
    block. The decorator then puts the new Block onto the costing package. For a unit
    model costing block named `blk`, this is symmatically equivalent to:

    .. testcode::
        :skipif: True

        blk.costing_package.parameter_block_name = Block(rule=build_rule)

    However, this decorator contains a few additional safeguards.

        1. If the wrapped unit model costing function is called more than once,
           the associated parameter block will only be added one time.

        2. Similarly, if a component named `parameter_block_name` already exists
           on the `costing_package`, this decorator will raise an error if the
           parameter block was built with a different `build_rule`. This prevents
           inadvertent name clashes between different parameter blocks.

        3. The decorator fixes any Pyomo Vars it finds on the constructed parameter
           block.

    Use of this decorator allows for defining unit-level costing parameters in the same
    module as the unit-level costing method. It also ensures that unit-level costing
    parameters are only added to the costing package if and only if that particular
    unit model is costed.

    Args:
        build_rule (callable): A function defining a "build rule" for a Pyomo `Block`.
        parameter_block_name (str): The name for the built Block on the costing package.
    """

    def register_costing_parameter_block_decorator(func):
        @functools.wraps(func)
        def add_costing_parameter_block(blk, *args, **kwargs):
            parameter_block = blk.costing_package.component(parameter_block_name)
            if parameter_block is None:
                parameter_block = pyo.Block(rule=build_rule)
                blk.costing_package.add_component(parameter_block_name, parameter_block)
                # fix the parameters in case the build_rule did not
                parameter_block.fix_all_vars()
            elif parameter_block._rule is None or not hasattr(
                parameter_block._rule, "_fcn"
            ):
                raise RuntimeError(
                    "Use the register_costing_parameter_block decorator for specifying"
                    "costing-package-level parameters"
                )
            elif parameter_block._rule._fcn is not build_rule:
                other_rule = parameter_block._rule._fcn
                raise RuntimeError(
                    "Attempting to add identically named costing parameter blocks with "
                    "different build rules to the costing package "
                    f"{blk.costing_package}. Parameter block named "
                    f"{parameter_block_name} was previously built by function "
                    f"{other_rule.__name__} from module {other_rule.__module__}"
                )
            # else parameter_block was constructed by build_rule previously
            return func(blk, *args, **kwargs)

        return add_costing_parameter_block

    return register_costing_parameter_block_decorator


def make_capital_cost_var(blk):
    """
    Generic function for adding a Var named `capital_cost`
    to the unit model costing block.

    Creates a pyomo Var named `capital_cost` with units
    `blk.costing_package.base_currency`.
    """
    blk.capital_cost = pyo.Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Unit capital cost",
    )


def make_fixed_operating_cost_var(blk):
    """
    Generic function for adding a Var named `fixed_operating_cost`
    to the unit model costing block.

    Creates a pyomo Var named `fixed_operating_cost` with units
    `blk.costing_package.base_currency / blk.costing_package.base_period`.
    """
    blk.fixed_operating_cost = pyo.Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency / blk.costing_package.base_period,
        doc="Unit fixed operating cost",
    )


def cost_membrane(blk, membrane_cost, factor_membrane_replacement):
    """
    Generic function for costing a membrane. Assumes the unit_model
    has an `area` variable or parameter.

    Args:
        membrane_cost: The cost of the membrane in currency per area
        factor_membrane_replacement: Membrane replacement factor
            (fraction of membrane replaced/year)
    """

    make_capital_cost_var(blk)
    make_fixed_operating_cost_var(blk)
    blk.membrane_cost = pyo.Expression(expr=membrane_cost)
    blk.factor_membrane_replacement = pyo.Expression(expr=factor_membrane_replacement)

    blk.costing_package.add_cost_factor(blk, "TIC")
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == blk.cost_factor
        * pyo.units.convert(
            blk.membrane_cost * blk.unit_model.area,
            to_units=blk.costing_package.base_currency,
        )
    )
    blk.fixed_operating_cost_constraint = pyo.Constraint(
        expr=blk.fixed_operating_cost
        == pyo.units.convert(
            blk.factor_membrane_replacement * blk.membrane_cost * blk.unit_model.area,
            to_units=blk.costing_package.base_currency
            / blk.costing_package.base_period,
        )
    )


def cost_rectifier(blk, power=100 * pyo.units.kW, ac_dc_conversion_efficiency=0.90):
    """
    Method to cost rectifiers for electrified process units that require direct current which must be converted
    from an alternating current source. Note that this should be used solely for units that require the conversion,
    and should not be used universally for electricity requirements.
    Assumes the unit_model has a `power` variable or parameter.

    Args:
        ac_dc_conversion_efficiency: Efficiency of the conversion from AC to DC current
    """

    # create variables on cost block
    make_capital_cost_var(blk)
    blk.ac_dc_conversion_efficiency = pyo.Expression(
        expr=ac_dc_conversion_efficiency,
        doc="fixing unit model vairable for upscaling required power considering "
        "the efficiency of converting alternating to direct current",
    )
    blk.ac_power = pyo.Var(
        initialize=100,
        domain=pyo.NonNegativeReals,
        units=pyo.units.kW,
        doc="Unit AC power",
    )

    # use unit.power variable in conversion with efficiency
    blk.power_conversion = pyo.Constraint(
        expr=blk.ac_power * blk.ac_dc_conversion_efficiency
        == pyo.units.convert(power, to_units=pyo.units.kW)
    )

    # USD_2021 embedded in equation
    rectifier_cost_coeff = {0: 508.6, 1: 2810}
    blk.rectifier_cost_coeff = pyo.Var(
        rectifier_cost_coeff.keys(),
        initialize=rectifier_cost_coeff,
        units=pyo.units.dimensionless,
        doc="Rectifier cost coefficients",
    )
    blk.capital_cost_rectifier = pyo.Var(
        initialize=100, units=blk.costing_package.base_currency
    )

    # refix variables to appropriate costing parameters
    for index, var in blk.rectifier_cost_coeff.items():
        var.fix(rectifier_cost_coeff[index])

    # calculate capital cost
    blk.capital_cost_rectifier_constraint = pyo.Constraint(
        expr=blk.capital_cost_rectifier
        == pyo.units.convert(
            pyo.units.USD_2021
            * (
                blk.rectifier_cost_coeff[1]
                + (blk.rectifier_cost_coeff[0] * (blk.ac_power * pyo.units.kW**-1))
            ),
            to_units=blk.costing_package.base_currency,
        )
    )

    # cost electricity flow
    blk.costing_package.cost_flow(blk.ac_power, "electricity")


def cost_by_flow_volume(blk, flow_cost, flow_to_cost):
    """
    Generic function for costing by flow volume.

    Args:
        flow_cost: The cost of the device in [currency]/([volume]/[time])
        flow_to_cost: The flow costed in [volume]/[time]
    """
    make_capital_cost_var(blk)
    blk.flow_cost = pyo.Expression(expr=flow_cost)
    blk.costing_package.add_cost_factor(blk, "TIC")
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == blk.cost_factor
        * pyo.units.convert(
            blk.flow_cost * flow_to_cost, to_units=blk.costing_package.base_currency
        )
    )


def cost_steam_flow(
    costing_package,
    steam_cost_type="steam",
    steam_volume_flow=None,
    steam_heat_duty=None,
    steam_mass_flow=None,
    steam_pressure=None,
):
    """
    Generic function for costing a steam flow.

    Args:
        steam_cost_type: name of the steam cost to use for costing flow, if not available will be created
        steam_volume_flow: the steam flow to cost in [volume]/[time]
        steam_mass_flow (optional): the steam flow to cost in [mass]/[time]
        steam_pressure (needed if mass flow is provided): the pressure of the steam flow in [pressure]
    """
    if steam_mass_flow is None:
        if steam_mass_flow is not None and steam_heat_duty is not None:
            raise ValueError(
                "If steam_volume_flow is not provided, either steam mass_flow or steam_heat_duty should be provided"
            )
        assert (
            steam_pressure is not None
        ), "steam pressure must be provided if steam mass flow or heat duty is provided"
        dh_vap, sp_vol = _compute_steam_properties(steam_pressure)
        if steam_volume_flow is not None:
            steam_mass_flow = steam_volume_flow / sp_vol
        elif steam_heat_duty is not None:
            steam_mass_flow = steam_heat_duty / dh_vap

    if costing_package.component(steam_cost_type) is None:
        costing_package.add_component(
            steam_cost_type,
            pyo.Var(initialize=0.004, units=pyo.units.USD_2018 / (pyo.units.kg)),
        )
        costing_package.find_component(steam_cost_type).fix()
        costing_package.register_flow_type(
            steam_cost_type, costing_package.find_component(steam_cost_type)
        )
    costing_package.cost_flow(
        pyo.units.convert(
            (steam_mass_flow),
            to_units=pyo.units.kg / pyo.units.s,
        ),
        steam_cost_type,
    )


def _compute_steam_properties(steam_pressure):
    """
    Function for computing saturated steam properties for thermal heating estimation.

    Args:
        steam_pressure:   Steam gauge pressure in bar

    Out:
        Steam thermal capacity (latent heat of condensation * density) in kJ/m3
    """
    pressure_sat = steam_pressure
    # 1. Compute saturation temperature of steam: computed from El-Dessouky expression
    tsat_constants = [
        42.6776 * pyo.units.K,
        -3892.7 * pyo.units.K,
        1000 * pyo.units.kPa,
        -9.48654 * pyo.units.dimensionless,
    ]
    psat = (
        pyo.units.convert(pressure_sat, to_units=pyo.units.kPa)
        + 101.325 * pyo.units.kPa
    )
    temperature_sat = tsat_constants[0] + tsat_constants[1] / (
        pyo.log(psat / tsat_constants[2]) + tsat_constants[3]
    )

    # 2. Compute latent heat of condensation/vaporization: computed from Sharqawy expression
    t = temperature_sat - 273.15 * pyo.units.K
    enth_mass_units = pyo.units.J / pyo.units.kg
    t_inv_units = pyo.units.K**-1
    dh_constants = [
        2.501e6 * enth_mass_units,
        -2.369e3 * enth_mass_units * t_inv_units**1,
        2.678e-1 * enth_mass_units * t_inv_units**2,
        -8.103e-3 * enth_mass_units * t_inv_units**3,
        -2.079e-5 * enth_mass_units * t_inv_units**4,
    ]
    dh_vap = (
        dh_constants[0]
        + dh_constants[1] * t
        + dh_constants[2] * t**2
        + dh_constants[3] * t**3
        + dh_constants[4] * t**4
    )
    dh_vap = pyo.units.convert(dh_vap, to_units=pyo.units.kJ / pyo.units.kg)

    # 3. Compute specific volume: computed from Affandi expression (Eq 5)
    t_critical = 647.096 * pyo.units.K
    t_red = temperature_sat / t_critical  # Reduced temperature
    sp_vol_constants = [
        -7.75883 * pyo.units.dimensionless,
        3.23753 * pyo.units.dimensionless,
        2.05755 * pyo.units.dimensionless,
        -0.06052 * pyo.units.dimensionless,
        0.00529 * pyo.units.dimensionless,
    ]
    log_sp_vol = (
        sp_vol_constants[0]
        + sp_vol_constants[1] * (pyo.log(1 / t_red)) ** 0.4
        + sp_vol_constants[2] / (t_red**2)
        + sp_vol_constants[3] / (t_red**4)
        + sp_vol_constants[4] / (t_red**5)
    )
    sp_vol = pyo.exp(log_sp_vol) * pyo.units.m**3 / pyo.units.kg

    # 4. Return specific energy: density and density
    return dh_vap, sp_vol


def plot_LCOW_breakdown(
    blk,
    log,
    lcow_name="LCOW",
    fs_name="fs",
    by="component",
    separate_flows=False,
    relative=False,
    colormap=None,
    save_as=None,
    close_fig=False,
    tol=1e-5,
    plot_tol=0,
    component_labels={},
    **kwargs,
):
    """
    This method creates a bar chart showing the breakdown of the LCOW calculation
    into its contributing components. The breakdown can be done by individual unit
    models or by unit model type, and flows can be shown separately or included as
    part of variable OPEX.

    Args:
        blk: a WaterTAP costing block
        log: logger object for logging messages
        lcow_name: name of the LCOW Expression to plot (default: "LCOW")
        fs_name: name of the flowsheet (default: "fs")
        by: breakdown by individual "component" unit or "aggregate" unit class (default: "component")
        separate_flows: whether to separate material/energy flows from OPEX portion (default: False)
        relative: whether to plot as fraction of total LCOW (default: False)
        colormap: colormap to use for the plot (default: None; a custom colormap of tab20b + tab20c is created)
        save_as: filename to save the plot (default: None)
        close_fig: whether to close the figure after creation (default: False)
        tol: tolerance for LCOW calculation check (default: 1e-5)
        plot_tol: minimum absolute value of SEC contribution to include in the plot (default: 0)
        component_labels: dictionary mapping flowsheet unit names to custom unit names (default: {})

    Returns:
        fig, ax: matplotlib figure and axis objects for the plot
    """

    def _add_opex(v):
        if pyo.value(v) > 0:
            pos_opex[k.replace(f"{fs_name}.", "")] += pyo.value(v)
        else:
            neg_opex[k.replace(f"{fs_name}.", "")] += pyo.value(v)

    def _get_label(u):
        if u not in component_labels:
            return u.split(".")[-1]
        return component_labels[u]

    hatch_dict = {
        "CAPEX": "",
        "OPEX": "\\\\",
    }
    flows_hatch = ".."

    if blk.find_component(lcow_name) is None:
        msg = f"{blk.name} does not have a component named {lcow_name}."
        msg += " Use the add_LCOW method to add an LCOW Expression to the costing"
        msg += " block before calling plot_LCOW_breakdown."
        raise RuntimeError(msg)

    if by not in ("component", "aggregate"):
        raise ValueError("'by' argument must be either 'component' or 'aggregate'.")

    lcow_actual = pyo.value(blk.find_component(lcow_name))

    capex_comp_exprs = [
        f"{lcow_name}_{by}_direct_capex",
        f"{lcow_name}_{by}_indirect_capex",
    ]

    # CAPEX portion of LCOW
    capex = defaultdict(float)

    # OPEX portions of LCOW
    pos_opex = defaultdict(float)
    neg_opex = defaultdict(float)
    flows = defaultdict(float)

    # CAPEX = direct + indirect
    for ce in capex_comp_exprs:
        expr = blk.find_component(ce)
        for k, v in expr.items():
            capex[k.replace(f"{fs_name}.", "")] += pyo.value(v)

    # OPEX will always include fixed OPEX,
    # but may or may not include variable OPEX depending on `separate_flows`
    for k, v in blk.find_component(f"{lcow_name}_{by}_fixed_opex").items():
        _add_opex(v)

    # All unique units contributing to LCOW
    all_units = set([k for _d in [pos_opex, neg_opex, capex] for k in _d])

    # NOTE: Nov. 2026: There are no WaterTAP unit models that have
    # a variable_operating_cost that is not a flow cost.
    # If that were to change in the future, it will break the LCOW
    # accounting that is done here.
    if separate_flows:
        # Aggregate flows are plotted separately from variable OPEX
        for k, v in blk.find_component(f"{lcow_name}_aggregate_variable_opex").items():
            if k not in blk.used_flows:
                continue
            flows[k] += pyo.value(v)
    else:
        # Aggregate flows are included as part of variable OPEX
        for k, v in blk.find_component(f"{lcow_name}_{by}_variable_opex").items():
            _add_opex(v)

    # All unique components contributing to LCOW
    all_components = list(flows.keys()) + list(all_units)

    if colormap is None:
        # Create a custom colormap with 40 unique colors
        cmap = plt.get_cmap("tab20b").colors + plt.get_cmap("tab20c").colors
        cmap = mcolors.ListedColormap(cmap)
    else:
        cmap = plt.get_cmap(colormap)

    colors = [cmap(i / len(all_components)) for i in range(len(all_components))]
    color_dict = dict(zip(all_components, colors))

    # Start plotting
    fig, ax = plt.subplots()
    x = 0  # location of bar
    width = 0.5  # bar width
    xmin = x - width
    xmax = x + width

    # Starting bottom needs to be at zero
    bottom = 0
    d = lcow_actual * 1e-2 if relative else 1

    # Create legend
    handles = list()
    labels = list()

    # For checking LCOW calculation
    lcow_check = 0
    # Sum of negative OPEX contributions
    lcow_neg = sum(neg_opex.values())
    # Track units to inclue in legend
    leg_units = list()

    # Starting from the bottom:
    # add negative OPEX, then flows, then positive OPEX/CAPEX for each unit
    if lcow_neg < 0:
        log.info(
            "Cannot plot relative values with negative LCOW contributions; setting 'relative' to False."
        )
        relative = False
        d = 1
        # Set bottom to the sum of negative OPEX contributions
        bottom = sum(v / d for v in neg_opex.values() if abs(v) > plot_tol)
        for u, v in neg_opex.items():
            neg_o = v / d
            if abs(v) > plot_tol:
                ax.bar(
                    [x],
                    [abs(neg_o)],
                    bottom=bottom,
                    hatch=hatch_dict["OPEX"],
                    color=color_dict[u],
                    edgecolor="red",
                    width=width,
                    zorder=3,
                )
                leg_units.append(u)
                bottom += abs(neg_o)

            lcow_check += neg_o

    for flow, v in flows.items():

        f = v / d
        if v > plot_tol:
            ax.bar(
                [x],
                [f],
                bottom=bottom,
                hatch=flows_hatch,
                color=color_dict[flow],
                width=width,
                edgecolor="black",
                zorder=3,
            )
            leg_units.append(flow)
            bottom += f
        lcow_check += f

    for u in all_units:

        o = pos_opex.get(u, 0) / d
        if pos_opex.get(u, 0) > plot_tol:
            ax.bar(
                [x],
                [o],
                bottom=bottom,
                hatch=hatch_dict["OPEX"],
                color=color_dict[u],
                edgecolor="black",
                width=width,
                zorder=3,
            )
            leg_units.append(u)
            bottom += o
        lcow_check += o

        c = capex.get(u, 0) / d
        if capex.get(u, 0) > plot_tol:
            ax.bar(
                [x],
                [c],
                bottom=bottom,
                hatch=hatch_dict["CAPEX"],
                color=color_dict[u],
                edgecolor="black",
                width=width,
                zorder=3,
            )
            leg_units.append(u)
            bottom += c
        lcow_check += c

    leg_units = list(dict.fromkeys(leg_units))
    for u in leg_units:
        hatch = flows_hatch if u in flows else None
        handles.append(Patch(facecolor=color_dict[u], edgecolor="k", hatch=hatch))
        labels.append(_get_label(u))

    # Reverse legend order to align with plot order (bottom to top)
    handles = handles[::-1]
    labels = labels[::-1]

    # Add hatch legends to top
    handles[0:0] = [
        Patch(facecolor="white", edgecolor="k", label=k, hatch=v)
        for k, v in hatch_dict.items()
    ]
    labels[0:0] = list(hatch_dict.keys())

    ax.set_axisbelow(True)
    ax.grid(visible=True, zorder=0)
    ax.set_xticks([])
    ax.set_xlim(xmin, xmax)
    ax.set_ylabel(f"LCOW (\\$/m$^3$)" if not relative else f"Relative LCOW (%)")

    if lcow_neg < 0:
        if relative:
            # This should not happen
            msg = "Cannot plot relative values with negative LCOW contributions. "
            msg += "Contact a WaterTAP developer."
            raise FrozenPipes(msg)
        # If there are negative LCOW contributions, the net LCOW
        # is plotted with a yellow star
        ax.scatter(
            x,
            lcow_actual,
            marker="*",
            color="yellow",
            edgecolor="black",
            zorder=100,
            s=150,
        )
        handles.insert(
            0,
            Line2D(
                [0],
                [0],
                marker="*",
                color="white",
                markeredgecolor="black",
                markerfacecolor="yellow",
                zorder=100,
                markersize=150**0.5,
            ),
        )
        labels.insert(0, f"Net LCOW: {lcow_actual:.2f}")
        ax.use_sticky_edges = False
        ax.hlines(0, xmin, xmax, color="k", lw=2, zorder=1)

    max_per_col = 16
    ncols = ceil(len(labels) / max_per_col)

    ax.legend(
        handles=handles,
        labels=labels,
        ncols=ncols,
        bbox_to_anchor=(1.05, 1.0),
        loc="upper left",
        borderaxespad=0.0,
    )
    ax.autoscale(axis="y")
    fig.tight_layout()
    fig.subplots_adjust(right=0.8)

    if save_as is not None:
        fig.savefig(f"{save_as}.png", dpi=300, bbox_inches="tight")

    log.debug(f"Calculated LCOW: {lcow_check * d:.4f}, Actual LCOW: {lcow_actual}")

    if not isclose(lcow_check, lcow_actual / d, rel_tol=tol):
        print(f"Calculated LCOW: {lcow_check * d:.6f}")
        print(f"Actual LCOW: {lcow_actual:.6f}")
        if relative:
            print(f"Calculated LCOW (% of actual): {lcow_check:.6f}")
            print(f"Actual LCOW (% of actual): {lcow_actual / d:.6f}")
        raise ValueError(
            f"Calculated LCOW and actual LCOW differ by more than tolerance ({tol})"
        )

    if close_fig:
        plt.close(fig)

    return fig, ax


def plot_SEC_breakdown(
    blk,
    log,
    sec_name="specific_energy_consumption",
    fs_name="fs",
    relative=False,
    colormap=None,
    save_as=None,
    close_fig=False,
    tol=1e-5,
    plot_tol=1e-4,
    component_labels={},
    **kwargs,
):
    """
    This method creates a bar chart showing the breakdown of the specific energy consumption calculation
    into its contributing components. The breakdown can be done by individual unit
    models or by unit model type.

    Args:
        blk: a WaterTAP costing block
        log: logger object for logging messages
        sec_name: name of the specific energy consumption Expression to plot (default: "specific_energy_consumption")
        fs_name: name of the flowsheet (default: "fs")
        relative: whether to plot as fraction of total SEC (default: False)
        colormap: colormap to use for the plot (default: None; a custom colormap of tab20b + tab20c is created)
        save_as: filename to save the plot (default: None)
        close_fig: whether to close the figure after creation (default: False)
        tol: tolerance for SEC calculation check (default: 1e-5)
        plot_tol: minimum absolute value of SEC contribution to include in the plot (default: 1e-4)
        component_labels: dictionary mapping flowsheet unit names to custom unit names (default: {})
    """

    def _get_label(u):
        if u not in component_labels:
            return u.split(".")[-1]
        return component_labels[u]

    if blk.find_component(sec_name) is None:
        msg = f"{blk.name} does not have a component named {sec_name}."
        msg += " Use the add_specific_energy_consumption method to add a specific energy consumption Expression to the costing"
        msg += " block before calling plot_SEC_breakdown."
        raise RuntimeError(msg)

    sec_actual = pyo.value(blk.find_component(sec_name))
    sec_expr = blk.find_component(f"{sec_name}_component")
    sec_comp = {
        k.replace(f"{fs_name}.", ""): pyo.value(v)
        for k, v in sec_expr.items()
        if pyo.value(v) > 0
    }
    neg_sec_comp = {
        k.replace(f"{fs_name}.", ""): pyo.value(v)
        for k, v in sec_expr.items()
        if pyo.value(v) < 0
    }

    all_units = set(
        list(k for k in sec_comp.keys() if sec_comp[k] > plot_tol)
        + list(k for k in neg_sec_comp.keys() if abs(neg_sec_comp[k]) > plot_tol)
    )

    if colormap is None:
        # Create a custom colormap with 40 unique colors
        cmap = plt.get_cmap("tab20b").colors + plt.get_cmap("tab20c").colors
        cmap = mcolors.ListedColormap(cmap)
    else:
        cmap = plt.get_cmap(colormap)

    colors = [cmap(i / len(all_units)) for i in range(len(all_units))]
    color_dict = dict(zip(all_units, colors))

    # Start plotting
    fig, ax = plt.subplots()
    x = 0  # location of bar
    width = 0.5  # bar width
    xmin = x - width
    xmax = x + width

    # Create legend
    handles = []
    labels = []

    # For checking SEC calculation
    sec_check = 0
    # Sum of negative SEC contributions
    sec_neg = sum(neg_sec_comp.values())

    # Starting bottom needs to be at zero
    bottom = 0
    d = sec_actual * 1e-2 if relative else 1

    leg_units = list()

    if sec_neg < 0:
        log.info(
            "Cannot plot relative values with negative SEC contributions; setting 'relative' to False."
        )
        relative = False
        d = 1
        # Set bottom to the sum of negative SEC contributions
        bottom = sum(v / d for v in neg_sec_comp.values())
        for u, v in neg_sec_comp.items():
            neg_ec = v / d
            if abs(v) > plot_tol:
                # Only plot negative contributions above tolerance
                ax.bar(
                    [x],
                    [abs(neg_ec)],
                    bottom=bottom,
                    color=color_dict[u],
                    edgecolor="red",
                    width=width,
                    zorder=3,
                )
                leg_units.append(u)
                bottom += abs(neg_ec)
            sec_check += neg_ec

    for u in sec_comp:

        ec = sec_comp.get(u, 0) / d

        if sec_comp.get(u, 0) > plot_tol:
            # Only plot absolute contributions above tolerance
            ax.bar(
                [x],
                [ec],
                bottom=bottom,
                color=color_dict[u],
                edgecolor="black",
                width=width,
                zorder=3,
            )
            leg_units.append(u)
            bottom += ec
        sec_check += ec

    leg_units = list(dict.fromkeys(leg_units))
    for u in leg_units:
        handles.append(Patch(facecolor=color_dict[u], edgecolor="k"))
        labels.append(_get_label(u))

    # Reverse order to align with plot order (bottom to top)
    handles = handles[::-1]
    labels = labels[::-1]

    ax.set_axisbelow(True)
    ax.grid(visible=True, zorder=0)
    ax.set_xticks([])
    ax.set_xlim(xmin, xmax)
    ax.set_ylabel(f"SEC (kWh/m$^3$)" if not relative else f"Relative SEC (%)")

    if sec_neg < 0:
        if relative:
            # This should not happen
            msg = "Cannot plot relative values with negative SEC contributions. "
            msg += "Contact a WaterTAP developer."
            raise FrozenPipes(msg)
        # If there are negative SEC contributions, the net SEC
        # is plotted with a yellow star
        ax.scatter(
            x,
            sec_actual,
            marker="*",
            color="yellow",
            edgecolor="black",
            zorder=100,
            s=150,
        )
        handles.insert(
            0,
            Line2D(
                [0],
                [0],
                marker="*",
                color="white",
                markeredgecolor="black",
                markerfacecolor="yellow",
                zorder=100,
                markersize=150**0.5,
            ),
        )
        labels.insert(0, f"Net SEC: {sec_actual:.2f}")
        ax.use_sticky_edges = False
        ax.hlines(0, xmin, xmax, color="k", lw=2, zorder=1)

    max_per_col = 16
    ncols = ceil(len(labels) / max_per_col)

    ax.legend(
        handles=handles,
        labels=labels,
        ncols=ncols,
        bbox_to_anchor=(1.05, 1.0),
        loc="upper left",
        borderaxespad=0.0,
    )
    ax.autoscale(axis="y")
    fig.tight_layout()
    fig.subplots_adjust(right=0.8)

    if save_as is not None:
        fig.savefig(f"{save_as}.png", dpi=300, bbox_inches="tight")

    log.debug(f"Calculated SEC: {sec_check * d:.4f}, Actual SEC: {sec_actual}")

    if not isclose(sec_check, sec_actual / d, rel_tol=tol):
        print(f"Calculated SEC: {sec_check}")
        print(f"Actual SEC: {sec_actual / d}")
        raise ValueError(
            f"Calculated SEC and actual SEC differ by more than tolerance ({tol})"
        )

    if close_fig:
        plt.close(fig)

    return fig, ax
