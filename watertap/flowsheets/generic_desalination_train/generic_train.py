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
    value,
    Var,
    Constraint,
    Objective,
    TransformationFactory,
    NonNegativeReals,
    units as pyunits,
)
from idaes.core import (
    FlowsheetBlock,
)

import idaes.core.util.scaling as iscale
from idaes.core.util.model_statistics import degrees_of_freedom


from watertap.core.solvers import get_solver

from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
    ActivityCoefficientModel,
    DensityCalculation,
    MaterialFlowBasis,
)

from idaes.core.util.initialization import propagate_state

from idaes.models.unit_models import (
    Product,
    Feed,
)

import watertap.flowsheets.generic_desalination_train.utils.scale_utils as scaleTools
from pyomo.util.check_units import assert_units_consistent
from watertap.flowsheets.generic_desalination_train.costing import (
    generic_costing,
    stream_costing,
)
from pyomo.network import Arc
from watertap.flowsheets.generic_desalination_train.unit_operations import (
    desalter,
    mixer,
    separator,
)

import logging

_logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
formatter = logging.Formatter(
    "gen_train %(asctime)s %(levelname)s: %(message)s", "%H:%M:%S"
)
handler.setFormatter(formatter)
_logger.addHandler(handler)
_logger.setLevel(logging.DEBUG)

__author__ = "Alexander V. Dudchenko"


def main():
    m = build()
    initialize(m)
    m.fs.Pretreatment.separator.component_removal_percent["X"].fix(50)
    m.fs.Pretreatment.separator.separation_cost["X"].fix(0.5)
    m.fs.Valorizer.separator.product_value["X"].fix(1)
    m.fs.Valorizer.separator.component_removal_percent["X"].fix(50)

    m.fs.Desal_1.desalter.water_recovery.fix(80)
    m.fs.Desal_2.desalter.water_recovery.fix(50)
    m.fs.Desal_2.desalter.recovery_cost.fix(0.01)
    m.fs.Desal_2.desalter.recovery_cost_offset.fix(35)
    m.fs.Desal_3.desalter.water_recovery.unfix()
    m.fs.Desal_3.desalter.brine_water_mass_percent.fix(80)
    result = solve(m)
    return result


def build(
    train_type="Pretreatment>Desal1>Desal2>Crystalizer>Valorizer",
    water_source="generic",
):
    train_orders = {
        "Pretreatment>Desal1>Desal2>Crystalizer>Valorizer": {
            0: {
                "process_type": "separator",
                "process_name": "Pretreatment",
                "default_kwargs": {
                    "base_cost": 0.1,
                    "additive_dose": 0,
                    "additive_cost": 0,
                },
            },
            1: {
                "process_type": "desalter",
                "process_name": "Desal_1",
                "default_kwargs": {
                    "base_cost": 0.3,
                    "recovery_cost": 0,
                    "min_recovery": 0,
                    "max_recovery": 85,
                    "default_recovery": 85,
                },
            },
            2: {
                "process_type": "desalter",
                "process_name": "Desal_2",
                "default_kwargs": {
                    "base_cost": 0.5,
                    "recovery_cost": 0,
                    "min_recovery": 0,
                    "max_recovery": 50,
                    "default_recovery": 50,
                },
            },
            3: {
                "process_type": "desalter",
                "process_name": "Desal_3",
                "default_kwargs": {
                    "base_cost": 10,
                    "recovery_cost": 0.0,
                    "min_recovery": 0,
                    "max_recovery": 99.9,
                    "default_recovery": 90,
                },
            },
            4: {
                "process_type": "valorizer",
                "process_name": "Valorizer",
                "default_kwargs": {
                    "base_cost": 0,
                    "additive_dose": 0.0,
                    "additive_cost": 0.0,
                },
            },
        },
    }

    train_order = train_orders[train_type]
    """ set up mcas as our default prop package for tracking ions etc."""

    mcas_props = {
        "activity_coefficient_model": ActivityCoefficientModel.ideal,
        "density_calculation": DensityCalculation.constant,
        "solute_list": ["TDS", "X"],
        "mw_data": {"H2O": 0.01801528, "TDS": 1, "X": 1},
        "charge": {"TDS": 0.0, "X": 0.0},
        "material_flow_basis": MaterialFlowBasis.mass,
    }

    # defines solids used for precipitator tracking

    tracked_solids = ["TDS"]
    """ build all the processes and connection"""
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(**mcas_props)
    m.fs.costing = generic_costing.add_generic_costing_block()

    m.working_name = water_source + "_" + train_type.replace(">", "_")

    m.fs.feed = Feed(property_package=m.fs.properties)
    """  make sure these are constructed."""
    m.fs.feed.properties[0].conc_mass_phase_comp[...]
    m.fs.feed.properties[0].flow_mass_phase_comp[...]
    m.fs.feed.properties[0].flow_vol_phase[...]

    m.fs.product = Product(property_package=m.fs.properties)
    m.fs.disposal = Product(property_package=m.fs.properties)

    stream_costing.cost_stream(
        m.fs.costing, m.fs.product, 0, opt_name="Product distribution"
    )
    stream_costing.cost_stream(
        m.fs.costing, m.fs.disposal, 1, opt_name="Waste disposal"
    )
    stream_costing.cost_stream(m.fs.costing, m.fs.feed, 0, opt_name="Feed sourcing")

    add_flowsheet_level_constraints(m, m.fs)

    """ build main processes and arcs"""

    m.fs.process_order = []
    m.arc_order = []
    m.end_point_order = []
    for item_number in range(len(train_order.keys())):
        process_name = train_order[item_number]["process_type"]
        if train_order[item_number].get("default_kwargs") is None:
            train_order[item_number]["default_kwargs"] = {}
        else:
            train_order[item_number]["default_kwargs"].update(
                {"tracked_solids": tracked_solids}
            )
        build_selected_processes(m, **train_order[item_number])
        if item_number == 0:
            arc_name = "Feed_to_{}".format(process_name)
            m.fs.add_component(
                arc_name,
                Arc(
                    source=m.fs.feed.outlet,
                    destination=m.fs.process_order[-1]["feed"],
                ),
            )
            m.arc_order.append(m.fs.find_component(arc_name))
        else:
            prior_process_name = m.fs.process_order[-2]["process_name"]
            arc_name = "{}_to_{}".format(prior_process_name, process_name)
            m.fs.add_component(
                arc_name,
                Arc(
                    source=m.fs.process_order[-2]["effluent"],
                    destination=m.fs.process_order[-1]["feed"],
                ),
            )
            m.arc_order.append(m.fs.find_component(arc_name))

    product_streams = []
    disposal_streams = []
    for i, process_dict in enumerate(m.fs.process_order):
        if process_dict.get("product") is not None:
            product_streams.append(process_dict)
        if process_dict.get("disposal") is not None or i == len(m.fs.process_order) - 1:
            if (
                process_dict.get("disposal") is None
                and process_dict.get("effluent") is not None
            ):
                process_dict["disposal"] = process_dict["effluent"]
            if process_dict.get("disposal") is not None:
                disposal_streams.append(process_dict)
    streams = {"product": product_streams, "disposal": disposal_streams}
    for end_point, streams in streams.items():
        if len(streams) > 1:
            mixer.build_mixer(
                m,
                "{}_mixer".format(end_point),
                [process_dict["process_name"] for process_dict in streams],
            )
            mixer_unit = m.fs.find_component("{}_mixer".format(end_point))
            for stream in streams:
                arc_name = "{}_to_{}_mixer".format(stream["process_name"], end_point)
                m.fs.add_component(
                    arc_name,
                    Arc(
                        source=stream[end_point],
                        destination=mixer_unit.find_component(stream["process_name"]),
                    ),
                )
                m.end_point_order.append(m.fs.find_component(arc_name))
            arc_name = "{}_mixer_to_{}".format(end_point, end_point)
            m.fs.add_component(
                arc_name,
                Arc(
                    source=mixer_unit.outlet,
                    destination=m.fs.find_component(end_point).inlet,
                ),
            )
            m.end_point_order.append(m.fs.find_component(arc_name))
        else:
            arc_name = "{}_to_{}".format(streams[0]["process_name"], end_point)
            m.fs.add_component(
                arc_name,
                Arc(
                    source=streams[0][end_point],
                    destination=m.fs.find_component(end_point).inlet,
                ),
            )
            m.end_point_order.append(m.fs.find_component(arc_name))

    generic_costing.cost_process(
        m.fs.costing, m.fs.product.properties[0].flow_vol_phase["Liq"]
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    """set in initial operating conditions"""
    set_feed(m)

    """ set scaling for feed_composition"""
    scaleTools.set_default_scaling(
        m.fs.feed.properties[0],
        m.fs.properties,
        "flow_mol_phase_comp",
        factor=1.0,
    )
    scaleTools.set_default_scaling(
        m.fs.feed.properties[0],
        m.fs.properties,
        "flow_mass_phase_comp",
        factor=1.0,
    )
    m.fs.cost_objective = Objective(expr=m.fs.costing.LCOW)
    assert_units_consistent(m)
    _logger.info(f"DOFS:{degrees_of_freedom(m)}")
    assert degrees_of_freedom(m) == 0
    return m


def set_feed(m):
    m.fs.feed.properties[0].temperature.fix(298.15)
    m.fs.feed.properties[0].pressure.fix(101325)
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(1)
    m.fs.feed.properties[0].conc_mass_phase_comp["Liq", "TDS"].fix(3.5)  # kg/m3 or g/L
    m.fs.feed.properties[0].conc_mass_phase_comp["Liq", "X"].fix(0.3)  # kg/m3 or g/L
    update_feed(m.fs)


def add_flowsheet_level_constraints(m, blk):
    blk.water_recovery = Var(
        initialize=50,
        bounds=(0, 99.5),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="System Water Recovery",
    )
    blk.water_recovery.setub(99.5)
    blk.water_recovery.setlb(0)
    blk.eq_water_recovery = Constraint(
        expr=sum(blk.feed.properties[0].flow_mass_phase_comp["Liq", :])
        * blk.water_recovery
        / 100
        == blk.product.properties[0].flow_mass_phase_comp["Liq", "H2O"]
    )
    iscale.set_scaling_factor(blk.water_recovery, 1 / 100)
    iscale.constraint_scaling_transform(blk.eq_water_recovery, 1 / 100)

    blk.disposal_mass_flow_feed_basis = Var(
        list(blk.feed.properties[0].flow_mass_phase_comp.keys()),
        units=pyunits.kg / pyunits.m**3,
        bounds=(None, None),
    )

    @blk.Constraint(list(blk.feed.properties[0].flow_mass_phase_comp.keys()))
    def eq_disposal_mass_flow_feed_basis(blk, liq, ion):
        return (
            blk.disposal_mass_flow_feed_basis[liq, ion]
            == blk.disposal.properties[0].flow_mass_phase_comp[liq, ion]
            / blk.feed.properties[0].flow_vol_phase["Liq"]
        )


def initialize(m, solver=None, **kwargs):
    if solver is None:
        solver = get_solver()
    m.fs.feed.initialize(optarg=solver.options)
    # loop over all the process and their arcs and initialize
    # right now only supporting linear trains with no recycle
    for i, proc_dict in enumerate(m.fs.process_order):
        propagate_state(m.arc_order[i])
        if proc_dict.get("init_args") is not None:
            proc_dict["init"](
                m, proc_dict["process_block"], solver, **proc_dict.get("init_args")
            )
        else:
            proc_dict["init"](m, proc_dict["process_block"], solver)
    for arg in m.end_point_order:
        propagate_state(arg)

    # Only init if it was constructed
    if m.fs.find_component("product_mixer") is not None:
        m.fs.product_mixer.initialize(optarg=solver.options)
    m.fs.product.initialize(optarg=solver.options)
    update_feed(m.fs, solver)
    fix_conc_feed(
        m.fs,
    )
    setup_optimization(m)
    solve(m, solver)


def build_selected_processes(m, process_name, process_type, default_kwargs=None):
    m.fs.add_component(process_name, FlowsheetBlock(dynamic=False))
    block = m.fs.find_component(process_name)
    block.process_name = process_name
    if default_kwargs == None:
        default_kwargs = {}
    if "desalter" in process_type:
        desalter.build(m, block, **default_kwargs)
        m.fs.process_order.append(
            {
                "process_type": process_type,
                "process_name": process_name,
                "process_block": block,
                "init": desalter.initialize,
                "feed": block.desalter.inlet,
                "product": block.desalter.product,
                "effluent": block.desalter.brine,
                "disposal": None,
                "post_scale_func": None,
                "fix_vars_func": None,
                "unfix_opt_vars": desalter.unfix_opt_vars,
                "opt_guess": None,
                "init_args": None,
                "display_func": desalter.display,
                "default_scale_func": None,
            }
        )
    elif "separator" in process_type:
        separator.build(m, block, **default_kwargs)
        m.fs.process_order.append(
            {
                "process_type": process_type,
                "process_name": process_name,
                "process_block": block,
                "init": separator.initialize,
                "feed": block.separator.inlet,
                "product": None,
                "effluent": block.separator.treated,
                "disposal": block.separator.product,
                "post_scale_func": None,
                "fix_vars_func": None,
                "unfix_opt_vars": None,
                "opt_guess": None,
                "init_args": None,
                "display_func": separator.display,
                "default_scale_func": None,
            }
        )
    elif "valorizer" in process_type:
        default_kwargs.update({"valorizer_costing": True})
        separator.build(m, block, **default_kwargs)
        m.fs.process_order.append(
            {
                "process_type": process_type,
                "process_name": process_name,
                "process_block": block,
                "init": separator.initialize,
                "feed": block.separator.inlet,
                "product": None,
                "effluent": block.separator.treated,
                "disposal": None,
                "post_scale_func": None,
                "fix_vars_func": None,
                "unfix_opt_vars": None,
                "opt_guess": None,
                "init_args": None,
                "display_func": separator.display,
                "default_scale_func": None,
            }
        )
    else:
        raise TypeError("Selected process type {} unavailable".format(process_type))


def setup_optimization(m):
    """unfix all vars in processes for optimization"""
    for i, process_dict in enumerate(m.fs.process_order):
        if process_dict.get("unfix_opt_vars") != None:
            process_dict["unfix_opt_vars"](m, process_dict["process_block"])
    """ enable objective function """
    # m.fs.water_recovery.fix(95)
    m.fs.cost_objective.activate()


def solve(m, solver=None, tee=None):
    if m.find_component("fs") is None:
        update_feed(m, solver)
    else:
        update_feed(m.fs, solver)
    _logger.info(f"DOFS:{degrees_of_freedom(m)}")
    if solver is None:
        solver = get_solver()

    result = solver.solve(m)
    if m.find_component("fs") is None:
        fix_conc_feed(m)
        display_processes(m)
    else:
        fix_conc_feed(m.fs)
        display_processes(m.fs)
    return result


def fix_conc_feed(blk):
    blk.feed.properties[0].flow_mass_phase_comp.unfix()
    blk.feed.properties[0].conc_mass_phase_comp.fix()
    blk.feed.properties[0].conc_mass_phase_comp["Liq", "H2O"].unfix()
    blk.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix()


def update_feed(blk, solver=None):
    """used to update feed composition based on what user enters
    in ui, which uses concentration, rather then mass flow basis"""
    if solver == None:
        solver = get_solver()
    _logger.info("solved feed")
    fix_conc_feed(blk)
    solver.solve(blk.feed.properties[0])
    blk.feed.properties[0].conc_mass_phase_comp.unfix()
    blk.feed.properties[0].flow_mass_phase_comp.fix()


def display_processes(blk):
    _logger.info("--------Display results for {}--------".format("Overall process"))
    _logger.info("LCOW {}  $/m3".format(value(blk.costing.LCOW)))
    for name, var in blk.costing.LCOW_unit.items():
        _logger.info(f"{name} LCOW {value(var)}")
    _logger.info("Water recovery {}%".format(blk.water_recovery.value))
    _logger.info(
        f"Feed flow (kg/s) {value(blk.feed.properties[0].flow_mass_phase_comp['Liq','H2O'])}"
    )
    _logger.info(
        f"Product flow (kg/s) {value(blk.product.properties[0].flow_mass_phase_comp['Liq','H2O'])}"
    )
    _logger.info(
        f"Disposal flow (kg/s) {value(blk.disposal.properties[0].flow_mass_phase_comp['Liq','H2O'])}"
    )

    _logger.info(f"Annual feed cost ($) {value(blk.feed.annual_cost)}")
    _logger.info(f"Annual product cost ($) {value(blk.product.annual_cost)}")
    _logger.info(f"Annual disposal cost ($) {value(blk.disposal.annual_cost)}")

    _logger.info("--------------------------------")
    for proc_dict in blk.process_order:
        if proc_dict.get("display_func") is not None:
            _logger.info(
                "--------Display results for {}--------".format(
                    proc_dict["process_name"]
                )
            )
            proc_dict.get("display_func")(blk, proc_dict["process_block"])
            _logger.info("--------------------------------")


if __name__ == "__main__":
    main()
