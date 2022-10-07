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

"""Desalination flowsheet components"""

from pyomo.environ import ConcreteModel, TransformationFactory, Constraint
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from idaes.models.unit_models import Mixer
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    set_scaling_factor,
    get_scaling_factor,
    constraint_scaling_transform,
)
from idaes.core.util.initialization import propagate_state
from watertap.unit_models.pressure_changer import Pump, EnergyRecoveryDevice
from watertap.examples.flowsheets.full_treatment_train.flowsheet_components import (
    feed_block,
)
from watertap.examples.flowsheets.full_treatment_train.model_components import (
    unit_separator,
    unit_0DRO,
    unit_1DRO,
    property_models,
)
from watertap.examples.flowsheets.full_treatment_train.util import (
    solve_block,
    check_dof,
)
from idaes.core.util.model_statistics import fixed_variables_generator


def build_desalination(
    m,
    has_desal_feed=False,
    is_twostage=False,
    has_ERD=False,
    RO_type="0D",
    RO_base="TDS",
    RO_level="simple",
):
    """
    Builds RO desalination including specified feed and auxiliary equipment.

    Args:
        has_desal_feed: True or False, default = False,
            if True a feed block is created and specified to the standard feed
        RO_type: 'Sep', '0D', or 1D, default = '0D'
        RO_level: 'simple' or 'detailed', default = 'simple'
        RO_base: 'TDS' only, default = 'ion'
    """
    desal_port = {}
    prop = property_models.get_prop(m, base=RO_base)

    if has_desal_feed:
        # build feed
        feed_block.build_feed(m, base=RO_base)

    # build RO
    if RO_type == "Sep":
        unit_separator.build_SepRO(m, base=RO_base)
    elif RO_type == "0D":
        unit_0DRO.build_RO(m, base="TDS", level=RO_level)
    elif RO_type == "1D":
        unit_1DRO.build_RO(m, base="TDS", level=RO_level)
    else:
        raise ValueError(
            "Unexpected model type {RO_type} provided to build_desalination"
            "".format(RO_type=RO_type)
        )
    if is_twostage:
        if RO_type == "0D":
            unit_0DRO.build_RO(m, base="TDS", level=RO_level, name_str="RO2")
        elif RO_type == "1D":
            unit_1DRO.build_RO(m, base="TDS", level=RO_level, name_str="RO2")
        else:
            raise ValueError(
                "Unexpected model type {RO_type} provided to build_desalination when is_twostage is True"
                "".format(RO_type=RO_type)
            )

    if has_ERD:
        if RO_type == "Sep":
            raise ValueError(
                "Unexpected model type {RO_type} provided to build_desalination when has_ERD is True"
                "".format(RO_type=RO_type)
            )
        m.fs.ERD = EnergyRecoveryDevice(property_package=prop)

    # auxiliary units
    if RO_type == "Sep":
        # build auxiliary units (none)

        # connect models
        if has_desal_feed:
            m.fs.s_desal_feed_RO = Arc(
                source=m.fs.feed.outlet, destination=m.fs.RO.inlet
            )

        # specify (RO already specified)

        # inlet/outlet ports for pretreatment
        if not has_desal_feed:
            desal_port["in"] = m.fs.RO.inlet

    elif RO_type == "0D" or RO_type == "1D":
        # build auxiliary units
        m.fs.pump_RO = Pump(property_package=prop)
        if is_twostage:
            m.fs.pump_RO2 = Pump(property_package=prop)
            m.fs.mixer_permeate = Mixer(property_package=prop, inlet_list=["RO", "RO2"])

        # connect models
        if has_desal_feed:
            m.fs.s_desal_feed_pumpRO = Arc(
                source=m.fs.feed.outlet, destination=m.fs.pump_RO.inlet
            )

        m.fs.s_desal_pumpRO_RO = Arc(
            source=m.fs.pump_RO.outlet, destination=m.fs.RO.inlet
        )

        if is_twostage:
            m.fs.s_desal_RO_pumpRO2 = Arc(
                source=m.fs.RO.retentate, destination=m.fs.pump_RO2.inlet
            )
            m.fs.s_desal_pumpRO2_RO2 = Arc(
                source=m.fs.pump_RO2.outlet, destination=m.fs.RO2.inlet
            )
            m.fs.s_desal_permeateRO_mixer = Arc(
                source=m.fs.RO.permeate, destination=m.fs.mixer_permeate.RO
            )
            m.fs.s_desal_permeateRO2_mixer = Arc(
                source=m.fs.RO2.permeate, destination=m.fs.mixer_permeate.RO2
            )

        if has_ERD:
            if is_twostage:
                m.fs.s_desal_RO2_ERD = Arc(
                    source=m.fs.RO2.retentate, destination=m.fs.ERD.inlet
                )
            else:
                m.fs.s_desal_RO_ERD = Arc(
                    source=m.fs.RO.retentate, destination=m.fs.ERD.inlet
                )

        # specify (RO already specified, Pump 2 DOF, ERD 2 DOF)
        m.fs.pump_RO.efficiency_pump.fix(0.80)
        m.fs.pump_RO.control_volume.properties_out[0].pressure.fix(50e5)
        if is_twostage:
            m.fs.pump_RO2.efficiency_pump.fix(0.80)
            m.fs.pump_RO2.control_volume.properties_out[0].pressure.fix(55e5)

        if has_ERD:
            m.fs.ERD.efficiency_pump.fix(0.95)
            m.fs.ERD.outlet.pressure[0].fix(101325)

        # inlet/outlet ports for pretreatment
        if not has_desal_feed:
            desal_port["in"] = m.fs.pump_RO.inlet

    if is_twostage:
        desal_port["out"] = m.fs.mixer_permeate.outlet
        if has_ERD:
            desal_port["waste"] = m.fs.ERD.outlet
        else:
            desal_port["waste"] = m.fs.RO2.retentate
    else:
        desal_port["out"] = m.fs.RO.permeate
        if has_ERD:
            desal_port["waste"] = m.fs.ERD.outlet
        else:
            desal_port["waste"] = m.fs.RO.retentate

    return desal_port


def scale_desalination(m, **kwargs):
    """
    Scales the model created by build_desalination.
    Arguments:
        m: pyomo concrete model with a built desalination train
        **kwargs: same keywords as provided to the build_desalination function
    """

    if kwargs["has_desal_feed"]:
        calculate_scaling_factors(m.fs.feed)

    if kwargs["RO_type"] == "1D":
        set_scaling_factor(m.fs.RO.feed_side.area, 1e3)
        set_scaling_factor(m.fs.RO.area, 1e-1)
        set_scaling_factor(m.fs.RO.width, 1)

    calculate_scaling_factors(m.fs.RO)

    if kwargs["is_twostage"]:
        if kwargs["RO_type"] == "1D":
            set_scaling_factor(m.fs.RO2.feed_side.area, 1e3)
            set_scaling_factor(m.fs.RO2.area, 1e-1)
            set_scaling_factor(m.fs.RO2.width, 1)
        calculate_scaling_factors(m.fs.RO2)

    if kwargs["RO_type"] == "0D" or kwargs["RO_type"] == "1D":
        set_scaling_factor(m.fs.pump_RO.control_volume.work, 1e-3)
        calculate_scaling_factors(m.fs.pump_RO)
        set_scaling_factor(
            m.fs.pump_RO.ratioP, 1
        )  # TODO: IDAES should have a default and link to the constraint

        if kwargs["is_twostage"]:
            set_scaling_factor(m.fs.pump_RO2.control_volume.work, 1e-3)
            calculate_scaling_factors(m.fs.pump_RO2)
            set_scaling_factor(
                m.fs.pump_RO2.ratioP, 1
            )  # TODO: IDAES should have a default and link to the constraint
            calculate_scaling_factors(m.fs.mixer_permeate)
            set_scaling_factor(
                m.fs.mixer_permeate.minimum_pressure,
                get_scaling_factor(m.fs.mixer_permeate.mixed_state[0].pressure),
            )  # TODO: IDAES should have a default and link to the constraint
            for c in [
                m.fs.mixer_permeate.minimum_pressure_constraint[0, 1],
                m.fs.mixer_permeate.minimum_pressure_constraint[0, 2],
                m.fs.mixer_permeate.mixture_pressure[0.0],
            ]:
                constraint_scaling_transform(
                    c, get_scaling_factor(m.fs.mixer_permeate.minimum_pressure)
                )

    if kwargs["has_ERD"]:
        set_scaling_factor(m.fs.ERD.control_volume.work, 1e-3)
        set_scaling_factor(
            m.fs.ERD.ratioP, 1
        )  # TODO: IDAES should have a default and link to the constraint
        calculate_scaling_factors(m.fs.ERD)


def initialize_desalination(m, **kwargs):
    """
    Initialized the model created by build_desalination.
    Arguments:
        m: pyomo concrete model with a built desalination train
        **kwargs: same keywords as provided to the build_desalination function
    """
    optarg = {"nlp_scaling_method": "user-scaling"}

    if kwargs["has_desal_feed"]:
        m.fs.feed.initialize(optarg=optarg)

    if kwargs["RO_type"] == "Sep":
        if kwargs["has_desal_feed"]:
            propagate_state(m.fs.s_desal_feed_RO)
        m.fs.RO.mixed_state[
            0
        ].mass_frac_phase_comp  # touch properties to have a constraint on stateblock
        m.fs.RO.permeate_state[0].mass_frac_phase_comp
        m.fs.RO.retentate_state[0].mass_frac_phase_comp
        # m.fs.RO.initialize(optarg=optarg)  # IDAES error on initializing separators, simple enough to not need it

    elif kwargs["RO_type"] == "0D" or kwargs["RO_type"] == "1D":
        if kwargs["has_desal_feed"]:
            propagate_state(m.fs.s_desal_feed_pumpRO)
        m.fs.pump_RO.initialize(optarg=optarg)
        propagate_state(m.fs.s_desal_pumpRO_RO)
        m.fs.RO.initialize(optarg=optarg)
        if kwargs["is_twostage"]:
            propagate_state(m.fs.s_desal_RO_pumpRO2)
            m.fs.pump_RO2.initialize(optarg=optarg)
            propagate_state(m.fs.s_desal_pumpRO2_RO2)
            m.fs.RO2.initialize(optarg=optarg)
            propagate_state(m.fs.s_desal_permeateRO_mixer)
            propagate_state(m.fs.s_desal_permeateRO2_mixer)
            m.fs.mixer_permeate.initialize(optarg=optarg)

    if kwargs["has_ERD"]:
        if kwargs["is_twostage"]:
            propagate_state(m.fs.s_desal_RO2_ERD)
        else:
            propagate_state(m.fs.s_desal_RO_ERD)
        m.fs.ERD.initialize(optarg=optarg)


def display_desalination(m, **kwargs):
    if kwargs["has_desal_feed"]:
        m.fs.feed.report()

    if kwargs["RO_type"] == "Sep":
        m.fs.RO.report()

    elif kwargs["RO_type"] == "0D" or kwargs["RO_type"] == "1D":
        m.fs.pump_RO.report()
        m.fs.RO.report()

        if kwargs["is_twostage"]:
            m.fs.pump_RO2.report()
            m.fs.RO2.report()
            m.fs.mixer_permeate.report()

    if kwargs["has_ERD"]:
        m.fs.ERD.report()


def solve_desalination(**kwargs):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    property_models.build_prop(m, base="TDS")
    build_desalination(m, **kwargs)
    TransformationFactory("network.expand_arcs").apply_to(m)

    scale_desalination(m, **kwargs)

    initialize_desalination(m, **kwargs)

    check_dof(m)

    solve_block(m)

    display_desalination(m, **kwargs)

    return m


if __name__ == "__main__":
    m = solve_desalination(
        has_desal_feed=True,
        is_twostage=True,
        has_ERD=True,
        RO_type="1D",
        RO_base="TDS",
        RO_level="detailed",
    )
