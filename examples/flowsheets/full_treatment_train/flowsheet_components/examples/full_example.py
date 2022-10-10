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

from pyomo.environ import (
    ConcreteModel,
    Objective,
    Expression,
    Constraint,
    TransformationFactory,
    value,
    Param,
    Block,
)
from pyomo.environ import units as pyunits
from pyomo.network import Arc
from pyomo.util import infeasible
from idaes.core import FlowsheetBlock
from idaes.core.util.scaling import calculate_scaling_factors
from idaes.core.util.initialization import propagate_state
from watertap.examples.flowsheets.full_treatment_train.flowsheet_components import (
    pretreatment_NF,
    desalination,
    translator_block,
    feed_block,
    gypsum_saturation_index,
    costing,
)
from watertap.examples.flowsheets.full_treatment_train.model_components import (
    property_models,
)
from watertap.examples.flowsheets.full_treatment_train.util import (
    solve_block,
    check_dof,
    check_build,
    check_scaling,
)

"""Flowsheet example that satisfy minimum viable product requirements"""


def build_flowsheet_mvp_NF(m, **kwargs):
    """
    Build a flowsheet with NF pretreatment and RO.
    """
    # set up keyword arguments for the sections of treatment train
    kwargs_pretreatment = {
        k: kwargs[k]
        for k in (
            "has_bypass",
            "NF_type",
            "NF_base",
        )
    }
    kwargs_desalination = {
        k: kwargs[k]
        for k in (
            "has_desal_feed",
            "is_twostage",
            "has_ERD",
            "RO_type",
            "RO_base",
            "RO_level",
        )
    }
    # build flowsheet
    property_models.build_prop(m, base="ion")
    pretrt_port = pretreatment_NF.build_pretreatment_NF(m, **kwargs_pretreatment)

    property_models.build_prop(m, base=kwargs["RO_base"])
    desal_port = desalination.build_desalination(m, **kwargs_desalination)

    property_models.build_prop(m, base="eNRTL")

    translator_block.build_tb(
        m,
        base_inlet=kwargs["NF_base"],
        base_outlet=kwargs["RO_base"],
        name_str="tb_pretrt_to_desal",
    )

    # set up Arcs between pretreatment and desalination
    m.fs.s_pretrt_tb = Arc(
        source=pretrt_port["out"], destination=m.fs.tb_pretrt_to_desal.inlet
    )
    m.fs.s_tb_desal = Arc(
        source=m.fs.tb_pretrt_to_desal.outlet, destination=desal_port["in"]
    )

    # add gypsum saturation index calculations
    gypsum_saturation_index.build(m, section="desalination", **kwargs_desalination)
    gypsum_saturation_index.build(m, section="pretreatment", **kwargs_desalination)

    # new initialization
    if kwargs["NF_type"] == "ZO":
        m.fs.NF.area.fix(175)
    if kwargs["has_bypass"]:
        m.fs.splitter.split_fraction[0, "bypass"].fix(0.50)
    m.fs.RO.area.fix(80)
    m.fs.pump_RO.control_volume.properties_out[0].pressure.fix(60e5)
    if kwargs["is_twostage"]:
        m.fs.RO2.area.fix(20)
        m.fs.pump_RO2.control_volume.properties_out[0].pressure.fix(90e5)

    # touch some properties used in optimization
    if kwargs["is_twostage"]:
        product_water_sb = m.fs.mixer_permeate.mixed_state[0]
        RO_waste_sb = m.fs.RO2.feed_side.properties[0, 1]
    else:
        product_water_sb = m.fs.RO.mixed_permeate[0]
        RO_waste_sb = m.fs.RO.feed_side.properties[0, 1]

    # NOTE: Building the costing here means it gets
    #       initialized during the simulation phase.
    #       This helps model stability.
    m.fs.feed.properties[0].flow_vol
    m.fs.feed.properties[0].conc_mol_phase_comp["Liq", "Ca"]

    m.fs.tb_pretrt_to_desal.properties_in[0].flow_vol
    m.fs.tb_pretrt_to_desal.properties_in[0].conc_mol_phase_comp["Liq", "Ca"]

    product_water_sb.flow_vol
    RO_waste_sb.flow_vol

    m.fs.system_recovery = Expression(
        expr=product_water_sb.flow_vol / m.fs.feed.properties[0].flow_vol
    )
    m.fs.total_work = Expression(
        expr=m.fs.pump_RO.work_mechanical[0]
        + (m.fs.pump_RO2.work_mechanical[0] if kwargs["is_twostage"] else 0.0)
    )

    # annual water production
    m.fs.treated_flow_vol = Expression(expr=product_water_sb.flow_vol)
    costing.build_costing(m, **kwargs)

    return m


def set_up_optimization(m, system_recovery=0.7, **kwargs_flowsheet):
    is_twostage = kwargs_flowsheet["is_twostage"]

    # scale
    calculate_scaling_factors(m)

    # unfix variables
    m.fs.splitter.split_fraction[0, "bypass"].unfix()
    m.fs.splitter.split_fraction[0, "bypass"].setlb(0.001)
    m.fs.splitter.split_fraction[0, "bypass"].setub(0.99)

    m.fs.NF.area.unfix()
    m.fs.NF.area.setlb(10)
    m.fs.NF.area.setub(1000)

    m.fs.pump_RO.control_volume.properties_out[0].pressure.unfix()
    m.fs.pump_RO.control_volume.properties_out[0].pressure.setlb(20e5)
    m.fs.pump_RO.control_volume.properties_out[0].pressure.setub(75e5)

    m.fs.RO.area.unfix()
    m.fs.RO.area.setlb(10)
    m.fs.RO.area.setub(300)

    # Set lower bound for water flux at the RO outlet, based on a minimum net driving pressure, NDPmin
    m.fs.RO.NDPmin = Param(initialize=1e5, mutable=True, units=pyunits.Pa)

    m.fs.RO.flux_mass_phase_comp[0, 1, "Liq", "H2O"].setlb(
        value(m.fs.RO.A_comp[0, "H2O"] * m.fs.RO.dens_solvent * m.fs.RO.NDPmin)
    )

    if is_twostage:
        m.fs.max_allowable_pressure = Param(
            initialize=120e5, mutable=True, units=pyunits.pascal
        )
        m.fs.pump_RO2.control_volume.properties_out[0].pressure.unfix()
        m.fs.pump_RO2.control_volume.properties_out[0].pressure.setlb(20e5)
        m.fs.pump_RO2.control_volume.properties_out[0].pressure.setub(
            m.fs.max_allowable_pressure
        )

        m.fs.RO2.area.unfix()
        m.fs.RO2.area.setlb(10)
        m.fs.RO2.area.setub(300)

        # Set lower bound for water flux at the RO outlet, based on a minimum net driving pressure, NDPmin
        m.fs.RO2.NDPmin = Param(initialize=1e5, mutable=True, units=pyunits.Pa)
        m.fs.RO2.flux_mass_phase_comp[0, 1, "Liq", "H2O"].setlb(
            value(m.fs.RO2.A_comp[0, "H2O"] * m.fs.RO2.dens_solvent * m.fs.RO2.NDPmin)
        )

    # add additional constraints
    # fixed system recovery
    m.fs.system_recovery_target = Param(initialize=system_recovery, mutable=True)
    m.fs.system_recovery_tol = Param(initialize=5e-3, mutable=True)
    m.fs.eq_system_recovery = Constraint(
        expr=(
            m.fs.system_recovery_target,
            m.fs.system_recovery,
            m.fs.system_recovery_target + m.fs.system_recovery_tol,
        )
    )

    # saturation index
    m.fs.max_saturation_index = Param(initialize=1.0, mutable=True)
    m.fs.eq_max_saturation_index_desal = Constraint(
        expr=m.fs.desal_saturation.saturation_index <= m.fs.max_saturation_index
    )

    m.fs.max_conc_factor_target = Param(initialize=3.5, mutable=True)
    m.fs.eq_max_conc_NF = Constraint(
        expr=m.fs.NF.feed_side.properties_out[0].mass_frac_phase_comp["Liq", "Ca"]
        <= m.fs.max_conc_factor_target
        * m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "Ca"]
    )

    # set objective
    m.fs.objective = Objective(expr=m.fs.costing.LCOW)

    # set additional constraints to limit local minima
    # NOTE: doesn't seem necessary with new objective
    if False:
        m.fs.inequality_RO_area = Constraint(expr=m.fs.RO.area >= m.fs.RO2.area)
        min_pressure_increase = 1e5
        m.fs.inequality_RO_pressure = Constraint(
            expr=m.fs.pump_RO.control_volume.properties_out[0].pressure
            + min_pressure_increase
            <= m.fs.pump_RO2.control_volume.properties_out[0].pressure
        )
        m.fs.inequality_RO_permeate = Constraint(
            expr=m.fs.RO.permeate_side.properties_mixed[0].flow_vol_phase["Liq"]
            >= m.fs.RO2.permeate_side.properties_mixed[0].flow_vol_phase["Liq"]
        )

    check_dof(m, dof_expected=6 if is_twostage else 4)
    # solve_block(m, tee=False, fail_flag=True)


def optimize(m):
    solve_block(m, tee=False, fail_flag=True)


def solve_flowsheet_mvp_NF(**kwargs):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    build_flowsheet_mvp_NF(m, **kwargs)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scale
    pretreatment_NF.scale_pretreatment_NF(m, **kwargs)
    calculate_scaling_factors(m.fs.tb_pretrt_to_desal)
    desalination.scale_desalination(m, **kwargs)
    calculate_scaling_factors(m)

    # initialize
    optarg = {"nlp_scaling_method": "user-scaling"}
    pretreatment_NF.initialize_pretreatment_NF(m, **kwargs)
    m.fs.pretrt_saturation.properties.initialize(optarg=optarg)
    propagate_state(m.fs.s_pretrt_tb)
    m.fs.tb_pretrt_to_desal.initialize(optarg=optarg)
    propagate_state(m.fs.s_tb_desal)
    desalination.initialize_desalination(m, **kwargs)
    m.fs.desal_saturation.properties.initialize()

    m.fs.costing.initialize()

    # check_build(m)
    # check_scaling(m)

    check_dof(m)
    solve_block(m, tee=False, fail_flag=True)

    pretreatment_NF.display_pretreatment_NF(m, **kwargs)
    m.fs.tb_pretrt_to_desal.report()
    desalination.display_desalination(m, **kwargs)
    print(
        "desalination solubility index:", value(m.fs.desal_saturation.saturation_index)
    )
    print(
        "pretreatment solubility index:", value(m.fs.pretrt_saturation.saturation_index)
    )
    print("water recovery:", value(m.fs.system_recovery))
    print("LCOW:", value(m.fs.costing.LCOW))
    print("CP modulus:", value(m.fs.desal_saturation.cp_modulus))

    return m


def solve_optimization(system_recovery=0.75, **kwargs_flowsheet):

    m = solve_flowsheet_mvp_NF(**kwargs_flowsheet)

    print("\n****** Optimization *****\n")
    set_up_optimization(m, system_recovery=system_recovery, **kwargs_flowsheet)
    optimize(m)

    pretreatment_NF.display_pretreatment_NF(m, **kwargs_flowsheet)
    m.fs.tb_pretrt_to_desal.report()
    desalination.display_desalination(m, **kwargs_flowsheet)
    costing.display_costing(m)
    print(
        "desalination saturation index:", value(m.fs.desal_saturation.saturation_index)
    )
    print(
        "pretreatment saturation index:", value(m.fs.pretrt_saturation.saturation_index)
    )
    print(
        "pretreatment Ca concentration factor:",
        value(
            m.fs.NF.feed_side.properties_out[0].mass_frac_phase_comp["Liq", "Ca"]
            / m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "Ca"]
        ),
    )
    print("water recovery:", value(m.fs.system_recovery))
    print("CP modulus:", value(m.fs.desal_saturation.cp_modulus))
    return m


if __name__ == "__main__":
    import sys

    kwargs_flowsheet = {
        "has_bypass": True,
        "has_desal_feed": False,
        "is_twostage": True,
        "has_ERD": True,
        "NF_type": "ZO",
        "NF_base": "ion",
        "RO_type": "0D",
        "RO_base": "TDS",
        "RO_level": "detailed",
    }
    if len(sys.argv) == 1:
        m = solve_flowsheet_mvp_NF(**kwargs_flowsheet)
    else:
        m = solve_optimization(system_recovery=float(sys.argv[1]), **kwargs_flowsheet)
