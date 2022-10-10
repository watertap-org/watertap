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
"""
mutable parameters for optimization:
    m.fs.system_recovery_target
    m.fs.max_allowable_pressure
"""
from pyomo.environ import ConcreteModel, TransformationFactory, Param, value
from pyomo.environ import units as pyunits

from idaes.core import FlowsheetBlock
from idaes.core.util.initialization import propagate_state
from idaes.core.util.scaling import calculate_scaling_factors

from watertap.examples.flowsheets.full_treatment_train.util import (
    solve_block,
    check_dof,
)

import watertap.examples.flowsheets.full_treatment_train.analysis.flowsheet_single_stage as single_stage
from watertap.examples.flowsheets.full_treatment_train.analysis.flowsheet_single_stage import (
    build,
    build_components,
    scale,
    initialize,
    report,
    optimize,
)


desal_kwargs = {
    "has_desal_feed": False,
    "is_twostage": True,
    "has_ERD": True,
    "RO_type": "0D",
    "RO_base": "TDS",
    "RO_level": "detailed",
}


def set_optimization_components(m, system_recovery, **kwargs):
    single_stage.set_optimization_components(m, system_recovery, **kwargs)

    m.fs.max_allowable_pressure = Param(
        initialize=300e5, mutable=True, units=pyunits.pascal
    )
    m.fs.pump_RO2.control_volume.properties_out[0].pressure.unfix()
    m.fs.pump_RO2.control_volume.properties_out[0].pressure.setlb(20e5)
    m.fs.pump_RO2.control_volume.properties_out[0].pressure.setub(
        m.fs.max_allowable_pressure
    )

    m.fs.RO2.area.unfix()
    m.fs.RO2.area.setlb(1)
    m.fs.RO2.area.setub(300)

    m.fs.RO2.N_Re[0, 0].unfix()

    # Set lower bound for water flux at the RO outlet, based on a minimum net driving pressure, NDPmin
    m.fs.RO2.NDPmin = Param(initialize=1e5, mutable=True, units=pyunits.Pa)
    m.fs.RO2.flux_mass_phase_comp[0, 1, "Liq", "H2O"].setlb(
        value(m.fs.RO2.A_comp[0, "H2O"] * m.fs.RO2.dens_solvent * m.fs.RO2.NDPmin)
    )


def set_up_optimization(m, system_recovery=0.50, **kwargs):
    set_optimization_components(m, system_recovery, **kwargs)
    check_dof(m, 6)


def solve_flowsheet(**desal_kwargs):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    build(m, **desal_kwargs)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scale
    scale(m, **desal_kwargs)
    calculate_scaling_factors(m)

    # initialize
    m.fs.feed.initialize()
    propagate_state(m.fs.s_pretrt_tb)
    initialize(m, **desal_kwargs)

    check_dof(m)
    solve_block(m, tee=False, fail_flag=True)

    # report
    print("===================================" "\n          Simulation          ")
    report(m, **desal_kwargs)

    return m


def optimize_flowsheet(system_recovery=0.50, **kwargs):
    m = solve_flowsheet(**kwargs)
    set_up_optimization(m, system_recovery=system_recovery, **kwargs)
    optimize(m)
    print("===================================" "\n       Optimization            ")
    report(m, **kwargs)

    return m


if __name__ == "__main__":
    import sys

    if len(sys.argv) == 1:
        m = solve_flowsheet(**desal_kwargs)
    else:
        m = optimize_flowsheet(system_recovery=float(sys.argv[1]), **desal_kwargs)
