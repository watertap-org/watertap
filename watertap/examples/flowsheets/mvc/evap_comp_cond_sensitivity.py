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
    assert_optimal_termination,
    value,
    TerminationCondition,
)
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
import numpy as np
import matplotlib.pyplot as plt

from watertap.unit_models.mvc.components import Evaporator, Condenser
import watertap.property_models.seawater_prop_pack as props_sw
import watertap.property_models.water_prop_pack as props_w
import evaporator_compressor_condenser as mvc_unit


def main():
    area_min = 300
    area_max = 600
    n = 5
    areas = np.linspace(area_min, area_max, n)
    rr_with_fixed_area = sweep_area(areas)
    # using the recoveries found, solve for the area - should be the same
    areas_with_fixed_rr = sweep_rr(rr_with_fixed_area)

    print("Area sweep:", areas)
    print("Recovery solved for: ", rr_with_fixed_area)
    print("Areas solved for with fixed recovery: ", areas_with_fixed_rr)
    plt.figure()
    plt.plot(areas, rr_with_fixed_area, "-bo", label=r"Fixed area")
    plt.plot(areas_with_fixed_rr, rr_with_fixed_area, "-ro", label=r"Fixed recovery")
    plt.xlabel(r"area [m^2]")
    plt.ylabel(r"recovery [-]")
    plt.legend()
    plt.show()


def sweep_area(area_sweep):
    # DOES NOT REBUILD THE MODEL
    recoveries = np.zeros(np.shape(area_sweep))
    n = np.size(area_sweep)
    # fix work
    for i in range(n):  # for rebuilding with each run
        m = mvc_unit.build()
        mvc_unit.set_operating_conditions(m)
        # for i in range(n): #for not rebuilding the model
        m.fs.evaporator.area.fix(area_sweep[i])
        mvc_unit.initialize_system(m)
        solver = get_solver()
        results = solver.solve(m, tee=False)

        recovery = m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp[
            "Vap", "H2O"
        ].value / (
            m.fs.evaporator.properties_feed[0].flow_mass_phase_comp["Liq", "TDS"].value
            + m.fs.evaporator.properties_feed[0]
            .flow_mass_phase_comp["Liq", "H2O"]
            .value
        )

        if results.solver.termination_condition == TerminationCondition.infeasible:
            recoveries[i] == None
        else:
            recoveries[i] = recovery
    return recoveries


def sweep_rr(recovery_sweep):
    # DOES NOT REBUILD THE MODEL
    areas = np.zeros(np.shape(recovery_sweep))
    n = np.size(recovery_sweep)
    # fix work
    for i in range(n):  # for rebuilding with each run
        m = mvc_unit.build()
        mvc_unit.set_operating_conditions(m)
        # for i in range(n): # for not rebuilding the model
        m.fs.evaporator.area.unfix()
        vap_flow = recovery_sweep[i] * (
            m.fs.evaporator.properties_feed[0].flow_mass_phase_comp["Liq", "TDS"].value
            + m.fs.evaporator.properties_feed[0]
            .flow_mass_phase_comp["Liq", "H2O"]
            .value
        )
        m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].fix(
            vap_flow
        )
        mvc_unit.initialize_system(m)

        solver = get_solver()
        results = solver.solve(m, tee=False)

        if results.solver.termination_condition == TerminationCondition.infeasible:
            areas[i] == None
        else:
            areas[i] = m.fs.evaporator.area.value

    return areas


def plot(x, y, x_label, y_label):
    plt.figure()
    plt.plot(x, y, "-o")
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.show()


if __name__ == "__main__":
    main()
