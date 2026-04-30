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
    TerminationCondition,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom

import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from watertap.core.solvers import get_solver
from idaes.core import UnitModelCostingBlock

from watertap.property_models.unit_specific import cryst_prop_pack as props
from watertap.unit_models.crystallizer import Crystallization
from watertap.costing import WaterTAPCosting, CrystallizerCostType


def main():
    ##########################################
    # # Case 1: Fix crystallizer temperature
    ##########################################
    m = build1()
    add_costing(m)
    solve(m)

    ##########################################
    # # Case 2: Fix crystallizer yield
    ##########################################
    build2(m)
    solve(m)

    # ##########################################
    # # # Case 3: Fix crystallizer solids outlet
    # ##########################################
    build3(m)
    solve(m)

    #########################################
    # Case 4: Fix magma density
    #########################################
    build4(m)
    solve(m)

    #########################################
    # Case 5: Fix heat addition
    #########################################
    build5(m)
    solve(m)

    return m


def build1():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    # attach property package
    m.fs.properties = props.NaClParameterBlock()
    # build the unit model
    m.fs.crystallizer = Crystallization(property_package=m.fs.properties)

    # now specify the model
    print("DOF before specifying:", degrees_of_freedom(m.fs))

    # Specify the Feed
    m.fs.crystallizer.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(10.5119)
    m.fs.crystallizer.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(38.9326)
    m.fs.crystallizer.inlet.flow_mass_phase_comp[0, "Sol", "NaCl"].fix(1e-6)
    m.fs.crystallizer.inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(1e-6)
    m.fs.crystallizer.inlet.pressure[0].fix(101325)
    m.fs.crystallizer.inlet.temperature[0].fix(273.15 + 20)

    print("DOF after specifying feed:", degrees_of_freedom(m.fs))
    print("\n--- Case 1 ---")
    m.fs.crystallizer.temperature_operating.fix(273.15 + 55)
    m.fs.crystallizer.solids.flow_mass_phase_comp[0, "Sol", "NaCl"].fix(5.556)

    # Fix
    m.fs.crystallizer.crystal_growth_rate.fix()
    m.fs.crystallizer.souders_brown_constant.fix()
    m.fs.crystallizer.crystal_median_length.fix()

    # # Scaling
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "NaCl")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e0, index=("Vap", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Sol", "NaCl")
    )
    iscale.calculate_scaling_factors(m.fs)

    # Initialize crystallizer
    m.fs.crystallizer.initialize(outlvl=idaeslog.WARNING)
    assert_units_consistent(m)  # check that units are consistent

    return m


def add_costing(m):
    m.fs.costing = WaterTAPCosting()
    m.fs.crystallizer.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method_arguments={"cost_type": CrystallizerCostType.mass_basis},
    )
    m.fs.costing.cost_process()

    return m


def solve(m, solver=None, tee=False):
    assert degrees_of_freedom(m) == 0

    if solver is None:
        solver = get_solver()
    results = solver.solve(m, tee=tee, symbolic_solver_labels=True)

    assert results.solver.termination_condition == TerminationCondition.optimal
    m.fs.crystallizer.report()

    return results


def build2(m):
    print("\n--- Case 2 ---")
    m.fs.crystallizer.solids.flow_mass_phase_comp[0, "Sol", "NaCl"].unfix()
    m.fs.crystallizer.crystallization_yield["NaCl"].fix(0.7)

    return m


def build3(m):
    print("\n--- Case 3 ---")
    m.fs.crystallizer.crystallization_yield["NaCl"].unfix()
    m.fs.crystallizer.product_volumetric_solids_fraction.fix(0.1182)

    return m


def build4(m):
    print("\n--- Case 4 ---")
    m.fs.crystallizer.product_volumetric_solids_fraction.unfix()
    m.fs.crystallizer.dens_mass_magma.fix(250)

    return m


def build5(m):
    print("\n--- Case 5 ---")
    m.fs.crystallizer.dens_mass_magma.unfix()
    m.fs.crystallizer.work_mechanical[0].fix(55000)

    return m


if __name__ == "__main__":
    m = main()
