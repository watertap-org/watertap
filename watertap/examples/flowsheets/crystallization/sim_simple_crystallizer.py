from pyomo.environ import (
    ConcreteModel,
    SolverFactory,
    TerminationCondition,
    value,
    Constraint,
    Var,
    Objective,
    Expression,
)
from pyomo.environ import units as pyunits
from pyomo.util.check_units import (
    assert_units_consistent,
    assert_units_equivalent,
    check_units_equivalent,
)
import pyomo.util.infeasible as infeas
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_activated_constraints,
    number_unfixed_variables_in_activated_equalities,
    number_activated_equalities,
    number_unused_variables,
)

import idaes.core.util.model_statistics as stats
from idaes.core.util.constants import Constants
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.core.solvers import get_solver
from idaes.core import UnitModelCostingBlock

from watertap.property_models import cryst_prop_pack as props
from watertap.unit_models.crystallizer import Crystallization
from watertap.costing import WaterTAPCosting
from watertap.costing.watertap_costing_package import CrystallizerCostType

from io import StringIO
from pyomo.util.infeasible import (
    log_active_constraints,
    log_close_to_bounds,
    log_infeasible_bounds,
    log_infeasible_constraints,
)
from pyomo.common.log import LoggingIntercept
import logging


if __name__ == "__main__":

    # create model, flowsheet
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    # attach property package
    m.fs.properties = props.NaClParameterBlock()
    m.fs.costing = WaterTAPCosting()
    # build the unit model
    m.fs.crystallizer = Crystallization(property_package=m.fs.properties)

    # now specify the model
    print("DOF before specifying:", degrees_of_freedom(m.fs))

    # Specify the Feed
    m.fs.crystallizer.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(10.5119)
    m.fs.crystallizer.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(38.9326)
    # m.fs.crystallizer.properties_in[0].flow_vol_phase['Liq'].fix(150/3600)
    # m.fs.crystallizer.properties_in[0].mass_frac_phase_comp['Liq', 'NaCl'].fix(0.2126)
    m.fs.crystallizer.inlet.flow_mass_phase_comp[0, "Sol", "NaCl"].fix(1e-6)
    m.fs.crystallizer.inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(1e-6)
    m.fs.crystallizer.inlet.pressure[0].fix(101325)
    m.fs.crystallizer.inlet.temperature[0].fix(273.15 + 20)
    print("DOF after specifying feed:", degrees_of_freedom(m.fs))  # Should be 2

    ##########################################
    # # Case 1: Fix crystallizer temperature
    ##########################################
    print("\n--- Case 1 ---")
    m.fs.crystallizer.temperature_operating.fix(273.15 + 55)
    #  m.fs.crystallizer.pressure_operating.fix(10000)
    m.fs.crystallizer.solids.flow_mass_phase_comp[0, "Sol", "NaCl"].fix(5.556)
    # m.fs.crystallizer.vapor.flow_mass_phase_comp[0, 'Vap', 'H2O'].fix(20.56)

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
    m.fs.crystallizer.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method_arguments={"cost_type": CrystallizerCostType.mass_basis},
    )
    m.fs.costing.cost_process()

    #  m.fs.crystallizer.k_param = 0.06
    # solving
    m.fs.crystallizer.initialize(outlvl=idaeslog.DEBUG)
    assert_units_consistent(m)  # check that units are consistent
    assert (
        degrees_of_freedom(m) == 0
    )  # check that the degrees of freedom are what we expect
    solver = get_solver()
    # solver = SolverFactory('ipopt')
    # solver.options = {'tol': 1e-8, 'nlp_scaling_method': 'user-scaling', 'halt_on_ampl_error': 'yes'}
    results = solver.solve(m, tee=True, symbolic_solver_labels=True)
    output = StringIO()
    with LoggingIntercept(output, "pyomo.util.infeasible", logging.INFO):
        log_infeasible_constraints(m)
    print(output.getvalue().splitlines())
    assert results.solver.termination_condition == TerminationCondition.optimal
    # m.fs.crystallizer.display()

    m.fs.crystallizer.report()
    m.fs.crystallizer.costing.capital_cost.pprint()

    #  assert False
    ##########################################
    # # Case 2: Fix crystallizer yield
    ##########################################
    print("\n--- Case 2 ---")
    # m.fs.crystallizer.T_crystallization.unfix()
    m.fs.crystallizer.solids.flow_mass_phase_comp[0, "Sol", "NaCl"].unfix()
    m.fs.crystallizer.crystallization_yield["NaCl"].fix(0.7)

    # solving
    assert_units_consistent(m)  # check that units are consistent
    assert (
        degrees_of_freedom(m) == 0
    )  # check that the degrees of freedom are what we expect
    # solver = SolverFactory('ipopt')
    # solver.options = {'tol': 1e-8, 'nlp_scaling_method': 'user-scaling'}
    solver = get_solver()
    results = solver.solve(m, tee=True)
    output = StringIO()
    with LoggingIntercept(output, "pyomo.util.infeasible", logging.INFO):
        log_infeasible_constraints(m)
    print(output.getvalue().splitlines())
    assert results.solver.termination_condition == TerminationCondition.optimal
    # m.fs.crystallizer.display()

    m.fs.crystallizer.report()

    # ##########################################
    # # # Case 3: Fix crystallizer solids outlet
    # ##########################################
    print("\n--- Case 3 ---")
    m.fs.crystallizer.crystallization_yield["NaCl"].unfix()
    m.fs.crystallizer.product_volumetric_solids_fraction.fix(0.1182)

    # solving
    assert_units_consistent(m)  # check that units are consistent
    assert (
        degrees_of_freedom(m) == 0
    )  # check that the degrees of freedom are what we expect
    solver = get_solver()
    # solver = SolverFactory('ipopt')
    # solver.options = {'tol': 1e-8, 'nlp_scaling_method': 'user-scaling'}
    results = solver.solve(m, tee=True)
    output = StringIO()
    with LoggingIntercept(output, "pyomo.util.infeasible", logging.INFO):
        log_infeasible_constraints(m)
    print(output.getvalue().splitlines())
    assert results.solver.termination_condition == TerminationCondition.optimal
    # m.fs.crystallizer.display()

    m.fs.crystallizer.report()

    #########################################
    # Case 4: Fix magma density
    #########################################
    print("\n--- Case 4 ---")
    m.fs.crystallizer.product_volumetric_solids_fraction.unfix()
    m.fs.crystallizer.dens_mass_magma.fix(250)

    # solving
    assert_units_consistent(m)  # check that units are consistent
    assert (
        degrees_of_freedom(m) == 0
    )  # check that the degrees of freedom are what we expect
    solver = get_solver()
    # solver = SolverFactory('ipopt')
    # solver.options = {'tol': 1e-8, 'nlp_scaling_method': 'user-scaling'}
    results = solver.solve(m, tee=True)
    output = StringIO()
    with LoggingIntercept(output, "pyomo.util.infeasible", logging.INFO):
        log_infeasible_constraints(m)
    print(output.getvalue().splitlines())
    assert results.solver.termination_condition == TerminationCondition.optimal
    # m.fs.crystallizer.display()

    m.fs.crystallizer.report()

    #########################################
    # Case 5: Fix heat addition
    #########################################
    print("\n--- Case 5 ---")
    m.fs.crystallizer.dens_mass_magma.unfix()
    m.fs.crystallizer.work_mechanical[0].fix(55000)

    # solving
    assert_units_consistent(m)  # check that units are consistent
    assert (
        degrees_of_freedom(m) == 0
    )  # check that the degrees of freedom are what we expect
    solver = get_solver()
    # solver = SolverFactory('ipopt')
    # solver.options = {'tol': 1e-8, 'nlp_scaling_method': 'user-scaling'}
    results = solver.solve(m, tee=True)
    output = StringIO()
    with LoggingIntercept(output, "pyomo.util.infeasible", logging.INFO):
        log_infeasible_constraints(m)
    print(output.getvalue().splitlines())
    assert results.solver.termination_condition == TerminationCondition.optimal
    # m.fs.crystallizer.display()

    m.fs.crystallizer.report()

    # print('number_variables:', number_variables(m.fs))
    # print('number_unused_variables:', number_unused_variables(m.fs))
    # print('number_total_constraints:', number_total_constraints(m.fs))
    # print('number_activated_constraints:', number_activated_constraints(m.fs))
    # print('number_activated_equalities', number_activated_equalities(m.fs))
