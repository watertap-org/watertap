##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Tests for 0D Nanofiltration unit model.
Author: Markus Drouven
"""

import pytest
from pyomo.environ import (ConcreteModel,
                           Constraint,
                           TerminationCondition,
                           SolverStatus,
                           units,
                           Var,
                           value)
from idaes.core import (FlowsheetBlock,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType)

# COMMENT: NEED TO UPDATE REFERENCE ONCE FINALIZED

from proteuslib.unit_models.nanofiltration import NF

# QUESTION: WHAT IS THE RIGHT SYNTAX HERE?

from proteuslib.property_models.NaCl_prop_pack import NaClParameterBlock
# import NaCl_prop_pack as props

from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_variables,
                                              number_total_constraints,
                                              number_unused_variables)
from idaes.core.util.testing import (get_default_solver,
                                     PhysicalParameterTestBlock,
                                     initialization_tester)
from pyomo.util.check_units import (assert_units_consistent,
                                    assert_units_equivalent)


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_default_solver()


# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    # QUESTION: WHY IS THIS GIVING ME AN ERROR?
    m.fs.properties = NaClParameterBlock()
    # m.fs.properties = props.NaClParameterBlock
    # m.fs.properties = PhysicalParameterTestBlock()
    
    m.fs.unit = NF(default={"property_package": m.fs.properties,"has_pressure_change": True})

    # Check unit config arguments
    assert len(m.fs.unit.config) == 8

    assert m.fs.unit.config.material_balance_type == \
        MaterialBalanceType.useDefault
    assert m.fs.unit.config.energy_balance_type == \
        EnergyBalanceType.useDefault
    assert m.fs.unit.config.momentum_balance_type == \
        MomentumBalanceType.pressureTotal
    assert m.fs.unit.config.has_pressure_change
    assert m.fs.unit.config.property_package is m.fs.properties



# -----------------------------------------------------------------------------
class TestNanofiltration(object):
    @pytest.fixture(scope="class")
    def defaultspec(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = NaClParameterBlock()
        
        m.fs.unit = NF(default={"property_package": m.fs.properties,"has_pressure_change": True})

        m.fs.unit.inlet.flow_mass_comp[0, 'NaCl'].fix(0.035)
        m.fs.unit.inlet.flow_mass_comp[0, 'H2O'].fix(0.965)
        
        m.fs.unit.inlet.pressure[0].fix(689476)
        m.fs.unit.inlet.temperature[0].fix(298.15)

        m.fs.unit.A.fix(3.77e-11)
        m.fs.unit.B.fix(4.724e-5)
        m.fs.unit.sigma.fix(0.28)
        m.fs.unit.deltaP.fix(-1e5)
        m.fs.unit.area.fix(2.5)

        m.fs.unit.permeate.pressure[0].fix(101325)

        return m


    @pytest.mark.unit
    def test_build(self, defaultspec):

        assert hasattr(defaultspec.fs.unit, "inlet")      
        assert hasattr(defaultspec.fs.unit.inlet, "flow_mass_comp")
        assert hasattr(defaultspec.fs.unit.inlet, "temperature")
        assert hasattr(defaultspec.fs.unit.inlet, "pressure")
        assert len(defaultspec.fs.unit.inlet.vars) == 3

        assert hasattr(defaultspec.fs.unit, "permeate")      
        assert hasattr(defaultspec.fs.unit.permeate, "flow_mass_comp")
        assert hasattr(defaultspec.fs.unit.permeate, "temperature")
        assert hasattr(defaultspec.fs.unit.permeate, "pressure")
        assert len(defaultspec.fs.unit.permeate.vars) == 3

        assert hasattr(defaultspec.fs.unit, "retentate")      
        assert hasattr(defaultspec.fs.unit.retentate, "flow_mass_comp")
        assert hasattr(defaultspec.fs.unit.retentate, "temperature")
        assert hasattr(defaultspec.fs.unit.retentate, "pressure")
        assert len(defaultspec.fs.unit.retentate.vars) == 3

        assert hasattr(defaultspec.fs.unit, "deltaP")

        assert number_variables(defaultspec) == 107
        assert number_total_constraints(defaultspec) == 42
        # QUESTION: WHY DO I HAVE SO MANY UNUSED VARIABLES IN HERE?
        assert number_unused_variables(defaultspec) == 13

    @pytest.mark.component
    def test_units(self, defaultspec):
        assert_units_consistent(defaultspec)
        assert_units_equivalent(defaultspec.fs.unit.deltaP[0], units.Pa)

    @pytest.mark.unit
    def test_dof(self, defaultspec):
        assert degrees_of_freedom(defaultspec) == 0

    @pytest.mark.component
    def test_initialize(self, defaultspec):
        initialization_tester(defaultspec)

    @pytest.mark.component
    def test_solve(self, defaultspec):
        results = solver.solve(defaultspec)

        # Check for optimal solution
        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, defaultspec):
        assert (pytest.approx(589476.0, abs=1e-2) ==
                value(defaultspec.fs.unit.retentate.pressure[0]))
        assert (pytest.approx(298.15, abs=1e-2) ==
                value(defaultspec.fs.unit.permeate.temperature[0]))
        assert (pytest.approx(0.03428220299, abs=1e-5) ==
                value(defaultspec.fs.unit.permeate.flow_mass_comp[0, "H2O"]))

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, defaultspec):
        assert (abs(value(defaultspec.fs.unit.inlet.flow_mass_comp[0, "NaCl"] - 
                          defaultspec.fs.unit.permeate.flow_mass_comp[0, "NaCl"] -
                          defaultspec.fs.unit.retentate.flow_mass_comp[0, "NaCl"]))
                <= 1e-6)



