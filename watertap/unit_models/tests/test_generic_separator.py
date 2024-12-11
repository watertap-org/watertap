#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
from pyomo.environ import ConcreteModel, value, assert_optimal_termination
import pytest
from watertap.core.solvers import get_solver

from idaes.core import FlowsheetBlock

import idaes.core.util.scaling as iscale


from watertap.unit_models.generic_separation import (
    GenericSeparation,
)
from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
    ActivityCoefficientModel,
    DensityCalculation,
    MaterialFlowBasis,
)

from idaes.core.util.model_statistics import degrees_of_freedom

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

# -----------------------------------------------------------------------------


def build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    mcas_props = {
        "activity_coefficient_model": ActivityCoefficientModel.ideal,
        "density_calculation": DensityCalculation.constant,
        "solute_list": ["TDS", "X"],
        "mw_data": {"H2O": 0.01801528, "TDS": 1, "X": 1},
        "charge": {"TDS": 0.0, "X": 0.0},
        "material_flow_basis": MaterialFlowBasis.mass,
    }

    m.fs.properties = MCASParameterBlock(**mcas_props)

    m.fs.unit = GenericSeparation(property_package=m.fs.properties)
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(1)
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(0.1)
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "X"].fix(0.01)
    m.fs.unit.additive_dose.fix(1)
    m.fs.unit.component_removal_percent["X"].fix(50)
    m.fs.unit.component_removal_percent["H2O"].fix(0)
    m.fs.unit.component_removal_percent["TDS"].fix(0)
    m.fs.unit.inlet.pressure[0].fix(1e5)
    m.fs.unit.inlet.temperature[0].fix(293.15)
    iscale.calculate_scaling_factors(m.fs.unit)

    return m


@pytest.mark.component
def test_solve():
    m = build()
    m.fs.unit.initialize()
    assert degrees_of_freedom(m) == 0
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)

    assert value(
        m.fs.unit.treated.flow_mass_phase_comp[0, "Liq", "H2O"]
    ) == pytest.approx(1, rel=1e-3)
    assert value(
        m.fs.unit.treated.flow_mass_phase_comp[0, "Liq", "TDS"]
    ) == pytest.approx(0.1, rel=1e-3)
    assert value(
        m.fs.unit.treated.flow_mass_phase_comp[0, "Liq", "X"]
    ) == pytest.approx(0.005, rel=1e-3)
