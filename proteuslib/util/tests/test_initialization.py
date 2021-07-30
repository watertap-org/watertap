###############################################################################
# ProteusLib Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/nawi-hub/proteuslib/"
#
###############################################################################

import pytest
from pyomo.environ import ConcreteModel

from idaes.core import (FlowsheetBlock,
                        ControlVolume1DBlock)
from proteuslib.util.initialization import check_solve, check_dof
from proteuslib.unit_models.reverse_osmosis_1D import (ReverseOsmosis1D,
                                                       ConcentrationPolarizationType,
                                                       MassTransferCoefficient,
                                                       PressureChangeType)
import proteuslib.property_models.NaCl_prop_pack as props
from idaes.core.util import get_solver
import idaes.logger as idaeslog


__author__ = "Adam Atia"

_log = idaeslog.getLogger(__name__)

# Set up solver
solver = get_solver()


@pytest.fixture(scope="class")
def model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.properties = props.NaClParameterBlock()

    m.fs.unit = ReverseOsmosis1D(default={
        "property_package": m.fs.properties,
        "has_pressure_change": True,
        "concentration_polarization_type": ConcentrationPolarizationType.calculated,
        "mass_transfer_coefficient": MassTransferCoefficient.calculated,
        "pressure_change_type": PressureChangeType.calculated,
        "transformation_scheme": "BACKWARD",
        "transformation_method": "dae.finite_difference",
        "finite_elements": 5,
        "has_full_reporting": True
    })


    return m

@pytest.mark.unit
def test_check_dof(model):
    m = model
    # check_dof should pass since fail_flag=False produces warning for DOF!=0
    check_dof(m, fail_flag=False)
    # Verify error since no variables were fixed
    with pytest.raises(ValueError, match="Non-zero degrees of freedom: Degrees of freedom on unknown = 11. "
                                         "Fix 11 more variable\(s\) or set keyword arg to ignore_dof=True"):
        check_dof(m, fail_flag=True)

    # fully specify system
    feed_flow_mass = 1000 / 3600
    feed_mass_frac_NaCl = 0.034283
    feed_pressure = 70e5

    feed_temperature = 273.15 + 25
    A = 4.2e-12
    B = 3.5e-8
    pressure_atmospheric = 1e5
    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl

    m.fs.unit.feed_inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl'].fix(
        feed_flow_mass * feed_mass_frac_NaCl)

    m.fs.unit.feed_inlet.flow_mass_phase_comp[0, 'Liq', 'H2O'].fix(
        feed_flow_mass * feed_mass_frac_H2O)

    m.fs.unit.feed_inlet.pressure[0].fix(feed_pressure)
    m.fs.unit.feed_inlet.temperature[0].fix(feed_temperature)
    m.fs.unit.A_comp.fix(A)
    m.fs.unit.B_comp.fix(B)
    m.fs.unit.permeate_outlet.pressure[0].fix(pressure_atmospheric)
    m.fs.unit.N_Re[0, 0].fix(400)
    m.fs.unit.recovery_mass_phase_comp[0, 'Liq', 'H2O'].fix(0.5)
    m.fs.unit.spacer_porosity.fix(0.97)
    m.fs.unit.channel_height.fix(0.001)

    # check should pass since DOF=0
    check_dof(m, fail_flag=True)


@pytest.mark.unit
def test_check_solve(model):
    solver.options = {'nlp_scaling_method': 'user-scaling'}

    results = solver.solve(model)
    # check_solve should pass since fail_flag=False and only warning will be produced
    check_solve(results, logger=_log, fail_flag=False)
    # Without calling calculate_scaling_factors() or initialization, expect the solve to fail and raise error
    with pytest.raises(ValueError, match="The solver failed to converge to an optimal solution. This suggests that the "
                                         "user provided infeasible inputs or that the model is poorly scaled."):
        check_solve(results, logger=_log, fail_flag=True)