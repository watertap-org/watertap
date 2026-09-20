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
    units as pyunits,
    TransformationFactory,
    assert_optimal_termination,
    value,
)
from pyomo.network import Port
from idaes.core.solvers import petsc

from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    StateBlock,
)
from idaes.core.util.model_statistics import (
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
import idaes.core.util.scaling as iscale
from watertap.core import (
    MembraneChannel0DBlock,
    FrictionFactor,
    ModuleType,
)
from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
)
from watertap.unit_models.reverse_osmosis_0D import (
    ReverseOsmosis0D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
from watertap.core.solvers import get_solver
from watertap.unit_models.reverse_osmosis_base import TransportModel

import watertap.property_models.NaCl_prop_pack as props

from watertap.unit_models.tests.unit_test_harness import UnitTestHarness
import pytest

from idaes.core.solvers import petsc
import numpy as np

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)

m.fs.properties = props.NaClParameterBlock()

m.fs.unit = ReverseOsmosis0D(
    property_package=m.fs.properties,
    has_pressure_change=True,
    concentration_polarization_type=ConcentrationPolarizationType.calculated,
    mass_transfer_coefficient=MassTransferCoefficient.fixed,
    pressure_change_type=PressureChangeType.calculated,
    transport_model=TransportModel.SKK,
)

m.fs.unit.list_vars_to_fix()

# fully specify system
feed_flow_mass = 1
feed_mass_frac_NaCl = 0.035
feed_pressure = 50e5
feed_temperature = 273.15 + 25
# membrane_pressure_drop = 3e5
# dP_dx = 3e5 / 6
membrane_area = 50
membrane_length = 6
# A = 4.2e-12
# B = 3.5e-8
alpha = 0.8
sigma = 0.95
pressure_atmospheric = 101325
# concentration_polarization_modulus = 1.1
K = 2.6259644014323275e-05

feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
    feed_flow_mass * feed_mass_frac_NaCl
)
m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
    feed_flow_mass * feed_mass_frac_H2O
)
m.fs.unit.inlet.pressure[0].fix(feed_pressure)
m.fs.unit.inlet.temperature[0].fix(feed_temperature)
# m.fs.unit.deltaP.fix(-membrane_pressure_drop)
# m.fs.unit.feed_side.dP_dx.fix(-dP_dx)
m.fs.unit.area.fix(membrane_area)
m.fs.unit.length.fix(membrane_length)
# m.fs.unit.A_comp.fix(A)
# m.fs.unit.B_comp.fix(B)
m.fs.unit.alpha.fix(alpha)
m.fs.unit.reflect_coeff.fix(sigma)
m.fs.unit.permeate.pressure[0].fix(pressure_atmospheric)
# m.fs.unit.feed_side.cp_modulus.fix(concentration_polarization_modulus)
m.fs.unit.feed_side.K.fix(K)
m.fs.unit.feed_side.channel_height.fix(0.002)
m.fs.unit.feed_side.spacer_porosity.fix(0.75)

# Set scaling factors for badly scaled variables
m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1e2, index=("Liq", "NaCl"))
iscale.calculate_scaling_factors(m.fs.unit)

m.fs.unit.list_vars_to_fix()


if __name__ == "__main__":
    results = solver.solve(m)
    print(value(m.fs.unit.length))
    print(value(m.fs.unit.width))
    print(value(m.fs.unit.feed_side.K[0, 0, "NaCl"]))  # 2.6259644014323275e-05
