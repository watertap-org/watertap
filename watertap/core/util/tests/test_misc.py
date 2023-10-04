#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import pytest

from pyomo.environ import ConcreteModel, Param, Var, units

from watertap.core.util.misc import is_constant_up_to_units


@pytest.mark.unit
def test_is_constant_up_to_units():

    # test float
    assert is_constant_up_to_units(42.0)

    # test int
    assert is_constant_up_to_units(42)

    # test unit-only expression
    assert is_constant_up_to_units(42.0 * units.m**2 / units.s)

    m = ConcreteModel()
    m.fixed_param = Param(initialize=42.0)

    # test non-mutable param
    assert is_constant_up_to_units(m.fixed_param)

    # test non-mutable param with unit-only expression
    assert is_constant_up_to_units(m.fixed_param * units.m**2 / units.s)

    m.fixed_param_2 = Param(initialize=6.28)

    # test expression of fixed params with units
    assert is_constant_up_to_units(
        (m.fixed_param**2) * m.fixed_param_2 * units.m**2
    )

    m.mutable_param = Param(initialize=42, mutable=True)

    # test mutable param
    assert not is_constant_up_to_units(m.mutable_param)

    m.mutable_param_units = Param(initialize=42, units=units.m**2 / units.s)
    # test mutable param with units specified
    assert not is_constant_up_to_units(m.mutable_param_units)

    m.fixed_variable = Var(initialize=42)
    m.fixed_variable.fix()
    m.fixed_variable_units = Var(initialize=42, units=units.m**2 / units.s)
    m.fixed_variable_units.fix()

    m.unfixed_variable = Var(initialize=6.28)
    m.unfixed_variable_units = Var(initialize=6.28, units=units.m**2 / units.s)

    # test variables
    assert not is_constant_up_to_units(m.fixed_variable)
    assert not is_constant_up_to_units(m.fixed_variable_units)
    assert not is_constant_up_to_units(m.unfixed_variable)
    assert not is_constant_up_to_units(m.unfixed_variable_units)

    # test combinations
    assert not is_constant_up_to_units(m.fixed_variable * m.fixed_param * units.m**2)
    assert not is_constant_up_to_units(m.fixed_variable * 42.0 * units.m**2)

    assert not is_constant_up_to_units(
        m.unfixed_variable * m.fixed_param * units.m**2
    )
    assert not is_constant_up_to_units(m.unfixed_variable * 42.0 * units.m**2)

    assert not is_constant_up_to_units(
        m.fixed_variable_units * m.fixed_param * units.m**2
    )
    assert not is_constant_up_to_units(m.fixed_variable_units * 42.0 * units.m**2)

    assert not is_constant_up_to_units(
        m.unfixed_variable_units * m.fixed_param * units.m**2
    )
    assert not is_constant_up_to_units(m.unfixed_variable_units * 42.0 * units.m**2)

    assert not is_constant_up_to_units(
        (1 / m.mutable_param_units) * m.fixed_param * units.m**2
    )
    assert not is_constant_up_to_units(
        (1 / m.mutable_param_units) * m.fixed_variable_units * 42.0 * units.m**2
    )
