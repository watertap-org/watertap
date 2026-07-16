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

import pytest

from pyomo.util.check_units import assert_units_consistent
import pyomo.environ as pyo
import idaes.core as idc

from watertap.costing.watertap_costing_package import WaterTAPCosting
import watertap.flowsheets.lsrro.lsrro as lsrro


@pytest.mark.component
def test_watertap_costing_package():
    m = pyo.ConcreteModel()
    m.fs = idc.FlowsheetBlock(dynamic=False)

    m.fs.costing = WaterTAPCosting()

    m.fs.electricity = pyo.Var(units=pyo.units.kW, initialize=1)

    m.fs.costing.cost_flow(m.fs.electricity, "electricity")

    assert "foo" not in m.fs.costing.flow_types
    with pytest.raises(
        ValueError,
        match="foo is not a recognized flow type. Please check "
        "your spelling and that the flow type has been registered with"
        " the FlowsheetCostingBlock.",
    ):
        m.fs.costing.cost_flow(m.fs.electricity, "foo")

    m.fs.costing.foo_cost = foo_cost = pyo.Var(
        initialize=42,
        doc="foo",
        units=pyo.units.USD_2020 / pyo.units.m / pyo.units.second,
    )

    m.fs.costing.register_flow_type("foo", m.fs.costing.foo_cost)

    # make sure the component was not replaced
    # by register_defined_flow
    assert foo_cost is m.fs.costing.foo_cost

    m.fs.foo = pyo.Var(units=pyo.units.m, initialize=10)

    m.fs.costing.cost_flow(m.fs.foo, "foo")

    m.fs.costing.bar_base_cost = pyo.Var(
        initialize=0.42,
        doc="bar",
        units=pyo.units.USD_2020 / pyo.units.g / pyo.units.hour,
    )
    m.fs.costing.bar_purity = pyo.Param(
        initialize=0.50, doc="bar purity", units=pyo.units.dimensionless
    )

    m.fs.costing.register_flow_type(
        "bar", m.fs.costing.bar_base_cost * m.fs.costing.bar_purity
    )

    bar_cost = m.fs.costing.bar_cost
    assert isinstance(bar_cost, pyo.Expression)
    assert pyo.value(bar_cost) == 0.21

    m.fs.costing.bar_base_cost.value = 1.5
    assert pyo.value(bar_cost) == 0.75

    m.fs.costing.baz_cost = pyo.Var(initialize=5)

    with pytest.raises(
        RuntimeError,
        match="Component baz_cost already exists on fs.costing but is not 42",
    ):
        m.fs.costing.register_flow_type(
            "baz", 42 * pyo.units.USD_2020 / pyo.units.m**2 / pyo.units.day
        )

    m.fs.costing.flow_types.remove("baz")

    m.fs.costing.register_flow_type(
        "ham", 42 * pyo.units.USD_2021 / pyo.units.kg / pyo.units.minute
    )

    assert isinstance(m.fs.costing.ham_cost, pyo.Var)

    m.fs.costing.cost_process()
    # no error, wacc, plant_lifetime fixed
    m.fs.costing.initialize()

    m.fs.costing.capital_recovery_factor.fix()
    with pytest.raises(
        RuntimeError,
        match="Exactly two of the variables fs.costing.plant_lifetime, "
        "fs.costing.wacc, fs.costing.capital_recovery_factor should be "
        "fixed and the other unfixed.",
    ):
        # error, capital_recovery_factor,  wacc, plant_lifetime all fixed
        m.fs.costing.initialize()
    m.fs.costing.wacc.unfix()

    # no error
    m.fs.costing.initialize()

    m.fs.costing.capital_recovery_factor.unfix()
    with pytest.raises(
        RuntimeError,
        match="Exactly two of the variables fs.costing.plant_lifetime, "
        "fs.costing.wacc, fs.costing.capital_recovery_factor should be "
        "fixed and the other unfixed.",
    ):
        # error, capital_recovery_factor, wacc, unfixed
        m.fs.costing.initialize()

    m.fs.costing.plant_lifetime.unfix()
    with pytest.raises(
        RuntimeError,
        match="Exactly two of the variables fs.costing.plant_lifetime, "
        "fs.costing.wacc, fs.costing.capital_recovery_factor should be "
        "fixed and the other unfixed.",
    ):
        # error, capital_recovery_factor, wacc, and plant_lifetime unfixed
        m.fs.costing.initialize()

    m.fs.costing.wacc.fix()
    m.fs.costing.capital_recovery_factor.fix()
    # no error, wacc, capital_recovery_factor fixed
    m.fs.costing.initialize()


@pytest.mark.component
def test_breakdowns():
    m = lsrro.build()

    m.fs.BoosterPumps[:].control_volume.work[0.0].value = 42e6
    m.fs.EnergyRecoveryDevices[:].control_volume.work[0.0].value = -42e6

    m.fs.costing.add_specific_electrical_carbon_intensity(
        m.fs.product.properties[0].flow_vol
    )
    # Tests for performance indices introduced from valorization_costing_block
    m.fs.costing.add_levelized_cost(
        sum(
            m.fs.product.properties[0].flow_mass_phase_comp["Liq", comp]
            for comp in m.fs.properties.component_list
        ),
        flow_basis="mass",
        name="LCOP",
    )
    m.fs.costing.add_levelized_cost(
        m.fs.product.properties[0].flow_vol,
        flow_basis="volumetric",
        name="LCOT",
    )
    m.fs.costing.add_specific_energy_consumption(
        sum(
            m.fs.product.properties[0].flow_mass_phase_comp["Liq", comp]
            for comp in m.fs.properties.component_list
        ),
        flow_basis="mass",
        name="specific_energy_consumption_with_product",
    )
    m.fs.costing.add_annual_total(
        sum(
            m.fs.product.properties[0].flow_mass_phase_comp["Liq", comp]
            for comp in m.fs.properties.component_list
        ),
        flow_basis="mass",
        name="annual_product_generation",
    )

    assert_units_consistent(m)

    m.fs.costing.initialize()

    total_LCOW = pyo.value(m.fs.costing.LCOW)

    summed_aggregates = pyo.value(
        sum(m.fs.costing.LCOW_aggregate_direct_capex.values())
        + sum(m.fs.costing.LCOW_aggregate_indirect_capex.values())
        + sum(m.fs.costing.LCOW_aggregate_fixed_opex.values())
        + sum(m.fs.costing.LCOW_aggregate_variable_opex.values())
        # electricity is counted both in the aggregate variable opex
        # per unit and per flow, so it is double-counted in this sum
        - m.fs.costing.LCOW_aggregate_variable_opex["electricity"]
    )
    assert pytest.approx(total_LCOW) == summed_aggregates

    summed_components = pyo.value(
        sum(m.fs.costing.LCOW_component_direct_capex.values())
        + sum(m.fs.costing.LCOW_component_indirect_capex.values())
        + sum(m.fs.costing.LCOW_component_fixed_opex.values())
        + sum(m.fs.costing.LCOW_component_variable_opex.values())
    )
    assert pytest.approx(total_LCOW) == summed_components

    sec = pyo.value(m.fs.costing.specific_energy_consumption)
    summed_sec = pyo.value(
        sum(m.fs.costing.specific_energy_consumption_component.values())
    )
    assert pytest.approx(sec) == summed_sec

    seci = pyo.value(m.fs.costing.specific_electrical_carbon_intensity)
    summed_seci = pyo.value(
        sum(m.fs.costing.specific_electrical_carbon_intensity_component.values())
    )
    assert pytest.approx(seci) == summed_seci


@pytest.mark.component
def test_flow_basis_mismatch():
    m = lsrro.build()

    comp = next(iter(m.fs.properties.component_list))

    # flow_vol supplied but flow_basis says "mass"
    with pytest.raises(
        ValueError,
        match=r"flow_basis was set to 'mass', but the supplied flow_rate "
        r".*flow_vol.* appears to be a 'volumetric' flow \(matched on "
        r"'flow_vol'\)\. Please check that flow_basis matches the flow_rate "
        r"provided\.",
    ):
        m.fs.costing.add_levelized_cost(
            m.fs.product.properties[0].flow_vol,
            flow_basis="mass",
            name="LCOW_mismatch",
        )

    # flow_mass_phase_comp supplied but flow_basis says "volumetric"
    with pytest.raises(
        ValueError,
        match=r"flow_basis was set to 'volumetric', but the supplied flow_rate "
        r".*flow_mass_phase_comp.* appears to be a 'mass' flow \(matched on "
        r"'flow_mass'\)\. Please check that flow_basis matches the flow_rate "
        r"provided\.",
    ):
        m.fs.costing.add_specific_energy_consumption(
            m.fs.product.properties[0].flow_mass_phase_comp["Liq", comp],
            flow_basis="volumetric",
            name="specific_energy_consumption_mismatch",
        )

    # flow_vol supplied but flow_basis says "energy"
    with pytest.raises(
        ValueError,
        match=r"flow_basis was set to 'energy', but the supplied flow_rate "
        r".*flow_vol.* appears to be a 'volumetric' flow \(matched on "
        r"'flow_vol'\)\. Please check that flow_basis matches the flow_rate "
        r"provided\.",
    ):
        m.fs.costing.add_annual_total(
            m.fs.product.properties[0].flow_vol,
            flow_basis="energy",
            name="annual_total_mismatch",
        )
