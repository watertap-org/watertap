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
from pyomo.core.base.units_container import InconsistentUnitsError
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
    # Tests for performance indices with flow_basis specified and inferred from flow_rate units
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
    m.fs.costing.add_process_throughput(
        sum(
            m.fs.product.properties[0].flow_mass_phase_comp["Liq", comp]
            for comp in m.fs.properties.component_list
        ),
        flow_basis="mass",
        name="annual_product_generation",
    )
    m.fs.costing.add_specific_electrical_carbon_intensity(
        sum(
            m.fs.product.properties[0].flow_mass_phase_comp["Liq", comp]
            for comp in m.fs.properties.component_list
        ),
        flow_basis="mass",
        name="specific_electrical_carbon_intensity_with_product",
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
def test_units_based_basis_validation():
    m = lsrro.build()

    comp = next(iter(m.fs.properties.component_list))
    mass_flow = m.fs.product.properties[0].flow_mass_phase_comp["Liq", comp]

    # Incompatible explicit flow_basis values now fail at unit conversion time.
    with pytest.raises(InconsistentUnitsError):
        m.fs.costing.add_levelized_cost(
            m.fs.product.properties[0].flow_vol,
            flow_basis="mass",
            name="LCOW_mass_basis",
        )

    # Incompatible output units are rejected.
    with pytest.raises(
        ValueError,
        match=r"Could not infer flow basis from units",
    ):
        m.fs.costing.add_levelized_cost(
            mass_flow,
            flow_basis_units=pyo.units.m,
            name="LCOP_invalid_units",
        )

    # Flow units that do not match the requested denominator conversion fail.
    with pytest.raises(
        InconsistentUnitsError,
    ):
        m.fs.costing.add_specific_energy_consumption(
            flow_basis="volumetric",
            flow_rate=mass_flow,
            name="specific_energy_consumption_mismatch",
        )

    # Matching units continue to work.
    m.fs.costing.add_process_throughput(
        mass_flow,
        flow_basis="mass",
        name="annual_throughput_ok",
    )


@pytest.mark.component
def test_flow_basis_units_inference_and_validation():
    m = lsrro.build()

    mass_flow = sum(
        m.fs.product.properties[0].flow_mass_phase_comp["Liq", comp]
        for comp in m.fs.properties.component_list
    )

    # infer basis from flow_rate units when neither flow_basis nor flow_basis_units is provided
    m.fs.costing.add_process_throughput(
        mass_flow,
        name="annual_throughput_inferred_mass",
    )

    # infer basis from explicit flow_basis_units
    m.fs.costing.add_specific_energy_consumption(
        mass_flow,
        flow_basis_units=pyo.units.kg,
        name="specific_energy_consumption_flow_basis_units_mass",
    )

    # flow_basis_units should be behaviorally equivalent to the matching flow_basis
    m.fs.costing.add_process_throughput(
        mass_flow,
        flow_basis="mass",
        name="annual_throughput_mass_basis",
    )
    m.fs.costing.add_process_throughput(
        mass_flow,
        flow_basis_units=pyo.units.kg,
        name="annual_throughput_flow_basis_units_mass",
    )

    # flow_basis must agree with flow_basis_units when both are provided
    with pytest.raises(
        ValueError,
        match=r"flow_basis 'volumetric' is inconsistent with flow_basis_units",
    ):
        m.fs.costing.add_levelized_cost(
            mass_flow,
            flow_basis="volumetric",
            flow_basis_units=pyo.units.kg,
            name="LCOP_invalid_basis_units",
        )

    assert_units_consistent(m)
    assert pytest.approx(
        pyo.value(m.fs.costing.annual_throughput_mass_basis)
    ) == pyo.value(m.fs.costing.annual_throughput_flow_basis_units_mass)


@pytest.mark.component
def test_flow_basis_units_accepts_convertible_units():
    m = lsrro.build()

    mass_flow = sum(
        m.fs.product.properties[0].flow_mass_phase_comp["Liq", comp]
        for comp in m.fs.properties.component_list
    )
    volumetric_flow = m.fs.product.properties[0].flow_vol

    m.fs.costing.add_process_throughput(
        volumetric_flow,
        flow_basis_units=pyo.units.m**3,
        name="annual_throughput_m3",
    )
    m.fs.costing.add_process_throughput(
        volumetric_flow,
        flow_basis_units=pyo.units.L,
        name="annual_throughput_L",
    )

    m.fs.costing.add_process_throughput(
        mass_flow,
        flow_basis_units=pyo.units.kg,
        name="annual_throughput_kg",
    )
    m.fs.costing.add_process_throughput(
        mass_flow,
        flow_basis_units=pyo.units.g,
        name="annual_throughput_g",
    )

    assert_units_consistent(m)
    assert pytest.approx(
        pyo.value(m.fs.costing.annual_throughput_L), rel=1e-10
    ) == 1000 * pyo.value(m.fs.costing.annual_throughput_m3)
    assert pytest.approx(
        pyo.value(m.fs.costing.annual_throughput_g), rel=1e-10
    ) == 1000 * pyo.value(m.fs.costing.annual_throughput_kg)


# Test for different period units for add_process_throughput
@pytest.mark.component
def test_flow_basis_units_period_units():
    m = lsrro.build()

    mass_flow = sum(
        m.fs.product.properties[0].flow_mass_phase_comp["Liq", comp]
        for comp in m.fs.properties.component_list
    )

    m.fs.costing.add_process_throughput(
        mass_flow,
        flow_basis_units=pyo.units.kg,
        period=pyo.units.day,
        name="annual_throughput_kg_per_day",
    )
    m.fs.costing.add_process_throughput(
        mass_flow,
        flow_basis_units=pyo.units.kg,
        period=pyo.units.year,
        name="annual_throughput_kg_per_year",
    )

    days_per_year = pyo.value(
        pyo.units.convert(
            1 * pyo.units.year,
            to_units=pyo.units.day,
        )
    )

    assert pytest.approx(
        pyo.value(m.fs.costing.annual_throughput_kg_per_year),
        rel=1e-10,
    ) == days_per_year * pyo.value(m.fs.costing.annual_throughput_kg_per_day)
