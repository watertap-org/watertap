import pandas as pd
import numpy as np
import pytest
from pyomo.environ import (
    ConcreteModel,
    TerminationCondition,
    SolverStatus,
    value,
    Var,
)
from pyomo.network import Port
from idaes.core import FlowsheetBlock
from pyomo.util.check_units import assert_units_consistent
from watertap.unit_models.crystallizer import Crystallization
import watertap.property_models.cryst_prop_pack as props

from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    badly_scaled_var_generator,
)
from idaes.core import UnitModelCostingBlock

from watertap.costing import WaterTAPCosting, CrystallizerCostType

solver = get_solver()


def build_multi_effect_crystallizer_fs(
    m=None,
    operating_pressure_eff1=0.35,  # bar
    operating_pressure_eff2=0.208,  # bar
    operating_pressure_eff3=0.095,  # bar
    feed_mass_frac_NaCl=0.2126,
    feed_pressure=101325,  # Pa
    feed_temperature=273.15 + 20,  # K
    crystallizer_yield=0.5,
):
    """
    This flowsheet depicts a 3-effect crystallizer, with brine fed in parallel
    to each effect, and the operating pressure is specfied individually.
    """

    if m is None:
        m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.props = props.NaClParameterBlock()

    eff_1 = m.fs.effect_1 = Crystallization(property_package=m.fs.props)
    eff_2 = m.fs.effect_2 = Crystallization(property_package=m.fs.props)
    eff_3 = m.fs.effect_3 = Crystallization(property_package=m.fs.props)

    # Specify the feed properties
    feed_mass_frac_NaCl = 0.2126
    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
    eps = 1e-6

    eff_1.inlet.pressure[0].fix(feed_pressure)
    eff_1.inlet.temperature[0].fix(feed_temperature)
    eff_1.inlet.flow_mass_phase_comp[0, "Sol", "NaCl"].fix(eps)
    eff_1.inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(eps)

    eff_2.inlet.pressure[0].fix(feed_pressure)
    eff_2.inlet.temperature[0].fix(feed_temperature)
    eff_2.inlet.flow_mass_phase_comp[0, "Sol", "NaCl"].fix(eps)
    eff_2.inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(eps)

    eff_3.inlet.pressure[0].fix(feed_pressure)
    eff_3.inlet.temperature[0].fix(feed_temperature)
    eff_3.inlet.flow_mass_phase_comp[0, "Sol", "NaCl"].fix(eps)
    eff_3.inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(eps)

    # Set up for effect 1
    # Operating pressure: 0.35 bar
    c1_pressure = 101325 * operating_pressure_eff1
    feed_flow_mass = 1

    eff_1.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )
    eff_1.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )

    eff_1.pressure_operating.fix(c1_pressure)
    eff_1.crystallization_yield["NaCl"].fix(crystallizer_yield)

    eff_1.crystal_growth_rate.fix()
    eff_1.souders_brown_constant.fix()
    eff_1.crystal_median_length.fix()

    # Set up for effect 2
    # Same crystallizer_yield
    # Operating pressure: 0.208 bar
    c2_pressure = 101325 * operating_pressure_eff2

    eff_2.pressure_operating.fix(c2_pressure)
    eff_2.crystallization_yield["NaCl"].fix(crystallizer_yield)

    eff_2.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )
    eff_2.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )

    eff_2.crystal_growth_rate.fix()
    eff_2.souders_brown_constant.fix()
    eff_2.crystal_median_length.fix()

    # Set up for effect 2
    # Same crystallizer_yield
    # Operating pressure: 0.095 bar
    c3_pressure = 101325 * operating_pressure_eff3

    eff_3.pressure_operating.fix(c3_pressure)
    eff_3.crystallization_yield["NaCl"].fix(crystallizer_yield)

    eff_3.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )
    eff_3.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )

    eff_3.crystal_growth_rate.fix()
    eff_3.souders_brown_constant.fix()
    eff_3.crystal_median_length.fix()

    assert degrees_of_freedom(m) == 0

    # Scale
    m.fs.props.set_default_scaling("flow_mass_phase_comp", 1e-1, index=("Liq", "H2O"))
    m.fs.props.set_default_scaling("flow_mass_phase_comp", 1e-1, index=("Liq", "NaCl"))
    m.fs.props.set_default_scaling("flow_mass_phase_comp", 1e-1, index=("Vap", "H2O"))
    m.fs.props.set_default_scaling("flow_mass_phase_comp", 1e-1, index=("Sol", "NaCl"))
    calculate_scaling_factors(m)

    return m


def initialize_and_unfix_dof(m):
    eff_1 = m.fs.effect_1
    eff_2 = m.fs.effect_2
    eff_3 = m.fs.effect_3

    # Initialize
    eff_1.initialize()
    eff_2.initialize()
    eff_3.initialize()

    # Unfix dof
    brine_salinity = eff_1.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"].value

    eff_2.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].unfix()
    eff_2.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].unfix()
    eff_2.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"].fix(brine_salinity)

    eff_3.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].unfix()
    eff_3.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].unfix()
    eff_3.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"].fix(brine_salinity)

    # Energy is provided from the previous effect
    @m.Constraint(doc="Energy supplied to the 2nd effect")
    def eqn_energy_from_eff1(b):
        return b.fs.effect_2.work_mechanical[0] == b.fs.effect_1.energy_from_vapor

    @m.Constraint(doc="Energy supplied to the 3rd effect")
    def eqn_energy_from_eff2(b):
        return b.fs.effect_3.work_mechanical[0] == b.fs.effect_2.energy_from_vapor

    # Solve
    results = solver.solve(m)
    assert results.solver.termination_condition == TerminationCondition.optimal


def get_model_performance(m):
    # Print result
    effs = [m.fs.effect_1, m.fs.effect_2, m.fs.effect_3]
    effect_names = ["Effect 1", "Effect 2", "Effect 3"]
    feed_salinities = [
        i.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"].value for i in effs
    ]
    feed_flow_rates = [
        sum(
            i.properties_in[0].flow_mass_phase_comp["Liq", j].value
            for j in ["H2O", "NaCl"]
        )
        for i in effs
    ]
    feed_vol_flow_rates = [i.properties_in[0].flow_vol_phase["Liq"].value for i in effs]
    temp_operating = [i.temperature_operating.value - 273.15 for i in effs]
    temp_vapor_cond = [
        i.properties_pure_water[0].temperature.value - 273.15 for i in effs
    ]
    p_operating = [i.pressure_operating.value / 101325 for i in effs]
    water_prod = [
        i.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].value for i in effs
    ]
    solid_prod = [
        i.properties_solids[0].flow_mass_phase_comp["Sol", "NaCl"].value for i in effs
    ]
    liquid_prod = [
        sum(
            i.properties_out[0].flow_mass_phase_comp["Liq", j].value
            for j in ["H2O", "NaCl"]
        )
        for i in effs
    ]
    power_required = [i.work_mechanical[0].value for i in effs]
    power_provided = [i.energy_from_vapor.value for i in effs]
    STEC = [
        i.work_mechanical[0].value
        / i.properties_in[0].flow_vol_phase["Liq"].value
        / 3600
        for i in effs
    ]

    overall_STEC = (
        m.fs.effect_1.work_mechanical[0].value
        / sum(i.properties_in[0].flow_vol_phase["Liq"].value for i in effs)
        / 3600
    )

    model_output = np.array(
        [
            feed_flow_rates,
            feed_vol_flow_rates,
            feed_salinities,
            temp_operating,
            temp_vapor_cond,
            p_operating,
            water_prod,
            solid_prod,
            liquid_prod,
            power_required,
            power_provided,
            STEC,
        ]
    )

    data_table = pd.DataFrame(
        data=model_output,
        columns=effect_names,
        index=[
            "Feed mass flow rate (kg/s)",
            "Feed volumetric flow rate (m3/s)",
            "Feed salinities (g/L)",
            "Operating temperature (C)",
            "Vapor condensation temperature (C)",
            "Operating pressure (bar)",
            "Water production (kg/s)",
            "Solid production (kg/s)",
            "Liquid waste (kg/s)",
            "Thermal energy requirement (kW)",
            "Thermal energy available from vapor (kW)",
            "STEC (kWh/m3 feed)",
        ],
    )

    overall_performance = {
        "Feed brine salinity (g/L)": m.fs.effect_1.properties_in[0]
        .conc_mass_phase_comp["Liq", "NaCl"]
        .value,
        "Total brine disposed (kg/s)": sum(feed_flow_rates),
        "Total water production (kg/s)": sum(water_prod),
        "Total solids collected (kg/s)": sum(solid_prod),
        "Total waste water remained (kg/s)": sum(liquid_prod),
        "Initial thermal energy consumption (kW)": m.fs.effect_1.work_mechanical[
            0
        ].value,
        "Overall STEC (kWh/m3 feed)": overall_STEC,
    }

    return data_table, overall_performance


if __name__ == "__main__":
    test_case = build_multi_effect_crystallizer_fs(
        operating_pressure_eff1=0.35,
        operating_pressure_eff2=0.208,
        operating_pressure_eff3=0.095,
        feed_mass_frac_NaCl=0.21,
        feed_pressure=101325,
        feed_temperature=273.15 + 20,
        crystallizer_yield=0.9,
    )

    initialize_and_unfix_dof(test_case)
    data_table, overall_performance = get_model_performance(test_case)

    # Print model results
    pd.set_option("display.precision", 3)
    print(data_table)
    print("")
    print("System overall performance")
    for (key, value) in overall_performance.items():
        print(key, round(value, 2))
