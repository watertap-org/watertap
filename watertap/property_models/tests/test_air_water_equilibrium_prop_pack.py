#################################################################################
# WaterTAP Copyright (c) 2020-2025, The Regents of the University of California,
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
from pyomo.environ import (
    ConcreteModel,
    value,
    Var,
    Param,
    Set,
    assert_optimal_termination,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.util.scaling import calculate_scaling_factors
from idaes.core.util.exceptions import ConfigurationError

from watertap.core.util.initialization import check_dof
from watertap.core.solvers import get_solver

from watertap.property_models.air_water_equilibrium_prop_pack import (
    AirWaterEq,
    AirWaterEqStateBlock,
    MolarVolumeCalculation,
    LiqDiffusivityCalculation,
    VapDiffusivityCalculation,
    SaturationVaporPressureCalculation,
    VaporPressureCalculation,
    RelativeHumidityCalculation,
    LatentHeatVaporizationCalculation,
    SpecificHeatWaterCalculation,
    DensityCalculation,
)

solver = get_solver()


@pytest.mark.unit
def test_rh_pressure_vap_config():

    props = dict(
        relative_humidity_calculation=RelativeHumidityCalculation.FromVaporPressureRatio,
        vapor_pressure_calculation=VaporPressureCalculation.FromRelativeHumidity,
    )
    m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)

    err_msg = "\nUsing VaporPressureCalculation.FromRelativeHumidity and\n"
    err_msg += "RelativeHumidityCalculation.FromVaporPressureRatio creates\n"
    err_msg += "a redundant constraint. Choose a different option for\n"
    err_msg += "either configuration argument."

    with pytest.raises(ConfigurationError, match=err_msg):
        m.fs.properties = AirWaterEq(**props)


@pytest.fixture(scope="module")
def m1():
    """
    Test NMSU case study 1 results
    """
    props = {
        "non_volatile_solute_list": ["TDS"],
        "volatile_solute_list": ["TCA"],
        "mw_data": {
            "TCA": 0.1334,
            "TDS": 31.4038218e-3,
        },
        "dynamic_viscosity_data": {"Liq": 0.00115, "Vap": 1.75e-5},
        "henry_constant_data": {"TCA": 0.725},  # salinity adjusted
        "standard_enthalpy_change_data": {"TCA": 28.7e3},
        "temperature_boiling_data": {"TCA": 347},
        "molar_volume_data": {"TCA": 9.81e-5},
        "critical_molar_volume_data": {"TCA": 2.94e-4},
        "density_calculation": DensityCalculation.calculated,
    }
    m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = AirWaterEq(**props)

    return m


@pytest.mark.unit
def test_parameter_block1(m1):
    m = m1

    assert m.fs.properties.config.temp_adjust_henry
    assert isinstance(m.fs.properties.enth_change_dissolution_comp, Var)

    assert isinstance(
        m.fs.properties.config.liq_diffus_calculation, LiqDiffusivityCalculation
    )
    assert (
        m.fs.properties.config.liq_diffus_calculation
        == LiqDiffusivityCalculation.HaydukLaudie
    )
    assert isinstance(
        m.fs.properties.config.vap_diffus_calculation, VapDiffusivityCalculation
    )
    assert (
        m.fs.properties.config.vap_diffus_calculation
        == VapDiffusivityCalculation.WilkeLee
    )
    assert isinstance(
        m.fs.properties.config.molar_volume_calculation, MolarVolumeCalculation
    )
    assert (
        m.fs.properties.config.molar_volume_calculation
        == MolarVolumeCalculation.TynCalus
    )
    assert isinstance(m.fs.properties.config.density_calculation, DensityCalculation)
    assert m.fs.properties.config.density_calculation == DensityCalculation.calculated
    assert (
        m.fs.properties.config.saturation_vapor_pressure_calculation
        == SaturationVaporPressureCalculation.ArdenBuck
    )
    assert (
        m.fs.properties.config.vapor_pressure_calculation
        == VaporPressureCalculation.FromRelativeHumidity
    )
    assert (
        m.fs.properties.config.relative_humidity_calculation
        == RelativeHumidityCalculation.none
    )
    assert (
        m.fs.properties.config.latent_heat_of_vaporization_calculation
        == LatentHeatVaporizationCalculation.Sharqawy
    )
    assert (
        m.fs.properties.config.specific_heat_of_water_calculation
        == SpecificHeatWaterCalculation.Sharqawy
    )


@pytest.mark.component
def test_properties1(m1):
    m = m1

    m.fs.stream = m.fs.properties.build_state_block([0], defined_state=True)
    stream = m.fs.stream[0]
    stream.flow_mass_phase_comp["Liq", "H2O"].fix(157.8657)
    stream.flow_mass_phase_comp["Liq", "TDS"].fix(0.157)
    stream.flow_mass_phase_comp["Liq", "TCA"].fix(2.61e-5)
    stream.flow_mass_phase_comp["Vap", "Air"].fix(1.34932)
    stream.flow_mass_phase_comp["Vap", "TCA"].fix(0)
    stream.flow_mass_phase_comp["Vap", "H2O"].fix(0.10)
    stream.temperature["Liq"].fix(283)
    stream.temperature["Vap"].fix(283)
    stream.pressure.fix(101325)

    stream.conc_mass_phase_comp[...]

    stream.flow_mol_phase_comp[...]
    stream.conc_mol_phase_comp[...]
    stream.mole_frac_phase_comp[...]

    stream.diffus_phase_comp[...]
    stream.flow_vol_phase[...]
    stream.flow_mass_phase[...]
    stream.henry_comp[...]

    stream.pressure_vap_sat[...]
    stream.pressure_vap[...]
    stream.relative_humidity[...]
    stream.dh_vap_mass_solvent[...]
    stream.cp_mass_solvent[...]

    stream_vars = [
        "pressure",
        "temperature",
        "flow_mass_phase_comp",
        "conc_mass_phase_comp",
        "mass_frac_phase_comp",
        "flow_mol_phase_comp",
        "conc_mol_phase_comp",
        "mole_frac_phase_comp",
        "diffus_phase_comp",
        "molar_volume_comp",
        "flow_vol_phase",
        "flow_mass_phase",
        "henry_comp",
        "pressure_vap_sat",
        "pressure_vap",
        "dh_vap_mass_solvent",
        "cp_mass_solvent",
        "dens_mass_phase",
        "dens_mass_solvent",
    ]

    for prop in stream_vars:
        assert hasattr(stream, prop)
        assert isinstance(getattr(stream, prop), Var)

    calculate_scaling_factors(m)

    assert_units_consistent(m)

    check_dof(m, fail_flag=True)
    m.fs.stream.initialize()

    results = solver.solve(m)
    assert_optimal_termination(results)

    stream_results = {
        "pressure": 101325.0,
        "temperature": {"Liq": 283.0, "Vap": 283.0},
        "flow_mass_phase_comp": {
            ("Liq", "TCA"): 2.61e-05,
            ("Liq", "TDS"): 0.157,
            ("Liq", "H2O"): 157.86,
            ("Vap", "TCA"): 0.0,
            ("Vap", "H2O"): 0.1,
            ("Vap", "Air"): 1.3493,
        },
        "conc_mass_phase_comp": {
            ("Liq", "TCA"): 0.0001652,
            ("Liq", "TDS"): 0.99382,
            ("Liq", "H2O"): 999.3,
            ("Vap", "TCA"): 0.0,
            ("Vap", "H2O"): 0.08608,
            ("Vap", "Air"): 1.1615,
        },
        "dens_mass_phase": {"Liq": 1000.3024, "Vap": 1.2476},
        "mass_frac_phase_comp": {
            ("Liq", "TCA"): 1.651e-07,
            ("Liq", "TDS"): 0.0009935,
            ("Liq", "H2O"): 0.99900630,
            ("Vap", "TCA"): 0.0,
            ("Vap", "H2O"): 0.068997,
            ("Vap", "Air"): 0.931,
        },
        "dens_mass_solvent": {"H2O": 999.5236, "Air": 1.2476},
        "flow_mol_phase_comp": {
            ("Liq", "TCA"): 0.0001956,
            ("Liq", "TDS"): 4.999,
            ("Liq", "H2O"): 8770.3,
            ("Vap", "TCA"): 0.0,
            ("Vap", "H2O"): 5.555,
            ("Vap", "Air"): 46.5282,
        },
        "conc_mol_phase_comp": {
            ("Liq", "TCA"): 0.001238,
            ("Liq", "TDS"): 31.64,
            ("Liq", "H2O"): 55517.1,
            ("Vap", "TCA"): 0.0,
            ("Vap", "H2O"): 4.782,
            ("Vap", "Air"): 40.054,
        },
        "mole_frac_phase_comp": {
            ("Liq", "TCA"): 2.2295e-08,
            ("Liq", "TDS"): 0.0005697,
            ("Liq", "H2O"): 0.99943,
            ("Vap", "TCA"): 0.0,
            ("Vap", "H2O"): 0.10666,
            ("Vap", "Air"): 0.89333,
        },
        "diffus_phase_comp": {("Vap", "TCA"): 7.673219e-06, ("Liq", "TCA"): 7.09e-10},
        "collision_molecular_separation": {"TCA": 0.46830},
        "collision_molecular_separation_comp": {"TCA": 0.56550},
        "molar_volume_comp": {"TCA": 0.0001100},
        "collision_function_comp": {"TCA": 0.59093},
        "collision_function_zeta_comp": {"TCA": -0.22846},
        "collision_function_ee_comp": {"TCA": 0.192517},
        "energy_molecular_attraction": {"TCA": 2.5081e-14},
        "energy_molecular_attraction_air": 1.0851e-14,
        "energy_molecular_attraction_comp": {"TCA": 5.7969e-14},
        "flow_vol_phase": {"Liq": 0.15797494, "Vap": 1.16161},
        "flow_mass_phase": {"Liq": 158.0, "Vap": 1.44932},
        "henry_comp": {"TCA": 0.392},
        "pressure_vap_sat": {"H2O": 1215.5},
        "arden_buck_exponential_term": 0.687,
        "pressure_vap": {"H2O": 303.8},
        "dh_vap_mass_solvent": 2477.6,
        "cp_mass_solvent": {"Liq": 4197.1, "Vap": 1861.6},
    }

    for prop, d in stream_results.items():
        sv = getattr(stream, prop)
        if isinstance(d, dict):
            for i, val in d.items():
                assert value(sv[i]) == pytest.approx(val, rel=1e-3)
        else:
            assert value(sv) == pytest.approx(d, rel=1e-3)


@pytest.fixture(scope="module")
def m2():
    """
    Test NMSU case study 2 results
    """
    props = {
        "volatile_solute_list": ["DCP"],
        "mw_data": {"DCP": 0.11298},
        "dynamic_viscosity_data": {"Liq": 0.001307, "Vap": 1.79e-5},
        "henry_constant_data": {"DCP": 0.146},  # salinity adjusted
        "standard_enthalpy_change_data": {"DCP": 31.1e3},
        "temperature_boiling_data": {"DCP": 369.1},
        "molar_volume_data": {"DCP": 0.00011007},
        "critical_molar_volume_data": {"DCP": 2.26e-4},
        "density_data": {"Liq": 999.7, "Vap": 1.247},
        "relative_humidity_data": 0.5,
        "saturation_vapor_pressure_calculation": SaturationVaporPressureCalculation.Antoine,
    }
    m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = AirWaterEq(**props)

    return m


@pytest.mark.unit
def test_parameter_block2(m2):
    m = m2

    assert isinstance(m.fs.properties.component_list, Set)
    assert len(m.fs.properties.component_list) == 3
    for j in m.fs.properties.component_list:
        assert j in ["H2O", "Air", "DCP"]

    assert isinstance(m.fs.properties.phase_list, Set)
    assert len(m.fs.properties.phase_list) == 2
    assert m.fs.properties.phase_list == ["Liq", "Vap"]

    ### Liquid components, Sets
    assert isinstance(m.fs.properties.liq_comps, Set)
    assert len(m.fs.properties.liq_comps) == 2
    for j in m.fs.properties.liq_comps:
        assert j in ["DCP", "H2O"]

    assert isinstance(m.fs.properties.liq_comp_set, Set)
    assert len(m.fs.properties.liq_comp_set) == 2
    assert len(m.fs.properties.liq_comp_set) == len(m.fs.properties.liq_comps)
    for p, j in m.fs.properties.liq_comp_set:
        assert j in ["DCP", "H2O"]
        assert p == "Liq"

    assert isinstance(m.fs.properties.non_volatile_comps, Set)
    assert len(m.fs.properties.non_volatile_comps) == 0
    assert len(m.fs.properties.non_volatile_comps) == len(
        m.fs.properties.config.non_volatile_solute_list
    )

    assert isinstance(m.fs.properties.volatile_comp_set, Set)
    assert len(m.fs.properties.volatile_comp_set) == 2
    assert (
        len(m.fs.properties.volatile_comp_set)
        == len(m.fs.properties.volatile_comps) * 2
    )

    ### Vapor components, Sets
    assert isinstance(m.fs.properties.vap_comps, Set)
    assert len(m.fs.properties.vap_comps) == 3
    for j in m.fs.properties.vap_comps:
        assert j in ["DCP", "H2O", "Air"]

    assert isinstance(m.fs.properties.vap_comp_set, Set)
    assert len(m.fs.properties.vap_comp_set) == 3
    assert len(m.fs.properties.vap_comp_set) == len(m.fs.properties.vap_comps)
    for p, j in m.fs.properties.vap_comp_set:
        assert j in ["DCP", "H2O", "Air"]
        assert p in ["Liq", "Vap"]

    assert isinstance(m.fs.properties.volatile_comps, Set)
    assert len(m.fs.properties.volatile_comps) == 1
    assert len(m.fs.properties.volatile_comps) == len(
        m.fs.properties.config.volatile_solute_list
    )
    for j in m.fs.properties.volatile_comps:
        assert j == "DCP"
        assert j == m.fs.properties.config.volatile_solute_list[0]

    assert isinstance(m.fs.properties.volatile_comp_set, Set)
    assert len(m.fs.properties.volatile_comp_set) == 2
    assert (
        len(m.fs.properties.volatile_comp_set)
        == len(m.fs.properties.volatile_comps) + 1
    )
    for p, j in m.fs.properties.volatile_comp_set:
        assert j == "DCP"
        assert p in ["Liq", "Vap"]

    assert isinstance(m.fs.properties.component_set, Set)
    assert len(m.fs.properties.component_set) == len(
        m.fs.properties.config.non_volatile_solute_list
    ) + len(m.fs.properties.config.volatile_solute_list) + len(["H2O", "Air"])
    assert len(m.fs.properties.component_set) == len(m.fs.properties.solute_set)

    assert m.fs.properties._state_block_class is AirWaterEqStateBlock
    assert isinstance(m.fs.properties.mw_comp, Param)
    assert m.fs.properties.mw_comp["H2O"].value == 18e-3
    assert m.fs.properties.mw_comp["Air"].value == 29e-3
    assert m.fs.properties.mw_comp["DCP"].value == 112.98e-3

    assert isinstance(m.fs.properties.visc_d_phase, Var)
    assert isinstance(m.fs.properties.molar_volume_comp, Var)
    assert isinstance(m.fs.properties.molar_volume_comp_crit, Var)
    assert isinstance(m.fs.properties.henry_comp, Var)
    assert isinstance(m.fs.properties.temperature_boiling_comp, Var)

    assert m.fs.properties.config.temp_adjust_henry
    assert isinstance(m.fs.properties.enth_change_dissolution_comp, Var)

    assert isinstance(
        m.fs.properties.config.liq_diffus_calculation, LiqDiffusivityCalculation
    )
    assert (
        m.fs.properties.config.liq_diffus_calculation
        == LiqDiffusivityCalculation.HaydukLaudie
    )
    assert isinstance(
        m.fs.properties.config.vap_diffus_calculation, VapDiffusivityCalculation
    )
    assert (
        m.fs.properties.config.vap_diffus_calculation
        == VapDiffusivityCalculation.WilkeLee
    )
    assert isinstance(
        m.fs.properties.config.molar_volume_calculation, MolarVolumeCalculation
    )
    assert (
        m.fs.properties.config.molar_volume_calculation
        == MolarVolumeCalculation.TynCalus
    )
    assert isinstance(m.fs.properties.config.density_calculation, DensityCalculation)
    assert m.fs.properties.config.density_calculation == DensityCalculation.constant
    assert (
        m.fs.properties.config.saturation_vapor_pressure_calculation
        == SaturationVaporPressureCalculation.Antoine
    )
    assert (
        m.fs.properties.config.vapor_pressure_calculation
        == VaporPressureCalculation.FromRelativeHumidity
    )
    assert (
        m.fs.properties.config.relative_humidity_calculation
        == RelativeHumidityCalculation.none
    )
    assert (
        m.fs.properties.config.latent_heat_of_vaporization_calculation
        == LatentHeatVaporizationCalculation.Sharqawy
    )
    assert (
        m.fs.properties.config.specific_heat_of_water_calculation
        == SpecificHeatWaterCalculation.Sharqawy
    )


@pytest.mark.component
def test_properties2(m2):
    m = m2
    m.fs.stream = m.fs.properties.build_state_block([0], defined_state=True)
    stream = m.fs.stream[0]

    stream.flow_mass_phase_comp["Liq", "H2O"].fix(99.97)
    stream.flow_mass_phase_comp["Liq", "DCP"].fix(1e-4)
    stream.flow_mass_phase_comp["Vap", "H2O"].fix(1e-3)
    stream.flow_mass_phase_comp["Vap", "Air"].fix(7.482)
    stream.flow_mass_phase_comp["Vap", "DCP"].fix(0)
    stream.temperature["Liq"].fix(283)
    stream.temperature["Vap"].fix(283)
    stream.pressure.fix(101325)

    stream.conc_mass_phase_comp[...]

    stream.flow_mol_phase_comp[...]
    stream.conc_mol_phase_comp[...]
    stream.mole_frac_phase_comp[...]

    stream.diffus_phase_comp[...]
    stream.flow_vol_phase[...]
    stream.flow_mass_phase[...]
    stream.henry_comp[...]

    stream.pressure_vap_sat[...]
    stream.pressure_vap[...]
    stream.relative_humidity[...]
    stream.dh_vap_mass_solvent[...]
    stream.cp_mass_solvent[...]

    stream_vars = [
        "pressure",
        "temperature",
        "flow_mass_phase_comp",
        "conc_mass_phase_comp",
        "mass_frac_phase_comp",
        "flow_mol_phase_comp",
        "conc_mol_phase_comp",
        "mole_frac_phase_comp",
        "diffus_phase_comp",
        "molar_volume_comp",
        "flow_vol_phase",
        "flow_mass_phase",
        "henry_comp",
        "pressure_vap_sat",
        "pressure_vap",
        "dh_vap_mass_solvent",
        "cp_mass_solvent",
    ]

    for prop in stream_vars:
        assert hasattr(stream, prop)
        assert isinstance(getattr(stream, prop), Var)

    calculate_scaling_factors(m)

    assert_units_consistent(m)

    check_dof(m, fail_flag=True)

    m.fs.stream.initialize()

    results = solver.solve(m)
    assert_optimal_termination(results)
    stream_results = {
        "pressure": 101325.0,
        "temperature": {"Liq": 283.0, "Vap": 283.0},
        "flow_mass_phase_comp": {
            ("Liq", "DCP"): 0.0001,
            ("Liq", "H2O"): 99.97,
            ("Vap", "DCP"): 0.0,
            ("Vap", "H2O"): 0.001,
            ("Vap", "Air"): 7.482,
        },
        "conc_mass_phase_comp": {
            ("Liq", "DCP"): 0.000999,
            ("Liq", "H2O"): 999.699,
            ("Vap", "DCP"): 0.0,
            ("Vap", "H2O"): 0.0001666,
            ("Vap", "Air"): 1.24683335,
        },
        "mass_frac_phase_comp": {
            ("Liq", "DCP"): 1.000299e-06,
            ("Liq", "H2O"): 0.99999899,
            ("Vap", "DCP"): 0.0,
            ("Vap", "H2O"): 0.00013363,
            ("Vap", "Air"): 0.99986636,
        },
        "flow_mol_phase_comp": {
            ("Liq", "DCP"): 0.0008851,
            ("Liq", "H2O"): 5553.8,
            ("Vap", "DCP"): 2e-12,
            ("Vap", "H2O"): 0.05555,
            ("Vap", "Air"): 258.0,
        },
        "conc_mol_phase_comp": {
            ("Liq", "DCP"): 0.0088511,
            ("Liq", "H2O"): 55538.8,
            ("Vap", "DCP"): 0.0,
            ("Vap", "H2O"): 0.0092580,
            ("Vap", "Air"): 42.9942,
        },
        "mole_frac_phase_comp": {
            ("Liq", "DCP"): 1.59368e-07,
            ("Liq", "H2O"): 0.9999,
            ("Vap", "DCP"): 0.0,
            ("Vap", "H2O"): 0.0002152,
            ("Vap", "Air"): 0.99978471,
        },
        "diffus_phase_comp": {("Vap", "DCP"): 8.579253e-06, ("Liq", "DCP"): 7.21e-10},
        "molar_volume_comp": {"DCP": 8.3550e-05},
        "flow_vol_phase": {"Liq": 0.1000, "Vap": 6.000},
        "flow_mass_phase": {"Liq": 99.9701, "Vap": 7.483},
        "henry_comp": {"DCP": 0.0750617},
        "pressure_vap_sat": {"H2O": 1208.812},
        "pressure_vap": {"H2O": 604.40},
        "dh_vap_mass_solvent": 2477.6833,
        "cp_mass_solvent": {"Liq": 4197.1362, "Vap": 1861.6326},
        "collision_molecular_separation": {"DCP": 0.44348019},
        "collision_molecular_separation_comp": {"DCP": 0.51586},
        "collision_function_comp": {"DCP": 0.59831308},
        "collision_function_zeta_comp": {"DCP": -0.2230715},
        "collision_function_ee_comp": {"DCP": 0.17911045},
        "energy_molecular_attraction": {"DCP": 2.586778e-14},
        "energy_molecular_attraction_comp": {"DCP": 6.1661e-14},
    }

    for prop, d in stream_results.items():
        sv = getattr(stream, prop)
        if isinstance(d, dict):
            for i, val in d.items():
                assert value(sv[i]) == pytest.approx(val, rel=1e-3)
        else:
            assert value(sv) == pytest.approx(d, rel=1e-3)


@pytest.fixture(scope="module")
def m3():
    """
    Test when molar volume, diffusivity, and Henry constant are all direct user input via CONFIG
    """
    props = {
        "volatile_solute_list": ["DCP"],
        "mw_data": {"DCP": 0.11298},
        "diffusivity_data": {("Liq", "DCP"): 7.21045e-10, ("Vap", "DCP"): 8.57928e-6},
        "dynamic_viscosity_data": {"Liq": 0.001307, "Vap": 1.79e-5},
        "henry_constant_data": {"DCP": 0.0525},  # salinity adjusted
        "standard_enthalpy_change_data": {"DCP": 31.1e3},
        "temperature_boiling_data": {"DCP": 369.1},
        "molar_volume_data": {"DCP": 8.355e-5},
        "density_data": {"Liq": 999.7, "Vap": 1.247},
        "temp_adjust_henry": False,
        "pressure_vap_sat_data": 1215.5768,
        "pressure_vap_data": 607.7884,
        "latent_heat_of_vaporization_data": 2477.6833,
        "specific_heat_of_water_data": {"Liq": 4197.1362, "Vap": 1861.6326},
        "relative_humidity_data": 0.5,
        "molar_volume_calculation": "none",
        "liq_diffus_calculation": "none",
        "vap_diffus_calculation": "none",
        "saturation_vapor_pressure_calculation": "none",
        "vapor_pressure_calculation": "none",
        "latent_heat_of_vaporization_calculation": "none",
        "specific_heat_of_water_calculation": "none",
        "relative_humidity_calculation": "none",
    }
    m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = AirWaterEq(**props)

    return m


@pytest.mark.unit
def test_parameter_block3(m3):
    m = m3

    assert isinstance(m.fs.properties.component_list, Set)
    assert len(m.fs.properties.component_list) == 3
    for j in m.fs.properties.component_list:
        assert j in ["H2O", "Air", "DCP"]

    assert isinstance(m.fs.properties.phase_list, Set)
    assert len(m.fs.properties.phase_list) == 2
    assert m.fs.properties.phase_list == ["Liq", "Vap"]

    ### Liquid components, Sets
    assert isinstance(m.fs.properties.liq_comps, Set)
    assert len(m.fs.properties.liq_comps) == 2
    for j in m.fs.properties.liq_comps:
        assert j in ["DCP", "H2O"]

    assert isinstance(m.fs.properties.liq_comp_set, Set)
    assert len(m.fs.properties.liq_comp_set) == 2
    assert len(m.fs.properties.liq_comp_set) == len(m.fs.properties.liq_comps)
    for p, j in m.fs.properties.liq_comp_set:
        assert j in ["DCP", "H2O"]
        assert p == "Liq"

    assert isinstance(m.fs.properties.non_volatile_comps, Set)
    assert len(m.fs.properties.non_volatile_comps) == 0
    assert len(m.fs.properties.non_volatile_comps) == len(
        m.fs.properties.config.non_volatile_solute_list
    )

    assert isinstance(m.fs.properties.volatile_comp_set, Set)
    assert len(m.fs.properties.volatile_comp_set) == 2
    assert (
        len(m.fs.properties.volatile_comp_set)
        == len(m.fs.properties.volatile_comps) * 2
    )

    ### Vapor components, Sets
    assert isinstance(m.fs.properties.vap_comps, Set)
    assert len(m.fs.properties.vap_comps) == 3
    for j in m.fs.properties.vap_comps:
        assert j in ["DCP", "H2O", "Air"]

    assert isinstance(m.fs.properties.vap_comp_set, Set)
    assert len(m.fs.properties.vap_comp_set) == 3
    assert len(m.fs.properties.vap_comp_set) == len(m.fs.properties.vap_comps)
    for p, j in m.fs.properties.vap_comp_set:
        assert j in ["DCP", "H2O", "Air"]
        assert p in ["Liq", "Vap"]

    assert isinstance(m.fs.properties.volatile_comps, Set)
    assert len(m.fs.properties.volatile_comps) == 1
    assert len(m.fs.properties.volatile_comps) == len(
        m.fs.properties.config.volatile_solute_list
    )
    for j in m.fs.properties.volatile_comps:
        assert j == "DCP"
        assert j == m.fs.properties.config.volatile_solute_list[0]

    assert isinstance(m.fs.properties.volatile_comp_set, Set)
    assert len(m.fs.properties.volatile_comp_set) == 2
    assert (
        len(m.fs.properties.volatile_comp_set)
        == len(m.fs.properties.volatile_comps) + 1
    )
    for p, j in m.fs.properties.volatile_comp_set:
        assert j == "DCP"
        assert p in ["Liq", "Vap"]

    assert isinstance(m.fs.properties.component_set, Set)
    assert len(m.fs.properties.component_set) == len(
        m.fs.properties.config.non_volatile_solute_list
    ) + len(m.fs.properties.config.volatile_solute_list) + len(["H2O", "Air"])
    assert len(m.fs.properties.component_set) == len(m.fs.properties.solute_set)

    assert m.fs.properties._state_block_class is AirWaterEqStateBlock
    assert isinstance(m.fs.properties.mw_comp, Param)
    assert m.fs.properties.mw_comp["H2O"].value == 18e-3
    assert m.fs.properties.mw_comp["Air"].value == 29e-3
    assert m.fs.properties.mw_comp["DCP"].value == 112.98e-3

    assert isinstance(m.fs.properties.visc_d_phase, Var)
    assert isinstance(m.fs.properties.molar_volume_comp, Var)
    assert isinstance(m.fs.properties.molar_volume_comp_crit, Var)
    assert isinstance(m.fs.properties.henry_comp, Var)
    assert isinstance(m.fs.properties.temperature_boiling_comp, Var)

    assert not m.fs.properties.config.temp_adjust_henry

    assert isinstance(
        m.fs.properties.config.liq_diffus_calculation, LiqDiffusivityCalculation
    )
    assert (
        m.fs.properties.config.liq_diffus_calculation == LiqDiffusivityCalculation.none
    )
    assert isinstance(
        m.fs.properties.config.vap_diffus_calculation, VapDiffusivityCalculation
    )
    assert (
        m.fs.properties.config.vap_diffus_calculation == VapDiffusivityCalculation.none
    )
    assert isinstance(
        m.fs.properties.config.molar_volume_calculation, MolarVolumeCalculation
    )
    assert (
        m.fs.properties.config.molar_volume_calculation == MolarVolumeCalculation.none
    )
    assert isinstance(m.fs.properties.config.density_calculation, DensityCalculation)
    assert m.fs.properties.config.density_calculation == DensityCalculation.constant
    assert (
        m.fs.properties.config.saturation_vapor_pressure_calculation
        == SaturationVaporPressureCalculation.none
    )
    assert (
        m.fs.properties.config.vapor_pressure_calculation
        == VaporPressureCalculation.none
    )
    assert (
        m.fs.properties.config.relative_humidity_calculation
        == RelativeHumidityCalculation.none
    )
    assert (
        m.fs.properties.config.latent_heat_of_vaporization_calculation
        == LatentHeatVaporizationCalculation.none
    )
    assert (
        m.fs.properties.config.specific_heat_of_water_calculation
        == SpecificHeatWaterCalculation.none
    )


@pytest.mark.component
def test_properties3(m3):
    m = m3

    m.fs.stream = m.fs.properties.build_state_block([0], defined_state=True)
    stream = m.fs.stream[0]
    stream.flow_mass_phase_comp["Liq", "H2O"].fix(99.97)
    stream.flow_mass_phase_comp["Liq", "DCP"].fix(1e-4)
    stream.flow_mass_phase_comp["Vap", "H2O"].fix(1e-3)
    stream.flow_mass_phase_comp["Vap", "Air"].fix(7.482)
    stream.flow_mass_phase_comp["Vap", "DCP"].fix(0)
    stream.temperature["Liq"].fix(283)
    stream.temperature["Vap"].fix(283)
    stream.pressure.fix(101325)

    stream.conc_mass_phase_comp[...]

    stream.flow_mol_phase_comp[...]
    stream.conc_mol_phase_comp[...]
    stream.mole_frac_phase_comp[...]

    stream.diffus_phase_comp[...]
    stream.flow_vol_phase[...]
    stream.flow_mass_phase[...]
    stream.henry_comp[...]

    stream.pressure_vap_sat[...]
    stream.pressure_vap[...]
    stream.relative_humidity[...]
    stream.dh_vap_mass_solvent[...]
    stream.cp_mass_solvent[...]

    stream_vars = [
        "pressure",
        "temperature",
        "flow_mass_phase_comp",
        "conc_mass_phase_comp",
        "mass_frac_phase_comp",
        "flow_mol_phase_comp",
        "conc_mol_phase_comp",
        "mole_frac_phase_comp",
        "diffus_phase_comp",
        "molar_volume_comp",
        "flow_vol_phase",
        "flow_mass_phase",
        "henry_comp",
        "pressure_vap_sat",
        "pressure_vap",
        "dh_vap_mass_solvent",
        "cp_mass_solvent",
    ]

    for prop in stream_vars:
        assert hasattr(stream, prop)
        assert isinstance(getattr(stream, prop), Var)

    calculate_scaling_factors(m)

    assert_units_consistent(m)

    check_dof(m, fail_flag=True)
    m.fs.stream.initialize()

    results = solver.solve(m)
    assert_optimal_termination(results)

    stream_results = {
        "pressure": 101325.0,
        "temperature": {"Liq": 283.0, "Vap": 283.0},
        "flow_mass_phase_comp": {
            ("Liq", "DCP"): 0.0001,
            ("Liq", "H2O"): 99.97,
            ("Vap", "DCP"): 0.0,
            ("Vap", "H2O"): 0.001,
            ("Vap", "Air"): 7.482,
        },
        "conc_mass_phase_comp": {
            ("Liq", "DCP"): 0.000999,
            ("Liq", "H2O"): 999.7,
            ("Vap", "DCP"): 0.0,
            ("Vap", "H2O"): 0.0001666,
            ("Vap", "Air"): 1.2468,
        },
        "mass_frac_phase_comp": {
            ("Liq", "DCP"): 1.00e-06,
            ("Liq", "H2O"): 0.99999,
            ("Vap", "DCP"): 0.0,
            ("Vap", "H2O"): 0.0001336,
            ("Vap", "Air"): 0.99986,
        },
        "flow_mol_phase_comp": {
            ("Liq", "DCP"): 0.0008851,
            ("Liq", "H2O"): 5553.8,
            ("Vap", "DCP"): 2e-12,
            ("Vap", "H2O"): 0.055555,
            ("Vap", "Air"): 258.0,
        },
        "conc_mol_phase_comp": {
            ("Liq", "DCP"): 0.008851,
            ("Liq", "H2O"): 55538.8,
            ("Vap", "DCP"): 0.0,
            ("Vap", "H2O"): 0.009258,
            ("Vap", "Air"): 42.9942,
        },
        "mole_frac_phase_comp": {
            ("Liq", "DCP"): 1.593e-07,
            ("Liq", "H2O"): 0.9999,
            ("Vap", "DCP"): 0.0,
            ("Vap", "H2O"): 0.0002152,
            ("Vap", "Air"): 0.999784,
        },
        "diffus_phase_comp": {("Vap", "DCP"): 8.579e-06, ("Liq", "DCP"): 7.21e-10},
        "flow_vol_phase": {"Liq": 0.10, "Vap": 6.0},
        "flow_mass_phase": {"Liq": 99.97, "Vap": 7.48},
    }
    for prop, d in stream_results.items():
        sv = getattr(stream, prop)
        if isinstance(d, dict):
            for i, val in d.items():
                assert value(sv[i]) == pytest.approx(val, rel=1e-3)
        else:
            assert value(sv) == pytest.approx(d, rel=1e-3)


@pytest.fixture(scope="module")
def m4():
    """
    Test air-water ONLY system
    """
    props = {
        "relative_humidity_data": 0.67,
        "vapor_pressure_calculation": VaporPressureCalculation.FromRelativeHumidity,
        "saturation_vapor_pressure_calculation": SaturationVaporPressureCalculation.Huang,
    }

    m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = AirWaterEq(**props)

    return m


@pytest.mark.unit
def test_parameter_block4(m4):
    m = m4

    assert len(m.fs.properties.config.volatile_solute_list) == 0
    assert len(m.fs.properties.config.non_volatile_solute_list) == 0
    assert isinstance(m.fs.properties.component_list, Set)
    assert len(m.fs.properties.component_list) == 2
    for j in m.fs.properties.component_list:
        assert j in ["H2O", "Air"]

    assert isinstance(m.fs.properties.phase_list, Set)
    assert len(m.fs.properties.phase_list) == 2
    assert m.fs.properties.phase_list == ["Liq", "Vap"]

    ### Liquid components, Sets
    assert isinstance(m.fs.properties.liq_comps, Set)
    assert len(m.fs.properties.liq_comps) == 1
    for j in m.fs.properties.liq_comps:
        assert j == "H2O"

    assert isinstance(m.fs.properties.liq_comp_set, Set)
    assert len(m.fs.properties.liq_comp_set) == 1
    assert len(m.fs.properties.liq_comp_set) == len(m.fs.properties.liq_comps)
    for p, j in m.fs.properties.liq_comp_set:
        assert j == "H2O"
        assert p == "Liq"

    assert isinstance(m.fs.properties.non_volatile_comps, Set)
    assert len(m.fs.properties.non_volatile_comps) == 0
    assert len(m.fs.properties.non_volatile_comps) == len(
        m.fs.properties.config.non_volatile_solute_list
    )

    assert isinstance(m.fs.properties.volatile_comp_set, Set)
    assert len(m.fs.properties.volatile_comp_set) == 0
    assert (
        len(m.fs.properties.volatile_comp_set)
        == len(m.fs.properties.volatile_comps) * 2
    )

    ### Vapor components, Sets
    assert isinstance(m.fs.properties.vap_comps, Set)
    assert len(m.fs.properties.vap_comps) == 2
    for j in m.fs.properties.vap_comps:
        assert j in ["H2O", "Air"]

    assert isinstance(m.fs.properties.vap_comp_set, Set)
    assert len(m.fs.properties.vap_comp_set) == 2
    assert len(m.fs.properties.vap_comp_set) == len(m.fs.properties.vap_comps)
    for p, j in m.fs.properties.vap_comp_set:
        assert j in ["H2O", "Air"]
        assert p in ["Liq", "Vap"]

    assert isinstance(m.fs.properties.volatile_comps, Set)
    assert len(m.fs.properties.volatile_comps) == 0
    assert len(m.fs.properties.volatile_comps) == len(
        m.fs.properties.config.volatile_solute_list
    )

    assert isinstance(m.fs.properties.volatile_comp_set, Set)
    assert len(m.fs.properties.volatile_comp_set) == 0
    assert (
        len(m.fs.properties.volatile_comp_set)
        == len(m.fs.properties.volatile_comps) * 2
    )

    assert isinstance(m.fs.properties.component_set, Set)
    assert len(m.fs.properties.component_set) == len(
        m.fs.properties.config.non_volatile_solute_list
    ) + len(m.fs.properties.config.volatile_solute_list) + len(["H2O", "Air"])
    assert len(m.fs.properties.component_set) == len(m.fs.properties.solute_set)

    assert m.fs.properties._state_block_class is AirWaterEqStateBlock
    assert isinstance(m.fs.properties.mw_comp, Param)
    assert m.fs.properties.mw_comp["H2O"].value == 18e-3
    assert m.fs.properties.mw_comp["Air"].value == 29e-3

    assert isinstance(m.fs.properties.visc_d_phase, Var)
    assert isinstance(m.fs.properties.molar_volume_comp, Var)
    assert isinstance(m.fs.properties.molar_volume_comp_crit, Var)
    assert isinstance(m.fs.properties.henry_comp, Var)
    assert isinstance(m.fs.properties.temperature_boiling_comp, Var)

    assert m.fs.properties.config.temp_adjust_henry
    assert isinstance(m.fs.properties.enth_change_dissolution_comp, Var)

    assert isinstance(
        m.fs.properties.config.liq_diffus_calculation, LiqDiffusivityCalculation
    )
    assert (
        m.fs.properties.config.liq_diffus_calculation
        == LiqDiffusivityCalculation.HaydukLaudie
    )
    assert isinstance(
        m.fs.properties.config.vap_diffus_calculation, VapDiffusivityCalculation
    )
    assert (
        m.fs.properties.config.vap_diffus_calculation
        == VapDiffusivityCalculation.WilkeLee
    )
    assert isinstance(
        m.fs.properties.config.molar_volume_calculation, MolarVolumeCalculation
    )
    assert (
        m.fs.properties.config.molar_volume_calculation
        == MolarVolumeCalculation.TynCalus
    )
    assert isinstance(m.fs.properties.config.density_calculation, DensityCalculation)
    assert m.fs.properties.config.density_calculation == DensityCalculation.constant
    assert (
        m.fs.properties.config.saturation_vapor_pressure_calculation
        == SaturationVaporPressureCalculation.Huang
    )
    assert (
        m.fs.properties.config.vapor_pressure_calculation
        == VaporPressureCalculation.FromRelativeHumidity
    )
    assert (
        m.fs.properties.config.relative_humidity_calculation
        == RelativeHumidityCalculation.none
    )
    assert (
        m.fs.properties.config.latent_heat_of_vaporization_calculation
        == LatentHeatVaporizationCalculation.Sharqawy
    )
    assert (
        m.fs.properties.config.specific_heat_of_water_calculation
        == SpecificHeatWaterCalculation.Sharqawy
    )


@pytest.mark.component
def test_properties4(m4):
    m = m4

    m.fs.stream = m.fs.properties.build_state_block([0], defined_state=True)
    stream = m.fs.stream[0]
    stream.flow_mass_phase_comp["Liq", "H2O"].fix(99.97)
    stream.flow_mass_phase_comp["Vap", "H2O"].fix(1e-3)
    stream.flow_mass_phase_comp["Vap", "Air"].fix(0.2)
    stream.temperature["Liq"].fix(322)
    stream.temperature["Vap"].fix(283)
    stream.pressure.fix(101325)

    stream.conc_mass_phase_comp[...]

    stream.flow_mol_phase_comp[...]
    stream.conc_mol_phase_comp[...]
    stream.mole_frac_phase_comp[...]

    stream.flow_vol_phase[...]
    stream.flow_mass_phase[...]

    stream.pressure_vap_sat[...]
    stream.pressure_vap[...]
    stream.relative_humidity[...]
    stream.dh_vap_mass_solvent[...]
    stream.cp_mass_solvent[...]

    stream_vars = [
        "pressure",
        "temperature",
        "flow_mass_phase_comp",
        "conc_mass_phase_comp",
        "mass_frac_phase_comp",
        "flow_mol_phase_comp",
        "conc_mol_phase_comp",
        "mole_frac_phase_comp",
        "molar_volume_comp",
        "flow_vol_phase",
        "flow_mass_phase",
        "pressure_vap_sat",
        "pressure_vap",
        "relative_humidity",
        "dh_vap_mass_solvent",
        "cp_mass_solvent",
    ]

    for prop in stream_vars:
        assert hasattr(stream, prop)
        assert isinstance(getattr(stream, prop), Var)

    calculate_scaling_factors(m)

    assert_units_consistent(m)

    check_dof(m, fail_flag=True)
    m.fs.stream.initialize()

    results = solver.solve(m)
    assert_optimal_termination(results)

    stream_results = {
        "pressure": 101325.0,
        "temperature": {"Liq": 322.0, "Vap": 283.0},
        "flow_mass_phase_comp": {
            ("Liq", "H2O"): 99.97,
            ("Vap", "H2O"): 0.001,
            ("Vap", "Air"): 0.2,
        },
        "conc_mass_phase_comp": {
            ("Liq", "H2O"): 998.2,
            ("Vap", "H2O"): 0.00599,
            ("Vap", "Air"): 1.198,
        },
        "mass_frac_phase_comp": {
            ("Liq", "H2O"): 1.0,
            ("Vap", "H2O"): 0.004975,
            ("Vap", "Air"): 0.995024,
        },
        "flow_mol_phase_comp": {
            ("Liq", "H2O"): 5553.8,
            ("Vap", "H2O"): 0.055555,
            ("Vap", "Air"): 6.896,
        },
        "conc_mol_phase_comp": {
            ("Liq", "H2O"): 55455.5,
            ("Vap", "H2O"): 0.33278,
            ("Vap", "Air"): 41.31,
        },
        "mole_frac_phase_comp": {
            ("Liq", "H2O"): 0.99999,
            ("Vap", "H2O"): 0.00799,
            ("Vap", "Air"): 0.99200,
        },
        "flow_vol_phase": {"Liq": 0.100150, "Vap": 0.1669},
        "flow_mass_phase": {"Liq": 99.97, "Vap": 0.201},
        "pressure_vap_sat": {"H2O": 11663.9},
        "pressure_vap": {"H2O": 7814.8},
        "relative_humidity": {"H2O": 0.67},
        "dh_vap_mass_solvent": 2384.8,
        "cp_mass_solvent": {"Liq": 4180.8, "Vap": 1870.8},
    }
    for prop, d in stream_results.items():
        sv = getattr(stream, prop)
        if isinstance(d, dict):
            for i, val in d.items():
                assert value(sv[i]) == pytest.approx(val, rel=1e-3)
        else:
            assert value(sv) == pytest.approx(d, rel=1e-3)
