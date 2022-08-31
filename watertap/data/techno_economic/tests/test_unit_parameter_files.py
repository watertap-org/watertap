###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################
"""
Tests for loading water source definitions
"""
import pytest
import os

from pyomo.environ import units
from pyomo.util.check_units import assert_units_equivalent

from watertap.core import Database

dbpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

db = Database()

exclude_files = ["water_sources.yaml", "component_list.yaml", "default_case_study.yaml"]

tech_list = []
for f in os.listdir(dbpath):
    filename = os.fsdecode(f)
    if filename.endswith(".yaml") and filename not in exclude_files:
        tech_list.append(filename[:-5])


@pytest.mark.integration
@pytest.mark.parametrize("tech", tech_list)
def test_unit_parameter_files(tech):
    data = db._get_technology(tech)

    # Check that data has as default key
    assert "default" in data

    # Iterate overall entries in tech data and check for expected contents
    # TODO : Need to check up on this once everything is done
    pass_through = [
        "blending_reservoir",
        "buffer_tank",
        "chemical_addition",
        "co2_addition",
        "cooling_supply",
        "energy_recovery",
        "feed_water_tank",
        "injection_well_disposal",
        "intrusion_mitigation",
        "landfill",
        "municipal_drinking",
        "municipal_wwtp",
        "pump_electricity",
        "pump",
        "smp",
        "static_mixer",
        "storage_tank",
        "surface_discharge",
        "sw_onshore_intake",
        "tramp_oil_tank",
        "water_pumping_station",
        "well_field",
    ]

    siso_full_recovery = ["uv_aop", "uv", "fixed_bed", "decarbonator", "chlorination"]

    no_energy_electric_flow_vol_inlet = [
        "anaerobic_mbr_mec",
        "autothermal_hydrothermal_liquefaction",
        "backwash_solids_handling",
        "brine_concentrator",
        "CANDO_P",
        "chemical_addition",
        "coag_and_floc",
        "cofermentation",
        "constructed_wetlands",
        "deep_well_injection",
        "electrochemical_nutrient_removal",
        "electrodialysis_reversal",
        "energy_recovery",
        "evaporation_pond",
        "filter_press",
        "gac",
        "gas_sparged_membrane",
        "hydrothermal_gasification",
        "ion_exchange",
        "iron_and_manganese_removal",
        "mbr",
        "metab",
        "municipal_drinking",
        "ozonation",
        "ozone_aop",
        "photothermal_membrane",
        "pump_electricity",
        "storage_tank",
        "supercritical_salt_precipitation",
        "surface_discharge",
        "sw_onshore_intake",
        "vfa_recovery",
        "water_pumping_station",
        "well_field",
    ]

    expected = ["recovery_frac_mass_H2O", "default_removal_frac_mass_comp"]

    for k in data.values():

        for e in expected:
            if tech not in pass_through and tech not in siso_full_recovery:
                assert e in k.keys()
                assert "units" in k[e].keys()
                assert_units_equivalent(k[e]["units"], units.dimensionless)
                assert "value" in k[e].keys()
                assert k[e]["value"] >= 0
                assert k[e]["value"] <= 1
            elif tech in siso_full_recovery:
                if e == "default_removal_frac_mass_comp":
                    assert e in k.keys()
            else:
                assert e not in k.keys()

        if tech not in no_energy_electric_flow_vol_inlet:
            e = "energy_electric_flow_vol_inlet"
            assert e in k.keys()
            assert "units" in k[e].keys()
            assert_units_equivalent(k[e]["units"], units.dimensionless)
            assert "value" in k[e].keys()
            assert k[e]["value"] >= 0

        # Check for specific removal fractions
        if "removal_frac_mass_comp" in k.keys():
            for (j, c_data) in k["removal_frac_mass_comp"].items():
                assert "units" in c_data.keys()
                assert_units_equivalent(c_data["units"], units.dimensionless)
                assert "value" in c_data.keys()
                assert c_data["value"] >= 0
                assert c_data["value"] <= 1
                assert j in db.component_list.keys()
