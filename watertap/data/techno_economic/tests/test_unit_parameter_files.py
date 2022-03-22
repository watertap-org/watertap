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

exclude_files = [
    "water_sources.yaml", "component_list.yaml", "default_case_study.yaml"]

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
    pass_through = ["blending_reservoir",
                    "buffer_tank",
                    "chemical_addition",
                    "cooling_supply",
                    "co2_addition",
                    "energy_recovery",
                    "feed_water_tank",
                    "injection_well_disposal",
                    "intrusion_mitigation",
                    "landfill",
                    "municipal_drinking",
                    "municipal_wwtp",
                    "pump",
                    "storage_tank",
                    "static_mixer",
                    "smp",
                    "sw_onshore_intake",
                    "tramp_oil_tank",
                   ]

    siso_full_recovery = ["uv_aop", "uv", "fixed_bed", "decarbonator", "chlorination"]

    no_energy_electric_flow_vol_inlet = ["constructed_wetlands",
                                         "energy_recovery",
                                         "mbr_denitrification",
                                         "mbr_nitrification",
                                         "multi_stage_bubble_aeration",
                                         "tri_media_filtration",
                                         "cartridge_filtration_with_backflush",
                                         "landfill",
                                         "well_field",
                                         "ion_exchange",
                                         "ozone_aop",
                                         "fixed_bed",
                                         "holding_tank",
                                         "heap_leaching",
                                         "nuclear_cooling_tower",
                                         "lime_softening",
                                         "ozonation",
                                         "cooling_tower",
                                         "gac",
                                         "tri_media_filtration_with_backflush",
                                         "sedimentation",
                                         "backwash_solids_handling",
                                         "ph_decrease",
                                         "ph_increase",
                                         "co2_addition",
                                         "coag_and_floc",
                                         "crystallizer",
                                         "iron_and_manganese_removal",
                                         "fluidized_bed",
                                         "surface_discharge",
                                         "solution_distribution_and_recovery_plant",
                                         "chemical_addition",
                                         "cartridge_filtration",
                                         "injection_well",
                                         "sw_onshore_intake",
                                         "filter_press",
                                         "municipal_drinking",
                                         "gac_pressure_30_min",
                                         "packed_tower_aeration",
                                         "treated_storage",
                                         "gac_gravity_60_min",
                                         "evaporation_pond",
                                         "lime_addition",
                                         "brine_concentrator",
                                         "agglom_stacking",
                                         "storage_tank"]

    expected = ["recovery_frac_mass_H2O",
                "default_removal_frac_mass_solute"]

    for k in data.values():

        for e in expected:
            if tech not in pass_through and tech not in siso_full_recovery:
                assert e in k.keys()
                assert "units" in k[e].keys()
                assert_units_equivalent(
                    k[e]["units"], units.dimensionless)
                assert "value" in k[e].keys()
                assert k[e]["value"] >= 0
                assert k[e]["value"] <= 1
            elif tech in siso_full_recovery:
                if e == "default_removal_frac_mass_solute":
                    assert e in k.keys()
            else:
                assert e not in k.keys()

        if tech not in no_energy_electric_flow_vol_inlet:
            e = "energy_electric_flow_vol_inlet"
            assert e in k.keys()
            assert "units" in k[e].keys()
            assert_units_equivalent(
                k[e]["units"], units.dimensionless)
            assert "value" in k[e].keys()
            assert k[e]["value"] >= 0

        # Check for specific removal fractions
        if "removal_frac_mass_solute" in k.keys():
            for (j, c_data) in k["removal_frac_mass_solute"].items():
                assert "units" in c_data.keys()
                assert_units_equivalent(
                    c_data["units"], units.dimensionless)
                assert "value" in c_data.keys()
                assert c_data["value"] >= 0
                assert c_data["value"] <= 1
                assert j in db.component_list.keys()
