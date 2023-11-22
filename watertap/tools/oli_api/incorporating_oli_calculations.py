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
#
#################################################################################

from pyomo.environ import units as pyunits

from numpy import linspace

from watertap.tools.oli_api.util.state_block_helper_functions import (
    create_state_block,
    extract_state_vars,
)

from watertap.tools.oli_api.credentials import CredentialManager
from watertap.tools.oli_api.client import OLIApi

from watertap.tools.oli_api.flash import Flash

source_water = {
    "temperature": 298.15,
    "pressure": 101325,
    "components": {
        "Cl_-": 870,
        "Na_+": 739,
        "SO4_2-": 1011,
        "Mg_2+": 90,
        "Ca_2+": 258,
        "K_+": 9,
        "HCO3_-": 385,
    },
    "units": {
        "temperature": pyunits.K,
        "pressure": pyunits.Pa,
        "components": pyunits.mg / pyunits.L,
    },
}

# use source water in a WaterTAP model
m = create_state_block(source_water)
state_block = m.fs.stream[0]
conc_var = state_block.conc_mass_phase_comp
state_vars = extract_state_vars(state_block, conc_var, source_water["units"])
print(state_vars)

survey_conditions = {"SO4_2-": linspace(0, 1e3, 10)}

water_analysis = WaterAnalysis(
    state_vars=state_vars, survey_conditions=survey_conditions
)

water_analysis.oli_water_analysis_properties["AllowSolidsToForm"] = True
props = {
    "scalingIndex": False,
    "prescalingTendencies": True,
    "prescalingTendenciesRigorous": True,
    "scalingTendencies": True,
    "MBGComposition": False,
    "materialBalanceGroup": False,
}
water_analysis.oli_optional_properties.update(props)

# uncomment line below and input credentials
#credential_manager = CredentialManager()

# will take 20-30 seconds to run

survey = water_analysis.build_composition_survey(survey_conditions)

solute_list = source_water["components"]
phases = ["liquid1", "solid"]

with OLIApi(credential_manager) as oliapi:
    dbs_file_id = oliapi.get_dbs_file_id(
        chemistry_source=solute_list, phases=phases, model_name="remote_file_from_dict"
    )

    water_analysis.run(oliapi=oliapi, dbs_file_id=dbs_file_id)


print(water_analysis.survey)
print("\nPhase Properties:")

extracted_properties = water_analysis.extract_basic_properties(
    phase="liquid1", properties=["osmoticPressure", "ph"]
)
print(extracted_properties)

print("\nScaling Tendencies:")

extracted_scaling_tendencies = water_analysis.extract_scaling_tendencies(
    scalants=["CACO3", "CASO4.2H2O"]
)
print(extracted_scaling_tendencies)

print(water_analysis.inputs_true)
