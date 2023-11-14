###############################################################################
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
# OLI Systems, Inc. Copyright © 2022, all rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation and/or
# other materials provided with the distribution.
#
# 3. Neither the name of OLI Systems, Inc. nor the names of any contributors to
# the software made available herein may be used to endorse or promote products derived
# from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
# SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
# OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
# TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# You are under no obligation whatsoever to provide any bug fixes, patches, or upgrades to the
# features, functionality or performance of the source code ("Enhancements") to anyone; however,
# if you choose to make your Enhancements available either publicly, or directly to OLI Systems, Inc.,
# without imposing a separate written license agreement for such Enhancements, then you hereby grant
# the following license: a non-exclusive, royalty-free perpetual license to install, use, modify, prepare
# derivative works, incorporate into other computer software, distribute, and sublicense such enhancements
# or derivative works thereof, in binary and source code form.
###############################################################################

__author__ = "Oluwamayowa Amusat, Paul Vecchiarelli"

from pyomo.environ import value, units as pyunits

import yaml
from copy import deepcopy
from itertools import product
from pandas import DataFrame, MultiIndex

from numpy import linspace

from watertap.tools.oli_api.util.state_block_helper_functions import (
    create_state_block,
    extract_state_vars,
)

from watertap.tools.oli_api.credentials import CredentialManager
from watertap.tools.oli_api.client import OLIApi

#from watertap.tools.oli_api.core.water_analysis import WaterAnalysis

from watertap.tools.oli_api.util.fixed_keys_dict import (
    default_oli_water_analysis_properties,
    default_oli_optional_properties,
    default_oli_unit_set_info,
)

from watertap.tools.oli_api.util.watertap_to_oli_helper_functions import (
    get_oli_name,
    get_charge,
    get_charge_group,
    get_molar_mass,
)

from watertap.tools.oli_api.util.survey import (
    build_survey)

# potential contents:
# flash cases
# helper functions
# extraction functions

# TODO: generate a base case using water analysis (maybe chemistry-info also works?)

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

f = Flash()
f.build_water_analysis_input(source_water)

f.water_analysis_properties["AllowSolidsToForm"] = True
props = {
    "scalingIndex": False,
    "prescalingTendencies": True,
    "prescalingTendenciesRigorous": True,
    "scalingTendencies": True,
    "MBGComposition": False,
    "materialBalanceGroup": False,
}
f.optional_properties.update(props)

water_analysis_inputs = {"params": {
    "waterAnalysisInputs": f.water_analysis_input},
    "optionalProperties": dict(f.optional_properties),
    "unitSetInfo": dict(f.unit_set_info)}

credential_manager = CredentialManager(encryption_key=encryption_key)

with OLIApi(credential_manager) as oliapi:
    dbs_file_id = oliapi.get_dbs_file_id(chemistry_source=source_water["components"],
                                         phases=["liquid1", "solid"],
                                         model_name="remote_file_from_dict")                        
        
    water_analysis_base_case = oliapi.call("wateranalysis", dbs_file_id, water_analysis_inputs)
    
# after generating a base case, we can modify it and then do a composition survey of water analyses
survey_conditions = {"SO4_2-": linspace(0, 1e3, 3),
                     "Cl_-": linspace(0, 1e4, 4),
                     "Na_+": linspace(0, 1e4, 4),
                     "Mg_2+": linspace(0, 1e3, 3),
                     "temperature": linspace(0, 1e2, 10)}
water_analysis_survey = build_survey(survey_conditions)

print(water_analysis_survey)

# we can do further flash analyses using the water analysis base case 

#f.extract_inflows(water_analysis_base_case)

# we can do a survey of flash analyses