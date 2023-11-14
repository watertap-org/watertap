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

# planned contents:
# water analysis input builder
# flash case input builder
# helper functions
# extraction functions

class Flash:
    def __init__(self,
                 optional_properties=default_oli_optional_properties,
                 unit_set_info=default_oli_unit_set_info,
                 water_analysis_properties=default_oli_water_analysis_properties):
        
        # set values based on inputs
        self.optional_properties = optional_properties
        self.unit_set_info = unit_set_info
        self.water_analysis_properties = water_analysis_properties
        
    def build_water_analysis_input(self, state_vars={}):
        """
        Creates input list for water-analysis flash method.

        :param state_vars: dictionary containing state variables
        """
        
        if not bool(state_vars): 
            raise IOError(
                " Provide a dictionary of state variables with values," + 
                " and units for each variable."
            )
            
        self.water_analysis_input = []
        
        self.water_analysis_input.append(
            {
                "group": "Properties",
                "name": "Temperature",
                "unit": str(state_vars["units"]["temperature"]),
                "value": float(state_vars["temperature"]),
            }
        )
        self.water_analysis_input.append(
            {
                "group": "Properties",
                "name": "Pressure",
                "unit": str(state_vars["units"]["pressure"]),
                "value": float(state_vars["pressure"]),
            }
        )
        
        for comp, conc in state_vars["components"].items():
            charge = get_charge(comp)
            self.water_analysis_input.append(
                {
                    "group": get_charge_group(charge),
                    "name": get_oli_name(comp),
                    # TODO: need this case output 
                    "unit": "mg/L", #str(state_vars["units"]["components"]),
                    "value": float(conc),
                    "charge": charge,
                }
            )

            electroneutrality_options = [
                "DominantIon",
                "ProrateCations",
                "ProrateAnions",
                "Prorate",
                "AutoNACL",
                "MakeupIon",
            ]
            self.water_analysis_properties._check_value("electroneutrality_value", electroneutrality_options)

            self.water_analysis_input.append(
                {
                    "group": "Electroneutrality Options",
                    "name": "ElectroNeutralityBalanceType",
                    "value": self.water_analysis_properties["electroneutrality_value"],
                }
            )
            if self.water_analysis_properties["electroneutrality_value"] == "MakeupIon":
                makeup_ion_options = [get_oli_name(comp) for comp in state_vars["components"]]
                self.water_analysis_properties._check_value("MakeupIonBaseTag", makeup_ion_options)    
                
                self.water_analysis_input.append(
                    {
                        "group": "Electroneutrality Options",
                        "name": "MakeupIonBaseTag",
                        "value": self.water_analysis_properties["MakeupIonBaseTag"],
                    }
                )

            reconciliation_options = [
                "EquilCalcOnly",
            ]
            # TODO: test additional reconciliation options
            """
                "ReconcilePh",
                "ReconcilePhAndAlkalinity",
                "ReconcilePhAndAlkalinityAndTic",
                "ReconcileCo2Gas",
            ]
            """
            
            self.water_analysis_properties._check_value("reconciliation_value", reconciliation_options)
            
            self.water_analysis_input.append(
                {
                    "group": "Calculation Options",
                    "name": "CalcType",
                    "value": self.water_analysis_properties["reconciliation_value"],
                }
            )
            # there is an option in the API to include/exclude specific solids if desired:
            # https://devdocs.olisystems.com/optional-inputs
            self.water_analysis_input.append(
                {
                    "group": "Calculation Options",
                    "name": "AllowSolidsToForm",
                    "value": bool(self.water_analysis_properties["AllowSolidsToForm"]),
                }
            )
            self.water_analysis_input.append(
                {
                    "group": "Calculation Options",
                    "name": "CalcAlkalnity",
                    "value": bool(self.water_analysis_properties["CalcAlkalnity"]),
                }
            )
            self.water_analysis_properties._check_value("PhAcidTitrant", ["HCL", "CO2"])
            self.water_analysis_input.append(
                {
                    "group": "Calculation Options",
                    "name": "PhAcidTitrant",
                    "value": self.water_analysis_properties["PhAcidTitrant"],
                }
            )
            self.water_analysis_properties._check_value("PhBaseTitrant", ["NAOH"])
            self.water_analysis_input.append(
                {
                    "group": "Calculation Options",
                    "name": "PhBaseTitrant",
                    "value": self.water_analysis_properties["PhBaseTitrant"],
                }
            )
            
            # reserved for further testing
            """
            if self.oli_input_dict["reconciliation_value"] == "ReconcileCo2Gas":
                self.oli_input_dict._check_value("CO2GasFraction", range(0, 101))
                self.inputs_true.append(
                    {
                        "group": "Properties",
                        "name": "CO2GasFraction",
                        "unit": self.oli_input_dict["gas_fraction_unit"],
                        "value": self.oli_input_dict["CO2GasFraction"],
                    }
                )
            else:
                if "Ph" in self.oli_input_dict["reconciliation_value"]:
                    self.oli_input_dict._check_value("pH", range(0, 15))
                    self.inputs_true.append(
                        {
                            "group": "Properties",
                            "name": "pH",
                            "value": self.oli_input_dict["pH"],
                        }
                    )

                if "Alkalinity" in self.oli_input_dict["reconciliation_value"]:
                    self.inputs_true.append(
                        {
                            "group": "Properties",
                            "name": "Alkalinity",
                            "unit": self.oli_input_dict["alkalinity_unit"],
                            "value": float(self.oli_input_dict["Alkalinity"]),
                        }
                    )
                    # TODO: which other titrants and pH endpoints may be used
                    self.oli_input_dict._check_value("AlkalinityPhTitrant", ["H2SO4"])
                    self.inputs_true.append(
                        {
                            "group": "Calculation Options",
                            "name": "AlkalinityPhTitrant",
                            "value": self.oli_input_dict["AlkalinityPhTitrant"],
                        }
                    )
                    self.oli_input_dict._check_value(
                        "AlkalinityTitrationEndpointPh", range(0, 15)
                    )
                    self.inputs_true.append(
                        {
                            "group": "Properties",
                            "name": "AlkalinityTitrationEndPointPh",
                            "value": self.oli_input_dict["AlkalinityTitrationEndPointPh"],
                        }
                    )

                if "Tic" in self.oli_input_dict["reconciliation_value"]:
                    self.inputs_true.append(
                        {
                            "group": "Properties",
                            "name": "TIC",
                            "unit": self.oli_input_dict["tic_unit"],
                            "value": float(self.oli_input_dict["TIC"]),
                        }
                    )        
            """
            
    def build_flash_calculation_input(self, method="", state_vars={}, additional_params={}, water_analysis_output=None):
        
        if not bool(water_analysis_output):
            raise IOError("Run wateranalysis flash to generate water_analysis_output data.")
        
        self.flash_analysis_inputs = {
            "params": {
                "temperature":  {
                    "unit": state_vars["units"]["temperature"],
                    "value": state_vars["temperature"]
                },
                "pressure": {
                    "unit": state_vars["units"]["pressure"],
                    "value": state_vars["pressure"]
                },
                "inflows": {self.extract_inflows(water_analysis_output)}
                }
            }
        # TODO: enable other flash functions by updating flash_analysis_inputs with additional_params required for flash
        #self.flash_analysis_inputs.update(self.extract_additional_params(method, additional_params))
        self.flash_analysis_inputs.update({"unitSetInfo": self.unit_set_info})
        
    def extract_additional_params(self, method, additional_params):
        return {}
    
    def extract_inflows(self, water_analysis_output):
        return water_analysis_output["result"]["total"]["molecularConcentration"]
    
