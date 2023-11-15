from numpy import linspace, array

import yaml
from copy import deepcopy
from itertools import product
from pandas import DataFrame, MultiIndex
from pyomo.environ import value, units as pyunits

from watertap.tools.oli_api.util.fixed_keys_dict import (
    default_oli_water_analysis_properties,
    default_oli_optional_properties,
    default_oli_unit_set_info,
)

from watertap.tools.oli_api.util.watertap_to_oli_helper_functions import (
    get_oli_name,
)

# TODO: add zero_species argument so users can track specific additives/scalants - not sure if necessary

def build_survey(survey_vars={}, get_oli_names=False, tee=False):
    """
    Builds a DataFrame used to modify flash calculation parameters.

    :param survey_vars: dictionary containing variables: arrays to survey
    :param get_oli_names: boolean switch to convert name into OLI form
    
    :return survey: DataFrame containing surveys
    """
    
    exclude_items = ["Temperature", "Pressure"]
    if bool(survey_vars):
        survey_vars = {(get_oli_name(k) if bool(get_oli_names) and (k not in exclude_items) else k): v for k, v in survey_vars.items()}
        survey_prod = list(product(*(survey_vars[key] for key in survey_vars)))
        survey = DataFrame(columns=survey_vars.keys(), index=range(len(survey_prod)), data=survey_prod)
        if bool(tee):
            print(f"Number of survey conditions: {len(survey)}.")
        return survey
    else:
        return DataFrame()

def modify_inputs(initial_flash_input, survey, flash_method=""):
    """
    Iterates over a survey to create modified clones of an initial flash analysis output.

    :param initial_flash_input: flash analysis input to copy
    :param survey: DataFrame containing modifications for each test
    :param flash_method: string name of flash method to use

    :return clones: dictionary containing modified state variables and survey index
    """

    clones = {}    
    for clone_index in survey.index:
        modified_clone = deepcopy(initial_flash_input)
        if flash_method == "wateranalysis":
            for param in modified_clone["params"]["waterAnalysisInputs"]:
                key = param["name"]
                if key in survey.columns:
                    param.update({"value": survey.loc[clone_index, key]})
        # TODO: figure out how to modify water analysis output (look at Mayo's work for simplicity)
        elif flash_method == "isothermal":
            for param in modified_clone:
                print(param in survey.columns)
        else:
            raise IOError(" Flash calculations besides 'wateranalysis' and 'isothermal' not yet implemented.")
        clones[clone_index] = modified_clone
    return clones

