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

def build_survey(survey_conditions={}, get_oli_name=False):
    """
    Builds a list of modified clone dictionaries.

    :param survey_conditions: dictionary containing variables: arrays to survey
    :param get_oli_name: boolean switch to convert name into OLI form
    
    :return survey: DataFrame containing parameters to survey
    """
    
    if bool(survey_conditions):
        if bool(get_oli_name):
            survey_conditions = {get_oli_name(k): v for k, v in survey_conditions.items()}
        survey_product = list(product(*(survey_conditions[key] for key in survey_conditions)))
        return array(survey_product)
    else:
        return array()
    
def modify_inputs(base_case, survey, flash_method=""):
    """
    Iterates over a dataframe to create a modified clone of an input dictionary.

    :param base_case: 
    :param survey: DataFrame containing modifications for each test
    :param flash_method: string name of flash method to use

    :return inputs_clone: dictionary containing modified state variables
    """

    # TODO: distinguish between wateranalysis and other flashes
    if flash_method == "wateranalysis":
        inputs_clone = deepcopy(base_case)
        
        for param in inputs_clone:
            key = param["name"]
            if key in survey:
                param.update({"value": survey[key].iloc[i]})
        return inputs_clone
    else:
        return {}
