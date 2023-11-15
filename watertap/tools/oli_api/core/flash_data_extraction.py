from pandas import MultiIndex, DataFrame

# TODO: Generalize for other flash calculations (currently tested for water analysis)
def extract_scaling_tendencies(raw_result=None, scalants=None, lower_bound=0):
    """
    Extracts scaling tendencies from OLI output for specific scalants.

    :param raw_result: dictionary containing raw data to extract from
    :param scalants: list containing names of scalants
    :param lower_bound: minimum scaling tendency to extract

    :return extracted_scaling_tendencies: copy of DataFrame containing extracted scaling tendencies
    """
       
    # TODO: make more informative (e.g., provide list of available scalants)
    if scalants is None:
        raise RuntimeError(
            f" Unable to find scaling tendency for species {scalants}."
        )
        
    header = MultiIndex.from_product(
        [scalants, ["prescaling", "eq. scaling"]], names=["species", "label"]
    )
    extracted_scaling_tendencies = DataFrame(
        columns=header, index=raw_result
    )
    for k in raw_result:
        root_path = raw_result[k]["result"]
        prescaling_path = root_path["additionalProperties"]["prescalingTendencies"][
            "values"
        ]
        eq_scaling_path = root_path["additionalProperties"]["scalingTendencies"][
            "values"
        ]
        for scalant in scalants:
            val = prescaling_path[scalant], eq_scaling_path[scalant]
            extracted_scaling_tendencies.loc[k, scalant] = val
    return extracted_scaling_tendencies

# TODO: probably condense these two methods into single method
def extract_basic_properties(raw_result=None, survey=None, phase="", properties=[]):
    """
    Extracts basic phase-specific properties from OLI output.

    :param phase: string name of phase to extract properties from
    :param properties: list containing string names of properties to extract from results

    :return extract: copy of DataFrame containing extracted properties
    """
        
    header = MultiIndex.from_product(
        [properties, ["value", "unit"]], names=["property", "label"]
    )
    extracted_properties = DataFrame(
        columns=header, index=raw_result
    )
    for k in raw_result:
        root_path = raw_result[k]["result"]
        for prop in properties:
            val = (
                root_path["phases"][phase]["properties"][prop]["value"],
                root_path["phases"][phase]["properties"][prop]["unit"],
            )
            extracted_properties.loc[k, prop] = val
    return extracted_properties