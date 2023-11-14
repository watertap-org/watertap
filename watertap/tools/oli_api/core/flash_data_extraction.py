# TODO: Make index labels more clear
# TODO: Generalize for other calls
def extract_scaling_tendencies(scalants=None, lower_bound=0):
    """
    Extracts scaling tendencies from OLI output for specific scalants.

    :param scalants: list containing names of scalants
    :param lower_bound: minimum scaling tendency to extract

    :return extract: copy of DataFrame containing extracted scaling tendencies
    """

    if scalants is None:
        raise RuntimeError(
            f" Unable to find scaling tendency for species {scalants}."
        )
    header = MultiIndex.from_product(
        [scalants, ["prescaling", "eq. scaling"]], names=["species", "label"]
    )
    self.extracted_scaling_tendencies = DataFrame(
        columns=header, index=range(len(self.results))
    )
    for i in range(len(self.results)):
        root_path = self.results[i]["result"]
        prescaling_path = root_path["additionalProperties"]["prescalingTendencies"][
            "values"
        ]
        eq_scaling_path = root_path["additionalProperties"]["scalingTendencies"][
            "values"
        ]
        for scalant in scalants:
            val = prescaling_path[scalant], eq_scaling_path[scalant]
            self.extracted_scaling_tendencies.loc[i, scalant] = val
    extract = deepcopy(self.extracted_scaling_tendencies)
    return extract

def extract_basic_properties(self, phase, properties):
    """
    Extracts basic phase-specific properties from OLI output.

    :param phase: string name of phase to extract properties from
    :param properties: list containing string names of properties to extract from results

    :return extract: copy of DataFrame containing extracted properties
    """

    header = MultiIndex.from_product(
        [properties, ["value", "unit"]], names=["property", "label"]
    )
    self.extracted_properties = DataFrame(
        columns=header, index=range(len(self.results))
    )
    for i in range(len(self.results)):
        root_path = self.results[i]["result"]
        for prop in properties:
            val = (
                root_path["phases"][phase]["properties"][prop]["value"],
                root_path["phases"][phase]["properties"][prop]["unit"],
            )
            self.extracted_properties.loc[i, prop] = val
    extract = deepcopy(self.extracted_properties)
    return extract
    
extract_scaling_tendencies()
