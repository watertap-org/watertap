from watertap.tools.parameter_sweep import (LinearSample,
                                            UniformSample,
                                            NormalSample,
                                            LatinHypercubeSample)
import yaml

def _split_pyomo_path(model_value):
    """ Splits a dot-separated string into a valid model path

    This function attempts to split a dot-separated string
    into a list where each element is a valid string representation of
    a flowsheet object or attribute.  Since floating point values are
    valid indices for indexed parameters, the function joins together
    instances in which the decimal point separates a floating point index.

    For example::

         a = _split_pyomo_path("m.fs.RO.flux_mass_io_phase_comp[0.0,'in','Liq','H2O']")
         a = ["m", "fs", "RO", "flux_mass_io_phase_comp[0.0,'in','Liq','H2O']"]
         len(a) >> 4

    Args:
        model_value (str):
            A single dot-separated string representing the path to the model attribute.

    Returns:
        healed_path (list(str)):
            A list of strings where each element is a valid model object/attribute

    """

    path = model_value.split('.')

    healed_path = [path[0]]

    for k in range(1, len(path)):
        if path[k-1][-1].isnumeric() and path[k][0].isnumeric():
            # We accidentally split between a floating point number,
            # reattach the extra bit of string to the preceding one
            healed_path[-1] += '.' + path[k]
        else:
            # We split correctly, append the new element
            healed_path.append(path[k])

    return healed_path


def _return_child_property(m, healed_path):
    """ This function returns the attribute from a dot-separated object path

    This function returns an innermost nested attribute specified by a
    dot-separated object path.  For example::

        m.fs.RO.feed_side.properties_interface_out[0.0].conc_mass_phase_comp['Liq','H2O']

    returns the ``conc_mass_phase_comp['Liq','H2O']`` attribute.  In instances
    where the attribute or parameter is indexed, an attempt is made to extract
    the correct index to continue the descent along the corrent branch.

    Args:
        m (flowsheet):
            The Pyomo model containing the flowsheet.

        healed_path (list(str)):
            The healed path along the dot-separated specification, string attribute
            names stored as a list.

    Returns:
        child (Pyomo attribute):
            The innermost attribute from the specified path, e.g., ``m.fs.[...].return_attribute``.

    """

    for k, attr in enumerate(healed_path):
        if k == 0:
            parent = m

        else:
            found_key = None

            if '[' in attr:
                # This is an indexed parameter, we need
                # to temporarily remove the index.
                s = attr.split('[')
                attr = s[0]
                found_key = s[1].split(']')[0]
                found_key = found_key.replace('\'', '')
                found_key = found_key.replace('\"', '')

            try:
                child = getattr(parent, attr)
            except:
                raise ValueError('Could not acccess attribute %s' % (attr))

            if child.is_indexed():
                if found_key is not None:
                    # We have an indexed child, and a bracketed key specifier.
                    # Convert the found key to the expected type and use it for child[key].

                    index_type = type(next(child.keys()))

                    # Convert the found key to the expected type
                    if index_type is float:
                        key = float(found_key)

                    elif index_type is tuple:
                        key = []
                        found_key = found_key.split(',')

                        for idx in found_key:
                            try:
                                key.append(float(idx))
                            except:
                                key.append(idx)

                        key = tuple(key)

                    else:
                        print('Child indexing type %s not recognized.' % (index_type))

                else:
                    # We have an indexed child, but no bracketed index was found in the string specifier.
                    # Iterate through the child's keys and raise an error if >1 options exist.
                    # If only one option exists, use it as the key to descend further along the path. 
                    ct = 0
                    for key in child.keys():
                        ct += 1

                    if ct > 1:
                        raise ValueError('No bracketed index found for %s path specifier but multiple options exist.')

                child = child[key]

            parent = child

    return child


def get_sweep_params_from_yaml(m, yaml_filename):
    """ Creates a dictionary of swept model parameters specified via yaml file

    This function creates a dictionary of the items to vary during a parameter
    sweep where the variable name, model attribute, and sweeping domain are
    specified in a YAML file.  The YAML file should have the following format::

        A_comp:
            type: NormalSample
            model_value: m.fs.RO.A_comp
            mean: 4.0e-12
            std: 0.5e-12

    where the top-level keyword can be any short, easily understood identifier
    for the parameter.  ``type`` must be one of ``LinearSample``, ``UniformSample``,
    ``NormalSample``, or ``LatinHypercubeSample``.  ``model_value`` must be a valid
    dot-sperated string path to the object attribute (in this case, an RO attribute
    on the flowsheet ``m``) that you wish to vary.  The remaining arguments are
    dependent on the sample type selected.  For ``NormalSample`` information about
    the mean and standard deviation is required.  Consult the ``parameter_sweep``
    help for more information on the different sample classes.

    Args:
        m (pyomo model):
            The flowsheet containing the model to deploy with the parameter sweep
            tool.
        yaml_filename (str): The path to the yaml file.

    Returns:
        sweep_params (dict):
            A dictionary containing different instances of parameter sweep samples

    """

    try:
        # Open the yaml file and import the contents into a 
        # dictionary with the same structure
        with open(yaml_filename) as fp:
            input_dict = yaml.load(fp, Loader = yaml.FullLoader)
    except:
        raise ValueError('Could not open file %s' % (yaml_filename))

    sweep_params = {}

    for param, values in input_dict.items():


        # Split the dot-separated object path into a list of
        # separate objects and attributes
        path = _split_pyomo_path(values['model_value'])

        # Follow the attribute path to the correct nested
        # property, should work with indexing!
        child = _return_child_property(m, path)

        # TODO: assign the value for child based on theq
        # value stored in the dictionary
        # e.g., child.fix(values['fixed_val'])
        # or,   child.set_value(values['fixed_val'])


        if values['type'] == 'LinearSample':
            sweep_params[param] = LinearSample(child,
                                                values['lower_limit'],
                                                values['upper_limit'],
                                                values['num_samples'])

        elif values['type'] == 'UniformSample':
            sweep_params[param] = UniformSample(child,
                                                values['lower_limit'],
                                                values['upper_limit'])

        elif values['type'] == 'NormalSample':
            sweep_params[param] = NormalSample(child,
                                               values['mean'],
                                               values['std'])

        elif values['type'] == 'LatinHypercubeSample':
            sweep_params[param] = LatinHypercubeSample(child,
                                                       values['lower_limit'],
                                                       values['upper_limit'])

    return sweep_params
