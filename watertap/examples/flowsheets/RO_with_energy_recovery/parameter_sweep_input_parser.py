from watertap.tools.parameter_sweep import (LinearSample,
                                            UniformSample,
                                            NormalSample,
                                            LatinHypercubeSample)
import yaml

def _split_pyomo_path(model_value):
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

def get_sweep_params_from_yaml(m, yaml_filename):
    try:
        # Open the yaml file and import the contents into a 
        # dictionary with the same structure
        with open(yaml_filename) as fp:
            input_dict = yaml.load(fp, Loader = yaml.FullLoader)
    except:
        raise ValueError('Could not open file %s' % (yaml_filename))

    sweep_params = {}

    for param, values in input_dict.items():

        path = _split_pyomo_path(values['model_value'])

        for k, pp in enumerate(path):
            if k == 0:
                child = m
            else:
                try:
                    child = getattr(child, pp)
                except:
                    raise ValueError('Could not acccess attribute %s' % (pp))

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
