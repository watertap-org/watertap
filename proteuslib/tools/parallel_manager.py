import numpy as np

def build_and_divide_combinations(d, rank, num_procs, return_global=False):

    param_values = []
    param_keys = []

    for k, v in d.items():
        # Build a vector of discrete values for this parameter
        # and record the parameter's name
        # v[0] = start, v[1] = stop, v[2] = resolution (number of elements)
        p = np.linspace(v[0], v[1], v[2])
        param_values.append(p)
        param_keys.append(k)

    num_var_params = len(param_values)

    # Form an array with every possible combination of parameter values
    full_combo_array = np.array(np.meshgrid(*param_values))
    full_combo_array = full_combo_array.T.reshape(-1, num_var_params)

    # Split the total list of combinations into NUM_PROCS chunks,
    # one per each of the MPI ranks
    divided_combo_array = np.array_split(full_combo_array, num_procs, axis=0)

    # Return only this rank's portion of the total workload
    local_combo_array = divided_combo_array[rank]

    return param_keys, local_combo_array, full_combo_array

def set_nested_attr(m, full_attr_path, value):
    
    parent = m
    attr_list = full_attr_path.split('.')
    stripped_list = []
    nn = len(attr_list)
    
    for attr in attr_list:
        if '[' in attr:
            attr_name = attr.split('[')[0]
            idx = int(attr.split('[')[1].split(']')[0])
        else:
            attr_name = attr
            idx = None
            
        stripped_list.append([attr_name, idx])
    
    for k, item in enumerate(stripped_list):
        attr_name = item[0]
        idx = item[1]
        child = getattr(parent, attr_name)

        if type(child) == list and idx is None:
            raise ValueError('The "%s" attribute is a list and must have an element specified with an integer index' % attr_name)
        
        if k < nn-1:
            if idx is not None:
                child = child[idx]
    
            parent = child

        else:
            if idx is not None:
                child[idx] = value
            else:
                setattr(parent, attr_name, value)

    return m
