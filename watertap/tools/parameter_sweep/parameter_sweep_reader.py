##############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################

from watertap.tools.parameter_sweep.sampling_types import (
    LinearSample,
    GeomSample,
    ReverseGeomSample,
    UniformSample,
    NormalSample,
    LatinHypercubeSample,
)
import yaml
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


class ParameterSweepReader:
    @staticmethod
    def _yaml_to_dict(yaml_filename):
        """Reads and stores a yaml file as a dictionary

        Args:
            yaml_filename (str):
                The filename of the yaml file to read.

        Returns:
            input_dict (dict):
                The result of reading the yaml file and translating
                its structure into a dictionary.

        """

        try:
            # Open the yaml file and import the contents into a
            # dictionary with the same structure
            with open(yaml_filename) as fp:
                input_dict = yaml.load(fp, Loader=yaml.FullLoader)

        except:
            raise ValueError("Could not open file %s" % (yaml_filename))

        return input_dict

    def get_sweep_params_from_yaml(self, m, yaml_filename):
        """Creates a dictionary of swept model parameters specified via yaml file

        This function creates a dictionary of the items to vary during a parameter
        sweep where the variable name, model attribute, and sweeping domain are
        specified in a YAML file.  The YAML file should have the following format::

            A_comp:
                type: NormalSample
                param: fs.RO.A_comp
                mean: 4.0e-12
                std: 0.5e-12

        where the top-level keyword can be any short, easily understood identifier
        for the parameter.  ``type`` must be one of ``LinearSample``, ``UniformSample``,
        ``NormalSample``, or ``LatinHypercubeSample``.  ``param`` must be a valid
        dot-sperated string path to the object attribute (in this case, an RO attribute
        on the flowsheet ``m``) that you wish to vary.  The remaining arguments are
        dependent on the sample type selected.  For ``NormalSample`` information about
        the mean and standard deviation is required.  Consult the ``parameter_sweep``
        help for more information on the different sample classes.

        Args:
            m (pyomo model):
                The flowsheet containing the model to deploy with the parameter sweep
                tool.
            yaml_filename (str):
                The path to the yaml file.

        Returns:
            sweep_params (dict):
                A dictionary containing different instances of parameter sweep samples

        """

        input_dict = self._yaml_to_dict(yaml_filename)
        return self._dict_to_params(m, input_dict)

    @staticmethod
    def _dict_to_params(m, input_dict):

        """Reads and stores a yaml file as a dictionary

        Args:
            dict (str):
                The dictionary of paramters that are turned into paramter sweep samples

        Returns:
            input_dict (dict):
                The result of reading the yaml file and translating
                its structure into a dictionary.

        """
        sweep_params = {}

        for param, values in input_dict.items():

            # Find the specified component on the model
            component = m.find_component(values["param"])

            if component is None:
                raise ValueError(f'Could not acccess attribute {values["param"]}')

            if values["type"] == "LinearSample":
                sweep_params[param] = LinearSample(
                    component,
                    values["lower_limit"],
                    values["upper_limit"],
                    values["num_samples"],
                )
            elif values["type"] == "GeomSample":
                sweep_params[param] = GeomSample(
                    component,
                    values["lower_limit"],
                    values["upper_limit"],
                    values["num_samples"],
                )
            elif values["type"] == "ReverseGeomSample":
                sweep_params[param] = ReverseGeomSample(
                    component,
                    values["lower_limit"],
                    values["upper_limit"],
                    values["num_samples"],
                )
            elif values["type"] == "UniformSample":
                sweep_params[param] = UniformSample(
                    component, values["lower_limit"], values["upper_limit"]
                )

            elif values["type"] == "NormalSample":
                sweep_params[param] = NormalSample(
                    component, values["mean"], values["std"]
                )

            elif values["type"] == "LatinHypercubeSample":
                sweep_params[param] = LatinHypercubeSample(
                    component, values["lower_limit"], values["upper_limit"]
                )

        return sweep_params

    @staticmethod
    def _set_value(component, key, default_value):
        _log.debug(f"Property: {key}")
        _log.debug(f"New Value: {default_value:.6e}, Old Value: {component.value:.6e}")
        component.value = default_value

    def set_defaults_from_yaml(self, m, yaml_filename, verbose=False):
        """Sets default model values using values stored in a yaml file

        This function reads a yaml file with the structure::

            fs.path.to.attribute_1: 0.123
            fs.path.to.attribute_2: 1.234
            ...

        and uses the (key, default_value) pairs to set default values
        for the attributes in model ``m``.

        Args:
            m (pyomo model):
                The flowsheet containing the model to set default values for
            yaml_filename (str):
                The path to the yaml file.

        Returns:
            N/A

        """

        input_dict = self._yaml_to_dict(yaml_filename)
        self._set_values_from_dict(m, input_dict, verbose)

    def _set_values_from_dict(self, m, input_dict, verbose=False):

        fail_count = 0

        for key, default_value in input_dict.items():
            # Find the specified component on the model
            component = m.find_component(key)

            if component is None:
                raise ValueError(f"Could not acccess component {key} on model {m}")

            if component.is_variable_type():
                self._set_value(component, key, default_value)
            elif component.is_parameter_type():
                if component.mutable:
                    self._set_value(component, key, default_value)
                else:
                    _log.warning(f"Cannot set value of non-mutable Param {component}")
                    fail_count += 1
            else:
                _log.warning(f"Cannot set value of component {component}")
                fail_count += 1

        number_defaults = len(input_dict)
        print(
            f"Set {number_defaults-fail_count} of {number_defaults} options to default values ({fail_count} failures, see warnings)"
        )


def get_sweep_params_from_yaml(m, yaml_filename):
    """Creates a dictionary of swept model parameters specified via yaml file

    This function creates a dictionary of the items to vary during a parameter
    sweep where the variable name, model attribute, and sweeping domain are
    specified in a YAML file.  The YAML file should have the following format::

        A_comp:
            type: NormalSample
            param: fs.RO.A_comp
            mean: 4.0e-12
            std: 0.5e-12

    where the top-level keyword can be any short, easily understood identifier
    for the parameter.  ``type`` must be one of ``LinearSample``, ``UniformSample``,
    ``NormalSample``, or ``LatinHypercubeSample``.  ``param`` must be a valid
    dot-sperated string path to the object attribute (in this case, an RO attribute
    on the flowsheet ``m``) that you wish to vary.  The remaining arguments are
    dependent on the sample type selected.  For ``NormalSample`` information about
    the mean and standard deviation is required.  Consult the ``parameter_sweep``
    help for more information on the different sample classes.

    Args:
        m (pyomo model):
            The flowsheet containing the model to deploy with the parameter sweep
            tool.
        yaml_filename (str):
            The path to the yaml file.

    Returns:
        sweep_params (dict):
            A dictionary containing different instances of parameter sweep samples

    """
    return ParameterSweepReader().get_sweep_params_from_yaml(m, yaml_filename)


def set_defaults_from_yaml(m, yaml_filename, verbose=False):
    """Sets default model values using values stored in a yaml file

    This function reads a yaml file with the structure::

        fs.path.to.attribute_1: 0.123
        fs.path.to.attribute_2: 1.234
        ...

    and uses the (key, default_value) pairs to set default values
    for the attributes in model ``m``.

    Args:
        m (pyomo model):
            The flowsheet containing the model to set default values for
        yaml_filename (str):
            The path to the yaml file.

    Returns:
        N/A

    """
    return ParameterSweepReader().set_defaults_from_yaml(m, yaml_filename, verbose)
