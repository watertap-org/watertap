#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

# sampling.py - This file contains all of the sampling classes

import numpy as np

from enum import Enum, auto
from abc import abstractmethod, ABC


class SamplingType(Enum):
    FIXED = auto()
    RANDOM = auto()
    RANDOM_LHS = auto()


class SetMode(Enum):
    FIX_VALUE = auto()
    SET_LB = auto()
    SET_UB = auto()
    SET_FIXED_STATE = auto()


class _Sample(ABC):
    def __init__(self, pyomo_object, *args, **kwargs):
        # Check for indexed with single value
        if pyomo_object.is_indexed() and len(pyomo_object) == 1:
            for _data_obj in pyomo_object.values():
                pyomo_object = _data_obj

        # Make sure we are a Var() or Param()
        if not (pyomo_object.is_parameter_type() or pyomo_object.is_variable_type()):
            raise ValueError(
                f"The sweep parameter needs to be a pyomo Param or Var but {type(pyomo_object)} was provided instead."
            )
        if pyomo_object.is_parameter_type() and not pyomo_object.mutable:
            raise ValueError(
                f"Parameter {pyomo_object} is not mutable, and so cannot be set by parameter_sweep"
            )
        # default is always FIX_VALUE MODE
        self.set_mode = SetMode.FIX_VALUE

        self.pyomo_object = pyomo_object

        self.setup(*args, **kwargs)

    @abstractmethod
    def sample(self, num_samples):
        pass

    @abstractmethod
    def setup(self, *args, **kwargs):
        pass

    def set_variable_update_mode(self, set_mode_type=None, default_fixed_value=None):
        if set_mode_type is None:
            self.set_mode = SetMode.FIX_VALUE
        else:
            self.set_mode = set_mode_type
        self.default_fixed_value = default_fixed_value
        if self.pyomo_object.is_parameter_type() and self.set_mode != SetMode.FIX_VALUE:
            raise ValueError(
                f"Set fixed state for {self.pyomo_object} is not supported, SET_FIXED_STATE is only supported by pyomo Var"
            )


class RandomSample(_Sample):
    sampling_type = SamplingType.RANDOM


class FixedSample(_Sample):
    sampling_type = SamplingType.FIXED


class LinearSample(FixedSample):
    def sample(self):
        return np.linspace(self.lower_limit, self.upper_limit, self.num_samples)

    def setup(self, lower_limit, upper_limit, num_samples):
        self.lower_limit = lower_limit
        self.upper_limit = upper_limit
        self.num_samples = num_samples


class GeomSample(FixedSample):
    def sample(self):
        return np.geomspace(
            self.lower_limit, self.upper_limit, self.num_samples, endpoint=True
        )

    def setup(self, lower_limit, upper_limit, num_samples):
        self.lower_limit = lower_limit
        self.upper_limit = upper_limit
        self.num_samples = num_samples


class ReverseGeomSample(FixedSample):
    def sample(self):
        return (
            (self.upper_limit + self.lower_limit)
            - np.geomspace(
                self.lower_limit, self.upper_limit, self.num_samples, endpoint=True
            )
        )[::-1]

    def setup(self, lower_limit, upper_limit, num_samples):
        self.lower_limit = lower_limit
        self.upper_limit = upper_limit
        self.num_samples = num_samples


class PredeterminedFixedSample(FixedSample):
    """
    Similar to other fixed sampling types except the setup function arguments.
    In this case a user needs to specify a numpy array (or a list) of
    predetermined values. For example:

    sample_obj = PredeterminedFixedSample(np.array([1,2,3,4]))
    """

    def sample(self):
        return self.values

    def setup(self, values):
        self.values = values
        self.num_samples = len(values)


class UniformSample(RandomSample):
    def sample(self):
        return np.random.uniform(self.lower_limit, self.upper_limit, self.num_samples)

    def setup(self, lower_limit, upper_limit, num_samples):
        self.lower_limit = lower_limit
        self.upper_limit = upper_limit
        self.num_samples = num_samples


class NormalSample(RandomSample):
    def sample(self):
        return np.random.normal(self.mean, self.sd, self.num_samples)

    def setup(self, mean, sd, num_samples):
        self.mean = mean
        self.sd = sd
        self.num_samples = num_samples


class PredeterminedRandomSample(RandomSample):
    """
    Similar to other fixed sampling types except the setup function arguments.
    In this case a user needs to specify a numpy array (or a list) of
    predetermined values. For example:

    sample_obj = PredeterminedRandomSample(np.array([1,2,3,4]))
    """

    def sample(self):
        return self.values

    def setup(self, values):
        self.values = values
        self.num_samples = len(values)


class LatinHypercubeSample(_Sample):
    sampling_type = SamplingType.RANDOM_LHS

    def sample(self):
        return [self.lower_limit, self.upper_limit]

    def setup(self, lower_limit, upper_limit, num_samples):
        self.lower_limit = lower_limit
        self.upper_limit = upper_limit
        self.num_samples = num_samples
