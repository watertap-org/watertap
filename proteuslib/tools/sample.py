###############################################################################
# ProteusLib Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/nawi-hub/proteuslib/"
#
###############################################################################
import numpy as np

from enum import Enum, auto
from abc import abstractmethod, ABC 


class SamplingType(Enum):
    FIXED = auto()
    RANDOM = auto()


class _Sample(ABC): 

    def __init__(self, pyomo_object, *args, **kwargs):
        # Check for indexed with single value
        if pyomo_object.is_indexed() and len(pyomo_object) == 1:
            pyomo_object = pyomo_object.values()[0]

        # Make sure we are a Var() or Param()
        if not (pyomo_object.is_parameter_type() or pyomo_object.is_variable_type()):
            raise ValueError(f"The sampled parameter needs to be a pyomo Param or Var but {type(pyomo_object)} was provided instead.")

        # Set the appropriate method for value updates
        if pyomo_object.is_parameter_type():
            self.pyo_set_value_method = pyomo_object.set_value
        else:
            self.pyo_set_value_method = pyomo_object.fix

        self._pyomo_object = pyomo_object 
        self.setup(*args, **kwargs)

    def set_value(self, val):
        self.pyo_set_value_method(val)

    @abstractmethod 
    def sample(self, num_samples): 
        pass 

    @abstractmethod 
    def setup(self, *args, **kwargs): 
        pass 


class RandomSample(_Sample):
    sampling_type = SamplingType.RANDOM


class FixedSample(_Sample):
    sampling_type = SamplingType.FIXED


class LinearSample(FixedSample):

    def sample(self, num_samples): 
        return np.linspace(self.lower_limit, self.upper_limit, self.num_samples)

    def setup(self, lower_limit, upper_limit, num_samples):
        self.lower_limit = lower_limit
        self.upper_limit = upper_limit
        self.num_samples = num_samples


class UniformSample(RandomSample):

    def sample(self, num_samples): 
        return np.random.uniform(self.lower_limit, self.upper_limit, num_samples)

    def setup(self, lower_limit, upper_limit):
        self.lower_limit = lower_limit
        self.upper_limit = upper_limit


class NormalSample(RandomSample):

    def sample(self, num_samples): 
        return np.random.normal(self.mean, self.sd, num_samples)

    def setup(self, mean, sd):
        self.mean = mean
        self.sd = sd
