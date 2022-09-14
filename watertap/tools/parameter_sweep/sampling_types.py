###############################################################################
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

# sampling.py - This file contains all of the sampling classes

import numpy as np

from enum import Enum, auto
from abc import abstractmethod, ABC


class SamplingType(Enum):
    FIXED = auto()
    RANDOM = auto()
    RANDOM_LHS = auto()


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
        self.pyomo_object = pyomo_object
        self.setup(*args, **kwargs)

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


class GeomSample(FixedSample):
    def sample(self, num_samples):
        return np.geomspace(
            self.lower_limit, self.upper_limit, self.num_samples, endpoint=True
        )

    def setup(self, lower_limit, upper_limit, num_samples):
        self.lower_limit = lower_limit
        self.upper_limit = upper_limit
        self.num_samples = num_samples


class ReverseGeomSample(FixedSample):
    def sample(self, num_samples):
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


class LatinHypercubeSample(_Sample):
    sampling_type = SamplingType.RANDOM_LHS

    def sample(self, num_samples):
        return [self.lower_limit, self.upper_limit]

    def setup(self, lower_limit, upper_limit):
        self.lower_limit = lower_limit
        self.upper_limit = upper_limit
