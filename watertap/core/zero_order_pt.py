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
"""
This module contains the base class for all zero order pass-through unit
models (i.e. units with a single inlet and single outlet where flow and
composition do not change, such as pumps).
"""
import idaes.logger as idaeslog
from idaes.core.util.tables import create_stream_table_dataframe

from watertap.core.zero_order_base import ZeroOrderBaseData

# Some more inforation about this module
__author__ = "Andrew Lee"

# Set up logger
_log = idaeslog.getLogger(__name__)


class PassThroughBaseData(ZeroOrderBaseData):
    """
    Standard base class for pass-through unit models.

    This class is intended to be used for creating derived model classes and
    should not be instantiated by itself.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        # Create state blocks for inlet and outlets
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["defined_state"] = True

        self.properties = self.config.property_package.build_state_block(
            self.flowsheet().time,
            doc="Material properties in unit",
            default=tmp_dict)

        # Create Ports
        self.add_port("inlet", self.properties, doc="Inlet port")
        self.add_port("outlet", self.properties, doc="Tutlet port")

    def initialize(blk, state_args=None, outlvl=idaeslog.NOTSET,
                   solver=None, optarg=None):
        '''
        Initialization routine for pass-through unit models.

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                           package(s) to provide an initial state for
                           initialization (see documentation of the specific
                           property package) (default = {}).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None, use default IDAES solver)

        Returns:
            None
        '''
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")

        # Get initial guesses for inlet if none provided
        if state_args is None:
            state_args = {}
            state_dict = (
                blk.properties[
                    blk.flowsheet().time.first()]
                .define_port_members())

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        # ---------------------------------------------------------------------
        # Initialize control volume block
        blk.properties.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
            hold_state=False
        )
        init_log.info('Initialization Complete.')

    def calculate_scaling_factors(self):
        # Get default scale factors and do calculations from base classes
        super().calculate_scaling_factors()

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {"Inlet": self.inlet,
             "Outlet": self.outlet},
            time_point=time_point)
