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

import idaes
from idaes.core.util.exceptions import BurntToast
from idaes.logger import solver_capture_off
from watertap.core.plugins.solvers import create_debug_solver_wrapper


class _ModelDebugMode:
    def __init__(self):
        self._prior_default_solver = None
        self._prior_idaes_logger_capture_solver = None

    def activate(self):
        _default_solver_config_value = idaes.cfg.get("default_solver")
        self._prior_default_solver = _default_solver_config_value.value()
        # create a debug solver around the current default solver
        debug_solver_name = create_debug_solver_wrapper(self._prior_default_solver)

        # reconfigure the default IDAES solver to use the debug wrapper
        _default_solver_config_value.set_default_value(debug_solver_name)
        if not _default_solver_config_value._userSet:
            _default_solver_config_value.reset()

        # disable solver log capturing so the resulting notebook
        # can use the whole terminal screen
        self._prior_idaes_logger_capture_solver = idaes.cfg.logger_capture_solver
        solver_capture_off()

    def deactivate(self):
        if self._prior_default_solver is None:
            if self._prior_idaes_logger_capture_solver is not None:
                raise BurntToast
            # TODO: should we raise an error instead?
            return

        # reconfigure the default IDAES solver to use the debug wrapper
        _default_solver_config_value = idaes.cfg.get("default_solver")
        _default_solver_config_value.set_default_value(self._prior_default_solver)
        if not _default_solver_config_value._userSet:
            _default_solver_config_value.reset()

        idaes.cfg.logger_capture_solver = self._prior_idaes_logger_capture_solver

        self._prior_idaes_logger_capture_solver = None
        self._prior_default_solver = None


_mdm = _ModelDebugMode()
activate = _mdm.activate
deactivate = _mdm.deactivate
