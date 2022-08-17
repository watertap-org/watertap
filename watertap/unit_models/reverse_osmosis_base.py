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

from copy import deepcopy
from enum import Enum, auto
from pyomo.environ import (
    Block,
    NonNegativeReals,
    Param,
    Suffix,
    Var,
    check_optimal_termination,
    exp,
    units as pyunits,
)
from idaes.core import UnitModelBlockData
from idaes.core.solvers import get_solver
from idaes.core.util.exceptions import ConfigurationError, InitializationError
from idaes.core.util.tables import create_stream_table_dataframe
import idaes.logger as idaeslog

from watertap.core.membrane_channel_base import validate_membrane_config_args, CONFIG_Template


class ReverseOsmosisBaseData(UnitModelBlockData):
    """
    Reverse Osmosis base class
    """

    CONFIG = CONFIG_Template() 

    def build(self):
        """
        Common variables and constraints for an RO unit model
        """
        super().build()

        if len(self.config.property_package.solvent_set) > 1:
            raise ConfigurationError(
                "Membrane models only support one solvent component,"
                "the provided property package has specified {} solvent components".format(
                    len(self.config.property_package.solvent_set)
                )
            )

        if len(self.config.property_package.phase_list) > 1 or "Liq" not in [
            p for p in self.config.property_package.phase_list
        ]:
            raise ConfigurationError(
                "Membrane models only support one liquid phase ['Liq'],"
                "the property package has specified the following phases {}".format(
                    [p for p in self.config.property_package.phase_list]
                )
            )

        validate_membrane_config_args(self)

        self._add_membrane_channel()

        self.membrane_channel.add_geometry()

        self.membrane_channel.add_state_blocks(has_phase_equilibrium=False)

        self.membrane_channel.add_material_balances(
            balance_type=self.config.material_balance_type, has_mass_transfer=True
        )

        self.membrane_channel.add_energy_balances(
            balance_type=self.config.energy_balance_type, has_enthalpy_transfer=True
        )

        self.membrane_channel.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            pressure_change_type=self.config.pressure_change_type,
            has_pressure_change=self.config.has_pressure_change,
        )

        self.membrane_channel.add_isothermal_conditions()

        self.membrane_channel.add_volumetric_flowrate_balance()

        self.membrane_channel.add_flux_balance()

        self.membrane_channel.add_mass_transfer()

        self.membrane_channel.add_concentration_polarization(
            concentration_polarization_type=self.config.concentration_polarization_type,
            mass_transfer_coefficient=self.config.mass_transfer_coefficient,
        )

        self.membrane_channel.add_recovery_vol_phase()

        self.membrane_channel.add_expressions()

        self.membrane_channel.apply_transformation()

        # Add Ports
        self.add_inlet_port(name="inlet", block=self.membrane_channel)
        self.add_outlet_port(name="retentate", block=self.membrane_channel)
        self.add_port(name="permeate", block=self.membrane_channel.mixed_permeate)

        # TODO: add references

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

    def initialize_build(
        self,
        initialize_guess=None,
        state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        General wrapper for RO initialization routines

        Keyword Arguments:

            initialize_guess : a dict of guesses for solvent_recovery, solute_recovery,
                               and cp_modulus. These guesses offset the initial values
                               for the retentate, permeate, and membrane interface
                               state blocks from the inlet feed
                               (default =
                               {'deltaP': -1e4,
                               'solvent_recovery': 0.5,
                               'solute_recovery': 0.01,
                               'cp_modulus': 1.1})
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for the inlet
                         feed side state block (see documentation of the specific
                         property package) (default = None).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None)
            solver : solver object or string indicating which solver to use during
                     initialization, if None provided the default solver will be used
                     (default = None)
        Returns:
            None
        """
        self.membrane_channel.initialize(state_args=state_args, outlvl=outlvl, optarg=optarg, solver=solver, initialize_guess=initialize_guess) 
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        # Create solver
        opt = get_solver(solver, optarg)
        # ---------------------------------------------------------------------
        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
            if not check_optimal_termination(res):
                raise InitializationError(f"The RO unit {self.name} failed to initialize")


    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Feed Inlet": self.inlet,
                "Feed Outlet": self.retentate,
                "Permeate Outlet": self.permeate,
            },
            time_point=time_point,
        )
