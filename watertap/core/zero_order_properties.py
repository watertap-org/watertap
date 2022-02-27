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
This module contains the general purpose property package for zero-order
unit models. Zero-order models do not track temperature and pressure, or any
form of energy flow.
"""
from idaes.core import (EnergyBalanceType,
                        MaterialBalanceType,
                        MaterialFlowBasis,
                        PhysicalParameterBlock,
                        StateBlock,
                        StateBlockData,
                        declare_process_block_class)
from idaes.core.components import Solvent, Solute
from idaes.core.phases import LiquidPhase
from idaes.core.util.initialization import fix_state_vars, revert_state_vars, solve_indexed_blocks
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util import get_solver

from pyomo.environ import (Expression,
                           Param,
                           PositiveReals,
                           units as pyunits,
                           Var,
                           check_optimal_termination,
                           value)
from pyomo.common.config import ConfigValue

# Some more information about this module
__author__ = "Andrew Lee"

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("WaterParameterBlock")
class WaterParameterBlockData(PhysicalParameterBlock):
    """
    Property Parameter Block Class

    Defines component and phase lists, along with base units and constant
    parameters.
    """
    CONFIG = PhysicalParameterBlock.CONFIG()

    CONFIG.declare('database', ConfigValue(
        description='An instance of a WaterTAP Database to use for parameters.'
        ))
    CONFIG.declare('water_source', ConfigValue(
        description=
        'Water source to use when looking up parameters from database.'))
    CONFIG.declare("solute_list", ConfigValue(
        domain=list,
        description="List of solute species of interest. If None, will use "
        "all species defined in the water_source provided."))

    def build(self):
        '''
        Callable method for Block construction.
        '''
        super().build()

        self._state_block_class = WaterStateBlock

        self.Liq = LiquidPhase()

        self.H2O = Solvent()

        # Get component set from database if provided
        comp_set = None
        if self.config.database is not None:
            comp_set = self.config.database.get_solute_set(
                self.config.water_source)

        # Check definition of solute list
        solute_list = self.config.solute_list
        if solute_list is None:
            # No user-provided solute list, look up list from database
            if comp_set is None:
                # No solute list in database and none provided.
                raise ConfigurationError(
                    f"{self.name} no solute_list or database was defined. "
                    f"Users must provide at least one of these arguments.")
            else:
                solute_list = comp_set
        elif self.config.database is not None:
            # User provided custom list and database - check that all
            # components are supported
            for j in solute_list:
                if j not in comp_set:
                    _log.info(f"{self.name} component {j} is not defined in "
                              f"the water_sources database file.")
        else:
            # User provided list but no database - assume they know what they
            # are doing
            pass

        for j in solute_list:
            self.add_component(str(j), Solute())

        # Define default value for mass density of solution
        self.dens_mass_default = 1000*pyunits.kg/pyunits.m**3
        # Define default value for dynamic viscosity of solution
        self.visc_d_default = 0.001*pyunits.kg/pyunits.m/pyunits.s


        # ---------------------------------------------------------------------
        # Set default scaling factors
        self.default_scaling_factor = {
            ("flow_vol"): 1e3,
            ("conc_mass_comp"): 1e2}

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units({
                'time': pyunits.s,
                'length': pyunits.m,
                'mass': pyunits.kg,
                'amount': pyunits.mol,
                'temperature': pyunits.K,
                })

        obj.add_properties(
            {'flow_mass_comp': {'method': None},
             'flow_vol': {'method': '_flow_vol'},
             'conc_mass_comp': {'method': '_conc_mass_comp'},
             'dens_mass': {'method': '_dens_mass'},
             'visc_d': {'method': '_visc_d'}})


class _WaterStateBlock(StateBlock):
    """
    This Class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """
    def initialize(blk,
                   state_args=None,
                   state_vars_fixed=False,
                   hold_state=False,
                   outlvl=idaeslog.NOTSET,
                   solver=None,
                   optarg=None):
        '''
        Initialization routine for property package.

        Keyword Arguments:
        state_args : Dictionary with initial guesses for the state vars
                     chosen. Note that if this method is triggered
                     through the control volume, and if initial guesses
                     were not provied at the unit model level, the
                     control volume passes the inlet values as initial
                     guess.The keys for the state_args dictionary are:

                     flow_mol_comp : value at which to initialize component
                                     flows (default=None)
                     pressure : value at which to initialize pressure
                                (default=None)
                     temperature : value at which to initialize temperature
                                  (default=None)
        outlvl : sets output level of initialization routine
        state_vars_fixed: Flag to denote if state vars have already been
                          fixed.
                          - True - states have already been fixed and
                                   initialization does not need to worry
                                   about fixing and unfixing variables.
                         - False - states have not been fixed. The state
                                   block will deal with fixing/unfixing.
        optarg : solver options dictionary object (default=None, use
                 default solver options)
        solver : str indicating which solver to use during
                 initialization (default = None, use default solver)
        hold_state : flag indicating whether the initialization routine
                     should unfix any state variables fixed during
                     initialization (default=False).
                     - True - states varaibles are not unfixed, and
                             a dict of returned containing flags for
                             which states were fixed during
                             initialization.
                    - False - state variables are unfixed after
                             initialization by calling the
                             relase_state method

        Returns:
            If hold_states is True, returns a dict containing flags for
            which states were fixed during initialization.
        '''
        # For now, there are no ocnstraints in the property package, so only
        # fix state variables if required
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="properties")

        init_log.info('Initialization Complete.')

        if hold_state is True:
            flags = fix_state_vars(blk, state_args)
            return flags
        else:
            return

    def release_state(blk, flags, outlvl=idaeslog.NOTSET):
        '''
        Method to release state variables fixed during initialization.

        Keyword Arguments:
            flags : dict containing information of which state variables
                    were fixed during initialization, and should now be
                    unfixed. This dict is returned by initialize if
                    hold_state=True.
            outlvl : sets output level of of logging
        '''
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="properties")

        if flags is None:
            return

        # Unfix state variables
        revert_state_vars(blk, flags)
        init_log.info('State Released.')

    def calculate_state(self, var_args=None, hold_state=False, outlvl=idaeslog.NOTSET,
                            solver=None, optarg=None):
        """
        Solves state blocks given a set of variables and their values. These variables can
        be state variables or properties. This method is typically used before
        initialization to solve for state variables because non-state variables (i.e. properties)
        cannot be fixed in initialization routines.
        Keyword Arguments:
            var_args : dictionary with variables and their values, they can be state variables or properties
                       {(VAR_NAME, INDEX): VALUE}
            hold_state : flag indicating whether all of the state variables should be fixed after calculate state.
                         True - State variables will be fixed.
                         False - State variables will remain unfixed, unless already fixed.
            outlvl : idaes logger object that sets output level of solve call (default=idaeslog.NOTSET)
            solver : solver name string if None is provided the default solver
                     for IDAES will be used (default = None)
            optarg : solver options dictionary object (default={})
        Returns:
            results object from state block solve
        """
        # Get logger
        solve_log = idaeslog.getSolveLogger(self.name, level=outlvl, tag="properties")

        # Initialize at current state values (not user provided)
        self.initialize(solver=solver, optarg=optarg, outlvl=outlvl)

        # Set solver and options
        opt = get_solver(solver, optarg)

        # Fix variables and check degrees of freedom
        flags = {}  # dictionary noting which variables were fixed and their previous state
        for k in self.keys():
            sb = self[k]
            for (v_name, ind), val in var_args.items():
                var = getattr(sb, v_name)
                if iscale.get_scaling_factor(var[ind]) is None:
                    _log.warning(
                            "While using the calculate_state method on {sb_name}, variable {v_name} "
                            "was provided as an argument in var_args, but it does not have a scaling "
                            "factor. This suggests that the calculate_scaling_factor method has not been "
                            "used or the variable was created on demand after the scaling factors were "
                            "calculated. It is recommended to touch all relevant variables (i.e. call "
                            "them or set an initial value) before using the calculate_scaling_factor "
                            "method.".format(v_name=v_name, sb_name=sb.name))
                if var[ind].is_fixed():
                    flags[(k, v_name, ind)] = True
                    if value(var[ind]) != val:
                        raise ConfigurationError(
                            "While using the calculate_state method on {sb_name}, {v_name} was "
                            "fixed to a value {val}, but it was already fixed to value {val_2}. "
                            "Unfix the variable before calling the calculate_state "
                            "method or update var_args."
                            "".format(sb_name=sb.name, v_name=var.name, val=val, val_2=value(var[ind])))
                else:
                    flags[(k, v_name, ind)] = False
                    var[ind].fix(val)

            if degrees_of_freedom(sb) != 0:
                raise RuntimeError("While using the calculate_state method on {sb_name}, the degrees "
                                   "of freedom were {dof}, but 0 is required. Check var_args and ensure "
                                   "the correct fixed variables are provided."
                                   "".format(sb_name=sb.name, dof=degrees_of_freedom(sb)))

        # Solve
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = solve_indexed_blocks(opt, [self], tee=slc.tee)
            solve_log.info_high("Calculate state: {}.".format(idaeslog.condition(results)))

        if not check_optimal_termination(results):
            _log.warning("While using the calculate_state method on {sb_name}, the solver failed "
                         "to converge to an optimal solution. This suggests that the user provided "
                         "infeasible inputs, or that the model is poorly scaled, poorly initialized, "
                         "or degenerate.")

        # unfix all variables fixed with var_args
        for (k, v_name, ind), previously_fixed in flags.items():
            if not previously_fixed:
                var = getattr(self[k], v_name)
                var[ind].unfix()

        # fix state variables if hold_state
        if hold_state:
            fix_state_vars(self)

        return results

@declare_process_block_class("WaterStateBlock",
                             block_class=_WaterStateBlock)
class WaterStateBlockData(StateBlockData):
    """
    General purpose StateBlock for Zero-Order unit models.
    """

    def build(self):
        super().build()

        # Create state variables
        self.flow_mass_comp = Var(self.component_list,
                                  initialize=1,
                                  domain=PositiveReals,
                                  doc='Mass flowrate of each component',
                                  units=pyunits.kg/pyunits.s)

    # -------------------------------------------------------------------------
    # Other properties
    def _conc_mass_comp(self):
        def rule_cmc(blk, j):
            return (blk.flow_mass_comp[j] /
                    sum(self.flow_mass_comp[k] for k in self.component_list) *
                    blk.dens_mass)
        self.conc_mass_comp = Expression(self.component_list,
                                         rule=rule_cmc)

    def _dens_mass(self):
        self.dens_mass = Param(initialize=self.params.dens_mass_default,
                               units=pyunits.kg/pyunits.m**3,
                               mutable=True,
                               doc="Mass density of flow")

    def _flow_vol(self):
        self.flow_vol = Expression(
            expr=sum(self.flow_mass_comp[j] for j in self.component_list) /
            self.dens_mass)

    def _visc_d(self):
        self.visc_d = Param(initialize=self.params.visc_d_default,
                            units=pyunits.kg/pyunits.m/pyunits.s,
                            mutable=True,
                            doc="Dynamic viscosity of solution")

    def get_material_flow_terms(blk, p, j):
        return blk.flow_mass_comp[j]

    def get_enthalpy_flow_terms(blk, p):
        raise NotImplementedError

    def get_material_density_terms(blk, p, j):
        return blk.conc_mass_comp[j]

    def get_energy_density_terms(blk, p):
        raise NotImplementedError

    def default_material_balance_type(self):
        return MaterialBalanceType.componentTotal

    def default_energy_balance_type(self):
        return EnergyBalanceType.none

    def define_state_vars(blk):
        return {"flow_mass_comp": blk.flow_mass_comp}

    def define_display_vars(blk):
        return {"Volumetric Flowrate": blk.flow_vol,
                "Mass Concentration": blk.conc_mass_comp}

    def get_material_flow_basis(blk):
        return MaterialFlowBasis.mass

    def calculate_scaling_factors(self):
        # Get default scale factors and do calculations from base classes
        super().calculate_scaling_factors()

        d_sf_Q = self.params.default_scaling_factor["flow_vol"]
        d_sf_c = self.params.default_scaling_factor["conc_mass_comp"]

        for j, v in self.flow_mass_comp.items():
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, d_sf_Q*d_sf_c)

        if self.is_property_constructed("flow_vol"):
            if iscale.get_scaling_factor(self.flow_vol) is None:
                iscale.set_scaling_factor(self.flow_vol, d_sf_Q)

        if self.is_property_constructed("conc_mass_comp"):
            for j, v in self.conc_mass_comp.items():
                sf_c = iscale.get_scaling_factor(self.conc_mass_comp[j])
                if sf_c is None:
                    try:
                        sf_c = self.params.default_scaling_factor[
                            ("conc_mass_comp", j)]
                    except KeyError:
                        sf_c = d_sf_c
                    iscale.set_scaling_factor(self.conc_mass_comp[j], sf_c)
