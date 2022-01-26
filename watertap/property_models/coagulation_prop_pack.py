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
Initial property package for water treatment via a Coagulation-Flocculation process
"""
from pyomo.environ import (Constraint,
                           Var,
                           Param,
                           Expression,
                           Reals,
                           NonNegativeReals,
                           Suffix)
from pyomo.environ import units as pyunits

# Import IDAES cores
from idaes.core import (declare_process_block_class,
                        MaterialFlowBasis,
                        PhysicalParameterBlock,
                        StateBlockData,
                        StateBlock,
                        MaterialBalanceType,
                        EnergyBalanceType)
from idaes.core.components import Component
from idaes.core.phases import LiquidPhase, SolidPhase, PhaseType
from idaes.core.util.initialization import fix_state_vars, revert_state_vars, solve_indexed_blocks
from idaes.core.util.model_statistics import degrees_of_freedom, number_unfixed_variables
from idaes.core.util.exceptions import PropertyPackageError
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.core.util import get_solver

__author__ = "Austin Ladshaw"

# Set up logger
_log = idaeslog.getLogger(__name__)


# Forward declaration of 'PropParameterBlock' from 'CoagulationParameterData'
@declare_process_block_class("PropParameterBlock")
class CoagulationParameterData(PhysicalParameterBlock):
    CONFIG = PhysicalParameterBlock.CONFIG()

    def build(self):
        """
        Callable method for Block construction.
        """
        super(CoagulationParameterData, self).build()

        self._state_block_class = PropStateBlock

        # phases
        self.Liq = LiquidPhase()

        # components
        self.H2O = Component(default={"valid_phase_types": PhaseType.liquidPhase})
        self.TSS = Component(default={"valid_phase_types": PhaseType.liquidPhase})

        # parameters
        self.cp = Param(mutable=False,
                        initialize=4.2e3,
                        units=pyunits.J / (pyunits.kg * pyunits.K))

        dens_mass_param_dict = {('Liq','H2O'): 1000,
                                ('Liq','TSS'): 1400}
        self.dens_mass_param = Param(
            dens_mass_param_dict.keys(),
            domain=Reals,
            initialize=dens_mass_param_dict,
            mutable=True,
            units=pyunits.kg / pyunits.m ** 3,
            doc='Mass density parameters')


        # ---default scaling---
        self.set_default_scaling('temperature', 1e-2)
        self.set_default_scaling('pressure', 1e-6)
        self.set_default_scaling('dens_mass_phase', 1e-3)

    @classmethod
    def define_metadata(cls, obj):
        """Define properties supported and units."""
        obj.add_properties(
            {'flow_mass_phase_comp': {'method': None},
             'temperature': {'method': None},
             'pressure': {'method': None},
             'mass_frac_phase_comp': {'method': '_mass_frac_phase_comp'},
             'dens_mass_phase': {'method': '_dens_mass_phase'},
             'flow_vol_phase': {'method': '_flow_vol_phase'},
             'conc_mass_phase_comp': {'method': '_conc_mass_phase_comp'},
             'enth_flow': {'method': '_enth_flow'},
            })

        obj.add_default_units({'time': pyunits.s,
                               'length': pyunits.m,
                               'mass': pyunits.kg,
                               'amount': pyunits.mol,
                               'temperature': pyunits.K})


# this second class contains methods that are applied to state blocks as a whole,
# rather than individual elements of indexed state blocks. Essentially it just includes
# the initialization routine
class _CoagulationStateBlock(StateBlock):
    def initialize(self, state_args=None, state_vars_fixed=False,
                   hold_state=False, outlvl=idaeslog.NOTSET,
                   solver=None, optarg=None):
        """
        Initialization routine for property package.
        Keyword Arguments:
            state_args : Dictionary with initial guesses for the state vars
                         chosen. Note that if this method is triggered
                         through the control volume, and if initial guesses
                         were not provided at the unit model level, the
                         control volume passes the inlet values as initial
                         guess.The keys for the state_args dictionary are:

                         flow_mass_phase_comp : value at which to initialize
                                               phase component flows
                         pressure : value at which to initialize pressure
                         temperature : value at which to initialize temperature

            state_vars_fixed: Flag to denote if state vars have already been
                              fixed.
                              - True - states have already been fixed by the
                                       control volume 1D. Control volume 0D
                                       does not fix the state vars, so will
                                       be False if this state block is used
                                       with 0D blocks.
                             - False - states have not been fixed. The state
                                       block will deal with fixing/unfixing.
            hold_state : flag indicating whether the initialization routine
                         should unfix any state variables fixed during
                         initialization (default=False).
                         - True - states variables are not unfixed, and
                                 a dict of returned containing flags for
                                 which states were fixed during
                                 initialization.
                        - False - state variables are unfixed after
                                 initialization by calling the
                                 release_state method
            outlvl : sets output level of initialization routine (default=idaeslog.NOTSET)
            solver : Solver object to use during initialization if None is provided
                     it will use the default solver for IDAES (default = None)
            optarg : solver options dictionary object (default=None)
        Returns:
            If hold_states is True, returns a dict containing flags for
            which states were fixed during initialization.
        """
        # Get loggers
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="properties")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="properties")

        # Set solver and options
        opt = get_solver(solver, optarg)

        # Fix state variables
        flags = fix_state_vars(self, state_args)
        # Check when the state vars are fixed already result in dof 0
        for k in self.keys():
            dof = degrees_of_freedom(self[k])
            if dof != 0:
                raise PropertyPackageError("\nWhile initializing {sb_name}, the degrees of freedom "
                                           "are {dof}, when zero is required. \nInitialization assumes "
                                           "that the state variables should be fixed and that no other "
                                           "variables are fixed. \nIf other properties have a "
                                           "predetermined value, use the calculate_state method "
                                           "before using initialize to determine the values for "
                                           "the state variables and avoid fixing the property variables."
                                           "".format(sb_name=self.name, dof=dof))

        # ---------------------------------------------------------------------
        skip_solve = True  # skip solve if only state variables are present
        for k in self.keys():
            if number_unfixed_variables(self[k]) != 0:
                skip_solve = False

        if not skip_solve:
            # Initialize properties
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                results = solve_indexed_blocks(opt, [self], tee=slc.tee)
            init_log.info_high("Property initialization: {}.".format(idaeslog.condition(results)))

        # ---------------------------------------------------------------------
        # If input block, return flags, else release state
        if state_vars_fixed is False:
            if hold_state is True:
                return flags
            else:
                self.release_state(flags)

    def release_state(self, flags, outlvl=idaeslog.NOTSET):
        '''
        Method to release state variables fixed during initialisation.

        Keyword Arguments:
            flags : dict containing information of which state variables
                    were fixed during initialization, and should now be
                    unfixed. This dict is returned by initialize if
                    hold_state=True.
            outlvl : sets output level of of logging
        '''
        # Unfix state variables
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="properties")
        revert_state_vars(self, flags)
        init_log.info_high('{} State Released.'.format(self.name))

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

# this third class, provides the instructions to create state blocks
@declare_process_block_class("PropStateBlock",
                             block_class=_CoagulationStateBlock)
class CoagulationStateBlockData(StateBlockData):
    def build(self):
        """Callable method for Block construction."""
        super(CoagulationStateBlockData, self).build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        # first we create the state variables
        # the following dictionary is used for providing initial values in line 229
        seawater_mass_frac_dict = {('Liq', 'TSS'): 50e-6,
                                   ('Liq', 'H2O'): 0.965}

        self.flow_mass_phase_comp = Var(
            self.params.phase_list,
            self.params.component_list,
            initialize=seawater_mass_frac_dict,
            bounds=(1e-8, 100),
            domain=NonNegativeReals,
            units=pyunits.kg/pyunits.s,
            doc='Mass flow rate')

        self.temperature = Var(
            initialize=298.15,
            bounds=(273.15, 1000),
            domain=NonNegativeReals,
            units=pyunits.degK,
            doc='State temperature')

        self.pressure = Var(
            initialize=101325,
            bounds=(1e5, 5e7),
            domain=NonNegativeReals,
            units=pyunits.Pa,
            doc='State pressure')

    # -----------------------------------------------------------------------------
    # Property Methods
    # these property methods build variables and constraints on demand (i.e. only when asked
    # for by a user or unit model)
    def _mass_frac_phase_comp(self):
        self.mass_frac_phase_comp = Var(
            self.params.phase_list,
            self.params.component_list,
            initialize=0.1,
            bounds=(1e-8, None),
            units=pyunits.dimensionless,
            doc='Mass fraction')

        def rule_mass_frac_phase_comp(b, p, j):
            return (b.mass_frac_phase_comp[p, j] == b.flow_mass_phase_comp[p, j] /
                    sum(b.flow_mass_phase_comp[p, i]
                        for i in self.params.component_list))
        self.eq_mass_frac_phase_comp = Constraint(self.params.phase_list, self.params.component_list, rule=rule_mass_frac_phase_comp)

    def _dens_mass_phase(self):
        self.dens_mass_phase = Var(
            self.params.phase_list,
            initialize=1e3,
            bounds=(5e2, 2e3),
            units=pyunits.kg * pyunits.m**-3,
            doc="Mass density")

        def rule_dens_mass_phase(b, p):
            return (b.dens_mass_phase[p] ==
                    sum(b.params.dens_mass_param[(p,i)] * b.mass_frac_phase_comp[p, i]
                        for i in self.params.component_list))
        self.eq_dens_mass_phase = Constraint(self.params.phase_list, rule=rule_dens_mass_phase)


    def _flow_vol_phase(self):
        self.flow_vol_phase = Var(
            self.params.phase_list,
            initialize=1e-3,
            bounds=(1e-8, None),
            units=pyunits.m**3/pyunits.s,
            doc="Volumetric flow rate")

        def rule_flow_vol_phase(b, p):
            return (b.flow_vol_phase[p]
                    == sum(b.flow_mass_phase_comp[p, j] for j in b.params.component_list)
                    / b.dens_mass_phase[p])
        self.eq_flow_vol_phase = Constraint(self.params.phase_list, rule=rule_flow_vol_phase)

    def _conc_mass_phase_comp(self):
        self.conc_mass_phase_comp = Var(
            self.params.phase_list,
            self.params.component_list,
            initialize=1,
            bounds=(1e-8, None),
            units=pyunits.kg/pyunits.m**3,
            doc="Molar flowrate")

        def rule_conc_mass_phase_comp(b, p, j):
            return (b.conc_mass_phase_comp[p, j] ==
                    b.mass_frac_phase_comp[p, j] * b.dens_mass_phase[p])
        self.eq_conc_mass_phase_comp = Constraint(self.params.phase_list, self.params.component_list, rule=rule_conc_mass_phase_comp)

    def _enth_flow(self):
        # enthalpy flow expression for get_enthalpy_flow_terms method
        temperature_ref = 273.15 * pyunits.K

        def rule_enth_flow(b):  # enthalpy flow [J/s]
            return (b.params.cp
                    * sum(b.flow_mass_phase_comp[p, j] for (p, j) in (b.params.phase_list, b.params.component_list))
                    * (b.temperature - temperature_ref))
        self.enth_flow = Expression(rule=rule_enth_flow)

    # -----------------------------------------------------------------------------
    # General Methods
    # NOTE: For scaling in the control volume to work properly, these methods must
    # return a pyomo Var or Expression

    def get_material_flow_terms(self, p, j):
        """Create material flow terms for control volume."""
        return self.flow_mass_phase_comp[p, j]

    def get_enthalpy_flow_terms(self, p):
        """Create enthalpy flow terms."""
        return self.enth_flow

    def default_material_balance_type(self):
        return MaterialBalanceType.componentTotal

    def default_energy_balance_type(self):
        return EnergyBalanceType.enthalpyTotal

    def get_material_flow_basis(b):
        return MaterialFlowBasis.mass

    def define_state_vars(self):
        """Define state vars."""
        return {"flow_mass_phase_comp": self.flow_mass_phase_comp,
                "temperature": self.temperature,
                "pressure": self.pressure}

    # -----------------------------------------------------------------------------
    # Scaling methods
    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # setting scaling factors for variables

        # default scaling factors have already been set with
        # idaes.core.property_base.calculate_scaling_factors()

        # scaling factors for parameters
        if iscale.get_scaling_factor(self.params.cp) is None:
            iscale.set_scaling_factor(self.params.cp, 1e-3)

        # these variables should have user input
        if iscale.get_scaling_factor(self.flow_mass_phase_comp['Liq', 'H2O']) is None:
            sf = iscale.get_scaling_factor(self.flow_mass_phase_comp['Liq', 'H2O'], default=1, warning=True)
            iscale.set_scaling_factor(self.flow_mass_phase_comp['Liq', 'H2O'], sf)

        if iscale.get_scaling_factor(self.flow_mass_phase_comp['Liq', 'TSS']) is None:
            sf = iscale.get_scaling_factor(self.flow_mass_phase_comp['Liq', 'TSS'], default=1e5, warning=True)
            iscale.set_scaling_factor(self.flow_mass_phase_comp['Liq', 'TSS'], sf)

        # these variables do not typically require user input,
        # will not override if the user does provide the scaling factor
        if self.is_property_constructed('mass_frac_phase_comp'):
            for j in self.params.component_list:
                if iscale.get_scaling_factor(self.mass_frac_phase_comp['Liq', j]) is None:
                    if j == 'H2O':
                        iscale.set_scaling_factor(self.mass_frac_phase_comp['Liq', j], 1)
                    else:
                        sf = (iscale.get_scaling_factor(self.flow_mass_phase_comp['Liq', j])
                              / iscale.get_scaling_factor(self.flow_mass_phase_comp['Liq', 'H2O']))
                        iscale.set_scaling_factor(self.mass_frac_phase_comp['Liq', j], sf)

        if self.is_property_constructed('flow_vol_phase'):
            sf = (iscale.get_scaling_factor(self.flow_mass_phase_comp['Liq', 'H2O'])
                  / iscale.get_scaling_factor(self.dens_mass_phase['Liq']))
            iscale.set_scaling_factor(self.flow_vol_phase, sf)

        if self.is_property_constructed('conc_mass_phase_comp'):
            for j in self.params.component_list:
                if iscale.get_scaling_factor(self.conc_mass_phase_comp['Liq', j]) is None:
                    sf = (iscale.get_scaling_factor(self.mass_frac_phase_comp['Liq', j])
                          * iscale.get_scaling_factor(self.dens_mass_phase['Liq']))
                    iscale.set_scaling_factor(self.conc_mass_phase_comp['Liq', j], sf)

        if self.is_property_constructed('enth_flow'):
            if iscale.get_scaling_factor(self.enth_flow) is None:
                sf = (iscale.get_scaling_factor(self.params.cp)
                      * iscale.get_scaling_factor(self.flow_mass_phase_comp['Liq', 'H2O'])
                      * 1e-1)  # temperature change on the order of 1e1
                iscale.set_scaling_factor(self.enth_flow, sf)

        if self.is_property_constructed('dens_mass_phase'):
            if iscale.get_scaling_factor(self.dens_mass_phase) is None:
                iscale.set_scaling_factor(self.dens_mass_phase, 1e-3)

        # # TODO:  FIX THIS GARBAGE
        '''
        # transforming constraints
        # property relationships with no index, simple constraint
        v_str_lst_simple = ['dens_mass_phase', 'flow_vol_phase']
        for v_str in v_str_lst_simple:
            if self.is_property_constructed(v_str):
                v = getattr(self, v_str)
                sf = iscale.get_scaling_factor(v, default=1, warning=True)
                print("---hher")
                print(sf)
                print()
                c = getattr(self, 'eq_' + v_str)
                print(c)
                c.pprint()
                #iscale.constraint_scaling_transform(c, sf)

        # property relationships indexed by component and phase
        v_str_lst_phase_comp = ['mass_frac_phase_comp', 'conc_mass_phase_comp']
        for v_str in v_str_lst_phase_comp:
            if self.is_property_constructed(v_str):
                v_comp = getattr(self, v_str)
                c_comp = getattr(self, 'eq_' + v_str)
                for p in self.params.phase_list:
                    for j in self.params.component_list:
                        sf = iscale.get_scaling_factor(v_comp[p, j], default=1, warning=True)
                        print("====NO")
                        c.pprint()
                        print()
                        iscale.constraint_scaling_transform(c[p], sf)
        #end
        '''
