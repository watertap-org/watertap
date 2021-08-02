##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Initial property package for H2O-NaCl-CaSO4 system
"""

# Import Python libraries
import idaes.logger as idaeslog

# Import Pyomo libraries
from pyomo.environ import Constraint, Reals, NonNegativeReals, \
    Var, Suffix
from pyomo.environ import units as pyunits
from pyomo.opt import SolverFactory

# Import IDAES cores
from idaes.core import (declare_process_block_class,
                        MaterialFlowBasis,
                        PhysicalParameterBlock,
                        StateBlockData,
                        StateBlock,
                        MaterialBalanceType,
                        EnergyBalanceType)
from idaes.core.components import Component, Solute, Solvent
from idaes.core.phases import LiquidPhase
from idaes.core.util.initialization import (fix_state_vars,
                                            revert_state_vars,
                                            solve_indexed_blocks)
from idaes.core.util.testing import get_default_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.exceptions import PropertyPackageError
import idaes.core.util.scaling as iscale

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("PropParameterBlock")
class PropParameterData(PhysicalParameterBlock):
    CONFIG = PhysicalParameterBlock.CONFIG()

    def build(self):
        '''
        Callable method for Block construction.
        '''
        super(PropParameterData, self).build()

        self._state_block_class = PropStateBlock

        # components
        self.H2O = Solvent()
        self.NaCl = Solute()
        self.CaSO4 = Solute()

        # phases
        self.Liq = LiquidPhase()

        # ---default scaling---
        self.set_default_scaling('temperature', 1e-2)
        self.set_default_scaling('pressure', 1e-6)

    @classmethod
    def define_metadata(cls, obj):
        """Define properties supported and units."""
        obj.add_properties(
            {'flow_mass_phase_comp': {'method': None},
             'temperature': {'method': None},
             'pressure': {'method': None},
             'mass_frac_phase_comp': {'method': '_mass_frac_phase_comp'}
            })

        obj.add_default_units({'time': pyunits.s,
                               'length': pyunits.m,
                               'mass': pyunits.kg,
                               'amount': pyunits.mol,
                               'temperature': pyunits.K})


class _PropStateBlock(StateBlock):
    """
    This Class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """

    def initialize(self, state_args={}, state_vars_fixed=False,
                   hold_state=False, outlvl=idaeslog.NOTSET,
                   solver=None, optarg={}):
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
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default={})
            state_vars_fixed: Flag to denote if state vars have already been
                              fixed.
                              - True - states have already been fixed by the
                                       control volume 1D. Control volume 0D
                                       does not fix the state vars, so will
                                       be False if this state block is used
                                       with 0D blocks.
                             - False - states have not been fixed. The state
                                       block will deal with fixing/unfixing.
            solver : Solver object to use during initialization if None is provided
                     it will use the default solver for IDAES (default = None)
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
        Returns:
            If hold_states is True, returns a dict containing flags for
            which states were fixed during initialization.
        """
        # Get loggers
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="properties")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="properties")

        # Set solver and options
        # TODO: clean up once IDAES new API for initialize solvers is released
        if isinstance(solver, str):
            opt = SolverFactory(solver)
            opt.options = optarg
        else:
            if solver is None:
                opt = get_default_solver()
            else:
                opt = solver
                opt.options = optarg

        # Fix state variables
        flags = fix_state_vars(self, state_args)
        # Check when the state vars are fixed already result in dof 0
        for k in self.keys():
            dof = degrees_of_freedom(self[k])
            if dof != 0:
                raise PropertyPackageError("State vars fixed but degrees of "
                                           "freedom for state block is not "
                                           "zero during initialization.")

        # ---------------------------------------------------------------------
        # Initialize properties
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = solve_indexed_blocks(opt, [self], tee=slc.tee)
        init_log.info("Property initialization: {}."
                      .format(idaeslog.condition(results)))

        # ---------------------------------------------------------------------
        # If input block, return flags, else release state
        if state_vars_fixed is False:
            if hold_state is True:
                return flags
            else:
                self.release_state(flags)

    def release_state(self, flags, outlvl=idaeslog.NOTSET):
        '''
        Method to relase state variables fixed during initialisation.
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
        init_log.info('{} State Released.'.format(self.name))

@declare_process_block_class("PropStateBlock",
                             block_class=_PropStateBlock)
class PropStateBlockData(StateBlockData):
    def build(self):
        """Callable method for Block construction."""
        super(PropStateBlockData, self).build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        # Add state variables
        self.flow_mass_phase_comp = Var(
            self.params.phase_list,
            self.params.component_list,
            initialize=1,
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
    def _mass_frac_phase_comp(self):
        self.mass_frac_phase_comp = Var(
            self.params.phase_list,
            self.params.component_list,
            initialize=0.1,
            bounds=(1e-8, None),
            units=pyunits.dimensionless,
            doc='Mass fraction')

        def rule_mass_frac_phase_comp(b, j):
            return (b.mass_frac_phase_comp['Liq', j] == b.flow_mass_phase_comp['Liq', j] /
                    sum(b.flow_mass_phase_comp['Liq', j]
                        for j in self.params.component_list))
        self.eq_mass_frac_phase_comp = Constraint(self.params.component_list, rule=rule_mass_frac_phase_comp)

    # -----------------------------------------------------------------------------
    # General Methods
    # NOTE: For scaling in the control volume to work properly, these methods must
    # return a pyomo Var or Expression

    def get_material_flow_terms(self, p, j):
        """Create material flow terms for control volume."""
        return self.flow_mass_phase_comp[p, j]

    def get_enthalpy_flow_terms(self, p):
        """Create enthalpy flow terms."""
        return 1

    # TODO: make property package compatible with dynamics
    # def get_material_density_terms(self, p, j):
    #     """Create material density terms."""

    # def get_enthalpy_density_terms(self, p):
    #     """Create enthalpy density terms."""

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
        # for the following variables: flow_mass_phase_comp, pressure,
        # temperature, dens_mass, visc_d, diffus, osm_coeff, and enth_mass

        # these variables should have user input
        if iscale.get_scaling_factor(self.flow_mass_phase_comp['Liq', 'H2O']) is None:
            sf = iscale.get_scaling_factor(self.flow_mass_phase_comp['Liq', 'H2O'], default=1e0, warning=True)
            iscale.set_scaling_factor(self.flow_mass_phase_comp['Liq', 'H2O'], sf)

        if iscale.get_scaling_factor(self.flow_mass_phase_comp['Liq', 'NaCl']) is None:
            sf = iscale.get_scaling_factor(self.flow_mass_phase_comp['Liq', 'NaCl'], default=1e2, warning=True)
            iscale.set_scaling_factor(self.flow_mass_phase_comp['Liq', 'NaCl'], sf)

        if iscale.get_scaling_factor(self.flow_mass_phase_comp['Liq', 'CaSO4']) is None:
            sf = iscale.get_scaling_factor(self.flow_mass_phase_comp['Liq', 'CaSO4'], default=1e4, warning=True)
            iscale.set_scaling_factor(self.flow_mass_phase_comp['Liq', 'CaSO4'], sf)


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

        # # transforming constraints
        # # property relationships with no index, simple constraint
        # v_str_lst_simple = []
        # for v_str in v_str_lst_simple:
        #     if self.is_property_constructed(v_str):
        #         v = getattr(self, v_str)
        #         sf = iscale.get_scaling_factor(v, default=1, warning=True)
        #         c = getattr(self, 'eq_' + v_str)
        #         iscale.constraint_scaling_transform(c, sf)
        #
        # # property relationships with phase index, but simple constraint
        # v_str_lst_phase = []
        # for v_str in v_str_lst_phase:
        #     if self.is_property_constructed(v_str):
        #         v = getattr(self, v_str)
        #         sf = iscale.get_scaling_factor(v['Liq'], default=1, warning=True)
        #         c = getattr(self, 'eq_' + v_str)
        #         iscale.constraint_scaling_transform(c, sf)
        #
        # # property relationship indexed by component
        # v_str_lst_comp = []
        # for v_str in v_str_lst_comp:
        #     if self.is_property_constructed(v_str):
        #         v_comp = getattr(self, v_str)
        #         c_comp = getattr(self, 'eq_' + v_str)
        #         for j, c in c_comp.items():
        #             sf = iscale.get_scaling_factor(v_comp[j], default=1, warning=True)
        #             iscale.constraint_scaling_transform(c, sf)

        # property relationships indexed by component and phase
        v_str_lst_phase_comp = ['mass_frac_phase_comp']
        for v_str in v_str_lst_phase_comp:
            if self.is_property_constructed(v_str):
                v_comp = getattr(self, v_str)
                c_comp = getattr(self, 'eq_' + v_str)
                for j, c in c_comp.items():
                    sf = iscale.get_scaling_factor(v_comp['Liq', j], default=1, warning=True)
                    iscale.constraint_scaling_transform(c, sf)