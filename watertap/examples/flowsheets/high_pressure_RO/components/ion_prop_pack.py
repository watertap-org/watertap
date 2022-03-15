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
Simple property package for solution represented with ions
"""

from pyomo.environ import (Constraint,
                           Var,
                           Param,
                           Expression,
                           Reals,
                           NonNegativeReals,
                           Suffix)
from pyomo.environ import units as pyunits
from pyomo.common.config import ConfigValue

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
from idaes.core.util.initialization import fix_state_vars, revert_state_vars, solve_indexed_blocks
from idaes.core.util.model_statistics import degrees_of_freedom, number_unfixed_variables
from idaes.core.util.exceptions import PropertyPackageError
from idaes.core.util.misc import extract_data
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.core.util import get_solver

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("PropParameterBlock")
class PropParameterData(PhysicalParameterBlock):
    CONFIG = PhysicalParameterBlock.CONFIG()
    CONFIG.declare("ion_list", ConfigValue(
        domain=list,
        description="List of ion species names"))
    # CONFIG.declare("mw_data", ConfigValue(
    #     default={},
    #     domain=dict,
    #     description="Dict of component names and molecular weight data"))

    def build(self):
        '''
        Callable method for Block construction.
        '''
        super(PropParameterData, self).build()

        self._state_block_class = PropStateBlock

        # phases
        self.Liq = LiquidPhase()

        # components
        self.H2O = Solvent()

        for j in self.config.ion_list:
            self.add_component(str(j), Solute())

        # # molecular weight
        # self.mw_comp = Param(
        #     self.component_list,
        #     mutable=True,
        #     default=18e-3,
        #     initialize=self.config.mw_data,
        #     units=pyunits.kg / pyunits.mol,
        #     doc="Molecular weight")

        # ---default scaling---
        self.set_default_scaling('flow_vol', 1e3)
        self.set_default_scaling('conc_mass_comp', 1e-1)
        self.set_default_scaling('temperature', 1e-2)
        self.set_default_scaling('pressure', 1e-6)

    @classmethod
    def define_metadata(cls, obj):
        """Define properties supported and units."""
        obj.add_properties(
            {'flow_vol': {'method': None},
             'conc_mass_comp': {'method': None},
             'temperature': {'method': None},
             'pressure': {'method': None}
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
            outlvl : sets output level of initialization routine (default=idaeslog.NOTSET)
            optarg : solver options dictionary object (default=None)
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

@declare_process_block_class("PropStateBlock",
                             block_class=_PropStateBlock)
class PropStateBlockData(StateBlockData):
    def build(self):
        """Callable method for Block construction."""
        super(PropStateBlockData, self).build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        # Add state variables
        self.flow_vol = Var(
            initialize=1e-3,
            bounds=(1e-8, 100),
            domain=NonNegativeReals,
            units=pyunits.m**3/pyunits.s,
            doc='Volumetric flow rate')

        self.conc_mass_comp = Var(
            self.params.component_list,
            initialize=1,
            bounds=(1e-6, 1e6),
            units=pyunits.kg/pyunits.m ** 3,
            doc="Concentration")

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

    def define_state_vars(self):
        """Define state vars."""
        return {"flow_vol": self.flow_vol,
                "conc_mass_comp": self.conc_mass_comp,
                "temperature": self.temperature,
                "pressure": self.pressure}

    # -----------------------------------------------------------------------------
    # Scaling methods
    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # default scaling factors have already been set with
        # idaes.core.property_base.calculate_scaling_factors()
