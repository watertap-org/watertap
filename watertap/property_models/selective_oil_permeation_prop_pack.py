#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

"""
Simple property package for selective oil permeation
"""

from pyomo.environ import Constraint, Var, Param, NonNegativeReals, Suffix
from pyomo.environ import units

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    MaterialFlowBasis,
    PhysicalParameterBlock,
    StateBlockData,
    StateBlock,
    MaterialBalanceType,
    EnergyBalanceType,
)
from idaes.core.base.components import Component
from idaes.core.base.phases import LiquidPhase
from idaes.core.util.initialization import (
    fix_state_vars,
    revert_state_vars,
    solve_indexed_blocks,
)
from idaes.core.util.misc import add_object_reference
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_unfixed_variables,
)
from idaes.core.util.exceptions import PropertyPackageError
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.core.solvers import get_solver
from watertap.core.util.scaling import transform_property_constraints

# Set up logger
_log = idaeslog.getLogger(__name__)


# this first class sets up the property package and creates a central parameter block
@declare_process_block_class("SopParameterBlock")
class SopParameterData(PhysicalParameterBlock):
    CONFIG = PhysicalParameterBlock.CONFIG()

    def build(self):
        """
        Callable method for Block construction.
        """
        super(SopParameterData, self).build()

        self._state_block_class = SopStateBlock

        # phases
        self.Liq = LiquidPhase()

        # components
        self.H2O = Component()
        self.oil = Component()

        # ---default scaling---
        self.set_default_scaling("temperature", 1e-2)
        self.set_default_scaling("pressure", 1e-6)
        self.set_default_scaling("dens_mass_phase", 1e-3, index="Liq")

        # parameters
        self.visc_d_phase_comp = Param(
            self.phase_list,
            self.component_list,
            mutable=True,
            initialize={("Liq", "H2O"): 1e-3, ("Liq", "oil"): 3.56e-3},
            units=units.kg * units.m**-1 * units.s**-1,
            doc="Dynamic viscosity",
        )

        self.dens_mass_phase_comp = Param(
            self.phase_list,
            self.component_list,
            mutable=True,
            initialize={("Liq", "H2O"): 1e3, ("Liq", "oil"): 780},
            units=units.kg * units.m**-3,
            doc="Component mass density",
        )

    @classmethod
    def define_metadata(cls, obj):
        """Define properties supported and units."""
        obj.add_properties(
            {
                "flow_mass_phase_comp": {"method": None},
                "temperature": {"method": None},
                "pressure": {"method": None},
                "visc_d_phase_comp": {"method": "_visc_d_phase_comp"},
                "dens_mass_phase_comp": {"method": "_dens_mass_phase_comp"},
                "mass_frac_phase_comp": {"method": "_mass_frac_phase_comp"},
                "flow_vol_phase_comp": {"method": "_flow_vol_phase_comp"},
                "flow_vol_phase": {"method": "_flow_vol_phase"},
                "flow_vol": {"method": "_flow_vol"},
                "flow_mass_phase": {"method": "_flow_mass_phase"},
                "dens_mass_phase": {"method": "_dens_mass_phase"},
                "vol_frac_phase_comp": {"method": "_vol_frac_phase_comp"},
                "conc_mass_phase_comp": {"method": "_conc_mass_phase_comp"},
            }
        )

        obj.add_default_units(
            {
                "time": units.s,
                "length": units.m,
                "mass": units.kg,
                "amount": units.mol,
                "temperature": units.K,
            }
        )


# this second class contains methods that are applied to state blocks as a whole,
# rather than individual elements of indexed state blocks. Essentially it just includes
# the initialization routine
class _SopStateBlock(StateBlock):
    def initialize(
        self,
        state_args=None,
        state_vars_fixed=False,
        hold_state=False,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
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
                raise PropertyPackageError(
                    "\nWhile initializing {sb_name}, the degrees of freedom "
                    "are {dof}, when zero is required. \nInitialization assumes "
                    "that the state variables should be fixed and that no other "
                    "variables are fixed. \nIf other properties have a "
                    "predetermined value, use the calculate_state method "
                    "before using initialize to determine the values for "
                    "the state variables and avoid fixing the property variables."
                    "".format(sb_name=self.name, dof=dof)
                )

        # ---------------------------------------------------------------------
        skip_solve = True  # skip solve if only state variables are present
        for k in self.keys():
            if number_unfixed_variables(self[k]) != 0:
                skip_solve = False

        if not skip_solve:
            # Initialize properties
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                results = solve_indexed_blocks(opt, [self], tee=slc.tee)
            init_log.info_high(
                "Property initialization: {}.".format(idaeslog.condition(results))
            )

        # ---------------------------------------------------------------------
        # If input block, return flags, else release state
        if state_vars_fixed is False:
            if hold_state is True:
                return flags
            else:
                self.release_state(flags)

    def release_state(self, flags, outlvl=idaeslog.NOTSET):
        """
        Method to release state variables fixed during initialisation.

        Keyword Arguments:
            flags : dict containing information of which state variables
                    were fixed during initialization, and should now be
                    unfixed. This dict is returned by initialize if
                    hold_state=True.
            outlvl : sets output level of of logging
        """
        # Unfix state variables
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="properties")
        revert_state_vars(self, flags)
        init_log.info_high("{} State Released.".format(self.name))


# this third class, provides the instructions to create state blocks
@declare_process_block_class("SopStateBlock", block_class=_SopStateBlock)
class SopStateBlockData(StateBlockData):
    def build(self):
        """Callable method for Block construction."""
        super(SopStateBlockData, self).build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        # first we create the state variables
        self.flow_mass_phase_comp = Var(
            self.params.phase_list,
            self.params.component_list,
            initialize={("Liq", "H2O"): 1, ("Liq", "oil"): 0.01},
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units.kg / units.s,
            doc="Mass flow rate",
        )

        self.temperature = Var(
            initialize=298.15,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units.degK,
            doc="State temperature",
        )

        self.pressure = Var(
            initialize=101325,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units.Pa,
            doc="State pressure",
        )

    # -----------------------------------------------------------------------------
    # Property Methods
    # these property methods build variables and constraints on demand (i.e. only when asked
    # for by a user or unit model)
    def _visc_d_phase_comp(self):
        add_object_reference(self, "visc_d_phase_comp", self.params.visc_d_phase_comp)

    def _dens_mass_phase_comp(self):
        add_object_reference(
            self, "dens_mass_phase_comp", self.params.dens_mass_phase_comp
        )

    def _mass_frac_phase_comp(self):
        self.mass_frac_phase_comp = Var(
            self.params.phase_list,
            self.params.component_list,
            initialize=0.1,
            bounds=(0, 1),
            domain=NonNegativeReals,
            units=units.dimensionless,
            doc="Mass fraction",
        )

        def rule_mass_frac_phase_comp(b, p, j):
            return (
                b.mass_frac_phase_comp[p, j] * b.flow_mass_phase[p]
                == b.flow_mass_phase_comp[p, j]
            )

        self.eq_mass_frac_phase_comp = Constraint(
            self.params.phase_list,
            self.params.component_list,
            rule=rule_mass_frac_phase_comp,
        )

    def _flow_vol_phase_comp(self):
        self.flow_vol_phase_comp = Var(
            self.params.phase_list,
            self.params.component_list,
            initialize=1e-6,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units.m**3 / units.s,
            doc="Component volumetric flowrate",
        )

        def rule_flow_vol_phase_comp(b, p, j):
            return (
                b.flow_vol_phase_comp[p, j] * b.dens_mass_phase_comp[p, j]
                == b.flow_mass_phase_comp[p, j]
            )

        self.eq_flow_vol_phase_comp = Constraint(
            self.params.phase_list,
            self.params.component_list,
            rule=rule_flow_vol_phase_comp,
        )

    def _flow_vol_phase(self):
        self.flow_vol_phase = Var(
            self.params.phase_list,
            initialize=1e-3,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units.m**3 / units.s,
            doc="Phase volumetric flow rate",
        )

        def rule_flow_vol_phase(b, p):
            return b.flow_vol_phase[p] == sum(
                b.flow_vol_phase_comp[p, j] for j in b.params.component_list
            )

        self.eq_flow_vol_phase = Constraint(
            self.params.phase_list, rule=rule_flow_vol_phase
        )

    def _flow_vol(self):
        self.flow_vol = Var(
            initialize=1e-3,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units.m**3 / units.s,
            doc="Total volumetric flow rate",
        )

        def rule_flow_vol(b):
            return b.flow_vol == sum(b.flow_vol_phase[p] for p in b.params.phase_list)

        self.eq_flow_vol = Constraint(rule=rule_flow_vol)

    def _flow_mass_phase(self):
        self.flow_mass_phase = Var(
            self.params.phase_list,
            initialize=1e-3,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units.kg / units.s,
            doc="Mass flow rate",
        )

        def rule_flow_mass_phase(b, p):
            return b.flow_mass_phase[p] == sum(
                b.flow_mass_phase_comp[p, j] for j in b.params.component_list
            )

        self.eq_flow_mass_phase = Constraint(
            self.params.phase_list, rule=rule_flow_mass_phase
        )

    def _dens_mass_phase(self):
        self.dens_mass_phase = Var(
            self.params.phase_list,
            initialize=1e3,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units.kg * units.m**-3,
            doc="Mass density",
        )

        def rule_dens_mass_phase(b, p):
            return b.dens_mass_phase[p] * b.flow_vol_phase[p] == b.flow_mass_phase[p]

        self.eq_dens_mass_phase = Constraint(
            self.params.phase_list, rule=rule_dens_mass_phase
        )

    def _vol_frac_phase_comp(self):
        self.vol_frac_phase_comp = Var(
            self.params.phase_list,
            self.params.component_list,
            initialize=0.1,
            bounds=(0, 1),
            domain=NonNegativeReals,
            units=units.dimensionless,
            doc="Volumetric fraction",
        )

        def rule_vol_frac_phase_comp(b, p, j):
            return (
                b.vol_frac_phase_comp[p, j] * b.flow_vol_phase[p]
                == b.flow_vol_phase_comp[p, j]
            )

        self.eq_vol_frac_phase_comp = Constraint(
            self.params.phase_list,
            self.params.component_list,
            rule=rule_vol_frac_phase_comp,
        )

    def _conc_mass_phase_comp(self):
        self.conc_mass_phase_comp = Var(
            self.params.phase_list,
            self.params.component_list,
            initialize=1,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units.kg / units.m**3,
            doc="Mass concentration",
        )

        def rule_conc_mass_phase_comp(b, p, j):
            return (
                b.conc_mass_phase_comp[p, j]
                == b.mass_frac_phase_comp[p, j] * b.dens_mass_phase[p]
            )

        self.eq_conc_mass_phase_comp = Constraint(
            self.params.phase_list,
            self.params.component_list,
            rule=rule_conc_mass_phase_comp,
        )

    # -----------------------------------------------------------------------------
    # General Methods
    # NOTE: For scaling in the control volume to work properly, these methods must
    # return a pyomo Var or Expression

    def get_material_flow_terms(self, p, j):
        """Create material flow terms for control volume."""
        return self.flow_mass_phase_comp[p, j]

    def default_material_balance_type(self):
        return MaterialBalanceType.componentTotal

    def default_energy_balance_type(self):
        return EnergyBalanceType.none

    def get_material_flow_basis(self):
        return MaterialFlowBasis.mass

    def define_state_vars(self):
        """Define state vars."""
        return {
            "flow_mass_phase_comp": self.flow_mass_phase_comp,
            "temperature": self.temperature,
            "pressure": self.pressure,
        }

    # -----------------------------------------------------------------------------
    # Scaling methods
    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # setting scaling factors for variables

        # default scaling factors have already been set with
        # idaes.core.property_base.calculate_scaling_factors()
        sf_dens_mass_phase = {}
        for p in self.params.phase_list:
            if self.is_property_constructed("dens_mass_phase"):
                sf_dens_mass_phase[p] = iscale.get_scaling_factor(
                    self.dens_mass_phase[p]
                )
            else:
                sf_dens_mass_phase[p] = self.params.get_default_scaling(
                    "dens_mass_phase", index=p
                )

        # these variables should have user input
        for p in self.params.phase_list:
            if (
                iscale.get_scaling_factor(
                    self.flow_mass_phase_comp[p, "H2O"], warning=True
                )
                is None
            ):
                iscale.set_scaling_factor(self.flow_mass_phase_comp[p, "H2O"], 1)

            if (
                iscale.get_scaling_factor(
                    self.flow_mass_phase_comp[p, "oil"], warning=True
                )
                is None
            ):
                iscale.set_scaling_factor(self.flow_mass_phase_comp[p, "oil"], 1e1)

        sf_flow_mass_phase = {}
        for p in self.params.phase_list:
            sf_flow_mass_phase[p] = 1 / (
                1 / iscale.get_scaling_factor(self.flow_mass_phase_comp[p, "H2O"])
                + 1 / iscale.get_scaling_factor(self.flow_mass_phase_comp[p, "oil"])
            )

        # these variables do not typically require user input
        if self.is_property_constructed("flow_mass_phase"):
            for p in self.params.phase_list:
                iscale.set_scaling_factor(
                    self.flow_mass_phase[p], sf_flow_mass_phase[p]
                )

        if self.is_property_constructed("mass_frac_phase_comp"):
            for p in self.params.phase_list:
                for j in self.params.component_list:
                    sf = (
                        iscale.get_scaling_factor(self.flow_mass_phase_comp[p, j])
                        / sf_flow_mass_phase[p]
                    )
                    iscale.set_scaling_factor(self.mass_frac_phase_comp[p, j], sf)

        if self.is_property_constructed("flow_vol_phase_comp"):
            for p in self.params.phase_list:
                for j in self.params.component_list:
                    sf = (
                        iscale.get_scaling_factor(self.flow_mass_phase_comp[p, j])
                        / sf_dens_mass_phase[p]
                    )
                    iscale.set_scaling_factor(self.flow_vol_phase_comp[p, j], sf)

        if self.is_property_constructed("flow_vol_phase"):
            for p in self.params.phase_list:
                sf = sf_flow_mass_phase[p] / sf_dens_mass_phase[p]
                iscale.set_scaling_factor(self.flow_vol_phase[p], sf)

        if self.is_property_constructed("flow_vol"):
            sf = sum(
                sf_flow_mass_phase[p] / sf_dens_mass_phase[p]
                for p in self.params.phase_list
            )
            iscale.set_scaling_factor(self.flow_vol, sf)

        if self.is_property_constructed("vol_frac_phase_comp"):
            for p in self.params.phase_list:
                for j in self.params.component_list:
                    sf = iscale.get_scaling_factor(
                        self.flow_vol_phase_comp[p, j]
                    ) / iscale.get_scaling_factor(self.flow_vol_phase[p])
                    iscale.set_scaling_factor(self.vol_frac_phase_comp[p, j], sf)

        if self.is_property_constructed("conc_mass_phase_comp"):
            for p in self.params.phase_list:
                for j in self.params.component_list:
                    sf = (
                        iscale.get_scaling_factor(self.mass_frac_phase_comp[p, j])
                        * sf_dens_mass_phase[p]
                    )
                    iscale.set_scaling_factor(self.conc_mass_phase_comp[p, j], sf)

        # transform property constraints
        transform_property_constraints(self)
