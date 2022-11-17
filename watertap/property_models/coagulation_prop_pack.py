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

from pyomo.environ import (
    Constraint,
    Var,
    Param,
    Expression,
    Reals,
    NonNegativeReals,
    Suffix,
    value,
    exp,
    assert_optimal_termination,
    check_optimal_termination,
)

from pyomo.environ import units as pyunits

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
from idaes.core.base.phases import LiquidPhase, SolidPhase, PhaseType
from idaes.core.util.initialization import (
    fix_state_vars,
    revert_state_vars,
    solve_indexed_blocks,
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_unfixed_variables,
)
from idaes.core.util.exceptions import (
    PropertyPackageError,
    InitializationError,
    ConfigurationError,
)
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.core.solvers import get_solver

__author__ = "Austin Ladshaw"

# Set up logger
_log = idaeslog.getLogger(__name__)


# Forward declaration of 'CoagulationParameterData'
@declare_process_block_class("CoagulationParameterBlock")
class CoagulationParameterData(PhysicalParameterBlock):
    CONFIG = PhysicalParameterBlock.CONFIG()

    def build(self):
        """
        Callable method for Block construction.
        """
        super(CoagulationParameterData, self).build()

        self._state_block_class = CoagulationStateBlock

        # phases
        self.Liq = LiquidPhase()

        # components
        self.H2O = Component()
        self.TSS = Component()
        self.TDS = Component()
        self.Sludge = Component()

        #   heat capacity of liquid
        self.cp = Param(
            mutable=False, initialize=4184, units=pyunits.J / (pyunits.kg * pyunits.K)
        )

        #   reference density of liquid
        self.ref_dens_liq = Param(
            domain=Reals,
            initialize=999.26,
            mutable=True,
            units=pyunits.kg / pyunits.m**3,
            doc="Reference water mass density parameter @ 0 oC and no salts",
        )

        #   change in liquid density with increasing mass fraction of salts/solids
        self.dens_slope = Param(
            domain=Reals,
            initialize=879.04,
            mutable=True,
            units=pyunits.kg / pyunits.m**3,
            doc="Relative increase in liquid density with mass fraction of salts",
        )

        #   adjustment parameters for density change with temperature
        #   Density calculation as a function of temperature and pressure
        #   --------------------------------------------------------------
        #   Engineering Toolbox. Water - Density, Specific Weight, and
        #   Thermal Expansion Coefficients. (2003) https://www.engineeringtoolbox.com/
        #   water-density-specific-weight-d_595.html [Accessed 02-01-2022]
        self.dens_param_A = Param(
            domain=Reals,
            initialize=-2.9335e-6,
            mutable=True,
            units=pyunits.K**-2,
            doc="Density correction parameter A for temperature variation",
        )

        self.dens_param_B = Param(
            domain=Reals,
            initialize=0.001529811,
            mutable=True,
            units=pyunits.K**-1,
            doc="Density correction parameter B for temperature variation",
        )

        self.dens_param_C = Param(
            domain=Reals,
            initialize=0.787973,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Density correction parameter C for temperature variation",
        )

        #   Correction factors for changes in density with changes in pressure
        self.ref_pressure_correction = Param(
            domain=Reals,
            initialize=1.0135,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Density reference correction parameter for changes in pressure",
        )

        self.ref_pressure_slope = Param(
            domain=Reals,
            initialize=4.9582e-10,
            mutable=True,
            units=pyunits.Pa**-1,
            doc="Slope of density change as a function of pressure",
        )

        #   Adjustment for dynamic viscosity as a function of temperature
        #   -------------------------------------------------------------
        #   D.S. Viswananth, G. Natarajan. Data Book on the Viscosity of
        #     Liquids. Hemisphere Publishing Corp. (1989)
        self.mu_A = Param(
            domain=Reals,
            initialize=2.939e-5,
            mutable=True,
            units=pyunits.kg / pyunits.m / pyunits.s,
            doc="Pre-exponential factor for viscosity calculation",
        )

        self.mu_B = Param(
            domain=Reals,
            initialize=507.88,
            mutable=True,
            units=pyunits.K,
            doc="Exponential numerator term for viscosity calculation",
        )

        self.mu_C = Param(
            domain=Reals,
            initialize=149.3,
            mutable=True,
            units=pyunits.K,
            doc="Exponential denominator term for viscosity calculation",
        )

        # ---default scaling---
        self.set_default_scaling("temperature", 1e-2)
        self.set_default_scaling("pressure", 1e-6)

    @classmethod
    def define_metadata(cls, obj):
        """Define properties supported and units."""
        obj.add_properties(
            {
                "flow_mass_phase_comp": {"method": None},
                "temperature": {"method": None},
                "pressure": {"method": None},
                "mass_frac_phase_comp": {"method": "_mass_frac_phase_comp"},
                "dens_mass_phase": {"method": "_dens_mass_phase"},
                "flow_vol_phase": {"method": "_flow_vol_phase"},
                "conc_mass_phase_comp": {"method": "_conc_mass_phase_comp"},
                "visc_d_phase": {"method": "_visc_d_phase"},
                "enth_flow": {"method": "_enth_flow"},
            }
        )

        obj.add_default_units(
            {
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            }
        )


# this second class contains methods that are applied to state blocks as a whole,
# rather than individual elements of indexed state blocks. Essentially it just includes
# the initialization routine
class _CoagulationStateBlock(StateBlock):
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
                if not check_optimal_termination(results):
                    raise InitializationError(
                        "The property package failed to solve during initialization"
                    )
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

    def calculate_state(
        self,
        var_args=None,
        hold_state=False,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
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
        flags = (
            {}
        )  # dictionary noting which variables were fixed and their previous state
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
                        "method.".format(v_name=v_name, sb_name=sb.name)
                    )
                if var[ind].is_fixed():
                    flags[(k, v_name, ind)] = True
                    if value(var[ind]) != val:
                        raise ConfigurationError(
                            "While using the calculate_state method on {sb_name}, {v_name} was "
                            "fixed to a value {val}, but it was already fixed to value {val_2}. "
                            "Unfix the variable before calling the calculate_state "
                            "method or update var_args."
                            "".format(
                                sb_name=sb.name,
                                v_name=var.name,
                                val=val,
                                val_2=value(var[ind]),
                            )
                        )
                else:
                    flags[(k, v_name, ind)] = False
                    var[ind].fix(val)

            if degrees_of_freedom(sb) != 0:
                raise RuntimeError(
                    "While using the calculate_state method on {sb_name}, the degrees "
                    "of freedom were {dof}, but 0 is required. Check var_args and ensure "
                    "the correct fixed variables are provided."
                    "".format(sb_name=sb.name, dof=degrees_of_freedom(sb))
                )

        # Solve
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = solve_indexed_blocks(opt, [self], tee=slc.tee)
            solve_log.info_high(
                "Calculate state: {}.".format(idaeslog.condition(results))
            )

        if not check_optimal_termination(results):
            _log.warning(
                "While using the calculate_state method on {sb_name}, the solver failed "
                "to converge to an optimal solution. This suggests that the user provided "
                "infeasible inputs, or that the model is poorly scaled, poorly initialized, "
                "or degenerate."
            )

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
@declare_process_block_class(
    "CoagulationStateBlock", block_class=_CoagulationStateBlock
)
class CoagulationStateBlockData(StateBlockData):
    def build(self):
        """Callable method for Block construction."""
        super(CoagulationStateBlockData, self).build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        # first we create the state variables
        #   note: the following dictionaries are used for providing initial values
        #   to state variables and on-demand variables
        self.seawater_mass_frac_dict = {
            ("Liq", "TSS"): 0.01,
            ("Liq", "H2O"): 1,
            ("Liq", "TDS"): 0.01,
            ("Liq", "Sludge"): 0.001,
        }

        self.seawater_mass_flow_dict = {
            ("Liq", "TSS"): 0.01,
            ("Liq", "H2O"): 1,
            ("Liq", "TDS"): 0.01,
            ("Liq", "Sludge"): 0.001,
        }

        self.seawater_mass_conc_dict = {
            ("Liq", "TSS"): 0.01 * 1000,
            ("Liq", "H2O"): 1 * 1000,
            ("Liq", "TDS"): 0.01 * 1000,
            ("Liq", "Sludge"): 0.001 * 1000,
        }

        self.flow_mass_phase_comp = Var(
            self.seawater_mass_frac_dict.keys(),
            initialize=self.seawater_mass_flow_dict,
            bounds=(0.0, None),
            domain=NonNegativeReals,
            units=pyunits.kg / pyunits.s,
            doc="Mass flow rate",
        )

        self.temperature = Var(
            initialize=298.15,
            bounds=(273.15, 700),
            domain=NonNegativeReals,
            units=pyunits.degK,
            doc="State temperature",
        )

        self.pressure = Var(
            initialize=101325,
            bounds=(0.9e5, 3e7),
            domain=NonNegativeReals,
            units=pyunits.Pa,
            doc="State pressure",
        )

    # -----------------------------------------------------------------------------
    # Property Methods
    # these property methods build variables and constraints on demand (i.e. only when asked
    # for by a user or unit model)
    def _mass_frac_phase_comp(self):
        self.mass_frac_phase_comp = Var(
            self.seawater_mass_frac_dict.keys(),
            initialize=self.seawater_mass_frac_dict,
            bounds=(0.0, 1.001),
            units=pyunits.dimensionless,
            doc="Mass fraction",
        )

        def rule_mass_frac_phase_comp(b, p, j):
            return b.mass_frac_phase_comp[p, j] == b.flow_mass_phase_comp[p, j] / (
                b.flow_mass_phase_comp[p, "H2O"]
                + b.flow_mass_phase_comp[p, "TDS"]
                + b.flow_mass_phase_comp[p, "TSS"]
                + b.flow_mass_phase_comp[p, "Sludge"]
            )

        self.eq_mass_frac_phase_comp = Constraint(
            self.seawater_mass_frac_dict.keys(), rule=rule_mass_frac_phase_comp
        )

    def _dens_mass_phase(self):
        self.dens_mass_phase = Var(
            self.params.phase_list,
            initialize={"Liq": 1000},
            bounds=(5e2, 2e4),
            units=pyunits.kg * pyunits.m**-3,
            doc="Mass density",
        )

        #   Density calculation as a function of temperature and pressure
        #   --------------------------------------------------------------
        #   Engineering Toolbox. Water - Density, Specific Weight, and
        #   Thermal Expansion Coefficients. (2003) https://www.engineeringtoolbox.com/
        #   water-density-specific-weight-d_595.html [Accessed 02-01-2022]
        def rule_dens_mass_phase(b, p):
            return b.dens_mass_phase[p] == (
                b.params.ref_dens_liq
                + b.params.dens_slope * b.mass_frac_phase_comp["Liq", "TDS"]
                + b.params.dens_slope * b.mass_frac_phase_comp["Liq", "TSS"]
                + b.params.dens_slope * b.mass_frac_phase_comp["Liq", "Sludge"]
            ) * (
                b.params.dens_param_A * b.temperature**2
                + b.params.dens_param_B * b.temperature
                + b.params.dens_param_C
            ) * (
                b.params.ref_pressure_correction
                + b.params.ref_pressure_slope * b.pressure
            )

        self.eq_dens_mass_phase = Constraint(
            self.params.phase_list, rule=rule_dens_mass_phase
        )

    def _flow_vol_phase(self):
        self.flow_vol_phase = Var(
            self.params.phase_list,
            initialize={"Liq": 1e-3},
            bounds=(0.0, None),
            units=pyunits.m**3 / pyunits.s,
            doc="Volumetric flow rate",
        )

        def rule_flow_vol_phase(b, p):
            return (
                b.flow_vol_phase[p]
                == (
                    b.flow_mass_phase_comp["Liq", "H2O"]
                    + b.flow_mass_phase_comp["Liq", "TDS"]
                    + b.flow_mass_phase_comp["Liq", "TSS"]
                    + b.flow_mass_phase_comp["Liq", "Sludge"]
                )
                / b.dens_mass_phase[p]
            )

        self.eq_flow_vol_phase = Constraint(
            self.params.phase_list, rule=rule_flow_vol_phase
        )

    def _conc_mass_phase_comp(self):
        self.conc_mass_phase_comp = Var(
            self.seawater_mass_frac_dict.keys(),
            initialize=self.seawater_mass_conc_dict,
            bounds=(0.0, None),
            units=pyunits.kg / pyunits.m**3,
            doc="Mass concentration",
        )

        def rule_conc_mass_phase_comp(b, p, j):
            return (
                b.conc_mass_phase_comp[p, j]
                == b.mass_frac_phase_comp[p, j] * b.dens_mass_phase[p]
            )

        self.eq_conc_mass_phase_comp = Constraint(
            self.seawater_mass_frac_dict.keys(), rule=rule_conc_mass_phase_comp
        )

    def _visc_d_phase(self):
        self.visc_d_phase = Var(
            self.params.phase_list,
            initialize={"Liq": 0.001},
            bounds=(0.0, 0.01),
            units=pyunits.kg / pyunits.m / pyunits.s,
            doc="Dynamic viscosity",
        )

        #   Adjustment for dynamic viscosity as a function of temperature
        #   -------------------------------------------------------------
        #   D.S. Viswananth, G. Natarajan. Data Book on the Viscosity of
        #     Liquids. Hemisphere Publishing Corp. (1989)
        def rule_visc_d_phase(b, p):
            return b.visc_d_phase[p] == (
                b.params.mu_A * exp(b.params.mu_B / (b.temperature - b.params.mu_C))
            )

        self.eq_visc_d_phase = Constraint(
            self.params.phase_list, rule=rule_visc_d_phase
        )

    def _enth_flow(self):
        # enthalpy flow expression for get_enthalpy_flow_terms method
        temperature_ref = 273.15 * pyunits.K

        def rule_enth_flow(b, p):  # enthalpy flow [J/s]
            return (
                b.params.cp
                * sum(b.flow_mass_phase_comp[pair] for pair in b.flow_mass_phase_comp)
                * (b.temperature - temperature_ref)
            )

        self.enth_flow = Expression(self.params.phase_list, rule=rule_enth_flow)

    # -----------------------------------------------------------------------------
    # General Methods
    # NOTE: For scaling in the control volume to work properly, these methods must
    # return a pyomo Var or Expression

    def get_material_flow_terms(self, p, j):
        """Create material flow terms for control volume."""
        return self.flow_mass_phase_comp[p, j]

    def get_enthalpy_flow_terms(self, p):
        """Create enthalpy flow terms."""
        return self.enth_flow[p]

    def default_material_balance_type(self):
        return MaterialBalanceType.componentPhase

    def default_energy_balance_type(self):
        return EnergyBalanceType.enthalpyTotal

    def get_material_flow_basis(b):
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

        # scaling factors for parameters
        if iscale.get_scaling_factor(self.params.cp) is None:
            iscale.set_scaling_factor(self.params.cp, 1e-3)

        # these variables should have user input
        if iscale.get_scaling_factor(self.flow_mass_phase_comp["Liq", "H2O"]) is None:
            sf = iscale.get_scaling_factor(
                self.flow_mass_phase_comp["Liq", "H2O"], default=1, warning=True
            )
            iscale.set_scaling_factor(self.flow_mass_phase_comp["Liq", "H2O"], sf)

        if iscale.get_scaling_factor(self.flow_mass_phase_comp["Liq", "TSS"]) is None:
            sf = iscale.get_scaling_factor(
                self.flow_mass_phase_comp["Liq", "TSS"], default=1e2, warning=True
            )
            iscale.set_scaling_factor(self.flow_mass_phase_comp["Liq", "TSS"], sf)

        if iscale.get_scaling_factor(self.flow_mass_phase_comp["Liq", "TDS"]) is None:
            sf = iscale.get_scaling_factor(
                self.flow_mass_phase_comp["Liq", "TDS"], default=1e2, warning=True
            )
            iscale.set_scaling_factor(self.flow_mass_phase_comp["Liq", "TDS"], sf)

        if (
            iscale.get_scaling_factor(self.flow_mass_phase_comp["Liq", "Sludge"])
            is None
        ):
            sf = iscale.get_scaling_factor(
                self.flow_mass_phase_comp["Liq", "Sludge"], default=1e3, warning=True
            )
            iscale.set_scaling_factor(self.flow_mass_phase_comp["Liq", "Sludge"], sf)

        # these variables do not typically require user input,
        # will not override if the user does provide the scaling factor

        if self.is_property_constructed("mass_frac_phase_comp"):
            # Apply variable scaling
            for pair in self.seawater_mass_frac_dict.keys():
                if iscale.get_scaling_factor(self.mass_frac_phase_comp[pair]) is None:
                    if pair[1] == "H2O":
                        iscale.set_scaling_factor(self.mass_frac_phase_comp[pair], 1)
                    else:
                        sf = iscale.get_scaling_factor(
                            self.flow_mass_phase_comp[pair]
                        ) / iscale.get_scaling_factor(
                            self.flow_mass_phase_comp["Liq", "H2O"]
                        )
                        iscale.set_scaling_factor(self.mass_frac_phase_comp[pair], sf)

            # Apply constraint scaling
            for pair in self.seawater_mass_frac_dict.keys():
                sf = iscale.get_scaling_factor(
                    self.mass_frac_phase_comp[pair], default=1, warning=True
                )
                iscale.constraint_scaling_transform(
                    self.eq_mass_frac_phase_comp[pair], sf
                )

        if self.is_property_constructed("dens_mass_phase"):
            if iscale.get_scaling_factor(self.dens_mass_phase) is None:
                iscale.set_scaling_factor(self.dens_mass_phase, 1e-3)

            # transforming constraints
            sf = iscale.get_scaling_factor(self.dens_mass_phase)
            iscale.constraint_scaling_transform(self.eq_dens_mass_phase["Liq"], sf)

        if self.is_property_constructed("flow_vol_phase"):
            sf = iscale.get_scaling_factor(
                self.flow_mass_phase_comp["Liq", "H2O"]
            ) / iscale.get_scaling_factor(self.dens_mass_phase["Liq"])
            iscale.set_scaling_factor(self.flow_vol_phase["Liq"], sf)
            # transforming constraints
            iscale.constraint_scaling_transform(self.eq_flow_vol_phase["Liq"], sf)

        if self.is_property_constructed("conc_mass_phase_comp"):
            # Apply variable scaling
            for pair in self.seawater_mass_frac_dict.keys():
                if iscale.get_scaling_factor(self.conc_mass_phase_comp[pair]) is None:
                    if pair[0] == "Liq":
                        sf = iscale.get_scaling_factor(
                            self.mass_frac_phase_comp[pair]
                        ) * iscale.get_scaling_factor(self.dens_mass_phase["Liq"])
                        iscale.set_scaling_factor(self.conc_mass_phase_comp[pair], sf)
                    else:
                        raise PropertyPackageError(
                            "Unsupported phase for CoagulationParameterData property package"
                        )

            # Apply constraint scaling
            for pair in self.seawater_mass_frac_dict.keys():
                sf = iscale.get_scaling_factor(
                    self.conc_mass_phase_comp[pair], default=1, warning=True
                )
                iscale.constraint_scaling_transform(
                    self.eq_conc_mass_phase_comp[pair], sf
                )

        if self.is_property_constructed("visc_d_phase"):
            if iscale.get_scaling_factor(self.visc_d_phase) is None:
                iscale.set_scaling_factor(self.visc_d_phase, 1e3)

            # transforming constraints
            sf = iscale.get_scaling_factor(self.visc_d_phase)
            iscale.constraint_scaling_transform(self.eq_visc_d_phase["Liq"], sf)

        if self.is_property_constructed("enth_flow"):
            if iscale.get_scaling_factor(self.enth_flow) is None:
                sf = (
                    iscale.get_scaling_factor(self.params.cp)
                    * iscale.get_scaling_factor(self.flow_mass_phase_comp["Liq", "H2O"])
                    * 1e-1
                )
                iscale.set_scaling_factor(self.enth_flow, sf)

    # ------------------------------------------------------------------
    # # TODO: Create functions for easier user access to setup intial values
    #       for all state vars and properties. This may be useful for improving
    #       how scaling factors are determined.
