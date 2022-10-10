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
Initial property package for H2O-NDMA system
"""

# Import Python libraries
import idaes.logger as idaeslog

# Import Pyomo libraries
from pyomo.environ import (
    Constraint,
    Expression,
    Reals,
    NonNegativeReals,
    Var,
    Param,
    Suffix,
    value,
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
from idaes.core.base.components import Component, Solute, Solvent
from idaes.core.base.phases import LiquidPhase
from idaes.core.util.constants import Constants
from idaes.core.util.initialization import (
    fix_state_vars,
    revert_state_vars,
    solve_indexed_blocks,
)
from idaes.core.util.misc import add_object_reference, extract_data
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_unfixed_variables,
)
from idaes.core.util.exceptions import (
    ConfigurationError,
    InitializationError,
    PropertyPackageError,
)
import idaes.core.util.scaling as iscale

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("NDMAParameterBlock")
class NDMAParameterData(PhysicalParameterBlock):
    CONFIG = PhysicalParameterBlock.CONFIG()

    def build(self):
        """
        Callable method for Block construction.
        """
        super(NDMAParameterData, self).build()

        self._state_block_class = NDMAStateBlock

        # components
        self.H2O = Solvent()
        self.NDMA = Solute()

        # phases
        self.Liq = LiquidPhase()

        # heat capacity
        self.cp = Param(
            mutable=False, initialize=4.2e3, units=pyunits.J / (pyunits.kg * pyunits.K)
        )

        # molecular weight
        mw_comp_data = {"H2O": 18.01528e-3, "NDMA": 74.0819e-3}
        self.mw_comp = Param(
            self.component_list,
            mutable=False,
            initialize=extract_data(mw_comp_data),
            units=pyunits.kg / pyunits.mol,
            doc="Molecular weight kg/mol",
        )

        # mass density of water
        self.dens_mass = Var(
            domain=Reals,
            initialize=997,
            units=pyunits.kg / pyunits.m**3,
            doc="Mass density kg/m^3",
        )

        # traditional parameters are the only Vars currently on the block and should be fixed
        for v in self.component_objects(Var):
            v.fix()

        # ---default scaling---
        self.set_default_scaling("temperature", 1e-2)
        self.set_default_scaling("pressure", 1e-6)
        self.set_default_scaling("dens_mass_phase", 1e-3, index="Liq")

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
                "flow_vol": {"method": "_flow_vol"},
                "conc_mass_phase_comp": {"method": "_conc_mass_phase_comp"},
                "flow_mol_phase_comp": {"method": "_flow_mol_phase_comp"},
                "mole_frac_phase_comp": {"method": "_mole_frac_phase_comp"},
                "molality_phase_comp": {"method": "_molality_phase_comp"},
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


class _NDMAStateBlock(StateBlock):
    """
    This Class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """

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

            if not check_optimal_termination(results):
                raise InitializationError(
                    f"{self.name} failed to initialize successfully. Please "
                    f"check the output logs for more information."
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
        # dictionary noting which variables were fixed and their previous state
        flags = {}
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


@declare_process_block_class("NDMAStateBlock", block_class=_NDMAStateBlock)
class NDMAStateBlockData(StateBlockData):
    def build(self):
        """Callable method for Block construction."""
        super(NDMAStateBlockData, self).build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        # Add state variables
        self.flow_mass_phase_comp = Var(
            self.params.phase_list,
            self.params.component_list,
            initialize={("Liq", "H2O"): 0.999999926, ("Liq", "NDMA"): 74e-9},
            bounds=(1e-18, None),
            domain=NonNegativeReals,
            units=pyunits.kg / pyunits.s,
            doc="Mass flow rate",
        )

        self.temperature = Var(
            initialize=298.15,
            bounds=(273.15, 373.15),
            domain=NonNegativeReals,
            units=pyunits.degK,
            doc="State temperature",
        )

        self.pressure = Var(
            initialize=101325,
            bounds=(1e4, 5e7),
            domain=NonNegativeReals,
            units=pyunits.Pa,
            doc="State pressure",
        )

    # -----------------------------------------------------------------------------
    # Property Methods
    def _mass_frac_phase_comp(self):
        self.mass_frac_phase_comp = Var(
            self.params.phase_list,
            self.params.component_list,
            initialize={("Liq", "H2O"): 0.999999926, ("Liq", "NDMA"): 74e-9},
            bounds=(
                1e-18,
                None,
            ),  # upper bound set to None because of stability benefits
            units=pyunits.dimensionless,
            doc="Mass fraction",
        )

        def rule_mass_frac_phase_comp(b, j):
            return b.mass_frac_phase_comp["Liq", j] == b.flow_mass_phase_comp[
                "Liq", j
            ] / sum(
                b.flow_mass_phase_comp["Liq", j] for j in self.params.component_list
            )

        self.eq_mass_frac_phase_comp = Constraint(
            self.params.component_list, rule=rule_mass_frac_phase_comp
        )

    def _dens_mass_phase(self):
        self.dens_mass_phase = Var(
            self.params.phase_list,
            initialize=997,
            bounds=(5e2, 2e3),
            units=pyunits.kg * pyunits.m**-3,
            doc="Mass density",
        )

        def rule_dens_mass_phase(b):
            # constant density, NDMA concentration is always minimal
            return b.dens_mass_phase["Liq"] == b.params.dens_mass

        self.eq_dens_mass_phase = Constraint(rule=rule_dens_mass_phase)

    def _flow_vol_phase(self):
        self.flow_vol_phase = Var(
            self.params.phase_list,
            initialize=1,
            bounds=(1e-18, None),
            units=pyunits.m**3 / pyunits.s,
            doc="Volumetric flow rate",
        )

        def rule_flow_vol_phase(b):
            return (
                b.flow_vol_phase["Liq"]
                == sum(
                    b.flow_mass_phase_comp["Liq", j] for j in self.params.component_list
                )
                / b.dens_mass_phase["Liq"]
            )

        self.eq_flow_vol_phase = Constraint(rule=rule_flow_vol_phase)

    def _flow_vol(self):
        def rule_flow_vol(b):
            return sum(b.flow_vol_phase[p] for p in self.params.phase_list)

        self.flow_vol = Expression(rule=rule_flow_vol)

    def _conc_mass_phase_comp(self):
        self.conc_mass_phase_comp = Var(
            self.params.phase_list,
            self.params.component_list,
            initialize=10,
            bounds=(1e-18, 2e3),
            units=pyunits.kg * pyunits.m**-3,
            doc="Mass concentration",
        )

        def rule_conc_mass_phase_comp(b, j):
            return (
                b.conc_mass_phase_comp["Liq", j]
                == b.dens_mass_phase["Liq"] * b.mass_frac_phase_comp["Liq", j]
            )

        self.eq_conc_mass_phase_comp = Constraint(
            self.params.component_list, rule=rule_conc_mass_phase_comp
        )

    def _flow_mol_phase_comp(self):
        self.flow_mol_phase_comp = Var(
            self.params.phase_list,
            self.params.component_list,
            initialize=100,
            bounds=(1e-18, None),
            units=pyunits.mol / pyunits.s,
            doc="Molar flowrate",
        )

        def rule_flow_mol_phase_comp(b, j):
            return (
                b.flow_mol_phase_comp["Liq", j]
                == b.flow_mass_phase_comp["Liq", j] / b.params.mw_comp[j]
            )

        self.eq_flow_mol_phase_comp = Constraint(
            self.params.component_list, rule=rule_flow_mol_phase_comp
        )

    def _mole_frac_phase_comp(self):
        self.mole_frac_phase_comp = Var(
            self.params.phase_list,
            self.params.component_list,
            initialize=0.1,
            bounds=(1e-18, None),
            units=pyunits.dimensionless,
            doc="Mole fraction",
        )

        def rule_mole_frac_phase_comp(b, j):
            return b.mole_frac_phase_comp["Liq", j] == b.flow_mol_phase_comp[
                "Liq", j
            ] / sum(b.flow_mol_phase_comp["Liq", j] for j in b.params.component_list)

        self.eq_mole_frac_phase_comp = Constraint(
            self.params.component_list, rule=rule_mole_frac_phase_comp
        )

    def _molality_phase_comp(self):
        self.molality_phase_comp = Var(
            self.params.phase_list,
            ["NDMA"],
            initialize=1,
            bounds=(1e-18, 10),
            units=pyunits.mole / pyunits.kg,
            doc="Molality",
        )

        def rule_molality_phase_comp(b, j):
            return (
                self.molality_phase_comp["Liq", j]
                == b.mass_frac_phase_comp["Liq", j]
                / (1 - b.mass_frac_phase_comp["Liq", j])
                / b.params.mw_comp[j]
            )

        self.eq_molality_phase_comp = Constraint(
            ["NDMA"], rule=rule_molality_phase_comp
        )

    # TODO: check enthalpy
    def _enth_flow(self):
        # enthalpy flow expression for get_enthalpy_flow_terms method
        temperature_ref = 273.15 * pyunits.K

        def rule_enth_flow(b):  # enthalpy flow [J/s]
            return (
                b.params.cp
                * sum(b.flow_mass_phase_comp["Liq", j] for j in b.params.component_list)
                * (b.temperature - temperature_ref)
            )

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

    # TODO: make property package compatible with dynamics
    # def get_material_density_terms(self, p, j):
    #     """Create material density terms."""

    # def get_enthalpy_density_terms(self, p):
    #     """Create enthalpy density terms."""

    def default_material_balance_type(self):
        return MaterialBalanceType.componentTotal

    def default_energy_balance_type(self):
        return EnergyBalanceType.enthalpyTotal

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
        # for the following variables: flow_mass_phase_comp, pressure,
        # temperature, dens_mass, and enth_mass

        if iscale.get_scaling_factor(self.params.cp) is None:
            iscale.set_scaling_factor(self.params.cp, 1e-3)

        # these variables should have user input
        if iscale.get_scaling_factor(self.flow_mass_phase_comp["Liq", "H2O"]) is None:
            sf = iscale.get_scaling_factor(
                self.flow_mass_phase_comp["Liq", "H2O"], default=1e0, warning=True
            )
            iscale.set_scaling_factor(self.flow_mass_phase_comp["Liq", "H2O"], sf)

        if iscale.get_scaling_factor(self.flow_mass_phase_comp["Liq", "NDMA"]) is None:
            sf = iscale.get_scaling_factor(
                self.flow_mass_phase_comp["Liq", "NDMA"], default=1e8, warning=True
            )
            iscale.set_scaling_factor(self.flow_mass_phase_comp["Liq", "NDMA"], sf)

        # scaling factors for parameters
        for j, v in self.params.mw_comp.items():
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(self.params.mw_comp, 1e2)

        # these variables do not typically require user input,
        # will not override if the user does provide the scaling factor
        if self.is_property_constructed("mass_frac_phase_comp"):
            for j in self.params.component_list:
                if (
                    iscale.get_scaling_factor(self.mass_frac_phase_comp["Liq", j])
                    is None
                ):
                    if j == "NDMA":
                        sf = iscale.get_scaling_factor(
                            self.flow_mass_phase_comp["Liq", j]
                        ) / iscale.get_scaling_factor(
                            self.flow_mass_phase_comp["Liq", "H2O"]
                        )
                        iscale.set_scaling_factor(
                            self.mass_frac_phase_comp["Liq", j], sf
                        )
                    elif j == "H2O":
                        iscale.set_scaling_factor(
                            self.mass_frac_phase_comp["Liq", j], 100
                        )

        if self.is_property_constructed("flow_vol_phase"):
            sf = iscale.get_scaling_factor(
                self.flow_mass_phase_comp["Liq", "H2O"]
            ) / iscale.get_scaling_factor(self.dens_mass_phase["Liq"])
            iscale.set_scaling_factor(self.flow_vol_phase, sf)

        if self.is_property_constructed("flow_vol"):
            sf = iscale.get_scaling_factor(self.flow_vol_phase)
            iscale.set_scaling_factor(self.flow_vol, sf)

        if self.is_property_constructed("conc_mass_phase_comp"):
            for j in self.params.component_list:
                sf_dens = iscale.get_scaling_factor(self.dens_mass_phase["Liq"])
                if (
                    iscale.get_scaling_factor(self.conc_mass_phase_comp["Liq", j])
                    is None
                ):
                    if j == "H2O":
                        # solvents typically have a mass fraction between 0.5-1
                        iscale.set_scaling_factor(
                            self.conc_mass_phase_comp["Liq", j], sf_dens
                        )
                    elif j == "NDMA":
                        iscale.set_scaling_factor(
                            self.conc_mass_phase_comp["Liq", j],
                            sf_dens
                            * iscale.get_scaling_factor(
                                self.mass_frac_phase_comp["Liq", j]
                            ),
                        )

        if self.is_property_constructed("flow_mol_phase_comp"):
            for j in self.params.component_list:
                if (
                    iscale.get_scaling_factor(self.flow_mol_phase_comp["Liq", j])
                    is None
                ):
                    sf = iscale.get_scaling_factor(self.flow_mass_phase_comp["Liq", j])
                    sf /= iscale.get_scaling_factor(self.params.mw_comp[j])
                    iscale.set_scaling_factor(self.flow_mol_phase_comp["Liq", j], sf)

        if self.is_property_constructed("mole_frac_phase_comp"):
            for j in self.params.component_list:
                if (
                    iscale.get_scaling_factor(self.mole_frac_phase_comp["Liq", j])
                    is None
                ):
                    if j == "NDMA":
                        sf = iscale.get_scaling_factor(
                            self.flow_mol_phase_comp["Liq", j]
                        ) / iscale.get_scaling_factor(
                            self.flow_mol_phase_comp["Liq", "H2O"]
                        )
                        iscale.set_scaling_factor(
                            self.mole_frac_phase_comp["Liq", j], sf
                        )
                    elif j == "H2O":
                        iscale.set_scaling_factor(
                            self.mole_frac_phase_comp["Liq", j], 1
                        )

        if self.is_property_constructed("molality_phase_comp"):
            for j in self.params.component_list:
                if isinstance(getattr(self.params, j), Solute):
                    if (
                        iscale.get_scaling_factor(self.molality_phase_comp["Liq", j])
                        is None
                    ):
                        sf = iscale.get_scaling_factor(
                            self.mass_frac_phase_comp["Liq", j]
                        )
                        sf /= iscale.get_scaling_factor(self.params.mw_comp[j])
                        iscale.set_scaling_factor(
                            self.molality_phase_comp["Liq", j], sf
                        )

        if self.is_property_constructed("enth_flow"):
            sf = (
                iscale.get_scaling_factor(self.params.cp)
                * iscale.get_scaling_factor(self.flow_mass_phase_comp["Liq", "H2O"])
                * 1e-1
            )  # temperature change on the order of 1e1
            iscale.set_scaling_factor(self.enth_flow, sf)

        # property relationships with phase index, but simple constraint
        v_str_lst_simple = ["dens_mass_phase", "flow_vol_phase"]
        for v_str in v_str_lst_simple:
            if self.is_property_constructed(v_str):
                v = getattr(self, v_str)
                sf = iscale.get_scaling_factor(v, default=1, warning=True)
                c = getattr(self, "eq_" + v_str)
                iscale.constraint_scaling_transform(c, sf)

        # property relationships indexed by component and phase
        for v_str in (
            "mass_frac_phase_comp",
            "conc_mass_phase_comp",
            "flow_mol_phase_comp",
            "mole_frac_phase_comp",
            "molality_phase_comp",
        ):
            if self.is_property_constructed(v_str):
                v_comp = self.component(v_str)
                for j, c in self.component("eq_" + v_str).items():
                    sf = iscale.get_scaling_factor(
                        v_comp["Liq", j], default=1, warning=True
                    )
                    iscale.constraint_scaling_transform(c, sf)
