#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
"""
Initial property package for seawater system
"""

# Import Pyomo libraries
from pyomo.environ import (
    check_optimal_termination,
    Constraint,
    Suffix,
    value,
    Var,
)
from pyomo.environ import units as pyunits
from pyomo.contrib.solver.common.util import NoSolutionError

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    StateBlock,
)
from idaes.core.base.components import Solute
from idaes.core.initialization.block_triangularization import (
    BlockTriangularizationInitializer,
)
from idaes.core.util.initialization import (
    fix_state_vars,
    revert_state_vars,
    solve_indexed_blocks,
)
import idaes.logger as idaeslog
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
)
from idaes.core.util.exceptions import (
    ConfigurationError,
    InitializationError,
    PropertyPackageError,
)
import idaes.core.util.scaling as iscale
from idaes.core.util.exceptions import BurntToast

# Import WaterTAP libraries
from watertap.core.solvers import get_solver
from watertap.core.util.scaling import transform_property_constraints

from watertap.property_models.seawater_prop_pack import (
    SeawaterParameterData,
    SeawaterStateBlockData,
)

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("SeawaterVTPCParameterBlock")
class SeawaterVTPCParameterData(SeawaterParameterData):
    """References
    Parameter block for a seawater property package. Two components: 'H2O' and 'TDS', and one phase: 'Liq'.

    This package was developed from the following references:
    - K.G.Nayar, M.H.Sharqawy, L.D.Banchik, and J.H.Lienhard V, "Thermophysical properties of seawater: A review and
    new correlations that include pressure dependence,"Desalination, Vol.390, pp.1 - 24, 2016.
    doi: 10.1016/j.desal.2016.02.024(preprint)
    - Mostafa H.Sharqawy, John H.Lienhard V, and Syed M.Zubair, "Thermophysical properties of seawater: A review of
    existing correlations and data,"Desalination and Water Treatment, Vol.16, pp.354 - 380, April 2010.
    (2017 corrections provided at http://web.mit.edu/seawater)
    Diffusivity for NaCl is being used temporarily based on
    Bartholomew & Mauter (2019) https://doi.org/10.1016/j.memsci.2018.11.067
    """

    CONFIG = SeawaterParameterData.CONFIG()

    def build(self):
        """
        Callable method for Block construction.
        """
        super(SeawaterVTPCParameterData, self).build()

        self._state_block_class = SeawaterVTPCStateBlock

        # ---default scaling---
        self.set_default_scaling("temperature", 1e-2)
        self.set_default_scaling("pressure", 1e-6)
        self.set_default_scaling("conc_mass_phase_comp", 1e-3, index=("Liq", "H2O"))
        self.set_default_scaling("conc_mass_phase_comp", 1e-1, index=("Liq", "TDS"))
        self.set_default_scaling("flow_vol_phase", 1e3, index="Liq")
        self.set_default_scaling("dens_mass_phase", 1e-3, index="Liq")
        self.set_default_scaling("dens_mass_solvent", 1e-3)
        self.set_default_scaling("visc_d_phase", 1e3, index="Liq")
        self.set_default_scaling("osm_coeff", 1e0)
        self.set_default_scaling("enth_mass_phase", 1e-5, index="Liq")
        self.set_default_scaling("pressure_sat", 1e-5)
        self.set_default_scaling("cp_mass_phase", 1e-3, index="Liq")
        self.set_default_scaling("therm_cond_phase", 1e0, index="Liq")
        self.set_default_scaling("dh_vap_mass", 1e-6)
        self.set_default_scaling("diffus_phase_comp", 1e9)
        self.set_default_scaling("boiling_point_elevation_phase", 1e0, index="Liq")

    @classmethod
    def define_metadata(cls, obj):
        """Define properties supported and units."""
        obj.add_properties(
            {
                "flow_mass_phase_comp": {"method": "_flow_mass_phase_comp"},
                "temperature": {"method": None},
                "pressure": {"method": None},
                "mass_frac_phase_comp": {"method": "_mass_frac_phase_comp"},
                "dens_mass_phase": {"method": "_dens_mass_phase"},
                "flow_vol_phase": {"method": None},
                "flow_vol": {"method": "_flow_vol"},
                "conc_mass_phase_comp": {"method": None},
                "flow_mol_phase_comp": {"method": "_flow_mol_phase_comp"},
                "mole_frac_phase_comp": {"method": "_mole_frac_phase_comp"},
                "molality_phase_comp": {"method": "_molality_phase_comp"},
                "visc_d_phase": {"method": "_visc_d_phase"},
                "pressure_osm_phase": {"method": "_pressure_osm_phase"},
                "enth_mass_phase": {"method": "_enth_mass_phase"},
                "pressure_sat": {"method": "_pressure_sat"},
                "cp_mass_phase": {"method": "_cp_mass_phase"},
                "therm_cond_phase": {"method": "_therm_cond_phase"},
                "diffus_phase_comp": {"method": "_diffus_phase_comp"},
            }
        )
        obj.define_custom_properties(
            {
                "dens_mass_solvent": {
                    "doc": "Mass Density of Pure Water",
                    "units": pyunits.kg * pyunits.m**-3,
                    "method": "_dens_mass_solvent",
                },
                "osm_coeff": {
                    "doc": "Osmotic Coefficient",
                    "units": pyunits.dimensionless,
                    "method": "_osm_coeff",
                },
                "enth_flow": {
                    "doc": "Enthalpy Flow",
                    "units": pyunits.J * pyunits.s**-1,
                    "method": "_enth_flow",
                },
                "dh_vap_mass": {
                    "doc": "Latent Heat of Vaporization",
                    "units": pyunits.J * pyunits.kg**-1,
                    "method": "_dh_vap_mass",
                },
                "boiling_point_elevation_phase": {
                    "doc": "Boiling Point Elevation",
                    "units": pyunits.K,
                    "method": "_boiling_point_elevation_phase",
                },
                "energy_density_phase": {
                    "doc": "Energy Density",
                    "units": pyunits.J * pyunits.m**-3,
                    "method": "_energy_density_phase",
                },
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


class _SeawaterVTPCStateBlock(StateBlock):
    """
    This Class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """

    def fix_initialization_states(self, state_args=None):
        """
        Fixes state variables for state blocks.

        Returns:
            None
        """
        # Fix state variables
        flags = fix_state_vars(self, state_args)

        # Constraint on water concentration at outlet - unfix in these cases
        for idx, b in self.items():
            if b.config.defined_state is False:
                b.conc_mass_phase_comp["Liq", "H2O"].unfix()
                # flags.pop((idx, "conc_mass_phase_comp", ("Liq", "H2O")))

        return flags

    def initialize(
        self,
        state_args=None,
        state_vars_fixed=False,
        hold_state=False,
        outlvl=idaeslog.NOTSET,
        solver="ipopt_v2",
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
        if solver is None:
            # Solver might get passed as None
            # from a parent initialization method
            solver = "ipopt_v2"

        # Fix state variables
        flags = self.fix_initialization_states(state_args)
        # Check when the state vars are fixed already result in dof 0
        for k in self.keys():
            dof = degrees_of_freedom(self[k])
            if dof != 0:
                raise PropertyPackageError(
                    "State vars fixed but degrees of "
                    "freedom for state block is not "
                    "zero during initialization."
                )

        init_obj = BlockTriangularizationInitializer(
            block_solver=solver, block_solver_options=optarg, output_level=outlvl
        )
        for blkdata in self.values():
            if not blkdata.active:
                continue
            try:
                init_obj.initialize(blkdata)
            except NoSolutionError as err:
                # The IDAES testing harness wants this exact error message
                raise InitializationError(
                    "fs.props failed to initialize successfully. "
                    "Please check the output logs for more information."
                ) from err

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
        init_log.info("{} State Released.".format(self.name))

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
            self.fix_initialization_states(self)

        return results


@declare_process_block_class(
    "SeawaterVTPCStateBlock", block_class=_SeawaterVTPCStateBlock
)
class SeawaterVTPCStateBlockData(SeawaterStateBlockData):
    """A seawater property package."""

    def build(self):
        """Callable method for Block construction."""
        # Do not call build from immediate parent class
        # but instead from its parent (StateBlockData)
        super(SeawaterStateBlockData, self).build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        # Add state variables
        self.flow_vol_phase = Var(
            self.params.phase_list,
            initialize=0.000977,
            bounds=(0.0, None),
            units=pyunits.m**3 / pyunits.s,
            doc="Volumetric flow rate",
        )

        self.conc_mass_phase_comp = Var(
            self.params.phase_list,
            self.params.component_list,
            initialize={
                ("Liq", "H2O"): 987.7,
                ("Liq", "TDS"): 35.82,
            },
            bounds=(0.0, 1e6),
            units=pyunits.kg * pyunits.m**-3,
            doc="Mass concentration",
        )

        self.temperature = Var(
            initialize=298.15,
            bounds=(273.15, 1000),
            units=pyunits.K,
            doc="Temperature",
        )

        self.pressure = Var(
            initialize=101325,
            bounds=(1e3, 5e7),
            units=pyunits.Pa,
            doc="Pressure",
        )

        if not self.config.defined_state:

            @self.Constraint(self.phase_list)
            def eq_sum_conc_mass_phase_comp(b, p):
                return (
                    sum(b.conc_mass_phase_comp[p, j] for j in b.component_list)
                    == b.dens_mass_phase[p]
                )

    # -----------------------------------------------------------------------------
    # Property Methods

    def _flow_vol_phase(self):
        # Overwrite method from base class.
        raise BurntToast(
            "This function should never be called because flow_vol_phase should "
            "have been created as a state variable. If it does not exist as a "
            "state variable, please report this bug to the WaterTAP developers."
        )

    def _conc_mass_phase_comp(self):
        # Overwrite method from base class.
        raise BurntToast(
            "This function should never be called because conc_mass_phase_comp should "
            "have been created as a state variable. If it does not exist as a "
            "state variable, please report this bug to the WaterTAP developers."
        )

    def _flow_mass_phase_comp(self):
        self.flow_mass_phase_comp = Var(
            self.params.phase_list,
            self.params.component_list,
            initialize={("Liq", "H2O"): 0.965, ("Liq", "TDS"): 0.035},
            bounds=(0.0, None),
            units=pyunits.kg / pyunits.s,
            doc="Mass flow rate",
        )

        def rule_flow_mass_phase_comp(b, p, j):
            return b.flow_mass_phase_comp[p, j] == (
                b.flow_vol_phase[p] * b.conc_mass_phase_comp[p, j]
            )

        self.eq_flow_mass_phase_comp = Constraint(
            self.params.phase_list,
            self.params.component_list,
            rule=rule_flow_mass_phase_comp,
        )

    def _mass_frac_phase_comp(self):
        self.mass_frac_phase_comp = Var(
            self.params.phase_list,
            self.params.component_list,
            initialize=0.1,
            bounds=(0.0, None),
            units=pyunits.dimensionless,
            doc="Mass fraction",
        )

        def rule_mass_frac_phase_comp(b, p, j):
            return b.mass_frac_phase_comp[p, j] == b.conc_mass_phase_comp[p, j] / sum(
                b.conc_mass_phase_comp[p, j] for j in b.params.component_list
            )

        self.eq_mass_frac_phase_comp = Constraint(
            self.params.phase_list,
            self.params.component_list,
            rule=rule_mass_frac_phase_comp,
        )

    def define_state_vars(self):
        """Define state vars."""
        return {
            "flow_vol_phase": self.flow_vol_phase,
            "conc_mass_phase_comp": self.conc_mass_phase_comp,
            "temperature": self.temperature,
            "pressure": self.pressure,
        }

    # -----------------------------------------------------------------------------
    # Scaling methods
    def calculate_scaling_factors(self):
        # Don't call method from immediate parent class, but the
        # one from that class's parent class (StateBlockData)
        super(SeawaterStateBlockData, self).calculate_scaling_factors()

        # setting scaling factors for variables

        # default scaling factors have already been set with
        # idaes.core.property_base.calculate_scaling_factors()
        # for the following variables: flow_vol_phase, conc_mass_phase_comp, pressure,
        # temperature, dens_mass_phase, visc_d_phase, osm_coeff, and enth_mass_phase

        if self.is_property_constructed("flow_mass_phase_comp"):
            for (p, j), vardata in self.flow_mass_phase_comp.items():
                # These scaling factors should have been set either by the
                # user directly or from the default scaling factors
                sf_C = iscale.get_scaling_factor(
                    self.conc_mass_phase_comp[p, j], warning=True
                )
                sf_V = iscale.get_scaling_factor(self.flow_vol_phase[p], warning=True)
                iscale.set_scaling_factor(vardata, sf_C * sf_V, overwrite=False)

        # scaling factors for parameters
        # TODO this should live on the parameter block.
        for j, v in self.params.mw_comp.items():
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(self.params.mw_comp, 1e2)

        # these variables do not typically require user input,
        # will not override if the user does provide the scaling factor
        if self.is_property_constructed("pressure_osm_phase"):
            if iscale.get_scaling_factor(self.pressure_osm_phase["Liq"]) is None:
                iscale.set_scaling_factor(
                    self.pressure_osm_phase["Liq"],
                    iscale.get_scaling_factor(self.pressure),
                )

        if self.is_property_constructed("mass_frac_phase_comp"):
            for j in self.params.component_list:
                if (
                    iscale.get_scaling_factor(self.mass_frac_phase_comp["Liq", j])
                    is None
                ):
                    if j == "TDS":
                        sf = iscale.get_scaling_factor(
                            self.conc_mass_phase_comp["Liq", j]
                        ) / iscale.get_scaling_factor(
                            self.conc_mass_phase_comp["Liq", "H2O"]
                        )
                        iscale.set_scaling_factor(
                            self.mass_frac_phase_comp["Liq", j], sf
                        )
                    elif j == "H2O":
                        iscale.set_scaling_factor(
                            self.mass_frac_phase_comp["Liq", j], 1
                        )

        if self.is_property_constructed("flow_vol"):
            sf = iscale.get_scaling_factor(self.flow_vol_phase["Liq"])
            iscale.set_scaling_factor(self.flow_vol, sf)

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
                    if j == "TDS":
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
            iscale.set_scaling_factor(
                self.enth_flow,
                iscale.get_scaling_factor(self.flow_mass_phase_comp["Liq", "H2O"])
                * iscale.get_scaling_factor(self.enth_mass_phase["Liq"]),
            )

        if self.is_property_constructed("boiling_point_elevation_phase"):
            iscale.set_scaling_factor(self.boiling_point_elevation_phase["Liq"], 1)

        # transforming constraints
        if self.is_property_constructed("eq_sum_conc_mass_phase_comp"):
            for p, condata in self.eq_sum_conc_mass_phase_comp.items():
                sf_rho = iscale.get_scaling_factor(
                    self.dens_mass_phase, default=1e-3, warning=True
                )
                iscale.constraint_scaling_transform(condata, sf_rho)
        # property relationships with no index, simple constraint
        v_str_lst_simple = [
            "dens_mass_solvent",
            "osm_coeff",
            "pressure_sat",
            "dh_vap_mass",
        ]
        for v_str in v_str_lst_simple:
            if self.is_property_constructed(v_str):
                v = getattr(self, v_str)
                sf = iscale.get_scaling_factor(v, default=1, warning=True)
                c = getattr(self, "eq_" + v_str)
                iscale.constraint_scaling_transform(c, sf)

        if self.is_property_constructed("pressure_osm_phase"):
            sf = iscale.get_scaling_factor(
                self.pressure_osm_phase["Liq"], default=1, warning=True
            )
            iscale.constraint_scaling_transform(self.eq_pressure_osm_phase["Liq"], sf)

        # transforming constraints
        transform_property_constraints(self)
