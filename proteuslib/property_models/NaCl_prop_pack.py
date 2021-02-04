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
Property package for H2O-NaCl solution
"""

# Import Pyomo modules
from pyomo.environ import Constraint, Expression, log, Reals, NonNegativeReals, \
    Var, Set, Param, sqrt, log10, TerminationCondition, Suffix
from pyomo.environ import units as pyunits
from pyomo.opt import SolverFactory

# Import IDAES modules
import idaes.logger as idaeslog
from idaes.core import (declare_process_block_class,
                        MaterialFlowBasis,
                        PhysicalParameterBlock,
                        StateBlockData,
                        StateBlock,
                        MaterialBalanceType,
                        EnergyBalanceType)
from idaes.core.components import Component, Solute, Solvent
from idaes.core.phases import LiquidPhase
from idaes.core.util.constants import Constants
from idaes.core.util.initialization import (fix_state_vars,
                                            revert_state_vars,
                                            solve_indexed_blocks)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.exceptions import PropertyPackageError
import idaes.core.util.scaling as iscale

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("NaClParameterBlock")
class NaClParameterData(PhysicalParameterBlock):
    CONFIG = PhysicalParameterBlock.CONFIG()

    def build(self):
        '''
        Callable method for Block construction.
        '''
        super(NaClParameterData, self).build()

        self._state_block_class = NaClStateBlock

        # components
        self.H2O = Solvent()
        self.NaCl = Solute()

        # phases
        self.Liq = LiquidPhase()

        # reference
        # this package is developed from Bartholomew & Mauter (2019) https://doi.org/10.1016/j.memsci.2018.11.067
        # the enthalpy calculations are from Sharqawy et al. (2010) http://dx.doi.org/10.5004/dwt.2010.1079

        # molecular weight
        mw_comp_data = {'H2O': 18.01528E-3,
                        'NaCl': 58.44E-3}
        self.mw_comp = Param(self.component_list,
                             mutable=False,
                             initialize=mw_comp_data,
                             units=pyunits.kg / pyunits.mol,
                             doc="Molecular weight")

        # mass density parameters, eq 4 in Bartholomew
        self.dens_mass_param_0 = Var(
            within=Reals, initialize=995, units=pyunits.kg / pyunits.m**3,
            doc='Mass density parameter 0')
        self.dens_mass_param_1 = Var(
            within=Reals, initialize=756, units=pyunits.kg / pyunits.m**3,
            doc='Mass density parameter 1')

        # dynamic viscosity parameters, eq 5 in Bartholomew
        self.visc_d_param_0 = Var(
            within=Reals, initialize=9.80E-4, units=pyunits.Pa * pyunits.s,
            doc='Dynamic viscosity parameter 0')
        self.visc_d_param_1 = Var(
            within=Reals, initialize=2.15E-3, units=pyunits.Pa * pyunits.s,
            doc='Dynamic viscosity parameter 1')

        # diffusivity parameters, eq 6 in Bartholomew
        self.diffus_param_0 = Var(
            within=Reals, initialize=1.51e-9, units=pyunits.m**2 / pyunits.s,
            doc='Dynamic viscosity parameter 0')
        self.diffus_param_1 = Var(
            within=Reals, initialize=-2.00e-9, units=pyunits.m**2 / pyunits.s,
            doc='Dynamic viscosity parameter 1')
        self.diffus_param_2 = Var(
            within=Reals, initialize=3.01e-8, units=pyunits.m**2 / pyunits.s,
            doc='Dynamic viscosity parameter 2')
        self.diffus_param_3 = Var(
            within=Reals, initialize=-1.22e-7, units=pyunits.m**2 / pyunits.s,
            doc='Dynamic viscosity parameter 3')
        self.diffus_param_4 = Var(
            within=Reals, initialize=1.53e-7, units=pyunits.m**2 / pyunits.s,
            doc='Dynamic viscosity parameter 4')

        # osmotic coefficient parameters, eq. 3b in Bartholomew
        self.osm_coeff_param_0 = Var(
            within=Reals, initialize=0.918, units=pyunits.dimensionless,
            doc='Osmotic coefficient parameter 0')
        self.osm_coeff_param_1 = Var(
            within=Reals, initialize=8.89e-2, units=pyunits.dimensionless,
            doc='Osmotic coefficient parameter 1')
        self.osm_coeff_param_2 = Var(
            within=Reals, initialize=4.92, units=pyunits.dimensionless,
            doc='Osmotic coefficient parameter 2')

        # TODO: update for NaCl solution, relationship from Sharqawy is for seawater
        # specific enthalpy parameters, eq. 55 and 43 in Sharqawy (2010)
        self.enth_mass_param_A1 = Var(
            within=Reals, initialize=124.790, units=pyunits.J / pyunits.kg,
            doc='Specific enthalpy parameter A1')
        self.enth_mass_param_A2 = Var(
            within=Reals, initialize=4203.075, units=(pyunits.J / pyunits.kg) * pyunits.K ** -1,
            doc='Specific enthalpy parameter A2')
        self.enth_mass_param_A3 = Var(
            within=Reals, initialize=-0.552, units=(pyunits.J / pyunits.kg) * pyunits.K ** -2,
            doc='Specific enthalpy parameter A3')
        self.enth_mass_param_A4 = Var(
            within=Reals, initialize=0.004, units=(pyunits.J / pyunits.kg) * pyunits.K ** -3,
            doc='Specific enthalpy parameter A4')
        self.enth_mass_param_B1 = Var(
            within=Reals, initialize=27062.623, units=pyunits.dimensionless,
            doc='Specific enthalpy parameter B1')
        self.enth_mass_param_B2 = Var(
            within=Reals, initialize=4835.675, units=pyunits.dimensionless,
            doc='Specific enthalpy parameter B2')

        # traditional parameters are the only Vars currently on the block and should be fixed
        for v in self.component_objects(Var):
            v.fix()

        # ---default scaling---
        self.set_default_scaling('temperature', 1e-2)
        self.set_default_scaling('pressure', 1e-6)
        self.set_default_scaling('dens_mass', 1e-3)
        self.set_default_scaling('visc_d', 1e3)
        self.set_default_scaling('diffus', 1e9)
        self.set_default_scaling('osm_coeff', 1e0)
        self.set_default_scaling('enth_mass', 1e-5)

    @classmethod
    def define_metadata(cls, obj):
        """Define properties supported and units."""
        obj.add_properties(
            {'flow_mass_comp': {'method': None},
             'temperature': {'method': None},
             'pressure': {'method': None},
             'mass_frac_comp': {'method': '_mass_frac_comp'},
             'dens_mass': {'method': '_dens_mass'},
             'flow_vol': {'method': '_flow_vol'},
             'visc_d': {'method': '_visc_d'},
             'diffus': {'method': '_diffus'},
             'conc_mass_comp': {'method': '_conc_mass_comp'},
             'osm_coeff': {'method': '_osm_coeff'},
             'pressure_osm': {'method': '_pressure_osm'},
             'enth_mass': {'method': '_enth_mass'},
             'enth_flow': {'method': '_enth_flow'}
             })

        obj.add_default_units({'time': pyunits.s,
                               'length': pyunits.m,
                               'mass': pyunits.kg,
                               'amount': pyunits.mol,
                               'temperature': pyunits.K})


class _NaClStateBlock(StateBlock):
    """
    This Class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """

    def initialize(self, state_args={}, state_vars_fixed=False,
                   hold_state=False, outlvl=idaeslog.NOTSET,
                   solver='ipopt', optarg={'tol': 1e-8}):
        """
        Initialization routine for property package.
        Keyword Arguments:
            state_args : Dictionary with initial guesses for the state vars
                         chosen. Note that if this method is triggered
                         through the control volume, and if initial guesses
                         were not provided at the unit model level, the
                         control volume passes the inlet values as initial
                         guess.The keys for the state_args dictionary are:

                         flow_mol_phase_comp : value at which to initialize
                                               phase component flows
                         pressure : value at which to initialize pressure
                         temperature : value at which to initialize temperature
            outlvl : sets output level of initialization routine
                     * 0 = no output (default)
                     * 1 = return solver state for each step in routine
                     * 2 = include solver output information (tee=True)
            optarg : solver options dictionary object (default={'tol': 1e-8})
            state_vars_fixed: Flag to denote if state vars have already been
                              fixed.
                              - True - states have already been fixed by the
                                       control volume 1D. Control volume 0D
                                       does not fix the state vars, so will
                                       be False if this state block is used
                                       with 0D blocks.
                             - False - states have not been fixed. The state
                                       block will deal with fixing/unfixing.
            solver : str indicating which solver to use during
                     initialization (default = 'ipopt')
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
        opt = SolverFactory(solver)
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

@declare_process_block_class("NaClStateBlock",
                             block_class=_NaClStateBlock)
class NaClStateBlockData(StateBlockData):
    def build(self):
        """Callable method for Block construction."""
        super(NaClStateBlockData, self).build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        # Add state variables
        self.flow_mass_comp = Var(
            self.params.component_list,
            initialize=1,
            bounds=(1e-8, 100),
            domain=NonNegativeReals,
            units=pyunits.kg/pyunits.s,
            doc='Mass flow rate [kg/s]')

        self.temperature = Var(
            initialize=298.15,
            bounds=(273.15, 1000),
            domain=NonNegativeReals,
            units=pyunits.degK,
            doc='State temperature [K]')

        self.pressure = Var(
            initialize=101325,
            bounds=(1e5, 5e7),
            domain=NonNegativeReals,
            units=pyunits.Pa,
            doc='State pressure [Pa]')

    # -----------------------------------------------------------------------------
    # Property Methods
    def _mass_frac_comp(self):
        self.mass_frac_comp = Var(
            self.params.component_list,
            initialize=0.1,
            bounds=(1e-8, 1),
            units=pyunits.dimensionless,
            doc='Mass fraction')

        def rule_mass_frac_comp(b, j):
            return (b.mass_frac_comp[j] == b.flow_mass_comp[j] /
                    sum(b.flow_mass_comp[j]
                        for j in self.params.component_list))
        self.eq_mass_frac_comp = Constraint(self.params.component_list, rule=rule_mass_frac_comp)

    def _dens_mass(self):
        self.dens_mass = Var(
            initialize=1e3,
            bounds=(1, 1e6),
            units=pyunits.kg * pyunits.m**-3,
            doc="Mass density")

        def rule_dens_mass(b):  # density, eq. 4 in Bartholomew
            return b.dens_mass == (b.params.dens_mass_param_1 * b.mass_frac_comp['NaCl']
                                   + b.params.dens_mass_param_0)
        self.eq_dens_mass = Constraint(rule=rule_dens_mass)

    def _flow_vol(self):
        self.flow_vol = Var(
            initialize=1,
            bounds=(1e-8, 1e8),
            units=pyunits.m**3 / pyunits.s,
            doc="Volumetric flowrate")

        def rule_flow_vol(b):
            return (b.flow_vol == sum(b.flow_mass_comp[j] for j in self.params.component_list)
                    / b.dens_mass)
        self.eq_flow_vol = Constraint(rule=rule_flow_vol)

    def _conc_mass_comp(self):
        self.conc_mass_comp = Var(
            self.params.component_list,
            initialize=10,
            bounds=(1e-6, 1e6),
            units=pyunits.kg * pyunits.m**-3,
            doc="Mass concentration")

        def rule_conc_mass_comp(b, j):
            return self.conc_mass_comp[j] == \
                   self.dens_mass * self.mass_frac_comp[j]
        self.eq_conc_mass_comp = Constraint(self.params.component_list, rule=rule_conc_mass_comp)

    def _visc_d(self):
        self.visc_d = Var(
            initialize=1e-3,
            bounds=(1e-8, 1),
            units=pyunits.Pa * pyunits.s,
            doc="Viscosity")

        def rule_visc_d(b):  # dynamic viscosity, eq 5 in Bartholomew
            return b.visc_d == (b.params.visc_d_param_1 * b.mass_frac_comp['NaCl']
                                + b.params.visc_d_param_0)
        self.eq_visc_d = Constraint(rule=rule_visc_d)

    def _diffus(self):
        self.diffus = Var(
            initialize=1e-9,
            bounds=(1e-12, 1e-6),
            units=pyunits.m**2 * pyunits.s**-1,
            doc="Diffusivity")

        def rule_diffus(b):  # diffusivity, eq 6 in Bartholomew
            return b.diffus == (b.params.diffus_param_4 * b.mass_frac_comp['NaCl'] ** 4
                                + b.params.diffus_param_3 * b.mass_frac_comp['NaCl'] ** 3
                                + b.params.diffus_param_2 * b.mass_frac_comp['NaCl'] ** 2
                                + b.params.diffus_param_1 * b.mass_frac_comp['NaCl']
                                + b.params.diffus_param_0)
        self.eq_diffus = Constraint(rule=rule_diffus)

    def _osm_coeff(self):
        self.osm_coeff = Var(
            initialize=1,
            bounds=(1e-8, 10),
            units=pyunits.dimensionless,
            doc="Osmotic coefficient")

        def rule_osm_coeff(b):
            return b.osm_coeff == (b.params.osm_coeff_param_2 * b.mass_frac_comp['NaCl'] ** 2
                                   + b.params.osm_coeff_param_1 * b.mass_frac_comp['NaCl']
                                   + b.params.osm_coeff_param_0)
        self.eq_osm_coeff = Constraint(rule=rule_osm_coeff)

    def _pressure_osm(self):
        self.pressure_osm = Var(
            initialize=1e6,
            bounds=(1, 1e8),
            units=pyunits.Pa,
            doc="Osmotic pressure")

        def rule_pressure_osm(b):
            i = 2  # number of ionic species
            molality = b.mass_frac_comp['NaCl'] / (1 - b.mass_frac_comp['NaCl']) / b.params.mw_comp['NaCl']
            rhow = 1000 * pyunits.kg/pyunits.m**3  # TODO: could make this variable based on temperature
            return (b.pressure_osm ==
                    i * b.osm_coeff * molality * rhow * Constants.gas_constant * b.temperature)
        self.eq_pressure_osm = Constraint(rule=rule_pressure_osm)

    def _enth_mass(self):
        self.enth_mass = Var(
            initialize=1e6,
            bounds=(1, 1e9),
            units=pyunits.J * pyunits.kg ** -1,
            doc="Specific enthalpy")

        def rule_enth_mass(b):  # specific enthalpy, eq. 55 and 43 in Sharqawy
            t = b.temperature - 273.15 * pyunits.K  # temperature in degC, but pyunits in K
            S = b.mass_frac_comp['NaCl']
            h_w = (b.params.enth_mass_param_A1
                   + b.params.enth_mass_param_A2 * t
                   + b.params.enth_mass_param_A3 * t ** 2
                   + b.params.enth_mass_param_A4 * t ** 3)
            # relationship requires dimensionless calculation and units added at end
            h_sw = (h_w -
                    (S * (b.params.enth_mass_param_B1 + S)
                     + S * (b.params.enth_mass_param_B2 + S) * t / pyunits.K)
                    * pyunits.J / pyunits.kg)
            return b.enth_mass == h_sw

        self.eq_enth_mass = Constraint(rule=rule_enth_mass)

    def _enth_flow(self):
        # enthalpy flow expression for get_enthalpy_flow_terms method

        def rule_enth_flow(b):  # enthalpy flow
            return sum(b.flow_mass_comp[j] for j in b.params.component_list) * b.enth_mass
        self.enth_flow = Expression(rule=rule_enth_flow)

    # TODO: add vapor pressure, specific heat, thermal conductivity,
    #   and heat of vaporization

    # -----------------------------------------------------------------------------
    # General Methods
    # NOTE: For scaling in the control volume to work properly, these methods must
    # return a pyomo Var or Expression

    def get_material_flow_terms(self, p, j):
        """Create material flow terms for control volume."""
        return self.flow_mass_comp[j]

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

    def get_material_flow_basis(b):
        return MaterialFlowBasis.mass

    def define_state_vars(self):
        """Define state vars."""
        return {"flow_mass_comp": self.flow_mass_comp,
                "temperature": self.temperature,
                "pressure": self.pressure}

    # -----------------------------------------------------------------------------
    # Scaling methods
    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # setting scaling factors for variables

        # default scaling factors have already been set with
        # idaes.core.property_base.calculate_scaling_factors()
        # for the following variables: flow_mass_comp, pressure,
        # temperature, dens_mass, visc_d, diffus, osm_coeff, and enth_mass

        # these variables should have user input
        if iscale.get_scaling_factor(self.flow_mass_comp['H2O']) is None:
            sf = iscale.get_scaling_factor(self.flow_mass_comp['H2O'], default=1e0, warning=True)
            iscale.set_scaling_factor(self.flow_mass_comp['H2O'], sf)

        if iscale.get_scaling_factor(self.flow_mass_comp['NaCl']) is None:
            sf = iscale.get_scaling_factor(self.flow_mass_comp['NaCl'], default=1e2, warning=True)
            iscale.set_scaling_factor(self.flow_mass_comp['NaCl'], sf)

        # these variables do not typically require user input,
        # will not override if the user does provide the scaling factor
        if self.is_property_constructed('pressure_osm'):
            if iscale.get_scaling_factor(self.pressure_osm) is None:
                iscale.set_scaling_factor(self.pressure_osm,
                                          iscale.get_scaling_factor(self.pressure))

        if self.is_property_constructed('mass_frac_comp'):
            for j in self.params.component_list:
                if iscale.get_scaling_factor(self.mass_frac_comp[j]) is None:
                    if j == 'H2O':
                        iscale.set_scaling_factor(self.mass_frac_comp[j], 1)
                    elif j == 'NaCl':
                        iscale.set_scaling_factor(
                            self.mass_frac_comp[j],
                            iscale.get_scaling_factor(self.flow_mass_comp['NaCl'])
                            / iscale.get_scaling_factor(self.flow_mass_comp['H2O']))

        if self.is_property_constructed('conc_mass_comp'):
            for j in self.params.component_list:
                sf_dens = iscale.get_scaling_factor(self.dens_mass)
                if iscale.get_scaling_factor(self.conc_mass_comp[j]) is None:
                    if j == 'H2O':
                        # solvents typically have a mass fraction between 0.5-1
                        iscale.set_scaling_factor(self.conc_mass_comp[j], sf_dens)
                    elif j == 'NaCl':
                        iscale.set_scaling_factor(
                            self.conc_mass_comp[j],
                            sf_dens * iscale.get_scaling_factor(self.mass_frac_comp[j]))

        if self.is_property_constructed('flow_vol'):
            iscale.set_scaling_factor(self.flow_vol,
                                      iscale.get_scaling_factor(self.flow_mass_comp['H2O'])
                                      / iscale.get_scaling_factor(self.dens_mass))

        if self.is_property_constructed('enth_flow'):
            iscale.set_scaling_factor(self.enth_flow,
                                      iscale.get_scaling_factor(self.flow_mass_comp['H2O'])
                                      * iscale.get_scaling_factor(self.enth_mass))

        # transforming constraints
        # property relationships with no index, simple constraint
        v_str_lst_simple = ['dens_mass', 'flow_vol', 'visc_d', 'diffus', 'osm_coeff', 'pressure_osm', 'enth_mass']
        for v_str in v_str_lst_simple:
            if self.is_property_constructed(v_str):
                v = getattr(self, v_str)
                sf = iscale.get_scaling_factor(v, default=1, warning=True)
                c = getattr(self, 'eq_' + v_str)
                iscale.constraint_scaling_transform(c, sf)

        # property relationships indexed by component
        v_str_lst_comp = ['mass_frac_comp', 'conc_mass_comp']
        for v_str in v_str_lst_comp:
            if self.is_property_constructed(v_str):
                v_comp = getattr(self, v_str)
                c_comp = getattr(self, 'eq_' + v_str)
                for j, c in c_comp.items():
                    sf = iscale.get_scaling_factor(v_comp[j], default=1, warning=True)
                    iscale.constraint_scaling_transform(c, sf)