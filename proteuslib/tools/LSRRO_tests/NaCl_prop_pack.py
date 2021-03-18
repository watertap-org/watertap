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
Initial property package for H2O-NaCl system
"""

# Import Python libraries
import idaes.logger as idaeslog

# Import Pyomo libraries
from pyomo.environ import Constraint, Expression, log, Reals, NonNegativeReals, \
    Var, Set, Param, sqrt, log10, TerminationCondition, Suffix
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
from idaes.core.util.constants import Constants
from idaes.core.util.initialization import (fix_state_vars,
                                            revert_state_vars,
                                            solve_indexed_blocks)
from idaes.core.util.misc import add_object_reference, extract_data
from idaes.core.util.model_statistics import degrees_of_freedom, \
    number_unfixed_variables
from idaes.core.util.exceptions import ConfigurationError, PropertyPackageError
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
        # self.component_list = Set(initialize=['H2O', 'NaCl'])  # deprecated, automatically added with comp objects
        self.H2O = Solvent()
        self.NaCl = Solute()

        # phases
        # self.phase_list = Set(initialize=['Liq'])  # deprecated, automatically added with phase objects
        self.Liq = LiquidPhase()

        # parameters
        # This file creates most parameters as fixed variables, so that they can be unfixed for parameter estimation.
        # Further, this file demonstrates two approaches that handle units differently. The first approach is to assign
        # units directly to each parameter. This approach may lead to numerous declarations for regression models and
        # is demonstrated for the dynamic viscosity parameters. The second approach is designed for regression models,
        # and involves declaring all regression parameters as unitless and including the final units for the property
        # in the regression equation. This second approach is demonstrated by all of the other properties.

        # molecular weight
        mw_comp_data = {'H2O': 18.01528E-3,
                        'NaCl': 58.44E-3}
        self.mw_comp = Param(self.component_list,
                             mutable=False,
                             initialize=extract_data(mw_comp_data),
                             units=pyunits.kg/pyunits.mol,
                             doc="Molecular weight kg/mol")

        # mass density parameters, unitless
        dens_mass_params_dic = {'A_1': 2, 'A_2': -200, 'A_3': 160,  # term A
                                'B_1': 2, 'B_2': -150, 'B_3': 150,  # term B
                                'F1': 0.5,  # term F1
                                'F3_1': 2, 'F3_2': -1,  # term F3
                                'F4_1': 4, 'F4_2': -3,  # term F4
                                'G1': 0.5,  # term G1
                                'G3_1': 2, 'G3_2': -1,  # term G3
                                'A1_1': 4.032, 'A1_2': 0.115, 'A1_3': 3.26E-4,  # term A1
                                'A2_1': -0.108, 'A2_2': 1.571e-3, 'A2_3': -4.23e-4,  # term A2
                                'A3_1': -0.012, 'A3_2': 1.74e-3, 'A3_3': -9e-6,  # term A3
                                'A4_1': 6.92e-4, 'A4_2': -8.7e-5, 'A4_3': -5.3e-5,  # term A4
                                }
        self.dens_mass_params = Var(
            [*dens_mass_params_dic],
            within=Reals, initialize=dens_mass_params_dic, units=None,
            doc='Mass density parameters')

        # dynamic viscosity parameters, with units
        self.visc_d_muw_A = Var(
            within=Reals, initialize=4.2844e-5, units=pyunits.Pa * pyunits.s,
            doc='Dynamic viscosity parameter A for pure water')
        self.visc_d_muw_B = Var(
            within=Reals, initialize=0.157, units=pyunits.degK**-2 * pyunits.Pa**-1 * pyunits.s**-1,
            doc='Dynamic viscosity parameter B for pure water')
        self.visc_d_muw_C = Var(
            within=Reals, initialize=64.993, units=pyunits.degK,
            doc='Dynamic viscosity parameter C for pure water')
        self.visc_d_muw_D = Var(
            within=Reals, initialize=91.296, units=pyunits.Pa**-1 * pyunits.s**-1,
            doc='Dynamic viscosity parameter D for pure water')
        self.visc_d_A_1 = Var(
            within=Reals, initialize=1.541, units=None,
            doc='Dynamic viscosity parameter 1 for term A')
        self.visc_d_A_2 = Var(
            within=Reals, initialize=1.998e-2, units=pyunits.degK**-1,
            doc='Dynamic viscosity parameter 2 for term A')
        self.visc_d_A_3 = Var(
            within=Reals, initialize=-9.52e-5, units=pyunits.degK**-2,
            doc='Dynamic viscosity parameter 3 for term A')
        self.visc_d_B_1 = Var(
            within=Reals, initialize=7.974, units=None,
            doc='Dynamic viscosity parameter 1 for term B')
        self.visc_d_B_2 = Var(
            within=Reals, initialize=-7.561e-2, units=pyunits.degK**-1,
            doc='Dynamic viscosity parameter 2 for term B')
        self.visc_d_B_3 = Var(
            within=Reals, initialize=4.724e-4, units=pyunits.degK**-2,
            doc='Dynamic viscosity parameter 3 for term B')

        # diffusivity parameters, unitless
        diffus_params_dic = {'A': 3.847e-4, 'B': -0.1984, 'C': 26.54}
        self.diffus_params = Var([*diffus_params_dic],
            within=Reals, initialize=diffus_params_dic, units=None,
            doc='Diffusivity parameters')

        # osmotic coefficient parameters, unitless
        osm_coeff_dic = {'1': 8.9453e-1,
                         '2': 4.1561e-4,
                         '3': -4.6262e-6,
                         '4': 2.2211e-11,
                         '5': -1.1445e-1,
                         '6': -1.4783e-3,
                         '7': -1.3526e-8,
                         '8': 7.0132,
                         '9': 5.696e-2,
                         '10': -2.8624e-4}
        self.osm_coeff_params = Var(
            [*osm_coeff_dic],
            within=Reals, initialize=osm_coeff_dic, units=None,
            doc='Osmotic coefficient parameters')

        # specific enthalpy parameters, unitless
        enth_mass_dic = {'A1': 124.790, 'A2': 4203.075, 'A3': -0.552, 'A4': 0.004,
                         'B1': 27062.623, 'B2': 4835.675}
        self.enth_mass_params = Var(
            [*enth_mass_dic],
            within=Reals, initialize=enth_mass_dic, units=None,
            doc='Specific enthalpy parameters')

        # traditional parameters are the only Vars on the block and should be fixed
        for v in self.component_objects(Var, descend_into=True):
            for i in v:
                if v[i].value is None:
                    raise ConfigurationError(
                        "{} parameter {} was not assigned"
                        " a value. Please check your configuration "
                        "arguments.".format(self.name, v.local_name))
                v[i].fix()

        # ---default scaling---
        self.set_default_scaling('temperature', 1e-2)
        self.set_default_scaling('dens_mass', 1e-3)
        self.set_default_scaling('visc_d', 1e3)
        self.set_default_scaling('diffus', 1e9)
        self.set_default_scaling('osm_coeff', 1e0)
        self.set_default_scaling('enth_mass', 1e-5)
        # TODO: require users to provide the following scaling factors
        self.set_default_scaling('flow_mass_comp', 1e0, 'H2O')
        self.set_default_scaling('flow_mass_comp', 1e2, 'NaCl')
        self.set_default_scaling('pressure', 1e-6)

    @classmethod
    def define_metadata(cls, obj):
        """Define properties supported and units."""
        obj.add_properties(
            {'flow_mass_comp': {'method': None, 'units': pyunits.kg/pyunits.s},
             'temperature': {'method': None, 'units': pyunits.degK},
             'pressure': {'method': None, 'units': pyunits.Pa},
             'mass_frac_comp': {'method': '_mass_frac_comp', 'units': None},
             'dens_mass': {'method': '_dens_mass', 'units': pyunits.kg/pyunits.m**3},
             'flow_vol': {'method': '_flow_vol', 'units': pyunits.m**3/pyunits.s},
             'visc_d': {'method': '_visc_d', 'units': pyunits.Pa*pyunits.s},
             'diffus': {'method': '_diffus', 'units': pyunits.m**2/pyunits.s},
             'conc_mass_comp': {'method': '_conc_mass_comp', 'units': pyunits.kg/pyunits.m**3},
             'osm_coeff': {'method': '_osm_coeff', 'units': None},
             'pressure_osm': {'method': '_pressure_osm', 'units': pyunits.Pa},
             'enth_mass': {'method': '_enth_mass', 'units': pyunits.J/pyunits.kg},
             'enth_flow': {'method': '_enth_flow', 'units': pyunits.J/pyunits.s}
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

    def initialize(blk, state_args={}, state_vars_fixed=False,
                   hold_state=False, outlvl=1,
                   solver='ipopt', optarg={'tol': 1e-8}):
        """
        Initialization routine for property package.
        Keyword Arguments:
            state_args : Dictionary with initial guesses for the state vars
                         chosen. Note that if this method is triggered
                         through the control volume, and if initial guesses
                         were not provied at the unit model level, the
                         control volume passes the inlet values as initial
                         guess.The keys for the state_args dictionary are:

                         flow_mol_phase_comp : value at which to initialize
                                               phase component flows
                         pressure : value at which to initialize pressure
                         temperature : value at which to initialize temperature
            outlvl : sets output level of initialization routine
                     * 0 = no output (default)
                     * 1 = return solver state for each step in routine
                     * 2 = include solver output infomation (tee=True)
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
            solver : str indicating whcih solver to use during
                     initialization (default = 'ipopt')
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
        """

        _log.info('Starting {} initialization'.format(blk.name))

        # Fix state variables if not already fixed
        if state_vars_fixed is False:
            flags = fix_state_vars(blk, state_args)

        else:
            # Check when the state vars are fixed already result in dof 0
            for k in blk.keys():
                if degrees_of_freedom(blk[k]) != 0:
                    raise Exception("State vars fixed but degrees of freedom "
                                    "for state block is not zero during "
                                    "initialization.")
        # Set solver options
        if outlvl > 1:
            stee = True
        else:
            stee = False

        if optarg is None:
            sopt = {'tol': 1e-8}
        else:
            sopt = optarg

        opt = SolverFactory('ipopt')
        opt.options = sopt
        # ---------------------------------------------------------------------
        # Initialize flow rates and compositions

        free_vars = 0
        for k in blk.keys():
            free_vars += number_unfixed_variables(blk[k])
        if free_vars > 0:
            try:
                results = solve_indexed_blocks(opt, [blk], tee=stee)
            except:
                results = None
        else:
            results = None

        if outlvl > 0:
            if results is None or results.solver.termination_condition \
                    == TerminationCondition.optimal:
                _log.info("Property initialization for "
                          "{} completed".format(blk.name))
            else:
                _log.warning("Property initialization for "
                             "{} failed".format(blk.name))

        # ---------------------------------------------------------------------
        # Return state to initial conditions
        if state_vars_fixed is False:
            if hold_state is True:
                return flags
            else:
                blk.release_state(flags)

        if outlvl > 0:
            _log.info("Initialization completed for {}".format(blk.name))

    def release_state(blk, flags, outlvl=0):
        '''
        Method to relase state variables fixed during initialization.
        Keyword Arguments:
            flags : dict containing information of which state variables
                    were fixed during initialization, and should now be
                    unfixed. This dict is returned by initialize if
                    hold_state=True.
            outlvl : sets output level of of logging
        '''
        if flags is None:
            return

        # Unfix state variables
        revert_state_vars(blk, flags)

        if outlvl > 0:
            if outlvl > 0:
                _log.info('{} states released.'.format(blk.name))

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
            #bounds=(1e-8, 0.26),
            # bounds=(1e-8, 0.999),
            units=None,
            doc='mass fraction [unitless]')

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
            doc="Mass density [kg/m3]")

        def rule_dens_mass(b):  # density [kg/m3]
            t = b.temperature / pyunits.degK - 273.15
            S = b.mass_frac_comp['NaCl'] * 1000
            params = b.params.dens_mass_params
            A = (params['A_1'] * t + params['A_2']) / params['A_3']
            B = (params['B_1'] * S + params['B_2']) / params['B_3']
            F1 = params['F1']
            F2 = A
            F3 = params['F3_1'] * A ** 2 + params['F3_2']
            F4 = params['F4_1'] * A ** 3 + params['F4_2'] * A
            G1 = params['G1']
            G2 = B
            G3 = params['G3_1'] * B ** 2 + params['G3_2']
            A1 = params['A1_1'] * G1 + params['A1_2'] * G2 + params['A1_3'] * G3
            A2 = params['A2_1'] * G1 + params['A2_2'] * G2 + params['A2_3'] * G3
            A3 = params['A3_1'] * G1 + params['A3_2'] * G2 + params['A3_3'] * G3
            A4 = params['A4_1'] * G1 + params['A4_2'] * G2 + params['A4_3'] * G3
            dens_mass = A1 * F1 + A2 * F2 + A3 * F3 + A4 * F4  # kg/L
            return b.dens_mass == dens_mass * 1e3 * pyunits.kg * pyunits.m**-3
        self.eq_dens_mass = Constraint(rule=rule_dens_mass)

    def _flow_vol(self):
        self.flow_vol = Var(
            initialize=1,
            bounds=(1e-8, 1e8),
            units=pyunits.m**3 / pyunits.s,
            doc="Volumetric flow rate [m3/s]")

        def rule_flow_vol(b):
            return (b.flow_vol == sum(b.flow_mass_comp[j] for j in self.params.component_list)
                    / b.dens_mass)
        self.eq_flow_vol = Constraint(rule=rule_flow_vol)

    def _visc_d(self):
        self.visc_d = Var(
            initialize=1e-3,
            bounds=(1e-8, 1),
            units=pyunits.Pa * pyunits.s,
            doc="Viscosity [Pa-s]")

        def rule_visc_d(b):  # dynamic viscosity [Pa-s]
            t = b.temperature - 273.15 * pyunits.degK
            s = b.mass_frac_comp['NaCl']
            mu_w = (b.params.visc_d_muw_A
                    + (b.params.visc_d_muw_B *
                       (t + b.params.visc_d_muw_C) ** 2
                       - b.params.visc_d_muw_D) ** -1)
            A = b.params.visc_d_A_1 + b.params.visc_d_A_2 * t + b.params.visc_d_A_3 * t ** 2
            B = b.params.visc_d_B_1 + b.params.visc_d_B_2 * t + b.params.visc_d_B_3 * t ** 2
            return b.visc_d == mu_w * (1 + A * s + B * s ** 2)
        self.eq_visc_d = Constraint(rule=rule_visc_d)

    def _diffus(self):
        self.diffus = Var(
            initialize=1e-9,
            bounds=(1e-12, 1e-6),
            units=pyunits.m**2 * pyunits.s**-1,
            doc="Diffusivity [m2/s]")

        def rule_diffus(b):  # diffusivity [m2/s]
            t = b.temperature / pyunits.degK
            params = b.params.diffus_params
            return (b.diffus == (params['A'] * t ** 2 + params['B'] * t
                                 + params['C']) * 1e-9 * pyunits.m**2 * pyunits.s**-1)
        self.eq_diffus = Constraint(rule=rule_diffus)

    def _conc_mass_comp(self):
        self.conc_mass_comp = Var(
            self.params.component_list,
            initialize=10,
            bounds=(1e-6, 1e6),
            units=pyunits.kg * pyunits.m**-3,
            doc="Mass concentration [kg/m3]")

        def rule_conc_mass_comp(b, j):
            return self.conc_mass_comp[j] == \
                   self.dens_mass * self.mass_frac_comp[j]
        self.eq_conc_mass_comp = Constraint(self.params.component_list,rule=rule_conc_mass_comp)


    def _osm_coeff(self):
        self.osm_coeff = Var(
            initialize=1,
            bounds=(1e-8, 10),
            units=None,
            doc="Osmotic coefficient [unitless]")

        def rule_osm_coeff(b):  # osmotic coefficient [-], eq. 49
            s = b.mass_frac_comp['NaCl']  # typo in Sharqawy, s is mass_frac
            t = b.temperature / pyunits.degK - 273.15
            params = b.params.osm_coeff_params
            osm_coeff = (params['1'] + params['2'] * t + params['3'] * t ** 2
                         + params['4'] * t ** 4 + params['5'] * s + params['6'] * s * t
                         + params['7'] * s * t ** 3 + params['8'] * s ** 2
                         + params['9'] * s ** 2 * t
                         + params['10'] * s ** 2 * t ** 2)
            return b.osm_coeff == osm_coeff
        self.eq_osm_coeff = Constraint(rule=rule_osm_coeff)

    def _pressure_osm(self):
        self.pressure_osm = Var(
            initialize=1e6,
            bounds=(1, 1e8),
            units=pyunits.Pa,
            doc="Osmotic pressure [Pa]")

        def rule_pressure_osm(b):  # osmotic pressure [Pa]
            i = 2  # number of ionic species
            return (b.pressure_osm ==
                    i * b.osm_coeff * b.conc_mass_comp['NaCl'] / b.params.mw_comp['NaCl']
                    * Constants.gas_constant * b.temperature)
        self.eq_pressure_osm = Constraint(rule=rule_pressure_osm)

    def _enth_mass(self):
        self.enth_mass = Var(
            initialize=1e6,
            bounds=(1, 1e9),
            units=pyunits.J * pyunits.kg**-1,
            doc="Specific enthalpy [J/kg]")

        def rule_enth_mass(b):  # specific enthalpy [J/kg]
            t = b.temperature / pyunits.degK - 273.15
            S = b.mass_frac_comp['NaCl']
            P = b.pressure
            params = b.params.enth_mass_params
            h_w = params['A1'] + params['A2'] * t + params['A3'] * t ** 2 + params['A4'] * t ** 3
            h_sw = (h_w - (S * (params['B1'] + S) + S * (params['B2'] + S) * t))
            return b.enth_mass == h_sw * pyunits.J * pyunits.kg ** -1
            # h_pterm = b.pressure / b.dens_mass
            # return b.enth_mass == h_sw * pyunits.J * pyunits.kg**-1 + h_pterm
        self.eq_enth_mass = Constraint(rule=rule_enth_mass)

    def _enth_flow(self):
        # enthalpy flow expression for get_enthalpy_flow_terms method

        def rule_enth_flow(b):  # enthalpy flow [J/s]
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
        if p is not 'Liq':
            PropertyPackageError("Property package {} does not support "
                                 "non-liquid phases.".format(self.name))
        return self.flow_mass_comp[j]

    def get_enthalpy_flow_terms(self, p):
        """Create enthalpy flow terms."""
        if p is not 'Liq':
            PropertyPackageError("Property package {} does not support "
                                 "non-liquid phases.".format(self.name))
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
        # temperature, dens_mass, isc_d, diffus, osm_coeff, and enth_mass

        # these variables should have user input
        if iscale.get_scaling_factor(self.flow_mass_comp['H2O']) is None:
            sf = iscale.get_scaling_factor(self.flow_mass_comp['H2O'], default=1e0, warning=True)
            iscale.set_scaling_factor(self.flow_mass_comp['H2O'], sf)

        if iscale.get_scaling_factor(self.flow_mass_comp['NaCl']) is None:
            sf = iscale.get_scaling_factor(self.flow_mass_comp['NaCl'], default=1e2, warning=True)
            iscale.set_scaling_factor(self.flow_mass_comp['NaCl'], sf)

        if iscale.get_scaling_factor(self.flow_mass_comp['H2O']) is None:
            sf = iscale.get_scaling_factor(self.pressure, default=1e-6, warning=True)
            iscale.set_scaling_factor(self.pressure, sf)

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
