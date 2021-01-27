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
Initial property package for seawater system
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
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.exceptions import PropertyPackageError
import idaes.core.util.scaling as iscale

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("SeawaterParameterBlock")
class SeawaterParameterData(PhysicalParameterBlock):
    CONFIG = PhysicalParameterBlock.CONFIG()

    def build(self):
        '''
        Callable method for Block construction.
        '''
        super(SeawaterParameterData, self).build()

        self._state_block_class = SeawaterStateBlock

        # components
        self.H2O = Solvent()
        self.TDS = Solute()

        # phases
        self.Liq = LiquidPhase()

        # reference
        # this package is developed from Sharqawy et al. (2010) http://dx.doi.org/10.5004/dwt.2010.1079

        # parameters
        # molecular weight
        mw_comp_data = {'H2O': 18.01528E-3,
                        'TDS': 58.44E-3}  # TODO: Confirm how Sharqawy converts TDS to moles
        self.mw_comp = Param(self.component_list,
                             mutable=False,
                             initialize=mw_comp_data,
                             units=pyunits.kg/pyunits.mol,
                             doc="Molecular weight")

        # mass density parameters, eq. 8 in Sharqawy
        self.dens_mass_param_A1 = Var(
            within=Reals, initialize=9.999e2, units=pyunits.kg/pyunits.m**3,
            doc='Mass density parameter A1')
        self.dens_mass_param_A2 = Var(
            within=Reals, initialize=2.034e-2, units=(pyunits.kg/pyunits.m**3) * pyunits.K**-1,
            doc='Mass density parameter A2')
        self.dens_mass_param_A3 = Var(
            within=Reals, initialize=-6.162e-3, units=(pyunits.kg/pyunits.m**3) * pyunits.K**-2,
            doc='Mass density parameter A3')
        self.dens_mass_param_A4 = Var(
            within=Reals, initialize=2.261e-5, units=(pyunits.kg/pyunits.m**3) * pyunits.K**-3,
            doc='Mass density parameter A4')
        self.dens_mass_param_A5 = Var(
            within=Reals, initialize=-4.657e-8, units=(pyunits.kg/pyunits.m**3) * pyunits.K**-4,
            doc='Mass density parameter A5')
        self.dens_mass_param_B1 = Var(
            within=Reals, initialize=8.020e2, units=pyunits.kg/pyunits.m ** 3,
            doc='Mass density parameter B1')
        self.dens_mass_param_B2 = Var(
            within=Reals, initialize=-2.001, units=(pyunits.kg/pyunits.m**3) * pyunits.K**-1,
            doc='Mass density parameter B2')
        self.dens_mass_param_B3 = Var(
            within=Reals, initialize=1.677e-2, units=(pyunits.kg/pyunits.m**3) * pyunits.K**-2,
            doc='Mass density parameter B3')
        self.dens_mass_param_B4 = Var(
            within=Reals, initialize=-3.060e-5, units=(pyunits.kg/pyunits.m**3) * pyunits.K**-3,
            doc='Mass density parameter B4')
        self.dens_mass_param_B5 = Var(
            within=Reals, initialize=-1.613e-5, units=(pyunits.kg/pyunits.m**3) * pyunits.K**-2,
            doc='Mass density parameter B5')
        self.dens_mass_param_dict = {
            'A1': self.dens_mass_param_A1,
            'A2': self.dens_mass_param_A2,
            'A3': self.dens_mass_param_A3,
            'A4': self.dens_mass_param_A4,
            'A5': self.dens_mass_param_A5,
            'B1': self.dens_mass_param_B1,
            'B2': self.dens_mass_param_B2,
            'B3': self.dens_mass_param_B3,
            'B4': self.dens_mass_param_B4,
            'B5': self.dens_mass_param_B5}

        # dynamic viscosity parameters, eq. 22 and 23 in Sharqawy
        self.visc_d_param_muw_A = Var(
            within=Reals, initialize=4.2844e-5, units=pyunits.Pa*pyunits.s,
            doc='Dynamic viscosity parameter A for pure water')
        self.visc_d_param_muw_B = Var(
            within=Reals, initialize=0.157, units=pyunits.degK**-2*pyunits.Pa**-1*pyunits.s**-1,
            doc='Dynamic viscosity parameter B for pure water')
        self.visc_d_param_muw_C = Var(
            within=Reals, initialize=64.993, units=pyunits.K,
            doc='Dynamic viscosity parameter C for pure water')
        self.visc_d_param_muw_D = Var(
            within=Reals, initialize=91.296, units=pyunits.Pa**-1*pyunits.s**-1,
            doc='Dynamic viscosity parameter D for pure water')
        self.visc_d_param_A_1 = Var(
            within=Reals, initialize=1.541, units=None,
            doc='Dynamic viscosity parameter 1 for term A')
        self.visc_d_param_A_2 = Var(
            within=Reals, initialize=1.998e-2, units=pyunits.K**-1,
            doc='Dynamic viscosity parameter 2 for term A')
        self.visc_d_param_A_3 = Var(
            within=Reals, initialize=-9.52e-5, units=pyunits.K**-2,
            doc='Dynamic viscosity parameter 3 for term A')
        self.visc_d_param_B_1 = Var(
            within=Reals, initialize=7.974, units=None,
            doc='Dynamic viscosity parameter 1 for term B')
        self.visc_d_param_B_2 = Var(
            within=Reals, initialize=-7.561e-2, units=pyunits.K**-1,
            doc='Dynamic viscosity parameter 2 for term B')
        self.visc_d_param_B_3 = Var(
            within=Reals, initialize=4.724e-4, units=pyunits.K**-2,
            doc='Dynamic viscosity parameter 3 for term B')
        self.visc_d_param_dict = {
            'muw_A': self.visc_d_param_muw_A,
            'muw_B': self.visc_d_param_muw_B,
            'muw_C': self.visc_d_param_muw_C,
            'muw_D': self.visc_d_param_muw_D,
            'A_1': self.visc_d_param_A_1,
            'A_2': self.visc_d_param_A_2,
            'A_3': self.visc_d_param_A_3,
            'B_1': self.visc_d_param_B_1,
            'B_2': self.visc_d_param_B_2,
            'B_3': self.visc_d_param_B_3}

        # osmotic coefficient parameters, eq. 49 in Sharqawy
        self.osm_coeff_param_1 = Var(
            within=Reals, initialize=8.9453e-1, units=pyunits.dimensionless,
            doc='Osmotic coefficient parameter 1')
        self.osm_coeff_param_2 = Var(
            within=Reals, initialize=4.1561e-4, units=pyunits.K**-1,
            doc='Osmotic coefficient parameter 2')
        self.osm_coeff_param_3 = Var(
            within=Reals, initialize=-4.6262e-6, units=pyunits.K**-2,
            doc='Osmotic coefficient parameter 3')
        self.osm_coeff_param_4 = Var(
            within=Reals, initialize=2.2211e-11, units=pyunits.K**-4,
            doc='Osmotic coefficient parameter 4')
        self.osm_coeff_param_5 = Var(
            within=Reals, initialize=-1.1445e-1, units=pyunits.dimensionless,
            doc='Osmotic coefficient parameter 5')
        self.osm_coeff_param_6 = Var(
            within=Reals, initialize=-1.4783e-3, units=pyunits.K**-1,
            doc='Osmotic coefficient parameter 6')
        self.osm_coeff_param_7 = Var(
            within=Reals, initialize=-1.3526e-8, units=pyunits.K**-3,
            doc='Osmotic coefficient parameter 7')
        self.osm_coeff_param_8 = Var(
            within=Reals, initialize=7.0132, units=pyunits.dimensionless,
            doc='Osmotic coefficient parameter 8')
        self.osm_coeff_param_9 = Var(
            within=Reals, initialize=5.696e-2, units=pyunits.K**-1,
            doc='Osmotic coefficient parameter 9')
        self.osm_coeff_param_10 = Var(
            within=Reals, initialize=-2.8624e-4, units=pyunits.K**-2,
            doc='Osmotic coefficient parameter 10')
        self.osm_coeff_param_dict = {
            '1': self.osm_coeff_param_1,
            '2': self.osm_coeff_param_2,
            '3': self.osm_coeff_param_3,
            '4': self.osm_coeff_param_4,
            '5': self.osm_coeff_param_5,
            '6': self.osm_coeff_param_6,
            '7': self.osm_coeff_param_7,
            '8': self.osm_coeff_param_8,
            '9': self.osm_coeff_param_9,
            '10': self.osm_coeff_param_10}

        # specific enthalpy parameters, eq. 55 and 43 in Sharqawy
        enth_mass_dic = {'A1': 124.790, 'A2': 4203.075, 'A3': -0.552, 'A4': 0.004,
                         'B1': 27062.623, 'B2': 4835.675}
        self.enth_mass_params = Var(
            [*enth_mass_dic],
            within=Reals, initialize=enth_mass_dic, units=None,
            doc='Specific enthalpy parameters')
        self.enth_mass_param_A1 = Var(
            within=Reals, initialize=124.790, units=pyunits.J/pyunits.kg,
            doc='Specific enthalpy parameter A1')
        self.enth_mass_param_A2 = Var(
            within=Reals, initialize=4203.075, units=(pyunits.J/pyunits.kg) * pyunits.K**-1,
            doc='Specific enthalpy parameter A2')
        self.enth_mass_param_A3 = Var(
            within=Reals, initialize=-0.552, units=(pyunits.J/pyunits.kg) * pyunits.K**-2,
            doc='Specific enthalpy parameter A3')
        self.enth_mass_param_A4 = Var(
            within=Reals, initialize=0.004, units=(pyunits.J/pyunits.kg) * pyunits.K**-3,
            doc='Specific enthalpy parameter A4')
        self.enth_mass_param_B1 = Var(
            within=Reals, initialize=27062.623, units=pyunits.dimensionless,
            doc='Specific enthalpy parameter B1')
        self.enth_mass_param_B2 = Var(
            within=Reals, initialize=4835.675, units=pyunits.dimensionless,
            doc='Specific enthalpy parameter B2')
        self.enth_mass_param_dict = {
            'A1': self.enth_mass_param_A1,
            'A2': self.enth_mass_param_A2,
            'A3': self.enth_mass_param_A3,
            'A4': self.enth_mass_param_A4,
            'B1': self.enth_mass_param_B1,
            'B2': self.enth_mass_param_B2}

        # traditional parameters are the only Vars currently on the block and should be fixed
        for v in self.component_objects(Var):
            v.fix()

        # ---default scaling---
        self.set_default_scaling('temperature', 1e-2)
        self.set_default_scaling('dens_mass', 1e-3)
        self.set_default_scaling('visc_d', 1e3)
        self.set_default_scaling('osm_coeff', 1e0)
        self.set_default_scaling('enth_mass', 1e-5)
        # TODO: require users to provide the following scaling factors
        self.set_default_scaling('flow_mass_comp', 1e0, 'H2O')
        self.set_default_scaling('flow_mass_comp', 1e2, 'TDS')
        self.set_default_scaling('pressure', 1e-6)

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


class _SeawaterStateBlock(StateBlock):
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

        # Fix state variables if not already fixed
        if state_vars_fixed is False:
            flags = fix_state_vars(self, state_args)
            # Check when the state vars are fixed already result in dof 0
            for k in self.keys():
                dof = degrees_of_freedom(self[k])
                if dof != 0:
                    init_log.error("Degrees of freedom not 0, ({})".format(dof))
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

@declare_process_block_class("SeawaterStateBlock",
                             block_class=_SeawaterStateBlock)
class SeawaterStateBlockData(StateBlockData):
    """A seawater property package."""
    def build(self):
        """Callable method for Block construction."""
        super(SeawaterStateBlockData, self).build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        # Add state variables
        self.flow_mass_comp = Var(
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
            units=pyunits.K,
            doc='Temperature')

        self.pressure = Var(
            initialize=101325,
            bounds=(1e5, 5e7),
            domain=NonNegativeReals,
            units=pyunits.Pa,
            doc='Pressure')

    # -----------------------------------------------------------------------------
    # Property Methods
    def _mass_frac_comp(self):
        self.mass_frac_comp = Var(
            self.params.component_list,
            initialize=0.1,
            bounds=(1e-8, 1),
            units=None,
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

        def rule_dens_mass(b):  # density, eq. 8 in Sharqawy
            t = b.temperature - 273.15 * pyunits.K
            s = b.mass_frac_comp['TDS']
            params = b.params.dens_mass_param_dict
            dens_mass = (params['A1'] + params['A2'] * t + params['A3'] * t ** 2
                         + params['A4'] * t ** 3 + params['A5'] * t ** 4
                         + params['B1'] * s + params['B2'] * s * t + params['B3'] * s * t ** 2
                         + params['B4'] * s * t ** 3 + params['B5'] * s ** 2 * t ** 2)
            return b.dens_mass == dens_mass
        self.eq_dens_mass = Constraint(rule=rule_dens_mass)

    def _flow_vol(self):
        self.flow_vol = Var(
            initialize=1,
            bounds=(1e-8, 1e8),
            units=pyunits.m**3 / pyunits.s,
            doc="Volumetric flow rate")

        def rule_flow_vol(b):
            return (b.flow_vol == sum(b.flow_mass_comp[j] for j in self.params.component_list)
                    / b.dens_mass)
        self.eq_flow_vol = Constraint(rule=rule_flow_vol)

    def _visc_d(self):
        self.visc_d = Var(
            initialize=1e-3,
            bounds=(1e-8, 1),
            units=pyunits.Pa * pyunits.s,
            doc="Viscosity")

        def rule_visc_d(b):  # dynamic viscosity, eq. 22 and 23 in Sharqawy
            t = b.temperature - 273.15 * pyunits.K  # temperature in degC, but pyunits are K
            s = b.mass_frac_comp['TDS']
            params = b.params.visc_d_param_dict
            mu_w = (params['muw_A']
                    + (params['muw_B'] *
                       (t + params['muw_C']) ** 2
                       - params['muw_D']) ** -1)
            A = (params['A_1'] + params['A_2'] * t + params['A_3'] * t ** 2)
            B = (params['B_1'] + params['B_2'] * t + params['B_3'] * t ** 2)
            return b.visc_d == mu_w * (1 + A * s + B * s ** 2)
        self.eq_visc_d = Constraint(rule=rule_visc_d)

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
        self.eq_conc_mass_comp = Constraint(self.params.component_list,rule=rule_conc_mass_comp)


    def _osm_coeff(self):
        self.osm_coeff = Var(
            initialize=1,
            bounds=(1e-8, 10),
            units=None,
            doc="Osmotic coefficient")

        def rule_osm_coeff(b):  # osmotic coefficient, eq. 49 in Sharqawy
            s = b.mass_frac_comp['TDS']
            t = b.temperature - 273.15 * pyunits.K  # temperature in degC, but pyunits are still K
            params = b.params.osm_coeff_param_dict
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
            doc="Osmotic pressure")

        def rule_pressure_osm(b):  # osmotic pressure, based on molality and assumes TDS is NaCl
            i = 2  # number of ionic species
            molality = b.mass_frac_comp['TDS'] / (1 - b.mass_frac_comp['TDS']) / b.params.mw_comp['TDS']
            rhow = 1000 * pyunits.kg/pyunits.m**3  # TODO: could make this variable based on temperature
            return (b.pressure_osm ==
                    i * b.osm_coeff * molality * rhow * Constants.gas_constant * b.temperature)
        self.eq_pressure_osm = Constraint(rule=rule_pressure_osm)

    def _enth_mass(self):
        self.enth_mass = Var(
            initialize=1e6,
            bounds=(1, 1e9),
            units=pyunits.J * pyunits.kg**-1,
            doc="Specific enthalpy")

        def rule_enth_mass(b):  # specific enthalpy, eq. 55 and 43 in Sharqawy
            t = b.temperature - 273.15 * pyunits.K  # temperature in degC, but pyunits in K
            S = b.mass_frac_comp['TDS']
            params = b.params.enth_mass_param_dict
            h_w = params['A1'] + params['A2'] * t + params['A3'] * t ** 2 + params['A4'] * t ** 3
            h_sw = (h_w -
                    (S * (params['B1'] + S)
                     + S * (params['B2'] + S) * t/pyunits.K)
                    * pyunits.J/pyunits.kg)  # must make relationship dimensionless and add units at end
            return b.enth_mass == h_sw
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
        # temperature, dens_mass, isc_d, diffus, osm_coeff, and enth_mass

        # these variables should have user input
        if iscale.get_scaling_factor(self.flow_mass_comp['H2O']) is None:
            sf = iscale.get_scaling_factor(self.flow_mass_comp['H2O'], default=1e0, warning=True)
            iscale.set_scaling_factor(self.flow_mass_comp['H2O'], sf)

        if iscale.get_scaling_factor(self.flow_mass_comp['TDS']) is None:
            sf = iscale.get_scaling_factor(self.flow_mass_comp['TDS'], default=1e2, warning=True)
            iscale.set_scaling_factor(self.flow_mass_comp['TDS'], sf)

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
                    elif j == 'TDS':
                        iscale.set_scaling_factor(
                            self.mass_frac_comp[j],
                            iscale.get_scaling_factor(self.flow_mass_comp['TDS'])
                            / iscale.get_scaling_factor(self.flow_mass_comp['H2O']))

        if self.is_property_constructed('conc_mass_comp'):
            for j in self.params.component_list:
                sf_dens = iscale.get_scaling_factor(self.dens_mass)
                if iscale.get_scaling_factor(self.conc_mass_comp[j]) is None:
                    if j == 'H2O':
                        # solvents typically have a mass fraction between 0.5-1
                        iscale.set_scaling_factor(self.conc_mass_comp[j], sf_dens)
                    elif j == 'TDS':
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