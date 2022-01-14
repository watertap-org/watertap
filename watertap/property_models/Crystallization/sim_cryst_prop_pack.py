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
Initial crystallization property package for H2O-NaCl system
"""

# Import Python libraries
import idaes.logger as idaeslog

# Import Pyomo libraries
from pyomo.environ import Constraint, Expression, Reals, NonNegativeReals, \
    Var, Param, exp, Suffix, value, check_optimal_termination
from pyomo.environ import units as pyunits

# Import IDAES cores
from idaes.core import (declare_process_block_class,
                        MaterialFlowBasis,
                        PhysicalParameterBlock,
                        StateBlockData,
                        StateBlock,
                        MaterialBalanceType,
                        EnergyBalanceType)
from idaes.core.components import Component, Solute, Solvent
from idaes.core.phases import LiquidPhase, VaporPhase, SolidPhase
from idaes.core.util.constants import Constants
from idaes.core.util.initialization import (fix_state_vars,
                                            revert_state_vars,
                                            solve_indexed_blocks)
from idaes.core.util.misc import add_object_reference, extract_data
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom, \
    number_unfixed_variables
from idaes.core.util.exceptions import ConfigurationError, PropertyPackageError
import idaes.core.util.scaling as iscale


@declare_process_block_class("NaClParameterBlock")
class NaClParameterData(PhysicalParameterBlock):

    def build(self):
        super().build()
        self._state_block_class = NaClStateBlock

        # Component
        self.H2O = Solvent()
        self.NaCl = Solute()

        # Phases
        self.Liq = LiquidPhase()
        self.Vap = VaporPhase()
        self.Sol = SolidPhase()

        ''' 
       References
        This package was developed from the following references:

        - K.G.Nayar, M.H.Sharqawy, L.D.Banchik, and J.H.Lienhard V, "Thermophysical properties of seawater: A review and
        new correlations that include pressure dependence,"Desalination, Vol.390, pp.1 - 24, 2016.
        doi: 10.1016/j.desal.2016.02.024(preprint)

        - Mostafa H.Sharqawy, John H.Lienhard V, and Syed M.Zubair, "Thermophysical properties of seawater: A review of 
        existing correlations and data,"Desalination and Water Treatment, Vol.16, pp.354 - 380, April 2010.
        (2017 corrections provided at http://web.mit.edu/seawater)

        - Laliberté, M. & Cooper, W. E. Model for Calculating the Density of Aqueous Electrolyte Solutions  
        Journal of Chemical & Engineering Data, American Chemical Society (ACS), 2004, 49, 1141-1151.

        - Laliberté, M. A Model for Calculating the Heat Capacity of Aqueous Solutions, with Updated Density and Viscosity Data 
        Journal of Chemical & Engineering Data, American Chemical Society (ACS), 2009, 54, 1725-1760 
        Liquid NaCl heat capacity and density parameters available https://pubs.acs.org/doi/10.1021/je8008123

        Chase, M.W., Jr., NIST-JANAF Themochemical Tables, Fourth Edition, J. Phys. Chem. Ref. Data, Monograph 9, 
        1998, 1-1951.

        - Tavare, N. S. Industrial Crystallization, Springer US, 2013.        
        '''

        ## TO-DO ##
        # 1. Add heat of crystalization
        # 2. Add enthalpy functions
        # 3. Add vapour pressure 
        # 4. Add scaling rules

        # Unit definitions
        dens_units = pyunits.kg/pyunits.m**3
        t_inv_units = pyunits.K**-1
        enth_mass_units = pyunits.J/pyunits.kg
        cp_units = pyunits.J/(pyunits.kg*pyunits.K)
        cp_units_2 = pyunits.kJ/(pyunits.kg*pyunits.K)

        # Solubility parameters from surrogate model
        self.sol_param_A1 = Param(initialize= 526.706475, units = pyunits.g / pyunits.L, doc=' Solubility parameter A1 [g/L] for NaCl surrogate')
        self.sol_param_A2 = Param(initialize= -1.326952, units = (pyunits.g / pyunits.L) * pyunits.K ** -1, doc=' Solubility parameter A2 [g/L] for NaCl surrogate')
        self.sol_param_A3 = Param(initialize= 0.002574, units = (pyunits.g / pyunits.L) * pyunits.K ** -2, doc=' Solubility parameter A3 [g/L] for NaCl surrogate')
        
        # Vapour mass density parameters (approximating using ideal gas)
        self.dens_mass_param_mw = Var(within=Reals, initialize=18.01528e-3, units=pyunits.kg/pyunits.mol,
            doc='Mass density parameter molecular weight')
        self.dens_mass_param_R = Var(within=Reals, initialize=8.31462618, units=pyunits.J/pyunits.mol/pyunits.K,
            doc='Mass density parameter universal gas constant')

        # Mass density and heat capacity values for NaCl crystals in solid phase: fixed for now at Tavare values - may not be accurate?
        self.dens_mass_param_NaCl = Param(initialize=2115, units=pyunits.kg/pyunits.m**3, doc='NaCl crystal density')
        self.cp_param_NaCl = Param(initialize=877, units=pyunits.J/(pyunits.kg*pyunits.K), doc='NaCl crystal density')

        # Mass density parameters for pure NaCl liquid based on Eq. 9 in Laliberte and Cooper (2004).
        self.dens_mass_param_NaCl_liq_C0 = Var(within=Reals, initialize=-0.00433, units=dens_units, doc='Mass density parameter C0 for liquid NaCl')
        self.dens_mass_param_NaCl_liq_C1 = Var(within=Reals, initialize=0.06471, units=dens_units, doc='Mass density parameter C1 for liquid NaCl')
        self.dens_mass_param_NaCl_liq_C2 = Var(within=Reals, initialize=1.01660, units=pyunits.dimensionless, doc='Mass density parameter C2 for liquid NaCl')
        self.dens_mass_param_NaCl_liq_C3 = Var(within=Reals, initialize=0.014624, units=t_inv_units, doc='Mass density parameter C3 for liquid NaCl')
        self.dens_mass_param_NaCl_liq_C4 = Var(within=Reals, initialize=3315.6, units=pyunits.K, doc='Mass density parameter C4 for liquid NaCl')

        # Mass density parameters for solvent in liquid phase, eq. 8 in Sharqawy et al. (2010)
        self.dens_mass_param_A1 = Var(within=Reals, initialize=9.999e2, units=dens_units, doc='Mass density parameter A1')
        self.dens_mass_param_A2 = Var(within=Reals, initialize=2.034e-2, units=dens_units * t_inv_units, doc='Mass density parameter A2')
        self.dens_mass_param_A3 = Var(within=Reals, initialize=-6.162e-3, units=dens_units * t_inv_units**2, doc='Mass density parameter A3')
        self.dens_mass_param_A4 = Var(within=Reals, initialize=2.261e-5, units=dens_units * t_inv_units**3, doc='Mass density parameter A4')
        self.dens_mass_param_A5 = Var(within=Reals, initialize=-4.657e-8, units=dens_units * t_inv_units**4, doc='Mass density parameter A5')

        # Latent heat of evaporation of pure water: Parameters from Sharqawy et al. (2010), eq. 54
        self.dh_vap_w_param_0 = Var(within=Reals, initialize=2.501e6, units=enth_mass_units, doc='Latent heat of pure water parameter 0')
        self.dh_vap_w_param_1 = Var(within=Reals, initialize=-2.369e3, units=enth_mass_units * t_inv_units**1, doc='Latent heat of pure water parameter 1')
        self.dh_vap_w_param_2 = Var(within=Reals, initialize=2.678e-1, units=enth_mass_units * t_inv_units**2, doc='Latent heat of pure water parameter 2')
        self.dh_vap_w_param_3 = Var(within=Reals, initialize=-8.103e-3, units=enth_mass_units * t_inv_units**3, doc='Latent heat of pure water parameter 3')
        self.dh_vap_w_param_4 = Var(within=Reals, initialize=-2.079e-5, units=enth_mass_units * t_inv_units**4, doc='Latent heat of pure water parameter 4')

        # Specific heat parameters for Cp vapor from NIST Webbook - Chase, M.W., Jr., NIST-JANAF Themochemical Tables
        self.cp_vap_param_A = Var(within=Reals, initialize=30.09200 / 18.01528e-3, units=cp_units, doc='Specific heat of water vapor parameter A')
        self.cp_vap_param_B = Var(within=Reals, initialize=6.832514 / 18.01528e-3, units=cp_units * t_inv_units, doc='Specific heat of water vapor parameter B')
        self.cp_vap_param_C = Var(within=Reals, initialize=6.793435 / 18.01528e-3, units=cp_units * t_inv_units**2, doc='Specific heat of water vapor parameter C')
        self.cp_vap_param_D = Var(within=Reals, initialize=-2.534480 / 18.01528e-3, units=cp_units * t_inv_units**3, doc='Specific heat of water vapor parameter D')
        self.cp_vap_param_E = Var(within=Reals, initialize=0.082139 / 18.01528e-3, units=cp_units * t_inv_units**-2, doc='Specific heat of water vapor parameter E')

        # Specific heat parameters for pure water from eq (9) in Sharqawy et al. (2010)
        self.cp_phase_param_A1 = Var(within=Reals, initialize=5.328, units=cp_units, doc='Specific heat of seawater parameter A1')
        self.cp_phase_param_B1 = Var(within=Reals, initialize=-6.913e-3, units=cp_units * t_inv_units, doc='Specific heat of seawater parameter B1')
        self.cp_phase_param_C1 = Var(within=Reals, initialize=9.6e-6, units=cp_units * t_inv_units**2, doc='Specific heat of seawater parameter C1')
        self.cp_phase_param_D1 = Var(within=Reals, initialize=2.5e-9, units=cp_units * t_inv_units**3, doc='Specific heat of seawater parameter D1')

        # Specific heat parameters for liquid NaCl for eqs. (11) & (12) in Laliberte (2009).
        self.cp_param_NaCl_liq_A1 = Var(within=Reals, initialize=-0.06936, units=cp_units_2, doc='Specific heat parameter A1 for liquid NaCl')
        self.cp_param_NaCl_liq_A2 = Var(within=Reals, initialize=-0.07821, units=t_inv_units, doc='Specific heat parameter A2 for liquid NaCl')
        self.cp_param_NaCl_liq_A3 = Var(within=Reals, initialize=3.8480, units=pyunits.dimensionless, doc='Specific heat parameter A3 for liquid NaCl')
        self.cp_param_NaCl_liq_A4 = Var(within=Reals, initialize=-11.2762, units=pyunits.dimensionless, doc='Specific heat parameter A4 for liquid NaCl')
        self.cp_param_NaCl_liq_A5 = Var(within=Reals, initialize=8.7319, units=cp_units_2, doc='Specific heat parameter A5 for liquid NaCl')
        self.cp_param_NaCl_liq_A6 = Var(within=Reals, initialize=1.8125, units=pyunits.dimensionless, doc='Specific heat parameter A6 for liquid NaCl')

        for v in self.component_objects(Var):
            v.fix()


    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units({'time': pyunits.s,
                               'length': pyunits.m,
                               'mass': pyunits.kg,
                               'amount': pyunits.mol,
                               'temperature': pyunits.K})


        obj.add_properties(
            {'flow_mass_phase_comp': {'method': None},
             'temperature': {'method': None},
             'pressure': {'method': None},
             'solubility_mass_phase_comp': {'method': '_solubility_mass_phase_comp'},
             'solubility_mass_frac_phase_comp': {'method': '_solubility_mass_frac_phase_comp'},
             'mass_frac_phase_comp': {'method': '_mass_frac_phase_comp'},
             'dens_mass_solvent': {'method': '_dens_mass_solvent'},
             'dens_mass_solute': {'method': '_dens_mass_solute'},
             'dens_mass_phase': {'method': '_dens_mass_phase'},
             'dh_vap_solvent': {'method': '_dh_vap_solvent'},
             'cp_solvent': {'method': '_cp_solvent'},
             'cp_solute': {'method': '_cp_solute'},
             'cp_phase': {'method': '_cp_phase'},
             'flow_vol_phase': {'method': '_flow_vol_phase'},
             'flow_vol': {'method': '_flow_vol'},
             # 'conc_mass_phase_comp': {'method': '_conc_mass_phase_comp'},
             # 'flow_mol_phase_comp': {'method': '_flow_mol_phase_comp'},
             # 'mole_frac_phase_comp': {'method': '_mole_frac_phase_comp'},
             # 'enth_mass_phase': {'method': '_enth_mass_phase'},
             # 'enth_flow': {'method': '_enth_flow'}
                })



class _NaClStateBlock(StateBlock):
    # Will build up later
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



@declare_process_block_class("NaClStateBlock", block_class=_NaClStateBlock)
class NaClStateBlockData(StateBlockData):
    def build(self):
        """Callable method for Block construction."""
        super(NaClStateBlockData, self).build()
        self._make_state_vars()

    def _make_state_vars(self):
        # Create state variables
        self.pressure = Var(domain=NonNegativeReals,
                            initialize=101325,
                            units = pyunits.Pa,
                            doc='State pressure [Pa]')

        self.temperature = Var(domain=NonNegativeReals,
                            initialize=298.15,
                            bounds = (273.15, 373.15),
                            units = pyunits.degK,
                            doc='State temperature [K]')

        self.flow_mass_phase_comp = Var(self.params.phase_list, self.params.component_list,
                            initialize={('Liq', 'H2O'): 0.965, ('Liq', 'NaCl'): 0.035},
                            bounds=(1e-8, None),
                            domain=NonNegativeReals,
                            units=pyunits.kg/pyunits.s,
                            doc='Mass flow rate')

    # Property Methods

    # 1 Mass fraction: From NaCl property package 
    def _mass_frac_phase_comp(self):
        self.mass_frac_phase_comp = Var(self.params.phase_list,
                                    self.params.component_list,
                                    initialize={('Liq', 'H2O'): 0.965, ('Liq', 'NaCl'): 0.035},
                                    bounds=(1e-6, None),  # upper bound set to None because of stability benefits
                                    units=pyunits.dimensionless,
                                    doc='Mass fraction')

        def rule_mass_frac_phase_comp(b, p, j):
            if p =='Liq':
                return (b.mass_frac_phase_comp[p, j] == b.flow_mass_phase_comp[p, j] /
                                        sum(b.flow_mass_phase_comp[p, j] for j in self.params.component_list))
            else:
                return Constraint.Skip ##MUST BE FIXED!##

        self.eq_mass_frac_phase_comp = Constraint(self.params.phase_list, self.params.component_list, rule=rule_mass_frac_phase_comp)


    # 2. Solubility: Based on regressed/surrogate model
    def _solubility_mass_phase_comp(self):
        self.solubility_mass_phase_comp = Var(['Liq'], ['NaCl'],
                                domain = NonNegativeReals,
                                bounds = (340, 400),
                                initialize = 356.5,
                                units = pyunits.g / pyunits.L,
                                doc="solubility of NaCl in water, g/L")

        def rule_solubility_mass_phase_comp(b, j):
            return (b.solubility_mass_phase_comp['Liq', j] == b.params.sol_param_A1 
                                        + b.params.sol_param_A2 * b.temperature 
                                        + b.params.sol_param_A3 * b.temperature ** 2
                                        )

        self.eq_solubility_mass_phase_comp = Constraint(['NaCl'], rule=rule_solubility_mass_phase_comp)


    # 3. Solubility as mass fraction
    def _solubility_mass_frac_phase_comp(self):
        self.solubility_mass_frac_phase_comp = Var(['Liq'], ['NaCl'],
                                domain = NonNegativeReals,
                                bounds = (0, 1),
                                initialize = 0.5,
                                units = pyunits.dimensionless,
                                doc="solubility (as mass fraction) of NaCl in water")

        def rule_solubility_mass_frac_phase_comp(b, j):
            return (b.solubility_mass_frac_phase_comp['Liq', j] 
                == b.solubility_mass_phase_comp['Liq', j] / (b.solubility_mass_phase_comp['Liq', j] + b.dens_mass_solvent['Liq'])
                                        )

        self.eq_solubility_mass_frac_phase_comp = Constraint(['NaCl'], rule=rule_solubility_mass_frac_phase_comp)


    # 4. Density of solvent (pure water in liquid and vapour phases)
    def _dens_mass_solvent(self):
        self.dens_mass_solvent = Var(['Liq', 'Vap'],
            initialize=1e3,
            bounds=(1e-4, 1e6),
            units=pyunits.kg*pyunits.m**-3,
            doc="Mass density of pure water")

        def rule_dens_mass_solvent(b, p):  
            if p == 'Liq': # density, eq. 8 in Sharqawy
                t = b.temperature - 273.15*pyunits.K
                dens_mass_w = (b.params.dens_mass_param_A1
                             + b.params.dens_mass_param_A2 * t
                             + b.params.dens_mass_param_A3 * t**2
                             + b.params.dens_mass_param_A4 * t**3
                             + b.params.dens_mass_param_A5 * t**4)
                return b.dens_mass_solvent[p] == dens_mass_w
            elif p =='Vap':
                return b.dens_mass_solvent[p] ==  (b.params.dens_mass_param_mw * b.pressure) / (b.params.dens_mass_param_R * b.temperature)

        self.eq_dens_mass_solvent = Constraint(['Liq', 'Vap'], rule=rule_dens_mass_solvent)


    # 5. Density of NaCl crystals and liquid
    def _dens_mass_solute(self):
        self.dens_mass_solute = Var(['Sol', 'Liq'],
            initialize=1e3,
            bounds=(1e-4, 1e6),
            units=pyunits.kg*pyunits.m**-3,
            doc="Mass density of solid NaCl crystals")

        def rule_dens_mass_solute(b, p):  
            if p == 'Sol':
                return b.dens_mass_solute[p] == b.params.dens_mass_param_NaCl
            elif p == 'Liq': # eq. 9 of Laliberte paper
                t = b.temperature - 273.15 * pyunits.K
                v_app = (b.mass_frac_phase_comp['Liq', 'NaCl'] + b.params.dens_mass_param_NaCl_liq_C2 + (b.params.dens_mass_param_NaCl_liq_C3 * t)) \
                    / ((b.mass_frac_phase_comp['Liq', 'NaCl'] * b.params.dens_mass_param_NaCl_liq_C0) + b.params.dens_mass_param_NaCl_liq_C1)
                v_app = v_app / exp(0.000001 * pyunits.K**-2 * (t + b.params.dens_mass_param_NaCl_liq_C4) ** 2) 
                return b.dens_mass_solute[p] == 1 / v_app

        self.eq_dens_mass_solute = Constraint(['Sol', 'Liq'], rule=rule_dens_mass_solute)


    # 6. Density of liquid solution (Water + NaCl)
    def _dens_mass_phase(self):
        self.dens_mass_phase = Var(
            ['Liq'],
            initialize=1e3,
            bounds=(5e2, 1e4),
            units=pyunits.kg * pyunits.m ** -3,
            doc="Mass density of liquid NaCl solution")

        def rule_dens_mass_phase(b):  # density, eq. 6 of Laliberte paper
            return (b.dens_mass_phase['Liq'] ==
                     1 / ((b.mass_frac_phase_comp['Liq', 'NaCl'] / b.dens_mass_solute['Liq'])  +  (b.mass_frac_phase_comp['Liq', 'H2O'] / b.dens_mass_solvent['Liq']))
                     )
        self.eq_dens_mass_phase = Constraint(rule=rule_dens_mass_phase)


    # 7. Latent heat of vapourization of pure water
    def _dh_vap_solvent(self):
        self.dh_vap_solvent = Var(
            initialize=2.4e3,
            bounds=(1, 1e9),
            units=pyunits.J/pyunits.kg,
            doc="Latent heat of vaporization of pure water")

        def rule_dh_vap_solvent(b):
            t = b.temperature - 273.15 * pyunits.K
            return b.dh_vap_solvent == b.params.dh_vap_w_param_0 + b.params.dh_vap_w_param_1 * t + b.params.dh_vap_w_param_2 * t**2 \
                       + b.params.dh_vap_w_param_3 * t**3 + b.params.dh_vap_w_param_4 * t**4

        self.eq_dh_vap_solvent = Constraint(rule=rule_dh_vap_solvent)


    # 8. Heat capacity of solvent (pure water in liquid and vapour phases)
    def _cp_solvent(self):
        self.cp_solvent = Var(
            ['Liq', 'Vap'],
            initialize=4e3,
            bounds=(1e-8, 1e8),
            units=pyunits.J / pyunits.kg / pyunits.K,
            doc="Specific heat capacity of pure solvent")

        def rule_cp_solvent(b, p):
            if p == 'Liq':
                # specific heat, eq. 9 in Sharqawy et al. (2010)
                # Convert T90 to T68, eq. 4 in Sharqawy et al. (2010); primary reference from Rusby (1991)
                t = (b.temperature - 0.00025 * 273.15 * pyunits.K) / (1 - 0.00025)
                A = b.params.cp_phase_param_A1
                B = b.params.cp_phase_param_B1
                C = b.params.cp_phase_param_C1
                D = b.params.cp_phase_param_D1
                return b.cp_solvent['Liq'] == (A + B * t + C * t ** 2 + D * t ** 3) * 1000
            elif p == 'Vap':
                t = b.temperature / 1000
                return b.cp_solvent['Vap'] == b.params.cp_vap_param_A \
                      + b.params.cp_vap_param_B * t \
                      + b.params.cp_vap_param_C * t ** 2 \
                      + b.params.cp_vap_param_D * t ** 3 \
                      + b.params.cp_vap_param_E / t ** 2

        self.eq_cp_solvent = Constraint(['Liq', 'Vap'], rule=rule_cp_solvent)


    # 9. Heat capacity of solid-phase NaCl crystals
    def _cp_solute(self):
        self.cp_solute = Var(['Liq', 'Sol'],
            initialize=1e3,
            bounds=(-1e4, 1e6),
            units=pyunits.J / pyunits.kg / pyunits.K,
            doc="Specific heat capacity of solid NaCl crystals")

        def rule_cp_solute(b, p):  
            if p == 'Sol':
                return b.cp_solute[p] == b.params.cp_param_NaCl
            if p == 'Liq': # density, eq. 11-12 of Laliberte (2009)
                t = b.temperature - 273.15 * pyunits.K
                alpha = (b.params.cp_param_NaCl_liq_A2 * t) + (b.params.cp_param_NaCl_liq_A4 * (1 - b.mass_frac_phase_comp['Liq', 'H2O'])) \
                        + (b.params.cp_param_NaCl_liq_A3 * exp(0.01 * pyunits.K**-1 * t))
                cp_nacl_liq = b.params.cp_param_NaCl_liq_A1 * exp(alpha) + \
                        b.params.cp_param_NaCl_liq_A5 * ((1 - b.mass_frac_phase_comp['Liq', 'H2O']) ** b.params.cp_param_NaCl_liq_A6)
                return b.cp_solute[p] == cp_nacl_liq * (1000 * pyunits.J / pyunits.kJ)

        self.eq_cp_solute = Constraint(['Liq', 'Sol'], rule=rule_cp_solute)


    # 10. cp of liquid solution (Water + NaCl)
    def _cp_phase(self):
        self.cp_phase = Var(
            ['Liq'],
            initialize=4e3,
            bounds=(1e-8, 1e8),
            units=pyunits.J/pyunits.kg/pyunits.K,
            doc="Specific heat capacity of liquid solution")

        def rule_cp_phase(b):  # heat capacity, eq. 10 of Laliberte (2009) paper
            return b.cp_phase['Liq'] == b.mass_frac_phase_comp['Liq', 'NaCl'] * b.cp_solute['Liq'] \
                                        + b.mass_frac_phase_comp['Liq', 'H2O'] * b.cp_solvent['Liq']
        self.eq_cp_phase = Constraint(rule=rule_cp_phase)

    # 11. Volumetric flow rate for each phase
    def _flow_vol_phase(self):
        self.flow_vol_phase = Var(
            self.params.phase_list,
            initialize=1,
            bounds=(1e-8, None),
            units=pyunits.m ** 3 / pyunits.s,
            doc="Volumetric flow rate")

        def rule_flow_vol_phase(b, p):
            if p == 'Liq':
                return (b.flow_vol_phase[p] == sum(b.flow_mass_phase_comp[p, j] for j in self.params.component_list)
                        / b.dens_mass_phase[p])
            elif p == 'Sol':
                return (b.flow_vol_phase[p] == sum(b.flow_mass_phase_comp[p, j] for j in self.params.component_list)
                        / b.dens_mass_solute['Sol'])
            elif p == 'Vap':
                return (b.flow_vol_phase[p] == sum(b.flow_mass_phase_comp[p, j] for j in self.params.component_list)
                        / b.dens_mass_solvent['Vap'])
        self.eq_flow_vol_phase = Constraint(self.params.phase_list, rule=rule_flow_vol_phase)


    # 12. Total volumetric flow rate
    def _flow_vol(self):

        def rule_flow_vol(b):
            return sum(b.flow_vol_phase[p] for p in self.params.phase_list)

        self.flow_vol = Expression(rule=rule_flow_vol)



    # Boilerplate Methods

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
        return {"flow_mass_phase": self.flow_mass_phase,
                "temperature": self.temperature,
                "pressure": self.pressure}







# def nacl_solubility_function():
#     import numpy as np 
#     from idaes.surrogate.pysmo import polynomial_regression as pr 
#     from matplotlib import pyplot as plt

#     # Data
#     sol_data = np.array([ [0, 356.5], [10, 357.2], [20, 358.9], [30, 360.9], [40, 363.7], [60, 370.4], [80, 379.3], [100, 389.9] ])
#     sol_data[:, 0] = sol_data[:, 0] + 273.15

#     # Model training
#     model = pr.PolynomialRegression(sol_data, sol_data, maximum_polynomial_order = 2, number_of_crossvalidations=10, training_split=0.9, overwrite=True, fname='nacl_sol.pickle')
#     model.get_feature_vector()
#     model.training()

#     # Model validation - eye test
#     xrange = np.linspace(273.15, 373.15, 101)
#     xrange = xrange.reshape(xrange.shape[0], 1)
#     y_pred = model.predict_output(xrange)
#     plt.plot(sol_data[:, 0], sol_data[:, 1], 'o', label='training points')
#     plt.plot(xrange, y_pred, label='Surrogate function')
#     plt.show()

#     return model