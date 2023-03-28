#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
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
Initial property package for H2O-NaCl system
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
from idaes.core.base.components import Solute, Solvent
from idaes.core.base.phases import LiquidPhase
from idaes.core.util.constants import Constants
from idaes.core.util.initialization import (
    fix_state_vars,
    revert_state_vars,
    solve_indexed_blocks,
)
from idaes.core.util.misc import extract_data
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

from watertap.core.util.scaling import transform_property_constraints

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("NaClParameterBlock")
class NaClParameterData(PhysicalParameterBlock):
    CONFIG = PhysicalParameterBlock.CONFIG()

    def build(self):
        """
        Callable method for Block construction.
        """
        super(NaClParameterData, self).build()

        self._state_block_class = NaClStateBlock

        # components
        self.H2O = Solvent()
        self.NaCl = Solute()

        # phases
        self.Liq = LiquidPhase()
        #self.Vap = VaporPhase()


        # molecular weight
        mw_comp_data = {"H2O": 18.01528e-3, "NaCl": 58.44e-3}
        self.mw_comp = Param(
            self.component_list,
            mutable=False,
            initialize=extract_data(mw_comp_data),
            units=pyunits.kg / pyunits.mol,
            doc="Molecular weight kg/mol",
        )

        # solubility, 0-450 C
        # Sparrow 2003, Eq 5: Xsat = param_0 + param_1 * T + param_2 * T**2
        solubility_param_dict = {"0": 0.2628, "1": 62.75e-6, "2": 1.084e-6}
        self.solubility_param = Var(
            solubility_param_dict.keys(),
            domain=Reals,
            initialize=solubility_param_dict,
            units=pyunits.dimensionless,
            doc="Mass density parameters",
        )

        # mass density parameters, 0-300 C
        # Sparrow 2003, Eq 7: density = A+BT+CT2+DT3+ET4, where A = ao + a1X + a2X2 + a3X3 + a4X4
        dens_mass_A_param_dict = {"0": 1.0001e3, "1": 0.7666e3,"2": -.0149e3, "3": 0.2663e3,"4": 0.8845e3}
        dens_mass_B_param_dict = {"0": -0.0214, "1": -3.496, "2": 10.02, "3": -6.56, "4": -31.37}
        dens_mass_C_param_dict = {"0": -5.263e-3, "1": 39.87e-3, "2": -176.2e-3, "3": 363.5e-3, "4": -7.784e-3}
        dens_mass_D_param_dict = {"0": 15.42e-6, "1": -167e-6, "2": 980.7e-6, "3": -2573e-6, "4": 876.6e-6}
        dens_mass_E_param_dict = {"0": -0.0276e-6, "1": 0.2978e-6, "2": -2.017e-6, "3": 6.345e-6, "4": -3.914e-6}

        self.dens_mass_param_A = Var(
            dens_mass_A_param_dict.keys(),
            domain=Reals,
            initialize=dens_mass_A_param_dict,
            units=pyunits.kg / pyunits.m**3,
            doc="Mass density parameter A",
        )
        self.dens_mass_param_B = Var(
            dens_mass_B_param_dict.keys(),
            domain=Reals,
            initialize=dens_mass_B_param_dict,
            units=pyunits.kg / pyunits.m ** 3,
            doc="Mass density parameter B",
        )
        self.dens_mass_param_C = Var(
            dens_mass_C_param_dict.keys(),
            domain=Reals,
            initialize=dens_mass_C_param_dict,
            units=pyunits.kg / pyunits.m ** 3,
            doc="Mass density parameter C",
        )
        self.dens_mass_param_D = Var(
            dens_mass_D_param_dict.keys(),
            domain=Reals,
            initialize=dens_mass_D_param_dict,
            units=pyunits.kg / pyunits.m ** 3,
            doc="Mass density parameter D",
        )
        self.dens_mass_param_E = Var(
            dens_mass_E_param_dict.keys(),
            domain=Reals,
            initialize=dens_mass_E_param_dict,
            units=pyunits.kg / pyunits.m ** 3,
            doc="Mass density parameter E",
        )

        # specific enthalpy parameters, 0-300 C
        # Sparrow 2003, Eq 8: h = A+BT+CT2+DT3+ET4, where A = ao + a1X + a2X2 + a3X3 + a4X4
        enth_A_param_dict = {"0": 0.0005e3, "1": 0.0378e3, "2": -0.3682e3, "3": -0.6529e3, "4": 2.89e3}
        enth_B_param_dict = {"0": 4.145, "1": -4.973, "2": 4.482, "3": 18.31, "4": -46.41}
        enth_C_param_dict = {"0": 0.0007, "1": -0.0059, "2": 0.0854, "3": -0.495, "4": 0.8255}
        enth_D_param_dict = {"0": -0.0048e-3, "1": 0.0639e-3, "2": -0.714e-3, "3": 3.273e-3, "4": -4.85e-3}
        enth_E_param_dict = {"0": 0.0202e-6, "1": -0.2432e-6, "2": 2.054e-6, "3": -8.211e-6, "4": 11.43e-6}

        self.enth_param_A = Var(
            enth_A_param_dict.keys(),
            domain=Reals,
            initialize=enth_A_param_dict,
            units=pyunits.kJ / pyunits.kg,
            doc="Specific enthalpy parameter A",
        )
        self.enth_param_B = Var(
            enth_B_param_dict.keys(),
            domain=Reals,
            initialize=enth_B_param_dict,
            units=pyunits.kJ / pyunits.kg,
            doc="Specific enthalpy parameter B",
        )
        self.enth_param_C = Var(
            enth_C_param_dict.keys(),
            domain=Reals,
            initialize=enth_C_param_dict,
            units=pyunits.kJ / pyunits.kg,
            doc="Specific enthalpy parameter C",
        )
        self.enth_param_D = Var(
            enth_D_param_dict.keys(),
            domain=Reals,
            initialize=enth_D_param_dict,
            units=pyunits.kJ / pyunits.kg,
            doc="Specific enthalpy parameter D",
        )
        self.enth_param_E = Var(
            enth_E_param_dict.keys(),
            domain=Reals,
            initialize=enth_E_param_dict,
            units=pyunits.kJ / pyunits.kg,
            doc="Specific enthalpy parameter E",
        )

        # vapor pressure parameters, 0-150 C
        # Sparrow 2003, Eq 6: Pvap = A+BT+CT2+DT3+ET4, where A = ao + a1X + a2X2 + a3X3 + a4X4
        vap_pressure_A1_param_dict = {"0": 0.9083e3, "1": -0.569e3, "2": 0.1945e3, "3": -3.736e3, "4": 2.82e3}
        vap_pressure_B1_param_dict = {"0": -0.0669e-3, "1": 0.0582e-3, "2": -0.1668e-3, "3": 0.676e-3, "4": -2.091e-3}
        vap_pressure_C1_param_dict = {"0": 7.541e-6, "1": -5.143e-6, "2": 6.482e-6, "3": -52.62e-6, "4": 115.7e-6}
        vap_pressure_D1_param_dict = {"0":-0.0922e-6, "1": 0.0649e-6, "2": -0.1313e-6, "3": 0.8024e-6, "4": -1.986e-6}
        vap_pressure_E1_param_dict = {"0": 1.237e-9, "1": -0.753e-9, "2": 0.1448e-9, "3": -6.964e-9, "4": 14.61e-9}

        # vapor pressure parameters, (150)-300 C
        vap_pressure_A2_param_dict = {"0": -3.248, "1": 7.081, "2": -49.93, "3": 219.6, "4": -308.5}
        vap_pressure_B2_param_dict = {"0": 0.0610, "1": -0.1185, "2": 0.7916, "3": -3.474, "4": 4.882}
        vap_pressure_C2_param_dict = {"0": -0.4109e-3, "1": 0.6789e-3, "2": -4.155e-3, "3": 18.34e-3, "4": -25.89e-3}
        vap_pressure_D2_param_dict = {"0": 1.13e-6, "1": -1.432e-6, "2": 7.169e-6, "3": -33.17e-6, "4": 47.45e-6}
        vap_pressure_E2_param_dict = {"0": 0, "1": 0, "2": 0, "3": 0, "4": 0}

        self.vap_pressure_A1_param = Var(
            vap_pressure_A1_param_dict.keys(),
            domain=Reals,
            initialize=vap_pressure_A1_param_dict,
            units=pyunits.Pa / 10**6 , #MPa
            doc="Vapor pressure parameters A1",
        )
        self.vap_pressure_A2_param = Var(
            vap_pressure_A2_param_dict.keys(),
            domain=Reals,
            initialize=vap_pressure_A2_param_dict,
            units=pyunits.Pa / 10 ** 6,  # MPa
            doc="Vapor pressure parameters A2",
        )
        self.vap_pressure_B1_param = Var(
            vap_pressure_B1_param_dict.keys(),
            domain=Reals,
            initialize=vap_pressure_B1_param_dict,
            units=pyunits.Pa / 10 ** 6,  # MPa
            doc="Vapor pressure parameters B1",
        )
        self.vap_pressure_B2_param = Var(
            vap_pressure_B2_param_dict.keys(),
            domain=Reals,
            initialize=vap_pressure_B2_param_dict,
            units=pyunits.Pa / 10 ** 6,  # MPa
            doc="Vapor pressure parameters B2",
        )
        self.vap_pressure_C1_param = Var(
            vap_pressure_C1_param_dict.keys(),
            domain=Reals,
            initialize=vap_pressure_C1_param_dict,
            units=pyunits.Pa / 10 ** 6,  # MPa
            doc="Vapor pressure parameters C1",
        )
        self.vap_pressure_C2_param = Var(
            vap_pressure_C2_param_dict.keys(),
            domain=Reals,
            initialize=vap_pressure_C2_param_dict,
            units=pyunits.Pa / 10 ** 6,  # MPa
            doc="Vapor pressure parameters C2",
        )
        self.vap_pressure_D1_param = Var(
            vap_pressure_D1_param_dict.keys(),
            domain=Reals,
            initialize=vap_pressure_D1_param_dict,
            units=pyunits.Pa / 10 ** 6,  # MPa
            doc="Vapor pressure parameters D1",
        )
        self.vap_pressure_D2_param = Var(
            vap_pressure_D2_param_dict.keys(),
            domain=Reals,
            initialize=vap_pressure_D1_param_dict,
            units=pyunits.Pa / 10 ** 6,  # MPa
            doc="Vapor pressure parameters D2",
        )
        self.vap_pressure_E1_param = Var(
            vap_pressure_E1_param_dict.keys(),
            domain=Reals,
            initialize=vap_pressure_E1_param_dict,
            units=pyunits.Pa / 10 ** 6,  # MPa
            doc="Vapor pressure parameters E1",
        )
        self.vap_pressure_E2_param = Var(
            vap_pressure_E2_param_dict.keys(),
            domain=Reals,
            initialize=vap_pressure_E2_param_dict,
            units=pyunits.Pa / 10 ** 6,  # MPa
            doc="Vapor pressure parameters E2",
        )

        # could add entropy if needed from Sparrow 2003, Eq. 9

        # thermal conductivity, 0.1-140 MPa, 0-200 C
        # Qasuem et al. 2021, k = k,H2O + k,ion-solvent, correlation model no. 1
        th_cond_param_water_dict = {'0': -0.8624, '1': 8.801e-3, '1P': 1.57e-6, "2": -1.5113e-5, "3": 7e-9}
        th_cond_param_Na_dict = {'1': 0, '2': 0}
        th_cond_param_Cl_dict = {'1': -0.360439,'2': 0.006076}

        self.th_cond_param_water = Var(
           th_cond_param_water_dict.keys(),
            domain=Reals,
            initialize=th_cond_param_water_dict,
            units=pyunits.J / pyunits.s * pyunits.m * pyunits.K,
            doc="Thermal conductivity parameters for pure water",
        )
        self.th_cond_param_Na = Var(
            th_cond_param_Na_dict.keys(),
            domain=Reals,
            initialize=th_cond_param_Na_dict,
            units=pyunits.J / pyunits.s * pyunits.m * pyunits.K,
            doc="Thermal conductivity parameters for Na",
        )
        self.th_cond_param_Cl = Var(
            th_cond_param_Cl_dict.keys(),
            domain=Reals,
            initialize=th_cond_param_Cl_dict,
            units=pyunits.J / pyunits.s * pyunits.m * pyunits.K,
            doc="Thermal conductivity parameters for Cl",
        )

        # viscosity, 0-150 salinity (g/kg), 20-180 T (C)
        # Qasuem et al. 2021 (Sharqawy et al. correlation) #######################################
        # temperature (T) in °C and practical salinity (Sp) in g/kg
        visc_param_water_dict = {'0': -0.00379418, '1': 0.604129, '2': 139.18}
        visc_A_param_NaCl_dict = {'0': 1.474e-3, '1': 1.5e-5, '2': -3.927e-8}
        visc_B_param_NaCl_dict = {'0': 1.073e-5, '1': -8.5e-8, '2': 2.230e-10}

        self.visc_d_param_water = Var(
           visc_param_water_dict.keys(),
            domain=Reals,
            initialize=visc_param_water_dict,
            units=pyunits.dimensionless, ############################
            doc="Viscosity parameters for water",
        )
        self.visc_d_param_A_NaCl = Var(
            visc_A_param_NaCl_dict.keys(),
            domain=Reals,
            initialize=visc_A_param_NaCl_dict,
            units=pyunits.dimensionless,  ############################
            doc="Viscosity parameters for NaCl, A",
        )
        self.visc_d_param_B_NaCl = Var(
            visc_B_param_NaCl_dict.keys(),
            domain=Reals,
            initialize=visc_B_param_NaCl_dict,
            units=pyunits.dimensionless,  ############################
            doc="Viscosity parameters for NaCl, B",
        )

        # specific heat, 0 - 180 Salinity (g/kg),  0-180 T (C)
        # Generous et al. 2020 (Ghallab and Elnahas correlation),
        # temperature(T) in  °C, pressure(P) in MPa and salinity(S) in g / kg
        cp_param_T0_dict = {'0': 4.21, '1': -6.66, '2': 14.38, '3': -40.4, '4': 379.7, '5': -1692, '6': 2920}
        cp_param_T1_dict = {'0': 1.18e-3, '1': 5.62e-2, '2': -2.7e-1, '3': 2.5e-1, '4': -1.38, '5': 2.434, '6': 0}
        cp_param_T2_dict = {'0': 1.5e-5, '1': -5.77e-4, '2': 2.45e-3, '3': -3e-4, '4':1.4e-3, '5': 0, '6': 0}
        cp_param_T3_dict = {'0': -5.73e-8, '1': -1.79e-6, '2': -7.29e-6, '3': -4.61e-7, '4': 0, '5': 0, '6': 0}
        cp_param_T4_dict = {'0': 5.28e-10, '1': -4.25e-10, '2': 2.47e-10, '3': 0, '4': 0, '5': 0, '6': 0}
        cp_param_T5_dict = {'0': -2.27e-12, '1': 6.72e-13, '2': 0, '3': 0, '4': 0, '5': 0, '6': 0}
        cp_param_T6_dict = {'0': 3.75e-15, '1': 0, '2': 0, '3': 0, '4': 0, '5': 0, '6': 0}

        self.cp_param_T0 = Var(
            cp_param_T0_dict.keys(),
            domain=Reals,
            initialize=cp_param_T0_dict,
            units=pyunits.dimensionless, # May need to update units
            doc="Specific heat parameters for T0",
        )
        self.cp_param_T1 = Var(
            cp_param_T1_dict.keys(),
            domain=Reals,
            initialize=cp_param_T1_dict,
            units=pyunits.dimensionless,  # May need to update units
            doc="Specific heat parameters for T1",
        )
        self.cp_param_T2 = Var(
            cp_param_T2_dict.keys(),
            domain=Reals,
            initialize=cp_param_T2_dict,
            units=pyunits.dimensionless,  # May need to update units
            doc="Specific heat parameters for T2",
        )
        self.cp_param_T3 = Var(
            cp_param_T3_dict.keys(),
            domain=Reals,
            initialize=cp_param_T3_dict,
            units=pyunits.dimensionless,  # May need to update units
            doc="Specific heat parameters for T3",
        )
        self.cp_param_T4 = Var(
            cp_param_T4_dict.keys(),
            domain=Reals,
            initialize=cp_param_T4_dict,
            units=pyunits.dimensionless,  # May need to update units
            doc="Specific heat parameters for T4",
        )
        self.cp_param_T5 = Var(
            cp_param_T5_dict.keys(),
            domain=Reals,
            initialize=cp_param_T5_dict,
            units=pyunits.dimensionless,  # May need to update units
            doc="Specific heat parameters for T5",
        )
        self.cp_param_T6 = Var(
            cp_param_T6_dict.keys(),
            domain=Reals,
            initialize=cp_param_T6_dict,
            units=pyunits.dimensionless,  # May need to update units
            doc="Specific heat parameters for T6",
        )


        # specific heat of vaporization, 0-200 C, 1 atm, 0-120 g/kg
        # Qasuem et al. 2021 (Valderrama et al. correlation)
        h_vaporization_param_dict = {'a1': 2501e6, 'a2': -2.369e3, 'a3': 2.678e-1, 'a4': -8.103e-3, 'a5': -2.079e-5}
        self.h_vaporization_param = Var(
            h_vaporization_param_dict.keys(),
            domain=Reals,
            initialize=h_vaporization_param_dict,
            units=pyunits.dimensionless, # may need to update units
            doc="Latent heat of vaporization parameters",
        )

        # water vapor specific quanities

            # specific enthalpy
            # density

        # diffusivity (solute)

        # diffusivity (vapor)

        # osmostic coef.

        ##################################From Existing Model, not T dep###############################################

        # diffusivity parameters, eq 6 in Bartholomew
        diffus_param_dict = {
            "0": 1.51e-9,
            "1": -2.00e-9,
            "2": 3.01e-8,
            "3": -1.22e-7,
            "4": 1.53e-7,
        }
        self.diffus_param = Var(
            diffus_param_dict.keys(),
            domain=Reals,
            initialize=diffus_param_dict,
            units=pyunits.m**2 / pyunits.s,
            doc="Dynamic viscosity parameters",
        )

        # osmotic coefficient parameters, eq. 3b in Bartholomew
        osm_coeff_param_dict = {"0": 0.918, "1": 8.89e-2, "2": 4.92}
        self.osm_coeff_param = Var(
            osm_coeff_param_dict.keys(),
            domain=Reals,
            initialize=osm_coeff_param_dict,
            units=pyunits.dimensionless,
            doc="Osmotic coefficient parameters",
        )

        #################################################################################

        # traditional parameters are the only Vars currently on the block and should be fixed
        for v in self.component_objects(Var):
            v.fix()

        # ---default scaling---
        self.set_default_scaling("temperature", 1e-2)
        self.set_default_scaling("pressure", 1e-6)
        self.set_default_scaling("dens_mass_phase", 1e-3, index="Liq")
        self.set_default_scaling("visc_d_phase", 1e3, index="Liq")
        self.set_default_scaling("diffus_phase_comp", 1e9, index=("Liq", "NaCl"))
        self.set_default_scaling("osm_coeff", 1e0)
        self.set_default_scaling("enth_mass_phase", 1e-4, index="Liq")

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
                "diffus_phase_comp": {"method": "_diffus_phase_comp"},
                "visc_d_phase": {"method": "_visc_d_phase"},
                "pressure_osm_phase": {"method": "_pressure_osm_phase"},
                "enth_mass_phase": {"method": "_enth_mass_phase"},
            }
        )

        #obj.define_custom_properties(
        #    {
        #        "osm_coeff": {"method": "_osm_coeff"},
        #       "enth_flow": {"method": "_enth_flow"},
        #   }
        #)

        obj.add_default_units(
            {
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            }
        )


class _NaClStateBlock(StateBlock):
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

        # If input block, return flags, else release state
        if state_vars_fixed is False:
            if hold_state is True:
                return flags
            else:
                self.release_state(flags)

        if (not skip_solve) and (not check_optimal_termination(results)):
            raise InitializationError(
                f"{self.name} failed to initialize successfully. Please "
                f"check the output logs for more information."
            )

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


@declare_process_block_class("NaClStateBlock", block_class=_NaClStateBlock)
class NaClStateBlockData(StateBlockData):
    def build(self):
        """Callable method for Block construction."""
        super(NaClStateBlockData, self).build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        # Add state variables
        self.flow_mass_phase_comp = Var(
            self.params.phase_list,
            self.params.component_list,
            initialize={("Liq", "H2O"): 0.965, ("Liq", "NaCl"): 0.035},
            bounds=(0.0, None),
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
            initialize={("Liq", "H2O"): 0.965, ("Liq", "NaCl"): 0.035},
            # upper bound set to None because of stability benefits
            bounds=(0.0, None),
            units=pyunits.dimensionless,
            doc="Mass fraction",
        )

        def rule_mass_frac_phase_comp(b, p, j):
            return b.mass_frac_phase_comp[p, j] == b.flow_mass_phase_comp[p, j] / sum(
                b.flow_mass_phase_comp[p, j] for j in self.params.component_list
            )

        self.eq_mass_frac_phase_comp = Constraint(
            self.params.phase_list,
            self.params.component_list,
            rule=rule_mass_frac_phase_comp,
        )

    def _dens_mass_phase(self):
        self.dens_mass_phase = Var(
            self.params.phase_list,
            initialize=1e3,
            bounds=(5e2, 2e3),
            units=pyunits.kg * pyunits.m**-3,
            doc="Mass density",
        )

        def rule_dens_mass_phase(b, p):  # density, eq. 7 in Sparrow 2003
            t = (b.temperature - 273.15 * pyunits.K) / pyunits.K
            param_vec = [b.params.dens_mass_param_A, b.params.dens_mass_param_B, b.params.dens_mass_param_C,
                       b.params.dens_mass_param_D, b.params.dens_mass_param_E]
            iter_param = {'A': 0, 'B': 0, 'C': 0, 'D': 0, 'E': 0}
            k = 0
            for key in iter_param:
                iter_param[key] = (param_vec[k]['0']
                + param_vec[k]['1'] * b.mass_frac_phase_comp[p, "NaCl"]
                + param_vec[k]['2'] * b.mass_frac_phase_comp[p, "NaCl"] ** 2
                + param_vec[k]['3'] * b.mass_frac_phase_comp[p, "NaCl"] ** 3
                + param_vec[k]['4'] * b.mass_frac_phase_comp[p, "NaCl"] ** 4)
                k += 1
            return (
                b.dens_mass_phase[p] ==
                iter_param['A'] + iter_param['B']*t + iter_param['C']*t**2 + iter_param['D']*t**3 + iter_param['E']*t**4
            )

        self.eq_dens_mass_phase = Constraint(
            self.params.phase_list, rule=rule_dens_mass_phase
        )

    def _flow_vol_phase(self):
        self.flow_vol_phase = Var(
            self.params.phase_list,
            initialize=1,
            bounds=(0.0, None),
            units=pyunits.m**3 / pyunits.s,
            doc="Volumetric flow rate",
        )

        def rule_flow_vol_phase(b, p):
            return (
                b.flow_vol_phase[p]
                == sum(b.flow_mass_phase_comp[p, j] for j in self.params.component_list)
                / b.dens_mass_phase[p]
            )

        self.eq_flow_vol_phase = Constraint(
            self.params.phase_list, rule=rule_flow_vol_phase
        )

    def _flow_vol(self):
        def rule_flow_vol(b):
            return sum(b.flow_vol_phase[p] for p in self.params.phase_list)

        self.flow_vol = Expression(rule=rule_flow_vol)

    def _conc_mass_phase_comp(self):
        self.conc_mass_phase_comp = Var(
            self.params.phase_list,
            self.params.component_list,
            initialize=10,
            bounds=(1e-3, 2e3),
            units=pyunits.kg * pyunits.m**-3,
            doc="Mass concentration",
        )

        def rule_conc_mass_phase_comp(b, p, j):
            return (
                b.conc_mass_phase_comp[p, j]
                == b.dens_mass_phase[p] * b.mass_frac_phase_comp[p, j]
            )

        self.eq_conc_mass_phase_comp = Constraint(
            self.params.phase_list,
            self.params.component_list,
            rule=rule_conc_mass_phase_comp,
        )

    def _flow_mol_phase_comp(self):
        self.flow_mol_phase_comp = Var(
            self.params.phase_list,
            self.params.component_list,
            initialize=100,
            bounds=(0.0, None),
            units=pyunits.mol / pyunits.s,
            doc="Molar flowrate",
        )

        def rule_flow_mol_phase_comp(b, p, j):
            return (
                b.flow_mol_phase_comp[p, j]
                == b.flow_mass_phase_comp[p, j] / b.params.mw_comp[j]
            )

        self.eq_flow_mol_phase_comp = Constraint(
            self.params.phase_list,
            self.params.component_list,
            rule=rule_flow_mol_phase_comp,
        )

    def _mole_frac_phase_comp(self):
        self.mole_frac_phase_comp = Var(
            self.params.phase_list,
            self.params.component_list,
            initialize=0.1,
            bounds=(0.0, None),
            units=pyunits.dimensionless,
            doc="Mole fraction",
        )

        def rule_mole_frac_phase_comp(b, p, j):
            return b.mole_frac_phase_comp[p, j] == b.flow_mol_phase_comp[p, j] / sum(
                b.flow_mol_phase_comp[p, j] for j in b.params.component_list
            )

        self.eq_mole_frac_phase_comp = Constraint(
            self.params.phase_list,
            self.params.component_list,
            rule=rule_mole_frac_phase_comp,
        )

    def _molality_phase_comp(self):
        self.molality_phase_comp = Var(
            self.params.phase_list,
            ["NaCl"],
            initialize=1,
            bounds=(1e-4, 10),
            units=pyunits.mole / pyunits.kg,
            doc="Molality",
        )

        def rule_molality_phase_comp(b, p, j):
            return (
                self.molality_phase_comp[p, j]
                == b.mass_frac_phase_comp[p, j]
                / (1 - b.mass_frac_phase_comp[p, j])
                / b.params.mw_comp[j]
            )

        self.eq_molality_phase_comp = Constraint(
            self.params.phase_list, ["NaCl"], rule=rule_molality_phase_comp
        )
    #######################################################
    def _visc_d_phase(self):
        self.visc_d_phase = Var(
            self.params.phase_list,
            initialize=1e-3,
            bounds=(1e-4, 1e-2),
            units=pyunits.Pa * pyunits.s,
            doc="Viscosity",
        )

        def rule_visc_d_phase(b, p):  # dynamic viscosity, Sharqawy et al. correlation
            #t = b.temperature - 273.15 * pyunits.K
            #s = b.flow_mass_phase_comp[p, 'NaCl']/ b.flow_mass_phase_comp[p, 'H2O'] * 1000 # kg/kg -> g/kg
            #visc_purewater = b.params.visc_d_param['pure water']['0'] + b.params.visc_d_param['pure water']['1']/(b.params.visc_d_param['pure water']['0'] + t)
            #A_param = b.params.visc_d_param['salt water']['A']['0'] + b.params.visc_d_param['salt water']['A']['1']*t + b.params.visc_d_param['salt water']['A']['2']*t**2
            #B_param = b.params.visc_d_param['salt water']['B']['0'] + b.params.visc_d_param['salt water']['B']['1']*t + b.params.visc_d_param['salt water']['B']['2']*t**2

            return b.visc_d_phase[p] == (
                    1 #visc_purewater * (1 + A_param * s + B_param * s ** 2) / 1000 #mPas to Pas
            )

        self.eq_visc_d_phase = Constraint(
            self.params.phase_list, rule=rule_visc_d_phase
        )

    def _diffus_phase_comp(self):
        self.diffus_phase_comp = Var(
            self.params.phase_list,
            ["NaCl"],
            initialize=1e-9,
            bounds=(1e-10, 1e-8),
            units=pyunits.m**2 * pyunits.s**-1,
            doc="Diffusivity",
        )

        def rule_diffus_phase_comp(b, p, j):  # diffusivity, eq 6 in Bartholomew
            return b.diffus_phase_comp[p, j] == (
                b.params.diffus_param["4"] * b.mass_frac_phase_comp[p, "NaCl"] ** 4
                + b.params.diffus_param["3"] * b.mass_frac_phase_comp[p, "NaCl"] ** 3
                + b.params.diffus_param["2"] * b.mass_frac_phase_comp[p, "NaCl"] ** 2
                + b.params.diffus_param["1"] * b.mass_frac_phase_comp[p, "NaCl"]
                + b.params.diffus_param["0"]
            )

        self.eq_diffus_phase_comp = Constraint(
            self.params.phase_list, ["NaCl"], rule=rule_diffus_phase_comp
        )

    def _osm_coeff(self):
        self.osm_coeff = Var(
            initialize=1,
            bounds=(0.5, 2),
            units=pyunits.dimensionless,
            doc="Osmotic coefficient",
        )

        def rule_osm_coeff(b):
            return b.osm_coeff == (
                b.params.osm_coeff_param["2"]
                * b.mass_frac_phase_comp["Liq", "NaCl"] ** 2
                + b.params.osm_coeff_param["1"] * b.mass_frac_phase_comp["Liq", "NaCl"]
                + b.params.osm_coeff_param["0"]
            )

        self.eq_osm_coeff = Constraint(rule=rule_osm_coeff)

    def _pressure_osm_phase(self):
        self.pressure_osm_phase = Var(
            self.params.phase_list,
            initialize=1e6,
            bounds=(5e2, 5e7),
            units=pyunits.Pa,
            doc="Osmotic pressure",
        )

        def rule_pressure_osm_phase(b, p):
            i = 2  # number of ionic species
            rhow = (
                1000 * pyunits.kg / pyunits.m**3
            )  # TODO: could make this variable based on temperature
            return (
                b.pressure_osm_phase[p]
                == i
                * b.osm_coeff
                * b.molality_phase_comp[p, "NaCl"]
                * rhow
                * Constants.gas_constant
                * b.temperature
            )

        self.eq_pressure_osm_phase = Constraint(
            self.params.phase_list, rule=rule_pressure_osm_phase
        )

    def _enth_mass_phase(self):
        self.enth_mass_phase = Var(
            self.params.phase_list,
            initialize=5e4,
            bounds=(1e4, 1e6),
            units=pyunits.kJ * pyunits.kg**-1,
            doc="Specific enthalpy",
        )

        def rule_enth_mass_phase(b, p):  # specific enthalpy, H' = Cp(T-Tref) + (P-Pref)/rho
            t = (b.temperature - 273.15 * pyunits.K) / pyunits.K #pyunits.convert(b.temperature, to_units=pyunits.degC )/pyunits.degC
            param_vec = [b.params.enth_param_A, b.params.enth_param_B, b.params.enth_param_C,
                         b.params.enth_param_D,b.params.enth_param_E]
            iter_param = {'A': 0, 'B': 0, 'C': 0, 'D': 0, 'E': 0}
            k = 0
            for key in iter_param:
                iter_param[key] = (param_vec[k]['0']
                + param_vec[k]['1'] * b.mass_frac_phase_comp[p, "NaCl"]
                + param_vec[k]['2'] * b.mass_frac_phase_comp[p, "NaCl"] ** 2
                + param_vec[k]['3'] * b.mass_frac_phase_comp[p, "NaCl"] ** 3
                + param_vec[k]['4'] * b.mass_frac_phase_comp[p, "NaCl"] ** 4)
                k += 1
            return (
                    b.enth_mass_phase[p] ==
                    iter_param['A'] + iter_param['B'] * t + iter_param['C'] * t ** 2 + iter_param['D'] * t ** 3 +
                    iter_param['E'] * t ** 4
            )

        self.eq_enth_mass_phase = Constraint(
            self.params.phase_list, rule=rule_enth_mass_phase
        )

    def _enth_flow(self):
        # enthalpy flow expression for get_enthalpy_flow_terms method

        def rule_enth_flow(b):  # enthalpy flow [J/s]
            return (
                sum(b.flow_mass_phase_comp["Liq", j] for j in b.params.component_list)
                * b.enth_mass_phase["Liq"]
            )

        self.enth_flow = Expression(rule=rule_enth_flow)

    # TODO: add vapor pressure, specific heat, thermal conductivity,
    #   and heat of vaporization

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
        # temperature, dens_mass, visc_d, diffus, osm_coeff, and enth_mass

        # these variables should have user input
        if iscale.get_scaling_factor(self.flow_mass_phase_comp["Liq", "H2O"]) is None:
            sf = iscale.get_scaling_factor(
                self.flow_mass_phase_comp["Liq", "H2O"], default=1e0, warning=True
            )
            iscale.set_scaling_factor(self.flow_mass_phase_comp["Liq", "H2O"], sf)

        if iscale.get_scaling_factor(self.flow_mass_phase_comp["Liq", "NaCl"]) is None:
            sf = iscale.get_scaling_factor(
                self.flow_mass_phase_comp["Liq", "NaCl"], default=1e2, warning=True
            )
            iscale.set_scaling_factor(self.flow_mass_phase_comp["Liq", "NaCl"], sf)

        # scaling factors for parameters
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
                    if j == "NaCl":
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
                    elif j == "NaCl":
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
                    if j == "NaCl":
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
        # transforming constraints
        if self.is_property_constructed("pressure_osm_phase"):
            sf = iscale.get_scaling_factor(
                self.pressure_osm_phase["Liq"], default=1, warning=True
            )
            iscale.constraint_scaling_transform(self.eq_pressure_osm_phase["Liq"], sf)
        if self.is_property_constructed("osm_coeff"):
            sf = iscale.get_scaling_factor(self.osm_coeff, default=1, warning=True)
            iscale.constraint_scaling_transform(self.eq_osm_coeff, sf)

        # transforming constraints
        transform_property_constraints(self)
