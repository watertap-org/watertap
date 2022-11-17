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

from enum import Enum, auto

# Import Pyomo libraries
from pyomo.environ import (
    Constraint,
    Expression,
    Reals,
    NonNegativeReals,
    Var,
    Param,
    exp,
    log,
    Suffix,
    value,
    check_optimal_termination,
)
from pyomo.environ import units as pyunits
from pyomo.common.config import ConfigValue, In

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
from idaes.core.base.phases import (
    LiquidPhase,
    VaporPhase,
    SolidPhase,
    Phase,
    PhaseType as PT,
)
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

__author__ = "Oluwamayowa Amusat"


class HeatOfCrystallizationModel(Enum):
    constant = auto()  # Use constant heat of crystallization
    zero = auto()  # Assume heat of crystallization is zero
    temp_dependent = auto()  # Use temperature-dependent heat of crystallization


@declare_process_block_class("NaClParameterBlock")
class NaClParameterData(PhysicalParameterBlock):
    CONFIG = PhysicalParameterBlock.CONFIG()

    CONFIG.declare(
        "heat_of_crystallization_model",
        ConfigValue(
            default=HeatOfCrystallizationModel.constant,
            domain=In(HeatOfCrystallizationModel),
            description="Heat of crystallization construction flag",
            doc="""
           Options to account for heat of crystallization value for NaCl.

           **default** - ``HeatOfCrystallizationModel.constant``

       .. csv-table::
           :header: "Configuration Options", "Description"

           "``HeatOfCrystallizationModel.constant``", "Fixed heat of crystallization for NaCl based on literature"
           "``HeatOfCrystallizationModel.zero``", "Zero heat of crystallization assumption"
           "``HeatOfCrystallizationModel.temp_dependent``", "Temperature-dependent heat of crystallization for NaCl"
       """,
        ),
    )

    def build(self):
        super().build()
        self._state_block_class = NaClStateBlock

        # Component
        self.H2O = Solvent(valid_phase_types=[PT.liquidPhase, PT.vaporPhase])
        self.NaCl = Solute(valid_phase_types=[PT.liquidPhase, PT.solidPhase])

        # Phases
        self.Liq = LiquidPhase(component_list=["H2O", "NaCl"])
        self.Vap = VaporPhase(component_list=["H2O"])
        self.Sol = SolidPhase(component_list=["NaCl"])

        """
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

        - Sparrow, B. S. Empirical equations for the thermodynamic properties of aqueous sodium chloride
        Desalination, Elsevier BV, 2003, 159, 161-170

        - Chase, M.W., Jr., NIST-JANAF Themochemical Tables, Fourth Edition, J. Phys. Chem. Ref. Data, Monograph 9,
        1998, 1-1951.

        - El-Dessouky, H. T. & Ettouney, H. M. (Eds.). Appendix A - Thermodynamic Properties.
        Fundamentals of Salt Water Desalination, Elsevier Science B.V., 2002, 525-563

        - Tavare, N. S. Industrial Crystallization, Springer US, 2013.
        """

        # Unit definitions
        dens_units = pyunits.kg / pyunits.m**3
        t_inv_units = pyunits.K**-1
        enth_mass_units = pyunits.J / pyunits.kg
        enth_mass_units_2 = pyunits.kJ / pyunits.kg
        cp_units = pyunits.J / (pyunits.kg * pyunits.K)
        cp_units_2 = pyunits.kJ / (pyunits.kg * pyunits.K)
        cp_units_3 = pyunits.J / (pyunits.mol * pyunits.K)

        # molecular weights of solute and solvent
        mw_comp_data = {"H2O": 18.01528e-3, "NaCl": 58.44e-3}
        self.mw_comp = Param(
            self.component_list,
            mutable=False,
            initialize=extract_data(mw_comp_data),
            units=pyunits.kg / pyunits.mol,
            doc="Molecular weight kg/mol",
        )

        # Solubility parameters from (1) surrogate model (2) Sparrow paper
        # # 1. Surrogate
        # self.sol_param_A1 = Param(initialize= 526.706475, units = pyunits.g / pyunits.L, doc=' Solubility parameter A1 [g/L] for NaCl surrogate')
        # self.sol_param_A2 = Param(initialize= -1.326952, units = (pyunits.g / pyunits.L) * pyunits.K ** -1, doc=' Solubility parameter A2 [g/L] for NaCl surrogate')
        # self.sol_param_A3 = Param(initialize= 0.002574, units = (pyunits.g / pyunits.L) * pyunits.K ** -2, doc=' Solubility parameter A3 [g/L] for NaCl surrogate')
        # 2. Sparrow
        self.sol_param_A1 = Param(
            initialize=0.2628,
            units=pyunits.dimensionless,
            doc=" Solubility parameter A1 for NaCl",
        )
        self.sol_param_A2 = Param(
            initialize=62.75e-6,
            units=pyunits.K**-1,
            doc=" Solubility parameter A2 for NaCl",
        )
        self.sol_param_A3 = Param(
            initialize=1.084e-6,
            units=pyunits.K**-2,
            doc=" Solubility parameter A3 for NaCl",
        )

        # Mass density value for NaCl crystals in solid phase: fixed for now at Tavare value - may not be accurate?
        self.dens_mass_param_NaCl = Param(
            initialize=2115,
            units=pyunits.kg / pyunits.m**3,
            doc="NaCl crystal density",
        )

        # Heat of crystallization parameter - fixed value based on heat of fusion from Perry (Table 2-147)
        self.dh_crystallization_param = Param(
            initialize=-520, units=enth_mass_units_2, doc="NaCl heat of crystallization"
        )

        # Mass density parameters for pure NaCl liquid based on Eq. 9 in Laliberte and Cooper (2004).
        self.dens_mass_param_NaCl_liq_C0 = Var(
            within=Reals,
            initialize=-0.00433,
            units=dens_units,
            doc="Mass density parameter C0 for liquid NaCl",
        )
        self.dens_mass_param_NaCl_liq_C1 = Var(
            within=Reals,
            initialize=0.06471,
            units=dens_units,
            doc="Mass density parameter C1 for liquid NaCl",
        )
        self.dens_mass_param_NaCl_liq_C2 = Var(
            within=Reals,
            initialize=1.01660,
            units=pyunits.dimensionless,
            doc="Mass density parameter C2 for liquid NaCl",
        )
        self.dens_mass_param_NaCl_liq_C3 = Var(
            within=Reals,
            initialize=0.014624,
            units=t_inv_units,
            doc="Mass density parameter C3 for liquid NaCl",
        )
        self.dens_mass_param_NaCl_liq_C4 = Var(
            within=Reals,
            initialize=3315.6,
            units=pyunits.K,
            doc="Mass density parameter C4 for liquid NaCl",
        )

        # Mass density parameters for solvent in liquid phase, eq. 8 in Sharqawy et al. (2010)
        self.dens_mass_param_A1 = Var(
            within=Reals,
            initialize=9.999e2,
            units=dens_units,
            doc="Mass density parameter A1",
        )
        self.dens_mass_param_A2 = Var(
            within=Reals,
            initialize=2.034e-2,
            units=dens_units * t_inv_units,
            doc="Mass density parameter A2",
        )
        self.dens_mass_param_A3 = Var(
            within=Reals,
            initialize=-6.162e-3,
            units=dens_units * t_inv_units**2,
            doc="Mass density parameter A3",
        )
        self.dens_mass_param_A4 = Var(
            within=Reals,
            initialize=2.261e-5,
            units=dens_units * t_inv_units**3,
            doc="Mass density parameter A4",
        )
        self.dens_mass_param_A5 = Var(
            within=Reals,
            initialize=-4.657e-8,
            units=dens_units * t_inv_units**4,
            doc="Mass density parameter A5",
        )

        # Latent heat of evaporation of pure water: Parameters from Sharqawy et al. (2010), eq. 54
        self.dh_vap_w_param_0 = Var(
            within=Reals,
            initialize=2.501e6,
            units=enth_mass_units,
            doc="Latent heat of pure water parameter 0",
        )
        self.dh_vap_w_param_1 = Var(
            within=Reals,
            initialize=-2.369e3,
            units=enth_mass_units * t_inv_units**1,
            doc="Latent heat of pure water parameter 1",
        )
        self.dh_vap_w_param_2 = Var(
            within=Reals,
            initialize=2.678e-1,
            units=enth_mass_units * t_inv_units**2,
            doc="Latent heat of pure water parameter 2",
        )
        self.dh_vap_w_param_3 = Var(
            within=Reals,
            initialize=-8.103e-3,
            units=enth_mass_units * t_inv_units**3,
            doc="Latent heat of pure water parameter 3",
        )
        self.dh_vap_w_param_4 = Var(
            within=Reals,
            initialize=-2.079e-5,
            units=enth_mass_units * t_inv_units**4,
            doc="Latent heat of pure water parameter 4",
        )

        # Specific heat parameters for Cp vapor from NIST Webbook - Chase, M.W., Jr., NIST-JANAF Themochemical Tables
        self.cp_vap_param_A = Var(
            within=Reals,
            initialize=30.09200 / 18.01528e-3,
            units=cp_units,
            doc="Specific heat of water vapor parameter A",
        )
        self.cp_vap_param_B = Var(
            within=Reals,
            initialize=6.832514 / 18.01528e-3,
            units=cp_units * t_inv_units,
            doc="Specific heat of water vapor parameter B",
        )
        self.cp_vap_param_C = Var(
            within=Reals,
            initialize=6.793435 / 18.01528e-3,
            units=cp_units * t_inv_units**2,
            doc="Specific heat of water vapor parameter C",
        )
        self.cp_vap_param_D = Var(
            within=Reals,
            initialize=-2.534480 / 18.01528e-3,
            units=cp_units * t_inv_units**3,
            doc="Specific heat of water vapor parameter D",
        )
        self.cp_vap_param_E = Var(
            within=Reals,
            initialize=0.082139 / 18.01528e-3,
            units=cp_units * t_inv_units**-2,
            doc="Specific heat of water vapor parameter E",
        )

        # Specific heat parameters for pure water from eq (9) in Sharqawy et al. (2010)
        self.cp_phase_param_A1 = Var(
            within=Reals,
            initialize=5.328,
            units=cp_units,
            doc="Specific heat of seawater parameter A1",
        )
        self.cp_phase_param_B1 = Var(
            within=Reals,
            initialize=-6.913e-3,
            units=cp_units * t_inv_units,
            doc="Specific heat of seawater parameter B1",
        )
        self.cp_phase_param_C1 = Var(
            within=Reals,
            initialize=9.6e-6,
            units=cp_units * t_inv_units**2,
            doc="Specific heat of seawater parameter C1",
        )
        self.cp_phase_param_D1 = Var(
            within=Reals,
            initialize=2.5e-9,
            units=cp_units * t_inv_units**3,
            doc="Specific heat of seawater parameter D1",
        )

        # Specific heat parameters for liquid NaCl from eqs. (11) & (12) in Laliberte (2009).
        self.cp_param_NaCl_liq_A1 = Var(
            within=Reals,
            initialize=-0.06936,
            units=cp_units_2,
            doc="Specific heat parameter A1 for liquid NaCl",
        )
        self.cp_param_NaCl_liq_A2 = Var(
            within=Reals,
            initialize=-0.07821,
            units=t_inv_units,
            doc="Specific heat parameter A2 for liquid NaCl",
        )
        self.cp_param_NaCl_liq_A3 = Var(
            within=Reals,
            initialize=3.8480,
            units=pyunits.dimensionless,
            doc="Specific heat parameter A3 for liquid NaCl",
        )
        self.cp_param_NaCl_liq_A4 = Var(
            within=Reals,
            initialize=-11.2762,
            units=pyunits.dimensionless,
            doc="Specific heat parameter A4 for liquid NaCl",
        )
        self.cp_param_NaCl_liq_A5 = Var(
            within=Reals,
            initialize=8.7319,
            units=cp_units_2,
            doc="Specific heat parameter A5 for liquid NaCl",
        )
        self.cp_param_NaCl_liq_A6 = Var(
            within=Reals,
            initialize=1.8125,
            units=pyunits.dimensionless,
            doc="Specific heat parameter A6 for liquid NaCl",
        )

        # Specific heat parameters for solid NaCl : Shomate equation from NIST webbook (https://webbook.nist.gov/cgi/cbook.cgi?ID=C7647145&Mask=6F).
        self.cp_param_NaCl_solid_A = Var(
            within=Reals,
            initialize=50.72389,
            units=cp_units_3,
            doc="Specific heat parameter A for solid NaCl",
        )
        self.cp_param_NaCl_solid_B = Var(
            within=Reals,
            initialize=6.672267,
            units=cp_units_3 / pyunits.K,
            doc="Specific heat parameter B for solid NaCl",
        )
        self.cp_param_NaCl_solid_C = Var(
            within=Reals,
            initialize=-2.517167,
            units=cp_units_3 / pyunits.K**2,
            doc="Specific heat parameter C for solid NaCl",
        )
        self.cp_param_NaCl_solid_D = Var(
            within=Reals,
            initialize=10.15934,
            units=cp_units_3 / pyunits.K**3,
            doc="Specific heat parameter D for solid NaCl",
        )
        self.cp_param_NaCl_solid_E = Var(
            within=Reals,
            initialize=-0.200675,
            units=cp_units_3 * pyunits.K**2,
            doc="Specific heat parameter E for solid NaCl",
        )
        self.cp_param_NaCl_solid_F = Var(
            within=Reals,
            initialize=-427.2115,
            units=cp_units_3 * pyunits.K,
            doc="Specific heat parameter F for solid NaCl",
        )
        # self.cp_param_NaCl_solid_G = Var(within=Reals, initialize=130.3973, units=cp_units_2, doc='Specific heat parameter G for solid NaCl')
        self.cp_param_NaCl_solid_H = Var(
            within=Reals,
            initialize=-411.1203,
            units=cp_units_3 * pyunits.K,
            doc="Specific heat parameter H for solid NaCl",
        )

        # Vapour pressure parameters for NaCl solution from Sparrow (2003): 0 < T < 150 degC
        self.pressure_sat_param_A1 = Var(
            within=Reals,
            initialize=0.9083e-3,
            units=pyunits.MPa,
            doc="Vapour pressure parameter A1",
        )
        self.pressure_sat_param_A2 = Var(
            within=Reals,
            initialize=-0.569e-3,
            units=pyunits.MPa,
            doc="Vapour pressure parameter A2",
        )
        self.pressure_sat_param_A3 = Var(
            within=Reals,
            initialize=0.1945e-3,
            units=pyunits.MPa,
            doc="Vapour pressure parameter A3",
        )
        self.pressure_sat_param_A4 = Var(
            within=Reals,
            initialize=-3.736e-3,
            units=pyunits.MPa,
            doc="Vapour pressure parameter A4",
        )
        self.pressure_sat_param_A5 = Var(
            within=Reals,
            initialize=2.82e-3,
            units=pyunits.MPa,
            doc="Vapour pressure parameter A5",
        )
        self.pressure_sat_param_B1 = Var(
            within=Reals,
            initialize=-0.0669e-3,
            units=pyunits.MPa / pyunits.K,
            doc="Vapour pressure parameter B1",
        )
        self.pressure_sat_param_B2 = Var(
            within=Reals,
            initialize=0.0582e-3,
            units=pyunits.MPa / pyunits.K,
            doc="Vapour pressure parameter B2",
        )
        self.pressure_sat_param_B3 = Var(
            within=Reals,
            initialize=-0.1668e-3,
            units=pyunits.MPa / pyunits.K,
            doc="Vapour pressure parameter B3",
        )
        self.pressure_sat_param_B4 = Var(
            within=Reals,
            initialize=0.6761e-3,
            units=pyunits.MPa / pyunits.K,
            doc="Vapour pressure parameter B4",
        )
        self.pressure_sat_param_B5 = Var(
            within=Reals,
            initialize=-2.091e-3,
            units=pyunits.MPa / pyunits.K,
            doc="Vapour pressure parameter B5",
        )
        self.pressure_sat_param_C1 = Var(
            within=Reals,
            initialize=7.541e-6,
            units=pyunits.MPa / pyunits.K**2,
            doc="Vapour pressure parameter C1",
        )
        self.pressure_sat_param_C2 = Var(
            within=Reals,
            initialize=-5.143e-6,
            units=pyunits.MPa / pyunits.K**2,
            doc="Vapour pressure parameter C2",
        )
        self.pressure_sat_param_C3 = Var(
            within=Reals,
            initialize=6.482e-6,
            units=pyunits.MPa / pyunits.K**2,
            doc="Vapour pressure parameter C3",
        )
        self.pressure_sat_param_C4 = Var(
            within=Reals,
            initialize=-52.62e-6,
            units=pyunits.MPa / pyunits.K**2,
            doc="Vapour pressure parameter C4",
        )
        self.pressure_sat_param_C5 = Var(
            within=Reals,
            initialize=115.7e-6,
            units=pyunits.MPa / pyunits.K**2,
            doc="Vapour pressure parameter C5",
        )
        self.pressure_sat_param_D1 = Var(
            within=Reals,
            initialize=-0.0922e-6,
            units=pyunits.MPa / pyunits.K**3,
            doc="Vapour pressure parameter D1",
        )
        self.pressure_sat_param_D2 = Var(
            within=Reals,
            initialize=0.0649e-6,
            units=pyunits.MPa / pyunits.K**3,
            doc="Vapour pressure parameter D2",
        )
        self.pressure_sat_param_D3 = Var(
            within=Reals,
            initialize=-0.1313e-6,
            units=pyunits.MPa / pyunits.K**3,
            doc="Vapour pressure parameter D3",
        )
        self.pressure_sat_param_D4 = Var(
            within=Reals,
            initialize=0.8024e-6,
            units=pyunits.MPa / pyunits.K**3,
            doc="Vapour pressure parameter D4",
        )
        self.pressure_sat_param_D5 = Var(
            within=Reals,
            initialize=-1.986e-6,
            units=pyunits.MPa / pyunits.K**3,
            doc="Vapour pressure parameter D5",
        )
        self.pressure_sat_param_E1 = Var(
            within=Reals,
            initialize=1.237e-9,
            units=pyunits.MPa / pyunits.K**4,
            doc="Vapour pressure parameter E1",
        )
        self.pressure_sat_param_E2 = Var(
            within=Reals,
            initialize=-0.753e-9,
            units=pyunits.MPa / pyunits.K**4,
            doc="Vapour pressure parameter E2",
        )
        self.pressure_sat_param_E3 = Var(
            within=Reals,
            initialize=0.1448e-9,
            units=pyunits.MPa / pyunits.K**4,
            doc="Vapour pressure parameter E3",
        )
        self.pressure_sat_param_E4 = Var(
            within=Reals,
            initialize=-6.964e-9,
            units=pyunits.MPa / pyunits.K**4,
            doc="Vapour pressure parameter E4",
        )
        self.pressure_sat_param_E5 = Var(
            within=Reals,
            initialize=14.61e-9,
            units=pyunits.MPa / pyunits.K**4,
            doc="Vapour pressure parameter E5",
        )

        # Parameters for saturation temperature of water vapour from eq. A.12 in El-Dessouky and Ettouney
        self.temp_sat_solvent_A1 = Var(
            within=Reals,
            initialize=42.6776,
            units=pyunits.K,
            doc="Water boiling point parameter A1",
        )
        self.temp_sat_solvent_A2 = Var(
            within=Reals,
            initialize=-3892.7,
            units=pyunits.K,
            doc="Water boiling point parameter A2",
        )
        self.temp_sat_solvent_A3 = Var(
            within=Reals,
            initialize=1000,
            units=pyunits.kPa,
            doc="Water boiling point parameter A3",
        )
        self.temp_sat_solvent_A4 = Var(
            within=Reals,
            initialize=-9.48654,
            units=pyunits.dimensionless,
            doc="Water boiling point parameter A4",
        )

        # Parameters for specific enthalpy of pure water in liquid phase from eq. 55 in Sharqawy et al. (2010)
        self.enth_mass_solvent_param_A1 = Var(
            within=Reals,
            initialize=141.355,
            units=enth_mass_units,
            doc="Specific enthalpy parameter A1",
        )
        self.enth_mass_solvent_param_A2 = Var(
            within=Reals,
            initialize=4202.07,
            units=enth_mass_units * t_inv_units,
            doc="Specific enthalpy parameter A2",
        )
        self.enth_mass_solvent_param_A3 = Var(
            within=Reals,
            initialize=-0.535,
            units=enth_mass_units * t_inv_units**2,
            doc="Specific enthalpy parameter A3",
        )
        self.enth_mass_solvent_param_A4 = Var(
            within=Reals,
            initialize=0.004,
            units=enth_mass_units * t_inv_units**3,
            doc="Specific enthalpy parameter A4",
        )

        # Enthalpy parameters for NaCl solution from Sparrow (2003): 0 < T < 300 degC
        self.enth_phase_param_A1 = Var(
            within=Reals,
            initialize=0.0005e3,
            units=enth_mass_units_2,
            doc="Solution enthalpy parameter A1",
        )
        self.enth_phase_param_A2 = Var(
            within=Reals,
            initialize=0.0378e3,
            units=enth_mass_units_2,
            doc="Solution enthalpy parameter A2",
        )
        self.enth_phase_param_A3 = Var(
            within=Reals,
            initialize=-0.3682e3,
            units=enth_mass_units_2,
            doc="Solution enthalpy parameter A3",
        )
        self.enth_phase_param_A4 = Var(
            within=Reals,
            initialize=-0.6529e3,
            units=enth_mass_units_2,
            doc="Solution enthalpy parameter A4",
        )
        self.enth_phase_param_A5 = Var(
            within=Reals,
            initialize=2.89e3,
            units=enth_mass_units_2,
            doc="Solution enthalpy parameter A5",
        )
        self.enth_phase_param_B1 = Var(
            within=Reals,
            initialize=4.145,
            units=enth_mass_units_2 / pyunits.K,
            doc="Solution enthalpy parameter B1",
        )
        self.enth_phase_param_B2 = Var(
            within=Reals,
            initialize=-4.973,
            units=enth_mass_units_2 / pyunits.K,
            doc="Solution enthalpy parameter B2",
        )
        self.enth_phase_param_B3 = Var(
            within=Reals,
            initialize=4.482,
            units=enth_mass_units_2 / pyunits.K,
            doc="Solution enthalpy parameter B3",
        )
        self.enth_phase_param_B4 = Var(
            within=Reals,
            initialize=18.31,
            units=enth_mass_units_2 / pyunits.K,
            doc="Solution enthalpy parameter B4",
        )
        self.enth_phase_param_B5 = Var(
            within=Reals,
            initialize=-46.41,
            units=enth_mass_units_2 / pyunits.K,
            doc="Solution enthalpy parameter B5",
        )
        self.enth_phase_param_C1 = Var(
            within=Reals,
            initialize=0.0007,
            units=enth_mass_units_2 / pyunits.K**2,
            doc="Solution enthalpy parameter C1",
        )
        self.enth_phase_param_C2 = Var(
            within=Reals,
            initialize=-0.0059,
            units=enth_mass_units_2 / pyunits.K**2,
            doc="Solution enthalpy parameter C2",
        )
        self.enth_phase_param_C3 = Var(
            within=Reals,
            initialize=0.0854,
            units=enth_mass_units_2 / pyunits.K**2,
            doc="Solution enthalpy parameter C3",
        )
        self.enth_phase_param_C4 = Var(
            within=Reals,
            initialize=-0.4951,
            units=enth_mass_units_2 / pyunits.K**2,
            doc="Solution enthalpy parameter C4",
        )
        self.enth_phase_param_C5 = Var(
            within=Reals,
            initialize=0.8255,
            units=enth_mass_units_2 / pyunits.K**2,
            doc="Solution enthalpy parameter C5",
        )
        self.enth_phase_param_D1 = Var(
            within=Reals,
            initialize=-0.0048e-3,
            units=enth_mass_units_2 / pyunits.K**3,
            doc="Solution enthalpy parameter D1",
        )
        self.enth_phase_param_D2 = Var(
            within=Reals,
            initialize=0.0639e-3,
            units=enth_mass_units_2 / pyunits.K**3,
            doc="Solution enthalpy parameter D2",
        )
        self.enth_phase_param_D3 = Var(
            within=Reals,
            initialize=-0.714e-3,
            units=enth_mass_units_2 / pyunits.K**3,
            doc="Solution enthalpy parameter D3",
        )
        self.enth_phase_param_D4 = Var(
            within=Reals,
            initialize=3.273e-3,
            units=enth_mass_units_2 / pyunits.K**3,
            doc="Solution enthalpy parameter D4",
        )
        self.enth_phase_param_D5 = Var(
            within=Reals,
            initialize=-4.85e-3,
            units=enth_mass_units_2 / pyunits.K**3,
            doc="Solution enthalpy parameter D5",
        )
        self.enth_phase_param_E1 = Var(
            within=Reals,
            initialize=0.0202e-6,
            units=enth_mass_units_2 / pyunits.K**4,
            doc="Solution enthalpy parameter E1",
        )
        self.enth_phase_param_E2 = Var(
            within=Reals,
            initialize=-0.2432e-6,
            units=enth_mass_units_2 / pyunits.K**4,
            doc="Solution enthalpy parameter E2",
        )
        self.enth_phase_param_E3 = Var(
            within=Reals,
            initialize=2.054e-6,
            units=enth_mass_units_2 / pyunits.K**4,
            doc="Solution enthalpy parameter E3",
        )
        self.enth_phase_param_E4 = Var(
            within=Reals,
            initialize=-8.211e-6,
            units=enth_mass_units_2 / pyunits.K**4,
            doc="Solution enthalpy parameter E4",
        )
        self.enth_phase_param_E5 = Var(
            within=Reals,
            initialize=11.43e-6,
            units=enth_mass_units_2 / pyunits.K**4,
            doc="Solution enthalpy parameter E5",
        )

        for v in self.component_objects(Var):
            v.fix()

        # ---default scaling---
        self.set_default_scaling("temperature", 1e-2)
        self.set_default_scaling("pressure", 1e-6)
        self.set_default_scaling("pressure_sat", 1e-5)
        self.set_default_scaling("dens_mass_solvent", 1e-3, index="Liq")
        self.set_default_scaling("dens_mass_solvent", 1, index="Vap")
        self.set_default_scaling("dens_mass_solute", 1e-3, index="Sol")
        self.set_default_scaling("dens_mass_solute", 1e-3, index="Liq")
        self.set_default_scaling("dens_mass_phase", 1e-3, index="Liq")
        self.set_default_scaling("enth_mass_solvent", 1e-2, index="Liq")
        self.set_default_scaling("enth_mass_solvent", 1e-3, index="Vap")
        self.set_default_scaling("cp_mass_phase", 1e-3, index="Liq")
        self.set_default_scaling("dh_vap_mass_solvent", 1e-3)
        self.set_default_scaling("dh_crystallization_mass_comp", 1e-2, index="NaCl")

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units(
            {
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            }
        )

        obj.add_properties(
            {
                "flow_mass_phase_comp": {"method": None},
                "temperature": {"method": None},
                "pressure": {"method": None},
                "solubility_mass_phase_comp": {"method": "_solubility_mass_phase_comp"},
                "solubility_mass_frac_phase_comp": {
                    "method": "_solubility_mass_frac_phase_comp"
                },
                "mass_frac_phase_comp": {"method": "_mass_frac_phase_comp"},
                "dens_mass_solvent": {"method": "_dens_mass_solvent"},
                "dens_mass_solute": {"method": "_dens_mass_solute"},
                "dens_mass_phase": {"method": "_dens_mass_phase"},
                "dh_vap_mass_solvent": {"method": "_dh_vap_mass_solvent"},
                "cp_mass_solvent": {"method": "_cp_mass_solvent"},
                "cp_mass_solute": {"method": "_cp_mass_solute"},
                "cp_mass_phase": {"method": "_cp_mass_phase"},
                "flow_vol_phase": {"method": "_flow_vol_phase"},
                "flow_vol": {"method": "_flow_vol"},
                "pressure_sat": {"method": "_pressure_sat"},
                "temperature_sat_solvent": {"method": "_temperature_sat_solvent"},
                "conc_mass_phase_comp": {"method": "_conc_mass_phase_comp"},
                "enth_mass_solvent": {"method": "_enth_mass_solvent"},
                "enth_mass_solute": {"method": "_enth_mass_solute"},
                "enth_mass_phase": {"method": "_enth_mass_phase"},
                "dh_crystallization_mass_comp": {
                    "method": "_dh_crystallization_mass_comp"
                },
                "enth_flow": {"method": "_enth_flow"},
                "flow_mol_phase_comp": {"method": "_flow_mol_phase_comp"},
                "mole_frac_phase_comp": {"method": "_mole_frac_phase_comp"},
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
        self._make_state_vars()

    def _make_state_vars(self):
        # Create state variables
        self.pressure = Var(
            domain=NonNegativeReals,
            initialize=101325,
            units=pyunits.Pa,
            doc="State pressure [Pa]",
        )

        self.temperature = Var(
            domain=NonNegativeReals,
            initialize=298.15,
            bounds=(273.15, 373.15),
            units=pyunits.degK,
            doc="State temperature [K]",
        )

        self.flow_mass_phase_comp = Var(
            self.phase_component_set,
            initialize={
                ("Liq", "H2O"): 0.965,
                ("Liq", "NaCl"): 0.035,
                ("Vap", "H2O"): 0,
                ("Sol", "NaCl"): 0,
            },
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.kg / pyunits.s,
            doc="Mass flow rate",
        )

    # Property Methods

    # 1 Mass fraction: From NaCl property package
    def _mass_frac_phase_comp(self):
        self.mass_frac_phase_comp = Var(
            self.phase_component_set,
            domain=NonNegativeReals,
            initialize={
                ("Liq", "H2O"): 0.965,
                ("Liq", "NaCl"): 0.035,
                ("Vap", "H2O"): 1.0,
                ("Sol", "NaCl"): 1.0,
            },
            bounds=(0, 1.0001),
            units=pyunits.dimensionless,
            doc="Mass fraction",
        )

        def rule_mass_frac_phase_comp(b, p, j):
            phase_comp_list = [
                (p, j)
                for j in self.params.component_list
                if (p, j) in b.phase_component_set
            ]
            if len(phase_comp_list) == 1:  # one component in this phase
                return b.mass_frac_phase_comp[p, j] == 1
            else:
                return b.mass_frac_phase_comp[p, j] == b.flow_mass_phase_comp[
                    p, j
                ] / sum(b.flow_mass_phase_comp[p_j] for p_j in phase_comp_list)

        self.eq_mass_frac_phase_comp = Constraint(
            self.phase_component_set, rule=rule_mass_frac_phase_comp
        )

    # 2. Solubility in g/L: calculated from solubility mass fraction
    def _solubility_mass_phase_comp(self):
        self.solubility_mass_phase_comp = Var(
            ["Liq"],
            ["NaCl"],
            domain=NonNegativeReals,
            bounds=(300, 1000),
            initialize=356.5,
            units=pyunits.g / pyunits.L,
            doc="solubility of NaCl in water, g/L",
        )

        def rule_solubility_mass_phase_comp(b, j):
            return b.solubility_mass_phase_comp[
                "Liq", j
            ] == b.solubility_mass_frac_phase_comp["Liq", j] * b.dens_mass_solvent[
                "Liq"
            ] / (
                1 - b.solubility_mass_frac_phase_comp["Liq", j]
            )

        self.eq_solubility_mass_phase_comp = Constraint(
            ["NaCl"], rule=rule_solubility_mass_phase_comp
        )

    # 3. Solubility as mass fraction
    def _solubility_mass_frac_phase_comp(self):
        self.solubility_mass_frac_phase_comp = Var(
            ["Liq"],
            ["NaCl"],
            domain=NonNegativeReals,
            bounds=(0, 1.0001),
            initialize=0.5,
            units=pyunits.dimensionless,
            doc="solubility (as mass fraction) of NaCl in water",
        )

        def rule_solubility_mass_frac_phase_comp(b, j):  # Sparrow (2003)
            t = b.temperature - 273.15 * pyunits.K
            return (
                b.solubility_mass_frac_phase_comp["Liq", j]
                == b.params.sol_param_A1
                + b.params.sol_param_A2 * t
                + b.params.sol_param_A3 * t**2
            )

        self.eq_solubility_mass_frac_phase_comp = Constraint(
            ["NaCl"], rule=rule_solubility_mass_frac_phase_comp
        )

    # 4. Density of solvent (pure water in liquid and vapour phases)
    def _dens_mass_solvent(self):
        self.dens_mass_solvent = Var(
            ["Liq", "Vap"],
            initialize=1e3,
            bounds=(1e-4, 1e4),
            units=pyunits.kg * pyunits.m**-3,
            doc="Mass density of pure water",
        )

        def rule_dens_mass_solvent(b, p):
            if p == "Liq":  # density, eq. 8 in Sharqawy
                t = b.temperature - 273.15 * pyunits.K
                dens_mass_w = (
                    b.params.dens_mass_param_A1
                    + b.params.dens_mass_param_A2 * t
                    + b.params.dens_mass_param_A3 * t**2
                    + b.params.dens_mass_param_A4 * t**3
                    + b.params.dens_mass_param_A5 * t**4
                )
                return b.dens_mass_solvent[p] == dens_mass_w
            elif p == "Vap":
                return b.dens_mass_solvent[p] == (
                    b.params.mw_comp["H2O"] * b.pressure
                ) / (Constants.gas_constant * b.temperature)

        self.eq_dens_mass_solvent = Constraint(
            ["Liq", "Vap"], rule=rule_dens_mass_solvent
        )

    # 5. Density of NaCl crystals and liquid
    def _dens_mass_solute(self):
        self.dens_mass_solute = Var(
            ["Sol", "Liq"],
            initialize=1e3,
            bounds=(1e-4, 1e4),
            units=pyunits.kg * pyunits.m**-3,
            doc="Mass density of solid NaCl crystals",
        )

        def rule_dens_mass_solute(b, p):
            if p == "Sol":
                return b.dens_mass_solute[p] == b.params.dens_mass_param_NaCl
            elif p == "Liq":  # Apparent density in eq. 9 of Laliberte paper
                t = b.temperature - 273.15 * pyunits.K
                v_app = (
                    b.mass_frac_phase_comp["Liq", "NaCl"]
                    + b.params.dens_mass_param_NaCl_liq_C2
                    + (b.params.dens_mass_param_NaCl_liq_C3 * t)
                ) / (
                    (
                        b.mass_frac_phase_comp["Liq", "NaCl"]
                        * b.params.dens_mass_param_NaCl_liq_C0
                    )
                    + b.params.dens_mass_param_NaCl_liq_C1
                )
                v_app = v_app / exp(
                    0.000001
                    * pyunits.K**-2
                    * (t + b.params.dens_mass_param_NaCl_liq_C4) ** 2
                )
                return b.dens_mass_solute[p] == 1 / v_app

        self.eq_dens_mass_solute = Constraint(
            ["Sol", "Liq"], rule=rule_dens_mass_solute
        )

    # 6. Density of liquid solution (Water + NaCl)
    def _dens_mass_phase(self):
        self.dens_mass_phase = Var(
            ["Liq"],
            initialize=1e3,
            bounds=(5e2, 1e4),
            units=pyunits.kg * pyunits.m**-3,
            doc="Mass density of liquid NaCl solution",
        )

        def rule_dens_mass_phase(b):  # density, eq. 6 of Laliberte paper
            return b.dens_mass_phase["Liq"] == 1 / (
                (b.mass_frac_phase_comp["Liq", "NaCl"] / b.dens_mass_solute["Liq"])
                + (b.mass_frac_phase_comp["Liq", "H2O"] / b.dens_mass_solvent["Liq"])
            )

        self.eq_dens_mass_phase = Constraint(rule=rule_dens_mass_phase)

    # 7. Latent heat of vapourization of pure water
    def _dh_vap_mass_solvent(self):
        self.dh_vap_mass_solvent = Var(
            initialize=2.4e3,
            bounds=(1, 1e5),
            units=pyunits.kJ / pyunits.kg,
            doc="Latent heat of vaporization of pure water",
        )

        def rule_dh_vap_mass_solvent(b):
            t = b.temperature - 273.15 * pyunits.K
            dh_vap_sol = (
                b.params.dh_vap_w_param_0
                + b.params.dh_vap_w_param_1 * t
                + b.params.dh_vap_w_param_2 * t**2
                + b.params.dh_vap_w_param_3 * t**3
                + b.params.dh_vap_w_param_4 * t**4
            )
            return b.dh_vap_mass_solvent == pyunits.convert(
                dh_vap_sol, to_units=pyunits.kJ / pyunits.kg
            )

        self.eq_dh_vap_mass_solvent = Constraint(rule=rule_dh_vap_mass_solvent)

    # 8. Heat capacity of solvent (pure water in liquid and vapour phases)
    def _cp_mass_solvent(self):
        self.cp_mass_solvent = Var(
            ["Liq", "Vap"],
            initialize=4e3,
            bounds=(1e-5, 1e5),
            units=pyunits.J / pyunits.kg / pyunits.K,
            doc="Specific heat capacity of pure solvent",
        )

        def rule_cp_mass_solvent(b, p):
            if p == "Liq":
                # specific heat, eq. 9 in Sharqawy et al. (2010)
                # Convert T90 to T68, eq. 4 in Sharqawy et al. (2010); primary reference from Rusby (1991)
                t = (b.temperature - 0.00025 * 273.15 * pyunits.K) / (1 - 0.00025)
                A = b.params.cp_phase_param_A1
                B = b.params.cp_phase_param_B1
                C = b.params.cp_phase_param_C1
                D = b.params.cp_phase_param_D1
                return (
                    b.cp_mass_solvent["Liq"]
                    == (A + B * t + C * t**2 + D * t**3) * 1000
                )
            elif p == "Vap":
                t = b.temperature / 1000
                return (
                    b.cp_mass_solvent["Vap"]
                    == b.params.cp_vap_param_A
                    + b.params.cp_vap_param_B * t
                    + b.params.cp_vap_param_C * t**2
                    + b.params.cp_vap_param_D * t**3
                    + b.params.cp_vap_param_E / t**2
                )

        self.eq_cp_mass_solvent = Constraint(["Liq", "Vap"], rule=rule_cp_mass_solvent)

    # 9. Heat capacity of solid-phase NaCl crystals
    def _cp_mass_solute(self):
        self.cp_mass_solute = Var(
            ["Liq", "Sol"],
            initialize=1e3,
            bounds=(-1e4, 1e5),
            units=pyunits.J / pyunits.kg / pyunits.K,
            doc="Specific heat capacity of solid NaCl crystals",
        )

        def rule_cp_mass_solute(b, p):
            if p == "Sol":  # Shomate equation for NaCl, NIST
                t = b.temperature / (1000 * pyunits.dimensionless)
                cp_mass_solute_mol = (
                    b.params.cp_param_NaCl_solid_A
                    + (b.params.cp_param_NaCl_solid_B * t)
                    + (b.params.cp_param_NaCl_solid_C * t**2)
                    + (b.params.cp_param_NaCl_solid_D * t**3)
                    + (b.params.cp_param_NaCl_solid_E / (t**2))
                )
                return (
                    b.cp_mass_solute[p] == cp_mass_solute_mol / b.params.mw_comp["NaCl"]
                )
            if (
                p == "Liq"
            ):  # NaCl liq. apparent specific heat capacity, eq. 11-12 of Laliberte (2009)
                t = b.temperature - 273.15 * pyunits.K
                alpha = (
                    (b.params.cp_param_NaCl_liq_A2 * t)
                    + (
                        b.params.cp_param_NaCl_liq_A4
                        * (1 - b.mass_frac_phase_comp["Liq", "H2O"])
                    )
                    + (b.params.cp_param_NaCl_liq_A3 * exp(0.01 * pyunits.K**-1 * t))
                )
                cp_nacl_liq = b.params.cp_param_NaCl_liq_A1 * exp(
                    alpha
                ) + b.params.cp_param_NaCl_liq_A5 * (
                    (1 - b.mass_frac_phase_comp["Liq", "H2O"])
                    ** b.params.cp_param_NaCl_liq_A6
                )
                return b.cp_mass_solute[p] == pyunits.convert(
                    cp_nacl_liq, to_units=pyunits.J / pyunits.kg / pyunits.K
                )

        self.eq_cp_mass_solute = Constraint(["Liq", "Sol"], rule=rule_cp_mass_solute)

    # 10. cp of liquid solution (Water + NaCl)
    def _cp_mass_phase(self):
        self.cp_mass_phase = Var(
            ["Liq"],
            initialize=4e3,
            bounds=(1e-4, 1e5),
            units=pyunits.J / pyunits.kg / pyunits.K,
            doc="Specific heat capacity of liquid solution",
        )

        def rule_cp_mass_phase(b):  # heat capacity, eq. 10 of Laliberte (2009) paper
            return (
                b.cp_mass_phase["Liq"]
                == b.mass_frac_phase_comp["Liq", "NaCl"] * b.cp_mass_solute["Liq"]
                + b.mass_frac_phase_comp["Liq", "H2O"] * b.cp_mass_solvent["Liq"]
            )

        self.eq_cp_mass_phase = Constraint(rule=rule_cp_mass_phase)

    # 11. Volumetric flow rate for each phase
    def _flow_vol_phase(self):
        self.flow_vol_phase = Var(
            self.params.phase_list,
            initialize=1,
            bounds=(0, None),
            units=pyunits.m**3 / pyunits.s,
            doc="Volumetric flow rate",
        )

        def rule_flow_vol_phase(b, p):
            if p == "Liq":
                return (
                    b.flow_vol_phase[p]
                    == sum(
                        b.flow_mass_phase_comp[p, j]
                        for j in self.params.component_list
                        if (p, j) in self.phase_component_set
                    )
                    / b.dens_mass_phase[p]
                )
            elif p == "Sol":
                return (
                    b.flow_vol_phase[p]
                    == sum(
                        b.flow_mass_phase_comp[p, j]
                        for j in self.params.component_list
                        if (p, j) in self.phase_component_set
                    )
                    / b.dens_mass_solute["Sol"]
                )
            elif p == "Vap":
                return (
                    b.flow_vol_phase[p]
                    == sum(
                        b.flow_mass_phase_comp[p, j]
                        for j in self.params.component_list
                        if (p, j) in self.phase_component_set
                    )
                    / b.dens_mass_solvent["Vap"]
                )

        self.eq_flow_vol_phase = Constraint(
            self.params.phase_list, rule=rule_flow_vol_phase
        )

    # 12. Total volumetric flow rate
    def _flow_vol(self):
        def rule_flow_vol(b):
            return sum(b.flow_vol_phase[p] for p in self.params.phase_list)

        self.flow_vol = Expression(rule=rule_flow_vol)

    # 13. Vapour pressure of the NaCl solution based on the boiling temperature
    def _pressure_sat(self):
        self.pressure_sat = Var(
            initialize=1e3,
            bounds=(0.001, 1e6),
            units=pyunits.Pa,
            doc="Vapor pressure of NaCl solution",
        )

        def rule_pressure_sat(b):  # vapor pressure, eq6 in Sparrow (2003)
            t = b.temperature - 273.15 * pyunits.K
            x = b.mass_frac_phase_comp["Liq", "NaCl"]

            ps_a = (
                b.params.pressure_sat_param_A1
                + (b.params.pressure_sat_param_A2 * x)
                + (b.params.pressure_sat_param_A3 * x**2)
                + (b.params.pressure_sat_param_A4 * x**3)
                + (b.params.pressure_sat_param_A5 * x**4)
            )

            ps_b = (
                b.params.pressure_sat_param_B1
                + (b.params.pressure_sat_param_B2 * x)
                + (b.params.pressure_sat_param_B3 * x**2)
                + (b.params.pressure_sat_param_B4 * x**3)
                + (b.params.pressure_sat_param_B5 * x**4)
            )

            ps_c = (
                b.params.pressure_sat_param_C1
                + (b.params.pressure_sat_param_C2 * x)
                + (b.params.pressure_sat_param_C3 * x**2)
                + (b.params.pressure_sat_param_C4 * x**3)
                + (b.params.pressure_sat_param_C5 * x**4)
            )

            ps_d = (
                b.params.pressure_sat_param_D1
                + (b.params.pressure_sat_param_D2 * x)
                + (b.params.pressure_sat_param_D3 * x**2)
                + (b.params.pressure_sat_param_D4 * x**3)
                + (b.params.pressure_sat_param_D5 * x**4)
            )

            ps_e = (
                b.params.pressure_sat_param_E1
                + (b.params.pressure_sat_param_E2 * x)
                + (b.params.pressure_sat_param_E3 * x**2)
                + (b.params.pressure_sat_param_E4 * x**3)
                + (b.params.pressure_sat_param_E5 * x**4)
            )

            p_sat = (
                ps_a + (ps_b * t) + (ps_c * t**2) + (ps_d * t**3) + (ps_e * t**4)
            )
            return b.pressure_sat == pyunits.convert(p_sat, to_units=pyunits.Pa)

        self.eq_pressure_sat = Constraint(rule=rule_pressure_sat)

    # 14. Saturation temperature for water vapour at calculated boiling pressure
    def _temperature_sat_solvent(self):
        self.temperature_sat_solvent = Var(
            initialize=298.15,
            bounds=(273.15, 1000.15),
            units=pyunits.K,
            doc="Vapour (saturation) temperature of pure solvent at boiling (i.e. crystallization) pressure",
        )

        def rule_temperature_sat_solvent(b):
            psat = pyunits.convert(b.pressure_sat, to_units=pyunits.kPa)
            return (
                b.temperature_sat_solvent
                == b.params.temp_sat_solvent_A1
                + b.params.temp_sat_solvent_A2
                / (
                    log(psat / b.params.temp_sat_solvent_A3)
                    + b.params.temp_sat_solvent_A4
                )
            )

        self.eq_temperature_sat_solvent = Constraint(rule=rule_temperature_sat_solvent)

    # 15. Mass concentration
    def _conc_mass_phase_comp(self):
        self.conc_mass_phase_comp = Var(
            ["Liq"],
            self.params.component_list,
            initialize=10,
            bounds=(0, 1e6),
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

    # 16. Specific enthalpy of solvent (pure water in liquid and vapour phases)
    def _enth_mass_solvent(self):
        self.enth_mass_solvent = Var(
            ["Liq", "Vap"],
            initialize=1e3,
            bounds=(1, 1e4),
            units=pyunits.kJ * pyunits.kg**-1,
            doc="Specific saturated enthalpy of pure solvent",
        )

        def rule_enth_mass_solvent(b, p):
            t = b.temperature - 273.15 * pyunits.K
            h_w = (
                b.params.enth_mass_solvent_param_A1
                + b.params.enth_mass_solvent_param_A2 * t
                + b.params.enth_mass_solvent_param_A3 * t**2
                + b.params.enth_mass_solvent_param_A4 * t**3
            )
            if p == "Liq":  # enthalpy, eq. 55 in Sharqawy
                return b.enth_mass_solvent[p] == pyunits.convert(
                    h_w, to_units=pyunits.kJ * pyunits.kg**-1
                )
            elif p == "Vap":

                return (
                    b.enth_mass_solvent[p]
                    == pyunits.convert(h_w, to_units=pyunits.kJ * pyunits.kg**-1)
                    + +b.dh_vap_mass_solvent
                )

        self.eq_enth_mass_solvent = Constraint(
            ["Liq", "Vap"], rule=rule_enth_mass_solvent
        )

    # 17. Specific enthalpy of NaCl solution
    def _enth_mass_phase(self):
        self.enth_mass_phase = Var(
            ["Liq"],
            initialize=500,
            bounds=(1, 1000),
            units=pyunits.kJ * pyunits.kg**-1,
            doc="Specific enthalpy of NaCl solution",
        )

        def rule_enth_mass_phase(
            b,
        ):  # specific enthalpy calculation based on Sparrow (2003).
            t = (
                b.temperature - 273.15 * pyunits.K
            )  # temperature in degC, but pyunits in K
            S = b.mass_frac_phase_comp["Liq", "NaCl"]

            enth_a = (
                b.params.enth_phase_param_A1
                + (b.params.enth_phase_param_A2 * S)
                + (b.params.enth_phase_param_A3 * S**2)
                + (b.params.enth_phase_param_A4 * S**3)
                + (b.params.enth_phase_param_A5 * S**4)
            )

            enth_b = (
                b.params.enth_phase_param_B1
                + (b.params.enth_phase_param_B2 * S)
                + (b.params.enth_phase_param_B3 * S**2)
                + (b.params.enth_phase_param_B4 * S**3)
                + (b.params.enth_phase_param_B5 * S**4)
            )

            enth_c = (
                b.params.enth_phase_param_C1
                + (b.params.enth_phase_param_C2 * S)
                + (b.params.enth_phase_param_C3 * S**2)
                + (b.params.enth_phase_param_C4 * S**3)
                + (b.params.enth_phase_param_C5 * S**4)
            )

            enth_d = (
                b.params.enth_phase_param_D1
                + (b.params.enth_phase_param_D2 * S)
                + (b.params.enth_phase_param_D3 * S**2)
                + (b.params.enth_phase_param_D4 * S**3)
                + (b.params.enth_phase_param_D5 * S**4)
            )

            enth_e = (
                b.params.enth_phase_param_E1
                + (b.params.enth_phase_param_E2 * S)
                + (b.params.enth_phase_param_E3 * S**2)
                + (b.params.enth_phase_param_E4 * S**3)
                + (b.params.enth_phase_param_E5 * S**4)
            )

            return b.enth_mass_phase["Liq"] == enth_a + (enth_b * t) + (
                enth_c * t**2
            ) + (enth_d * t**3) + (enth_e * t**4)

        self.eq_enth_mass_phase = Constraint(rule=rule_enth_mass_phase)

    # 18. Heat of crystallization
    def _dh_crystallization_mass_comp(self):
        self.dh_crystallization_mass_comp = Var(
            ["NaCl"],
            initialize=1,
            bounds=(-1e3, 1e3),
            units=pyunits.kJ / pyunits.kg,
            doc="NaCl heat of crystallization",
        )

        def rule_dh_crystallization_mass_comp(b):
            if (
                b.params.config.heat_of_crystallization_model
                == HeatOfCrystallizationModel.constant
            ):
                return (
                    b.dh_crystallization_mass_comp["NaCl"]
                    == b.params.dh_crystallization_param
                )
            elif (
                b.params.config.heat_of_crystallization_model
                == HeatOfCrystallizationModel.zero
            ):
                return (
                    b.dh_crystallization_mass_comp["NaCl"]
                    == 0 * pyunits.kJ / pyunits.kg
                )
            elif (
                b.params.config.heat_of_crystallization_model
                == HeatOfCrystallizationModel.temp_dependent
            ):
                raise NotImplementedError(
                    f"Temperature-dependent heat of crystallization model has not been implemented yet."
                )

        self.eq_dh_crystallization_mass_comp = Constraint(
            rule=rule_dh_crystallization_mass_comp
        )

    # 19. Heat capacity of solid-phase NaCl crystals
    def _enth_mass_solute(self):
        self.enth_mass_solute = Var(
            ["Sol"],
            initialize=1e3,
            bounds=(1e-3, 1e4),
            units=pyunits.kJ / pyunits.kg,
            doc="Specific enthalpy of solid NaCl crystals",
        )

        def rule_enth_mass_solute(b, p):
            ############################
            # Shomate equation for molar enthalpy ofNaCl, NIST
            # Note: Tref is 298 K, so changing the Tref to 273 K to match IAPWS is necessary.
            # Computation formula for reference temperature change:
            #    Enthalpy at T relative to 273 K = Enthalpy change relative to 298 K + (Enthalpy at 298 K - Enthalpy at 273 K)
            ############################

            # (i) Enthalpy at original reference temperature (298 K)
            enth_mass_solute_mol_tref_1 = 0  # Enthalpy at original temperature (298 K)

            # (ii) Enthalpy at new reference temperature (273.15 K)
            tref_2 = 273.15 * pyunits.K / 1000
            enth_mass_solute_mol_tref_2 = (
                (b.params.cp_param_NaCl_solid_A * tref_2)
                + ((1 / 2) * b.params.cp_param_NaCl_solid_B * tref_2**2)
                + ((1 / 3) * b.params.cp_param_NaCl_solid_C * tref_2**3)
                + ((1 / 4) * b.params.cp_param_NaCl_solid_D * tref_2**4)
                - (b.params.cp_param_NaCl_solid_E / tref_2)
                + b.params.cp_param_NaCl_solid_F
                - b.params.cp_param_NaCl_solid_H
            )

            # (iii) Compute enthalpy change at temperature t relative to old reference temperature t_ref_1
            t = b.temperature / (1000 * pyunits.dimensionless)
            dh_mass_solute_mol_tref_1 = (
                (b.params.cp_param_NaCl_solid_A * t)
                + ((1 / 2) * b.params.cp_param_NaCl_solid_B * t**2)
                + ((1 / 3) * b.params.cp_param_NaCl_solid_C * t**3)
                + ((1 / 4) * b.params.cp_param_NaCl_solid_D * t**4)
                - (b.params.cp_param_NaCl_solid_E / t)
                + b.params.cp_param_NaCl_solid_F
                - b.params.cp_param_NaCl_solid_H
            )

            # (iv) Compute enthalpy change at temperature t relative to new reference temperature t_ref_2
            enth_mass_solute_mol = dh_mass_solute_mol_tref_1 + (
                enth_mass_solute_mol_tref_1 - enth_mass_solute_mol_tref_2
            )

            # (v) Convert from molar enthalpy to mass enthalpy
            enth_mass_solute_mol = enth_mass_solute_mol / b.params.mw_comp["NaCl"]
            return b.enth_mass_solute[p] == enth_mass_solute_mol * (
                pyunits.kJ / pyunits.J
            )

        self.eq_enth_mass_solute = Constraint(["Sol"], rule=rule_enth_mass_solute)

    # 20. Total enthalpy flow for any stream: adds up the enthalpies for the solid, liquid and vapour phases
    # Assumes no NaCl is vapour stream or water in crystals
    def _enth_flow(self):
        # enthalpy flow expression for get_enthalpy_flow_terms method

        def rule_enth_flow(b):  # enthalpy flow [J/s]
            return (
                sum(b.flow_mass_phase_comp["Liq", j] for j in b.params.component_list)
                * b.enth_mass_phase["Liq"]
                + b.flow_mass_phase_comp["Vap", "H2O"] * b.enth_mass_solvent["Vap"]
                + b.flow_mass_phase_comp["Sol", "NaCl"] * b.enth_mass_solute["Sol"]
            )

        self.enth_flow = Expression(rule=rule_enth_flow)

    # 21. Molar flows
    def _flow_mol_phase_comp(self):
        self.flow_mol_phase_comp = Var(
            self.phase_component_set,
            initialize=100,
            bounds=(None, None),
            domain=NonNegativeReals,
            units=pyunits.mol / pyunits.s,
            doc="Molar flowrate",
        )

        def rule_flow_mol_phase_comp(b, p, j):
            return (
                b.flow_mol_phase_comp[p, j]
                == b.flow_mass_phase_comp[p, j] / b.params.mw_comp[j]
            )

        self.eq_flow_mol_phase_comp = Constraint(
            self.phase_component_set, rule=rule_flow_mol_phase_comp
        )

    # 22. Mole fractions
    def _mole_frac_phase_comp(self):
        self.mole_frac_phase_comp = Var(
            self.phase_component_set,
            initialize=0.1,
            bounds=(0, 1.0001),
            units=pyunits.dimensionless,
            doc="Mole fraction",
        )

        def rule_mole_frac_phase_comp(b, p, j):
            phase_comp_list = [
                (p, j)
                for j in self.params.component_list
                if (p, j) in b.phase_component_set
            ]
            if len(phase_comp_list) == 1:  # one component in this phase
                return b.mole_frac_phase_comp[p, j] == 1
            else:
                return b.mole_frac_phase_comp[p, j] == b.flow_mol_phase_comp[
                    p, j
                ] / sum(b.flow_mol_phase_comp[p_j] for (p_j) in phase_comp_list)

        self.eq_mole_frac_phase_comp = Constraint(
            self.phase_component_set, rule=rule_mole_frac_phase_comp
        )

    # -----------------------------------------------------------------------------
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
        return {
            "flow_mass_phase_comp": self.flow_mass_phase_comp,
            "temperature": self.temperature,
            "pressure": self.pressure,
        }

    # -----------------------------------------------------------------------------
    # Scaling methods
    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # default scaling factors have already been set with idaes.core.property_base.calculate_scaling_factors()
        # for the following variables: flow_mass_phase_comp, pressure, temperature, dens_mass_phase, enth_mass_phase

        # These variables should have user input
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

        if iscale.get_scaling_factor(self.flow_mass_phase_comp["Sol", "NaCl"]) is None:
            sf = iscale.get_scaling_factor(
                self.flow_mass_phase_comp["Sol", "NaCl"], default=1e2, warning=True
            )
            iscale.set_scaling_factor(self.flow_mass_phase_comp["Sol", "NaCl"], sf)

        if iscale.get_scaling_factor(self.flow_mass_phase_comp["Vap", "H2O"]) is None:
            sf = iscale.get_scaling_factor(
                self.flow_mass_phase_comp["Vap", "H2O"], default=1e0, warning=True
            )
            iscale.set_scaling_factor(self.flow_mass_phase_comp["Vap", "H2O"], sf)

        # scaling factors for molecular weights
        for j, v in self.params.mw_comp.items():
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(self.params.mw_comp, 1e3)

        # Scaling for solubility (g/L) parameters. Values typically about 300-500, so scale by 1e-3.
        if self.is_property_constructed("solubility_mass_phase_comp"):
            if iscale.get_scaling_factor(self.solubility_mass_phase_comp) is None:
                iscale.set_scaling_factor(self.solubility_mass_phase_comp, 1e-3)

        # Scaling for solubility mass fraction. Values typically about 0-1, so scale by 1e0.
        if self.is_property_constructed("solubility_mass_frac_phase_comp"):
            if iscale.get_scaling_factor(self.solubility_mass_frac_phase_comp) is None:
                iscale.set_scaling_factor(self.solubility_mass_frac_phase_comp, 1e0)

        # Scaling for flow_vol_phase: scaled as scale of dominant component in phase / density of phase
        if self.is_property_constructed("flow_vol_phase"):
            for p in self.params.phase_list:
                if p == "Liq":
                    if iscale.get_scaling_factor(self.flow_vol_phase[p]) is None:
                        sf = iscale.get_scaling_factor(
                            self.flow_mass_phase_comp[p, "H2O"]
                        ) / iscale.get_scaling_factor(self.dens_mass_phase[p])
                        iscale.set_scaling_factor(self.flow_vol_phase[p], sf)
                elif p == "Vap":
                    if iscale.get_scaling_factor(self.flow_vol_phase[p]) is None:
                        sf = iscale.get_scaling_factor(
                            self.flow_mass_phase_comp[p, "H2O"]
                        ) / iscale.get_scaling_factor(self.dens_mass_solvent[p])
                        iscale.set_scaling_factor(self.flow_vol_phase[p], sf)
                elif p == "Sol":
                    if iscale.get_scaling_factor(self.flow_vol_phase[p]) is None:
                        sf = iscale.get_scaling_factor(
                            self.flow_mass_phase_comp[p, "NaCl"]
                        ) / iscale.get_scaling_factor(self.dens_mass_solute[p])
                        iscale.set_scaling_factor(self.flow_vol_phase[p], sf)

        # Scaling for flow volume
        if self.is_property_constructed("flow_vol"):
            sf = iscale.get_scaling_factor(self.flow_vol_phase)
            iscale.set_scaling_factor(self.flow_vol, sf)

        # Scaling material heat capacities
        if self.is_property_constructed("cp_mass_solute"):
            for p in ["Sol", "Liq"]:
                if iscale.get_scaling_factor(self.cp_mass_solute[p]) is None:
                    iscale.set_scaling_factor(
                        self.cp_mass_solute[p], 1e-3
                    )  # same as scaling factor of .cp_mass_phase['Liq']

        if self.is_property_constructed("cp_mass_solvent"):
            for p in ["Liq", "Vap"]:
                if iscale.get_scaling_factor(self.cp_mass_solvent[p]) is None:
                    iscale.set_scaling_factor(
                        self.cp_mass_solvent[p], 1e-3
                    )  # same as scaling factor of .cp_mass_phase['Liq']

        # Scaling saturation temperature
        if self.is_property_constructed("temperature_sat_solvent"):
            if iscale.get_scaling_factor(self.temperature_sat_solvent) is None:
                iscale.set_scaling_factor(
                    self.temperature_sat_solvent,
                    iscale.get_scaling_factor(self.temperature),
                )

        # Scaling solute and solvent enthalpies
        if self.is_property_constructed("enth_mass_solute"):
            if iscale.get_scaling_factor(self.enth_mass_solute) is None:
                iscale.set_scaling_factor(
                    self.enth_mass_solute["Sol"],
                    iscale.get_scaling_factor(self.enth_mass_solvent["Liq"]),
                )

        if self.is_property_constructed("enth_mass_phase"):
            if iscale.get_scaling_factor(self.enth_mass_phase) is None:
                iscale.set_scaling_factor(
                    self.enth_mass_phase["Liq"],
                    iscale.get_scaling_factor(self.enth_mass_solvent["Liq"]),
                )

        # Scaling enthapy flow - not sure about this one
        if self.is_property_constructed("enth_flow"):
            iscale.set_scaling_factor(
                self.enth_flow,
                iscale.get_scaling_factor(self.flow_mass_phase_comp["Liq", "H2O"])
                * iscale.get_scaling_factor(self.enth_mass_phase["Liq"]),
            )

        # Scaling molar flows - derived from flow_mass
        if self.is_property_constructed("flow_mol_phase_comp"):
            for p, j in self.phase_component_set:
                if iscale.get_scaling_factor(self.flow_mol_phase_comp[p, j]) is None:
                    sf = iscale.get_scaling_factor(self.flow_mass_phase_comp[p, j])
                    sf /= iscale.get_scaling_factor(self.params.mw_comp[j])
                    iscale.set_scaling_factor(self.flow_mol_phase_comp[p, j], sf)

        ######################################################
        # Scaling for mass fractions - needs verification!
        if self.is_property_constructed("mass_frac_phase_comp"):
            # Option 1:
            for p, j in self.phase_component_set:
                if iscale.get_scaling_factor(self.mass_frac_phase_comp[p, j]) is None:
                    if p == "Sol":
                        iscale.set_scaling_factor(self.mass_frac_phase_comp[p, j], 1e0)
                    else:
                        if j == "NaCl":
                            sf = iscale.get_scaling_factor(
                                self.flow_mass_phase_comp[p, j]
                            ) / iscale.get_scaling_factor(
                                self.flow_mass_phase_comp[p, "H2O"]
                            )
                            iscale.set_scaling_factor(
                                self.mass_frac_phase_comp[p, j], sf
                            )
                        elif j == "H2O":
                            iscale.set_scaling_factor(
                                self.mass_frac_phase_comp[p, j], 1e0
                            )

        # Scaling for mole fractions - same approach as mass fractions - needs verification!
        # Appears to make things worse!
        if self.is_property_constructed("mole_frac_phase_comp"):
            # Option 1:
            for p, j in self.phase_component_set:
                if iscale.get_scaling_factor(self.mole_frac_phase_comp[p, j]) is None:
                    if p == "Sol":
                        iscale.set_scaling_factor(self.mole_frac_phase_comp[p, j], 1e-1)
                    else:
                        if j == "NaCl":
                            sf = iscale.get_scaling_factor(
                                self.flow_mol_phase_comp[p, j]
                            ) / iscale.get_scaling_factor(
                                self.flow_mol_phase_comp[p, "H2O"]
                            )
                            iscale.set_scaling_factor(
                                self.mole_frac_phase_comp[p, j], sf
                            )
                        elif j == "H2O":
                            iscale.set_scaling_factor(
                                self.mole_frac_phase_comp[p, j], 1e0
                            )

        #     ########################################################

        # Scaling for mass concentrations
        if self.is_property_constructed("conc_mass_phase_comp"):
            for j in self.params.component_list:
                sf_dens = iscale.get_scaling_factor(self.dens_mass_phase["Liq"])
                if (
                    iscale.get_scaling_factor(self.conc_mass_phase_comp["Liq", j])
                    is None
                ):
                    if j == "H2O":
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

        # Transforming constraints
        # property relationships with no index, simple constraint
        v_str_lst_simple = [
            "pressure_sat",
            "dh_vap_mass_solvent",
            "temperature_sat_solvent",
        ]
        for v_str in v_str_lst_simple:
            if self.is_property_constructed(v_str):
                v = getattr(self, v_str)
                sf = iscale.get_scaling_factor(v, default=1, warning=True)
                c = getattr(self, "eq_" + v_str)
                iscale.constraint_scaling_transform(c, sf)

        # Property relationships with phase index, but simple constraint
        v_str_lst_phase = ["dens_mass_phase", "enth_mass_phase", "cp_mass_phase"]
        for v_str in v_str_lst_phase:
            if self.is_property_constructed(v_str):
                v = getattr(self, v_str)
                sf = iscale.get_scaling_factor(v["Liq"], default=1, warning=True)
                c = getattr(self, "eq_" + v_str)
                iscale.constraint_scaling_transform(c, sf)

        # Property relationship indexed by component
        v_str_lst_comp = [
            "solubility_mass_phase_comp",
            "solubility_mass_frac_phase_comp",
            "conc_mass_phase_comp",
        ]
        for v_str in v_str_lst_comp:
            if self.is_property_constructed(v_str):
                v_comp = getattr(self, v_str)
                c_comp = getattr(self, "eq_" + v_str)
                for j, c in c_comp.items():
                    sf = iscale.get_scaling_factor(
                        v_comp["Liq", j], default=1, warning=True
                    )
                    iscale.constraint_scaling_transform(c, sf)

        # Property relationship indexed by single component
        if self.is_property_constructed("dh_crystallization_mass_comp"):
            sf = iscale.get_scaling_factor(self.dh_crystallization_mass_comp["NaCl"])
            iscale.constraint_scaling_transform(
                self.eq_dh_crystallization_mass_comp, sf
            )

        # Property relationships with phase index and indexed constraints
        v_str_lst_phase = [
            "dens_mass_solvent",
            "dens_mass_solute",
            "flow_vol_phase",
            "enth_mass_solvent",
            "enth_mass_solute",
            "cp_mass_solvent",
            "cp_mass_solute",
        ]
        for v_str in v_str_lst_phase:
            if self.is_property_constructed(v_str):
                v = getattr(self, v_str)
                c_phase = getattr(self, "eq_" + v_str)
                for (ind, c) in c_phase.items():
                    sf = iscale.get_scaling_factor(v[ind], default=1, warning=True)
                    iscale.constraint_scaling_transform(c, sf)

        # Property relationships indexed by component and phase
        v_str_lst_phase_comp = [
            "mass_frac_phase_comp",
            "flow_mol_phase_comp",
            "mole_frac_phase_comp",
        ]
        for v_str in v_str_lst_phase_comp:
            if self.is_property_constructed(v_str):
                v_comp = getattr(self, v_str)
                c_comp = getattr(self, "eq_" + v_str)
                for j, c in c_comp.items():
                    sf = iscale.get_scaling_factor(v_comp[j], default=1, warning=True)
                    iscale.constraint_scaling_transform(c, sf)
