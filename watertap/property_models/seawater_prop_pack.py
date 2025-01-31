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
Initial property package for seawater system
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
    log,
    log10,
    exp,
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
from watertap.core.solvers import get_solver
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


@declare_process_block_class("SeawaterParameterBlock")
class SeawaterParameterData(PhysicalParameterBlock):
    """Parameter block for a seawater property package."""

    CONFIG = PhysicalParameterBlock.CONFIG()

    def build(self):
        """
        Callable method for Block construction.
        """
        super(SeawaterParameterData, self).build()

        self._state_block_class = SeawaterStateBlock

        # components
        self.H2O = Solvent()
        self.TDS = Solute()

        # phases
        self.Liq = LiquidPhase()

        """ References
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
        # TODO: Edit comment above about diffusivity when/if relationship is changed

        # parameters
        # molecular weight
        mw_comp_data = {
            "H2O": 18.01528e-3,
            "TDS": 31.4038218e-3,
        }
        # molecular weight of TDS is taken as the average atomic weight of sea salt, based on
        # "Reference-Composition Salinity Scale" in Millero et al. (2008) and cited by Sharqawy et al. (2010)

        self.mw_comp = Param(
            self.component_list,
            initialize=mw_comp_data,
            units=pyunits.kg / pyunits.mol,
            doc="Molecular weight",
        )

        # mass density parameters, 0-180 C, 0-150 g/kg, 0-12 MPa
        # eq. 8 in Sharqawy et al. (2010)
        dens_units = pyunits.kg / pyunits.m**3
        t_inv_units = pyunits.K**-1
        s_inv_units = pyunits.kg / pyunits.g

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
        self.dens_mass_param_B1 = Var(
            within=Reals,
            initialize=8.020e2,
            units=dens_units,
            doc="Mass density parameter B1",
        )
        self.dens_mass_param_B2 = Var(
            within=Reals,
            initialize=-2.001,
            units=dens_units * t_inv_units,
            doc="Mass density parameter B2",
        )
        self.dens_mass_param_B3 = Var(
            within=Reals,
            initialize=1.677e-2,
            units=dens_units * t_inv_units**2,
            doc="Mass density parameter B3",
        )
        self.dens_mass_param_B4 = Var(
            within=Reals,
            initialize=-3.060e-5,
            units=dens_units * t_inv_units**3,
            doc="Mass density parameter B4",
        )
        self.dens_mass_param_B5 = Var(
            within=Reals,
            initialize=-1.613e-5,
            units=dens_units * t_inv_units**2,
            doc="Mass density parameter B5",
        )

        visc_d_units = pyunits.Pa * pyunits.s
        # dynamic viscosity parameters, 0-180 C, 0-150 g/kg
        # eq. 22 and 23 in Sharqawy et al. (2010)
        self.visc_d_param_muw_A = Var(
            within=Reals,
            initialize=4.2844e-5,
            units=visc_d_units,
            doc="Dynamic viscosity parameter A for pure water",
        )
        self.visc_d_param_muw_B = Var(
            within=Reals,
            initialize=0.157,
            units=t_inv_units**2 * visc_d_units**-1,
            doc="Dynamic viscosity parameter B for pure water",
        )
        self.visc_d_param_muw_C = Var(
            within=Reals,
            initialize=64.993,
            units=pyunits.K,
            doc="Dynamic viscosity parameter C for pure water",
        )
        self.visc_d_param_muw_D = Var(
            within=Reals,
            initialize=91.296,
            units=visc_d_units**-1,
            doc="Dynamic viscosity parameter D for pure water",
        )
        self.visc_d_param_A_1 = Var(
            within=Reals,
            initialize=1.541,
            units=pyunits.dimensionless,
            doc="Dynamic viscosity parameter 1 for term A",
        )
        self.visc_d_param_A_2 = Var(
            within=Reals,
            initialize=1.998e-2,
            units=t_inv_units,
            doc="Dynamic viscosity parameter 2 for term A",
        )
        self.visc_d_param_A_3 = Var(
            within=Reals,
            initialize=-9.52e-5,
            units=t_inv_units**2,
            doc="Dynamic viscosity parameter 3 for term A",
        )
        self.visc_d_param_B_1 = Var(
            within=Reals,
            initialize=7.974,
            units=pyunits.dimensionless,
            doc="Dynamic viscosity parameter 1 for term B",
        )
        self.visc_d_param_B_2 = Var(
            within=Reals,
            initialize=-7.561e-2,
            units=t_inv_units,
            doc="Dynamic viscosity parameter 2 for term B",
        )
        self.visc_d_param_B_3 = Var(
            within=Reals,
            initialize=4.724e-4,
            units=t_inv_units**2,
            doc="Dynamic viscosity parameter 3 for term B",
        )

        # diffusivity parameters, 25 C
        # eq. 6 in Bartholomew & Mauter (2019)
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
            doc="Diffusivity parameters",
        )

        # osmotic coefficient parameters, 0-200 C, 0-120 g/kg
        # eq. 49 in Sharqawy et al. (2010)
        self.osm_coeff_param_1 = Var(
            within=Reals,
            initialize=8.9453e-1,
            units=pyunits.dimensionless,
            doc="Osmotic coefficient parameter 1",
        )
        self.osm_coeff_param_2 = Var(
            within=Reals,
            initialize=4.1561e-4,
            units=t_inv_units,
            doc="Osmotic coefficient parameter 2",
        )
        self.osm_coeff_param_3 = Var(
            within=Reals,
            initialize=-4.6262e-6,
            units=t_inv_units**2,
            doc="Osmotic coefficient parameter 3",
        )
        self.osm_coeff_param_4 = Var(
            within=Reals,
            initialize=2.2211e-11,
            units=t_inv_units**4,
            doc="Osmotic coefficient parameter 4",
        )
        self.osm_coeff_param_5 = Var(
            within=Reals,
            initialize=-1.1445e-1,
            units=pyunits.dimensionless,
            doc="Osmotic coefficient parameter 5",
        )
        self.osm_coeff_param_6 = Var(
            within=Reals,
            initialize=-1.4783e-3,
            units=t_inv_units,
            doc="Osmotic coefficient parameter 6",
        )
        self.osm_coeff_param_7 = Var(
            within=Reals,
            initialize=-1.3526e-8,
            units=t_inv_units**3,
            doc="Osmotic coefficient parameter 7",
        )
        self.osm_coeff_param_8 = Var(
            within=Reals,
            initialize=7.0132,
            units=pyunits.dimensionless,
            doc="Osmotic coefficient parameter 8",
        )
        self.osm_coeff_param_9 = Var(
            within=Reals,
            initialize=5.696e-2,
            units=t_inv_units,
            doc="Osmotic coefficient parameter 9",
        )
        self.osm_coeff_param_10 = Var(
            within=Reals,
            initialize=-2.8624e-4,
            units=t_inv_units**2,
            doc="Osmotic coefficient parameter 10",
        )

        # specific enthalpy parameters, 10-120 C, 0-120 g/kg, 0-12 MPa
        # Table 9 in Nayar et al. (2016)
        enth_mass_units = pyunits.J / pyunits.kg
        P_inv_units = pyunits.MPa**-1

        self.enth_mass_param_A1 = Var(
            within=Reals,
            initialize=996.7767,
            units=enth_mass_units * P_inv_units,
            doc="Specific enthalpy parameter A1",
        )
        self.enth_mass_param_A2 = Var(
            within=Reals,
            initialize=-3.2406,
            units=enth_mass_units * P_inv_units * t_inv_units,
            doc="Specific enthalpy parameter A2",
        )
        self.enth_mass_param_A3 = Var(
            within=Reals,
            initialize=0.0127,
            units=enth_mass_units * P_inv_units * t_inv_units**2,
            doc="Specific enthalpy parameter A3",
        )
        self.enth_mass_param_A4 = Var(
            within=Reals,
            initialize=-4.7723e-5,
            units=enth_mass_units * P_inv_units * t_inv_units**3,
            doc="Specific enthalpy parameter A4",
        )
        self.enth_mass_param_A5 = Var(
            within=Reals,
            initialize=-1.1748,
            units=enth_mass_units * P_inv_units,
            doc="Specific enthalpy parameter A5",
        )
        self.enth_mass_param_A6 = Var(
            within=Reals,
            initialize=0.01169,
            units=enth_mass_units * P_inv_units * t_inv_units,
            doc="Specific enthalpy parameter A6",
        )
        self.enth_mass_param_A7 = Var(
            within=Reals,
            initialize=-2.6185e-5,
            units=enth_mass_units * P_inv_units * t_inv_units**2,
            doc="Specific enthalpy parameter A7",
        )
        self.enth_mass_param_A8 = Var(
            within=Reals,
            initialize=7.0661e-8,
            units=enth_mass_units * P_inv_units * t_inv_units**3,
            doc="Specific enthalpy parameter A8",
        )
        self.enth_mass_param_B1 = Var(
            within=Reals,
            initialize=-2.34825e4,
            units=enth_mass_units,
            doc="Specific enthalpy parameter B1",
        )
        self.enth_mass_param_B2 = Var(
            within=Reals,
            initialize=3.15183e5,
            units=enth_mass_units,
            doc="Specific enthalpy parameter B2",
        )
        self.enth_mass_param_B3 = Var(
            within=Reals,
            initialize=2.80269e6,
            units=enth_mass_units,
            doc="Specific enthalpy parameter B3",
        )
        self.enth_mass_param_B4 = Var(
            within=Reals,
            initialize=-1.44606e7,
            units=enth_mass_units,
            doc="Specific enthalpy parameter B4",
        )
        self.enth_mass_param_B5 = Var(
            within=Reals,
            initialize=7.82607e3,
            units=enth_mass_units * t_inv_units,
            doc="Specific enthalpy parameter B5",
        )
        self.enth_mass_param_B6 = Var(
            within=Reals,
            initialize=-4.41733,
            units=enth_mass_units * t_inv_units**2,
            doc="Specific enthalpy parameter B6",
        )
        self.enth_mass_param_B7 = Var(
            within=Reals,
            initialize=2.1394e-1,
            units=enth_mass_units * t_inv_units**3,
            doc="Specific enthalpy parameter B7",
        )
        self.enth_mass_param_B8 = Var(
            within=Reals,
            initialize=-1.99108e4,
            units=enth_mass_units * t_inv_units,
            doc="Specific enthalpy parameter B8",
        )
        self.enth_mass_param_B9 = Var(
            within=Reals,
            initialize=2.77846e4,
            units=enth_mass_units * t_inv_units,
            doc="Specific enthalpy parameter B9",
        )
        self.enth_mass_param_B10 = Var(
            within=Reals,
            initialize=9.72801,
            units=enth_mass_units * t_inv_units**2,
            doc="Specific enthalpy parameter B10",
        )
        self.enth_mass_param_C1 = Var(
            within=Reals,
            initialize=141.355,
            units=enth_mass_units,
            doc="Specific enthalpy parameter C1",
        )
        self.enth_mass_param_C2 = Var(
            within=Reals,
            initialize=4202.07,
            units=enth_mass_units * t_inv_units,
            doc="Specific enthalpy parameter C2",
        )
        self.enth_mass_param_C3 = Var(
            within=Reals,
            initialize=-0.535,
            units=enth_mass_units * t_inv_units**2,
            doc="Specific enthalpy parameter C3",
        )
        self.enth_mass_param_C4 = Var(
            within=Reals,
            initialize=0.004,
            units=enth_mass_units * t_inv_units**3,
            doc="Specific enthalpy parameter C4",
        )

        # vapor pressure parameters,  0-180 C, 0-160 g/kg
        # eq. 5 and 6 in Nayar et al.(2016)
        self.pressure_sat_param_psatw_A1 = Var(
            within=Reals,
            initialize=-5.8002206e3,
            units=pyunits.K,
            doc="Vapor pressure of pure water parameter A1",
        )
        self.pressure_sat_param_psatw_A2 = Var(
            within=Reals,
            initialize=1.3914993,
            units=pyunits.dimensionless,
            doc="Vapor pressure of pure water parameter A2",
        )
        self.pressure_sat_param_psatw_A3 = Var(
            within=Reals,
            initialize=-4.8640239e-2,
            units=t_inv_units,
            doc="Vapor pressure of pure water parameter A3",
        )
        self.pressure_sat_param_psatw_A4 = Var(
            within=Reals,
            initialize=4.1764768e-5,
            units=t_inv_units**2,
            doc="Vapor pressure of pure water parameter A4",
        )
        self.pressure_sat_param_psatw_A5 = Var(
            within=Reals,
            initialize=-1.4452093e-8,
            units=t_inv_units**3,
            doc="Vapor pressure of pure water parameter A5",
        )
        self.pressure_sat_param_psatw_A6 = Var(
            within=Reals,
            initialize=6.5459673,
            units=pyunits.dimensionless,
            doc="Vapor pressure of pure water parameter A6",
        )
        self.pressure_sat_param_B1 = Var(
            within=Reals,
            initialize=-4.5818e-4,
            units=s_inv_units,
            doc="Vapor pressure of seawater parameter B1",
        )
        self.pressure_sat_param_B2 = Var(
            within=Reals,
            initialize=-2.0443e-6,
            units=s_inv_units**2,
            doc="Vapor pressure of seawater parameter B2",
        )

        # specific heat parameters, 0-180 C, 0-180 g/kg, 0-12 MPa
        # eq. 9 in Sharqawy et al. (2010)
        cp_units = pyunits.J / (pyunits.kg * pyunits.K)
        self.cp_phase_param_A1 = Var(
            within=Reals,
            initialize=5.328,
            units=cp_units,
            doc="Specific heat of seawater parameter A1",
        )
        self.cp_phase_param_A2 = Var(
            within=Reals,
            initialize=-9.76e-2,
            units=cp_units * s_inv_units,
            doc="Specific heat of seawater parameter A2",
        )
        self.cp_phase_param_A3 = Var(
            within=Reals,
            initialize=4.04e-4,
            units=cp_units * s_inv_units**2,
            doc="Specific heat of seawater parameter A3",
        )
        self.cp_phase_param_B1 = Var(
            within=Reals,
            initialize=-6.913e-3,
            units=cp_units * t_inv_units,
            doc="Specific heat of seawater parameter B1",
        )
        self.cp_phase_param_B2 = Var(
            within=Reals,
            initialize=7.351e-4,
            units=cp_units * s_inv_units * t_inv_units,
            doc="Specific heat of seawater parameter B2",
        )
        self.cp_phase_param_B3 = Var(
            within=Reals,
            initialize=-3.15e-6,
            units=cp_units * s_inv_units**2 * t_inv_units,
            doc="Specific heat of seawater parameter B3",
        )
        self.cp_phase_param_C1 = Var(
            within=Reals,
            initialize=9.6e-6,
            units=cp_units * t_inv_units**2,
            doc="Specific heat of seawater parameter C1",
        )
        self.cp_phase_param_C2 = Var(
            within=Reals,
            initialize=-1.927e-6,
            units=cp_units * s_inv_units * t_inv_units**2,
            doc="Specific heat of seawater parameter C2",
        )
        self.cp_phase_param_C3 = Var(
            within=Reals,
            initialize=8.23e-9,
            units=cp_units * s_inv_units**2 * t_inv_units**2,
            doc="Specific heat of seawater parameter C3",
        )
        self.cp_phase_param_D1 = Var(
            within=Reals,
            initialize=2.5e-9,
            units=cp_units * t_inv_units**3,
            doc="Specific heat of seawater parameter D1",
        )
        self.cp_phase_param_D2 = Var(
            within=Reals,
            initialize=1.666e-9,
            units=cp_units * s_inv_units * t_inv_units**3,
            doc="Specific heat of seawater parameter D2",
        )
        self.cp_phase_param_D3 = Var(
            within=Reals,
            initialize=-7.125e-12,
            units=cp_units * s_inv_units**2 * t_inv_units**3,
            doc="Specific heat of seawater parameter D3",
        )

        # thermal conductivity parameters, 0-180 C, 0-160 g/kg
        # eq. 13 in Sharqawy et al. (2010)
        self.therm_cond_phase_param_1 = Var(
            within=Reals,
            initialize=240,
            units=pyunits.dimensionless,
            doc="Thermal conductivity of seawater parameter 1",
        )
        self.therm_cond_phase_param_2 = Var(
            within=Reals,
            initialize=0.0002,
            units=s_inv_units,
            doc="Thermal conductivity of seawater parameter 2",
        )
        self.therm_cond_phase_param_3 = Var(
            within=Reals,
            initialize=0.434,
            units=pyunits.dimensionless,
            doc="Thermal conductivity of seawater parameter 3",
        )
        self.therm_cond_phase_param_4 = Var(
            within=Reals,
            initialize=2.3,
            units=pyunits.dimensionless,
            doc="Thermal conductivity of seawater parameter 4",
        )
        self.therm_cond_phase_param_5 = Var(
            within=Reals,
            initialize=343.5,
            units=t_inv_units**-1,
            doc="Thermal conductivity of seawater parameter 5",
        )
        self.therm_cond_phase_param_6 = Var(
            within=Reals,
            initialize=0.037,
            units=s_inv_units * t_inv_units**-1,
            doc="Thermal conductivity of seawater parameter 6",
        )
        self.therm_cond_phase_param_7 = Var(
            within=Reals,
            initialize=647,
            units=t_inv_units**-1,
            doc="Thermal conductivity of seawater parameter 7",
        )
        self.therm_cond_phase_param_8 = Var(
            within=Reals,
            initialize=0.03,
            units=t_inv_units**-1 * s_inv_units,
            doc="Thermal conductivity of seawater parameter 8",
        )

        # latent heat of pure water parameters, 0-200 C
        # eq. 54 in Sharqawy et al. (2010)
        self.dh_vap_w_param_0 = Var(
            within=Reals,
            initialize=2.501e6,
            units=enth_mass_units,
            doc="Latent heat of pure water parameter 0",
        )
        self.dh_vap_w_param_1 = Var(
            within=Reals,
            initialize=-2.369e3,
            units=cp_units,
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

        # Boiling point elevation parameters, 0-200 C, 0-120 g/kg
        # eq. 36 in Sharqawy et al. (2010)
        self.bpe_A0 = Var(
            within=Reals,
            initialize=17.95,
            units=pyunits.K,
            doc="Boiling point parameter A0",
        )
        self.bpe_A1 = Var(
            within=Reals,
            initialize=2.823e-1,
            units=pyunits.dimensionless,
            doc="Boiling point parameter A1",
        )
        self.bpe_A2 = Var(
            within=Reals,
            initialize=-4.584e-4,
            units=t_inv_units,
            doc="Boiling point parameter A2",
        )
        self.bpe_B0 = Var(
            within=Reals,
            initialize=6.56,
            units=pyunits.K,
            doc="Boiling point parameter B0",
        )
        self.bpe_B1 = Var(
            within=Reals,
            initialize=5.267e-2,
            units=pyunits.dimensionless,
            doc="Boiling point parameter B1",
        )
        self.bpe_B2 = Var(
            within=Reals,
            initialize=1.536e-4,
            units=t_inv_units,
            doc="Boiling point parameter B2",
        )

        # traditional parameters are the only Vars currently on the block and should be fixed
        for v in self.component_objects(Var):
            v.fix()

        # ---default scaling---
        self.set_default_scaling("temperature", 1e-2)
        self.set_default_scaling("pressure", 1e-6)
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
                "dens_mass_solvent": {"method": "_dens_mass_solvent"},
                "osm_coeff": {"method": "_osm_coeff"},
                "enth_flow": {"method": "_enth_flow"},
                "dh_vap_mass": {"method": "_dh_vap_mass"},
                "boiling_point_elevation_phase": {
                    "method": "_boiling_point_elevation_phase"
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


class _SeawaterStateBlock(StateBlock):
    """
    This Class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """

    def fix_initialization_states(self):
        """
        Fixes state variables for state blocks.

        Returns:
            None
        """
        # Fix state variables
        fix_state_vars(self)

        # Constraint on water concentration at outlet - unfix in these cases
        for b in self.values():
            if b.config.defined_state is False:
                b.conc_mol_comp["H2O"].unfix()

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
                    "State vars fixed but degrees of "
                    "freedom for state block is not "
                    "zero during initialization."
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
            fix_state_vars(self)

        return results


@declare_process_block_class("SeawaterStateBlock", block_class=_SeawaterStateBlock)
class SeawaterStateBlockData(StateBlockData):
    """A seawater property package."""

    def build(self):
        """Callable method for Block construction."""
        super().build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        # Add state variables
        self.flow_mass_phase_comp = Var(
            self.params.phase_list,
            self.params.component_list,
            initialize={("Liq", "H2O"): 0.965, ("Liq", "TDS"): 0.035},
            bounds=(0.0, None),
            domain=NonNegativeReals,
            units=pyunits.kg / pyunits.s,
            doc="Mass flow rate",
        )

        self.temperature = Var(
            initialize=298.15,
            bounds=(273.15, 1000),
            domain=NonNegativeReals,
            units=pyunits.K,
            doc="Temperature",
        )

        self.pressure = Var(
            initialize=101325,
            bounds=(1e3, 5e7),
            domain=NonNegativeReals,
            units=pyunits.Pa,
            doc="Pressure",
        )

    # -----------------------------------------------------------------------------
    # Property Methods
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
            return b.flow_mass_phase_comp[p, j] == b.mass_frac_phase_comp[p, j] * sum(
                b.flow_mass_phase_comp[p, j] for j in b.params.component_list
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
            bounds=(1, 1e6),
            units=pyunits.kg * pyunits.m**-3,
            doc="Mass density of seawater",
        )

        # Sharqawy et al. (2010), eq. 8, 0-180 C, 0-150 g/kg, 0-12 MPa
        def rule_dens_mass_phase(b, p):
            t = b.temperature - 273.15 * pyunits.K
            s = b.mass_frac_phase_comp[p, "TDS"]
            dens_mass = (
                b.dens_mass_solvent
                + b.params.dens_mass_param_B1 * s
                + b.params.dens_mass_param_B2 * s * t
                + b.params.dens_mass_param_B3 * s * t**2
                + b.params.dens_mass_param_B4 * s * t**3
                + b.params.dens_mass_param_B5 * s**2 * t**2
            )
            return b.dens_mass_phase[p] == dens_mass

        self.eq_dens_mass_phase = Constraint(
            self.params.phase_list, rule=rule_dens_mass_phase
        )

    def _dens_mass_solvent(self):
        self.dens_mass_solvent = Var(
            initialize=1e3,
            bounds=(1, 1e6),
            units=pyunits.kg * pyunits.m**-3,
            doc="Mass density of pure water",
        )

        # Sharqawy et al. (2010), eq. 8, 0-180 C
        def rule_dens_mass_solvent(b):
            t = b.temperature - 273.15 * pyunits.K
            dens_mass_w = (
                b.params.dens_mass_param_A1
                + b.params.dens_mass_param_A2 * t
                + b.params.dens_mass_param_A3 * t**2
                + b.params.dens_mass_param_A4 * t**3
                + b.params.dens_mass_param_A5 * t**4
            )
            return b.dens_mass_solvent == dens_mass_w

        self.eq_dens_mass_solvent = Constraint(rule=rule_dens_mass_solvent)

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
                == sum(b.flow_mass_phase_comp[p, j] for j in b.params.component_list)
                / b.dens_mass_phase[p]
            )

        self.eq_flow_vol_phase = Constraint(
            self.params.phase_list, rule=rule_flow_vol_phase
        )

    def _flow_vol(self):
        def rule_flow_vol(b):
            return sum(b.flow_vol_phase[p] for p in b.params.phase_list)

        self.flow_vol = Expression(rule=rule_flow_vol)

    def _conc_mass_phase_comp(self):
        self.conc_mass_phase_comp = Var(
            self.params.phase_list,
            self.params.component_list,
            initialize=10,
            bounds=(0.0, 1e6),
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
            return b.flow_mol_phase_comp[p, j] == b.mole_frac_phase_comp[p, j] * sum(
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
            ["TDS"],
            initialize=1,
            bounds=(0.0, 1e6),
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
            self.params.phase_list, ["TDS"], rule=rule_molality_phase_comp
        )

    def _visc_d_phase(self):
        self.visc_d_phase = Var(
            self.params.phase_list,
            initialize=1e-3,
            bounds=(0.0, 1),
            units=pyunits.Pa * pyunits.s,
            doc="Viscosity",
        )

        # Sharqawy et al. (2010), eq. 22 and 23, 0-180 C, 0-150 g/kg
        def rule_visc_d_phase(b, p):
            # temperature in degC, but pyunits are K
            t = b.temperature - 273.15 * pyunits.K
            s = b.mass_frac_phase_comp[p, "TDS"]
            mu_w = (
                b.params.visc_d_param_muw_A
                + (
                    b.params.visc_d_param_muw_B * (t + b.params.visc_d_param_muw_C) ** 2
                    - b.params.visc_d_param_muw_D
                )
                ** -1
            )
            A = (
                b.params.visc_d_param_A_1
                + b.params.visc_d_param_A_2 * t
                + b.params.visc_d_param_A_3 * t**2
            )
            B = (
                b.params.visc_d_param_B_1
                + b.params.visc_d_param_B_2 * t
                + b.params.visc_d_param_B_3 * t**2
            )
            return b.visc_d_phase[p] == mu_w * (1 + A * s + B * s**2)

        self.eq_visc_d_phase = Constraint(
            self.params.phase_list, rule=rule_visc_d_phase
        )

    # TODO: diffusivity from NaCl prop model used temporarily--reconsider this
    def _diffus_phase_comp(self):
        self.diffus_phase_comp = Var(
            self.params.phase_list,
            ["TDS"],
            initialize=1e-9,
            bounds=(1e-10, 1e-8),
            units=pyunits.m**2 * pyunits.s**-1,
            doc="Diffusivity",
        )

        # Bartholomew & Mauter (2019), eq. 6 (substituting NaCl w/ TDS), 25 C
        def rule_diffus_phase_comp(b, p, j):
            return b.diffus_phase_comp[p, j] == (
                b.params.diffus_param["4"] * b.mass_frac_phase_comp[p, j] ** 4
                + b.params.diffus_param["3"] * b.mass_frac_phase_comp[p, j] ** 3
                + b.params.diffus_param["2"] * b.mass_frac_phase_comp[p, j] ** 2
                + b.params.diffus_param["1"] * b.mass_frac_phase_comp[p, j]
                + b.params.diffus_param["0"]
            )

        self.eq_diffus_phase_comp = Constraint(
            self.params.phase_list, ["TDS"], rule=rule_diffus_phase_comp
        )

    def _osm_coeff(self):
        self.osm_coeff = Var(
            initialize=1,
            bounds=(0.0, 10),
            units=pyunits.dimensionless,
            doc="Osmotic coefficient",
        )

        # Sharqawy et al. (2010), eq. 49, 0-200 C, 0-120 g/kg
        def rule_osm_coeff(b):
            s = b.mass_frac_phase_comp["Liq", "TDS"]
            # temperature in degC, but pyunits are still K
            t = b.temperature - 273.15 * pyunits.K
            osm_coeff = (
                b.params.osm_coeff_param_1
                + b.params.osm_coeff_param_2 * t
                + b.params.osm_coeff_param_3 * t**2
                + b.params.osm_coeff_param_4 * t**4
                + b.params.osm_coeff_param_5 * s
                + b.params.osm_coeff_param_6 * s * t
                + b.params.osm_coeff_param_7 * s * t**3
                + b.params.osm_coeff_param_8 * s**2
                + b.params.osm_coeff_param_9 * s**2 * t
                + b.params.osm_coeff_param_10 * s**2 * t**2
            )
            return b.osm_coeff == osm_coeff

        self.eq_osm_coeff = Constraint(rule=rule_osm_coeff)

    def _pressure_osm_phase(self):
        self.pressure_osm_phase = Var(
            self.params.phase_list,
            initialize=1e6,
            bounds=(1, 1e8),
            units=pyunits.Pa,
            doc="Osmotic pressure",
        )

        # Nayar et al. (2016), eq. 48, 0-200 C, 0-120 g/kg
        def rule_pressure_osm_phase(b, p):
            i = 2  # number of ionic species
            rhow = b.dens_mass_solvent
            return (
                b.pressure_osm_phase[p]
                == b.osm_coeff
                * b.molality_phase_comp[p, "TDS"]
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
            initialize=1e6,
            bounds=(1, 1e9),
            units=pyunits.J * pyunits.kg**-1,
            doc="Specific enthalpy",
        )

        # Nayar et al. (2016), eq. 25 and 26, 10-120 C, 0-120 g/kg, 0-12 MPa
        def rule_enth_mass_phase(b, p):
            # temperature in degC, but pyunits in K
            t = b.temperature - 273.15 * pyunits.K
            S_kg_kg = b.mass_frac_phase_comp[p, "TDS"]
            S_g_kg = S_kg_kg * 1000
            P = b.pressure - 101325 * pyunits.Pa
            P_MPa = pyunits.convert(P, to_units=pyunits.MPa)

            h_w = (
                b.params.enth_mass_param_C1
                + b.params.enth_mass_param_C2 * t
                + b.params.enth_mass_param_C3 * t**2
                + b.params.enth_mass_param_C4 * t**3
            )
            h_sw0 = h_w - S_kg_kg * (
                b.params.enth_mass_param_B1
                + b.params.enth_mass_param_B2 * S_kg_kg
                + b.params.enth_mass_param_B3 * S_kg_kg**2
                + b.params.enth_mass_param_B4 * S_kg_kg**3
                + b.params.enth_mass_param_B5 * t
                + b.params.enth_mass_param_B6 * t**2
                + b.params.enth_mass_param_B7 * t**3
                + b.params.enth_mass_param_B8 * S_kg_kg * t
                + b.params.enth_mass_param_B9 * S_kg_kg**2 * t
                + b.params.enth_mass_param_B10 * S_kg_kg * t**2
            )
            h_sw = h_sw0 + P_MPa * (
                b.params.enth_mass_param_A1
                + b.params.enth_mass_param_A2 * t
                + b.params.enth_mass_param_A3 * t**2
                + b.params.enth_mass_param_A4 * t**3
                + S_g_kg
                * (
                    +b.params.enth_mass_param_A5
                    + b.params.enth_mass_param_A6 * t
                    + b.params.enth_mass_param_A7 * t**2
                    + b.params.enth_mass_param_A8 * t**3
                )
            )
            return b.enth_mass_phase[p] == h_sw

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

    def _pressure_sat(self):
        self.pressure_sat = Var(
            initialize=1e3, bounds=(1, 1e8), units=pyunits.Pa, doc="Vapor pressure"
        )

        # Nayar et al.(2016), eq. 5 and 6, 0-180 C, 0-160 g/kg
        def rule_pressure_sat(b):
            t = b.temperature
            s = b.mass_frac_phase_comp["Liq", "TDS"] * 1000 * pyunits.g / pyunits.kg
            psatw = (
                exp(
                    b.params.pressure_sat_param_psatw_A1 * t**-1
                    + b.params.pressure_sat_param_psatw_A2
                    + b.params.pressure_sat_param_psatw_A3 * t
                    + b.params.pressure_sat_param_psatw_A4 * t**2
                    + b.params.pressure_sat_param_psatw_A5 * t**3
                    + b.params.pressure_sat_param_psatw_A6 * log(t / pyunits.K)
                )
                * pyunits.Pa
            )
            return b.pressure_sat == psatw * exp(
                b.params.pressure_sat_param_B1 * s
                + b.params.pressure_sat_param_B2 * s**2
            )

        self.eq_pressure_sat = Constraint(rule=rule_pressure_sat)

    def _cp_mass_phase(self):
        self.cp_mass_phase = Var(
            self.params.phase_list,
            initialize=4e3,
            bounds=(0.0, 1e8),
            units=pyunits.J / pyunits.kg / pyunits.K,
            doc="Specific heat capacity",
        )

        # Sharqawy et al. (2010), eq. 9, 0-180 C, 0-180 g/kg, 0-12 MPa
        def rule_cp_mass_phase(b, p):
            # Convert T90 to T68, eq. 4 in Sharqawy et al. (2010); primary reference from Rusby (1991)
            t = (b.temperature - 0.00025 * 273.15 * pyunits.K) / (1 - 0.00025)
            s = b.mass_frac_phase_comp[p, "TDS"] * 1000 * pyunits.g / pyunits.kg
            A = (
                b.params.cp_phase_param_A1
                + b.params.cp_phase_param_A2 * s
                + b.params.cp_phase_param_A3 * s**2
            )
            B = (
                b.params.cp_phase_param_B1
                + b.params.cp_phase_param_B2 * s
                + b.params.cp_phase_param_B3 * s**2
            )
            C = (
                b.params.cp_phase_param_C1
                + b.params.cp_phase_param_C2 * s
                + b.params.cp_phase_param_C3 * s**2
            )
            D = (
                b.params.cp_phase_param_D1
                + b.params.cp_phase_param_D2 * s
                + b.params.cp_phase_param_D3 * s**2
            )
            return b.cp_mass_phase[p] == (A + B * t + C * t**2 + D * t**3) * 1000

        self.eq_cp_mass_phase = Constraint(
            self.params.phase_list, rule=rule_cp_mass_phase
        )

    def _therm_cond_phase(self):
        self.therm_cond_phase = Var(
            self.params.phase_list,
            initialize=0.6,
            bounds=(0.0, 1),
            units=pyunits.W / pyunits.m / pyunits.K,
            doc="Thermal conductivity",
        )

        # Sharqawy  et al. (2010), eq. 13, 0-180 C, 0-160 g/kg
        def rule_therm_cond_phase(b, p):
            # Convert T90 to T68, eq. 4 in Sharqawy et al. (2010); primary reference from Rusby (1991)
            t = (b.temperature - 0.00025 * 273.15 * pyunits.K) / (1 - 0.00025)
            s = b.mass_frac_phase_comp[p, "TDS"] * 1000 * pyunits.g / pyunits.kg
            log10_ksw = log10(
                b.params.therm_cond_phase_param_1
                + b.params.therm_cond_phase_param_2 * s
            ) + b.params.therm_cond_phase_param_3 * (
                b.params.therm_cond_phase_param_4
                - (
                    b.params.therm_cond_phase_param_5
                    + b.params.therm_cond_phase_param_6 * s
                )
                / t
            ) * (
                1
                - t
                / (
                    b.params.therm_cond_phase_param_7
                    + b.params.therm_cond_phase_param_8 * s
                )
            ) ** (
                1 / 3
            )
            return (
                b.therm_cond_phase[p]
                == 10**log10_ksw * 1e-3 * pyunits.W / pyunits.m / pyunits.K
            )

        self.eq_therm_cond_phase = Constraint(
            self.params.phase_list, rule=rule_therm_cond_phase
        )

    def _dh_vap_mass(self):
        self.dh_vap_mass = Var(
            initialize=2.4e3,
            bounds=(1, 1e9),
            units=pyunits.J / pyunits.kg,
            doc="Latent heat of vaporization",
        )

        # Sharqawy et al. (2010), eq. 37 and 54, 0-200 C, 0-240 g/kg
        def rule_dh_vap_mass(b):
            t = b.temperature - 273.15 * pyunits.K
            s = b.mass_frac_phase_comp["Liq", "TDS"]
            dh_vap_mass_w = (
                b.params.dh_vap_w_param_0
                + b.params.dh_vap_w_param_1 * t
                + b.params.dh_vap_w_param_2 * t**2
                + b.params.dh_vap_w_param_3 * t**3
                + b.params.dh_vap_w_param_4 * t**4
            )
            return b.dh_vap_mass == dh_vap_mass_w * (1 - s)

        self.eq_dh_vap_mass = Constraint(rule=rule_dh_vap_mass)

    def _boiling_point_elevation_phase(self):
        self.boiling_point_elevation_phase = Var(
            self.params.phase_list,
            initialize=5e-1,
            bounds=(0, 1e9),
            units=pyunits.K,
            doc="Boiling point elevation",
        )

        # Sharqawy et al. (2010), eq. 36, 0-200 C, 0-120 g/kg
        def rule_boiling_point_elevation_phase(b, p):
            t = b.temperature - 273.15 * pyunits.K
            s = b.mass_frac_phase_comp["Liq", "TDS"]
            A = b.params.bpe_A0 + b.params.bpe_A1 * t + b.params.bpe_A2 * t**2
            B = b.params.bpe_B0 + b.params.bpe_B1 * t + b.params.bpe_B2 * t**2
            return b.boiling_point_elevation_phase[p] == A * s**2 + B * s

        self.eq_boiling_point_elevation_phase = Constraint(
            self.params.phase_list, rule=rule_boiling_point_elevation_phase
        )

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
        # temperature, dens_mass_phase, visc_d_phase, osm_coeff, and enth_mass_phase

        # these variables should have user input
        if iscale.get_scaling_factor(self.flow_mass_phase_comp["Liq", "H2O"]) is None:
            sf = iscale.get_scaling_factor(
                self.flow_mass_phase_comp["Liq", "H2O"], default=1e0, warning=True
            )
            iscale.set_scaling_factor(self.flow_mass_phase_comp["Liq", "H2O"], sf)

        if iscale.get_scaling_factor(self.flow_mass_phase_comp["Liq", "TDS"]) is None:
            sf = iscale.get_scaling_factor(
                self.flow_mass_phase_comp["Liq", "TDS"], default=1e2, warning=True
            )
            iscale.set_scaling_factor(self.flow_mass_phase_comp["Liq", "TDS"], sf)

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
                    if j == "TDS":
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
                            self.mass_frac_phase_comp["Liq", j], 1
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
                    elif j == "TDS":
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
