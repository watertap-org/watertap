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
Initial property package for H2O-NaCl system with temperature dependence
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

        # molecular weight
        mw_comp_data = {"H2O": 18.01528e-3, "NaCl": 58.44e-3}
        self.mw_comp = Param(
            self.component_list,
            initialize=extract_data(mw_comp_data),
            units=pyunits.kg / pyunits.mol,
            doc="Molecular weight kg/mol",
        )

        # solubility_comp, 0-450 C
        # Sparrow 2003, Eq 5: Xsat = param_0 + param_1 * T + param_2 * T**2
        solubility_comp_param_dict = {"0": 0.2628, "1": 62.75e-6, "2": 1.084e-6}
        self.solubility_comp_param = Var(
            solubility_comp_param_dict.keys(),
            domain=Reals,
            initialize=solubility_comp_param_dict,
            units=pyunits.dimensionless,
            doc="solubility_comp parameters",
        )

        # mass density parameters, 0-300 C
        # Sparrow 2003, Eq 7: density = A+BT+CT2+DT3+ET4, where A = ao + a1X + a2X2 + a3X3 + a4X4
        dens_mass_A_param_dict = {
            "0": 1.0001e3,
            "1": 0.7666e3,
            "2": -0.0149e3,
            "3": 0.2663e3,
            "4": 0.8845e3,
        }
        dens_mass_B_param_dict = {
            "0": -0.0214,
            "1": -3.496,
            "2": 10.02,
            "3": -6.56,
            "4": -31.37,
        }
        dens_mass_C_param_dict = {
            "0": -5.263e-3,
            "1": 39.87e-3,
            "2": -176.2e-3,
            "3": 363.5e-3,
            "4": -7.784e-3,
        }
        dens_mass_D_param_dict = {
            "0": 15.42e-6,
            "1": -167e-6,
            "2": 980.7e-6,
            "3": -2573e-6,
            "4": 876.6e-6,
        }
        dens_mass_E_param_dict = {
            "0": -0.0276e-6,
            "1": 0.2978e-6,
            "2": -2.017e-6,
            "3": 6.345e-6,
            "4": -3.914e-6,
        }

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
            units=pyunits.kg / pyunits.m**3,
            doc="Mass density parameter B",
        )
        self.dens_mass_param_C = Var(
            dens_mass_C_param_dict.keys(),
            domain=Reals,
            initialize=dens_mass_C_param_dict,
            units=pyunits.kg / pyunits.m**3,
            doc="Mass density parameter C",
        )
        self.dens_mass_param_D = Var(
            dens_mass_D_param_dict.keys(),
            domain=Reals,
            initialize=dens_mass_D_param_dict,
            units=pyunits.kg / pyunits.m**3,
            doc="Mass density parameter D",
        )
        self.dens_mass_param_E = Var(
            dens_mass_E_param_dict.keys(),
            domain=Reals,
            initialize=dens_mass_E_param_dict,
            units=pyunits.kg / pyunits.m**3,
            doc="Mass density parameter E",
        )

        # specific enthalpy parameters, 0-300 C
        # Sparrow 2003, Eq 8: h = A+BT+CT2+DT3+ET4, where A = ao + a1X + a2X2 + a3X3 + a4X4
        enth_A_param_dict = {
            "0": 0.0005e3,
            "1": 0.0378e3,
            "2": -0.3682e3,
            "3": -0.6529e3,
            "4": 2.89e3,
        }
        enth_B_param_dict = {
            "0": 4.145,
            "1": -4.973,
            "2": 4.482,
            "3": 18.31,
            "4": -46.41,
        }
        enth_C_param_dict = {
            "0": 0.0007,
            "1": -0.0059,
            "2": 0.0854,
            "3": -0.495,
            "4": 0.8255,
        }
        enth_D_param_dict = {
            "0": -0.0048e-3,
            "1": 0.0639e-3,
            "2": -0.714e-3,
            "3": 3.273e-3,
            "4": -4.85e-3,
        }
        enth_E_param_dict = {
            "0": 0.0202e-6,
            "1": -0.2432e-6,
            "2": 2.054e-6,
            "3": -8.211e-6,
            "4": 11.43e-6,
        }

        self.enth_param_A = Var(
            enth_A_param_dict.keys(),
            domain=Reals,
            initialize=enth_A_param_dict,
            units=pyunits.J / pyunits.kg,
            doc="Specific enthalpy parameter A",
        )
        self.enth_param_B = Var(
            enth_B_param_dict.keys(),
            domain=Reals,
            initialize=enth_B_param_dict,
            units=pyunits.J / pyunits.kg,
            doc="Specific enthalpy parameter B",
        )
        self.enth_param_C = Var(
            enth_C_param_dict.keys(),
            domain=Reals,
            initialize=enth_C_param_dict,
            units=pyunits.J / pyunits.kg,
            doc="Specific enthalpy parameter C",
        )
        self.enth_param_D = Var(
            enth_D_param_dict.keys(),
            domain=Reals,
            initialize=enth_D_param_dict,
            units=pyunits.J / pyunits.kg,
            doc="Specific enthalpy parameter D",
        )
        self.enth_param_E = Var(
            enth_E_param_dict.keys(),
            domain=Reals,
            initialize=enth_E_param_dict,
            units=pyunits.J / pyunits.kg,
            doc="Specific enthalpy parameter E",
        )

        # vapor pressure parameters, 0-150 C
        # Sparrow 2003, Eq 6: Pvap = A+BT+CT2+DT3+ET4, where A = ao + a1X + a2X2 + a3X3 + a4X4
        vap_pressure_A1_param_dict = {
            "0": 0.9083e-3,
            "1": -0.569e-3,
            "2": 0.1945e-3,
            "3": -3.736e-3,
            "4": 2.82e-3,
        }
        vap_pressure_B1_param_dict = {
            "0": -0.0669e-3,
            "1": 0.0582e-3,
            "2": -0.1668e-3,
            "3": 0.676e-3,
            "4": -2.091e-3,
        }
        vap_pressure_C1_param_dict = {
            "0": 7.541e-6,
            "1": -5.143e-6,
            "2": 6.482e-6,
            "3": -52.62e-6,
            "4": 115.7e-6,
        }
        vap_pressure_D1_param_dict = {
            "0": -0.0922e-6,
            "1": 0.0649e-6,
            "2": -0.1313e-6,
            "3": 0.8024e-6,
            "4": -1.986e-6,
        }
        vap_pressure_E1_param_dict = {
            "0": 1.237e-9,
            "1": -0.753e-9,
            "2": 0.1448e-9,
            "3": -6.964e-9,
            "4": 14.61e-9,
        }

        self.vap_pressure_A1_param = Var(
            vap_pressure_A1_param_dict.keys(),
            domain=Reals,
            initialize=vap_pressure_A1_param_dict,
            units=pyunits.Pa,
            doc="Vapor pressure parameters A1",
        )
        self.vap_pressure_B1_param = Var(
            vap_pressure_B1_param_dict.keys(),
            domain=Reals,
            initialize=vap_pressure_B1_param_dict,
            units=pyunits.Pa,
            doc="Vapor pressure parameters B1",
        )
        self.vap_pressure_C1_param = Var(
            vap_pressure_C1_param_dict.keys(),
            domain=Reals,
            initialize=vap_pressure_C1_param_dict,
            units=pyunits.Pa,
            doc="Vapor pressure parameters C1",
        )
        self.vap_pressure_D1_param = Var(
            vap_pressure_D1_param_dict.keys(),
            domain=Reals,
            initialize=vap_pressure_D1_param_dict,
            units=pyunits.Pa,
            doc="Vapor pressure parameters D1",
        )
        self.vap_pressure_E1_param = Var(
            vap_pressure_E1_param_dict.keys(),
            domain=Reals,
            initialize=vap_pressure_E1_param_dict,
            units=pyunits.Pa,
            doc="Vapor pressure parameters E1",
        )

        # TODO: could add entropy if needed from Sparrow 2003, Eq. 9

        # thermal conductivity, 0-155 C
        # Regressed from Zaytsev & Aseev (1992):
        # th cond = A+BT+CT2+DT3, where A = ao + a1X + a2X2 + a3X3
        therm_cond_A_param_dict = {
            "0": 0.5424026,
            "1": 0.01283929,
            "2": -0.587953,
            "3": 1.090895,
        }
        therm_cond_B_param_dict = {
            "0": 0.002909031,
            "1": -0.001817648,
            "2": 0.007804725,
            "3": -0.01199839,
        }
        therm_cond_C_param_dict = {
            "0": -0.00002129933,
            "1": 0.0000275758,
            "2": -0.0001439831,
            "3": 0.000237931,
        }
        therm_cond_D_param_dict = {
            "0": 0.00000005486099,
            "1": -0.0000001044598,
            "2": 0.0000005747034,
            "3": -0.0000009645982,
        }

        self.therm_cond_param_A = Var(
            therm_cond_A_param_dict.keys(),
            domain=Reals,
            initialize=therm_cond_A_param_dict,
            units=pyunits.W / (pyunits.m * pyunits.K),
            doc="Thermal conductivity parameter A",
        )
        self.therm_cond_param_B = Var(
            therm_cond_B_param_dict.keys(),
            domain=Reals,
            initialize=therm_cond_B_param_dict,
            units=pyunits.W / (pyunits.m * pyunits.K),
            doc="Thermal conductivity parameter B",
        )
        self.therm_cond_param_C = Var(
            therm_cond_C_param_dict.keys(),
            domain=Reals,
            initialize=therm_cond_C_param_dict,
            units=pyunits.W / (pyunits.m * pyunits.K),
            doc="Thermal conductivity parameter C",
        )
        self.therm_cond_param_D = Var(
            therm_cond_D_param_dict.keys(),
            domain=Reals,
            initialize=therm_cond_D_param_dict,
            units=pyunits.W / (pyunits.m * pyunits.K),
            doc="Thermal conductivity parameter D",
        )

        # viscosity, 0-200 C,
        # Regressed from Zaytsev & Aseev (1992):
        # vis = A+BT+CT2+DT3+ET4, where A = ao + a1X + a2X2 + a3X3 + a4X4
        visc_A_param_dict = {
            "0": 1.740036,
            "1": 0.2437716,
            "2": 33.97162,
            "3": -140.441,
            "4": 381.4262,
        }
        visc_B_param_dict = {
            "0": -0.04347415,
            "1": 0.01794253,
            "2": -0.9344071,
            "3": 3.607937,
            "4": -9.966248,
        }
        visc_C_param_dict = {
            "0": 0.0005076655,
            "1": -0.0001613445,
            "2": 0.009094354,
            "3": -0.03051292,
            "4": 0.09635153,
        }
        visc_D_param_dict = {
            "0": -0.00000271643,
            "1": 0.0000004705967,
            "2": -0.00003883917,
            "3": 0.0001014729,
            "4": -0.0004023515,
        }
        visc_E_param_dict = {
            "0": 0.000000005339076,
            "1": -0.0000000002519368,
            "2": 0.00000006098774,
            "3": -0.000000102327,
            "4": 0.0000006067517,
        }

        self.visc_param_A = Var(
            visc_A_param_dict.keys(),
            domain=Reals,
            initialize=visc_A_param_dict,
            units=pyunits.Pa * pyunits.s,
            doc="Dynamic viscosity parameter A",
        )
        self.visc_param_B = Var(
            visc_B_param_dict.keys(),
            domain=Reals,
            initialize=visc_B_param_dict,
            units=pyunits.Pa * pyunits.s,
            doc="Dynamic viscosity parameter B",
        )
        self.visc_param_C = Var(
            visc_C_param_dict.keys(),
            domain=Reals,
            initialize=visc_C_param_dict,
            units=pyunits.Pa * pyunits.s,
            doc="Dynamic viscosity parameter C",
        )
        self.visc_param_D = Var(
            visc_D_param_dict.keys(),
            domain=Reals,
            initialize=visc_D_param_dict,
            units=pyunits.Pa * pyunits.s,
            doc="Dynamic viscosity parameter D",
        )
        self.visc_param_E = Var(
            visc_E_param_dict.keys(),
            domain=Reals,
            initialize=visc_E_param_dict,
            units=pyunits.Pa * pyunits.s,
            doc="Dynamic viscosity parameter E",
        )

        # heat capacity, 0-200 C
        # Regressed from Zaytsev & Aseev (1992):
        # cp = A+BT+CT2+DT3, where A = ao + a1X + a2X2 + a3X3
        cp_A_param_dict = {
            "0": 4174.74838,
            "1": -5533.00792,
            "2": 9564.358017,
            "3": -11918.208084,
        }
        cp_B_param_dict = {
            "0": -0.115581,
            "1": 17.432724,
            "2": -157.366583,
            "3": 392.883082,
        }
        cp_C_param_dict = {"0": 0.001698, "1": -0.137938, "2": 1.541653, "3": -3.799609}
        cp_D_param_dict = {"0": 0.000034, "1": 0.000376, "2": -0.004463, "3": 0.010904}
        self.cp_param_A = Var(
            cp_A_param_dict.keys(),
            domain=Reals,
            initialize=cp_A_param_dict,
            units=pyunits.J / (pyunits.kg * pyunits.K),
            doc="Heat Capacity parameter A",
        )
        self.cp_param_B = Var(
            cp_B_param_dict.keys(),
            domain=Reals,
            initialize=cp_B_param_dict,
            units=pyunits.J / (pyunits.kg * pyunits.K),
            doc="Heat Capacity parameter B",
        )
        self.cp_param_C = Var(
            cp_C_param_dict.keys(),
            domain=Reals,
            initialize=cp_C_param_dict,
            units=pyunits.J / (pyunits.kg * pyunits.K),
            doc="Heat Capacity parameter C",
        )
        self.cp_param_D = Var(
            cp_D_param_dict.keys(),
            domain=Reals,
            initialize=cp_D_param_dict,
            units=pyunits.J / (pyunits.kg * pyunits.K),
            doc="Heat Capacity parameter D",
        )

        # diffusivity parameters, 0-60 C
        # Regressed from Zaytsev & Aseev (1992):
        # diffus = A+BT+CT2+DT3, where A = ao + a1X + a2X2 + a3X3
        diffus_aq_A_param_dict = {
            "0": 0.4131329195225919,
            "1": -7.742251927618488,
            "2": 103.85062756812115,
            "3": -237.04211192526648,
        }
        diffus_aq_B_param_dict = {
            "0": 0.04632401864393698,
            "1": 0.8936760261966561,
            "2": -9.939718876008943,
            "3": 22.522378774378765,
        }
        diffus_aq_C_param_dict = {
            "0": -7.116300356324601e-06,
            "1": -0.03354460908477397,
            "2": 0.33333186888927646,
            "3": -0.7541027909550477,
        }
        diffus_aq_D_param_dict = {
            "0": -1.1147902601038595e-08,
            "1": 0.00027118535808978717,
            "2": -0.0026528149775375542,
            "3": 0.00603727069050608,
        }
        self.diffus_aq_param_A = Var(
            diffus_aq_A_param_dict.keys(),
            domain=Reals,
            initialize=diffus_aq_A_param_dict,
            units=pyunits.m**2 / pyunits.s,
            doc="Diffusivity (solution) parameter A",
        )
        self.diffus_aq_param_B = Var(
            diffus_aq_B_param_dict.keys(),
            domain=Reals,
            initialize=diffus_aq_B_param_dict,
            units=pyunits.m**2 / pyunits.s,
            doc="Diffusivity (solution) parameter B",
        )
        self.diffus_aq_param_C = Var(
            diffus_aq_C_param_dict.keys(),
            domain=Reals,
            initialize=diffus_aq_C_param_dict,
            units=pyunits.m**2 / pyunits.s,
            doc="Diffusivity (solution) parameter C",
        )
        self.diffus_aq_param_D = Var(
            diffus_aq_D_param_dict.keys(),
            domain=Reals,
            initialize=diffus_aq_D_param_dict,
            units=pyunits.m**2 / pyunits.s,
            doc="Diffusivity (solution) parameter D",
        )

        # osmotic coefficient parameters, 0-300 C
        # Regressed from Pitzer et. al. (1984):
        # osm coeff = A+BT+CT2+DT3, where A = ao + a1X + a2X2 + a3X3
        osm_coeff_A_param_dict = {
            "0": 0.9399062962108361,
            "1": -0.0189882837141575,
            "2": 0.019107595372900864,
            "3": -0.0011117486936487648,
        }
        osm_coeff_B_param_dict = {
            "0": -0.0007061546095086372,
            "1": 0.0008291452376952367,
            "2": -6.51004044513984e-05,
            "3": -8.040912167732404e-06,
        }
        osm_coeff_C_param_dict = {
            "0": 2.5295381526926904e-06,
            "1": -4.995407097616723e-06,
            "2": 3.626080914442103e-07,
            "3": 2.643648394495091e-08,
        }
        osm_coeff_D_param_dict = {
            "0": -4.6659633952665815e-09,
            "1": 3.766553344020095e-09,
            "2": 1.761632037484241e-10,
            "3": -7.493701294810002e-11,
        }
        self.osm_coeff_param_A = Var(
            osm_coeff_A_param_dict.keys(),
            domain=Reals,
            initialize=osm_coeff_A_param_dict,
            units=pyunits.dimensionless,
            doc="Osmotic coefficient parameter A",
        )
        self.osm_coeff_param_B = Var(
            osm_coeff_B_param_dict.keys(),
            domain=Reals,
            initialize=osm_coeff_B_param_dict,
            units=pyunits.dimensionless,
            doc="Osmotic coefficient parameter B",
        )
        self.osm_coeff_param_C = Var(
            osm_coeff_C_param_dict.keys(),
            domain=Reals,
            initialize=osm_coeff_C_param_dict,
            units=pyunits.dimensionless,
            doc="Osmotic coefficient parameter C",
        )
        self.osm_coeff_param_D = Var(
            osm_coeff_D_param_dict.keys(),
            domain=Reals,
            initialize=osm_coeff_D_param_dict,
            units=pyunits.dimensionless,
            doc="Osmotic coefficient parameter D",
        )

        # water density parameters from: water_prop_pack for liq water density
        dens_units = pyunits.kg / pyunits.m**3
        t_inv_units = pyunits.K**-1

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

        # TODO: add vapor phase and according properties (if deemed necessary)
        # including: specific enthalpy, diffusivity(water-air), density, specific heat, heat of vaporization

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
        self.set_default_scaling("enth_mass_phase", 1e-5, index="Liq")
        self.set_default_scaling("cp_mass_phase", 1e-4, index="Liq")
        self.set_default_scaling("pressure_sat", 1e-5)
        self.set_default_scaling("therm_cond_phase", 1e0, index="Liq")
        self.set_default_scaling("solubility_comp", 1e0)

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
                "cp_mass_phase": {"method": "_cp_mass_phase"},
                "pressure_sat": {"method": "_pressure_sat"},
            }
        )
        obj.define_custom_properties(
            {
                "osm_coeff": {"method": "_osm_coeff"},
                "enth_flow": {"method": "_enth_flow"},
                "solubility_comp": {"method": "_solubility_comp"},
                "therm_cond_phase": {"method": "_therm_cond_phase"},
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


class _NaClStateBlock(StateBlock):
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
        Method to release state variables fixed during initialization.

        Keyword Arguments:
            flags : dict containing information of which state variables
                    were fixed during initialization, and should now be
                    unfixed. This dict is returned by initialize if
                    hold_state=True.
            outlvl : sets output level of logging
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
            bounds=(273.15, 1000),
            domain=NonNegativeReals,
            units=pyunits.degK,
            doc="State temperature",
        )

        self.pressure = Var(
            initialize=101325,
            bounds=(1e3, 5e7),
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
            bounds=(1, 1e6),
            units=pyunits.kg * pyunits.m**-3,
            doc="Mass density",
        )
        # Sparrow 2003, Eq. 7, 0-300 C
        def rule_dens_mass_phase(b, p):
            t = (b.temperature - 273.15 * pyunits.K) / pyunits.K
            param_vec = [
                b.params.dens_mass_param_A,
                b.params.dens_mass_param_B,
                b.params.dens_mass_param_C,
                b.params.dens_mass_param_D,
                b.params.dens_mass_param_E,
            ]
            iter_param = {"A": 0, "B": 0, "C": 0, "D": 0, "E": 0}
            k = 0
            for key in iter_param:
                iter_param[key] = (
                    param_vec[k]["0"]
                    + param_vec[k]["1"] * b.mass_frac_phase_comp[p, "NaCl"]
                    + param_vec[k]["2"] * b.mass_frac_phase_comp[p, "NaCl"] ** 2
                    + param_vec[k]["3"] * b.mass_frac_phase_comp[p, "NaCl"] ** 3
                    + param_vec[k]["4"] * b.mass_frac_phase_comp[p, "NaCl"] ** 4
                )
                k += 1
            return (
                b.dens_mass_phase[p]
                == iter_param["A"]
                + iter_param["B"] * t
                + iter_param["C"] * t**2
                + iter_param["D"] * t**3
                + iter_param["E"] * t**4
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
            self.params.phase_list, ["NaCl"], rule=rule_molality_phase_comp
        )

    def _solubility_comp(self):
        self.solubility_comp = Var(
            ["NaCl"],
            initialize=0.5,
            bounds=(0.0, 1.01),
            units=pyunits.dimensionless,
            doc="solubility_comp",
        )
        # Sparrow 2003, Eq 5,  0-450 C
        def rule_solubility_comp(b, j):
            t = (b.temperature - 273.15 * pyunits.K) / pyunits.K
            return (
                self.solubility_comp[j]
                == b.params.solubility_comp_param["0"]
                + b.params.solubility_comp_param["1"] * t
                + b.params.solubility_comp_param["2"] * t**2
            )

        self.eq_solubility_comp = Constraint(["NaCl"], rule=rule_solubility_comp)

    def _pressure_sat(self):
        self.pressure_sat = Var(
            initialize=1e3,
            bounds=(1, 1e8),
            units=pyunits.Pa,
            doc="Vapor Pressure",
        )
        # Sparrow 2003, Eq. 6, 0-150 C
        def rule_pressure_sat(b):
            t = (b.temperature - 273.15 * pyunits.K) / pyunits.K
            x = b.mass_frac_phase_comp["Liq", "NaCl"]
            param_vec = [
                b.params.vap_pressure_A1_param,
                b.params.vap_pressure_B1_param,
                b.params.vap_pressure_C1_param,
                b.params.vap_pressure_D1_param,
                b.params.vap_pressure_E1_param,
            ]
            iter_param = {"A": 0, "B": 0, "C": 0, "D": 0, "E": 0}
            k = 0
            for key in iter_param:
                iter_param[key] = (
                    param_vec[k]["0"]
                    + param_vec[k]["1"] * x
                    + param_vec[k]["2"] * x**2
                    + param_vec[k]["3"] * x**3
                    + param_vec[k]["4"] * x**4
                )
                k += 1
            return (
                b.pressure_sat
                == (
                    iter_param["A"]
                    + iter_param["B"] * t
                    + iter_param["C"] * t**2
                    + iter_param["D"] * t**3
                    + iter_param["E"] * t**4
                )
                * 1e6  # convert MPa to Pa
            )

        self.eq_pressure_sat = Constraint(rule=rule_pressure_sat)

    def _visc_d_phase(self):
        self.visc_d_phase = Var(
            self.params.phase_list,
            initialize=1e-3,
            bounds=(0.0, 1),
            units=pyunits.Pa * pyunits.s,
            doc="Viscosity",
        )
        # Regressed from Zaytsev & Aseev (1992), 0-200 C
        def rule_visc_d_phase(b, p):
            t = (b.temperature - 273.15 * pyunits.K) / pyunits.K
            param_vec = [
                b.params.visc_param_A,
                b.params.visc_param_B,
                b.params.visc_param_C,
                b.params.visc_param_D,
                b.params.visc_param_E,
            ]
            iter_param = {"A": 0, "B": 0, "C": 0, "D": 0, "E": 0}
            k = 0
            for key in iter_param:
                iter_param[key] = (
                    param_vec[k]["0"]
                    + param_vec[k]["1"] * b.mass_frac_phase_comp[p, "NaCl"]
                    + param_vec[k]["2"] * b.mass_frac_phase_comp[p, "NaCl"] ** 2
                    + param_vec[k]["3"] * b.mass_frac_phase_comp[p, "NaCl"] ** 3
                    + param_vec[k]["4"] * b.mass_frac_phase_comp[p, "NaCl"] ** 4
                )
                k += 1
            return b.visc_d_phase[p] == (
                iter_param["A"]
                + iter_param["B"] * t
                + iter_param["C"] * t**2
                + iter_param["D"] * t**3
                + iter_param["E"] * t**4
            ) * 10 ** (-3)

        self.eq_visc_d_phase = Constraint(
            self.params.phase_list, rule=rule_visc_d_phase
        )

    def _therm_cond_phase(self):
        self.therm_cond_phase = Var(
            self.params.phase_list,
            initialize=0.6,
            bounds=(0.0, 2),
            units=pyunits.W / (pyunits.m * pyunits.K),
            doc="Thermal Conductivity",
        )
        # Regressed from Zaytsev & Aseev (1992), 0-155 C
        def rule_therm_cond_phase(b, p):
            t = (b.temperature - 273.15 * pyunits.K) / pyunits.K
            param_vec = [
                b.params.therm_cond_param_A,
                b.params.therm_cond_param_B,
                b.params.therm_cond_param_C,
                b.params.therm_cond_param_D,
            ]
            iter_param = {"A": 0, "B": 0, "C": 0, "D": 0}
            k = 0
            for key in iter_param:
                iter_param[key] = (
                    param_vec[k]["0"]
                    + param_vec[k]["1"] * b.mass_frac_phase_comp[p, "NaCl"]
                    + param_vec[k]["2"] * b.mass_frac_phase_comp[p, "NaCl"] ** 2
                    + param_vec[k]["3"] * b.mass_frac_phase_comp[p, "NaCl"] ** 3
                )
                k += 1
            return (
                b.therm_cond_phase[p]
                == iter_param["A"]
                + iter_param["B"] * t
                + iter_param["C"] * t**2
                + iter_param["D"] * t**3
            )

        self.eq_therm_cond_phase = Constraint(
            self.params.phase_list, rule=rule_therm_cond_phase
        )

    def _cp_mass_phase(self):
        self.cp_mass_phase = Var(
            self.params.phase_list,
            initialize=4e3,
            bounds=(0.0, 1e6),
            units=pyunits.J * pyunits.kg**-1 * pyunits.K**-1,
            doc="specific heat",
        )
        # Regressed from Zaytsev & Aseev (1992), 0-200 C
        def rule_cp_mass_phase(b, p):
            t = (b.temperature - 273.15 * pyunits.K) / pyunits.K
            param_vec = [
                b.params.cp_param_A,
                b.params.cp_param_B,
                b.params.cp_param_C,
                b.params.cp_param_D,
            ]
            iter_param = {"A": 0, "B": 0, "C": 0, "D": 0}
            k = 0
            for key in iter_param:
                iter_param[key] = (
                    param_vec[k]["0"]
                    + param_vec[k]["1"] * b.mass_frac_phase_comp[p, "NaCl"]
                    + param_vec[k]["2"] * b.mass_frac_phase_comp[p, "NaCl"] ** 2
                    + param_vec[k]["3"] * b.mass_frac_phase_comp[p, "NaCl"] ** 3
                )
                k += 1
            return (
                b.cp_mass_phase[p]
                == iter_param["A"]
                + iter_param["B"] * t
                + iter_param["C"] * t**2
                + iter_param["D"] * t**3
            )

        self.eq_cp_mass_phase = Constraint(
            self.params.phase_list, rule=rule_cp_mass_phase
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
        # Regressed from Zaytsev & Aseev (1992), 0-60 C
        def rule_diffus_phase_comp(b, p, j):
            t = (b.temperature - 273.15 * pyunits.K) / pyunits.K
            param_vec = [
                b.params.diffus_aq_param_A,
                b.params.diffus_aq_param_B,
                b.params.diffus_aq_param_C,
                b.params.diffus_aq_param_D,
            ]
            iter_param = {"A": 0, "B": 0, "C": 0, "D": 0}
            k = 0
            for key in iter_param:
                iter_param[key] = (
                    param_vec[k]["0"]
                    + param_vec[k]["1"] * b.mass_frac_phase_comp[p, "NaCl"]
                    + param_vec[k]["2"] * b.mass_frac_phase_comp[p, "NaCl"] ** 2
                    + param_vec[k]["3"] * b.mass_frac_phase_comp[p, "NaCl"] ** 3
                )
                k += 1
            return b.diffus_phase_comp[p, j] == (
                (
                    iter_param["A"]
                    + iter_param["B"] * t
                    + iter_param["C"] * t**2
                    + iter_param["D"] * t**3
                )
                * 1e-9
            )

        self.eq_diffus_phase_comp = Constraint(
            self.params.phase_list, ["NaCl"], rule=rule_diffus_phase_comp
        )

    def _osm_coeff(self):
        self.osm_coeff = Var(
            initialize=1,
            bounds=(0.0, 10),
            units=pyunits.dimensionless,
            doc="Osmotic coefficient",
        )
        # Regressed from Pitzer et. al. (1984), 0-300 C
        def rule_osm_coeff(b):
            t = (b.temperature - 273.15 * pyunits.K) / pyunits.K
            m = (b.molality_phase_comp["Liq", "NaCl"]) / (pyunits.mole / pyunits.kg)
            param_vec = [
                b.params.osm_coeff_param_A,
                b.params.osm_coeff_param_B,
                b.params.osm_coeff_param_C,
                b.params.osm_coeff_param_D,
            ]
            iter_param = {"A": 0, "B": 0, "C": 0, "D": 0}
            k = 0
            for key in iter_param:
                iter_param[key] = (
                    param_vec[k]["0"]
                    + param_vec[k]["1"] * m
                    + param_vec[k]["2"] * m**2
                    + param_vec[k]["3"] * m**3
                )
                k += 1
            return (
                b.osm_coeff
                == iter_param["A"]
                + iter_param["B"] * t
                + iter_param["C"] * t**2
                + iter_param["D"] * t**3
            )

        self.eq_osm_coeff = Constraint(rule=rule_osm_coeff)

    def _pressure_osm_phase(self):
        self.pressure_osm_phase = Var(
            self.params.phase_list,
            initialize=1e6,
            bounds=(1, 1e8),
            units=pyunits.Pa,
            doc="Osmotic pressure",
        )

        def rule_pressure_osm_phase(b, p):
            i = 2  # number of ionic species
            t = b.temperature - 273.15 * pyunits.K
            dens_mass = (
                b.params.dens_mass_param_A1
                + b.params.dens_mass_param_A2 * t
                + b.params.dens_mass_param_A3 * t**2
                + b.params.dens_mass_param_A4 * t**3
                + b.params.dens_mass_param_A5 * t**4
            )
            return (
                b.pressure_osm_phase[p]
                == i
                * b.osm_coeff
                * b.molality_phase_comp[p, "NaCl"]
                * dens_mass
                * Constants.gas_constant
                * b.temperature
            )

        self.eq_pressure_osm_phase = Constraint(
            self.params.phase_list, rule=rule_pressure_osm_phase
        )

    def _enth_mass_phase(self):
        self.enth_mass_phase = Var(
            self.params.phase_list,
            initialize=1e2,
            bounds=(0, 1e8),
            units=pyunits.J * pyunits.kg**-1,
            doc="Specific enthalpy",
        )
        # Sparrow 2003, Eq 8, 0-300 C
        def rule_enth_mass_phase(b, p):
            t = (b.temperature - 273.15 * pyunits.K) / pyunits.K
            param_vec = [
                b.params.enth_param_A,
                b.params.enth_param_B,
                b.params.enth_param_C,
                b.params.enth_param_D,
                b.params.enth_param_E,
            ]
            iter_param = {"A": 0, "B": 0, "C": 0, "D": 0, "E": 0}
            k = 0
            for key in iter_param:
                iter_param[key] = (
                    param_vec[k]["0"]
                    + param_vec[k]["1"] * b.mass_frac_phase_comp[p, "NaCl"]
                    + param_vec[k]["2"] * b.mass_frac_phase_comp[p, "NaCl"] ** 2
                    + param_vec[k]["3"] * b.mass_frac_phase_comp[p, "NaCl"] ** 3
                    + param_vec[k]["4"] * b.mass_frac_phase_comp[p, "NaCl"] ** 4
                )
                k += 1
            return b.enth_mass_phase[p] == (
                (
                    iter_param["A"]
                    + iter_param["B"] * t
                    + iter_param["C"] * t**2
                    + iter_param["D"] * t**3
                    + iter_param["E"] * t**4
                )
                * 1000
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
