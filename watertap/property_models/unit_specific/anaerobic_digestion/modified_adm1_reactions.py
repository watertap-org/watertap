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
Modified ADM1 reaction package.

Reference:

X. Flores-Alsina, K. Solon, C.K. Mbamba, S. Tait, K.V. Gernaey, U. Jeppsson, D.J. Batstone,
Modelling phosphorus (P), sulfur (S) and iron (Fe) interactions fordynamic simulations of anaerobic digestion processes,
Water Research. 95 (2016) 370-382. https://www.sciencedirect.com/science/article/pii/S0043135416301397

"""

# Import Pyomo libraries
import pyomo.environ as pyo
from pyomo.environ import Suffix

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    MaterialFlowBasis,
    ReactionParameterBlock,
    ReactionBlockDataBase,
    ReactionBlockBase,
)
from idaes.core.util.constants import Constants
from idaes.core.util.misc import add_object_reference
from idaes.core.util.exceptions import BurntToast
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale
from idaes.core.util.math import smooth_max
from idaes.core.scaling import CustomScalerBase, ConstraintScalingScheme

# Some more information about this module
__author__ = "Chenyu Wang, Marcus Holly, Xinhong Liu"

# Set up logger
_log = idaeslog.getLogger(__name__)

mw_n = 14 * pyo.units.kg / pyo.units.kmol
mw_c = 12 * pyo.units.kg / pyo.units.kmol
mw_p = 31 * pyo.units.kg / pyo.units.kmol


@declare_process_block_class("ModifiedADM1ReactionParameterBlock")
class ModifiedADM1ReactionParameterData(ReactionParameterBlock):
    """
    Property Parameter Block Class
    """

    def build(self):
        """
        Callable method for Block construction.
        """
        super().build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        self._reaction_block_class = ModifiedADM1ReactionBlock

        # Reaction Index
        # Reaction names based on standard numbering in ADM1
        # R1:  Hydrolysis of carbohydrates
        # R2:  Hydrolysis of proteins
        # R3:  Hydrolysis of lipids
        # R4:  Uptake of sugars
        # R5:  Uptake of amino acids
        # R6:  Uptake of long chain fatty acids (LCFAs)
        # R7:  Uptake of valerate
        # R8:  Uptake of butyrate
        # R9:  Uptake of propionate
        # R10: Uptake of acetate
        # R11: Uptake of hydrogen
        # R12: Decay of X_su
        # R13: Decay of X_aa
        # R14: Decay of X_fa
        # R15: Decay of X_c4
        # R16: Decay of X_pro
        # R17: Decay of X_ac
        # R18: Decay of X_h2
        # R19: Storage of S_va in X_PHA
        # R20: Storage of S_bu in X_PHA
        # R21: Storage of S_pro in X_PHA
        # R22: Storage of S_ac in X_PHA
        # R23: Lysis of X_PAO
        # R24: Lysis of X_PP
        # R25: Lysis of X_PHA

        self.rate_reaction_idx = pyo.Set(
            initialize=[
                "R1",
                "R2",
                "R3",
                "R4",
                "R5",
                "R6",
                "R7",
                "R8",
                "R9",
                "R10",
                "R11",
                "R12",
                "R13",
                "R14",
                "R15",
                "R16",
                "R17",
                "R18",
                "R19",
                "R20",
                "R21",
                "R22",
                "R23",
                "R24",
                "R25",
            ]
        )

        # Carbon content
        Ci_dict = {
            "S_su": 0.031250000,
            "S_aa": 0.030741549,
            "S_fa": 0.021404110,
            "S_va": 0.024038462,
            "S_bu": 0.02500000,
            "S_pro": 0.026785714,
            "S_ac": 0.0312500000,
            "S_ch4": 0.015625000,
            "S_I": 0.030148202,
            "X_ch": 0.031250000,
            "X_pr": 0.030741549,
            "X_li": 0.021925591,
            "X_su": 0.030509792,
            "X_aa": 0.030509792,
            "X_fa": 0.030509792,
            "X_c4": 0.030509792,
            "X_pro": 0.030509792,
            "X_ac": 0.030509792,
            "X_h2": 0.030509792,
            "X_I": 0.030148202,
            "X_PHA": 0.02500000,
            "X_PAO": 0.030509792,
        }

        self.Ci = pyo.Var(
            Ci_dict.keys(),
            initialize=Ci_dict,
            units=pyo.units.kmol / pyo.units.kg,
            domain=pyo.PositiveReals,
            doc="Carbon content of component [kmole C/kg COD]",
        )

        # Nitrogen content
        Ni_dict = {
            "S_aa": 0.0079034,
            "S_I": 0.0042876,
            "X_pr": 0.0079034,
            "X_su": 0.0061532,
            "X_aa": 0.0061532,
            "X_fa": 0.0061532,
            "X_c4": 0.0061532,
            "X_pro": 0.0061532,
            "X_ac": 0.0061532,
            "X_h2": 0.0061532,
            "X_I": 0.0042876,
            "X_PAO": 0.0061532,
        }

        self.Ni = pyo.Var(
            Ni_dict.keys(),
            initialize=Ni_dict,
            units=pyo.units.kmol / pyo.units.kg,
            domain=pyo.PositiveReals,
            doc="Nitrogen content of component [kmole N/kg COD]",
        )

        # Phosphorus content
        Pi_dict = {
            "S_I": 0.0002093322,
            "X_li": 0.0003440808,
            "X_su": 0.0006947201,
            "X_aa": 0.0006947201,
            "X_fa": 0.0006947201,
            "X_c4": 0.0006947201,
            "X_pro": 0.0006947201,
            "X_ac": 0.0006947201,
            "X_h2": 0.0006947201,
            "X_I": 0.0002093322,
            "X_PP": 1,
            "X_PAO": 0.0006947201,
        }

        self.Pi = pyo.Var(
            Pi_dict.keys(),
            initialize=Pi_dict,
            units=pyo.units.kmol / pyo.units.kg,
            domain=pyo.PositiveReals,
            doc="Phosphorus content of component [kmole P/kg COD]",
        )

        # Common parameters for both translator blocks
        self.f_sI_xc = pyo.Param(
            initialize=0,
            units=pyo.units.dimensionless,
            mutable=True,
            doc="Soluble inerts from composites",
        )
        self.f_xI_xc = pyo.Param(
            initialize=0.1,
            units=pyo.units.dimensionless,
            mutable=True,
            doc="Particulate inerts from composites",
        )
        self.f_ch_xc = pyo.Param(
            initialize=0.275,
            units=pyo.units.dimensionless,
            mutable=True,
            doc="Carbohydrates from composites",
        )
        self.f_pr_xc = pyo.Param(
            initialize=0.275,
            units=pyo.units.dimensionless,
            mutable=True,
            doc="Proteins from composites",
        )
        self.f_li_xc = pyo.Param(
            initialize=0.35,
            units=pyo.units.dimensionless,
            mutable=True,
            doc="Lipids from composites",
        )

        self.P_ch = pyo.Param(
            initialize=0,
            units=pyo.units.kmol / pyo.units.kg,
            mutable=True,
            doc="P content of X_ch",
        )

        # TODO: Consider inheriting these parameters from ADM1 such that there is less repeated code

        # Stoichiometric Parameters (Table 1.1 and 2.1 in Flores-Alsina et al., 2016)
        self.Z_h2s = pyo.Param(
            within=pyo.NonNegativeReals,
            mutable=True,
            default=0,
            doc="Reference component mass concentration of hydrogen sulfide",
            units=pyo.units.kg / pyo.units.m**3,
        )
        self.f_xi_xb = pyo.Var(
            initialize=0.1,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Fraction of inert particulate organics from biomass",
        )
        self.f_ch_xb = pyo.Var(
            initialize=0.275,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Fraction of carbohydrates from biomass",
        )
        self.f_li_xb = pyo.Var(
            initialize=0.350,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Fraction of lipids from biomass",
        )
        self.f_pr_xb = pyo.Var(
            initialize=0.275,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Fraction of proteins from biomass",
        )
        self.f_si_xb = pyo.Var(
            initialize=0,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Fraction of soluble inerts from biomass",
        )
        self.f_fa_li = pyo.Var(
            initialize=0.95,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Fatty acids from lipids",
        )
        self.f_h2_su = pyo.Var(
            initialize=0.1906,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Hydrogen from sugars",
        )
        self.f_bu_su = pyo.Var(
            initialize=0.1328,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Butyrate from sugars",
        )
        self.f_pro_su = pyo.Var(
            initialize=0.2691,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Propionate from sugars",
        )
        self.f_ac_su = pyo.Var(
            initialize=0.4076,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Acetate from sugars",
        )
        self.f_h2_aa = pyo.Var(
            initialize=0.06,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Hydrogen from amino acids",
        )

        self.f_va_aa = pyo.Var(
            initialize=0.23,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Valerate from amino acids",
        )
        self.f_bu_aa = pyo.Var(
            initialize=0.26,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Butyrate from amino acids",
        )
        self.f_pro_aa = pyo.Var(
            initialize=0.05,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Propionate from amino acids",
        )
        self.f_ac_aa = pyo.Var(
            initialize=0.40,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Acetate from amino acids",
        )
        self.Y_su = pyo.Var(
            initialize=0.10,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Yield of biomass on sugar substrate [kg COD X/ kg COD S]",
        )
        self.Y_aa = pyo.Var(
            initialize=0.08,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Yield of biomass on amino acid substrate [kg COD X/ kg COD S]",
        )
        self.Y_fa = pyo.Var(
            initialize=0.06,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Yield of biomass on fatty acid substrate [kg COD X/ kg COD S]",
        )
        self.Y_c4 = pyo.Var(
            initialize=0.06,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Yield of biomass on valerate and butyrate substrate [kg COD X/ kg COD S]",
        )
        self.Y_pro = pyo.Var(
            initialize=0.04,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Yield of biomass on propionate substrate [kg COD X/ kg COD S]",
        )
        self.Y_ac = pyo.Var(
            initialize=0.05,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Yield of biomass on acetate substrate [kg COD X/ kg COD S]",
        )
        self.Y_h2 = pyo.Var(
            initialize=0.06,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Yield of hydrogen per biomass [kg COD S/ kg COD X]",
        )
        # Biochemical Parameters
        self.k_hyd_ch = pyo.Var(
            initialize=10,
            units=pyo.units.day**-1,
            domain=pyo.PositiveReals,
            doc="First-order kinetic parameter for hydrolysis of carbohydrates",
        )
        self.k_hyd_pr = pyo.Var(
            initialize=10,
            units=pyo.units.day**-1,
            domain=pyo.PositiveReals,
            doc="First-order kinetic parameter for hydrolysis of proteins",
        )
        self.k_hyd_li = pyo.Var(
            initialize=10,
            units=pyo.units.day**-1,
            domain=pyo.PositiveReals,
            doc="First-order kinetic parameter for hydrolysis of lipids",
        )
        self.K_S_IN = pyo.Var(
            initialize=1e-4,
            units=pyo.units.kmol * pyo.units.m**-3,
            domain=pyo.PositiveReals,
            doc="Inhibition parameter for inorganic nitrogen",
        )
        self.k_m_su = pyo.Var(
            initialize=30,
            units=pyo.units.day**-1,
            domain=pyo.PositiveReals,
            doc="Monod maximum specific uptake rate of sugars",
        )
        self.K_S_su = pyo.Var(
            initialize=0.5,
            units=pyo.units.kg * pyo.units.m**-3,
            domain=pyo.PositiveReals,
            doc="Half saturation value for uptake of sugars",
        )
        self.pH_UL_aa = pyo.Var(
            initialize=5.5,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Upper limit of pH for uptake rate of amino acids",
        )
        self.pH_LL_aa = pyo.Var(
            initialize=4,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Lower limit of pH for uptake rate of amino acids",
        )
        self.k_m_aa = pyo.Var(
            initialize=50,
            units=pyo.units.day**-1,
            domain=pyo.PositiveReals,
            doc="Monod maximum specific uptake rate of amino acids",
        )

        self.K_S_aa = pyo.Var(
            initialize=0.3,
            units=pyo.units.kg * pyo.units.m**-3,
            domain=pyo.PositiveReals,
            doc="Half saturation value for uptake of amino acids",
        )
        self.k_m_fa = pyo.Var(
            initialize=6,
            units=pyo.units.day**-1,
            domain=pyo.PositiveReals,
            doc="Monod maximum specific uptake rate of fatty acids",
        )
        self.K_S_fa = pyo.Var(
            initialize=0.4,
            units=pyo.units.kg * pyo.units.m**-3,
            domain=pyo.PositiveReals,
            doc="Half saturation value for uptake of fatty acids",
        )
        self.K_I_h2_fa = pyo.Var(
            initialize=5e-6,
            units=pyo.units.kg * pyo.units.m**-3,
            domain=pyo.PositiveReals,
            doc="Inhibition parameter for hydrogen during uptake of fatty acids",
        )
        self.k_m_c4 = pyo.Var(
            initialize=20,
            units=pyo.units.day**-1,
            domain=pyo.PositiveReals,
            doc="Monod maximum specific uptake rate of valerate and butyrate",
        )
        self.K_S_c4 = pyo.Var(
            initialize=0.2,
            units=pyo.units.kg * pyo.units.m**-3,
            domain=pyo.PositiveReals,
            doc="Half saturation value for uptake of valerate and butyrate",
        )
        self.K_I_h2_c4 = pyo.Var(
            initialize=1e-5,
            units=pyo.units.kg * pyo.units.m**-3,
            domain=pyo.PositiveReals,
            doc="Inhibition parameter for hydrogen during uptake of valerate and butyrate",
        )
        self.k_m_pro = pyo.Var(
            initialize=13,
            units=pyo.units.day**-1,
            domain=pyo.PositiveReals,
            doc="Monod maximum specific uptake rate of propionate",
        )
        self.K_S_pro = pyo.Var(
            initialize=0.1,
            units=pyo.units.kg * pyo.units.m**-3,
            domain=pyo.PositiveReals,
            doc="Half saturation value for uptake of propionate",
        )
        self.K_I_h2_pro = pyo.Var(
            initialize=3.5e-6,
            units=pyo.units.kg * pyo.units.m**-3,
            domain=pyo.PositiveReals,
            doc="Inhibition parameter for hydrogen during uptake of propionate",
        )
        self.k_m_ac = pyo.Var(
            initialize=8,
            units=pyo.units.day**-1,
            domain=pyo.PositiveReals,
            doc="Monod maximum specific uptake rate of acetate",
        )
        self.K_S_ac = pyo.Var(
            initialize=0.15,
            units=pyo.units.kg * pyo.units.m**-3,
            domain=pyo.PositiveReals,
            doc="Half saturation value for uptake of acetate",
        )
        self.K_I_nh3 = pyo.Var(
            initialize=0.0018,
            units=pyo.units.kmol * pyo.units.m**-3,
            domain=pyo.PositiveReals,
            doc="Inhibition parameter for ammonia during uptake of acetate",
        )
        self.pH_UL_ac = pyo.Var(
            initialize=7,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Upper limit of pH for uptake rate of acetate",
        )
        self.pH_LL_ac = pyo.Var(
            initialize=6,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Lower limit of pH for uptake rate of acetate",
        )
        self.k_m_h2 = pyo.Var(
            initialize=35,
            units=pyo.units.day**-1,
            domain=pyo.PositiveReals,
            doc="Monod maximum specific uptake rate of hydrogen",
        )
        self.K_S_h2 = pyo.Var(
            initialize=7e-6,
            units=pyo.units.kg * pyo.units.m**-3,
            domain=pyo.PositiveReals,
            doc="Half saturation value for uptake of hydrogen",
        )
        self.pH_UL_h2 = pyo.Var(
            initialize=6,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Upper limit of pH for uptake rate of hydrogen",
        )
        self.pH_LL_h2 = pyo.Var(
            initialize=5,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Lower limit of pH for uptake rate of hydrogen",
        )
        self.k_dec_X_su = pyo.Var(
            initialize=0.02,
            units=pyo.units.day**-1,
            domain=pyo.PositiveReals,
            doc="First-order decay rate for X_su",
        )
        self.k_dec_X_aa = pyo.Var(
            initialize=0.02,
            units=pyo.units.day**-1,
            domain=pyo.PositiveReals,
            doc="First-order decay rate for X_aa",
        )
        self.k_dec_X_fa = pyo.Var(
            initialize=0.02,
            units=pyo.units.day**-1,
            domain=pyo.PositiveReals,
            doc="First-order decay rate for X_fa",
        )
        self.k_dec_X_c4 = pyo.Var(
            initialize=0.02,
            units=pyo.units.day**-1,
            domain=pyo.PositiveReals,
            doc="First-order decay rate for X_c4",
        )
        self.k_dec_X_pro = pyo.Var(
            initialize=0.02,
            units=pyo.units.day**-1,
            domain=pyo.PositiveReals,
            doc="First-order decay rate for X_pro",
        )
        self.k_dec_X_ac = pyo.Var(
            initialize=0.02,
            units=pyo.units.day**-1,
            domain=pyo.PositiveReals,
            doc="First-order decay rate for X_ac",
        )
        self.k_dec_X_h2 = pyo.Var(
            initialize=0.02,
            units=pyo.units.day**-1,
            domain=pyo.PositiveReals,
            doc="First-order decay rate for X_h2",
        )
        self.K_a_va = pyo.Var(
            initialize=10 ** (-4.86),
            units=pyo.units.kmol / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Valerate acid-base equilibrium constant",
        )
        self.K_a_bu = pyo.Var(
            initialize=10 ** (-4.82),
            units=pyo.units.kmol / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Butyrate acid-base equilibrium constant",
        )
        self.K_a_pro = pyo.Var(
            initialize=10 ** (-4.88),
            units=pyo.units.kmol / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Propionate acid-base equilibrium constant",
        )
        self.K_a_ac = pyo.Var(
            initialize=10 ** (-4.76),
            units=pyo.units.kmol / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Acetate acid-base equilibrium constant",
        )
        self.K_I_h2s_ac = pyo.Var(
            initialize=460e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="50% inhibitory concentration of H2S on acetogens",
        )
        self.K_I_h2s_c4 = pyo.Var(
            initialize=481e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="50% inhibitory concentration of H2S on c4 degraders",
        )
        self.K_I_h2s_h2 = pyo.Var(
            initialize=400e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="50% inhibitory concentration of H2S on hydrogenotrophic methanogens",
        )
        self.K_I_h2s_pro = pyo.Var(
            initialize=481e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="50% inhibitory concentration of propionate degraders",
        )
        self.K_S_IP = pyo.Var(
            initialize=2e-5,
            units=pyo.units.kmol / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="P limitation for inorganic phosphorus",
        )
        self.b_PAO = pyo.Var(
            initialize=0.2,
            units=pyo.units.day**-1,
            domain=pyo.PositiveReals,
            doc="Lysis rate of phosphorus accumulating organisms",
        )
        self.b_PHA = pyo.Var(
            initialize=0.2,
            units=pyo.units.day**-1,
            domain=pyo.PositiveReals,
            doc="Lysis rate of polyhydroxyalkanoates",
        )
        self.b_PP = pyo.Var(
            initialize=0.2,
            units=pyo.units.day**-1,
            domain=pyo.PositiveReals,
            doc="Lysis rate of polyphosphates",
        )
        self.f_ac_PHA = pyo.Var(
            initialize=0.4,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Yield of acetate on polyhydroxyalkanoates",
        )
        self.f_bu_PHA = pyo.Var(
            initialize=0.1,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Yield of butyrate on polyhydroxyalkanoates",
        )
        self.f_pro_PHA = pyo.Var(
            initialize=0.4,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Yield of propionate on polyhydroxyalkanoates",
        )
        self.f_va_PHA = pyo.Var(
            initialize=0.1,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Yield of valerate on polyhydroxyalkanoates",
        )
        self.K_A = pyo.Var(
            initialize=4e-3,
            units=pyo.units.kg * pyo.units.m**-3,
            domain=pyo.PositiveReals,
            doc="Saturation coefficient for acetate",
        )
        self.K_PP = pyo.Var(
            initialize=0.32e-3,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Saturation coefficient for polyphosphate",
        )
        self.q_PHA = pyo.Var(
            initialize=3.0,
            units=pyo.units.day**-1,
            domain=pyo.PositiveReals,
            doc="Rate constant for storage of polyhydroxyalkanoates",
        )
        self.Y_PO4 = pyo.Var(
            initialize=12.903e-3,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Yield of biomass on phosphate (kmol P/kg COD)",
        )
        self.K_XPP = pyo.Var(
            initialize=1 / 3,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Potassium coefficient for polyphosphates",
        )
        self.Mg_XPP = pyo.Var(
            initialize=1 / 3,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Magnesium coefficient for polyphosphates",
        )
        self.temperature_ref = pyo.Param(
            within=pyo.PositiveReals,
            mutable=True,
            default=298.15,
            doc="Reference temperature",
            units=pyo.units.K,
        )

        # Reaction Stoichiometry
        # This is the stoichiometric part of the Peterson matrix in dict form.
        # See Table 1.1 and 2.1 in Flores-Alsina et al., 2016.

        # Exclude non-zero stoichiometric coefficients for S_IC initially since they depend on other stoichiometric coefficients.
        self.rate_reaction_stoichiometry = {
            # R1: Hydrolysis of carbohydrates
            ("R1", "Liq", "H2O"): 0,
            ("R1", "Liq", "S_su"): 1,
            ("R1", "Liq", "S_aa"): 0,
            ("R1", "Liq", "S_fa"): 0,
            ("R1", "Liq", "S_va"): 0,
            ("R1", "Liq", "S_bu"): 0,
            ("R1", "Liq", "S_pro"): 0,
            ("R1", "Liq", "S_ac"): 0,
            ("R1", "Liq", "S_h2"): 0,
            ("R1", "Liq", "S_ch4"): 0,
            ("R1", "Liq", "S_IC"): -(self.Ci["S_su"] - self.Ci["X_ch"]) * mw_c,
            ("R1", "Liq", "S_IN"): 0,
            ("R1", "Liq", "S_IP"): 0,
            ("R1", "Liq", "S_I"): 0,
            ("R1", "Liq", "X_ch"): -1,
            ("R1", "Liq", "X_pr"): 0,
            ("R1", "Liq", "X_li"): 0,
            ("R1", "Liq", "X_su"): 0,
            ("R1", "Liq", "X_aa"): 0,
            ("R1", "Liq", "X_fa"): 0,
            ("R1", "Liq", "X_c4"): 0,
            ("R1", "Liq", "X_pro"): 0,
            ("R1", "Liq", "X_ac"): 0,
            ("R1", "Liq", "X_h2"): 0,
            ("R1", "Liq", "X_I"): 0,
            ("R1", "Liq", "X_PHA"): 0,
            ("R1", "Liq", "X_PP"): 0,
            ("R1", "Liq", "X_PAO"): 0,
            ("R1", "Liq", "S_K"): 0,
            ("R1", "Liq", "S_Mg"): 0,
            # R2:  Hydrolysis of proteins
            ("R2", "Liq", "H2O"): 0,
            ("R2", "Liq", "S_su"): 0,
            ("R2", "Liq", "S_aa"): 1,
            ("R2", "Liq", "S_fa"): 0,
            ("R2", "Liq", "S_va"): 0,
            ("R2", "Liq", "S_bu"): 0,
            ("R2", "Liq", "S_pro"): 0,
            ("R2", "Liq", "S_ac"): 0,
            ("R2", "Liq", "S_h2"): 0,
            ("R2", "Liq", "S_ch4"): 0,
            ("R2", "Liq", "S_IC"): -(self.Ci["S_aa"] - self.Ci["X_pr"]) * mw_c,
            ("R2", "Liq", "S_IN"): -(self.Ni["S_aa"] - self.Ni["X_pr"]) * mw_n,
            ("R2", "Liq", "S_IP"): 0,
            ("R2", "Liq", "S_I"): 0,
            ("R2", "Liq", "X_ch"): 0,
            ("R2", "Liq", "X_pr"): -1,
            ("R2", "Liq", "X_li"): 0,
            ("R2", "Liq", "X_su"): 0,
            ("R2", "Liq", "X_aa"): 0,
            ("R2", "Liq", "X_fa"): 0,
            ("R2", "Liq", "X_c4"): 0,
            ("R2", "Liq", "X_pro"): 0,
            ("R2", "Liq", "X_ac"): 0,
            ("R2", "Liq", "X_h2"): 0,
            ("R2", "Liq", "X_I"): 0,
            ("R2", "Liq", "X_PHA"): 0,
            ("R2", "Liq", "X_PP"): 0,
            ("R2", "Liq", "X_PAO"): 0,
            ("R2", "Liq", "S_K"): 0,
            ("R2", "Liq", "S_Mg"): 0,
            # R3:  Hydrolysis of lipids
            ("R3", "Liq", "H2O"): 0,
            ("R3", "Liq", "S_su"): 1 - self.f_fa_li,
            ("R3", "Liq", "S_aa"): 0,
            ("R3", "Liq", "S_fa"): self.f_fa_li,
            ("R3", "Liq", "S_va"): 0,
            ("R3", "Liq", "S_bu"): 0,
            ("R3", "Liq", "S_pro"): 0,
            ("R3", "Liq", "S_ac"): 0,
            ("R3", "Liq", "S_h2"): 0,
            ("R3", "Liq", "S_ch4"): 0,
            ("R3", "Liq", "S_IC"): (
                self.Ci["X_li"]
                - (1 - self.f_fa_li) * self.Ci["S_su"]
                - self.f_fa_li * self.Ci["S_fa"]
            )
            * mw_c,
            ("R3", "Liq", "S_IN"): 0,
            ("R3", "Liq", "S_IP"): self.Pi["X_li"] * mw_p,
            ("R3", "Liq", "S_I"): 0,
            ("R3", "Liq", "X_ch"): 0,
            ("R3", "Liq", "X_pr"): 0,
            ("R3", "Liq", "X_li"): -1,
            ("R3", "Liq", "X_su"): 0,
            ("R3", "Liq", "X_aa"): 0,
            ("R3", "Liq", "X_fa"): 0,
            ("R3", "Liq", "X_c4"): 0,
            ("R3", "Liq", "X_pro"): 0,
            ("R3", "Liq", "X_ac"): 0,
            ("R3", "Liq", "X_h2"): 0,
            ("R3", "Liq", "X_I"): 0,
            ("R3", "Liq", "X_PHA"): 0,
            ("R3", "Liq", "X_PP"): 0,
            ("R3", "Liq", "X_PAO"): 0,
            ("R3", "Liq", "S_K"): 0,
            ("R3", "Liq", "S_Mg"): 0,
            # R4:  Uptake of sugars
            ("R4", "Liq", "H2O"): 0,
            ("R4", "Liq", "S_su"): -1,
            ("R4", "Liq", "S_aa"): 0,
            ("R4", "Liq", "S_fa"): 0,
            ("R4", "Liq", "S_va"): 0,
            ("R4", "Liq", "S_bu"): (1 - self.Y_su) * self.f_bu_su,
            ("R4", "Liq", "S_pro"): (1 - self.Y_su) * self.f_pro_su,
            ("R4", "Liq", "S_ac"): (1 - self.Y_su) * self.f_ac_su,
            ("R4", "Liq", "S_h2"): (1 - self.Y_su) * self.f_h2_su,
            ("R4", "Liq", "S_ch4"): 0,
            ("R4", "Liq", "S_IC"): (
                self.Ci["S_su"]
                - (1 - self.Y_su)
                * (
                    self.f_bu_su * self.Ci["S_bu"]
                    + self.f_pro_su * self.Ci["S_pro"]
                    + self.f_ac_su * self.Ci["S_ac"]
                )
                - self.Y_su * self.Ci["X_su"]
            )
            * mw_c,
            ("R4", "Liq", "S_IN"): (-self.Y_su * self.Ni["X_su"]) * mw_n,
            ("R4", "Liq", "S_IP"): (-self.Y_su * self.Pi["X_su"]) * mw_p,
            ("R4", "Liq", "S_I"): 0,
            ("R4", "Liq", "X_ch"): 0,
            ("R4", "Liq", "X_pr"): 0,
            ("R4", "Liq", "X_li"): 0,
            ("R4", "Liq", "X_su"): self.Y_su,
            ("R4", "Liq", "X_aa"): 0,
            ("R4", "Liq", "X_fa"): 0,
            ("R4", "Liq", "X_c4"): 0,
            ("R4", "Liq", "X_pro"): 0,
            ("R4", "Liq", "X_ac"): 0,
            ("R4", "Liq", "X_h2"): 0,
            ("R4", "Liq", "X_I"): 0,
            ("R4", "Liq", "X_PHA"): 0,
            ("R4", "Liq", "X_PP"): 0,
            ("R4", "Liq", "X_PAO"): 0,
            ("R4", "Liq", "S_K"): 0,
            ("R4", "Liq", "S_Mg"): 0,
            # R5:  Uptake of amino acids
            ("R5", "Liq", "H2O"): 0,
            ("R5", "Liq", "S_su"): 0,
            ("R5", "Liq", "S_aa"): -1,
            ("R5", "Liq", "S_fa"): 0,
            ("R5", "Liq", "S_va"): (1 - self.Y_aa) * self.f_va_aa,
            ("R5", "Liq", "S_bu"): (1 - self.Y_aa) * self.f_bu_aa,
            ("R5", "Liq", "S_pro"): (1 - self.Y_aa) * self.f_pro_aa,
            ("R5", "Liq", "S_ac"): (1 - self.Y_aa) * self.f_ac_aa,
            ("R5", "Liq", "S_h2"): (1 - self.Y_aa) * self.f_h2_aa,
            ("R5", "Liq", "S_ch4"): 0,
            ("R5", "Liq", "S_IC"): (
                self.Ci["S_aa"]
                - (1 - self.Y_aa)
                * (
                    self.f_va_aa * self.Ci["S_va"]
                    + self.f_bu_aa * self.Ci["S_bu"]
                    + self.f_pro_aa * self.Ci["S_pro"]
                    + self.f_ac_aa * self.Ci["S_ac"]
                )
                - self.Y_aa * self.Ci["X_aa"]
            )
            * mw_c,
            ("R5", "Liq", "S_IN"): -(-self.Ni["S_aa"] + self.Y_aa * self.Ni["X_aa"])
            * mw_n,
            ("R5", "Liq", "S_IP"): -(self.Y_aa * self.Pi["X_aa"]) * mw_p,
            ("R5", "Liq", "S_I"): 0,
            ("R5", "Liq", "X_ch"): 0,
            ("R5", "Liq", "X_pr"): 0,
            ("R5", "Liq", "X_li"): 0,
            ("R5", "Liq", "X_su"): 0,
            ("R5", "Liq", "X_aa"): self.Y_aa,
            ("R5", "Liq", "X_fa"): 0,
            ("R5", "Liq", "X_c4"): 0,
            ("R5", "Liq", "X_pro"): 0,
            ("R5", "Liq", "X_ac"): 0,
            ("R5", "Liq", "X_h2"): 0,
            ("R5", "Liq", "X_I"): 0,
            ("R5", "Liq", "X_PHA"): 0,
            ("R5", "Liq", "X_PP"): 0,
            ("R5", "Liq", "X_PAO"): 0,
            ("R5", "Liq", "S_K"): 0,
            ("R5", "Liq", "S_Mg"): 0,
            # R6:  Uptake of long chain fatty acids (LCFAs)
            ("R6", "Liq", "H2O"): 0,
            ("R6", "Liq", "S_su"): 0,
            ("R6", "Liq", "S_aa"): 0,
            ("R6", "Liq", "S_fa"): -1,
            ("R6", "Liq", "S_va"): 0,
            ("R6", "Liq", "S_bu"): 0,
            ("R6", "Liq", "S_pro"): 0,
            ("R6", "Liq", "S_ac"): (1 - self.Y_fa) * 0.7,
            ("R6", "Liq", "S_h2"): (1 - self.Y_fa) * 0.3,
            ("R6", "Liq", "S_ch4"): 0,
            ("R6", "Liq", "S_IC"): (
                self.Ci["S_fa"]
                - (1 - self.Y_fa) * 0.7 * self.Ci["S_ac"]
                - self.Y_fa * self.Ci["X_fa"]
            )
            * mw_c,
            ("R6", "Liq", "S_IN"): (-self.Y_fa * self.Ni["X_fa"]) * mw_n,
            ("R6", "Liq", "S_IP"): (-self.Y_fa * self.Pi["X_fa"]) * mw_p,
            ("R6", "Liq", "S_I"): 0,
            ("R6", "Liq", "X_ch"): 0,
            ("R6", "Liq", "X_pr"): 0,
            ("R6", "Liq", "X_li"): 0,
            ("R6", "Liq", "X_su"): 0,
            ("R6", "Liq", "X_aa"): 0,
            ("R6", "Liq", "X_fa"): self.Y_fa,
            ("R6", "Liq", "X_c4"): 0,
            ("R6", "Liq", "X_pro"): 0,
            ("R6", "Liq", "X_ac"): 0,
            ("R6", "Liq", "X_h2"): 0,
            ("R6", "Liq", "X_I"): 0,
            ("R6", "Liq", "X_PHA"): 0,
            ("R6", "Liq", "X_PP"): 0,
            ("R6", "Liq", "X_PAO"): 0,
            ("R6", "Liq", "S_K"): 0,
            ("R6", "Liq", "S_Mg"): 0,
            # R7:  Uptake of valerate
            ("R7", "Liq", "H2O"): 0,
            ("R7", "Liq", "S_su"): 0,
            ("R7", "Liq", "S_aa"): 0,
            ("R7", "Liq", "S_fa"): 0,
            ("R7", "Liq", "S_va"): -1,
            ("R7", "Liq", "S_bu"): 0,
            ("R7", "Liq", "S_pro"): (1 - self.Y_c4) * 0.54,
            ("R7", "Liq", "S_ac"): (1 - self.Y_c4) * 0.31,
            ("R7", "Liq", "S_h2"): (1 - self.Y_c4) * 0.15,
            ("R7", "Liq", "S_ch4"): 0,
            ("R7", "Liq", "S_IC"): (
                self.Ci["S_va"]
                - (1 - self.Y_c4) * 0.54 * self.Ci["S_pro"]
                - (1 - self.Y_c4) * 0.31 * self.Ci["S_ac"]
                - self.Y_c4 * self.Ci["X_c4"]
            )
            * mw_c,
            ("R7", "Liq", "S_IN"): (-self.Y_c4 * self.Ni["X_c4"]) * mw_n,
            ("R7", "Liq", "S_IP"): (-self.Y_c4 * self.Pi["X_c4"]) * mw_p,
            ("R7", "Liq", "S_I"): 0,
            ("R7", "Liq", "X_ch"): 0,
            ("R7", "Liq", "X_pr"): 0,
            ("R7", "Liq", "X_li"): 0,
            ("R7", "Liq", "X_su"): 0,
            ("R7", "Liq", "X_aa"): 0,
            ("R7", "Liq", "X_fa"): 0,
            ("R7", "Liq", "X_c4"): self.Y_c4,
            ("R7", "Liq", "X_pro"): 0,
            ("R7", "Liq", "X_ac"): 0,
            ("R7", "Liq", "X_h2"): 0,
            ("R7", "Liq", "X_I"): 0,
            ("R7", "Liq", "X_PHA"): 0,
            ("R7", "Liq", "X_PP"): 0,
            ("R7", "Liq", "X_PAO"): 0,
            ("R7", "Liq", "S_K"): 0,
            ("R7", "Liq", "S_Mg"): 0,
            # R8:  Uptake of butyrate
            ("R8", "Liq", "H2O"): 0,
            ("R8", "Liq", "S_su"): 0,
            ("R8", "Liq", "S_aa"): 0,
            ("R8", "Liq", "S_fa"): 0,
            ("R8", "Liq", "S_va"): 0,
            ("R8", "Liq", "S_bu"): -1,
            ("R8", "Liq", "S_pro"): 0,
            ("R8", "Liq", "S_ac"): (1 - self.Y_c4) * 0.8,
            ("R8", "Liq", "S_h2"): (1 - self.Y_c4) * 0.2,
            ("R8", "Liq", "S_ch4"): 0,
            ("R8", "Liq", "S_IC"): (
                self.Ci["S_bu"]
                - (1 - self.Y_c4) * 0.8 * self.Ci["S_ac"]
                - self.Y_c4 * self.Ci["X_c4"]
            )
            * mw_c,
            ("R8", "Liq", "S_IN"): (-self.Y_c4 * self.Ni["X_c4"]) * mw_n,
            ("R8", "Liq", "S_IP"): (-self.Y_c4 * self.Pi["X_c4"]) * mw_p,
            ("R8", "Liq", "S_I"): 0,
            ("R8", "Liq", "X_ch"): 0,
            ("R8", "Liq", "X_pr"): 0,
            ("R8", "Liq", "X_li"): 0,
            ("R8", "Liq", "X_su"): 0,
            ("R8", "Liq", "X_aa"): 0,
            ("R8", "Liq", "X_fa"): 0,
            ("R8", "Liq", "X_c4"): self.Y_c4,
            ("R8", "Liq", "X_pro"): 0,
            ("R8", "Liq", "X_ac"): 0,
            ("R8", "Liq", "X_h2"): 0,
            ("R8", "Liq", "X_I"): 0,
            ("R8", "Liq", "X_PHA"): 0,
            ("R8", "Liq", "X_PP"): 0,
            ("R8", "Liq", "X_PAO"): 0,
            ("R8", "Liq", "S_K"): 0,
            ("R8", "Liq", "S_Mg"): 0,
            # R9: Uptake of propionate
            ("R9", "Liq", "H2O"): 0,
            ("R9", "Liq", "S_su"): 0,
            ("R9", "Liq", "S_aa"): 0,
            ("R9", "Liq", "S_fa"): 0,
            ("R9", "Liq", "S_va"): 0,
            ("R9", "Liq", "S_bu"): 0,
            ("R9", "Liq", "S_pro"): -1,
            ("R9", "Liq", "S_ac"): (1 - self.Y_pro) * 0.57,
            ("R9", "Liq", "S_h2"): (1 - self.Y_pro) * 0.43,
            ("R9", "Liq", "S_ch4"): 0,
            ("R9", "Liq", "S_IC"): (
                self.Ci["S_pro"]
                - (1 - self.Y_pro) * 0.57 * self.Ci["S_ac"]
                - self.Y_pro * self.Ci["X_pro"]
            )
            * mw_c,
            ("R9", "Liq", "S_IN"): (-self.Y_pro * self.Ni["X_pro"]) * mw_n,
            ("R9", "Liq", "S_IP"): (-self.Y_pro * self.Pi["X_pro"]) * mw_p,
            ("R9", "Liq", "S_I"): 0,
            ("R9", "Liq", "X_ch"): 0,
            ("R9", "Liq", "X_pr"): 0,
            ("R9", "Liq", "X_li"): 0,
            ("R9", "Liq", "X_su"): 0,
            ("R9", "Liq", "X_aa"): 0,
            ("R9", "Liq", "X_fa"): 0,
            ("R9", "Liq", "X_c4"): 0,
            ("R9", "Liq", "X_pro"): self.Y_pro,
            ("R9", "Liq", "X_ac"): 0,
            ("R9", "Liq", "X_h2"): 0,
            ("R9", "Liq", "X_I"): 0,
            ("R9", "Liq", "X_PHA"): 0,
            ("R9", "Liq", "X_PP"): 0,
            ("R9", "Liq", "X_PAO"): 0,
            ("R9", "Liq", "S_K"): 0,
            ("R9", "Liq", "S_Mg"): 0,
            # R10: Uptake of acetate
            ("R10", "Liq", "H2O"): 0,
            ("R10", "Liq", "S_su"): 0,
            ("R10", "Liq", "S_aa"): 0,
            ("R10", "Liq", "S_fa"): 0,
            ("R10", "Liq", "S_va"): 0,
            ("R10", "Liq", "S_bu"): 0,
            ("R10", "Liq", "S_pro"): 0,
            ("R10", "Liq", "S_ac"): -1,
            ("R10", "Liq", "S_h2"): 0,
            ("R10", "Liq", "S_ch4"): 1 - self.Y_ac,
            ("R10", "Liq", "S_IC"): (
                self.Ci["S_ac"]
                - (1 - self.Y_ac) * self.Ci["S_ch4"]
                - self.Y_ac * self.Ci["X_ac"]
            )
            * mw_c,
            ("R10", "Liq", "S_IN"): (-self.Y_ac * self.Ni["X_ac"]) * mw_n,
            ("R10", "Liq", "S_IP"): (-self.Y_ac * self.Pi["X_ac"]) * mw_p,
            ("R10", "Liq", "S_I"): 0,
            ("R10", "Liq", "X_ch"): 0,
            ("R10", "Liq", "X_pr"): 0,
            ("R10", "Liq", "X_li"): 0,
            ("R10", "Liq", "X_su"): 0,
            ("R10", "Liq", "X_aa"): 0,
            ("R10", "Liq", "X_fa"): 0,
            ("R10", "Liq", "X_c4"): 0,
            ("R10", "Liq", "X_pro"): 0,
            ("R10", "Liq", "X_ac"): self.Y_ac,
            ("R10", "Liq", "X_h2"): 0,
            ("R10", "Liq", "X_I"): 0,
            ("R10", "Liq", "X_PHA"): 0,
            ("R10", "Liq", "X_PP"): 0,
            ("R10", "Liq", "X_PAO"): 0,
            ("R10", "Liq", "S_K"): 0,
            ("R10", "Liq", "S_Mg"): 0,
            # R11: Uptake of hydrogen
            ("R11", "Liq", "H2O"): 0,
            ("R11", "Liq", "S_su"): 0,
            ("R11", "Liq", "S_aa"): 0,
            ("R11", "Liq", "S_fa"): 0,
            ("R11", "Liq", "S_va"): 0,
            ("R11", "Liq", "S_bu"): 0,
            ("R11", "Liq", "S_pro"): 0,
            ("R11", "Liq", "S_ac"): 0,
            ("R11", "Liq", "S_h2"): -1,
            ("R11", "Liq", "S_ch4"): 1 - self.Y_h2,
            ("R11", "Liq", "S_IC"): (
                -(1 - self.Y_h2) * self.Ci["S_ch4"] - self.Y_h2 * self.Ci["X_h2"]
            )
            * mw_c,
            ("R11", "Liq", "S_IN"): (-self.Y_h2 * self.Ni["X_h2"]) * mw_n,
            ("R11", "Liq", "S_IP"): (-self.Y_h2 * self.Pi["X_h2"]) * mw_p,
            ("R11", "Liq", "S_I"): 0,
            ("R11", "Liq", "X_ch"): 0,
            ("R11", "Liq", "X_pr"): 0,
            ("R11", "Liq", "X_li"): 0,
            ("R11", "Liq", "X_su"): 0,
            ("R11", "Liq", "X_aa"): 0,
            ("R11", "Liq", "X_fa"): 0,
            ("R11", "Liq", "X_c4"): 0,
            ("R11", "Liq", "X_pro"): 0,
            ("R11", "Liq", "X_ac"): 0,
            ("R11", "Liq", "X_h2"): self.Y_h2,
            ("R11", "Liq", "X_I"): 0,
            ("R11", "Liq", "X_PHA"): 0,
            ("R11", "Liq", "X_PP"): 0,
            ("R11", "Liq", "X_PAO"): 0,
            ("R11", "Liq", "S_K"): 0,
            ("R11", "Liq", "S_Mg"): 0,
            # R12: Decay of X_su
            ("R12", "Liq", "H2O"): 0,
            ("R12", "Liq", "S_su"): 0,
            ("R12", "Liq", "S_aa"): 0,
            ("R12", "Liq", "S_fa"): 0,
            ("R12", "Liq", "S_va"): 0,
            ("R12", "Liq", "S_bu"): 0,
            ("R12", "Liq", "S_pro"): 0,
            ("R12", "Liq", "S_ac"): 0,
            ("R12", "Liq", "S_h2"): 0,
            ("R12", "Liq", "S_ch4"): 0,
            ("R12", "Liq", "S_IC"): (
                self.Ci["X_su"]
                - self.f_ch_xb * self.Ci["X_ch"]
                - self.f_pr_xb * self.Ci["X_pr"]
                - self.f_li_xb * self.Ci["X_li"]
                - self.f_xi_xb * self.Ci["X_I"]
            )
            * mw_c,
            ("R12", "Liq", "S_IN"): (
                self.Ni["X_su"]
                - self.f_pr_xb * self.Ni["X_pr"]
                - self.f_xi_xb * self.Ni["X_I"]
            )
            * mw_n,
            ("R12", "Liq", "S_IP"): (
                self.Pi["X_su"]
                - self.f_li_xb * self.Pi["X_li"]
                - self.f_xi_xb * self.Pi["X_I"]
            )
            * mw_p,
            ("R12", "Liq", "S_I"): self.f_si_xb,
            ("R12", "Liq", "X_ch"): self.f_ch_xb,
            ("R12", "Liq", "X_pr"): self.f_pr_xb,
            ("R12", "Liq", "X_li"): self.f_li_xb,
            ("R12", "Liq", "X_su"): -1,
            ("R12", "Liq", "X_aa"): 0,
            ("R12", "Liq", "X_fa"): 0,
            ("R12", "Liq", "X_c4"): 0,
            ("R12", "Liq", "X_pro"): 0,
            ("R12", "Liq", "X_ac"): 0,
            ("R12", "Liq", "X_h2"): 0,
            ("R12", "Liq", "X_I"): self.f_xi_xb,
            ("R12", "Liq", "X_PHA"): 0,
            ("R12", "Liq", "X_PP"): 0,
            ("R12", "Liq", "X_PAO"): 0,
            ("R12", "Liq", "S_K"): 0,
            ("R12", "Liq", "S_Mg"): 0,
            # R13: Decay of X_aa
            ("R13", "Liq", "H2O"): 0,
            ("R13", "Liq", "S_su"): 0,
            ("R13", "Liq", "S_aa"): 0,
            ("R13", "Liq", "S_fa"): 0,
            ("R13", "Liq", "S_va"): 0,
            ("R13", "Liq", "S_bu"): 0,
            ("R13", "Liq", "S_pro"): 0,
            ("R13", "Liq", "S_ac"): 0,
            ("R13", "Liq", "S_h2"): 0,
            ("R13", "Liq", "S_ch4"): 0,
            ("R13", "Liq", "S_IC"): (
                self.Ci["X_aa"]
                - self.f_ch_xb * self.Ci["X_ch"]
                - self.f_pr_xb * self.Ci["X_pr"]
                - self.f_li_xb * self.Ci["X_li"]
                - self.f_xi_xb * self.Ci["X_I"]
            )
            * mw_c,
            ("R13", "Liq", "S_IN"): (
                self.Ni["X_aa"]
                - self.f_pr_xb * self.Ni["X_pr"]
                - self.f_xi_xb * self.Ni["X_I"]
            )
            * mw_n,
            ("R13", "Liq", "S_IP"): (
                self.Pi["X_aa"]
                - self.f_li_xb * self.Pi["X_li"]
                - self.f_xi_xb * self.Pi["X_I"]
            )
            * mw_p,
            ("R13", "Liq", "S_I"): self.f_si_xb,
            ("R13", "Liq", "X_ch"): self.f_ch_xb,
            ("R13", "Liq", "X_pr"): self.f_pr_xb,
            ("R13", "Liq", "X_li"): self.f_li_xb,
            ("R13", "Liq", "X_su"): 0,
            ("R13", "Liq", "X_aa"): -1,
            ("R13", "Liq", "X_fa"): 0,
            ("R13", "Liq", "X_c4"): 0,
            ("R13", "Liq", "X_pro"): 0,
            ("R13", "Liq", "X_ac"): 0,
            ("R13", "Liq", "X_h2"): 0,
            ("R13", "Liq", "X_I"): self.f_xi_xb,
            ("R13", "Liq", "X_PHA"): 0,
            ("R13", "Liq", "X_PP"): 0,
            ("R13", "Liq", "X_PAO"): 0,
            ("R13", "Liq", "S_K"): 0,
            ("R13", "Liq", "S_Mg"): 0,
            # R14: Decay of X_fa
            ("R14", "Liq", "H2O"): 0,
            ("R14", "Liq", "S_su"): 0,
            ("R14", "Liq", "S_aa"): 0,
            ("R14", "Liq", "S_fa"): 0,
            ("R14", "Liq", "S_va"): 0,
            ("R14", "Liq", "S_bu"): 0,
            ("R14", "Liq", "S_pro"): 0,
            ("R14", "Liq", "S_ac"): 0,
            ("R14", "Liq", "S_h2"): 0,
            ("R14", "Liq", "S_ch4"): 0,
            ("R14", "Liq", "S_IC"): (
                self.Ci["X_fa"]
                - self.f_ch_xb * self.Ci["X_ch"]
                - self.f_pr_xb * self.Ci["X_pr"]
                - self.f_li_xb * self.Ci["X_li"]
                - self.f_xi_xb * self.Ci["X_I"]
            )
            * mw_c,
            ("R14", "Liq", "S_IN"): (
                self.Ni["X_fa"]
                - self.f_pr_xb * self.Ni["X_pr"]
                - self.f_xi_xb * self.Ni["X_I"]
            )
            * mw_n,
            ("R14", "Liq", "S_IP"): (
                self.Pi["X_fa"]
                - self.f_li_xb * self.Pi["X_li"]
                - self.f_xi_xb * self.Pi["X_I"]
            )
            * mw_p,
            ("R14", "Liq", "S_I"): self.f_si_xb,
            ("R14", "Liq", "X_ch"): self.f_ch_xb,
            ("R14", "Liq", "X_pr"): self.f_pr_xb,
            ("R14", "Liq", "X_li"): self.f_li_xb,
            ("R14", "Liq", "X_su"): 0,
            ("R14", "Liq", "X_aa"): 0,
            ("R14", "Liq", "X_fa"): -1,
            ("R14", "Liq", "X_c4"): 0,
            ("R14", "Liq", "X_pro"): 0,
            ("R14", "Liq", "X_ac"): 0,
            ("R14", "Liq", "X_h2"): 0,
            ("R14", "Liq", "X_I"): self.f_xi_xb,
            ("R14", "Liq", "X_PHA"): 0,
            ("R14", "Liq", "X_PP"): 0,
            ("R14", "Liq", "X_PAO"): 0,
            ("R14", "Liq", "S_K"): 0,
            ("R14", "Liq", "S_Mg"): 0,
            # R15: Decay of X_c4
            ("R15", "Liq", "H2O"): 0,
            ("R15", "Liq", "S_su"): 0,
            ("R15", "Liq", "S_aa"): 0,
            ("R15", "Liq", "S_fa"): 0,
            ("R15", "Liq", "S_va"): 0,
            ("R15", "Liq", "S_bu"): 0,
            ("R15", "Liq", "S_pro"): 0,
            ("R15", "Liq", "S_ac"): 0,
            ("R15", "Liq", "S_h2"): 0,
            ("R15", "Liq", "S_ch4"): 0,
            ("R15", "Liq", "S_IC"): (
                self.Ci["X_c4"]
                - self.f_ch_xb * self.Ci["X_ch"]
                - self.f_pr_xb * self.Ci["X_pr"]
                - self.f_li_xb * self.Ci["X_li"]
                - self.f_xi_xb * self.Ci["X_I"]
            )
            * mw_c,
            ("R15", "Liq", "S_IN"): (
                self.Ni["X_c4"]
                - self.f_pr_xb * self.Ni["X_pr"]
                - self.f_xi_xb * self.Ni["X_I"]
            )
            * mw_n,
            ("R15", "Liq", "S_IP"): (
                self.Pi["X_c4"]
                - self.f_li_xb * self.Pi["X_li"]
                - self.f_xi_xb * self.Pi["X_I"]
            )
            * mw_p,
            ("R15", "Liq", "S_I"): self.f_si_xb,
            ("R15", "Liq", "X_ch"): self.f_ch_xb,
            ("R15", "Liq", "X_pr"): self.f_pr_xb,
            ("R15", "Liq", "X_li"): self.f_li_xb,
            ("R15", "Liq", "X_su"): 0,
            ("R15", "Liq", "X_aa"): 0,
            ("R15", "Liq", "X_fa"): 0,
            ("R15", "Liq", "X_c4"): -1,
            ("R15", "Liq", "X_pro"): 0,
            ("R15", "Liq", "X_ac"): 0,
            ("R15", "Liq", "X_h2"): 0,
            ("R15", "Liq", "X_I"): self.f_xi_xb,
            ("R15", "Liq", "X_PHA"): 0,
            ("R15", "Liq", "X_PP"): 0,
            ("R15", "Liq", "X_PAO"): 0,
            ("R15", "Liq", "S_K"): 0,
            ("R15", "Liq", "S_Mg"): 0,
            # R16: Decay of X_pro
            ("R16", "Liq", "H2O"): 0,
            ("R16", "Liq", "S_su"): 0,
            ("R16", "Liq", "S_aa"): 0,
            ("R16", "Liq", "S_fa"): 0,
            ("R16", "Liq", "S_va"): 0,
            ("R16", "Liq", "S_bu"): 0,
            ("R16", "Liq", "S_pro"): 0,
            ("R16", "Liq", "S_ac"): 0,
            ("R16", "Liq", "S_h2"): 0,
            ("R16", "Liq", "S_ch4"): 0,
            ("R16", "Liq", "S_IC"): (
                self.Ci["X_pro"]
                - self.f_ch_xb * self.Ci["X_ch"]
                - self.f_pr_xb * self.Ci["X_pr"]
                - self.f_li_xb * self.Ci["X_li"]
                - self.f_xi_xb * self.Ci["X_I"]
            )
            * mw_c,
            ("R16", "Liq", "S_IN"): (
                self.Ni["X_pro"]
                - self.f_pr_xb * self.Ni["X_pr"]
                - self.f_xi_xb * self.Ni["X_I"]
            )
            * mw_n,
            ("R16", "Liq", "S_IP"): (
                self.Pi["X_pro"]
                - self.f_li_xb * self.Pi["X_li"]
                - self.f_xi_xb * self.Pi["X_I"]
            )
            * mw_p,
            ("R16", "Liq", "S_I"): self.f_si_xb,
            ("R16", "Liq", "X_ch"): self.f_ch_xb,
            ("R16", "Liq", "X_pr"): self.f_pr_xb,
            ("R16", "Liq", "X_li"): self.f_li_xb,
            ("R16", "Liq", "X_su"): 0,
            ("R16", "Liq", "X_aa"): 0,
            ("R16", "Liq", "X_fa"): 0,
            ("R16", "Liq", "X_c4"): 0,
            ("R16", "Liq", "X_pro"): -1,
            ("R16", "Liq", "X_ac"): 0,
            ("R16", "Liq", "X_h2"): 0,
            ("R16", "Liq", "X_I"): self.f_xi_xb,
            ("R16", "Liq", "X_PHA"): 0,
            ("R16", "Liq", "X_PP"): 0,
            ("R16", "Liq", "X_PAO"): 0,
            ("R16", "Liq", "S_K"): 0,
            ("R16", "Liq", "S_Mg"): 0,
            # R17: Decay of X_ac
            ("R17", "Liq", "H2O"): 0,
            ("R17", "Liq", "S_su"): 0,
            ("R17", "Liq", "S_aa"): 0,
            ("R17", "Liq", "S_fa"): 0,
            ("R17", "Liq", "S_va"): 0,
            ("R17", "Liq", "S_bu"): 0,
            ("R17", "Liq", "S_pro"): 0,
            ("R17", "Liq", "S_ac"): 0,
            ("R17", "Liq", "S_h2"): 0,
            ("R17", "Liq", "S_ch4"): 0,
            ("R17", "Liq", "S_IC"): (
                self.Ci["X_ac"]
                - self.f_ch_xb * self.Ci["X_ch"]
                - self.f_pr_xb * self.Ci["X_pr"]
                - self.f_li_xb * self.Ci["X_li"]
                - self.f_xi_xb * self.Ci["X_I"]
            )
            * mw_c,
            ("R17", "Liq", "S_IN"): (
                self.Ni["X_ac"]
                - self.f_pr_xb * self.Ni["X_pr"]
                - self.f_xi_xb * self.Ni["X_I"]
            )
            * mw_n,
            ("R17", "Liq", "S_IP"): (
                self.Pi["X_ac"]
                - self.f_li_xb * self.Pi["X_li"]
                - self.f_xi_xb * self.Pi["X_I"]
            )
            * mw_p,
            ("R17", "Liq", "S_I"): self.f_si_xb,
            ("R17", "Liq", "X_ch"): self.f_ch_xb,
            ("R17", "Liq", "X_pr"): self.f_pr_xb,
            ("R17", "Liq", "X_li"): self.f_li_xb,
            ("R17", "Liq", "X_su"): 0,
            ("R17", "Liq", "X_aa"): 0,
            ("R17", "Liq", "X_fa"): 0,
            ("R17", "Liq", "X_c4"): 0,
            ("R17", "Liq", "X_pro"): 0,
            ("R17", "Liq", "X_ac"): -1,
            ("R17", "Liq", "X_h2"): 0,
            ("R17", "Liq", "X_I"): self.f_xi_xb,
            ("R17", "Liq", "X_PHA"): 0,
            ("R17", "Liq", "X_PP"): 0,
            ("R17", "Liq", "X_PAO"): 0,
            ("R17", "Liq", "S_K"): 0,
            ("R17", "Liq", "S_Mg"): 0,
            # R18: Decay of X_h2
            ("R18", "Liq", "H2O"): 0,
            ("R18", "Liq", "S_su"): 0,
            ("R18", "Liq", "S_aa"): 0,
            ("R18", "Liq", "S_fa"): 0,
            ("R18", "Liq", "S_va"): 0,
            ("R18", "Liq", "S_bu"): 0,
            ("R18", "Liq", "S_pro"): 0,
            ("R18", "Liq", "S_ac"): 0,
            ("R18", "Liq", "S_h2"): 0,
            ("R18", "Liq", "S_ch4"): 0,
            ("R18", "Liq", "S_IC"): (
                self.Ci["X_h2"]
                - self.f_ch_xb * self.Ci["X_ch"]
                - self.f_pr_xb * self.Ci["X_pr"]
                - self.f_li_xb * self.Ci["X_li"]
                - self.f_xi_xb * self.Ci["X_I"]
            )
            * mw_c,
            ("R18", "Liq", "S_IN"): (
                self.Ni["X_h2"]
                - self.f_pr_xb * self.Ni["X_pr"]
                - self.f_xi_xb * self.Ni["X_I"]
            )
            * mw_n,
            ("R18", "Liq", "S_IP"): (
                self.Pi["X_h2"]
                - self.f_li_xb * self.Pi["X_li"]
                - self.f_xi_xb * self.Pi["X_I"]
            )
            * mw_p,
            ("R18", "Liq", "S_I"): self.f_si_xb,
            ("R18", "Liq", "X_ch"): self.f_ch_xb,
            ("R18", "Liq", "X_pr"): self.f_pr_xb,
            ("R18", "Liq", "X_li"): self.f_li_xb,
            ("R18", "Liq", "X_su"): 0,
            ("R18", "Liq", "X_aa"): 0,
            ("R18", "Liq", "X_fa"): 0,
            ("R18", "Liq", "X_c4"): 0,
            ("R18", "Liq", "X_pro"): 0,
            ("R18", "Liq", "X_ac"): 0,
            ("R18", "Liq", "X_h2"): -1,
            ("R18", "Liq", "X_I"): self.f_xi_xb,
            ("R18", "Liq", "X_PHA"): 0,
            ("R18", "Liq", "X_PP"): 0,
            ("R18", "Liq", "X_PAO"): 0,
            ("R18", "Liq", "S_K"): 0,
            ("R18", "Liq", "S_Mg"): 0,
            # R19: Storage of S_va in X_PHA
            ("R19", "Liq", "H2O"): 0,
            ("R19", "Liq", "S_su"): 0,
            ("R19", "Liq", "S_aa"): 0,
            ("R19", "Liq", "S_fa"): 0,
            ("R19", "Liq", "S_va"): -1,
            ("R19", "Liq", "S_bu"): 0,
            ("R19", "Liq", "S_pro"): 0,
            ("R19", "Liq", "S_ac"): 0,
            ("R19", "Liq", "S_h2"): 0,
            ("R19", "Liq", "S_ch4"): 0,
            ("R19", "Liq", "S_IC"): (self.Ci["S_va"] - self.Ci["X_PHA"]) * mw_c,
            ("R19", "Liq", "S_IN"): 0,
            ("R19", "Liq", "S_IP"): (self.Y_PO4 * self.Pi["X_PP"]) * mw_p,
            ("R19", "Liq", "S_I"): 0,
            ("R19", "Liq", "X_ch"): 0,
            ("R19", "Liq", "X_pr"): 0,
            ("R19", "Liq", "X_li"): 0,
            ("R19", "Liq", "X_su"): 0,
            ("R19", "Liq", "X_aa"): 0,
            ("R19", "Liq", "X_fa"): 0,
            ("R19", "Liq", "X_c4"): 0,
            ("R19", "Liq", "X_pro"): 0,
            ("R19", "Liq", "X_ac"): 0,
            ("R19", "Liq", "X_h2"): 0,
            ("R19", "Liq", "X_I"): 0,
            ("R19", "Liq", "X_PHA"): 1,
            ("R19", "Liq", "X_PP"): -self.Y_PO4,
            ("R19", "Liq", "X_PAO"): 0,
            ("R19", "Liq", "S_K"): self.Y_PO4 * self.K_XPP,
            ("R19", "Liq", "S_Mg"): self.Y_PO4 * self.Mg_XPP,
            # R20: Storage of S_bu in X_PHA
            ("R20", "Liq", "H2O"): 0,
            ("R20", "Liq", "S_su"): 0,
            ("R20", "Liq", "S_aa"): 0,
            ("R20", "Liq", "S_fa"): 0,
            ("R20", "Liq", "S_va"): 0,
            ("R20", "Liq", "S_bu"): -1,
            ("R20", "Liq", "S_pro"): 0,
            ("R20", "Liq", "S_ac"): 0,
            ("R20", "Liq", "S_h2"): 0,
            ("R20", "Liq", "S_ch4"): 0,
            ("R20", "Liq", "S_IC"): (self.Ci["S_bu"] - self.Ci["X_PHA"]) * mw_c,
            ("R20", "Liq", "S_IN"): 0,
            ("R20", "Liq", "S_IP"): (self.Y_PO4 * self.Pi["X_PP"]) * mw_p,
            ("R20", "Liq", "S_I"): 0,
            ("R20", "Liq", "X_ch"): 0,
            ("R20", "Liq", "X_pr"): 0,
            ("R20", "Liq", "X_li"): 0,
            ("R20", "Liq", "X_su"): 0,
            ("R20", "Liq", "X_aa"): 0,
            ("R20", "Liq", "X_fa"): 0,
            ("R20", "Liq", "X_c4"): 0,
            ("R20", "Liq", "X_pro"): 0,
            ("R20", "Liq", "X_ac"): 0,
            ("R20", "Liq", "X_h2"): 0,
            ("R20", "Liq", "X_I"): 0,
            ("R20", "Liq", "X_PHA"): 1,
            ("R20", "Liq", "X_PP"): -self.Y_PO4,
            ("R20", "Liq", "X_PAO"): 0,
            ("R20", "Liq", "S_K"): self.Y_PO4 * self.K_XPP,
            ("R20", "Liq", "S_Mg"): self.Y_PO4 * self.Mg_XPP,
            # R21: Storage of S_pro in X_PHA
            ("R21", "Liq", "H2O"): 0,
            ("R21", "Liq", "S_su"): 0,
            ("R21", "Liq", "S_aa"): 0,
            ("R21", "Liq", "S_fa"): 0,
            ("R21", "Liq", "S_va"): 0,
            ("R21", "Liq", "S_bu"): 0,
            ("R21", "Liq", "S_pro"): -1,
            ("R21", "Liq", "S_ac"): 0,
            ("R21", "Liq", "S_h2"): 0,
            ("R21", "Liq", "S_ch4"): 0,
            ("R21", "Liq", "S_IC"): (self.Ci["S_pro"] - self.Ci["X_PHA"]) * mw_c,
            ("R21", "Liq", "S_IN"): 0,
            ("R21", "Liq", "S_IP"): (self.Y_PO4 * self.Pi["X_PP"]) * mw_p,
            ("R21", "Liq", "S_I"): 0,
            ("R21", "Liq", "X_ch"): 0,
            ("R21", "Liq", "X_pr"): 0,
            ("R21", "Liq", "X_li"): 0,
            ("R21", "Liq", "X_su"): 0,
            ("R21", "Liq", "X_aa"): 0,
            ("R21", "Liq", "X_fa"): 0,
            ("R21", "Liq", "X_c4"): 0,
            ("R21", "Liq", "X_pro"): 0,
            ("R21", "Liq", "X_ac"): 0,
            ("R21", "Liq", "X_h2"): 0,
            ("R21", "Liq", "X_I"): 0,
            ("R21", "Liq", "X_PHA"): 1,
            ("R21", "Liq", "X_PP"): -self.Y_PO4,
            ("R21", "Liq", "X_PAO"): 0,
            ("R21", "Liq", "S_K"): self.Y_PO4 * self.K_XPP,
            ("R21", "Liq", "S_Mg"): self.Y_PO4 * self.Mg_XPP,
            # R22: Storage of S_ac in X_PHA
            ("R22", "Liq", "H2O"): 0,
            ("R22", "Liq", "S_su"): 0,
            ("R22", "Liq", "S_aa"): 0,
            ("R22", "Liq", "S_fa"): 0,
            ("R22", "Liq", "S_va"): 0,
            ("R22", "Liq", "S_bu"): 0,
            ("R22", "Liq", "S_pro"): 0,
            ("R22", "Liq", "S_ac"): -1,
            ("R22", "Liq", "S_h2"): 0,
            ("R22", "Liq", "S_ch4"): 0,
            ("R22", "Liq", "S_IC"): (self.Ci["S_ac"] - self.Ci["X_PHA"]) * mw_c,
            ("R22", "Liq", "S_IN"): 0,
            ("R22", "Liq", "S_IP"): (self.Y_PO4 * self.Pi["X_PP"]) * mw_p,
            ("R22", "Liq", "S_I"): 0,
            ("R22", "Liq", "X_ch"): 0,
            ("R22", "Liq", "X_pr"): 0,
            ("R22", "Liq", "X_li"): 0,
            ("R22", "Liq", "X_su"): 0,
            ("R22", "Liq", "X_aa"): 0,
            ("R22", "Liq", "X_fa"): 0,
            ("R22", "Liq", "X_c4"): 0,
            ("R22", "Liq", "X_pro"): 0,
            ("R22", "Liq", "X_ac"): 0,
            ("R22", "Liq", "X_h2"): 0,
            ("R22", "Liq", "X_I"): 0,
            ("R22", "Liq", "X_PHA"): 1,
            ("R22", "Liq", "X_PP"): -self.Y_PO4,
            ("R22", "Liq", "X_PAO"): 0,
            ("R22", "Liq", "S_K"): self.Y_PO4 * self.K_XPP,
            ("R22", "Liq", "S_Mg"): self.Y_PO4 * self.Mg_XPP,
            # R23: Lysis of X_PAO
            ("R23", "Liq", "H2O"): 0,
            ("R23", "Liq", "S_su"): 0,
            ("R23", "Liq", "S_aa"): 0,
            ("R23", "Liq", "S_fa"): 0,
            ("R23", "Liq", "S_va"): 0,
            ("R23", "Liq", "S_bu"): 0,
            ("R23", "Liq", "S_pro"): 0,
            ("R23", "Liq", "S_ac"): 0,
            ("R23", "Liq", "S_h2"): 0,
            ("R23", "Liq", "S_ch4"): 0,
            ("R23", "Liq", "S_IC"): (
                self.Ci["X_PAO"]
                - self.f_ch_xb * self.Ci["X_ch"]
                - self.f_pr_xb * self.Ci["X_pr"]
                - self.f_li_xb * self.Ci["X_li"]
                - self.f_xi_xb * self.Ci["X_I"]
            )
            * mw_c,
            ("R23", "Liq", "S_IN"): (
                self.Ni["X_PAO"]
                - self.f_pr_xb * self.Ni["X_pr"]
                - self.f_xi_xb * self.Ni["X_I"]
            )
            * mw_n,
            ("R23", "Liq", "S_IP"): (
                self.Pi["X_PAO"]
                - self.f_li_xb * self.Pi["X_li"]
                - self.f_xi_xb * self.Pi["X_I"]
            )
            * mw_p,
            ("R23", "Liq", "S_I"): self.f_si_xb,
            ("R23", "Liq", "X_ch"): self.f_ch_xb,
            ("R23", "Liq", "X_pr"): self.f_pr_xb,
            ("R23", "Liq", "X_li"): self.f_li_xb,
            ("R23", "Liq", "X_su"): 0,
            ("R23", "Liq", "X_aa"): 0,
            ("R23", "Liq", "X_fa"): 0,
            ("R23", "Liq", "X_c4"): 0,
            ("R23", "Liq", "X_pro"): 0,
            ("R23", "Liq", "X_ac"): 0,
            ("R23", "Liq", "X_h2"): 0,
            ("R23", "Liq", "X_I"): self.f_xi_xb,
            ("R23", "Liq", "X_PHA"): 0,
            ("R23", "Liq", "X_PP"): 0,
            ("R23", "Liq", "X_PAO"): -1,
            ("R23", "Liq", "S_K"): 0,
            ("R23", "Liq", "S_Mg"): 0,
            # R24: Lysis of X_PP
            ("R24", "Liq", "H2O"): 0,
            ("R24", "Liq", "S_su"): 0,
            ("R24", "Liq", "S_aa"): 0,
            ("R24", "Liq", "S_fa"): 0,
            ("R24", "Liq", "S_va"): 0,
            ("R24", "Liq", "S_bu"): 0,
            ("R24", "Liq", "S_pro"): 0,
            ("R24", "Liq", "S_ac"): 0,
            ("R24", "Liq", "S_h2"): 0,
            ("R24", "Liq", "S_ch4"): 0,
            ("R24", "Liq", "S_IC"): 0,
            ("R24", "Liq", "S_IN"): 0,
            ("R24", "Liq", "S_IP"): self.Pi["X_PP"] * mw_p,
            ("R24", "Liq", "S_I"): 0,
            ("R24", "Liq", "X_ch"): 0,
            ("R24", "Liq", "X_pr"): 0,
            ("R24", "Liq", "X_li"): 0,
            ("R24", "Liq", "X_su"): 0,
            ("R24", "Liq", "X_aa"): 0,
            ("R24", "Liq", "X_fa"): 0,
            ("R24", "Liq", "X_c4"): 0,
            ("R24", "Liq", "X_pro"): 0,
            ("R24", "Liq", "X_ac"): 0,
            ("R24", "Liq", "X_h2"): 0,
            ("R24", "Liq", "X_I"): 0,
            ("R24", "Liq", "X_PHA"): 0,
            ("R24", "Liq", "X_PP"): -1,
            ("R24", "Liq", "X_PAO"): 0,
            ("R24", "Liq", "S_K"): self.K_XPP,
            ("R24", "Liq", "S_Mg"): self.Mg_XPP,
            # R25: Lysis of X_PHA
            ("R25", "Liq", "H2O"): 0,
            ("R25", "Liq", "S_su"): 0,
            ("R25", "Liq", "S_aa"): 0,
            ("R25", "Liq", "S_fa"): 0,
            ("R25", "Liq", "S_va"): self.f_va_PHA,
            ("R25", "Liq", "S_bu"): self.f_bu_PHA,
            ("R25", "Liq", "S_pro"): self.f_pro_PHA,
            ("R25", "Liq", "S_ac"): self.f_ac_PHA,
            ("R25", "Liq", "S_h2"): 0,
            ("R25", "Liq", "S_ch4"): 0,
            ("R25", "Liq", "S_IC"): (
                self.Ci["X_PHA"]
                - self.f_va_PHA * self.Ci["S_va"]
                - self.f_bu_PHA * self.Ci["S_bu"]
                - self.f_pro_PHA * self.Ci["S_pro"]
                - self.f_ac_PHA * self.Ci["S_ac"]
            )
            * mw_c,
            ("R25", "Liq", "S_IN"): 0,
            ("R25", "Liq", "S_IP"): 0,
            ("R25", "Liq", "S_I"): 0,
            ("R25", "Liq", "X_ch"): 0,
            ("R25", "Liq", "X_pr"): 0,
            ("R25", "Liq", "X_li"): 0,
            ("R25", "Liq", "X_su"): 0,
            ("R25", "Liq", "X_aa"): 0,
            ("R25", "Liq", "X_fa"): 0,
            ("R25", "Liq", "X_c4"): 0,
            ("R25", "Liq", "X_pro"): 0,
            ("R25", "Liq", "X_ac"): 0,
            ("R25", "Liq", "X_h2"): 0,
            ("R25", "Liq", "X_I"): 0,
            ("R25", "Liq", "X_PHA"): -1,
            ("R25", "Liq", "X_PP"): 0,
            ("R25", "Liq", "X_PAO"): 0,
            ("R25", "Liq", "S_K"): 0,
            ("R25", "Liq", "S_Mg"): 0,
        }

        for R in self.rate_reaction_idx:
            self.rate_reaction_stoichiometry[R, "Liq", "S_cat"] = 0
            self.rate_reaction_stoichiometry[R, "Liq", "S_an"] = 0

        # Fix all the variables we just created
        for v in self.component_objects(pyo.Var, descend_into=False):
            v.fix()

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties(
            {
                "reaction_rate": {"method": "_rxn_rate"},
            }
        )
        obj.define_custom_properties(
            {
                "conc_mol_co2": {"method": "_rxn_rate"},
                "I": {"method": "_I"},
            }
        )
        obj.add_default_units(
            {
                "time": pyo.units.s,
                "length": pyo.units.m,
                "mass": pyo.units.kg,
                "amount": pyo.units.kmol,
                "temperature": pyo.units.K,
            }
        )


class ModifiedADM1ReactionScaler(CustomScalerBase):
    """
    Scaler for the Modified Anaerobic Digestion Model No.1 reaction package.
    Variables are scaled by their default scaling factor (if no user input provided), and constraints
    are scaled using the inverse maximum scheme.
    """

    # TODO: Revisit this scaling factor
    DEFAULT_SCALING_FACTORS = {
        "reaction_rate": 1e2,
        "I": 1e1,
    }

    def variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        for r in model.params.rate_reaction_idx:
            self.scale_variable_by_default(model.I[r], overwrite=overwrite)

        if model.is_property_constructed("reaction_rate"):
            for j in model.reaction_rate.values():
                self.scale_variable_by_default(j, overwrite=overwrite)

    def constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        # TODO: Revisit these scaling methodologies
        # Consider other schemes, scale_constraint_by_default, or scale_constraints_by_jacobian_norm
        if model.is_property_constructed("rate_expression"):
            for j in model.rate_expression.values():
                self.scale_constraint_by_nominal_value(
                    j,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )
        if model.is_property_constructed("Dissociation"):
            self.scale_constraint_by_nominal_value(
                model.Dissociation,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )
        if model.is_property_constructed("CO2_acid_base_equilibrium"):
            self.scale_constraint_by_nominal_value(
                model.CO2_acid_base_equilibrium,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )
        if model.is_property_constructed("IN_acid_base_equilibrium"):
            self.scale_constraint_by_nominal_value(
                model.IN_acid_base_equilibrium,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )
        if model.is_property_constructed("pH_calc"):
            self.scale_constraint_by_nominal_value(
                model.pH_calc,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )
        if model.is_property_constructed("concentration_of_va"):
            self.scale_constraint_by_nominal_value(
                model.concentration_of_va,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )
        if model.is_property_constructed("concentration_of_bu"):
            self.scale_constraint_by_nominal_value(
                model.concentration_of_bu,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )
        if model.is_property_constructed("concentration_of_pro"):
            self.scale_constraint_by_nominal_value(
                model.concentration_of_pro,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )
        if model.is_property_constructed("concentration_of_ac"):
            self.scale_constraint_by_nominal_value(
                model.concentration_of_ac,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )
        if model.is_property_constructed("concentration_of_hco3"):
            self.scale_constraint_by_nominal_value(
                model.concentration_of_hco3,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )
        if model.is_property_constructed("concentration_of_nh3"):
            self.scale_constraint_by_nominal_value(
                model.concentration_of_nh3,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )
        if model.is_property_constructed("concentration_of_co2"):
            self.scale_constraint_by_nominal_value(
                model.concentration_of_co2,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )
        if model.is_property_constructed("concentration_of_nh4"):
            self.scale_constraint_by_nominal_value(
                model.concentration_of_nh4,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )
        if model.is_property_constructed("concentration_of_Mg"):
            self.scale_constraint_by_nominal_value(
                model.concentration_of_Mg,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )
        if model.is_property_constructed("concentration_of_K"):
            self.scale_constraint_by_nominal_value(
                model.concentration_of_K,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )
        if model.is_property_constructed("S_H_cons"):
            self.scale_constraint_by_nominal_value(
                model.S_H_cons,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )
        if model.is_property_constructed("I_fun"):
            for r in model.params.rate_reaction_idx:
                self.scale_constraint_by_nominal_value(
                    model.I_fun[r],
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )


class _ModifiedADM1ReactionBlock(ReactionBlockBase):
    """
    This Class contains methods which should be applied to Reaction Blocks as a
    whole, rather than individual elements of indexed Reaction Blocks.
    """

    default_scaler = ModifiedADM1ReactionScaler

    def initialize(self, outlvl=idaeslog.NOTSET, **kwargs):
        """
        Initialization routine for reaction package.

        Keyword Arguments:
            outlvl : sets output level of initialization routine

        Returns:
            None
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="properties")
        init_log.info("Initialization Complete.")


@declare_process_block_class(
    "ModifiedADM1ReactionBlock", block_class=_ModifiedADM1ReactionBlock
)
class ModifiedADM1ReactionBlockData(ReactionBlockDataBase):
    """
    ReactionBlock for ADM1.
    """

    def build(self):
        """
        Callable method for Block construction
        """
        super().build()

        # Create references to state vars
        # Concentration
        add_object_reference(self, "conc_mass_comp_ref", self.state_ref.conc_mass_comp)
        add_object_reference(self, "temperature", self.state_ref.temperature)

        # Initial values of rates of reaction derived from Flores-Alsina 2016 GitHub
        self.rates = {
            "R1": 1.651e-04,
            "R2": 1.723e-04,
            "R3": 2.290e-04,
            "R4": 5.312e-06,
            "R5": 5.176e-06,
            "R6": 6.436e-06,
            "R7": 1.074e-06,
            "R8": 1.790e-06,
            "R9": 2.355e-06,
            "R10": 5.531e-06,
            "R11": 4.472e-06,
            "R12": 1.294e-07,
            "R13": 1.009e-07,
            "R14": 9.406e-08,
            "R15": 4.185e-08,
            "R16": 2.294e-08,
            "R17": 1.259e-07,
            "R18": 6.535e-08,
            "R19": 1.507e-06,
            "R20": 2.078e-06,
            "R21": 2.630e-06,
            "R22": 1.195e-05,
            "R23": 1.901e-06,
            "R24": 1.481e-09,
            "R25": 1.425e-06,
        }

    # Rate of reaction method
    def _rxn_rate(self):
        self.reaction_rate = pyo.Var(
            self.params.rate_reaction_idx,
            initialize=self.rates,
            domain=pyo.NonNegativeReals,
            doc="Rate of reaction",
            units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
        )
        self.I = pyo.Var(
            self.params.rate_reaction_idx,
            initialize=1,
            bounds=(1e-8, 10),
            doc="Process inhibition term",
            units=pyo.units.dimensionless,
        )
        self.pKW = pyo.Var(
            initialize=14,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Water dissociation constant",
        )
        self.pK_a_co2 = pyo.Var(
            initialize=6.35,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Carbon dioxide acid-base equilibrium constant",
        )
        self.pK_a_IN = pyo.Var(
            initialize=9.25,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Inorganic nitrogen acid-base equilibrium constant",
        )
        self.conc_mass_va = pyo.Var(
            initialize=0.01159624,
            domain=pyo.NonNegativeReals,
            doc="mass concentration of va-",
            units=pyo.units.kg / pyo.units.m**3,
        )
        self.conc_mass_bu = pyo.Var(
            initialize=0.0132208,
            domain=pyo.NonNegativeReals,
            doc="mass concentration of bu-",
            units=pyo.units.kg / pyo.units.m**3,
        )
        self.conc_mass_pro = pyo.Var(
            initialize=0.015742,
            domain=pyo.NonNegativeReals,
            doc="mass concentration of pro-",
            units=pyo.units.kg / pyo.units.m**3,
        )
        self.conc_mass_ac = pyo.Var(
            initialize=0.1972,
            domain=pyo.NonNegativeReals,
            doc="mass concentration of ac-",
            units=pyo.units.kg / pyo.units.m**3,
        )
        self.conc_mol_hco3 = pyo.Var(
            initialize=0.142777,
            domain=pyo.NonNegativeReals,
            doc="molar concentration of hco3",
            units=pyo.units.kmol / pyo.units.m**3,
        )
        self.conc_mol_nh3 = pyo.Var(
            initialize=0.004,
            domain=pyo.NonNegativeReals,
            doc="molar concentration of nh3",
            units=pyo.units.kmol / pyo.units.m**3,
        )
        self.conc_mol_co2 = pyo.Var(
            initialize=0.0099,
            domain=pyo.NonNegativeReals,
            doc="molar concentration of co2",
            units=pyo.units.kmol / pyo.units.m**3,
        )
        self.conc_mol_nh4 = pyo.Var(
            initialize=0.1261,
            domain=pyo.NonNegativeReals,
            doc="molar concentration of nh4",
            units=pyo.units.kmol / pyo.units.m**3,
        )
        self.conc_mol_Mg = pyo.Var(
            initialize=4.5822e-05,
            domain=pyo.NonNegativeReals,
            doc="molar concentration of Mg+2",
            units=pyo.units.kmol / pyo.units.m**3,
        )
        self.conc_mol_K = pyo.Var(
            initialize=0.010934,
            domain=pyo.NonNegativeReals,
            doc="molar concentration of K+",
            units=pyo.units.kmol / pyo.units.m**3,
        )
        self.S_H = pyo.Var(
            initialize=3.4e-8,
            domain=pyo.NonNegativeReals,
            doc="molar concentration of H",
            units=pyo.units.kmol / pyo.units.m**3,
        )
        self.pH = pyo.Var(
            initialize=7,
            bounds=(0, 14),
            domain=pyo.PositiveReals,
            doc="pH of solution",
            units=pyo.units.dimensionless,
        )

        def Dissociation_rule(self, t):
            return pyo.log(10**-self.pKW) == (
                pyo.log(1e-14)
                + 55900
                / pyo.units.mole
                * pyo.units.joule
                / (Constants.gas_constant)
                * ((1 / self.params.temperature_ref) - (1 / self.temperature))
            )

        self.Dissociation = pyo.Constraint(
            rule=Dissociation_rule,
            doc="Water dissociation constant constraint",
        )

        def CO2_acid_base_equilibrium_rule(self, t):
            return pyo.log(10**-self.pK_a_co2) == (
                pyo.log(10**-6.35)
                + 7646
                / pyo.units.mole
                * pyo.units.joule
                / (Constants.gas_constant)
                * ((1 / self.params.temperature_ref) - (1 / self.temperature))
            )

        self.CO2_acid_base_equilibrium = pyo.Constraint(
            rule=CO2_acid_base_equilibrium_rule,
            doc="Carbon dioxide acid-base equilibrium constraint",
        )

        def IN_acid_base_equilibrium_rule(self, t):
            return pyo.log(10**-self.pK_a_IN) == (
                pyo.log(10**-9.25)
                + 51965
                / pyo.units.mole
                * pyo.units.joule
                / (Constants.gas_constant)
                * ((1 / self.params.temperature_ref) - (1 / self.temperature))
            )

        self.IN_acid_base_equilibrium = pyo.Constraint(
            rule=IN_acid_base_equilibrium_rule,
            doc="Nitrogen acid-base equilibrium constraint",
        )

        def rule_pH(self):
            return self.pH == -pyo.log10(self.S_H / (pyo.units.kmole / pyo.units.m**3))

        self.pH_calc = pyo.Constraint(rule=rule_pH, doc="pH of solution")

        def concentration_of_va_rule(self):
            return (
                self.conc_mass_va * (1 + self.S_H / self.params.K_a_va)
                == self.conc_mass_comp_ref["S_va"]
            )

        self.concentration_of_va = pyo.Constraint(
            rule=concentration_of_va_rule,
            doc="constraint concentration of va-",
        )

        def concentration_of_bu_rule(self):
            return (
                self.conc_mass_bu * (1 + self.S_H / self.params.K_a_bu)
                == self.conc_mass_comp_ref["S_bu"]
            )

        self.concentration_of_bu = pyo.Constraint(
            rule=concentration_of_bu_rule,
            doc="constraint concentration of bu-",
        )

        def concentration_of_pro_rule(self):
            return (
                self.conc_mass_pro * (1 + self.S_H / self.params.K_a_pro)
                == self.conc_mass_comp_ref["S_pro"]
            )

        self.concentration_of_pro = pyo.Constraint(
            rule=concentration_of_pro_rule,
            doc="constraint concentration of pro-",
        )

        def concentration_of_ac_rule(self):
            return (
                self.conc_mass_ac * (1 + self.S_H / self.params.K_a_ac)
                == self.conc_mass_comp_ref["S_ac"]
            )

        self.concentration_of_ac = pyo.Constraint(
            rule=concentration_of_ac_rule,
            doc="constraint concentration of ac-",
        )

        def concentration_of_hco3_rule(self):
            return (
                self.pK_a_co2
                == pyo.log10(self.conc_mol_co2 / (pyo.units.kmole / pyo.units.m**3))
                - pyo.log10(self.conc_mol_hco3 / (pyo.units.kmole / pyo.units.m**3))
                + self.pH
            )

        self.concentration_of_hco3 = pyo.Constraint(
            rule=concentration_of_hco3_rule,
            doc="constraint concentration of hco3",
        )

        def concentration_of_nh3_rule(self):
            return (
                self.pK_a_IN
                == pyo.log10(self.conc_mol_nh4 / (pyo.units.kmole / pyo.units.m**3))
                - pyo.log10(self.conc_mol_nh3 / (pyo.units.kmole / pyo.units.m**3))
                + self.pH
            )

        self.concentration_of_nh3 = pyo.Constraint(
            rule=concentration_of_nh3_rule,
            doc="constraint concentration of nh3",
        )

        # TO DO: use correct conversion number
        def concentration_of_co2_rule(self):
            return (
                self.conc_mol_co2
                == self.conc_mass_comp_ref["S_IC"] / mw_c - self.conc_mol_hco3
            )

        self.concentration_of_co2 = pyo.Constraint(
            rule=concentration_of_co2_rule,
            doc="constraint concentration of co2",
        )

        def concentration_of_nh4_rule(self):
            return (
                self.conc_mol_nh4
                == self.conc_mass_comp_ref["S_IN"] / mw_n - self.conc_mol_nh3
            )

        self.concentration_of_nh4 = pyo.Constraint(
            rule=concentration_of_nh4_rule,
            doc="constraint concentration of nh4",
        )

        def concentration_of_Mg_rule(self):
            return self.conc_mol_Mg == self.conc_mass_comp_ref["X_PP"] / (
                300.41 * pyo.units.kg / pyo.units.kmol
            )

        self.concentration_of_Mg = pyo.Constraint(
            rule=concentration_of_Mg_rule,
            doc="constraint concentration of Mg",
        )

        def concentration_of_K_rule(self):
            return self.conc_mol_K == self.conc_mass_comp_ref["X_PP"] / (
                300.41 * pyo.units.kg / pyo.units.kmol
            )

        self.concentration_of_K = pyo.Constraint(
            rule=concentration_of_K_rule,
            doc="constraint concentration of K",
        )

        def S_H_rule(self):
            return (
                self.state_ref.cations
                + self.conc_mol_nh4
                + self.conc_mol_Mg
                + self.conc_mol_K
                + self.S_H
                - self.conc_mol_hco3
                - self.conc_mass_ac / (64 * pyo.units.kg / pyo.units.kmol)
                - self.conc_mass_pro / (112 * pyo.units.kg / pyo.units.kmol)
                - self.conc_mass_bu / (160 * pyo.units.kg / pyo.units.kmol)
                - self.conc_mass_va / (208 * pyo.units.kg / pyo.units.kmol)
                - 10 ** (self.pH - self.pKW) * (pyo.units.kmole / pyo.units.m**3)
                - self.state_ref.anions
                == 0
            )

        self.S_H_cons = pyo.Constraint(
            rule=S_H_rule,
            doc="constraint concentration of H",
        )

        def rule_I_IN_lim(self):
            return 1 / (
                1 + self.params.K_S_IN / (self.conc_mass_comp_ref["S_IN"] / mw_n)
            )

        self.I_IN_lim = pyo.Expression(
            rule=rule_I_IN_lim,
            doc="Inhibition function related to secondary substrate; inhibit uptake when inorganic nitrogen S_IN~ 0",
        )

        def rule_I_IP_lim(self):
            return 1 / (
                1 + self.params.K_S_IP / (self.conc_mass_comp_ref["S_IP"] / mw_p)
            )

        self.I_IP_lim = pyo.Expression(
            rule=rule_I_IP_lim,
            doc="Inhibition function related to secondary substrate; inhibit uptake when inorganic phosphorus S_IP~ 0",
        )

        def rule_I_h2_fa(self):
            return 1 / (1 + self.conc_mass_comp_ref["S_h2"] / self.params.K_I_h2_fa)

        self.I_h2_fa = pyo.Expression(
            rule=rule_I_h2_fa,
            doc="hydrogen inhibition attributed to long chain fatty acids",
        )

        def rule_I_h2_c4(self):
            return 1 / (1 + self.conc_mass_comp_ref["S_h2"] / self.params.K_I_h2_c4)

        self.I_h2_c4 = pyo.Expression(
            rule=rule_I_h2_c4,
            doc="hydrogen inhibition attributed to valerate and butyrate uptake",
        )

        def rule_I_h2_pro(self):
            return 1 / (1 + self.conc_mass_comp_ref["S_h2"] / self.params.K_I_h2_pro)

        self.I_h2_pro = pyo.Expression(
            rule=rule_I_h2_pro,
            doc="hydrogen inhibition attributed to propionate uptake",
        )

        # TODO: revisit Z_h2s value if we have ref state for S_h2s (currently assumed to be 0)
        def rule_I_h2s_ac(self):
            return 1 / (1 + self.params.Z_h2s / self.params.K_I_h2s_ac)

        self.I_h2s_ac = pyo.Expression(
            rule=rule_I_h2s_ac,
            doc="hydrogen sulfide inhibition attributed to acetate uptake",
        )

        def rule_I_h2s_c4(self):
            return 1 / (1 + self.params.Z_h2s / self.params.K_I_h2s_c4)

        self.I_h2s_c4 = pyo.Expression(
            rule=rule_I_h2s_c4,
            doc="hydrogen sulfide inhibition attributed to valerate and butyrate uptake",
        )

        def rule_I_h2s_h2(self):
            return 1 / (1 + self.params.Z_h2s / self.params.K_I_h2s_h2)

        self.I_h2s_h2 = pyo.Expression(
            rule=rule_I_h2s_h2,
            doc="hydrogen sulfide inhibition attributed to hydrogen uptake",
        )

        def rule_I_h2s_pro(self):
            return 1 / (1 + self.params.Z_h2s / self.params.K_I_h2s_pro)

        self.I_h2s_pro = pyo.Expression(
            rule=rule_I_h2s_pro,
            doc="hydrogen sulfide inhibition attributed to propionate uptake",
        )

        def rule_I_nh3(self):
            return 1 / (1 + self.conc_mol_nh3 / self.params.K_I_nh3)

        self.I_nh3 = pyo.Expression(
            rule=rule_I_nh3, doc="ammonia inibition attributed to acetate uptake"
        )

        def rule_I_pH_aa(self):
            return (
                -3
                * (
                    smooth_max(0, self.params.pH_UL_aa - self.pH, eps=1e-8)
                    / (self.params.pH_UL_aa - self.params.pH_LL_aa)
                )
                ** 2
            )

        self.I_pH_aa = pyo.Expression(
            rule=rule_I_pH_aa,
            doc="pH inhibition of amino-acid-utilizing microorganisms",
        )

        def rule_I_pH_ac(self):
            return (
                -3
                * (
                    smooth_max(0, self.params.pH_UL_ac - self.pH, eps=1e-8)
                    / (self.params.pH_UL_ac - self.params.pH_LL_ac)
                )
                ** 2
            )

        self.I_pH_ac = pyo.Expression(
            rule=rule_I_pH_ac, doc="pH inhibition of acetate-utilizing microorganisms"
        )

        def rule_I_pH_h2(self):
            return (
                -3
                * (
                    smooth_max(0, self.params.pH_UL_h2 - self.pH, eps=1e-8)
                    / (self.params.pH_UL_h2 - self.params.pH_LL_h2)
                )
                ** 2
            )

        self.I_pH_h2 = pyo.Expression(
            rule=rule_I_pH_h2, doc="pH inhibition of hydrogen-utilizing microorganisms"
        )

        def rule_I(self, r):
            if r == "R4" or r == "R5":
                return (
                    self.I[r] == pyo.exp(self.I_pH_aa) * self.I_IN_lim * self.I_IP_lim
                )
            elif r == "R6":
                return (
                    self.I[r]
                    == pyo.exp(self.I_pH_aa)
                    * self.I_IN_lim
                    * self.I_h2_fa
                    * self.I_IP_lim
                )
            elif r == "R7" or r == "R8":
                return (
                    self.I[r]
                    == pyo.exp(self.I_pH_aa)
                    * self.I_IN_lim
                    * self.I_h2_c4
                    * self.I_IP_lim
                )
            elif r == "R9":
                return (
                    self.I[r]
                    == pyo.exp(self.I_pH_aa)
                    * self.I_IN_lim
                    * self.I_h2_pro
                    * self.I_IP_lim
                )
            elif r == "R10":
                return (
                    self.I[r]
                    == pyo.exp(self.I_pH_ac)
                    * self.I_IN_lim
                    * self.I_nh3
                    * self.I_IP_lim
                )
            elif r == "R11":
                return (
                    self.I[r] == pyo.exp(self.I_pH_h2) * self.I_IN_lim * self.I_IP_lim
                )
            else:
                return self.I[r] == 1.0

        self.I_fun = pyo.Constraint(
            self.params.rate_reaction_idx,
            rule=rule_I,
            doc="Process inhibition functions",
        )

        try:

            def rate_expression_rule(b, r):
                if r == "R1":
                    # R1: Hydrolysis of carbohydrates
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_hyd_ch * b.conc_mass_comp_ref["X_ch"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R2":
                    # R2: Hydrolysis of proteins
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_hyd_pr * b.conc_mass_comp_ref["X_pr"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R3":
                    # R3: Hydrolysis of lipids
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_hyd_li * b.conc_mass_comp_ref["X_li"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R4":
                    # R4: Uptake of sugars
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_m_su
                        * b.conc_mass_comp_ref["S_su"]
                        / (b.params.K_S_su + b.conc_mass_comp_ref["S_su"])
                        * b.conc_mass_comp_ref["X_su"]
                        * b.I[r],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R5":
                    # R5: Uptake of amino acids
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_m_aa
                        * b.conc_mass_comp_ref["S_aa"]
                        / (b.params.K_S_aa + b.conc_mass_comp_ref["S_aa"])
                        * b.conc_mass_comp_ref["X_aa"]
                        * b.I[r],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R6":
                    # R6: Uptake of long chain fatty acids (LCFAs)
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_m_fa
                        * b.conc_mass_comp_ref["S_fa"]
                        / (b.params.K_S_fa + b.conc_mass_comp_ref["S_fa"])
                        * b.conc_mass_comp_ref["X_fa"]
                        * b.I[r],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R7":
                    # R7: Uptake of valerate
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_m_c4
                        * b.conc_mass_comp_ref["S_va"]
                        / (b.params.K_S_c4 + b.conc_mass_comp_ref["S_va"])
                        * b.conc_mass_comp_ref["X_c4"]
                        * (
                            b.conc_mass_comp_ref["S_va"]
                            / (
                                b.conc_mass_comp_ref["S_va"]
                                + b.conc_mass_comp_ref["S_bu"]
                                + 1e-10 * pyo.units.kg / pyo.units.m**3
                            )
                        )
                        * b.I[r]
                        * b.I_h2s_c4,
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R8":
                    # R8:  Uptake of butyrate
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_m_c4
                        * b.conc_mass_comp_ref["S_bu"]
                        / (b.params.K_S_c4 + b.conc_mass_comp_ref["S_bu"])
                        * b.conc_mass_comp_ref["X_c4"]
                        * (
                            b.conc_mass_comp_ref["S_bu"]
                            / (
                                b.conc_mass_comp_ref["S_va"]
                                + b.conc_mass_comp_ref["S_bu"]
                                + 1e-10 * pyo.units.kg / pyo.units.m**3
                            )
                        )
                        * b.I[r]
                        * b.I_h2s_c4,
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R9":
                    # R9: Uptake of propionate
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_m_pro
                        * b.conc_mass_comp_ref["S_pro"]
                        / (b.params.K_S_pro + b.conc_mass_comp_ref["S_pro"])
                        * b.conc_mass_comp_ref["X_pro"]
                        * b.I[r]
                        * b.I_h2s_pro,
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R10":
                    # R10: Uptake of acetate
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_m_ac
                        * b.conc_mass_comp_ref["S_ac"]
                        / (b.params.K_S_ac + b.conc_mass_comp_ref["S_ac"])
                        * b.conc_mass_comp_ref["X_ac"]
                        * b.I[r]
                        * b.I_h2s_ac,
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R11":
                    # R11: Uptake of hydrogen
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_m_h2
                        * b.conc_mass_comp_ref["S_h2"]
                        / (b.params.K_S_h2 + b.conc_mass_comp_ref["S_h2"])
                        * b.conc_mass_comp_ref["X_h2"]
                        * b.I[r]
                        * b.I_h2s_h2,
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R12":
                    # R12: Decay of X_su
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_dec_X_su * b.conc_mass_comp_ref["X_su"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R13":
                    # R13: Decay of X_aa
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_dec_X_aa * b.conc_mass_comp_ref["X_aa"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R14":
                    # R14: Decay of X_fa
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_dec_X_fa * b.conc_mass_comp_ref["X_fa"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R15":
                    # R15: Decay of X_c4
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_dec_X_c4 * b.conc_mass_comp_ref["X_c4"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R16":
                    # R16: Decay of X_pro
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_dec_X_pro * b.conc_mass_comp_ref["X_pro"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R17":
                    # R17: Decay of X_ac
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_dec_X_ac * b.conc_mass_comp_ref["X_ac"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R18":
                    # R18: Decay of X_h2
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_dec_X_h2 * b.conc_mass_comp_ref["X_h2"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R19":
                    # R19: Storage of S_va in X_PHA
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.q_PHA
                        * b.conc_mass_comp_ref["S_va"]
                        / (b.params.K_A + b.conc_mass_comp_ref["S_va"])
                        * b.conc_mass_comp_ref["X_PP"]
                        / (
                            b.params.K_PP * b.conc_mass_comp_ref["X_PAO"]
                            + b.conc_mass_comp_ref["X_PP"]
                        )
                        * b.conc_mass_comp_ref["X_PAO"]
                        * b.conc_mass_comp_ref["S_va"]
                        / (
                            b.conc_mass_comp_ref["S_va"]
                            + b.conc_mass_comp_ref["S_bu"]
                            + b.conc_mass_comp_ref["S_pro"]
                            + b.conc_mass_comp_ref["S_ac"]
                            + 1e-10 * pyo.units.kg / pyo.units.m**3
                        ),
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R20":
                    # R20: Storage of S_bu in X_PHA
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.q_PHA
                        * b.conc_mass_comp_ref["S_bu"]
                        / (b.params.K_A + b.conc_mass_comp_ref["S_bu"])
                        * b.conc_mass_comp_ref["X_PP"]
                        / (
                            b.params.K_PP * b.conc_mass_comp_ref["X_PAO"]
                            + b.conc_mass_comp_ref["X_PP"]
                        )
                        * b.conc_mass_comp_ref["X_PAO"]
                        * b.conc_mass_comp_ref["S_bu"]
                        / (
                            b.conc_mass_comp_ref["S_va"]
                            + b.conc_mass_comp_ref["S_bu"]
                            + b.conc_mass_comp_ref["S_pro"]
                            + b.conc_mass_comp_ref["S_ac"]
                            + 1e-10 * pyo.units.kg / pyo.units.m**3
                        ),
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R21":
                    # R21: Storage of S_pro in X_PHA
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.q_PHA
                        * b.conc_mass_comp_ref["S_pro"]
                        / (b.params.K_A + b.conc_mass_comp_ref["S_pro"])
                        * b.conc_mass_comp_ref["X_PP"]
                        / (
                            b.params.K_PP * b.conc_mass_comp_ref["X_PAO"]
                            + b.conc_mass_comp_ref["X_PP"]
                        )
                        * b.conc_mass_comp_ref["X_PAO"]
                        * b.conc_mass_comp_ref["S_pro"]
                        / (
                            b.conc_mass_comp_ref["S_va"]
                            + b.conc_mass_comp_ref["S_bu"]
                            + b.conc_mass_comp_ref["S_pro"]
                            + b.conc_mass_comp_ref["S_ac"]
                            + 1e-10 * pyo.units.kg / pyo.units.m**3
                        ),
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R22":
                    # R22: Storage of S_ac in X_PHA
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.q_PHA
                        * b.conc_mass_comp_ref["S_ac"]
                        / (b.params.K_A + b.conc_mass_comp_ref["S_ac"])
                        * b.conc_mass_comp_ref["X_PP"]
                        / (
                            b.params.K_PP * b.conc_mass_comp_ref["X_PAO"]
                            + b.conc_mass_comp_ref["X_PP"]
                        )
                        * b.conc_mass_comp_ref["X_PAO"]
                        * b.conc_mass_comp_ref["S_ac"]
                        / (
                            b.conc_mass_comp_ref["S_va"]
                            + b.conc_mass_comp_ref["S_bu"]
                            + b.conc_mass_comp_ref["S_pro"]
                            + b.conc_mass_comp_ref["S_ac"]
                            + 1e-10 * pyo.units.kg / pyo.units.m**3
                        ),
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R23":
                    # R23: Lysis of X_PAO
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.b_PAO * b.conc_mass_comp_ref["X_PAO"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R24":
                    # R24: Lysis of X_PP
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.b_PP * b.conc_mass_comp_ref["X_PP"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R25":
                    # R25: Lysis of X_PHA
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.b_PHA * b.conc_mass_comp_ref["X_PHA"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                else:
                    raise BurntToast()

            self.rate_expression = pyo.Constraint(
                self.params.rate_reaction_idx,
                rule=rate_expression_rule,
                doc="ADM1 rate expressions",
            )

        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.reaction_rate)
            self.del_component(self.rate_expression)
            raise

    def get_reaction_rate_basis(self):
        return MaterialFlowBasis.mass

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        for i, c in self.rates.items():
            iscale.set_scaling_factor(self.reaction_rate[i], 1 / c)

        iscale.set_scaling_factor(self.I, 1e1)
        iscale.set_scaling_factor(self.conc_mass_va, 1e2)
        iscale.set_scaling_factor(self.conc_mass_bu, 1e2)
        iscale.set_scaling_factor(self.conc_mass_pro, 1e2)
        iscale.set_scaling_factor(self.conc_mass_ac, 1e1)
        iscale.set_scaling_factor(self.conc_mol_hco3, 1e1)
        iscale.set_scaling_factor(self.conc_mol_nh3, 1e3)
        iscale.set_scaling_factor(self.conc_mol_co2, 1e3)
        iscale.set_scaling_factor(self.conc_mol_nh4, 1e1)
        iscale.set_scaling_factor(self.conc_mol_Mg, 1e3)
        iscale.set_scaling_factor(self.conc_mol_K, 1e3)
        iscale.set_scaling_factor(self.S_H, 1e5)
        iscale.set_scaling_factor(self.pKW, 1e0)
        iscale.set_scaling_factor(self.pK_a_co2, 1e0)
        iscale.set_scaling_factor(self.pK_a_IN, 1e0)
        iscale.set_scaling_factor(self.pH, 1e0)

        for i, c in self.rate_expression.items():
            iscale.constraint_scaling_transform(
                c,
                iscale.get_scaling_factor(
                    self.reaction_rate[i], default=1, warning=True
                ),
                overwrite=True,
            )
