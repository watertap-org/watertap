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
ADM1 reaction package.

References:

[1] D.J. Batstone, J. Keller, I. Angelidaki, S.V. Kalyuzhnyi, S.G. Pavlostathis,
A. Rozzi, W.T.M. Sanders, H. Siegrist, V.A. Vavilin,
The IWA Anaerobic Digestion Model No 1 (ADM1), Water Science and Technology.
45 (2002) 65–73. https://doi.org/10.2166/wst.2002.0292.

[2] C. Rosén, U. Jeppsson, Aspects on ADM1 Implementation within the BSM2 Framework,
Department of Industrial Electrical Engineering and Automation, Lund University, Lund, Sweden. (2006) 1–35.

"""

# Import Pyomo libraries
import pyomo.environ as pyo

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
__author__ = "Adam Atia, Alejandro Garciadiego, Xinhong Liu"

# Set up logger
_log = idaeslog.getLogger(__name__)

mw_n = 14 * pyo.units.kg / pyo.units.kmol
mw_c = 12 * pyo.units.kg / pyo.units.kmol


@declare_process_block_class("ADM1ReactionParameterBlock")
class ADM1ReactionParameterData(ReactionParameterBlock):
    """
    Property Parameter Block Class
    """

    def build(self):
        """
        Callable method for Block construction.
        """
        super().build()

        self._reaction_block_class = ADM1ReactionBlock

        # Reaction Index
        # Reaction names based on standard numbering in ADM1
        # R1:  Disintegration
        # R2:  Hydrolysis of carbohydrates
        # R3:  Hydrolysis of proteins
        # R4:  Hydrolysis of lipids
        # R5:  Uptake of sugars
        # R6:  Uptake of amino acids
        # R7:  Uptake of long chain fatty acids (LCFAs)
        # R8:  Uptake of valerate
        # R9:  Uptake of butyrate
        # R10: Uptake of propionate
        # R11: Uptake of acetate
        # R12: Uptake of hydrogen
        # R13: Decay of X_su
        # R14: Decay of X_aa
        # R15: Decay of X_fa
        # R16: Decay of X_c4
        # R17: Decay of X_pro
        # R18: Decay of X_ac
        # R19: Decay of X_h2

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
            ]
        )

        # Carbon content
        Ci_dict = {
            "S_su": 0.0313,
            "S_aa": 0.03,  # varies
            "S_fa": 0.0217,
            "S_va": 0.024,
            "S_bu": 0.025,
            "S_pro": 0.0268,
            "S_ac": 0.0313,
            "S_ch4": 0.0156,
            "S_I": 0.03,  # varies
            "X_c": 0.02786,  # varies
            "X_ch": 0.0313,
            "X_pr": 0.03,  # varies
            "X_li": 0.022,
            "X_su": 0.0313,
            "X_aa": 0.0313,
            "X_fa": 0.0313,
            "X_c4": 0.0313,
            "X_pro": 0.0313,
            "X_ac": 0.0313,
            "X_h2": 0.0313,
            "X_I": 0.03,  # varies
        }

        self.Ci = pyo.Var(
            Ci_dict.keys(),
            initialize=Ci_dict,
            units=pyo.units.kmol / pyo.units.kg,
            domain=pyo.PositiveReals,
            doc="Carbon content of component [kmole C/kg COD]",
        )

        # Stoichiometric Parameters (Table 6.1 in Batstone et al., 2002)
        self.f_sI_xc = pyo.Var(
            initialize=0.1,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Soluble inerts from composites",
        )
        self.f_xI_xc = pyo.Var(
            initialize=0.20,  # replacing 0.25 with 0.2 in accordance with Rosen & Jeppsson, 2006
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Particulate inerts from composites",
        )
        self.f_ch_xc = pyo.Var(
            initialize=0.2,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Carbohydrates from composites",
        )
        self.f_pr_xc = pyo.Var(
            initialize=0.2,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Proteins from composites",
        )
        self.f_li_xc = pyo.Var(
            initialize=0.30,  # replacing 0.25 with 0.3 in accordance with Rosen & Jeppsson, 2006
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Lipids from composites",
        )
        self.N_xc = pyo.Var(
            initialize=0.0376
            / 14,  # change from 0.002 to 0.0376/14 based on Rosen & Jeppsson, 2006
            units=pyo.units.kmol * pyo.units.kg**-1,
            domain=pyo.PositiveReals,
            doc="Nitrogen content of composites [kmole N/kg COD]",
        )
        self.N_I = pyo.Var(
            initialize=0.06
            / 14,  # change from 0.002 to 0.06/14 based on Rosen & Jeppsson, 2006
            units=pyo.units.kmol * pyo.units.kg**-1,
            domain=pyo.PositiveReals,
            doc="Nitrogen content of inerts [kmole N/kg COD]",
        )
        self.N_aa = pyo.Var(
            initialize=0.007,
            units=pyo.units.kmol * pyo.units.kg**-1,
            domain=pyo.PositiveReals,
            doc="Nitrogen in amino acids and proteins [kmole N/kg COD]",
        )
        self.N_bac = pyo.Var(
            initialize=0.08 / 14,
            units=pyo.units.kmol * pyo.units.kg**-1,
            # units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Nitrogen content in bacteria [kmole N/kg COD]",
        )
        self.f_fa_li = pyo.Var(
            initialize=0.95,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Fatty acids from lipids",
        )
        self.f_h2_su = pyo.Var(
            initialize=0.19,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Hydrogen from sugars",
        )
        self.f_bu_su = pyo.Var(
            initialize=0.13,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Butyrate from sugars",
        )
        self.f_pro_su = pyo.Var(
            initialize=0.27,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Propionate from sugars",
        )
        self.f_ac_su = pyo.Var(
            initialize=0.41,
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
        self.k_dis = pyo.Var(
            initialize=0.5,
            units=pyo.units.day**-1,
            domain=pyo.PositiveReals,
            doc="First-order kinetic parameter for disintegration",
        )
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
        self.temperature_ref = pyo.Param(
            within=pyo.PositiveReals,
            mutable=True,
            default=298.15,
            doc="Reference temperature",
            units=pyo.units.K,
        )

        # Reaction Stoichiometry
        # This is the stoichiometric part of the Peterson matrix in dict form.
        # See Table 3.1 in Batstone et al., 2002.

        # Exclude non-zero stoichiometric coefficients for S_IC initially since they depend on other stoichiometric coefficients.
        self.rate_reaction_stoichiometry = {
            # R1: Disintegration
            ("R1", "Liq", "H2O"): 0,
            ("R1", "Liq", "S_su"): 0,
            ("R1", "Liq", "S_aa"): 0,
            ("R1", "Liq", "S_fa"): 0,
            ("R1", "Liq", "S_va"): 0,
            ("R1", "Liq", "S_bu"): 0,
            ("R1", "Liq", "S_pro"): 0,
            ("R1", "Liq", "S_ac"): 0,
            ("R1", "Liq", "S_h2"): 0,
            ("R1", "Liq", "S_ch4"): 0,
            ("R1", "Liq", "S_IC"): (
                self.Ci["X_c"]
                - self.f_sI_xc * self.Ci["S_I"]
                - self.f_ch_xc * self.Ci["X_ch"]
                - self.f_pr_xc * self.Ci["X_pr"]
                - self.f_li_xc * self.Ci["X_li"]
                - self.f_xI_xc * self.Ci["X_I"]
            )
            * mw_c,
            ("R1", "Liq", "S_IN"): (
                self.N_xc
                - self.f_xI_xc * self.N_I
                - self.f_sI_xc * self.N_I
                - self.f_pr_xc * self.N_aa
            )
            * mw_n,
            ("R1", "Liq", "S_I"): self.f_sI_xc,
            ("R1", "Liq", "X_c"): -1,
            ("R1", "Liq", "X_ch"): self.f_ch_xc,
            ("R1", "Liq", "X_pr"): self.f_pr_xc,
            ("R1", "Liq", "X_li"): self.f_li_xc,
            ("R1", "Liq", "X_su"): 0,
            ("R1", "Liq", "X_aa"): 0,
            ("R1", "Liq", "X_fa"): 0,
            ("R1", "Liq", "X_c4"): 0,
            ("R1", "Liq", "X_pro"): 0,
            ("R1", "Liq", "X_ac"): 0,
            ("R1", "Liq", "X_h2"): 0,
            ("R1", "Liq", "X_I"): self.f_xI_xc,
            # R2: Hydrolysis of carbohydrates
            ("R2", "Liq", "H2O"): 0,
            ("R2", "Liq", "S_su"): 1,
            ("R2", "Liq", "S_aa"): 0,
            ("R2", "Liq", "S_fa"): 0,
            ("R2", "Liq", "S_va"): 0,
            ("R2", "Liq", "S_bu"): 0,
            ("R2", "Liq", "S_pro"): 0,
            ("R2", "Liq", "S_ac"): 0,
            ("R2", "Liq", "S_h2"): 0,
            ("R2", "Liq", "S_ch4"): 0,
            ("R2", "Liq", "S_IC"): 0,
            ("R2", "Liq", "S_IN"): 0,
            ("R2", "Liq", "S_I"): 0,
            ("R2", "Liq", "X_c"): 0,
            ("R2", "Liq", "X_ch"): -1,
            ("R2", "Liq", "X_pr"): 0,
            ("R2", "Liq", "X_li"): 0,
            ("R2", "Liq", "X_su"): 0,
            ("R2", "Liq", "X_aa"): 0,
            ("R2", "Liq", "X_fa"): 0,
            ("R2", "Liq", "X_c4"): 0,
            ("R2", "Liq", "X_pro"): 0,
            ("R2", "Liq", "X_ac"): 0,
            ("R2", "Liq", "X_h2"): 0,
            ("R2", "Liq", "X_I"): 0,
            # R3:  Hydrolysis of proteins
            ("R3", "Liq", "H2O"): 0,
            ("R3", "Liq", "S_su"): 0,
            ("R3", "Liq", "S_aa"): 1,
            ("R3", "Liq", "S_fa"): 0,
            ("R3", "Liq", "S_va"): 0,
            ("R3", "Liq", "S_bu"): 0,
            ("R3", "Liq", "S_pro"): 0,
            ("R3", "Liq", "S_ac"): 0,
            ("R3", "Liq", "S_h2"): 0,
            ("R3", "Liq", "S_ch4"): 0,
            ("R3", "Liq", "S_IC"): 0,
            ("R3", "Liq", "S_IN"): 0,
            ("R3", "Liq", "S_I"): 0,
            ("R3", "Liq", "X_c"): 0,
            ("R3", "Liq", "X_ch"): 0,
            ("R3", "Liq", "X_pr"): -1,
            ("R3", "Liq", "X_li"): 0,
            ("R3", "Liq", "X_su"): 0,
            ("R3", "Liq", "X_aa"): 0,
            ("R3", "Liq", "X_fa"): 0,
            ("R3", "Liq", "X_c4"): 0,
            ("R3", "Liq", "X_pro"): 0,
            ("R3", "Liq", "X_ac"): 0,
            ("R3", "Liq", "X_h2"): 0,
            ("R3", "Liq", "X_I"): 0,
            # R4:  Hydrolysis of lipids
            ("R4", "Liq", "H2O"): 0,
            ("R4", "Liq", "S_su"): 1 - self.f_fa_li,
            ("R4", "Liq", "S_aa"): 0,
            ("R4", "Liq", "S_fa"): self.f_fa_li,
            ("R4", "Liq", "S_va"): 0,
            ("R4", "Liq", "S_bu"): 0,
            ("R4", "Liq", "S_pro"): 0,
            ("R4", "Liq", "S_ac"): 0,
            ("R4", "Liq", "S_h2"): 0,
            ("R4", "Liq", "S_ch4"): 0,
            ("R4", "Liq", "S_IC"): (
                -(1 - self.f_fa_li) * self.Ci["S_su"]
                - self.f_fa_li * self.Ci["S_fa"]
                + self.Ci["X_li"]
            )
            * mw_c,
            ("R4", "Liq", "S_IN"): 0,
            ("R4", "Liq", "S_I"): 0,
            ("R4", "Liq", "X_c"): 0,
            ("R4", "Liq", "X_ch"): 0,
            ("R4", "Liq", "X_pr"): 0,
            ("R4", "Liq", "X_li"): -1,
            ("R4", "Liq", "X_su"): 0,
            ("R4", "Liq", "X_aa"): 0,
            ("R4", "Liq", "X_fa"): 0,
            ("R4", "Liq", "X_c4"): 0,
            ("R4", "Liq", "X_pro"): 0,
            ("R4", "Liq", "X_ac"): 0,
            ("R4", "Liq", "X_h2"): 0,
            ("R4", "Liq", "X_I"): 0,
            # R5:  Uptake of sugars
            ("R5", "Liq", "H2O"): 0,
            ("R5", "Liq", "S_su"): -1,
            ("R5", "Liq", "S_aa"): 0,
            ("R5", "Liq", "S_fa"): 0,
            ("R5", "Liq", "S_va"): 0,
            ("R5", "Liq", "S_bu"): (1 - self.Y_su) * self.f_bu_su,
            ("R5", "Liq", "S_pro"): (1 - self.Y_su) * self.f_pro_su,
            ("R5", "Liq", "S_ac"): (1 - self.Y_su) * self.f_ac_su,
            ("R5", "Liq", "S_h2"): (1 - self.Y_su) * self.f_h2_su,
            ("R5", "Liq", "S_ch4"): 0,
            ("R5", "Liq", "S_IC"): (
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
            ("R5", "Liq", "S_IN"): (-self.Y_su * self.N_bac) * mw_n,
            ("R5", "Liq", "S_I"): 0,
            ("R5", "Liq", "X_c"): 0,
            ("R5", "Liq", "X_ch"): 0,
            ("R5", "Liq", "X_pr"): 0,
            ("R5", "Liq", "X_li"): 0,
            ("R5", "Liq", "X_su"): self.Y_su,
            ("R5", "Liq", "X_aa"): 0,
            ("R5", "Liq", "X_fa"): 0,
            ("R5", "Liq", "X_c4"): 0,
            ("R5", "Liq", "X_pro"): 0,
            ("R5", "Liq", "X_ac"): 0,
            ("R5", "Liq", "X_h2"): 0,
            ("R5", "Liq", "X_I"): 0,
            # R6:  Uptake of amino acids
            ("R6", "Liq", "H2O"): 0,
            ("R6", "Liq", "S_su"): 0,
            ("R6", "Liq", "S_aa"): -1,
            ("R6", "Liq", "S_fa"): 0,
            ("R6", "Liq", "S_va"): (1 - self.Y_aa) * self.f_va_aa,
            ("R6", "Liq", "S_bu"): (1 - self.Y_aa) * self.f_bu_aa,
            ("R6", "Liq", "S_pro"): (1 - self.Y_aa) * self.f_pro_aa,
            ("R6", "Liq", "S_ac"): (1 - self.Y_aa) * self.f_ac_aa,
            ("R6", "Liq", "S_h2"): (1 - self.Y_aa) * self.f_h2_aa,
            ("R6", "Liq", "S_ch4"): 0,
            ("R6", "Liq", "S_IC"): (
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
            ("R6", "Liq", "S_IN"): (self.N_aa - self.Y_aa * self.N_bac) * mw_n,
            ("R6", "Liq", "S_I"): 0,
            ("R6", "Liq", "X_c"): 0,
            ("R6", "Liq", "X_ch"): 0,
            ("R6", "Liq", "X_pr"): 0,
            ("R6", "Liq", "X_li"): 0,
            ("R6", "Liq", "X_su"): 0,
            ("R6", "Liq", "X_aa"): self.Y_aa,
            ("R6", "Liq", "X_fa"): 0,
            ("R6", "Liq", "X_c4"): 0,
            ("R6", "Liq", "X_pro"): 0,
            ("R6", "Liq", "X_ac"): 0,
            ("R6", "Liq", "X_h2"): 0,
            ("R6", "Liq", "X_I"): 0,
            # R7:  Uptake of long chain fatty acids (LCFAs)
            ("R7", "Liq", "H2O"): 0,
            ("R7", "Liq", "S_su"): 0,
            ("R7", "Liq", "S_aa"): 0,
            ("R7", "Liq", "S_fa"): -1,
            ("R7", "Liq", "S_va"): 0,
            ("R7", "Liq", "S_bu"): 0,
            ("R7", "Liq", "S_pro"): 0,
            ("R7", "Liq", "S_ac"): (1 - self.Y_fa) * 0.7,
            ("R7", "Liq", "S_h2"): (1 - self.Y_fa) * 0.3,
            ("R7", "Liq", "S_ch4"): 0,
            ("R7", "Liq", "S_IC"): (
                self.Ci["S_fa"]
                - (1 - self.Y_fa) * 0.7 * self.Ci["S_ac"]
                - self.Y_fa * self.Ci["X_fa"]
            )
            * mw_c,
            ("R7", "Liq", "S_IN"): (-self.Y_fa * self.N_bac) * mw_n,
            ("R7", "Liq", "S_I"): 0,
            ("R7", "Liq", "X_c"): 0,
            ("R7", "Liq", "X_ch"): 0,
            ("R7", "Liq", "X_pr"): 0,
            ("R7", "Liq", "X_li"): 0,
            ("R7", "Liq", "X_su"): 0,
            ("R7", "Liq", "X_aa"): 0,
            ("R7", "Liq", "X_fa"): self.Y_fa,
            ("R7", "Liq", "X_c4"): 0,
            ("R7", "Liq", "X_pro"): 0,
            ("R7", "Liq", "X_ac"): 0,
            ("R7", "Liq", "X_h2"): 0,
            ("R7", "Liq", "X_I"): 0,
            # R8:  Uptake of valerate
            ("R8", "Liq", "H2O"): 0,
            ("R8", "Liq", "S_su"): 0,
            ("R8", "Liq", "S_aa"): 0,
            ("R8", "Liq", "S_fa"): 0,
            ("R8", "Liq", "S_va"): -1,
            ("R8", "Liq", "S_bu"): 0,
            ("R8", "Liq", "S_pro"): (1 - self.Y_c4) * 0.54,
            ("R8", "Liq", "S_ac"): (1 - self.Y_c4) * 0.31,
            ("R8", "Liq", "S_h2"): (1 - self.Y_c4) * 0.15,
            ("R8", "Liq", "S_ch4"): 0,
            ("R8", "Liq", "S_IC"): (
                self.Ci["S_va"]
                - (1 - self.Y_c4) * 0.54 * self.Ci["S_pro"]
                - (1 - self.Y_c4) * 0.31 * self.Ci["S_ac"]
                - self.Y_c4 * self.Ci["X_c4"]
            )
            * mw_c,
            ("R8", "Liq", "S_IN"): (-self.Y_c4 * self.N_bac) * mw_n,
            ("R8", "Liq", "S_I"): 0,
            ("R8", "Liq", "X_c"): 0,
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
            # R9:  Uptake of butyrate
            ("R9", "Liq", "H2O"): 0,
            ("R9", "Liq", "S_su"): 0,
            ("R9", "Liq", "S_aa"): 0,
            ("R9", "Liq", "S_fa"): 0,
            ("R9", "Liq", "S_va"): 0,
            ("R9", "Liq", "S_bu"): -1,
            ("R9", "Liq", "S_pro"): 0,
            ("R9", "Liq", "S_ac"): (1 - self.Y_c4) * 0.8,
            ("R9", "Liq", "S_h2"): (1 - self.Y_c4) * 0.2,
            ("R9", "Liq", "S_ch4"): 0,
            ("R9", "Liq", "S_IC"): (
                self.Ci["S_bu"]
                - (1 - self.Y_c4) * 0.8 * self.Ci["S_ac"]
                - self.Y_c4 * self.Ci["X_c4"]
            )
            * mw_c,
            ("R9", "Liq", "S_IN"): (-self.Y_c4 * self.N_bac) * mw_n,
            ("R9", "Liq", "S_I"): 0,
            ("R9", "Liq", "X_c"): 0,
            ("R9", "Liq", "X_ch"): 0,
            ("R9", "Liq", "X_pr"): 0,
            ("R9", "Liq", "X_li"): 0,
            ("R9", "Liq", "X_su"): 0,
            ("R9", "Liq", "X_aa"): 0,
            ("R9", "Liq", "X_fa"): 0,
            ("R9", "Liq", "X_c4"): self.Y_c4,
            ("R9", "Liq", "X_pro"): 0,
            ("R9", "Liq", "X_ac"): 0,
            ("R9", "Liq", "X_h2"): 0,
            ("R9", "Liq", "X_I"): 0,
            # R10: Uptake of propionate
            ("R10", "Liq", "H2O"): 0,
            ("R10", "Liq", "S_su"): 0,
            ("R10", "Liq", "S_aa"): 0,
            ("R10", "Liq", "S_fa"): 0,
            ("R10", "Liq", "S_va"): 0,
            ("R10", "Liq", "S_bu"): 0,
            ("R10", "Liq", "S_pro"): -1,
            ("R10", "Liq", "S_ac"): (1 - self.Y_pro) * 0.57,
            ("R10", "Liq", "S_h2"): (1 - self.Y_pro) * 0.43,
            ("R10", "Liq", "S_ch4"): 0,
            ("R10", "Liq", "S_IC"): (
                self.Ci["S_pro"]
                - (1 - self.Y_pro) * 0.57 * self.Ci["S_ac"]
                - self.Y_pro * self.Ci["X_pro"]
            )
            * mw_c,
            ("R10", "Liq", "S_IN"): (-self.Y_pro * self.N_bac) * mw_n,
            ("R10", "Liq", "S_I"): 0,
            ("R10", "Liq", "X_c"): 0,
            ("R10", "Liq", "X_ch"): 0,
            ("R10", "Liq", "X_pr"): 0,
            ("R10", "Liq", "X_li"): 0,
            ("R10", "Liq", "X_su"): 0,
            ("R10", "Liq", "X_aa"): 0,
            ("R10", "Liq", "X_fa"): 0,
            ("R10", "Liq", "X_c4"): 0,
            ("R10", "Liq", "X_pro"): self.Y_pro,
            ("R10", "Liq", "X_ac"): 0,
            ("R10", "Liq", "X_h2"): 0,
            ("R10", "Liq", "X_I"): 0,
            # R11: Uptake of acetate
            ("R11", "Liq", "H2O"): 0,
            ("R11", "Liq", "S_su"): 0,
            ("R11", "Liq", "S_aa"): 0,
            ("R11", "Liq", "S_fa"): 0,
            ("R11", "Liq", "S_va"): 0,
            ("R11", "Liq", "S_bu"): 0,
            ("R11", "Liq", "S_pro"): 0,
            ("R11", "Liq", "S_ac"): -1,
            ("R11", "Liq", "S_h2"): 0,
            ("R11", "Liq", "S_ch4"): 1 - self.Y_ac,
            ("R11", "Liq", "S_IC"): (
                self.Ci["S_ac"]
                - (1 - self.Y_ac) * self.Ci["S_ch4"]
                - self.Y_ac * self.Ci["X_ac"]
            )
            * mw_c,
            ("R11", "Liq", "S_IN"): (-self.Y_ac * self.N_bac) * mw_n,
            ("R11", "Liq", "S_I"): 0,
            ("R11", "Liq", "X_c"): 0,
            ("R11", "Liq", "X_ch"): 0,
            ("R11", "Liq", "X_pr"): 0,
            ("R11", "Liq", "X_li"): 0,
            ("R11", "Liq", "X_su"): 0,
            ("R11", "Liq", "X_aa"): 0,
            ("R11", "Liq", "X_fa"): 0,
            ("R11", "Liq", "X_c4"): 0,
            ("R11", "Liq", "X_pro"): 0,
            ("R11", "Liq", "X_ac"): self.Y_ac,
            ("R11", "Liq", "X_h2"): 0,
            ("R11", "Liq", "X_I"): 0,
            # R12: Uptake of hydrogen
            ("R12", "Liq", "H2O"): 0,
            ("R12", "Liq", "S_su"): 0,
            ("R12", "Liq", "S_aa"): 0,
            ("R12", "Liq", "S_fa"): 0,
            ("R12", "Liq", "S_va"): 0,
            ("R12", "Liq", "S_bu"): 0,
            ("R12", "Liq", "S_pro"): 0,
            ("R12", "Liq", "S_ac"): 0,
            ("R12", "Liq", "S_h2"): -1,
            ("R12", "Liq", "S_ch4"): (1 - self.Y_h2),
            ("R12", "Liq", "S_IC"): (
                -(1 - self.Y_h2) * self.Ci["S_ch4"] - self.Y_h2 * self.Ci["X_h2"]
            )
            * mw_c,
            ("R12", "Liq", "S_IN"): (-self.Y_h2 * self.N_bac) * mw_n,
            ("R12", "Liq", "S_I"): 0,
            ("R12", "Liq", "X_c"): 0,
            ("R12", "Liq", "X_ch"): 0,
            ("R12", "Liq", "X_pr"): 0,
            ("R12", "Liq", "X_li"): 0,
            ("R12", "Liq", "X_su"): 0,
            ("R12", "Liq", "X_aa"): 0,
            ("R12", "Liq", "X_fa"): 0,
            ("R12", "Liq", "X_c4"): 0,
            ("R12", "Liq", "X_pro"): 0,
            ("R12", "Liq", "X_ac"): 0,
            ("R12", "Liq", "X_h2"): self.Y_h2,
            ("R12", "Liq", "X_I"): 0,
            # R13: Decay of X_su
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
            ("R13", "Liq", "S_IC"): (-self.Ci["X_c"] + self.Ci["X_ac"]) * mw_c,
            ("R13", "Liq", "S_IN"): (self.N_bac - self.N_xc) * mw_n,
            ("R13", "Liq", "S_I"): 0,
            ("R13", "Liq", "X_c"): 1,
            ("R13", "Liq", "X_ch"): 0,
            ("R13", "Liq", "X_pr"): 0,
            ("R13", "Liq", "X_li"): 0,
            ("R13", "Liq", "X_su"): -1,
            ("R13", "Liq", "X_aa"): 0,
            ("R13", "Liq", "X_fa"): 0,
            ("R13", "Liq", "X_c4"): 0,
            ("R13", "Liq", "X_pro"): 0,
            ("R13", "Liq", "X_ac"): 0,
            ("R13", "Liq", "X_h2"): 0,
            ("R13", "Liq", "X_I"): 0,
            # R14: Decay of X_aa
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
            ("R14", "Liq", "S_IC"): (-self.Ci["X_c"] + self.Ci["X_ac"]) * mw_c,
            ("R14", "Liq", "S_IN"): (self.N_bac - self.N_xc) * mw_n,
            ("R14", "Liq", "S_I"): 0,
            ("R14", "Liq", "X_c"): 1,
            ("R14", "Liq", "X_ch"): 0,
            ("R14", "Liq", "X_pr"): 0,
            ("R14", "Liq", "X_li"): 0,
            ("R14", "Liq", "X_su"): 0,
            ("R14", "Liq", "X_aa"): -1,
            ("R14", "Liq", "X_fa"): 0,
            ("R14", "Liq", "X_c4"): 0,
            ("R14", "Liq", "X_pro"): 0,
            ("R14", "Liq", "X_ac"): 0,
            ("R14", "Liq", "X_h2"): 0,
            ("R14", "Liq", "X_I"): 0,
            # R15: Decay of X_fa
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
            ("R15", "Liq", "S_IC"): (-self.Ci["X_c"] + self.Ci["X_ac"]) * mw_c,
            ("R15", "Liq", "S_IN"): (self.N_bac - self.N_xc) * mw_n,
            ("R15", "Liq", "S_I"): 0,
            ("R15", "Liq", "X_c"): 1,
            ("R15", "Liq", "X_ch"): 0,
            ("R15", "Liq", "X_pr"): 0,
            ("R15", "Liq", "X_li"): 0,
            ("R15", "Liq", "X_su"): 0,
            ("R15", "Liq", "X_aa"): 0,
            ("R15", "Liq", "X_fa"): -1,
            ("R15", "Liq", "X_c4"): 0,
            ("R15", "Liq", "X_pro"): 0,
            ("R15", "Liq", "X_ac"): 0,
            ("R15", "Liq", "X_h2"): 0,
            ("R15", "Liq", "X_I"): 0,
            # R16: Decay of X_c4
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
            ("R16", "Liq", "S_IC"): (-self.Ci["X_c"] + self.Ci["X_ac"]) * mw_c,
            ("R16", "Liq", "S_IN"): (self.N_bac - self.N_xc) * mw_n,
            ("R16", "Liq", "S_I"): 0,
            ("R16", "Liq", "X_c"): 1,
            ("R16", "Liq", "X_ch"): 0,
            ("R16", "Liq", "X_pr"): 0,
            ("R16", "Liq", "X_li"): 0,
            ("R16", "Liq", "X_su"): 0,
            ("R16", "Liq", "X_aa"): 0,
            ("R16", "Liq", "X_fa"): 0,
            ("R16", "Liq", "X_c4"): -1,
            ("R16", "Liq", "X_pro"): 0,
            ("R16", "Liq", "X_ac"): 0,
            ("R16", "Liq", "X_h2"): 0,
            ("R16", "Liq", "X_I"): 0,
            # R17: Decay of X_pro
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
            ("R17", "Liq", "S_IC"): (-self.Ci["X_c"] + self.Ci["X_ac"]) * mw_c,
            ("R17", "Liq", "S_IN"): (self.N_bac - self.N_xc) * mw_n,
            ("R17", "Liq", "S_I"): 0,
            ("R17", "Liq", "X_c"): 1,
            ("R17", "Liq", "X_ch"): 0,
            ("R17", "Liq", "X_pr"): 0,
            ("R17", "Liq", "X_li"): 0,
            ("R17", "Liq", "X_su"): 0,
            ("R17", "Liq", "X_aa"): 0,
            ("R17", "Liq", "X_fa"): 0,
            ("R17", "Liq", "X_c4"): 0,
            ("R17", "Liq", "X_pro"): -1,
            ("R17", "Liq", "X_ac"): 0,
            ("R17", "Liq", "X_h2"): 0,
            ("R17", "Liq", "X_I"): 0,
            # R18: Decay of X_ac
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
            ("R18", "Liq", "S_IC"): (-self.Ci["X_c"] + self.Ci["X_ac"]) * mw_c,
            ("R18", "Liq", "S_IN"): (self.N_bac - self.N_xc) * mw_n,
            ("R18", "Liq", "S_I"): 0,
            ("R18", "Liq", "X_c"): 1,
            ("R18", "Liq", "X_ch"): 0,
            ("R18", "Liq", "X_pr"): 0,
            ("R18", "Liq", "X_li"): 0,
            ("R18", "Liq", "X_su"): 0,
            ("R18", "Liq", "X_aa"): 0,
            ("R18", "Liq", "X_fa"): 0,
            ("R18", "Liq", "X_c4"): 0,
            ("R18", "Liq", "X_pro"): 0,
            ("R18", "Liq", "X_ac"): -1,
            ("R18", "Liq", "X_h2"): 0,
            ("R18", "Liq", "X_I"): 0,
            # R19: Decay of X_h2
            ("R19", "Liq", "H2O"): 0,
            ("R19", "Liq", "S_su"): 0,
            ("R19", "Liq", "S_aa"): 0,
            ("R19", "Liq", "S_fa"): 0,
            ("R19", "Liq", "S_va"): 0,
            ("R19", "Liq", "S_bu"): 0,
            ("R19", "Liq", "S_pro"): 0,
            ("R19", "Liq", "S_ac"): 0,
            ("R19", "Liq", "S_h2"): 0,
            ("R19", "Liq", "S_ch4"): 0,
            ("R19", "Liq", "S_IC"): (-self.Ci["X_c"] + self.Ci["X_ac"]) * mw_c,
            ("R19", "Liq", "S_IN"): (self.N_bac - self.N_xc) * mw_n,
            ("R19", "Liq", "S_I"): 0,
            ("R19", "Liq", "X_c"): 1,
            ("R19", "Liq", "X_ch"): 0,
            ("R19", "Liq", "X_pr"): 0,
            ("R19", "Liq", "X_li"): 0,
            ("R19", "Liq", "X_su"): 0,
            ("R19", "Liq", "X_aa"): 0,
            ("R19", "Liq", "X_fa"): 0,
            ("R19", "Liq", "X_c4"): 0,
            ("R19", "Liq", "X_pro"): 0,
            ("R19", "Liq", "X_ac"): 0,
            ("R19", "Liq", "X_h2"): -1,
            ("R19", "Liq", "X_I"): 0,
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


class ADM1ReactionScaler(CustomScalerBase):
    """
    Scaler for the Anaerobic Digestion Model No.1 reaction package.

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


class _ADM1ReactionBlock(ReactionBlockBase):
    """
    This Class contains methods which should be applied to Reaction Blocks as a
    whole, rather than individual elements of indexed Reaction Blocks.
    """

    default_scaler = ADM1ReactionScaler

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


@declare_process_block_class("ADM1ReactionBlock", block_class=_ADM1ReactionBlock)
class ADM1ReactionBlockData(ReactionBlockDataBase):
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

        # Initial values of rates of reaction [2]
        self.rates = {
            "R1": 1.786e-06,
            "R2": 3.235e-06,
            "R3": 1.187e-05,
            "R4": 3.412e-06,
            "R5": 3.404e-06,
            "R6": 1.187e-05,
            "R7": 3.185e-06,
            "R8": 2.505e-06,
            "R9": 3.230e-06,
            "R10": 2.636e-06,
            "R11": 1.220e-05,
            "R12": 4.184e-06,
            "R13": 9.736e-08,
            "R14": 2.730e-07,
            "R15": 5.627e-08,
            "R16": 9.999e-08,
            "R17": 3.179e-08,
            "R18": 1.761e-07,
            "R19": 7.339e-08,
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
            initialize=0.011,
            domain=pyo.NonNegativeReals,
            doc="mass concentration of va-",
            units=pyo.units.kg / pyo.units.m**3,
        )
        self.conc_mass_bu = pyo.Var(
            initialize=0.013,
            domain=pyo.NonNegativeReals,
            doc="mass concentration of bu-",
            units=pyo.units.kg / pyo.units.m**3,
        )
        self.conc_mass_pro = pyo.Var(
            initialize=0.016,
            domain=pyo.NonNegativeReals,
            doc="mass concentration of pro-",
            units=pyo.units.kg / pyo.units.m**3,
        )
        self.conc_mass_ac = pyo.Var(
            initialize=0.2,
            domain=pyo.NonNegativeReals,
            doc="mass concentration of ac-",
            units=pyo.units.kg / pyo.units.m**3,
        )
        self.conc_mol_hco3 = pyo.Var(
            initialize=0.14,
            domain=pyo.NonNegativeReals,
            doc="molar concentration of hco3",
            units=pyo.units.kmol / pyo.units.m**3,
        )
        self.conc_mol_nh3 = pyo.Var(
            initialize=0.0041,
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
        self.S_H = pyo.Var(
            initialize=1e-7,
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

        # Equation from [2]
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

        # Equation from [2]
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

        # Equation from [2]
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

        def S_H_rule(self):
            return (
                self.state_ref.cations
                + self.conc_mol_nh4
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
            if r == "R5" or r == "R6":
                return self.I[r] == pyo.exp(self.I_pH_aa) * self.I_IN_lim
            elif r == "R7":
                return self.I[r] == pyo.exp(self.I_pH_aa) * self.I_IN_lim * self.I_h2_fa
            elif r == "R8" or r == "R9":
                return self.I[r] == pyo.exp(self.I_pH_aa) * self.I_IN_lim * self.I_h2_c4
            elif r == "R10":
                return (
                    self.I[r] == pyo.exp(self.I_pH_aa) * self.I_IN_lim * self.I_h2_pro
                )
            elif r == "R11":
                return self.I[r] == pyo.exp(self.I_pH_ac) * self.I_IN_lim * self.I_nh3
            elif r == "R12":
                return self.I[r] == pyo.exp(self.I_pH_h2) * self.I_IN_lim
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
                    # R1:  Disintegration
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_dis * b.conc_mass_comp_ref["X_c"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R2":
                    # R2: Hydrolysis of carbohydrates
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_hyd_ch * b.conc_mass_comp_ref["X_ch"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R3":
                    # R3: Hydrolysis of proteins
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_hyd_pr * b.conc_mass_comp_ref["X_pr"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R4":
                    # R4: Hydrolysis of lipids
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_hyd_li * b.conc_mass_comp_ref["X_li"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R5":
                    # R5: Uptake of sugars
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_m_su
                        * b.conc_mass_comp_ref["S_su"]
                        / (b.params.K_S_su + b.conc_mass_comp_ref["S_su"])
                        * b.conc_mass_comp_ref["X_su"]
                        * b.I[r],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R6":
                    # R6: Uptake of amino acids
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_m_aa
                        * b.conc_mass_comp_ref["S_aa"]
                        / (b.params.K_S_aa + b.conc_mass_comp_ref["S_aa"])
                        * b.conc_mass_comp_ref["X_aa"]
                        * b.I[r],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R7":
                    # R7: Uptake of long chain fatty acids (LCFAs)
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_m_fa
                        * b.conc_mass_comp_ref["S_fa"]
                        / (b.params.K_S_fa + b.conc_mass_comp_ref["S_fa"])
                        * b.conc_mass_comp_ref["X_fa"]
                        * b.I[r],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R8":
                    # R8: Uptake of valerate
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
                        * b.I[r],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R9":
                    # R9:  Uptake of butyrate
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
                        * b.I[r],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R10":
                    # R10: Uptake of propionate
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_m_pro
                        * b.conc_mass_comp_ref["S_pro"]
                        / (b.params.K_S_pro + b.conc_mass_comp_ref["S_pro"])
                        * b.conc_mass_comp_ref["X_pro"]
                        * b.I[r],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R11":
                    # R11: Uptake of acetate
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_m_ac
                        * b.conc_mass_comp_ref["S_ac"]
                        / (b.params.K_S_ac + b.conc_mass_comp_ref["S_ac"])
                        * b.conc_mass_comp_ref["X_ac"]
                        * b.I[r],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R12":
                    # R12: Uptake of hydrogen
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_m_h2
                        * b.conc_mass_comp_ref["S_h2"]
                        / (b.params.K_S_h2 + b.conc_mass_comp_ref["S_h2"])
                        * b.conc_mass_comp_ref["X_h2"]
                        * b.I[r],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R13":
                    # R13: Decay of X_su
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_dec_X_su * b.conc_mass_comp_ref["X_su"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R14":
                    # R14: Decay of X_aa
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_dec_X_aa * b.conc_mass_comp_ref["X_aa"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R15":
                    # R15: Decay of X_fa
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_dec_X_fa * b.conc_mass_comp_ref["X_fa"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R16":
                    # R16: Decay of X_c4
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_dec_X_c4 * b.conc_mass_comp_ref["X_c4"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R17":
                    # R17: Decay of X_pro
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_dec_X_pro * b.conc_mass_comp_ref["X_pro"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R18":
                    # R18: Decay of X_ac
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_dec_X_ac * b.conc_mass_comp_ref["X_ac"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R19":
                    # R19: Decay of X_h2
                    return b.reaction_rate[r] == (
                        pyo.units.convert(b.params.k_dec_X_h2, to_units=1 / pyo.units.s)
                        * b.conc_mass_comp_ref["X_h2"]
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

        iscale.set_scaling_factor(self.I, 1e1)
        iscale.set_scaling_factor(self.conc_mass_va, 1e2)
        iscale.set_scaling_factor(self.conc_mass_bu, 1e2)
        iscale.set_scaling_factor(self.conc_mass_pro, 1e2)
        iscale.set_scaling_factor(self.conc_mass_ac, 1e1)
        iscale.set_scaling_factor(self.conc_mol_hco3, 1e1)
        iscale.set_scaling_factor(self.conc_mol_nh3, 1e3)
        iscale.set_scaling_factor(self.conc_mol_co2, 1e3)
        iscale.set_scaling_factor(self.conc_mol_nh4, 1e1)
        iscale.set_scaling_factor(self.S_H, 1e5)
        iscale.set_scaling_factor(self.pKW, 1e0)
        iscale.set_scaling_factor(self.pK_a_co2, 1e0)
        iscale.set_scaling_factor(self.pK_a_IN, 1e0)
        iscale.set_scaling_factor(self.pH, 1e0)

        for i, c in self.rates.items():
            iscale.set_scaling_factor(self.reaction_rate[i], 1 / c)

        for i, c in self.rate_expression.items():
            iscale.constraint_scaling_transform(
                c,
                iscale.get_scaling_factor(
                    self.reaction_rate[i], default=1, warning=True
                ),
                overwrite=True,
            )
