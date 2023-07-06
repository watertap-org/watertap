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
Translator block representing the ASM2d/ADM1 interface.

Assumptions:
     * Steady-state only

Model formulated from:

Flores-Alsina, X., Solon, K., Mbamba, C.K., Tait, S., Gernaey, K.V., Jeppsson, U. and Batstone, D.J., 2016.
Modelling phosphorus (P), sulfur (S) and iron (Fe) interactions for dynamic simulations of anaerobic digestion processes.
Water Research, 95, pp.370-382. https://github.com/wwtmodels/Plant-Wide-Models
"""

# Import Pyomo libraries
from pyomo.common.config import ConfigBlock, ConfigValue

# Import IDAES cores
from idaes.core import declare_process_block_class
from idaes.models.unit_models.translator import TranslatorData
from idaes.core.util.config import (
    is_reaction_parameter_block,
)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.solvers import get_solver
import idaes.logger as idaeslog

from pyomo.environ import (
    Param,
    PositiveReals,
    Var,
    sqrt,
    units as pyunits,
    check_optimal_termination,
    Set,
    Expr_if,
)

from idaes.core.util.exceptions import InitializationError

from watertap.property_models.activated_sludge.modified_asm2d_reactions import (
    ModifiedASM2dReactionParameterBlock as ASM2d,
)
from watertap.property_models.activated_sludge.modified_asm2d_reactions import (
    DecaySwitch,
)
from watertap.property_models.anaerobic_digestion.modified_adm1_reactions import (
    ModifiedADM1ReactionParameterBlock as ADM1,
)

__author__ = "Marcus Holly"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("Translator_ASM2d_ADM1")
class TranslatorDataASM2dADM1(TranslatorData):
    """
    Translator block representing the ASM2d/ADM1 interface
    """

    CONFIG = TranslatorData.CONFIG()
    CONFIG.declare(
        "inlet_reaction_package",
        ConfigValue(
            default=None,
            domain=is_reaction_parameter_block,
            description="Reaction package to use for control volume",
            doc="""Reaction parameter object used to define reaction calculations,
**default** - None.
**Valid values:** {
**None** - no reaction package,
**ReactionParameterBlock** - a ReactionParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "inlet_reaction_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing reaction packages",
            doc="""A ConfigBlock with arguments to be passed to a reaction block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see reaction package for documentation.}""",
        ),
    )
    CONFIG.declare(
        "outlet_reaction_package",
        ConfigValue(
            default=None,
            domain=is_reaction_parameter_block,
            description="Reaction package to use for control volume",
            doc="""Reaction parameter object used to define reaction calculations,
    **default** - None.
    **Valid values:** {
    **None** - no reaction package,
    **ReactionParameterBlock** - a ReactionParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "outlet_reaction_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing reaction packages",
            doc="""A ConfigBlock with arguments to be passed to a reaction block(s)
    and used when constructing these,
    **default** - None.
    **Valid values:** {
    see reaction package for documentation.}""",
        ),
    )

    def build(self):
        """
        Begin building model.
        Args:
            None
        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super(TranslatorDataASM2dADM1, self).build()

        mw_p = 31 * pyunits.kg / pyunits.kmol
        mw_n = 14 * pyunits.kg / pyunits.kmol
        mw_c = 12 * pyunits.kg / pyunits.kmol
        mw_k = 39 * pyunits.kg / pyunits.kmol
        mw_mg = 24.3 * pyunits.kg / pyunits.kmol

        @self.Constraint(
            self.flowsheet().time,
            doc="Equality volumetric flow equation",
        )
        def eq_flow_vol_rule(blk, t):
            return blk.properties_out[t].flow_vol == blk.properties_in[t].flow_vol

        @self.Constraint(
            self.flowsheet().time,
            doc="Equality temperature equation",
        )
        def eq_temperature_rule(blk, t):
            return blk.properties_out[t].temperature == blk.properties_in[t].temperature

        @self.Constraint(
            self.flowsheet().time,
            doc="Equality pressure equation",
        )
        def eq_pressure_rule(blk, t):
            return blk.properties_out[t].pressure == blk.properties_in[t].pressure

        # -------------------------------------------Step 1----------------------------------------------------------------#
        self.Y_H = Var(
            initialize=0.625,
            units=pyunits.dimensionless,
            domain=PositiveReals,
            doc="Yield coefficient for heterotrophic biomass, [kg COD/kg COD]",
        )

        self.COD_SO2 = Var(
            self.flowsheet().time,
            initialize=self.properties_in[0].conc_mass_comp["S_O2"]
            / ((1 - self.Y_H) / self.Y_H),
            units=pyunits.kg / pyunits.m**3,
            domain=PositiveReals,
            doc="Concentration of COD demanding compounds in S_O2",
        )

        @self.Expression(self.flowsheet().time, doc="S_O2 concentration at step 1")
        def SO2_AS1(blk, t):
            return (
                blk.properties_in[t].conc_mass_comp["S_O2"]
                - (1 - blk.config.inlet_reaction_package.Y_H)
                / blk.config.inlet_reaction_package.Y_H
                * self.COD_SO2[t]
            )

        @self.Expression(self.flowsheet().time, doc="S_A concentration at step 1")
        def SA_AS1(blk, t):
            return (
                blk.properties_in[t].conc_mass_comp["S_A"]
                - (1 / blk.config.inlet_reaction_package.Y_H) * self.COD_SO2[t]
            )

        @self.Expression(self.flowsheet().time, doc="S_NH4 concentration at step 1")
        def SNH4_AS1(blk, t):
            return (
                blk.properties_in[t].conc_mass_comp["S_NH4"]
                - blk.config.inlet_reaction_package.i_NBM * self.COD_SO2[t]
            )

        @self.Expression(self.flowsheet().time, doc="S_PO4 concentration at step 1")
        def SPO4_AS1(blk, t):
            return (
                blk.properties_in[t].conc_mass_comp["S_PO4"]
                - blk.config.inlet_reaction_package.i_PBM * self.COD_SO2[t]
            )

        @self.Expression(self.flowsheet().time, doc="S_IC concentration at step 1")
        def SIC_AS1(blk, t):
            return (
                blk.properties_in[t].conc_mass_comp["S_IC"]
                - (-1 / blk.config.inlet_reaction_package.Y_H)
                * self.COD_SO2[t]
                * blk.config.inlet_reaction_package.i_CSA
                + self.COD_SO2[t] * blk.config.inlet_reaction_package.i_CXB
            )

        @self.Expression(self.flowsheet().time, doc="X_H concentration at step 1")
        def XH_AS1(blk, t):
            return blk.properties_in[t].conc_mass_comp["X_H"] + self.COD_SO2[t]

        # -------------------------------------------Step 2----------------------------------------------------------------#
        self.i_NOx_N2 = Var(
            initialize=2.8571,
            units=pyunits.dimensionless,
            domain=PositiveReals,
            doc="Nitrogen oxide coefficient for N2",
        )
        self.COD_SNO3 = Var(
            self.flowsheet().time,
            initialize=self.properties_in[0].conc_mass_comp["S_NO3"]
            / ((1 - self.Y_H) / (self.i_NOx_N2 * self.Y_H)),
            units=pyunits.kg / pyunits.m**3,
            domain=PositiveReals,
            doc="Concentration of COD demanding compounds in S_NO3",
        )

        @self.Expression(self.flowsheet().time, doc="S_A concentration at step 2")
        def SA_AS2(blk, t):
            return (
                blk.SA_AS1[t]
                - (1 / blk.config.inlet_reaction_package.Y_H) * self.COD_SNO3[t]
            )

        @self.Expression(self.flowsheet().time, doc="S_NH4 concentration at step 2")
        def SNH4_AS2(blk, t):
            return (
                blk.SNH4_AS1[t]
                - self.COD_SNO3[t] * blk.config.inlet_reaction_package.i_NBM
            )

        @self.Expression(self.flowsheet().time, doc="S_N2 concentration at step 2")
        def SN2_AS2(blk, t):
            return (
                blk.properties_in[t].conc_mass_comp["S_N2"]
                + (
                    (1 - blk.config.inlet_reaction_package.Y_H)
                    / (
                        blk.config.inlet_reaction_package.i_NOx_N2
                        * blk.config.inlet_reaction_package.Y_H
                    )
                )
                * self.COD_SNO3[t]
            )

        @self.Expression(self.flowsheet().time, doc="S_NO3 concentration at step 2")
        def SNO3_AS2(blk, t):
            return (
                blk.properties_in[t].conc_mass_comp["S_NO3"]
                - (
                    (1 - blk.config.inlet_reaction_package.Y_H)
                    / (
                        blk.config.inlet_reaction_package.i_NOx_N2
                        * blk.config.inlet_reaction_package.Y_H
                    )
                )
                * self.COD_SNO3[t]
            )

        @self.Expression(self.flowsheet().time, doc="S_PO4 concentration at step 2")
        def SPO4_AS2(blk, t):
            return (
                blk.SPO4_AS1[t]
                - self.COD_SNO3[t] * blk.config.inlet_reaction_package.i_PBM
            )

        @self.Expression(self.flowsheet().time, doc="S_IC concentration at step 2")
        def SIC_AS2(blk, t):
            return (
                blk.SIC_AS1[t]
                - (-1 / blk.config.inlet_reaction_package.Y_H)
                * self.COD_SO2[t]
                * blk.config.inlet_reaction_package.i_CSA
                + self.COD_SNO3[t] * blk.config.inlet_reaction_package.i_CXB
            )

        @self.Expression(self.flowsheet().time, doc="X_H concentration at step 2")
        def XH_AS2(blk, t):
            return blk.XH_AS1[t] + self.COD_SNO3[t]

        # -------------------------------------------Step 3----------------------------------------------------------------#
        self.i_NSF = Var(
            initialize=0.03352,
            units=pyunits.dimensionless,
            domain=PositiveReals,
            doc="N content of fermentable substrate, S_F, [kg N/kg COD]",
        )
        self.i_PSF = Var(
            initialize=0.00559,
            units=pyunits.dimensionless,
            domain=PositiveReals,
            doc="P content of fermentable substrate, S_F, [kg P/kg COD]",
        )
        self.N_aa = Var(
            initialize=0.0079034,
            units=pyunits.kmol / pyunits.kg,
            domain=PositiveReals,
            doc="Nitrogen content of S_aa",
        )
        self.S_ND = Var(
            initialize=self.properties_in[0].conc_mass_comp["S_F"] * self.i_NSF,
            units=pyunits.kg / pyunits.m**3,
            domain=PositiveReals,
            doc="Nitrogen demand for soluble inerts",
        )
        self.S_PD = Var(
            initialize=self.properties_in[0].conc_mass_comp["S_F"] * self.i_PSF,
            units=pyunits.kg / pyunits.m**3,
            domain=PositiveReals,
            doc="Phosphorus demand for soluble inerts",
        )
        self.SN_org = Var(
            self.flowsheet().time,
            initialize=self.S_ND / (self.N_aa * mw_n),
            units=pyunits.kg / pyunits.m**3,
            domain=PositiveReals,
            doc="Organic nitrogen from soluble inerts",
        )

        @self.Constraint(
            self.flowsheet().time,
            doc="Nitrogen demand for soluble inerts divided by nitrogen content in amino acids",
        )
        def ReqCOD_S(blk, t):
            return blk.SN_org[t] == self.S_ND / (self.N_aa * mw_n)

        # TODO: Double check if 1000 should be here and throughout the rest (just them going from g/m3 to kg/m3?)
        @self.Expression(self.flowsheet().time, doc="Monosaccharides mapping")
        def Ssu_mapping(blk, t):
            return Expr_if(
                blk.SN_org[t] >= self.properties_in[0].conc_mass_comp["S_F"],
                1e-9 * pyunits.kg / pyunits.m**3,
                (self.properties_in[0].conc_mass_comp["S_F"] - self.SN_org[t]) / 1000,
            )

        @self.Expression(self.flowsheet().time, doc="Amino acids mapping")
        def Saa_mapping(blk, t):
            return Expr_if(
                blk.SN_org[t] >= self.properties_in[0].conc_mass_comp["S_F"],
                self.properties_in[0].conc_mass_comp["S_F"] / 1000,
                self.SN_org[t] / 1000,
            )

        @self.Constraint(
            self.flowsheet().time,
            doc="Monosaccharides output",
        )
        def Ssu_output(blk, t):
            return blk.properties_out[t].conc_mass_comp["S_su"] == blk.Ssu_mapping[t]

        @self.Constraint(
            self.flowsheet().time,
            doc="Amino acids output",
        )
        def Saa_output(blk, t):
            return blk.properties_out[t].conc_mass_comp["S_aa"] == blk.Saa_mapping[t]

        @self.Expression(self.flowsheet().time, doc="S_F concentration at step 3")
        def SF_AS3(blk, t):
            return (
                self.properties_in[0].conc_mass_comp["S_F"]
                - blk.Ssu_mapping[t] * 1000
                - blk.Saa_mapping[t] * 1000
            )

        @self.Expression(self.flowsheet().time, doc="S_NH4 concentration at step 3")
        def SNH4_AS3(blk, t):
            return (
                blk.SNH4_AS2[t]
                + self.properties_in[0].conc_mass_comp["S_F"]
                * blk.config.inlet_reaction_package.i_NSF
                - blk.Saa_mapping[t] * self.N_aa * 1000 * mw_n
            )

        @self.Expression(self.flowsheet().time, doc="S_PO4 concentration at step 3")
        def SPO4_AS3(blk, t):
            return (
                blk.SPO4_AS2[t]
                + self.properties_in[0].conc_mass_comp["S_F"]
                * blk.config.inlet_reaction_package.i_PSF
            )

        @self.Expression(self.flowsheet().time, doc="S_IC concentration at step 3")
        def SIC_AS3(blk, t):
            return (
                blk.SIC_AS2[t]
                + self.properties_in[0].conc_mass_comp["S_F"]
                * blk.config.inlet_reaction_package.i_CSF
                - blk.Ssu_mapping[t]
                * blk.config.outlet_reaction_package.Ci["S_su"]
                * 1000
                * mw_c
                - blk.Saa_mapping[t]
                * blk.config.outlet_reaction_package.Ci["S_aa"]
                * 1000
                * mw_c
            )

        # -------------------------------------------Step 4----------------------------------------------------------------#
        self.f_sI_xc = Var(
            initialize=1e-9,
            units=pyunits.dimensionless,
            domain=PositiveReals,
            doc="Soluble inerts from composites",
        )
        self.f_xI_xc = Var(
            initialize=0.1,
            units=pyunits.dimensionless,
            domain=PositiveReals,
            doc="Particulate inerts from composites",
        )
        self.f_ch_xc = Var(
            initialize=0.275,
            units=pyunits.dimensionless,
            domain=PositiveReals,
            doc="Carbohydrates from composites",
        )
        self.f_pr_xc = Var(
            initialize=0.275,
            units=pyunits.dimensionless,
            domain=PositiveReals,
            doc="Proteins from composites",
        )
        self.f_li_xc = Var(
            initialize=0.35,
            units=pyunits.dimensionless,
            domain=PositiveReals,
            doc="Proteins from composites",
        )
        if self.config.inlet_reaction_package.config.decay_switch == DecaySwitch.on:

            @self.Expression(
                self.flowsheet().time, doc="Biomass concentration at step 4"
            )
            # TODO: include sulfur (X_SO) in this expression after implementing sulfur-extension
            def biomass(blk, t):
                return (
                    blk.XH_AS2[t]
                    + self.properties_in[0].conc_mass_comp["X_PAO"]
                    + self.properties_in[0].conc_mass_comp["X_AUT"]
                )

            @self.Expression(self.flowsheet().time, doc="S_I concentration at step 4")
            def SI_AS4(blk, t):
                return (
                    self.properties_in[0].conc_mass_comp["S_I"]
                    + blk.biomass[t] * self.f_sI_xc
                )

            @self.Constraint(self.flowsheet().time, doc="S_I concentration output")
            def SI_output(blk, t):
                return (
                    blk.properties_out[t].conc_mass_comp["S_I"] == blk.SI_AS4[t] / 1000
                )

            @self.Expression(self.flowsheet().time, doc="S_NH4 concentration at step 4")
            def SNH4_AS4(blk, t):
                return (
                    blk.SNH4_AS3[t]
                    + blk.biomass[t] * blk.config.inlet_reaction_package.i_NBM
                    - blk.biomass[t]
                    * self.f_sI_xc
                    * blk.config.inlet_reaction_package.i_NSI
                    - blk.biomass[t]
                    * self.f_xI_xc
                    * blk.config.inlet_reaction_package.i_NSI
                    - blk.biomass[t]
                    * self.f_pr_xc
                    * blk.config.outlet_reaction_package.Ni["X_pr"]
                    * mw_n
                )

            @self.Expression(self.flowsheet().time, doc="S_PO4 concentration at step 4")
            def SPO4_AS4(blk, t):
                return (
                    blk.SPO4_AS3[t]
                    + blk.biomass[t] * blk.config.inlet_reaction_package.i_PBM
                    - blk.biomass[t]
                    * self.f_sI_xc
                    * blk.config.inlet_reaction_package.i_PSI
                    - blk.biomass[t]
                    * self.f_xI_xc
                    * blk.config.inlet_reaction_package.i_PXI
                    - blk.biomass[t]
                    * self.f_ch_xc
                    * 1e-9
                    * pyunits.kmol
                    / pyunits.kg
                    * mw_p
                    - blk.biomass[t]
                    * self.f_li_xc
                    * blk.config.outlet_reaction_package.Pi["X_li"]
                    * mw_p
                )

            @self.Expression(self.flowsheet().time, doc="S_PO4 concentration at step 4")
            def SIC_AS4(blk, t):
                return (
                    blk.SIC_AS3[t]
                    + blk.biomass[t] * blk.config.inlet_reaction_package.i_CXB
                    - blk.biomass[t]
                    * self.f_sI_xc
                    * blk.config.inlet_reaction_package.i_CSI
                    - blk.biomass[t]
                    * self.f_xI_xc
                    * blk.config.inlet_reaction_package.i_CXI
                    - blk.biomass[t]
                    * self.f_pr_xc
                    * blk.config.outlet_reaction_package.Ci["X_pr"]
                    * mw_c
                    - blk.biomass[t]
                    * self.f_ch_xc
                    * blk.config.outlet_reaction_package.Ci["X_ch"]
                    * mw_c
                    - blk.biomass[t]
                    * self.f_li_xc
                    * blk.config.outlet_reaction_package.Ci["X_li"]
                    * mw_c
                )

            @self.Expression(self.flowsheet().time, doc="X_I concentration at step 4")
            def XI_AS4(blk, t):
                return (
                    self.properties_in[0].conc_mass_comp["X_I"]
                    + blk.biomass[t] * self.f_xI_xc
                )

            @self.Constraint(self.flowsheet().time, doc="X_I concentration output")
            def XI_output(blk, t):
                return (
                    blk.properties_out[t].conc_mass_comp["X_I"] == blk.XI_AS4[t] / 1000
                )

            # @self.Expression(doc="X_H concentration at step 4")
            # def XH_AS4():
            #     return 1e-9 * pyunits.kg / pyunits.m**3
            # TODO: Consider removing these cases where vars set conc to 0 since they are usually unused
            # TODO: These were left b/c this first draft follows Flores-Alsina as close as possible
            self.XH_AS4 = Var(
                initialize=1e-9,
                units=pyunits.kg / pyunits.m**3,
                domain=PositiveReals,
                doc="X_H concentration at step 4",
            )

            # @self.Expression(self.flowsheet().time, doc="X_PAO concentration at step 4")
            # def XPAO_AS4():
            #     return 1e-9 * pyunits.kg / pyunits.m**3
            self.XPAO_AS4 = Var(
                initialize=1e-9,
                units=pyunits.kg / pyunits.m**3,
                domain=PositiveReals,
                doc="X_PAO concentration at step 4",
            )

            @self.Constraint(self.flowsheet().time, doc="X_PAO concentration output")
            def XPAO_output(blk, t):
                return (
                    blk.properties_out[t].conc_mass_comp["X_PAO"]
                    == self.XPAO_AS4 / 1000
                )

            # @self.Expression(self.flowsheet().time, doc="X_AUT concentration at step 4")
            # def XAUT_AS4():
            #     return 1e-9 * pyunits.kg / pyunits.m**3
            self.XAUT_AS4 = Var(
                initialize=1e-9,
                units=pyunits.kg / pyunits.m**3,
                domain=PositiveReals,
                doc="X_AUT concentration at step 4",
            )

        elif self.config.inlet_reaction_package.config.decay_switch == DecaySwitch.off:

            @self.Expression(
                self.flowsheet().time, doc="Biomass concentration at step 4"
            )
            # TODO: include sulfur (X_SO) in this expression after implementing sulfur-extension
            def biomass(blk, t):
                return blk.XH_AS2[t] + self.properties_in[0].conc_mass_comp["X_AUT"]

            @self.Expression(self.flowsheet().time, doc="S_I concentration at step 4")
            def SI_AS4(blk, t):
                return (
                    self.properties_in[0].conc_mass_comp["S_I"]
                    + blk.biomass[t] * self.f_sI_xc
                )

            @self.Constraint(self.flowsheet().time, doc="S_I concentration output")
            def SI_output(blk, t):
                return (
                    blk.properties_out[t].conc_mass_comp["S_I"] == blk.SI_AS4[t] / 1000
                )

            @self.Expression(self.flowsheet().time, doc="S_NH4 concentration at step 4")
            def SNH4_AS4(blk, t):
                return (
                    blk.SNH4_AS3[t]
                    + blk.biomass[t] * blk.config.inlet_reaction_package.i_NBM
                    - blk.biomass[t]
                    * self.f_sI_xc
                    * blk.config.inlet_reaction_package.i_NSI
                    - blk.biomass[t]
                    * self.f_xI_xc
                    * blk.config.inlet_reaction_package.i_NSI
                    - blk.biomass[t]
                    * self.f_pr_xc
                    * blk.config.outlet_reaction_package.Ni["X_pr"]
                    * mw_n
                )

            @self.Expression(self.flowsheet().time, doc="S_PO4 concentration at step 4")
            def SPO4_AS4(blk, t):
                return (
                    blk.SPO4_AS3[t]
                    + blk.biomass[t] * blk.config.inlet_reaction_package.i_PBM
                    - blk.biomass[t]
                    * self.f_sI_xc
                    * blk.config.inlet_reaction_package.i_PSI
                    - blk.biomass[t]
                    * self.f_xI_xc
                    * blk.config.inlet_reaction_package.i_PXI
                    - blk.biomass[t]
                    * self.f_ch_xc
                    * 1e-9
                    * pyunits.kmol
                    / pyunits.kg
                    * mw_p
                    - blk.biomass[t]
                    * self.f_li_xc
                    * blk.config.outlet_reaction_package.Pi["X_li"]
                    * mw_p
                )

            @self.Expression(self.flowsheet().time, doc="S_PO4 concentration at step 4")
            def SIC_AS4(blk, t):
                return (
                    blk.SIC_AS3[t]
                    + blk.biomass[t] * blk.config.inlet_reaction_package.i_CXB
                    - blk.biomass[t]
                    * self.f_sI_xc
                    * blk.config.inlet_reaction_package.i_CSI
                    - blk.biomass[t]
                    * self.f_xI_xc
                    * blk.config.inlet_reaction_package.i_CXI
                    - blk.biomass[t]
                    * self.f_pr_xc
                    * blk.config.outlet_reaction_package.Ci["X_pr"]
                    * mw_c
                    - blk.biomass[t]
                    * self.f_ch_xc
                    * blk.config.outlet_reaction_package.Ci["X_ch"]
                    * mw_c
                    - blk.biomass[t]
                    * self.f_li_xc
                    * blk.config.outlet_reaction_package.Ci["X_li"]
                    * mw_c
                )

            @self.Expression(self.flowsheet().time, doc="X_I concentration at step 4")
            def XI_AS4(blk, t):
                return (
                    self.properties_in[0].conc_mass_comp["X_I"]
                    + blk.biomass[t] * self.f_xI_xc
                )

            @self.Constraint(self.flowsheet().time, doc="X_I concentration output")
            def XI_output(blk, t):
                return (
                    blk.properties_out[t].conc_mass_comp["X_I"] == blk.XI_AS4[t] / 1000
                )

            # @self.Expression(doc="X_H concentration at step 4")
            # def XH_AS4():
            #     return 1e-9 * pyunits.kg / pyunits.m**3
            self.XH_AS4 = Var(
                initialize=1e-9,
                units=pyunits.kg / pyunits.m**3,
                domain=PositiveReals,
                doc="X_H concentration at step 4",
            )

            @self.Expression(self.flowsheet().time, doc="X_PAO concentration at step 4")
            def XPAO_AS4():
                return self.properties_in[0].conc_mass_comp["X_PAO"]

            @self.Expression(self.flowsheet().time, doc="X_PP concentration at step 4")
            def XPP_AS4():
                return self.properties_in[0].conc_mass_comp["X_PP"]

            @self.Expression(self.flowsheet().time, doc="X_PHA concentration at step 4")
            def XPHA_AS4():
                return self.properties_in[0].conc_mass_comp["X_PHA"]

            # @self.Expression(self.flowsheet().time, doc="X_AUT concentration at step 4")
            # def XAUT_AS4():
            #     return 1e-9 * pyunits.kg / pyunits.m**3
            self.XAUT_AS4 = Var(
                initialize=1e-9,
                units=pyunits.kg / pyunits.m**3,
                domain=PositiveReals,
                doc="X_AUT concentration at step 4",
            )

        # -------------------------------------------Step 5----------------------------------------------------------------#
        self.i_NXS = Var(
            initialize=0.03352,
            units=pyunits.dimensionless,
            domain=PositiveReals,
            doc="N content of slowly biodegradable substrate X_S, [kg N/kg COD]",
        )
        self.i_PXS = Var(
            initialize=0.00559,
            units=pyunits.dimensionless,
            domain=PositiveReals,
            doc="P content of slowly biodegradable substrate X_S, [kg P/kg COD]",
        )
        self.N_pr = Var(
            initialize=0.0079034,
            units=pyunits.kmol / pyunits.kg,
            domain=PositiveReals,
            doc="Nitrogen content of X_pr",
        )
        self.X_ND = Var(
            initialize=self.properties_in[0].conc_mass_comp["X_S"] * self.i_NXS,
            units=pyunits.kg / pyunits.m**3,
            domain=PositiveReals,
            doc="Nitrogen demand for particulate inerts",
        )
        self.X_PD = Var(
            initialize=self.properties_in[0].conc_mass_comp["X_S"] * self.i_PXS,
            units=pyunits.kg / pyunits.m**3,
            domain=PositiveReals,
            doc="Phosphorus demand for particulate inerts",
        )
        self.XN_org = Var(
            self.flowsheet().time,
            initialize=self.X_ND / (self.N_pr * mw_n),
            units=pyunits.kg / pyunits.m**3,
            domain=PositiveReals,
            doc="Organic nitrogen from particulate inerts",
        )

        @self.Constraint(
            self.flowsheet().time,
            doc="Nitrogen demand for soluble inerts divided by nitrogen content in amino acids",
        )
        def ReqCOD_X(blk, t):
            return blk.XN_org[t] == self.X_ND - (
                blk.config.outlet_reaction_package.Ni["X_pr"] * mw_n
            )

        # TODO: Double check if 1000 should be here and throughout the rest (just them going from g/m3 to kg/m3?)
        @self.Expression(self.flowsheet().time, doc="Carbohydrates mapping")
        def Xch_mapping(blk, t):
            return Expr_if(
                blk.XN_org[t] >= self.properties_in[0].conc_mass_comp["X_S"],
                1e-9 * pyunits.kg / pyunits.m**3,
                (self.properties_in[0].conc_mass_comp["X_S"] - self.XN_org[t])
                * 0.4
                / 1000,
            )

        @self.Expression(self.flowsheet().time, doc="Protein mapping")
        def Xpr_mapping(blk, t):
            return Expr_if(
                blk.XN_org[t] >= self.properties_in[0].conc_mass_comp["X_S"],
                blk.SF_AS3[t] / 1000,
                self.XN_org[t] / 1000,
            )

        @self.Expression(self.flowsheet().time, doc="Lipids mapping")
        def Xli_mapping(blk, t):
            return Expr_if(
                blk.XN_org[t] >= self.properties_in[0].conc_mass_comp["X_S"],
                1e-9 * pyunits.kg / pyunits.m**3,
                (self.properties_in[0].conc_mass_comp["X_S"] - self.XN_org[t])
                * 0.6
                / 1000,
            )

        if self.config.inlet_reaction_package.config.decay_switch == DecaySwitch.on:

            @self.Constraint(
                self.flowsheet().time,
                doc="Carbohydrates output",
            )
            def Xch_output(blk, t):
                return (
                    blk.properties_out[t].conc_mass_comp["X_ch"]
                    == blk.Xch_mapping[t] + blk.biomass[t] * self.f_ch_xc / 1000
                )

            @self.Constraint(
                self.flowsheet().time,
                doc="Protein output",
            )
            def Xpr_output(blk, t):
                return (
                    blk.properties_out[t].conc_mass_comp["X_pr"]
                    == blk.Xpr_mapping[t] + blk.biomass[t] * self.f_pr_xc / 1000
                )

            @self.Constraint(
                self.flowsheet().time,
                doc="Lipids output",
            )
            def Xli_output(blk, t):
                return (
                    blk.properties_out[t].conc_mass_comp["X_li"]
                    == blk.Xli_mapping[t] + blk.biomass[t] * self.f_li_xc / 1000
                )

            @self.Expression(self.flowsheet().time, doc="S_NH4 concentration at step 5")
            def SNH4_AS5(blk, t):
                return (
                    blk.SNH4_AS4[t]
                    + self.properties_in[0].conc_mass_comp["X_S"]
                    * blk.config.inlet_reaction_package.i_NXS
                    - blk.Xpr_mapping[t]
                    * blk.config.outlet_reaction_package.Ni["X_pr"]
                    * 1000
                    * mw_n
                )

            # TODO: Default implementation outputs this as kmol/m3 - I assume we want kg/m3?
            @self.Constraint(self.flowsheet().time, doc="S_NH4 concentration output")
            def SIN_output(blk, t):
                return (
                    # blk.properties_out[t].conc_mol_comp["S_IN"]
                    # == blk.SNH4_AS5[t] / (mw_n * 1000)
                    blk.properties_out[t].conc_mass_comp["S_IN"]
                    == blk.SNH4_AS5[t] / 1000
                )

            @self.Expression(self.flowsheet().time, doc="S_PO4 concentration at step 5")
            def SPO4_AS5(blk, t):
                return (
                    blk.SPO4_AS4[t]
                    + self.properties_in[0].conc_mass_comp["X_S"]
                    * blk.config.inlet_reaction_package.i_PXS
                    - blk.Xch_mapping[t]
                    * 1e-9
                    * pyunits.kmol
                    / pyunits.kg
                    * 1000
                    * mw_p
                    - blk.Xli_mapping[t]
                    * blk.config.outlet_reaction_package.Pi["X_li"]
                    * 1000
                    * mw_p
                )

            @self.Expression(self.flowsheet().time, doc="S_IC concentration at step 5")
            def SIC_AS5(blk, t):
                return (
                    blk.SIC_AS4[t]
                    + self.properties_in[0].conc_mass_comp["S_F"]
                    * blk.config.inlet_reaction_package.i_CXS
                    - blk.Xch_mapping[t]
                    * blk.config.outlet_reaction_package.Ci["X_ch"]
                    * 1000
                    * mw_c
                    - blk.Xpr_mapping[t]
                    * blk.config.outlet_reaction_package.Ci["X_pr"]
                    * 1000
                    * mw_c
                    - blk.Xli_mapping[t]
                    * blk.config.outlet_reaction_package.Ci["X_li"]
                    * 1000
                    * mw_c
                )

            # @self.Expression(self.flowsheet().time, doc="X_S concentration at step 5")
            # def XS_AS5():
            #     return 1e-9 * pyunits.kg / pyunits.m**3
            self.XS_AS5 = Var(
                initialize=1e-9,
                units=pyunits.kg / pyunits.m**3,
                domain=PositiveReals,
                doc="X_S concentration at step 5",
            )

        elif self.config.inlet_reaction_package.config.decay_switch == DecaySwitch.off:

            @self.Constraint(
                self.flowsheet().time,
                doc="Carbohydrates output",
            )
            def Xch_output(blk, t):
                return (
                    blk.properties_out[t].conc_mass_comp["X_ch"]
                    == blk.Xch_mapping[t] + blk.biomass[t] * self.f_ch_xc / 1000
                )

            @self.Constraint(
                self.flowsheet().time,
                doc="Protein output",
            )
            def Xpr_output(blk, t):
                return (
                    blk.properties_out[t].conc_mass_comp["X_pr"]
                    == blk.Xpr_mapping[t] + blk.biomass[t] * self.f_pr_xc / 1000
                )

            @self.Constraint(
                self.flowsheet().time,
                doc="Lipids output",
            )
            def Xli_output(blk, t):
                return (
                    blk.properties_out[t].conc_mass_comp["X_li"]
                    == blk.Xli_mapping[t] + blk.biomass[t] * self.f_li_xc / 1000
                )

            @self.Expression(self.flowsheet().time, doc="S_NH4 concentration at step 5")
            def SNH4_AS5(blk, t):
                return (
                    blk.SNH4_AS4[t]
                    + self.properties_in[0].conc_mass_comp["X_S"]
                    * blk.config.inlet_reaction_package.i_NXS
                    - blk.Xpr_mapping[t]
                    * blk.config.outlet_reaction_package.Ni["X_pr"]
                    * 1000
                    * mw_n
                )

            # TODO: Default implementation outputs this as kmol/m3 - I assume we want kg/m3?
            @self.Constraint(self.flowsheet().time, doc="S_NH4 concentration output")
            def SIN_output(blk, t):
                return (
                    # blk.properties_out[t].conc_mol_comp["S_IN"]
                    # == blk.SNH4_AS5[t] / (mw_n * 1000)
                    blk.properties_out[t].conc_mass_comp["S_IN"]
                    == blk.SNH4_AS5[t] / 1000
                )

            @self.Expression(self.flowsheet().time, doc="S_PO4 concentration at step 5")
            def SPO4_AS5(blk, t):
                return (
                    blk.SPO4_AS4[t]
                    + self.properties_in[0].conc_mass_comp["X_S"]
                    * blk.config.inlet_reaction_package.i_PXS
                    - blk.Xch_mapping[t]
                    * 1e-9
                    * pyunits.kmol
                    / pyunits.kg
                    * 1000
                    * mw_p
                    - blk.Xli_mapping[t]
                    * blk.config.outlet_reaction_package.Pi["X_li"]
                    * 1000
                    * mw_p
                )

            # TODO: Default implementation outputs this as kmol/m3 - I assume we want kg/m3?
            @self.Constraint(self.flowsheet().time, doc="S_IP concentration output")
            def SIP_output(blk, t):
                return (
                    # blk.properties_out[t].conc_mol_comp["S_IP"]
                    # == blk.SPO4_AS5[t] / (mw_p * 1000)
                    blk.properties_out[t].conc_mass_comp["S_IP"]
                    == blk.SPO4_AS5[t] / 1000
                )

            @self.Expression(self.flowsheet().time, doc="S_IC concentration at step 5")
            def SIC_AS5(blk, t):
                return (
                    blk.SIC_AS4[t]
                    + self.properties_in[0].conc_mass_comp["S_F"]
                    * blk.config.inlet_reaction_package.i_CXS
                    - blk.Xch_mapping[t]
                    * blk.config.outlet_reaction_package.Ci["X_ch"]
                    * 1000
                    * mw_c
                    - blk.Xpr_mapping[t]
                    * blk.config.outlet_reaction_package.Ci["X_pr"]
                    * 1000
                    * mw_c
                    - blk.Xli_mapping[t]
                    * blk.config.outlet_reaction_package.Ci["X_li"]
                    * 1000
                    * mw_c
                )

            # TODO: Default implementation outputs this as kmol/m3 - I assume we want kg/m3?
            @self.Constraint(self.flowsheet().time, doc="S_IC concentration output")
            def SIC_output(blk, t):
                return (
                    # blk.properties_out[t].conc_mol_comp["S_IC"]
                    # == blk.SIC_AS5[t] / (mw_c * 1000)
                    blk.properties_out[t].conc_mass_comp["S_IC"]
                    == blk.SIC_AS5[t] / 1000
                )

            # @self.Expression(self.flowsheet().time, doc="X_S concentration at step 5")
            # def XS_AS5():
            #     return 1e-9 * pyunits.kg / pyunits.m**3
            self.XS_AS5 = Var(
                initialize=1e-9,
                units=pyunits.kg / pyunits.m**3,
                domain=PositiveReals,
                doc="X_S concentration at step 5",
            )

            # @self.Expression(self.flowsheet().time, doc="X_PAO concentration at step 5")
            # def XPAO_AS5():
            #     return self.properties_in[0].conc_mass_comp["X_PAO"]
            self.XPAO_AS5 = Var(
                initialize=1e-9,
                units=pyunits.kg / pyunits.m**3,
                domain=PositiveReals,
                doc="X_PAO concentration at step 5",
            )

            # @self.Expression(self.flowsheet().time, doc="X_PP concentration at step 4")
            # def XPP_AS5():
            #     return self.properties_in[0].conc_mass_comp["X_PP"]
            self.XPP_AS5 = Var(
                initialize=1e-9,
                units=pyunits.kg / pyunits.m**3,
                domain=PositiveReals,
                doc="X_PP concentration at step 5",
            )

            # @self.Expression(self.flowsheet().time, doc="X_PHA concentration at step 4")
            # def XPHA_AS5():
            #     return self.properties_in[0].conc_mass_comp["X_PHA"]
            self.XPHA_AS5 = Var(
                initialize=1e-9,
                units=pyunits.kg / pyunits.m**3,
                domain=PositiveReals,
                doc="X_PHA concentration at step 5",
            )

        # -------------------------------------------Step 6----------------------------------------------------------------#
        self.f_XPHA_Sva = Var(
            initialize=0.1,
            units=pyunits.dimensionless,
            domain=PositiveReals,
            doc="Valerate from Polyhydroxyalkanoates",
        )

        self.f_XPHA_Sbu = Var(
            initialize=0.1,
            units=pyunits.dimensionless,
            domain=PositiveReals,
            doc="Butyrate from Polyhydroxyalkanoates",
        )

        self.f_XPHA_Spro = Var(
            initialize=0.4,
            units=pyunits.dimensionless,
            domain=PositiveReals,
            doc="Propionate from Polyhydroxyalkanoates",
        )

        self.f_XPHA_Sac = Var(
            initialize=0.4,
            units=pyunits.dimensionless,
            domain=PositiveReals,
            doc="Acetate from Polyhydroxyalkanoates",
        )

        if self.config.inlet_reaction_package.config.decay_switch == DecaySwitch.on:

            # @self.Expression(self.flowsheet().time, doc="X_PP concentration at step 6")
            # def XPP_AS6():
            #     return 1e-9 * pyunits.kg / pyunits.m**3
            self.XPP_AS6 = Var(
                initialize=1e-9,
                units=pyunits.kg / pyunits.m**3,
                domain=PositiveReals,
                doc="X_PP concentration at step 6",
            )

            # TODO: Default implementation outputs this as kmol/m3 - I assume we want kg/m3?
            @self.Constraint(self.flowsheet().time, doc="X_PP concentration output")
            def XPP_output(blk, t):
                # return blk.properties_out[t].conc_mol_comp["X_PP"] == self.XPP_AS6 / (1000 * mw_p)
                return (
                    blk.properties_out[t].conc_mass_comp["X_PP"] == self.XPP_AS6 / 1000
                )

            # @self.Expression(self.flowsheet().time, doc="X_PP concentration at step 6")
            # def XPHA_AS6():
            #     return 1e-9 * pyunits.kg / pyunits.m**3
            self.XPHA_AS6 = Var(
                initialize=1e-9,
                units=pyunits.kg / pyunits.m**3,
                domain=PositiveReals,
                doc="X_PHA concentration at step 6",
            )

            @self.Constraint(self.flowsheet().time, doc="X_PHA concentration output")
            def XPHA_output(blk, t):
                return (
                    blk.properties_out[t].conc_mass_comp["X_PHA"]
                    == self.XPHA_AS6 / 1000
                )

            @self.Expression(
                self.flowsheet().time,
                doc="Total valerate concentration",
            )
            def Sva_AS6(blk):
                return blk.XPHA_AS6 * self.f_XPHA_Sva / 1000

            @self.Constraint(
                self.flowsheet().time,
                doc="Total valerate concentration output",
            )
            def Sva_output(blk, t):
                return blk.properties_out[t].conc_mass_comp["S_va"] == blk.Sva_AS6

            @self.Expression(
                self.flowsheet().time,
                doc="Total butyrate concentration",
            )
            def Sbu_AS6(blk):
                return blk.XPHA_AS6 * self.f_XPHA_Sbu / 1000

            @self.Constraint(
                self.flowsheet().time,
                doc="Total butyrate concentration output",
            )
            def Sbu_output(blk, t):
                return blk.properties_out[t].conc_mass_comp["S_bu"] == blk.Sbu_AS6

            @self.Expression(
                self.flowsheet().time,
                doc="Total propionate concentration",
            )
            def Spro_AS6(blk):
                return blk.XPHA_AS6 * self.f_XPHA_Spro / 1000

            @self.Constraint(
                self.flowsheet().time,
                doc="Total propionate concentration output",
            )
            def Spro_output(blk, t):
                return blk.properties_out[t].conc_mass_comp["S_pro"] == blk.Spro_AS6

            @self.Expression(
                self.flowsheet().time,
                doc="Total acetate concentration",
            )
            def Sac_AS6(blk):
                return blk.XPHA_AS6 * self.f_XPHA_Sac / 1000

            @self.Constraint(
                self.flowsheet().time,
                doc="Total acetate concentration output",
            )
            def Sac_output(blk, t):
                return (
                    blk.properties_out[t].conc_mass_comp["S_ac"]
                    == blk.Sac_AS6 + blk.SA_AS2[t] / 1000
                )

            @self.Expression(self.flowsheet().time, doc="S_PO4 concentration at step 6")
            def SPO4_AS6(blk, t):
                return blk.SPO4_AS5[t] + self.properties_in[0].conc_mass_comp["X_PP"]

            # TODO: Default implementation outputs this as kmol/m3 - I assume we want kg/m3?
            @self.Constraint(self.flowsheet().time, doc="S_IP concentration output")
            def SIP_output(blk, t):
                return (
                    # blk.properties_out[t].conc_mol_comp["S_IP"]
                    # == blk.SPO4_AS6[t] / (mw_p * 1000)
                    blk.properties_out[t].conc_mass_comp["S_IP"]
                    == blk.SPO4_AS6[t] / 1000
                )

            @self.Expression(self.flowsheet().time, doc="S_IC concentration at step 6")
            def SIC_AS6(blk, t):
                return (
                    blk.SIC_AS5[t]
                    + self.properties_in[0].conc_mass_comp["X_PHA"] * 0.3
                    - blk.Sva_AS6
                    * blk.config.outlet_reaction_package.Ci["S_va"]
                    * mw_c
                    * 1000
                    - blk.Sbu_AS6
                    * blk.config.outlet_reaction_package.Ci["S_bu"]
                    * mw_c
                    * 1000
                    - blk.Spro_AS6
                    * blk.config.outlet_reaction_package.Ci["S_pro"]
                    * mw_c
                    * 1000
                    - blk.Sac_AS6
                    * blk.config.outlet_reaction_package.Ci["S_ac"]
                    * mw_c
                    * 1000
                )

            # TODO: Default implementation outputs this as kmol/m3 - I assume we want kg/m3?
            @self.Constraint(self.flowsheet().time, doc="S_IC concentration output")
            def SIC_output(blk, t):
                return (
                    # blk.properties_out[t].conc_mol_comp["S_IC"]
                    # == blk.SIC_AS6[t] / (mw_c * 1000)
                    blk.properties_out[t].conc_mass_comp["S_IC"]
                    == blk.SIC_AS6[t] / 1000
                )

            @self.Expression(self.flowsheet().time, doc="S_K concentration at step 6")
            def SK_AS6(blk):
                return (
                    self.properties_in[0].conc_mass_comp["S_K"]
                    + blk.config.inlet_reaction_package.K_XPP
                    * self.properties_in[0].conc_mass_comp["X_PP"]
                )

            # TODO: Default implementation outputs this as kmol/m3 - I assume we want kg/m3?
            @self.Constraint(
                self.flowsheet().time,
                doc="Total potassium concentration output",
            )
            def SK_output(blk, t):
                return (
                    # blk.properties_out[t].conc_mass_comp["S_K"]
                    # == blk.SK_AS6[t] / (1000 * mw_k)
                    blk.properties_out[t].conc_mass_comp["S_K"]
                    == blk.SK_AS6[t] / 1000
                )

            @self.Expression(self.flowsheet().time, doc="S_Mg concentration at step 6")
            def SMg_AS6(blk):
                return (
                    self.properties_in[0].conc_mass_comp["S_Mg"]
                    + blk.config.inlet_reaction_package.Mg_XPP
                    * self.properties_in[0].conc_mass_comp["X_PP"]
                )

            # TODO: Default implementation outputs this as kmol/m3 - I assume we want kg/m3?
            @self.Constraint(
                self.flowsheet().time,
                doc="Total magnesium concentration output",
            )
            def SMg_output(blk, t):
                return (
                    # blk.properties_out[t].conc_mass_comp["S_Mg"]
                    # == blk.SMg_AS6[t] / (1000 * mw_mg)
                    blk.properties_out[t].conc_mass_comp["S_Mg"]
                    == blk.SMg_AS6[t] / 1000
                )

        elif self.config.inlet_reaction_package.config.decay_switch == DecaySwitch.off:

            @self.Constraint(
                self.flowsheet().time,
                doc="Total valerate concentration",
            )
            def Sva_output(blk, t):
                return (
                    blk.properties_out[t].conc_mass_comp["S_va"]
                    == 1e-9 * pyunits.kg / pyunits.m**3
                )

            @self.Constraint(
                self.flowsheet().time,
                doc="Total butyrate concentration",
            )
            def Sbu_output(blk, t):
                return (
                    blk.properties_out[t].conc_mass_comp["S_bu"]
                    == 1e-9 * pyunits.kg / pyunits.m**3
                )

            @self.Constraint(
                self.flowsheet().time,
                doc="Total propionate concentration",
            )
            def Spro_output(blk, t):
                return (
                    blk.properties_out[t].conc_mass_comp["S_pro"]
                    == 1e-9 * pyunits.kg / pyunits.m**3
                )

            @self.Constraint(
                self.flowsheet().time,
                doc="Total acetate concentration",
            )
            def Sac_output(blk, t):
                return (
                    blk.properties_out[t].conc_mass_comp["S_ac"]
                    == 1e-9 * pyunits.kg / pyunits.m**3
                )

            # @self.Expression(self.flowsheet().time, doc="X_PAO concentration at step 6")
            # def XPAO_AS6():
            #     return self.properties_in[0].conc_mass_comp["X_PAO"]
            self.XPAO_AS6 = Var(
                initialize=self.properties_in[0].conc_mass_comp["X_PAO"],
                units=pyunits.kg / pyunits.m**3,
                domain=PositiveReals,
                doc="X_PAO concentration at step 6",
            )

            @self.Constraint(self.flowsheet().time, doc="X_PAO concentration output")
            def XPAO_output(blk, t):
                return (
                    blk.properties_out[t].conc_mass_comp["X_PAO"]
                    == self.XPAO_AS6 / 1000
                )

            # @self.Expression(self.flowsheet().time, doc="X_PP concentration at step 6")
            # def XPP_AS6():
            #     return self.properties_in[0].conc_mass_comp["X_PP"]
            self.XPP_AS6 = Var(
                initialize=self.properties_in[0].conc_mass_comp["X_PP"],
                units=pyunits.kg / pyunits.m**3,
                domain=PositiveReals,
                doc="X_PP concentration at step 6",
            )

            # TODO: Default implementation outputs this as kmol/m3 - I assume we want kg/m3?
            @self.Constraint(self.flowsheet().time, doc="X_PP concentration output")
            def XPP_output(blk, t):
                # return blk.properties_out[t].conc_mol_comp["X_PP"] == blk.XPP_AS6[
                #     t
                # ] / (1000 * mw_p)
                return (
                    blk.properties_out[t].conc_mass_comp["X_PP"] == self.XPP_AS6 / 1000
                )

            # @self.Expression(self.flowsheet().time, doc="X_PHA concentration at step 6")
            # def XPHA_AS6():
            #     return self.properties_in[0].conc_mass_comp["X_PHA"]
            self.XPHA_AS6 = Var(
                initialize=self.properties_in[0].conc_mass_comp["X_PHA"],
                units=pyunits.kg / pyunits.m**3,
                domain=PositiveReals,
                doc="X_PHA concentration at step 6",
            )

            @self.Constraint(self.flowsheet().time, doc="X_PHA concentration output")
            def XPHA_output(blk, t):
                return (
                    blk.properties_out[t].conc_mass_comp["X_PHA"]
                    == self.XPHA_AS6 / 1000
                )

            # TODO: Default implementation outputs this as kmol/m3 - I assume we want kg/m3?
            @self.Constraint(
                self.flowsheet().time,
                doc="Total magnesium concentration output",
            )
            def SK_output(blk, t):
                return (
                    # blk.properties_out[t].conc_mass_comp["S_K"]
                    # == self.properties_in[0].conc_mass_comp["S_K"] / (1000 * mw_k)
                    blk.properties_out[t].conc_mass_comp["S_K"]
                    == self.properties_in[0].conc_mass_comp["S_K"] / 1000
                )

            # TODO: Default implementation outputs this as kmol/m3 - I assume we want kg/m3?
            @self.Constraint(
                self.flowsheet().time,
                doc="Total magnesium concentration output",
            )
            def SMg_output(blk, t):
                return (
                    # blk.properties_out[t].conc_mass_comp["S_Mg"]
                    # == self.properties_in[0].conc_mass_comp["S_Mg"] / (1000 * mw_mg)
                    blk.properties_out[t].conc_mass_comp["S_Mg"]
                    == self.properties_in[0].conc_mass_comp["S_Mg"] / 1000
                )

        self.zero_flow_components = Set(
            initialize=[
                "S_fa",
                "S_h2",
                "S_ch4",
                "X_su",
                "X_aa",
                "X_fa",
                "X_c4",
                "X_pro",
                "X_ac",
                "X_h2",
            ]
        )

        @self.Constraint(
            self.flowsheet().time,
            self.zero_flow_components,
            doc="Components with no flow equation",
        )
        def return_zero_flow_comp(blk, t, i):
            return (
                blk.properties_out[t].conc_mass_comp[i]
                == 1e-9 * pyunits.kg / pyunits.m**3
            )

        # TODO: Relationship carried over from ASM1/ADM1 - see perf_plant_AD_ss.m for Flores-Alsina interpretation (Scat = S_K + S_Mg)
        @self.Constraint(
            self.flowsheet().time,
            doc="Anions balance",
        )
        def return_an(blk, t):
            return (
                blk.properties_out[t].anions
                == blk.properties_out[t].conc_mass_comp["S_IN"] / mw_n
            )

        @self.Constraint(
            self.flowsheet().time,
            doc="Cations balance",
        )
        def return_cat(blk, t):
            return (
                blk.properties_out[t].cations
                == blk.properties_out[t].conc_mass_comp["S_IC"] / mw_c
            )

    def initialize_build(
        self,
        state_args_in=None,
        state_args_out=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        This method calls the initialization method of the state blocks.

        Keyword Arguments:
            state_args_in : a dict of arguments to be passed to the inlet
                property package (to provide an initial state for
                initialization (see documentation of the specific
                property package) (default = None).
            state_args_out : a dict of arguments to be passed to the outlet
                property package (to provide an initial state for
                initialization (see documentation of the specific
                property package) (default = None).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None, use default solver)

        Returns:
            None
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")

        # Create solver
        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        # Initialize state block
        flags = self.properties_in.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_in,
            hold_state=True,
        )

        self.properties_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_out,
        )

        if degrees_of_freedom(self) != 0:
            raise Exception(
                f"{self.name} degrees of freedom were not 0 at the beginning "
                f"of initialization. DoF = {degrees_of_freedom(self)}"
            )

        with idaeslog.solver_log(init_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)

        self.properties_in.release_state(flags=flags, outlvl=outlvl)

        init_log.info(f"Initialization Complete: {idaeslog.condition(res)}")

        if not check_optimal_termination(res):
            raise InitializationError(
                f"{self.name} failed to initialize successfully. Please check "
                f"the output logs for more information."
            )
