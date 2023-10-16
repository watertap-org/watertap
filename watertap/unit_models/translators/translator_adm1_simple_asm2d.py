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
Translator block representing the ADM1/Simple ASM2d interface.
This is copied from the Generic template for a translator block.

Assumptions:
     * Steady-state only

Model formulated from:

Flores-Alsina, X., Solon, K., Mbamba, C.K., Tait, S., Gernaey, K.V., Jeppsson, U. and Batstone, D.J., 2016.
Modelling phosphorus (P), sulfur (S) and iron (Fe) interactions for dynamic simulations of anaerobic digestion processes.
Water Research, 95, pp.370-382.
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
import idaes.core.util.scaling as iscale

from idaes.core.util.exceptions import InitializationError

from pyomo.environ import (
    units as pyunits,
    check_optimal_termination,
    Set,
)

__author__ = "Chenyu Wang, Marcus Holly, Xinhong Liu"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("Translator_ADM1_Simple_ASM2D")
class TranslatorDataADM1SimpleASM2D(TranslatorData):
    """
    Translator block representing the ADM1/Simple ASM2D interface
    """

    CONFIG = TranslatorData.CONFIG()
    CONFIG.declare(
        "reaction_package",
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
        "reaction_package_args",
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
        super(TranslatorDataADM1SimpleASM2D, self).build()

        mw_c = 12 * pyunits.kg / pyunits.kmol

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

        self.readily_biodegradable = Set(initialize=["S_su", "S_aa", "S_fa"])

        @self.Constraint(
            self.flowsheet().time,
            doc="Equality S_F equation",
        )
        def eq_SF_conc(blk, t):
            return blk.properties_out[t].conc_mass_comp["S_F"] == sum(
                blk.properties_in[t].conc_mass_comp[i]
                for i in blk.readily_biodegradable
            )

        self.readily_biodegradable2 = Set(initialize=["S_va", "S_bu", "S_pro", "S_ac"])

        @self.Constraint(
            self.flowsheet().time,
            doc="Equality S_A equation",
        )
        def eq_SA_conc(blk, t):
            return blk.properties_out[t].conc_mass_comp["S_A"] == sum(
                blk.properties_in[t].conc_mass_comp[i]
                for i in blk.readily_biodegradable2
            )

        self.unchanged_component = Set(initialize=["S_I", "X_I", "X_PP", "X_PHA"])

        @self.Constraint(
            self.flowsheet().time,
            self.unchanged_component,
            doc="Equality equation for unchanged components",
        )
        def eq_unchanged_conc(blk, t, i):
            return (
                blk.properties_out[t].conc_mass_comp[i]
                == blk.properties_in[t].conc_mass_comp[i]
            )

        @self.Constraint(
            self.flowsheet().time,
            doc="Equality S_NH4 equation",
        )
        def eq_SNH4_conc(blk, t):
            return (
                blk.properties_out[t].conc_mass_comp["S_NH4"]
                == blk.properties_in[t].conc_mass_comp["S_IN"]
            )

        @self.Constraint(
            self.flowsheet().time,
            doc="Equality S_PO4 equation",
        )
        def eq_SPO4_conc(blk, t):
            return (
                blk.properties_out[t].conc_mass_comp["S_PO4"]
                == blk.properties_in[t].conc_mass_comp["S_IP"]
            )

        @self.Constraint(
            self.flowsheet().time,
            doc="Equality alkalinity equation",
        )
        def return_Salk(blk, t):
            return (
                blk.properties_out[t].alkalinity
                == blk.properties_in[t].conc_mass_comp["S_IC"] / mw_c
            )

        self.slowly_biodegradable = Set(
            initialize=[
                "X_ch",
                "X_pr",
                "X_li",
            ]
        )

        @self.Constraint(
            self.flowsheet().time,
            doc="Equality X_S equation",
        )
        def eq_XS_conc(blk, t):
            return blk.properties_out[t].conc_mass_comp["X_S"] == sum(
                blk.properties_in[t].conc_mass_comp[i] for i in blk.slowly_biodegradable
            )

        self.zero_flow_components = Set(
            initialize=[
                "S_N2",
                "S_NO3",
                "S_O2",
                "X_AUT",
                "X_H",
                "X_PAO",
                "X_TSS",
                "X_MeOH",
                "X_MeP",
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
                == 1e-10 * pyunits.kg / pyunits.m**3
            )

        if (
            self.config.outlet_property_package.config.additional_solute_list
            is not None
        ):

            @self.Constraint(
                self.flowsheet().time,
                self.config.outlet_property_package.config.additional_solute_list,
                doc="Equality ASM2D additional solute equation",
            )
            def eq_ASM2D_additional_conc(blk, t, i):
                return (
                    blk.properties_out[t].conc_mass_comp[i]
                    == blk.properties_in[t].conc_mass_comp[i]
                )

        iscale.set_scaling_factor(self.properties_out[0].flow_vol, 1e5)

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
