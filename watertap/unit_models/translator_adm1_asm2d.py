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
Translator block representing the ADM1/ASM2d interface.
This is copied from the Generic template for a translator block.

Assumptions:
     * Steady-state only

Model formulated from:

Flores-Alsina, X., Solon, K., Mbamba, C.K., Tait, S., Gernaey, K.V., Jeppsson, U. and Batstone, D.J., 2016.
Modelling phosphorus (P), sulfur (S) and iron (Fe) interactions for dynamic simulations of anaerobic digestion processes.
Water Research, 95, pp.370-382.
"""

# Import Pyomo libraries
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool

# Import IDAES cores
from idaes.core import declare_process_block_class, UnitModelBlockData
from idaes.core.util.config import (
    is_physical_parameter_block,
    is_reaction_parameter_block,
)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.solvers import get_solver
import idaes.logger as idaeslog

from pyomo.environ import (
    Reference,
    Var,
    value,
    Constraint,
    Param,
    units as pyunits,
    check_optimal_termination,
    exp,
    Set,
    PositiveReals,
)

__author__ = "Chenyu Wang, Marcus Holly"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("Translator_ADM1_ASM2D")
class TranslatorData(UnitModelBlockData):
    """
    Standard Translator Block Class
    """

    CONFIG = ConfigBlock()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Translator blocks are always steady-state.""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Translator blocks do not contain holdup.""",
        ),
    )
    CONFIG.declare(
        "has_phase_equilibrium",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Phase equilibrium construction flag",
            doc="""Indicates whether terms for phase equilibrium should be
    constructed,
    **default** = False.
    **Valid values:** {
    **True** - include phase equilibrium terms
    **False** - exclude phase equilibrium terms.}""",
        ),
    )
    CONFIG.declare(
        "outlet_state_defined",
        ConfigValue(
            default=True,
            domain=Bool,
            description="Indicated whether outlet state will be fully defined",
            doc="""Indicates whether unit model will fully define outlet state.
    If False, the outlet property package will enforce constraints such as sum
    of mole fractions and phase equilibrium.
    **default** - True.
    **Valid values:** {
    **True** - outlet state will be fully defined,
    **False** - outlet property package should enforce sumation and equilibrium
    constraints.}""",
        ),
    )
    CONFIG.declare(
        "inlet_property_package",
        ConfigValue(
            default=None,
            domain=is_physical_parameter_block,
            description="Property package to use for incoming stream",
            doc="""Property parameter object used to define property
    calculations for the incoming stream,
    **default** - None.
    **Valid values:** {
    **PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "inlet_property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property package "
            "of the incoming stream",
            doc="""A ConfigBlock with arguments to be passed to the property
    block associated with the incoming stream,
    **default** - None.
    **Valid values:** {
    see property package for documentation.}""",
        ),
    )
    CONFIG.declare(
        "outlet_property_package",
        ConfigValue(
            default=None,
            domain=is_physical_parameter_block,
            description="Property package to use for outgoing stream",
            doc="""Property parameter object used to define property
    calculations for the outgoing stream,
    **default** - None.
    **Valid values:** {
    **PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "outlet_property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property package "
            "of the outgoing stream",
            doc="""A ConfigBlock with arguments to be passed to the property
    block associated with the outgoing stream,
    **default** - None.
    **Valid values:** {
    see property package for documentation.}""",
        ),
    )
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
        super(TranslatorData, self).build()

        # TODO: check this parameter later
        # self.i_ec = Param(
        #     initialize=0.06,
        #     units=pyunits.dimensionless,
        #     mutable=True,
        #     doc="Nitrogen inert content",
        # )

        # Add State Blocks
        self.properties_in = self.config.inlet_property_package.build_state_block(
            self.flowsheet().time,
            doc="Material properties in incoming stream",
            defined_state=True,
            has_phase_equilibrium=False,
            **self.config.inlet_property_package_args
        )

        self.properties_out = self.config.outlet_property_package.build_state_block(
            self.flowsheet().time,
            doc="Material properties in outgoing stream",
            defined_state=self.config.outlet_state_defined,
            has_phase_equilibrium=False,
            **self.config.outlet_property_package_args
        )

        # Add ports
        self.add_port(name="inlet", block=self.properties_in, doc="Inlet Port")
        self.add_port(name="outlet", block=self.properties_out, doc="Outlet Port")

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

        @self.Constraint(
            self.flowsheet().time,
            doc="Equality S_I equation",
        )
        def eq_SI_conc(blk, t):
            return (
                blk.properties_out[t].conc_mass_comp["S_I"]
                == blk.properties_in[t].conc_mass_comp["S_I"]
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

        # TODO: ADM1 does not have S_IP
        # @self.Constraint(
        #     self.flowsheet().time,
        #     doc="Equality S_PO4 equation",
        # )
        # def eq_SPO4_conc(blk, t):
        #     return (
        #         blk.properties_out[t].conc_mass_comp["S_IP"]
        #         == blk.properties_in[t].conc_mass_comp["S_PO4"]
        #     )

        # TODO: No S_IC in ASM2D
        # @self.Constraint(
        #     self.flowsheet().time,
        #     doc="Equality S_IC equation",
        # )
        # def eq_SIC_conc(blk, t):
        #     return (
        #         blk.properties_out[t].conc_mass_comp["S_IC"]
        #         == blk.properties_in[t].conc_mass_comp["S_IC"]
        #     )

        @self.Constraint(
            self.flowsheet().time,
            doc="Equality X_I equation",
        )
        def eq_XI_conc(blk, t):
            return (
                blk.properties_out[t].conc_mass_comp["X_I"]
                == blk.properties_in[t].conc_mass_comp["X_I"]
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

        # TODO: ADM1 model does not has X_PP as poly-phosphates
        # Assume both sides are 0 - Page 373
        # @self.Constraint(
        #     self.flowsheet().time,
        #     doc="Equality X_PP equation",
        # )
        # def eq_XPP_conc(blk, t):
        #     return (
        #         blk.properties_out[t].conc_mass_comp["X_PP"]
        #         == blk.properties_in[t].conc_mass_comp["X_PP"]
        #     )

        # TODO: ADM1 model does not has X_PHA
        # @self.Constraint(
        #     self.flowsheet().time,
        #     doc="Equality X_PHA equation",
        # )
        # def eq_XPHA_conc(blk, t):
        #     return (
        #         blk.properties_out[t].conc_mass_comp["X_PHA"]
        #         == blk.properties_in[t].conc_mass_comp["X_PHA"]
        #     )

        # TODO: check if we track S_SO4, S_Na, S_K, S_Cl, S_Ca, S_Mg, X_Ca2(PO4)3, X_MgNH4PO4

        @self.Constraint(
            self.flowsheet().time,
            doc="Equality alkalinity equation",
        )
        def return_Salk(blk, t):
            return (
                blk.properties_out[t].alkalinity
                == blk.properties_in[t].conc_mass_comp["S_IC"] / mw_c
            )

        # TODO: check S_ALK
        self.zero_flow_components = Set(
            initialize=[
                "S_N2",
                "S_NO3",
                "S_O2",
                "S_PO4",
                "X_AUT",
                "X_H",
                "X_MeOH",
                "X_MeP",
                "X_PAO",
                "X_PHA",
                "X_PP",
                "X_TSS",
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
                == 1e-6 * pyunits.kg / pyunits.m**3
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

        if degrees_of_freedom(self) == 0:
            with idaeslog.solver_log(init_log, idaeslog.DEBUG) as slc:
                res = opt.solve(self, tee=slc.tee)

            init_log.info("Initialization Complete {}.".format(idaeslog.condition(res)))
        else:
            init_log.warning(
                "Initialization incomplete. Degrees of freedom "
                "were not zero. Please provide sufficient number "
                "of constraints linking the state variables "
                "between the two state blocks."
            )

        self.properties_in.release_state(flags=flags, outlvl=outlvl)
