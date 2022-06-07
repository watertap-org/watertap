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

from pyomo.environ import (
    Block,
    Set,
    Var,
    Param,
    Suffix,
    NonNegativeReals,
    PositiveReals,
    Reference,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.dae import ContinuousSet, DerivativeVar

from enum import Enum, auto

from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    MaterialBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    useDefault,
    MaterialFlowBasis,
)
from idaes.core.util import get_solver
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

__author__ = "Hunter Barber"

_log = idaeslog.getLogger(__name__)

# ---------------------------------------------------------------------
# TODO: Start of adding user options for specific variables to reduce required input
"""
class SurfaceDiffusionCoeffVal(Enum):
    specified = auto()  # surface diffusion coefficient is a user-specified value
    calculated = auto() # pressure drop across membrane channel is calculated
"""

# ---------------------------------------------------------------------
@declare_process_block_class("GAC")
class GACData(UnitModelBlockData):
    """
    Initial Granular Activated Carbon Model
    """

    # CONFIG are options for the unit model, this simple model only has the mandatory config options
    CONFIG = ConfigBlock()

    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
    **default** = False. The filtration unit does not support dynamic
    behavior, thus this must be False.""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Indicates whether holdup terms should be constructed or not.
    **default** - False. The filtration unit does not have defined volume, thus
    this must be False.""",
        ),
    )
    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property calculations,
    **default** - useDefault.
    **Valid values:** {
    **useDefault** - use default package from parent model or flowsheet,
    **PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property packages",
            doc="""A ConfigBlock with arguments to be passed to a property block(s)
    and used when constructing these,
    **default** - None.
    **Valid values:** {
    see property package for documentation.}""",
        ),
    )
    CONFIG.declare(
        "material_balance_type",
        ConfigValue(
            default=MaterialBalanceType.useDefault,
            domain=In(MaterialBalanceType),
            description="Material balance construction flag",
            doc="""Indicates what type of mass balance should be constructed,
    **default** - MaterialBalanceType.useDefault.
    **Valid values:** {
    **MaterialBalanceType.useDefault - refer to property package for default
    balance type
    **MaterialBalanceType.none** - exclude material balances,
    **MaterialBalanceType.componentPhase** - use phase component balances,
    **MaterialBalanceType.componentTotal** - use total component balances,
    **MaterialBalanceType.elementTotal** - use total element balances,
    **MaterialBalanceType.total** - use total material balance.}""",
        ),
    )
    CONFIG.declare(
        "momentum_balance_type",
        ConfigValue(
            default=MomentumBalanceType.pressureTotal,
            domain=In(MomentumBalanceType),
            description="Momentum balance construction flag",
            doc="""Indicates what type of momentum balance should be constructed,
        **default** - MomentumBalanceType.pressureTotal.
        **Valid values:** {
        **MomentumBalanceType.none** - exclude momentum balances,
        **MomentumBalanceType.pressureTotal** - single pressure balance for material,
        **MomentumBalanceType.pressurePhase** - pressure balances for each phase,
        **MomentumBalanceType.momentumTotal** - single momentum balance for material,
        **MomentumBalanceType.momentumPhase** - momentum balances for each phase.}""",
        ),
    )
    CONFIG.declare(
        "has_pressure_change",
        ConfigValue(
            default=False,
            domain=In([False]),  # TODO: domain=In([True, False])
            description="Pressure change term construction flag",
            doc="""Indicates whether terms for pressure change should be
        constructed,
        **default** - False.
        **Valid values:** {
        **True** - include pressure change terms,
        **False** - exclude pressure change terms.}""",
        ),
    )

    # ---------------------------------------------------------------------
    def build(self):

        super().build()
        # create blank scaling factors to be populated later
        self.scaling_factor = Suffix(direction=Suffix.EXPORT)
        units_meta = self.config.property_package.get_metadata().get_derived_units

        # build control volume
        self.treatwater = ControlVolume0DBlock(
            default={
                "dynamic": False,
                "has_holdup": False,
                "property_package": self.config.property_package,
                "property_package_args": self.config.property_package_args,
            }
        )
        self.treatwater.add_state_blocks(has_phase_equilibrium=False)
        self.treatwater.add_material_balances(
            balance_type=self.config.material_balance_type, has_mass_transfer=True
        )
        self.treatwater.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=False,
        )

        @self.treatwater.Constraint(
            self.flowsheet().config.time, doc="isothermal assumption for water flow"
        )
        def eq_isothermal(b, t):
            return b.properties_in[t].temperature == b.properties_out[t].temperature

        # add block for spent GAC
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = False  # permeate block is not an inlet
        self.removal = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of spent gac",
            default=tmp_dict,
        )

        @self.Constraint(
            self.flowsheet().config.time, doc="Isothermal assumption for spent GAC"
        )
        def eq_isothermal_adsorbate(b, t):
            return b.treatwater.properties_in[t].temperature == b.removal[t].temperature

        @self.Constraint(
            self.flowsheet().config.time, doc="Isobaric assumption for spent GAC"
        )
        def eq_isobaric_adsorbate(b, t):
            return b.treatwater.properties_in[t].pressure == b.removal[t].pressure

        # Add ports
        self.add_inlet_port(name="inlet", block=self.treatwater)
        self.add_outlet_port(name="outlet", block=self.treatwater)
        self.add_port(name="spent_gac", block=self.removal)

        # ---------------------------------------------------------------------
        # variable declaration
        # mass transfer TODO: Add model capacity for multiple solutes
        self.contam_removal_frac = Var(
            self.config.property_package.solute_set,
            initialize=0.9,
            bounds=(0, 1),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Desired removal fraction of key contaminant",
        )

        self.freund_ninv = Var(
            initialize=0.5,
            bounds=(0, 5),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Freundlich 1/n parameter",
        )

        # TODO: Figure out non-fixed exponent
        # TODO: Determine whether unit correction procedure interferes with regressed k values
        self.freund_k = Var(
            initialize=10,
            bounds=(0, 1000),
            domain=NonNegativeReals,
            # TODO: Correct variable units
            units=pyunits.dimensionless,  # ((units_meta("length") ** 3) * units_meta("mass") ** -1) ** 0.8316,
            doc="Freundlich k parameter",
        )

        self.eps_bed = Var(
            initialize=0.5,
            bounds=(0, 1),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="GAC bed void fraction",
        )

        self.ebct = Var(
            initialize=500,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("time"),
            doc="Empty bed contact time",
        )

        self.thru = Var(
            initialize=10,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Mass throughput",
        )

        self.tau = Var(
            initialize=100,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("time"),
            doc="Packed bed contact time",
        )

        self.replace_time = Var(
            initialize=1e8,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("time"),
            doc="Replace time for saturated GAC",
        )

        # ---------------------------------------------------------------------
        # Steady state assumption variables
        self.num_active_beds = Var(
            initialize=5,
            bounds=(1, None),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Number of active beds in parallel, one additional bed offline for replacement",
        )

        # ---------------------------------------------------------------------
        # GAC particle properties
        self.replace_saturation_frac = Var(
            self.config.property_package.solute_set,
            initialize=0.9,
            bounds=(0, 1),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Fraction of carbon saturation before replacement",
        )

        self.particle_dens_app = Var(
            initialize=1000,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("mass") * units_meta("length") ** -3,
            doc="GAC particle apparent density",
        )

        self.particle_dens_bulk = Var(
            initialize=500,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("mass") * units_meta("length") ** -3,
            doc="GAC particle apparent density",
        )

        self.particle_dp = Var(
            initialize=0.001,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("length"),
            doc="GAC particle diameter",
        )

        self.kf = Var(
            initialize=1e-5,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("length") * units_meta("time") ** -1,
            doc="Liquid phase film transfer rate",
        )

        self.ds = Var(
            initialize=1e-12,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("length") ** 2 * units_meta("time") ** -1,
            doc="Surface diffusion coefficient",
        )

        # ---------------------------------------------------------------------
        # Minimum conditions to achieve a constant pattern solution
        self.min_st = Var(
            initialize=10,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Minimum Stanton number",
        )

        self.min_ebct = Var(
            initialize=1000,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("time"),
            doc="Minimum empty bed contact time",
        )

        self.min_tau = Var(
            initialize=1000,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("time"),
            doc="Minimum packed bed contact time",
        )

        self.min_time = Var(
            initialize=1e8,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("time"),
            doc="Minimum elapsed time",
        )

        # ---------------------------------------------------------------------
        # Constants in regressed equations
        self.a0 = Var(
            initialize=1,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Stanton equation parameter",
        )

        self.a1 = Var(
            initialize=1,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Stanton equation parameter",
        )

        self.b0 = Var(
            initialize=0.1,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Throughput equation parameter",
        )

        self.b1 = Var(
            initialize=0.1,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Throughput equation parameter",
        )

        self.b2 = Var(
            initialize=0.1,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Throughput equation parameter",
        )

        self.b3 = Var(
            initialize=0.1,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Throughput equation parameter",
        )

        self.b4 = Var(
            initialize=0.1,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Throughput equation parameter",
        )

        # ---------------------------------------------------------------------
        # Intermediate variables
        self.dg = Var(
            initialize=1e5,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Solute distribution parameter",
        )

        self.bi = Var(
            initialize=10,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Biot number",
        )

        # ---------------------------------------------------------------------
        # Performance Variables
        self.particle_mass_req = Var(
            initialize=1e-2,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("mass") * units_meta("time") ** -1,
            doc="GAC usage and replacement rate",
        )

        # ---------------------------------------------------------------------
        # TODO: Add support mole or mass based property packs
        # TODO: Ensure other phases in addition to 'Liq' can be handled

        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.solute_set,
            doc="Mass transfer term for solutes",
        )
        def eq_mass_transfer_solute(b, t, j):
            return (
                b.contam_removal_frac[j]
                * b.treatwater.properties_in[t].get_material_flow_terms("Liq", j)
                == -b.treatwater.mass_transfer_term[t, "Liq", j]
            )

        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.solvent_set,
            doc="No mass transfer of solvents",
        )
        def eq_mass_transfer_solvent(b, t, j):
            return b.treatwater.mass_transfer_term[t, "Liq", j] == 0
            # will set to flow_mol_phase_comp.lb of 1e-8

        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.component_list,
            doc="Contaminant absorbed",
        )
        def eq_mass_transfer_absorbed(b, t, j):
            return (
                b.removal[t].get_material_flow_terms("Liq", j)
                == -b.treatwater.mass_transfer_term[t, "Liq", j]
            )

        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.solute_set,
            doc="Solute distribution parameter",
        )
        def eq_dg(b, t, j):
            freund_k_units = (
                (units_meta("length") ** 3) * units_meta("mass") ** -1
            ) ** b.freund_ninv
            return b.dg * b.eps_bed * b.treatwater.properties_in[
                t
            ].conc_mass_phase_comp[
                "Liq", j
            ] == b.particle_dens_app * b.freund_k * freund_k_units * (
                b.treatwater.properties_in[t].conc_mass_phase_comp["Liq", j]
                ** b.freund_ninv
            ) * (
                1 - b.eps_bed
            )

        @self.Constraint(doc="Biot number")
        def eq_bi(b):
            return b.bi * b.ds * b.dg * b.eps_bed == b.kf * (b.particle_dp / 2) * (
                1 - b.eps_bed
            )

        @self.Constraint(
            doc="Minimum Stanton number to achieve constant pattern solution"
        )
        def eq_min_st_cps(b):
            return b.min_st == b.a0 * b.bi + b.a1

        @self.Constraint(
            doc="Minimum empty bed contact time to achieve constant pattern solution"
        )
        def eq_min_ebct_cps(b):
            return b.min_ebct * (1 - b.eps_bed) * b.kf == b.min_st * (b.particle_dp / 2)

        @self.Constraint(
            doc="Minimum packed bed contact time to achieve constant pattern solution"
        )
        def eq_min_tau_cps(b):
            return b.min_tau == b.eps_bed * b.min_ebct

        @self.Constraint(doc="residence time")
        def eq_tau(b):
            return b.tau == b.eps_bed * b.ebct

        # ---------------------------------------------------------------------
        # loop for beds in parallel
        @self.Constraint(
            self.config.property_package.solute_set,
            doc="Throughput based on 5-parameter regression",
        )
        def eq_thru(b, j):
            return b.thru == b.b0 + b.b1 * (
                (b.replace_saturation_frac[j]) ** b.b2
            ) + b.b3 / (1.01 - ((b.replace_saturation_frac[j]) ** b.b4))

        @self.Constraint(doc="Minimum elapsed time for constant pattern solution")
        def eq_min_time_cps(b):
            return b.min_time == b.min_tau * (b.dg + 1) * b.thru

        @self.Constraint(doc="Bed replacement time")
        def eq_replacement_time(b):
            return b.replace_time == b.min_time + (b.tau - b.min_tau) * (b.dg + 1)

        # ---------------------------------------------------------------------

        @self.Constraint(doc="Relate void fraction and GAC densities")
        def eq_bed_void_fraction(b):
            return b.eps_bed == 1 - (b.particle_dens_bulk / b.particle_dens_app)

        @self.Constraint(doc="Bed replacement mass required")
        def eq_replacement_mass(b):
            return (
                b.particle_mass_req * b.replace_time
                == b.treatwater.properties_in[0].flow_vol_phase["Liq"]
                * b.particle_dens_bulk
                * b.tau
            )

    # ---------------------------------------------------------------------
    # TODO: Correct initialization procedure from starting at an infeasible point, add robustness for intialization and scaling
    # initialize method
    def initialize_build(
        blk, state_args=None, outlvl=idaeslog.NOTSET, solver=None, optarg=None
    ):
        """
        General wrapper for initialization routines

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = {}).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None)
            solver : str indicating which solver to use during
                     initialization (default = None)

        Returns: None
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")
        # Set solver options
        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        # Set state_args from inlet state
        if state_args is None:
            state_args = {}
            state_dict = blk.treatwater.properties_in[
                blk.flowsheet().config.time.first()
            ].define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        # TODO: Determine if scale of the flowrate matters for this initialization
        # specify conditions to solve at a feasible point during initialization
        for j in blk.config.property_package.solute_set:
            state_args["flow_mol_phase_comp"][("Liq", j)] = 1e-5

        for j in blk.config.property_package.solvent_set:
            state_args["flow_mol_phase_comp"][("Liq", j)] = 100

        # Initialize control volume
        flags = blk.treatwater.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )

        init_log.info_high("Initialization Step 1 Complete.")
        # ---------------------------------------------------------------------
        # Initialize adsorbate port

        blk.removal.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )

        init_log.info_high("Initialization Step 2 Complete.")
        # --------------------------------------------------------------------
        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))
        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.treatwater.release_state(flags, outlvl + 1)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

    # ---------------------------------------------------------------------

    # TODO: def _get_performance_contents(self, time_point=0):
    # TODO: def _get_stream_table_contents(self, time_point=0):
    # TODO: def get_costing(self, module=None, **kwargs):

    # ---------------------------------------------------------------------
    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # overwrite default scaling

        # scaling for gac created variables
        for j in self.contam_removal_frac.keys():
            if iscale.get_scaling_factor(self.contam_removal_frac[j]) is None:
                iscale.set_scaling_factor(self.contam_removal_frac[j], 1)

        if iscale.get_scaling_factor(self.freund_ninv) is None:
            iscale.set_scaling_factor(self.freund_ninv, 1)

        if iscale.get_scaling_factor(self.freund_k) is None:
            iscale.set_scaling_factor(self.freund_k, 1)

        if iscale.get_scaling_factor(self.eps_bed) is None:
            iscale.set_scaling_factor(self.eps_bed, 1)

        if iscale.get_scaling_factor(self.ebct) is None:
            iscale.set_scaling_factor(self.ebct, 1e-2)

        if iscale.get_scaling_factor(self.thru) is None:
            iscale.set_scaling_factor(self.thru, 1)

        if iscale.get_scaling_factor(self.tau) is None:
            iscale.set_scaling_factor(self.tau, 1e-2)

        if iscale.get_scaling_factor(self.replace_time) is None:
            iscale.set_scaling_factor(self.replace_time, 1e-6)

        if iscale.get_scaling_factor(self.num_active_beds) is None:
            iscale.set_scaling_factor(self.num_active_beds, 1)

        for j in self.contam_removal_frac.keys():
            if iscale.get_scaling_factor(self.replace_saturation_frac[j]) is None:
                iscale.set_scaling_factor(self.replace_saturation_frac[j], 1)

        if iscale.get_scaling_factor(self.particle_dens_app) is None:
            iscale.set_scaling_factor(self.particle_dens_app, 1e-2)

        if iscale.get_scaling_factor(self.particle_dens_bulk) is None:
            iscale.set_scaling_factor(self.particle_dens_bulk, 1e-2)

        if iscale.get_scaling_factor(self.particle_dp) is None:
            iscale.set_scaling_factor(self.particle_dp, 1e3)

        if iscale.get_scaling_factor(self.kf) is None:
            iscale.set_scaling_factor(self.kf, 1e5)

        if iscale.get_scaling_factor(self.ds) is None:
            iscale.set_scaling_factor(self.ds, 1e13)

        if iscale.get_scaling_factor(self.min_st) is None:
            iscale.set_scaling_factor(self.min_st, 1e-1)

        if iscale.get_scaling_factor(self.min_ebct) is None:
            iscale.set_scaling_factor(self.min_ebct, 1e-3)

        if iscale.get_scaling_factor(self.min_tau) is None:
            iscale.set_scaling_factor(self.min_tau, 1e-2)

        if iscale.get_scaling_factor(self.min_time) is None:
            iscale.set_scaling_factor(self.min_time, 1e-7)

        iscale.set_scaling_factor(self.a0, 1)
        iscale.set_scaling_factor(self.a1, 1)
        iscale.set_scaling_factor(self.b0, 1)
        iscale.set_scaling_factor(self.b1, 1)
        iscale.set_scaling_factor(self.b2, 1)
        iscale.set_scaling_factor(self.b3, 10)
        iscale.set_scaling_factor(self.b4, 1)

        if iscale.get_scaling_factor(self.dg) is None:
            iscale.set_scaling_factor(self.dg, 1e-4)

        if iscale.get_scaling_factor(self.bi) is None:
            iscale.set_scaling_factor(self.bi, 1)

        if iscale.get_scaling_factor(self.particle_mass_req) is None:
            iscale.set_scaling_factor(self.particle_mass_req, 1e2)

        # (optional) transforming constraints
        """
        for (t, j), c in self.eq_mass_transfer_absorbed.items():
            if j in self.config.property_package.solvent_set:
                sf = 1e-12
                iscale.constraint_scaling_transform(c, sf)


        for (t, j), c in self.eq_mass_transfer_solute.items():
            sf = iscale.get_scaling_factor(
                self.treatwater.properties_in[t].flow_mol_phase_comp["Liq", j]
            )
            iscale.constraint_scaling_transform(c, sf)

        for (t, j), c in self.eq_mass_transfer_solvent.items():
            sf = 1
            iscale.constraint_scaling_transform(c, sf)

        for (t, j), c in self.eq_mass_transfer_absorbed.items():
            if j in self.config.property_package.solute_set:
                sf = iscale.get_scaling_factor(
                    self.treatwater.properties_in[t].flow_mol_phase_comp["Liq", j]
                )
                iscale.constraint_scaling_transform(c, sf)
            if j in self.config.property_package.solvent_set:
                sf = 1
                iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.treatwater.eq_isothermal.items():
            sf = iscale.get_scaling_factor(self.treatwater.properties_in[0].temperature)
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_isothermal_adsorbate.items():
            sf = iscale.get_scaling_factor(self.treatwater.properties_in[0].temperature)
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_isobaric_adsorbate.items():
            sf = iscale.get_scaling_factor(self.treatwater.properties_in[0].pressure)
            iscale.constraint_scaling_transform(c, sf)
        """
