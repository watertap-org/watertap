#################################################################################
# WaterTAP Copyright (c) 2020-2025, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

# Import Pyomo libraries
from pyomo.environ import (
    Param,
    Set,
    Var,
    log,
    check_optimal_termination,
    units as pyunits,
)
from pyomo.common.config import ConfigValue

# Import IDAES cores
from idaes.core import declare_process_block_class
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.core.util.exceptions import InitializationError, ConfigurationError

from watertap.core.solvers import get_solver
from watertap.unit_models.ion_exchange.ion_exchange_base import (
    IonExchangeBaseData,
    IonExchangeType,
    add_ss_approximation,
)

__author__ = "Kurban Sitterley"


@declare_process_block_class("IonExchangeClark")
class IonExchangeClarkData(IonExchangeBaseData):
    """
    Ion exchange Clark model
    """

    CONFIG = IonExchangeBaseData.CONFIG()

    CONFIG.declare(
        "reactive_components",
        ConfigValue(
            default=list(),
            domain=list,
            description="Designates other reactive species",
        ),
    )
    CONFIG.declare(
        "number_traps",
        ConfigValue(
            default=5,
            domain=int,
            description="Number of trapezoids to use for steady-state effluent concentration estimation",
        ),
    )

    def build(self):
        super().build()

        prop_in = self.process_flow.properties_in[0]
        regen = self.regeneration_stream[0]
        comps = self.config.property_package.component_list
        target_component = self.config.target_component
        reactive_components = self.config.reactive_components

        if not target_component in reactive_components:
            # Raise warning?
            reactive_components.append(target_component)

        self.reactive_component_set = Set(initialize=reactive_components)

        inerts = comps - self.reactive_component_set
        self.inert_set = Set(initialize=inerts)

        if len(self.target_component_set) > 1:
            raise ConfigurationError(
                f"IonExchangeClark can only accept a single target ion but {len(self.target_component_set)} were provided."
            )
        if self.config.property_package.charge_comp[target_component].value > 0:
            self.ion_exchange_type = IonExchangeType.cation
        elif self.config.property_package.charge_comp[target_component].value < 0:
            self.ion_exchange_type = IonExchangeType.anion
        else:
            raise ConfigurationError("Target ion must have non-zero charge.")

        for j in inerts:
            self.process_flow.mass_transfer_term[:, "Liq", j].fix(0)
            if j != "H2O":
                regen.get_material_flow_terms("Liq", j).fix(0)

        self.flow_basis = self.process_flow.properties_in[
            self.flowsheet().config.time.first()
        ].get_material_flow_basis()
        self.eps = Param(initialize=1e-4, mutable=True)

        self.c_norm = Var(
            self.reactive_component_set,
            initialize=0.5,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Normalized (relative) effluent concentration of reactive components",
        )

        self.freundlich_n = Var(
            self.reactive_component_set,
            initialize=1.5,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Freundlich isotherm exponent",
        )

        self.mass_transfer_coeff = Var(
            self.reactive_component_set,  # k_T
            initialize=0.001,
            units=pyunits.s**-1,
            bounds=(0, None),
            doc="Mass transfer coefficient for Clark model (kT)",
        )

        self.bv_50 = Var(
            self.reactive_component_set,  # BV_50
            initialize=2e5,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Bed volumes of feed at 50 percent breakthrough",
        )

        @self.Constraint(doc="Bed volumes at breakthrough")
        def eq_bv(b):
            return b.breakthrough_time * b.loading_rate == b.bv * b.bed_depth

        @self.Constraint(
            self.reactive_component_set, doc="Clark equation with fundamental constants"
        )  # Croll et al (2023), Eq.9
        def eq_clark(b, j):
            left_side = (
                (b.mass_transfer_coeff[j] * b.bed_depth * (b.freundlich_n[j] - 1))
                / (b.bv_50[j] * b.loading_rate)
            ) * (b.bv_50[j] - b.bv)

            right_side = log(
                ((1 / b.c_norm[j]) ** (b.freundlich_n[j] - 1) - 1)
                / (2 ** (b.freundlich_n[j] - 1) - 1)
            )
            return left_side - right_side == 0

        @self.Constraint(doc="Total bed volume required based on volumetric flow")
        def eq_bed_volume_required(b):
            return b.bed_volume == pyunits.convert(
                b.flow_per_column * b.ebct, to_units=pyunits.m**3
            )

        if self.config.add_steady_state_approximation:
            add_ss_approximation(self, ix_model_type="clark")
        else:
            # If not using steady state approximation, need mass transfer terms;
            # assume effluent concentration is equal to concentration at breakthrough

            @self.Constraint(
                self.reactive_component_set,
                doc="CV mass transfer term",
            )
            def eq_mass_transfer_term(b, j):
                return (1 - b.c_norm[j]) * prop_in.get_material_flow_terms(
                    "Liq", j
                ) == -b.process_flow.mass_transfer_term[0, "Liq", j]

            @self.Constraint(
                self.reactive_component_set,
                doc="Target component mass transfer for regeneration stream",
            )
            def eq_mass_transfer_regen(b, j):
                return (
                    regen.get_material_flow_terms("Liq", j)
                    == -b.process_flow.mass_transfer_term[0, "Liq", j]
                )

            @self.Constraint(doc="Regeneration stream volumetric flow rate")
            def eq_regen_flow_rate(b):
                return regen.flow_vol_phase["Liq"] == pyunits.convert(
                    b.rinse_flow_rate * (b.rinse_time / b.cycle_time)
                    + b.backwash_flow_rate * (b.backwash_time / b.cycle_time)
                    + b.regen_flow_rate * (b.regeneration_time / b.cycle_time),
                    to_units=pyunits.m**3 / pyunits.s,
                )

    def initialize_build(
        self,
        state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        opt = get_solver(solver, optarg)

        flags = self.process_flow.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
            hold_state=True,
        )
        init_log.info("Initialization Step 1a Complete.")

        # Initialize other state blocks
        # Set state_args from inlet state
        if state_args is None:
            self.state_args = state_args = {}
            state_dict = self.process_flow.properties_in[
                self.flowsheet().config.time.first()
            ].define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        self.regeneration_stream.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )

        init_log.info("Initialization Step 1b Complete.")
        # interval_initializer(self)

        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
            if not check_optimal_termination(res):
                init_log.warning(
                    f"Trouble solving unit model {self.name}, trying one more time"
                )
                res = opt.solve(self, tee=slc.tee)

        init_log.info("Initialization Step 2 {}.".format(idaeslog.condition(res)))

        # Release Inlet state
        self.process_flow.properties_in.release_state(flags, outlvl=outlvl)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {self.name} failed to initialize.")

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        if iscale.get_scaling_factor(self.c_norm) is None:
            iscale.set_scaling_factor(self.c_norm, 10)

        if iscale.get_scaling_factor(self.freundlich_n) is None:
            iscale.set_scaling_factor(self.freundlich_n, 0.1)

        if iscale.get_scaling_factor(self.mass_transfer_coeff) is None:
            iscale.set_scaling_factor(self.mass_transfer_coeff, 10)

        if iscale.get_scaling_factor(self.bv_50) is None:
            iscale.set_scaling_factor(self.bv_50, 1e-5)
