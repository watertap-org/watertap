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
    Var,
    Param,
    check_optimal_termination,
    units as pyunits,
)

# Import IDAES cores
import idaes.logger as idaeslog
from idaes.core import declare_process_block_class
import idaes.core.util.scaling as iscale
from idaes.core.util.exceptions import InitializationError

from watertap.unit_models.ion_exchange.ion_exchange_base import (
    IonExchangeBaseData,
    IonExchangeType,
)
from watertap.core.solvers import get_solver
from watertap.core.util.initialization import interval_initializer

__author__ = "Kurban Sitterley"


@declare_process_block_class("IonExchangeDemin")
class IonExchangeDeminData(IonExchangeBaseData):
    """
    Ion exchange demineralization
    """

    # For comparison to EPA-WBS
    # self.bv = BV_regen
    # self.breakthrough_time = regen_days
    # self.bed_volume_total = resin_volume_total
    # self.loading_rate = comm_load_rate
    # self.bed_volume_total = resin_volume_total

    def build(self):
        super().build()

        self.ion_exchange_type = IonExchangeType.demineralize
        self.config.add_steady_state_approximation = False

        prop_in = self.process_flow.properties_in[0]
        regen = self.regeneration_stream[0]

        solutes = (
            self.config.property_package.anion_set
            | self.config.property_package.cation_set
        )

        for j in self.config.property_package.neutral_set:
            self.process_flow.mass_transfer_term[:, "Liq", j].fix(0)
            regen.get_material_flow_terms("Liq", j).fix(0)

        self.process_flow.mass_transfer_term[:, "Liq", "H2O"].fix(0)

        # Set EPA-WBS bounds for EBCT
        self.ebct.setlb(60)
        self.ebct.setub(480)

        @self.Expression(doc="Equivalent flow rate for cations in feed")
        def flow_equiv_cation(b):
            return sum(
                pyunits.convert(
                    prop_in.flow_equiv_phase_comp["Liq", c]
                    * pyunits.s
                    * pyunits.mol**-1,
                    to_units=pyunits.dimensionless,
                )
                for c in self.config.property_package.cation_set
            )

        @self.Expression(doc="Equivalent flow rate for anions in feed")
        def flow_equiv_anion(b):
            return sum(
                pyunits.convert(
                    prop_in.flow_equiv_phase_comp["Liq", c]
                    * pyunits.s
                    * pyunits.mol**-1,
                    to_units=pyunits.dimensionless,
                )
                for c in self.config.property_package.anion_set
            )

        @self.Expression(doc="Equivalent flow rate for total in feed")
        def flow_equiv_total(b):
            return pyunits.convert(
                b.flow_equiv_cation + b.flow_equiv_anion,
                to_units=pyunits.dimensionless,
            )

        @self.Expression(doc="Charge ratio for cation exchange")
        def charge_ratio_cx(b):
            return b.flow_equiv_cation / b.flow_equiv_total

        @self.Expression(doc="Charge ratio for anion exchange")
        def charge_ratio_ax(b):
            return b.flow_equiv_anion / b.flow_equiv_total

        @self.Expression(doc="Waste to influent flow ratio")
        def frac_waste_to_inlet_flow(b):
            return pyunits.convert(
                regen.flow_vol_phase["Liq"] / prop_in.flow_vol_phase["Liq"],
                to_units=pyunits.dimensionless,
            )

        self.removal_efficiency = Param(
            solutes,
            initialize=0.99,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Removal efficiency",
        )

        # typical total capacity for synthetic resin is between 0.14-2.3 eq/L
        # typical operating/usable capacity for synthetic resin is between 0.69-1.15 eq/L
        # Wachinski, chap. 3

        self.resin_capacity_ax = Var(
            initialize=0.75,
            bounds=(0.14, 2.3),
            units=pyunits.mol / pyunits.liter,
            doc="Operating capacity of the anion exchange resin in equivalents per liter resin",
        )

        self.resin_capacity_cx = Var(
            initialize=0.75,
            bounds=(0.14, 2.3),
            units=pyunits.mol / pyunits.liter,
            doc="Operating capacity of the cation exchange resin in equivalents per liter resin",
        )

        self.resin_capacity_op = Var(
            initialize=0.75,
            bounds=(0.14, 2.3),
            units=pyunits.mol / pyunits.liter,
            doc="Operating capacity of the mixed bed resin",
        )

        @self.Constraint(doc="Effective resin capacity")
        def eq_resin_capacity_op(b):
            return (
                b.resin_capacity_op
                == b.resin_capacity_cx * b.charge_ratio_cx
                + b.resin_capacity_ax * b.charge_ratio_ax
            )

        @self.Constraint(doc="Bed volumes at breakthrough point")
        def eq_bv(b):
            return b.bv == pyunits.convert(
                (prop_in.flow_vol_phase["Liq"] * b.resin_capacity_op)
                / sum(
                    prop_in.flow_equiv_phase_comp["Liq", j] * b.removal_efficiency[j]
                    for j in solutes
                ),
                to_units=pyunits.dimensionless,
            )

        @self.Constraint()
        def eq_breakthrough_time(b):
            return b.breakthrough_time == pyunits.convert(
                (b.bed_depth * b.bv) / b.loading_rate, to_units=pyunits.s
            )

        @self.Constraint(solutes, doc="Mass transfer term for control volume")
        def eq_mass_transfer_term(b, j):
            return (
                b.process_flow.mass_transfer_term[0, "Liq", j]
                == -1 * prop_in.flow_mass_phase_comp["Liq", j] * b.removal_efficiency[j]
            )

        @self.Constraint(solutes, doc="Regeneration stream mass flow")
        def eq_mass_transfer_regen(b, j):
            return (
                regen.get_material_flow_terms("Liq", j)
                == -b.process_flow.mass_transfer_term[0, "Liq", j]
            )

        @self.Constraint(doc="Total bed volume required based on volumetric flow")
        def eq_bed_volume_required(b):
            return b.bed_volume == pyunits.convert(
                b.flow_per_column * b.ebct, to_units=pyunits.m**3
            )

        @self.Constraint(doc="Regeneration stream flow rate")
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

        state_args_regen = dict()
        state_dict_regen = self.regeneration_stream[
            self.flowsheet().config.time.first()
        ].define_port_members()

        for k in state_dict_regen.keys():
            if state_dict_regen[k].is_indexed():
                state_args_regen[k] = {}
                for m in state_dict_regen[k].keys():
                    state_args_regen[k][m] = state_dict_regen[k][m].value
            else:
                state_args_regen[k] = state_dict_regen[k].value

        self.regeneration_stream.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_regen,
        )

        init_log.info("Initialization Step 1b Complete.")
        interval_initializer(self)

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

        if iscale.get_scaling_factor(self.resin_capacity_ax) is None:
            iscale.set_scaling_factor(self.resin_capacity_ax, 1)

        if iscale.get_scaling_factor(self.resin_capacity_cx) is None:
            iscale.set_scaling_factor(self.resin_capacity_cx, 1)

        if iscale.get_scaling_factor(self.resin_capacity_op) is None:
            iscale.set_scaling_factor(self.resin_capacity_op, 1)
