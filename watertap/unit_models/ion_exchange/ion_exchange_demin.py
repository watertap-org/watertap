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

from copy import deepcopy

# Import Pyomo libraries
from pyomo.environ import (
    Set,
    Var,
    check_optimal_termination,
    Param,
    Suffix,
    log,
    units as pyunits,
)

# Import IDAES cores
from idaes.core import declare_process_block_class
from watertap.core.solvers import get_solver
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.constants import Constants
from idaes.core.util.exceptions import InitializationError, ConfigurationError

import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog


from watertap.unit_models.ion_exchange.ion_exchange_base import IonExchangeBaseData

__author__ = "Kurban Sitterley"


@declare_process_block_class("IonExchangeDemin")
class IonExchangeDeminData(IonExchangeBaseData):
    """
    Ion exchange demineralization
    """

    def build(self):
        super().build()

        prop_in = self.process_flow.properties_in[0]
        prop_out = self.process_flow.properties_out[0]
        regen = self.regeneration_stream[0]
        target = self.config.target_component
        comps = self.config.property_package.component_list
        solutes = comps - ["H2O"]
        self.ebct.setlb(60)
        self.ebct.setub(480)
        for j in solutes:
            self.process_flow.mass_transfer_term[:, "Liq", j].fix(0)
            regen.get_material_flow_terms("Liq", j).fix(0)

        # typical operating capacity for resin is between 30-60 kg/m3; Wachinski, chap. 3
        self.resin_capacity_total = Var(
            initialize=1,
            bounds=(0.1374, 2.29),
            units=pyunits.mol / pyunits.liter,
            doc="Total capacity of the resin in mol /L",
        )

        self.flow_equiv_total_in = Var(
            initialize=1, bounds=(0, None), units=pyunits.mol / pyunits.s
        )

        self.flow_equiv_total_out = Var(
            initialize=1, bounds=(0, None), units=pyunits.mol / pyunits.s
        )

        self.resin_capacity_op = Var(
            initialize=0.75,
            bounds=(0.687, 2.1),
            units=pyunits.mol / pyunits.liter,
            doc="Operational (usable) capacity of the resin in mol /L",
        )

        self.removal_efficiency = Var(
            # solutes,
            initialize=0.01,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Removal efficiency",
        )

        self.bv = Var(
            initialize=300,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="BV per cycle",
        )

        @self.Expression()
        def mass_removed_per_day(b):
            return pyunits.convert(
                b.flow_equiv_total_in * (1 * pyunits.day), to_units=pyunits.mol
            ) / (1 * pyunits.day)

        # @self.Constraint()
        # def eq_bed_volume_total(b):
        #     return b.bed_volume_total == pyunits.convert(b.mass_to_be_removed / b.resin_capacity_op,to_units=pyunits.m**3)

        @self.Constraint()
        def eq_flow_equiv_total_in(b):
            return b.flow_equiv_total_in == pyunits.convert(
                prop_in.flow_vol_phase["Liq"]
                * sum(prop_in.conc_equiv_phase_comp["Liq", j] for j in solutes),
                to_units=pyunits.mol / pyunits.s,
            )

        @self.Constraint()
        def eq_bv(b):
            expr = (prop_in.flow_vol_phase["Liq"] * b.resin_capacity_op) / (
                b.flow_equiv_total_in
                # b.mass_removed_per_day
            )
            return b.bv == pyunits.convert(expr, to_units=pyunits.dimensionless)

        @self.Constraint()
        def eq_breakthrough_time(b):
            return b.breakthrough_time == pyunits.convert(
                (b.bed_depth * b.bv) / b.loading_rate, to_units=pyunits.s
            )

        @self.Constraint()
        def eq_bed_vol_tot(b):
            return b.bed_volume_total == b.bed_volume * b.number_columns

        @self.Constraint()
        def eq_bed_vol(b):
            return b.bed_volume_total == pyunits.convert(
                prop_in.flow_vol_phase["Liq"] * b.ebct, to_units=pyunits.m**3
            )

        # @self.Constraint()
        # def eq_flow_equiv_total_out(b):
        #     return b.flow_equiv_total_out == b.removal_efficiency * b.flow_equiv_total_in

        # @self.Expression(doc="Mass of target to be removed per cycle")
        # def mass_removed_per_cycle(b):
        #     return pyunits.convert(
        #         prop_out.flow_mass_phase_comp["Liq", target] * b.breakthrough_time,
        #         to_units=pyunits.kg,
        #     )

        # @self.Constraint(doc="Effluent concentration")
        # def eq_effluent(b):
        #     return (
        #         prop_out.conc_mass_phase_comp["Liq", target]
        #         == b.c_norm[target] * prop_in.conc_mass_phase_comp["Liq", target]
        #     )

        # @self.Constraint(doc="Regeneration stream")

        # @self.Constraint(doc="Total bed volume required")
        # def bed_volume_total_constraint(b):
        #     return b.bed_volume_total == pyunits.convert(
        #         b.mass_removed_per_cycle / b.resin_capacity, to_units=pyunits.m**3
        #     )
