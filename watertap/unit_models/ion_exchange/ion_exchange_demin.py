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
    Var,
    Param,
    value,
    units as pyunits,
)

# Import IDAES cores
from idaes.core import declare_process_block_class

import idaes.core.util.scaling as iscale

from watertap.unit_models.ion_exchange.ion_exchange_base import (
    IonExchangeBaseData,
    IonExchangeType,
)

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

    def build(self):
        super().build()

        prop_in = self.process_flow.properties_in[0]
        regen = self.regeneration_stream[0]

        solutes = (
            self.config.property_package.anion_set
            | self.config.property_package.cation_set
        )
        for j in self.config.property_package.neutral_set:
            self.process_flow.mass_transfer_term[:, "Liq", j].fix(0)
            regen.get_material_flow_terms("Liq", j).fix(0)

        # Set EPA-WBS bounds for EBCT
        self.ebct.setlb(60)
        self.ebct.setub(480)

        self.flow_equiv_cation = sum(
            value(prop_in.flow_equiv_phase_comp["Liq", c])
            for c in self.config.property_package.cation_set
        )
        self.flow_equiv_anion = sum(
            value(prop_in.flow_equiv_phase_comp["Liq", c])
            for c in self.config.property_package.anion_set
        )
        self.flow_equiv_total = self.flow_equiv_cation + self.flow_equiv_anion

        self.ion_exchange_type = IonExchangeType.demineralize

        assert (
            self.flow_equiv_cation + self.flow_equiv_anion
        ) / self.flow_equiv_total == 1
        self.charge_ratio_cx = (
            self.flow_equiv_cation / self.flow_equiv_total * pyunits.dimensionless
        )
        self.charge_ratio_ax = (
            self.flow_equiv_anion / self.flow_equiv_total * pyunits.dimensionless
        )

        self.process_flow.mass_transfer_term[:, "Liq", "H2O"].fix(0)

        self.removal_efficiency = Param(
            solutes,
            initialize=0.9999,
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
            return b.bed_volume_total == pyunits.convert(
                prop_in.flow_vol_phase["Liq"] * b.ebct, to_units=pyunits.m**3
            )

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        if iscale.get_scaling_factor(self.resin_capacity_ax) is None:
            iscale.set_scaling_factor(self.resin_capacity_ax, 1)

        if iscale.get_scaling_factor(self.resin_capacity_cx) is None:
            iscale.set_scaling_factor(self.resin_capacity_cx, 1)

        if iscale.get_scaling_factor(self.resin_capacity_op) is None:
            iscale.set_scaling_factor(self.resin_capacity_op, 1)
