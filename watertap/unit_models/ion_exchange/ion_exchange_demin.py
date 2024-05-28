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
    value,
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


from watertap.unit_models.ion_exchange.ion_exchange_base import IonExchangeBaseData, IonExchangeType

__author__ = "Kurban Sitterley"


@declare_process_block_class("IonExchangeDemin")
class IonExchangeDeminData(IonExchangeBaseData):
    """
    Ion exchange demineralization
    """

    def build(self):
        super().build()

        prop_in = self.process_flow.properties_in[0]
        solutes = [x for x in self.config.property_package.config.solute_list if x != "H2O"]
        self.ebct.setlb(60)
        self.ebct.setub(480)

        self.process_flow.mass_transfer_term[:, "Liq", "H2O"].fix(0)
        
        self.removal_efficiency = Param(
            solutes,
            initialize=0.9999,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Removal efficiency",
        )

        self.resin_capacity_op = Var(
            initialize=0.75,
            bounds=(0.687, 2.1),
            units=pyunits.mol / pyunits.liter,
            doc="Operating (usable) capacity of the resin in equivalents per liter resin",
        ) # typical operating capacity for resin is between 30-60 kg/m3; Wachinski, chap. 3

        self.bv = Var(
            initialize=300,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Bed volumes per cycle",
        )

        @self.Constraint(doc="Bed volumes at breakthrough point")
        def eq_bv(b):
            return b.bv == pyunits.convert(
                (prop_in.flow_vol_phase["Liq"] * b.resin_capacity_op)
                / sum(prop_in.flow_equiv_phase_comp["Liq", j] * b.removal_efficiency[j] for j in solutes),
                to_units=pyunits.dimensionless,
            )

        @self.Constraint()
        def eq_breakthrough_time(b):
            return b.breakthrough_time == pyunits.convert(
                (b.bed_depth * b.bv) / b.loading_rate, to_units=pyunits.s
            )

        # @self.Constraint()
        # def eq_bed_vol_tot(b):
        #     return b.bed_volume_total == b.bed_volume * b.number_columns

        @self.Constraint(solutes)
        def eq_mass_transfer_term(b, j):
            return (
                b.process_flow.mass_transfer_term[0, "Liq", j]
                == -1 * prop_in.flow_mass_phase_comp["Liq", j] * b.removal_efficiency[j]
            )

        @self.Constraint()
        def eq_bed_vol(b):
            return b.bed_volume_total == pyunits.convert(
                prop_in.flow_vol_phase["Liq"] * b.ebct, to_units=pyunits.m**3
            )
