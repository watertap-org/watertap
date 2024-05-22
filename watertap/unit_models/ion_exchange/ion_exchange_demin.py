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

        # typical operating capacity for resin is between 30-60 kg/m3; Wachinski, chap. 3
        self.resin_capacity = Var(
            initialize=42,
            bounds=(
                25,
                65,
            ),
            units=pyunits.kg / pyunits.m**3,
            doc="Operational capacity of the resin in kg/m3 as CaCO3",
        )

        @self.Expression(doc="Mass of target to be removed per cycle")
        def mass_removed_per_cycle(b):
            return pyunits.convert(
                prop_out.flow_mass_phase_comp["Liq", target] * b.breakthrough_time,
                to_units=pyunits.kg,
            )

        @self.Constraint(doc="Effluent concentration")
        def eq_effluent(b):
            return (
                prop_out.conc_mass_phase_comp["Liq", target]
                == b.c_norm[target] * prop_in.conc_mass_phase_comp["Liq", target]
            )

        # @self.Constraint(doc="Regeneration stream")

        @self.Constraint(doc="Total bed volume required")
        def bed_volume_total_constraint(b):
            return b.bed_volume_total == pyunits.convert(
                b.mass_removed_per_cycle / b.resin_capacity, to_units=pyunits.m**3
            )
