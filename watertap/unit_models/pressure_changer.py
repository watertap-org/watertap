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

from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.environ import Var, units as pyunits, Expr_if, value

from enum import Enum, auto

# Import IDAES cores
from idaes.models.unit_models.pressure_changer import PumpData
from idaes.core import declare_process_block_class
from idaes.core.util.constants import Constants
import idaes.core.util.scaling as iscale

import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


class VariableEfficiency(Enum):
    none = auto()  # default is constant efficiency
    flow = auto()  # flow-only correlation
    flow_head = auto()  # flow and head correlation


@declare_process_block_class("Pump")
class PumpIsothermalData(PumpData):
    """
    Standard Isothermal Pump Unit Model Class
    """

    CONFIG = PumpData.CONFIG()
    CONFIG.declare(
        "variable_efficiency",
        ConfigValue(
            default=VariableEfficiency.none,
            domain=In(VariableEfficiency),
            description="Variable pump efficiency flag",
            doc="""Indicates the relationship used to define pump efficiency
    **VariableEfficiency.none** - uses default pump efficiency at BEP
    **VariableEfficiency.flow** - uses an efficiency correlation scaled to the BEP flow rate
    **VariableEfficiency.flow_head** - uses an efficiency correlation scaled to the BEP flow rate and head
    """,
        ),
    )

    def build(self):
        super().build()

        # ---------------------------------------
        # Isothermal pump set up
        # ---------------------------------------
        if hasattr(self.control_volume, "enthalpy_balances"):
            self.control_volume.del_component(self.control_volume.enthalpy_balances)

        @self.control_volume.Constraint(
            self.flowsheet().config.time, doc="Isothermal constraint"
        )
        def isothermal_balance(b, t):
            return b.properties_in[t].temperature == b.properties_out[t].temperature

        # ---------------------------------------
        # Variable efficiency pump set-up
        # ---------------------------------------

        if self.config.variable_efficiency is not VariableEfficiency.none:
            # create additional pyomo variables
            self.bep_flow = Var(
                initialize=1.0,
                doc="Best efficiency point flowrate of the centrifugal pump",
                units=pyunits.m**3 / pyunits.s,
            )

            self.bep_eta = Var(
                initialize=0.8,
                doc="Best efficiency of the centrifugal pump",
                units=pyunits.dimensionless,
            )

            self.flow_ratio = Var(
                self.flowsheet().time,
                initialize=1.0,
                doc="Ratio of pump flowrate to best efficiency point flowrate",
                units=pyunits.dimensionless,
            )

            # add constraints
            @self.Constraint(self.flowsheet().time, doc="Pump flow ratio")
            def flow_ratio_constraint(b, t):
                return (
                    b.flow_ratio[t] * b.bep_flow
                    == b.control_volume.properties_in[t].flow_vol
                )

        if self.config.variable_efficiency is VariableEfficiency.flow:

            @self.Expression(
                self.flowsheet().time,
                doc="Expression for variable pump efficiency based on flow only",
            )
            def eta_ratio(b, t):
                return Expr_if(
                    b.flow_ratio[t] < 0.6,
                    0.4,
                    Expr_if(
                        b.flow_ratio[t] > 1.4,
                        0.4,
                        -0.995 * b.flow_ratio[t] ** 2 + 1.977 * b.flow_ratio[t] + 0.018,
                    ),
                )

        elif self.config.variable_efficiency is VariableEfficiency.flow_head:
            raise NotImplementedError(
                "Config option 'VariableEfficiency.flow_head' is not fully implemented yet"
            )
            # TODO - Implement pump efficiency expression based on flow and head (bep_head, head_ratio)
        else:
            pass

        if self.config.variable_efficiency is not VariableEfficiency.none:
            # replace the constant efficiency assumption using eta_ratio
            # must be done after the eta_ratio expression is created.
            @self.Constraint(self.flowsheet().time, doc="Actual pump efficiency")
            def eta_constraint(b, t):
                return b.efficiency_pump[t] == (b.bep_eta * b.eta_ratio[t])

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        for ind, c in self.control_volume.isothermal_balance.items():
            sf = iscale.get_scaling_factor(
                self.control_volume.properties_in[0].temperature
            )
            iscale.constraint_scaling_transform(c, sf)

        if hasattr(self, "bep_flow"):
            if iscale.get_scaling_factor(self.bep_flow) is None:
                sf = value(self.bep_flow) ** -1
                iscale.set_scaling_factor(self.bep_flow, sf)

        if hasattr(self, "bep_head"):
            if iscale.get_scaling_factor(self.bep_head) is None:
                sf = value(self.bep_head) ** -1
                iscale.set_scaling_factor(self.bep_head, sf)

        if hasattr(self, "bep_eta"):
            if iscale.get_scaling_factor(self.bep_eta) is None:
                iscale.set_scaling_factor(self.bep_eta, 1)

        for t in self.flowsheet().time:
            if hasattr(self, "flow_ratio"):
                if iscale.get_scaling_factor(self.flow_ratio[t]) is None:
                    iscale.set_scaling_factor(self.flow_ratio[t], 1)

            if hasattr(self, "efficiency_pump"):
                if iscale.get_scaling_factor(self.efficiency_pump[t]) is None:
                    iscale.set_scaling_factor(self.efficiency_pump[t], 1)

            # scale constraints

            if hasattr(self, "flow_ratio_constraint"):
                if iscale.get_scaling_factor(self.flow_ratio_constraint[t]) is None:
                    iscale.set_scaling_factor(self.flow_ratio_constraint[t], 1)

            if hasattr(self, "eta_constraint"):
                if iscale.get_scaling_factor(self.eta_constraint[t]) is None:
                    iscale.set_scaling_factor(self.eta_constraint[t], 1)


@declare_process_block_class("EnergyRecoveryDevice")
class EnergyRecoveryDeviceData(PumpIsothermalData):
    """
    Turbine-type isothermal energy recovery device
    """

    # switch compressor to False
    CONFIG = PumpIsothermalData.CONFIG()
    CONFIG.get("compressor")._default = False
    CONFIG.get("compressor")._domain = In([False])
    CONFIG.compressor = False
