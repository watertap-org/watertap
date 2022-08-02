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

from pyomo.common.config import In
from pyomo.environ import Var, units as pyunits

# Import IDAES cores
from idaes.models.unit_models.pressure_changer import PumpData
from idaes.core import declare_process_block_class
from idaes.core.util.constants import Constants
import idaes.core.util.scaling as iscale

import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


@declare_process_block_class("Pump")
class PumpIsothermalData(PumpData):
    """
    Standard Isothermal Pump Unit Model Class
    """

    def build(self):
        super().build()

        self.control_volume.del_component(self.control_volume.enthalpy_balances)

        @self.control_volume.Constraint(
            self.flowsheet().config.time, doc="Isothermal constraint"
        )
        def isothermal_balance(b, t):
            return b.properties_in[t].temperature == b.properties_out[t].temperature

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        for ind, c in self.control_volume.isothermal_balance.items():
            sf = iscale.get_scaling_factor(
                self.control_volume.properties_in[0].temperature
            )
            iscale.constraint_scaling_transform(c, sf)


@declare_process_block_class("PumpVar")
class PumpVariableData(PumpData):
    """
    Variable efficiency Pump Unit Model Class
    Uses pump affinity laws for industrial scale centrifugal pumps.
    Atia et al. (2019) https://doi.org/10.1016/j.desal.2019.07.002
    """

    def build(self):
        super().build()

        self.bep_flow = Var(
            initialize=1.0,
            doc="Best efficiency point flowrate of the centrifugal pump",
            units=pyunits.m**3 / pyunits.s,
        )
        self.bep_head = Var(
            initialize=1.0,
            doc="Best efficiency point head of the centrifugal pump",
            units=pyunits.m,
        )
        self.bep_eta = Var(
            initialize=0.85,
            doc="Best efficiency of the centrifugal pump",
            units=pyunits.dimensionless,
        )
        self.head_ratio = Var(
            self.flowsheet().time,
            initialize=1.0,
            doc="Ratio of pump head to best efficiency point head",
            units=pyunits.dimensionless,
        )
        self.flow_ratio = Var(
            self.flowsheet().time,
            initialize=1.0,
            doc="Ratio of pump flowrate to best efficiency point flowrate",
            units=pyunits.dimensionless,
        )
        self.eta_ratio = Var(
            self.flowsheet().time,
            initialize=1.0,
            doc="Pump efficiency change due to variable operation",
            units=pyunits.dimensionless,
        )

        if hasattr(self, "efficiency_pump"):
            self.efficiency_pump.unfix()

        @self.Constraint(self.flowsheet().time, doc="Pump head ratio")
        def head_ratio_constraint(b, t):
            return b.head_ratio[t] * b.bep_head * Constants.acceleration_gravity == (
                (
                    b.control_volume.properties_out[t].pressure
                    / b.control_volume.properties_out[t].density
                )
                - (
                    b.control_volume.properties_in[t].pressure
                    / b.control_volume.properties_in[t].density
                )
            )

        @self.Constraint(self.flowsheet().time, doc="Pump flow ratio")
        def flow_ratio_constraint(b, t):
            return b.flow_ratio[t] * b.bep_flow == (
                b.control_volume.properties_out[t].pressure
            )

        @self.Constraint(self.flowsheet().time, doc="Pump efficiency ratio")
        def eta_ratio_constraint(b, t):
            return b.eta_ratio[t] == (
                (b.flow_ratio[t] ** 1.5) * (b.head_ratio[t] ** 0.75)
            )

        @self.Constraint(self.flowsheet().time, doc="Actual pump efficiency")
        def eta_constraint(b, t):
            return b.efficiency_pump[t] == (b.eta_bep * b.eta_ratio[t])


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
