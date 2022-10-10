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
This module contains a zero-order representation of a gas-sparged membrane unit.
"""
from idaes.core import declare_process_block_class

from watertap.core import pump_electricity, ZeroOrderBaseData
import idaes.logger as idaeslog
from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale
from pyomo.environ import NonNegativeReals, Var, units as pyunits, Reference
from watertap.core.zero_order_sido import initialize_sido

# Some more information about this module
__author__ = "Adam Atia"

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("GasSpargedMembraneZO")
class GasSpargedMembraneZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a gas-sparged membrane.
    This unit is similar to a SIDO, but there is technically a third outlet for gas extraction.
    Three StateBlocks are added with corresponding Ports:

        * properties_inlet
        * properties_treated
        * properties_byproduct

    Two additional variables are added:

        * recovery_vol (indexed by time)
        * removal_frac_mass_comp (indexed by time and component)

    Four additional constraints are added to represent the material balances, with modifications
    to account for gas extraction.

        * water_recovery_equation (indexed by time)
        * flow_balance (indexed by time)
        * solute_removal_equation (indexed by time and solute)
        * solute_treated_equation (indexed by time and solute)

    The build method also sets private attributes on the unit model with references
    to the appropriate initialization and scaling methods to use and to return
    the inlet volumetric flow rate.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "gas_sparged_membrane"

        self._has_recovery_removal = True
        self._initialize = initialize_sido
        self._scaling = calculate_scaling_factors_gas_extraction

        # Create state blocks for inlet and outlets
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["defined_state"] = True

        self.properties_in = self.config.property_package.build_state_block(
            self.flowsheet().time, doc="Material properties at inlet", **tmp_dict
        )

        tmp_dict_2 = dict(**tmp_dict)
        tmp_dict_2["defined_state"] = False

        self.properties_treated = self.config.property_package.build_state_block(
            self.flowsheet().time,
            doc="Material properties of treated water",
            **tmp_dict_2
        )
        self.properties_byproduct = self.config.property_package.build_state_block(
            self.flowsheet().time,
            doc="Material properties of byproduct stream",
            **tmp_dict_2
        )

        # Create Ports
        self.add_port("inlet", self.properties_in, doc="Inlet port")
        self.add_port(
            "treated", self.properties_treated, doc="Treated water outlet port"
        )
        self.add_port(
            "byproduct", self.properties_byproduct, doc="Byproduct outlet port"
        )

        # Add performance variables
        self.recovery_frac_mass_H2O = Var(
            self.flowsheet().time,
            initialize=0.8,
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            bounds=(0.0, 1.0000001),
            doc="Mass recovery fraction of water in the treated stream",
        )
        self.removal_frac_mass_comp = Var(
            self.flowsheet().time,
            self.config.property_package.solute_set,
            domain=NonNegativeReals,
            initialize=0.01,
            units=pyunits.dimensionless,
            doc="Solute removal fraction on a mass basis",
        )
        self.gas_mass_influent_ratio = Var(
            self.flowsheet().time,
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Mass flow of gas extracted per mass flow of influent",
        )
        self.flow_mass_gas_extraction = Var(
            self.flowsheet().time,
            domain=NonNegativeReals,
            units=pyunits.kg / pyunits.s,
            doc="Mass flow of hydrogen extracted",
        )
        self._fixed_perf_vars.append(self.gas_mass_influent_ratio)
        self._perf_var_dict[
            "Mass of gas extracted per mass flow of influent(kg/d/(kg/d)"
        ] = self.gas_mass_influent_ratio
        self._perf_var_dict[
            "Mass flow of gas extracted (kg/s))"
        ] = self.flow_mass_gas_extraction

        # Add performance constraints
        # Water recovery
        @self.Constraint(self.flowsheet().time, doc="Water recovery equation")
        def water_recovery_equation(b, t):
            return (
                b.recovery_frac_mass_H2O[t] * b.properties_in[t].flow_mass_comp["H2O"]
                == b.properties_treated[t].flow_mass_comp["H2O"]
            )

        # Flow balance
        @self.Constraint(self.flowsheet().time, doc="Overall flow balance")
        def mass_balance(b, t):
            return (
                sum(
                    b.properties_in[t].flow_mass_comp[j]
                    for j in self.config.property_package.component_list
                )
                == sum(
                    b.properties_treated[t].flow_mass_comp[j]
                    for j in self.config.property_package.component_list
                )
                + sum(
                    b.properties_byproduct[t].flow_mass_comp[j]
                    for j in self.config.property_package.component_list
                )
                + b.flow_mass_gas_extraction[t]
            )

        # Gas extraction
        @self.Constraint(self.flowsheet().time, doc="Gas extraction equation")
        def mass_gas_extraction_equation(b, t):
            return b.flow_mass_gas_extraction[t] == b.gas_mass_influent_ratio[t] * sum(
                b.properties_in[t].flow_mass_comp[j]
                for j in self.config.property_package.component_list
            )

        # Solute removal
        @self.Constraint(
            self.flowsheet().time,
            self.config.property_package.solute_set,
            doc="Solute removal equations",
        )
        def solute_removal_equation(b, t, j):
            return (
                b.removal_frac_mass_comp[t, j] * b.properties_in[t].flow_mass_comp[j]
                == b.properties_byproduct[t].flow_mass_comp[j]
            )

        # Solute concentration of treated stream
        @self.Constraint(
            self.flowsheet().time,
            self.config.property_package.solute_set,
            doc="Constraint for solute concentration in treated " "stream.",
        )
        def solute_treated_equation(b, t, j):
            return (1 - b.removal_frac_mass_comp[t, j]) * b.properties_in[
                t
            ].flow_mass_comp[j] == b.properties_treated[t].flow_mass_comp[j]

        self._stream_table_dict = {
            "Inlet": self.inlet,
            "Treated": self.treated,
            "Byproduct": self.byproduct,
        }

        self._perf_var_dict["Water Recovery"] = self.recovery_frac_mass_H2O
        self._perf_var_dict["Solute Removal"] = self.removal_frac_mass_comp

        self._get_Q = _get_Q_gas_extraction

        self._Q = Reference(self.properties_in[:].flow_vol)
        pump_electricity(self, self._Q)


def calculate_scaling_factors_gas_extraction(self):
    # Get default scale factors and do calculations from base classes
    for t, v in self.water_recovery_equation.items():
        iscale.constraint_scaling_transform(
            v,
            iscale.get_scaling_factor(
                self.properties_in[t].flow_mass_comp["H2O"],
                default=1,
                warning=True,
                hint=" for water recovery",
            ),
        )

    for t, v in self.mass_balance.items():
        iscale.constraint_scaling_transform(
            v,
            iscale.get_scaling_factor(
                self.properties_in[t].flow_mass_comp["H2O"], default=1, warning=False
            ),
        )  # would just be a duplicate of above

    for (t, j), v in self.solute_removal_equation.items():
        iscale.constraint_scaling_transform(
            v,
            iscale.get_scaling_factor(
                self.properties_in[t].flow_mass_comp[j],
                default=1,
                warning=True,
                hint=" for solute removal",
            ),
        )

    for (t, j), v in self.solute_treated_equation.items():
        iscale.constraint_scaling_transform(
            v,
            iscale.get_scaling_factor(
                self.properties_in[t].flow_mass_comp[j], default=1, warning=False
            ),
        )  # would just be a duplicate of above
    for t, v in self.mass_gas_extraction_equation.items():
        iscale.constraint_scaling_transform(v, 1e3)


def _get_Q_gas_extraction(self, t):
    return self.properties_in[t].flow_vol
