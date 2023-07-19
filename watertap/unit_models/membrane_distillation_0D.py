#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
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
from pyomo.environ import Var, NonNegativeReals, value, Constraint

from pyomo.common.config import Bool, ConfigDict, ConfigValue, In, ConfigBlock

from idaes.core import declare_process_block_class, FlowDirection
from idaes.core.util import scaling as iscale
from watertap.core import (
    MembraneChannel0DBlock,
    MassTransferCoefficient,
    PressureChangeType,
)
from watertap.core.temperature_polarization_mixn import CONFIG_Template
from watertap.unit_models.membrane_distillation_base import (
    MembraneDistillationBaseData,
    _add_has_full_reporting,
)

__author__ = "Elmira Shamlou"





@declare_process_block_class("MembraneDistillation0D")
class MembraneDistillationData(MembraneDistillationBaseData):
    """
    Standard DCMD Unit Model Class:
    - zero dimensional model
    - steady state only
    """

    CONFIG = ConfigBlock()


    def create_config_block(CONFIG_Template):
        config_block = ConfigBlock(implicit=True)
        config_block.update(CONFIG_Template)
        return config_block


    CONFIG.declare(
        "hot_side",
        create_config_block(CONFIG_Template),
        description="Config block for feed side",
        doc="""A config block used to construct the feed side control volume.
        This config can be given by the feed side name instead of feed_side.""",
    )

    CONFIG.declare(
        "cold_side",
        create_config_block(CONFIG_Template),
        description="Config block for permeate side",
        doc="""A config block used to construct the permeate side control volume.
        This config can be given by the permeate side name instead of permeate_side.""",
    )

    def _add_membrane_channel_and_geometry(
        self, config, channel="evaporator", flow_direction=FlowDirection.forward
    ):

        # Build membrane channel control volume
        setattr(
            self,
            channel,
            MembraneChannel0DBlock(
                dynamic=False,
                has_holdup=False,
                property_package=config.property_package,
                property_package_args=self.config.property_package_args,
            ),
        )
        mem_side = getattr(self, channel)

        if (self.config.pressure_change_type != PressureChangeType.fixed_per_stage) or (
            self.config.mass_transfer_coefficient == MassTransferCoefficient.calculated
        ):
            if not hasattr(self, "length") and not hasattr(self, "width"):
                self._add_length_and_width()
            mem_side.add_geometry(
                length_var=self.length,
                width_var=self.width,
                flow_direction=flow_direction,
            )
            if not hasattr(self, "eq_area"):
                add_eq_area = True
            else:
                add_eq_area = False
            self._add_area(include_constraint=add_eq_area)
        else:
            mem_side.add_geometry(
                length_var=None, width_var=None, flow_direction=flow_direction
            )
            self._add_area(include_constraint=False)

            #############flux and mass transfer section

    def _add_mass_transfer(self):

        units_meta = self.config.property_package.get_metadata().get_derived_units

        # mass transfer

        def mass_transfer_phase_comp_initialize(b, t, p, j):
            return value(
                self.hot_side.properties_in[t].get_material_flow_terms("Liq", j)
                * self.recovery_mass_phase_comp[t, "Liq", j]
            )

        self.mass_transfer_phase_comp = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            initialize=mass_transfer_phase_comp_initialize,
            bounds=(0.0, 1e6),
            domain=NonNegativeReals,
            units=units_meta("mass") * units_meta("time") ** -1,
            doc="Mass transfer to cold side",
        )

        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Mass transfer term",
        )
        def eq_mass_transfer_term(self, t, p, j):
            return (
                self.mass_transfer_phase_comp[t, p, j]
                == -self.hot_side.mass_transfer_term[t, p, j]
            )

        # Feed and permeate-side connection
        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Mass transfer from feed to permeate",
        )
        def eq_connect_mass_transfer(b, t, p, j):
            return (
                b.cold_side.mass_transfer_term[t, p, j]
                == -b.hot_side.mass_transfer_term[t, p, j]
            )

        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Permeate production",
        )
        def eq_permeate_production(b, t, p, j):
            if j == "H2O":
                return b.cold_side.mass_transfer_term[
                    t, p, j
                ] == b.area * b.flux_mass_avg(b, t)
            else:
                b.cold_side.mass_transfer_term[t, p, j].fix(
                    0
                )  # no mass transfer for components other than water
                return (
                    Constraint.Skip
                )  # skip the constraint for components other than water

    def _add_heat_transfer(self):

        # Conductive heat transfer for the hot side
        @self.Constraint(
            self.flowsheet().config.time,
            doc="Conductive heat transfer from hot side",
        )
        def eq_conductive_heat_transfer_hot(b, t):
            return b.hot_side.heat[t] == -b.area * b.flux_conduction_heat_avg[t]

        # hot side and cold side conduction heat connection
        @self.Constraint(
            self.flowsheet().config.time,
            doc="Conductive heat transfer to cold side",
        )
        def eq_conductive_heat_transfer_cold(b, t):
            return b.cold_side.heat[t] == -b.hot_side.heat[t]

        # Enthalpy transfer for the hot side
        @self.Constraint(
            self.flowsheet().config.time,
            doc="Enthalpy heat transfer from the hot side",
        )
        def eq_enthalpy_transfer_hot(b, t):
            return b.hot_side.enthalpy_transfer[t] == -b.area * b.flux_enth_hot_avg[t]

        # Enthalpy heat transfer for the cold side
        @self.Constraint(
            self.flowsheet().config.time,
            doc="Enthalpy heat transfer to the cold side",
        )
        def eq_enthalpy_transfer_cold(b, t):
            return b.cold_side.enthalpy_transfer[t] == b.area * b.flux_enth_cold_avg[t]

    def calculate_scaling_factors(self):
        # First, we'll call the `calculate_scaling_factors` method of the super class
        super().calculate_scaling_factors()

        # Next, we'll define the scaling factors for the variables and constraints
        # that we've defined in this class.

        for (t, p, j), v in self.mass_transfer_phase_comp.items():
            sf = iscale.get_scaling_factor(
                self.hot_side.properties_in[t].get_material_flow_terms(p, j)
            )
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, sf)
            v = self.hot_side.mass_transfer_term[t, p, j]
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, sf)
            v = self.cold_side.mass_transfer_term[t, p, j]
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, sf)

        # Get scaling factor from vapor enthalpy flow
        sf_vap = iscale.get_scaling_factor(
            self.properties_vapor[0].enth_flow_phase["Vap"], default=1
        )

        # Scaling factors for heat transfer
        for t in self.flowsheet().config.time:
            # If scaling factor has not been set for heat transfer, set it to sf_vap
            if iscale.get_scaling_factor(self.hot_side.heat[t]) is None:
                iscale.set_scaling_factor(self.hot_side.heat[t], sf_vap)
            if iscale.get_scaling_factor(self.hot_side.enthalpy_transfer[t]) is None:
                iscale.set_scaling_factor(self.hot_side.enthalpy_transfer[t], sf_vap)
            if iscale.get_scaling_factor(self.cold_side.enthalpy_transfer[t]) is None:
                iscale.set_scaling_factor(self.cold_side.enthalpy_transfer[t], sf_vap)

        # Scaling factors for area, length and width if they are defined
        if hasattr(self, "area"):
            if iscale.get_scaling_factor(self.area) is None:
                iscale.set_scaling_factor(self.area, 1)

        if hasattr(self, "length"):
            if iscale.get_scaling_factor(self.length) is None:
                iscale.set_scaling_factor(self.length, 1)

        if hasattr(self, "width"):
            if iscale.get_scaling_factor(self.width) is None:
                iscale.set_scaling_factor(self.width, 1)
