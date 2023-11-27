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
from pyomo.environ import (
    Var,
    NonNegativeReals,
    value,
)

from idaes.core import declare_process_block_class, FlowDirection
from idaes.core.util import scaling as iscale
from idaes.core.util.misc import add_object_reference
from watertap.core import (
    MembraneChannel0DBlock,
    MassTransferCoefficient,
    PressureChangeType,
)
from watertap.core.membrane_channel0d import CONFIG_Template
from watertap.unit_models.osmotically_assisted_reverse_osmosis_base import (
    OsmoticallyAssistedReverseOsmosisBaseData,
    _add_has_full_reporting,
)


__author__ = "Adam Atia, Chenyu Wang"


@declare_process_block_class("OsmoticallyAssistedReverseOsmosis0D")
class OsmoticallyAssistedReverseOsmosisData(OsmoticallyAssistedReverseOsmosisBaseData):
    """
    Standard OARO Unit Model Class:
    - zero dimensional model
    - steady state only
    - single liquid phase only
    """

    CONFIG = CONFIG_Template()

    _add_has_full_reporting(CONFIG)

    def _add_membrane_channels_and_geometry(self):

        # Build membrane channel control volume
        channel_kwargs = {
            "dynamic": False,
            "has_holdup": False,
            "property_package": self.config.property_package,
            "property_package_args": self.config.property_package_args,
        }
        self.feed_side = MembraneChannel0DBlock(**channel_kwargs)
        self.permeate_side = MembraneChannel0DBlock(**channel_kwargs)

        if (self.config.pressure_change_type != PressureChangeType.fixed_per_stage) or (
            self.config.mass_transfer_coefficient == MassTransferCoefficient.calculated
        ):
            self._add_length_and_width()
            add_geometry_kwargs = {
                "length_var": self.length,
                "width_var": self.width,
            }
            self.feed_side.add_geometry(
                flow_direction=FlowDirection.forward, **add_geometry_kwargs
            )
            self.permeate_side.add_geometry(
                flow_direction=FlowDirection.backward, **add_geometry_kwargs
            )
            self._add_area(include_constraint=True)
        else:
            add_geometry_kwargs = {
                "length_var": None,
                "width_var": None,
            }
            self.feed_side.add_geometry(
                flow_direction=FlowDirection.forward, **add_geometry_kwargs
            )
            self.permeate_side.add_geometry(
                flow_direction=FlowDirection.backward, **add_geometry_kwargs
            )
            self._add_area(include_constraint=False)

    def _add_deltaP(self, side):
        mem_side = self.component(side)
        add_object_reference(mem_side, "deltaP_stage", mem_side.deltaP)

    def _add_mass_transfer(self):

        units_meta = self.config.property_package.get_metadata().get_derived_units

        # mass transfer
        def mass_transfer_phase_comp_initialize(b, t, p, j):
            return value(
                self.feed_side.properties_in[t].get_material_flow_terms("Liq", j)
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
            doc="Mass transfer to permeate",
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
                == -self.feed_side.mass_transfer_term[t, p, j]
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
                b.permeate_side.mass_transfer_term[t, p, j]
                == -b.feed_side.mass_transfer_term[t, p, j]
            )

        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Permeate production",
        )
        def eq_permeate_production(b, t, p, j):
            return (
                b.permeate_side.mass_transfer_term[t, p, j]
                == b.area * b.flux_mass_phase_comp_avg[t, p, j]
            )

    def calculate_scaling_factors(self):

        super().calculate_scaling_factors()

        for (t, p, j), v in self.mass_transfer_phase_comp.items():
            sf = iscale.get_scaling_factor(
                self.feed_side.properties_in[t].get_material_flow_terms(p, j)
            )
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, sf)
            v = self.feed_side.mass_transfer_term[t, p, j]
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, sf)
            v = self.permeate_side.mass_transfer_term[t, p, j]
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, sf)
