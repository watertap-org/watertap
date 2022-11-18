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

# Import Pyomo libraries
from pyomo.environ import (
    Var,
    NonNegativeReals,
    NegativeReals,
    value,
)

from idaes.core import declare_process_block_class, FlowDirection, useDefault
from idaes.core.util import scaling as iscale
from watertap.core import (
    MembraneChannel1DBlock,
    MassTransferCoefficient,
    PressureChangeType,
)
from watertap.core.membrane_channel1d import CONFIG_Template
from watertap.unit_models.osmotically_assisted_reverse_osmosis_base import (
    OsmoticallyAssistedReverseOsmosisBaseData,
    _add_has_full_reporting,
)
import idaes.logger as idaeslog


__author__ = "Adam Atia, Chenyu Wang"

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("OsmoticallyAssistedReverseOsmosis1D")
class OsmoticallyAssistedReverseOsmosis1DData(
    OsmoticallyAssistedReverseOsmosisBaseData
):
    """
    Standard 1D OARO Unit Model Class:
    - one dimensional model
    - steady state only
    - single liquid phase only
    """

    CONFIG = CONFIG_Template()

    _add_has_full_reporting(CONFIG)

    def _process_config(self):
        if self.config.transformation_method is useDefault:
            _log.warning(
                "Discretization method was "
                "not specified for the "
                "reverse osmosis module. "
                "Defaulting to finite "
                "difference method."
            )
            self.config.transformation_method = "dae.finite_difference"

        if self.config.transformation_scheme is useDefault:
            _log.warning(
                "Discretization scheme was "
                "not specified for the "
                "reverse osmosis module."
                "Defaulting to backward finite "
                "difference."
            )
            self.config.transformation_scheme = "BACKWARD"

    def _add_membrane_channel_and_geometry(
        self, side="feed_side", flow_direction=FlowDirection.forward
    ):
        # Check configuration errors
        self._process_config()

        if not isinstance(side, str):
            raise TypeError(
                f"{side} is not a string. Please provide a string for the side argument."
            )

        # Build membrane channel control volume
        setattr(
            self,
            side,
            MembraneChannel1DBlock(
                dynamic=self.config.dynamic,
                has_holdup=self.config.has_holdup,
                area_definition=self.config.area_definition,
                property_package=self.config.property_package,
                property_package_args=self.config.property_package_args,
                transformation_method=self.config.transformation_method,
                transformation_scheme=self.config.transformation_scheme,
                finite_elements=self.config.finite_elements,
                collocation_points=self.config.collocation_points,
            ),
        )
        mem_side = getattr(self, side)

        # if (self.config.pressure_change_type != PressureChangeType.fixed_per_stage) or (
        #     self.config.mass_transfer_coefficient == MassTransferCoefficient.calculated
        # ):
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

    def _add_deltaP(self, side="feed_side"):
        if not isinstance(side, str):
            raise TypeError(
                f"{side} is not a string. Please provide a string for the side argument."
            )

        units_meta = self.config.property_package.get_metadata().get_derived_units
        mem_side = getattr(self, side)
        setattr(
            mem_side,
            "deltaP_stage",
            Var(
                self.flowsheet().config.time,
                initialize=-1e5,
                bounds=(-1e6, 0),
                domain=NegativeReals,
                units=units_meta("pressure"),
                doc=f"Pressure drop across {side}",
            ),
        )
        if self.config.pressure_change_type == PressureChangeType.fixed_per_stage:

            @mem_side.Constraint(
                self.flowsheet().config.time,
                self.length_domain,
                doc="Fixed pressure drop across unit",
            )
            def eq_pressure_drop(b, t, x):
                return mem_side.deltaP_stage[t] == b.length * mem_side.dP_dx[t, x]

        else:

            @mem_side.Constraint(
                self.flowsheet().config.time, doc="Pressure drop across unit"
            )
            def eq_pressure_drop(b, t):
                return mem_side.deltaP_stage[t] == sum(
                    mem_side.dP_dx[t, x] * b.length / b.nfe
                    for x in b.difference_elements
                )

    def _add_mass_transfer(self):

        units_meta = self.config.property_package.get_metadata().get_derived_units

        # mass transfer
        def mass_transfer_phase_comp_initialize(b, t, x, p, j):
            return value(
                self.feed_side.properties[t, x].get_material_flow_terms("Liq", j)
                * self.recovery_mass_phase_comp[t, "Liq", j]
            )

        self.mass_transfer_phase_comp = Var(
            self.flowsheet().config.time,
            self.difference_elements,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            initialize=mass_transfer_phase_comp_initialize,
            bounds=(0.0, 1e6),
            domain=NonNegativeReals,
            units=units_meta("mass")
            * units_meta("time") ** -1
            * units_meta("length") ** -1,
            doc="Mass transfer to permeate",
        )

        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Mass transfer term",
        )
        def eq_mass_transfer_term(self, t, x, p, j):
            return (
                self.mass_transfer_phase_comp[t, x, p, j]
                == -self.feed_side.mass_transfer_term[t, x, p, j]
            )

        # Feed and permeate-side connection
        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Mass transfer from feed to permeate",
        )
        def eq_connect_mass_transfer(b, t, x, p, j):
            return (
                b.permeate_side.mass_transfer_term[t, x, p, j]
                == -b.feed_side.mass_transfer_term[t, x, p, j]
            )

        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Permeate production",
        )
        def eq_permeate_production(b, t, x, p, j):
            return (
                -b.feed_side.mass_transfer_term[t, x, p, j]
                == b.width * b.flux_mass_phase_comp[t, x, p, j]
            )

    def calculate_scaling_factors(self):
        if iscale.get_scaling_factor(self.dens_solvent) is None:
            sf = iscale.get_scaling_factor(
                self.feed_side.properties[0, 0].dens_mass_phase["Liq"]
            )
            iscale.set_scaling_factor(self.dens_solvent, sf)

        super().calculate_scaling_factors()

        for (t, x, p, j), v in self.mass_transfer_phase_comp.items():
            sf = iscale.get_scaling_factor(
                self.feed_side.properties[t, x].get_material_flow_terms(p, j)
            )
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, sf)
            v = self.feed_side.mass_transfer_term[t, x, p, j]
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, sf)
            v = self.permeate_side.mass_transfer_term[t, x, p, j]
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, sf)

        if hasattr(self, "length"):
            if iscale.get_scaling_factor(self.length) is None:
                iscale.set_scaling_factor(self.length, 1)

        if hasattr(self, "width"):
            if iscale.get_scaling_factor(self.width) is None:
                iscale.set_scaling_factor(self.width, 1)
