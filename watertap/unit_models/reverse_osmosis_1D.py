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
    check_optimal_termination,
    Set,
)
from pyomo.common.config import ConfigValue, In

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    useDefault,
)
from idaes.core.util import scaling as iscale
from idaes.core.util.misc import add_object_reference
import idaes.logger as idaeslog

from watertap.core import (
    MembraneChannel1DBlock,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
from watertap.core.membrane_channel1d import CONFIG_Template
from watertap.unit_models.reverse_osmosis_base import (
    ReverseOsmosisBaseData,
    _add_has_full_reporting,
    _add_object_reference_if_exists,
)

__author__ = "Adam Atia"

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("ReverseOsmosis1D")
class ReverseOsmosis1DData(ReverseOsmosisBaseData):
    """Standard 1D Reverse Osmosis Unit Model Class."""

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

    def _add_feed_side_membrane_channel_and_geometry(self):
        # Check configuration errors
        self._process_config()

        # Build 1D Membrane Channel
        self.feed_side = MembraneChannel1DBlock(
            dynamic=self.config.dynamic,
            has_holdup=self.config.has_holdup,
            area_definition=self.config.area_definition,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
            transformation_method=self.config.transformation_method,
            transformation_scheme=self.config.transformation_scheme,
            finite_elements=self.config.finite_elements,
            collocation_points=self.config.collocation_points,
        )
        self._add_length_and_width()
        self.feed_side.add_geometry(length_var=self.length, width_var=self.width)
        self._add_area(include_constraint=True)

    def _add_deltaP(self):
        units_meta = self.config.property_package.get_metadata().get_derived_units
        self.deltaP = Var(
            self.flowsheet().config.time,
            initialize=-1e5,
            bounds=(-1e6, 0),
            domain=NegativeReals,
            units=units_meta("pressure"),
            doc="Pressure drop across unit",
        )
        if self.config.pressure_change_type == PressureChangeType.fixed_per_stage:

            @self.Constraint(
                self.flowsheet().config.time,
                self.length_domain,
                doc="Fixed pressure drop across unit",
            )
            def eq_pressure_drop(b, t, x):
                return b.deltaP[t] == b.length * b.dP_dx[t, x]

        else:

            @self.Constraint(
                self.flowsheet().config.time, doc="Pressure drop across unit"
            )
            def eq_pressure_drop(b, t):
                return b.deltaP[t] == sum(
                    b.dP_dx[t, x] * b.length / b.nfe for x in b.difference_elements
                )

    def _add_mass_transfer(self):

        units_meta = self.config.property_package.get_metadata().get_derived_units

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

        # ==========================================================================
        # Mass transfer term equation
        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Mass transfer term",
        )
        def eq_mass_transfer_term(b, t, x, p, j):
            return (
                b.mass_transfer_phase_comp[t, x, p, j]
                == -b.feed_side.mass_transfer_term[t, x, p, j]
            )

        # ==========================================================================
        # Feed and permeate-side mass transfer connection --> Mp,j = Mf,transfer = Jj * W * L/n
        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Mass transfer from feed to permeate",
        )
        def eq_connect_mass_transfer(b, t, x, p, j):
            return (
                b.permeate_side[t, x].get_material_flow_terms(p, j)
                == -b.feed_side.mass_transfer_term[t, x, p, j] * b.length / b.nfe
            )

        # ==========================================================================
        # Mass flux = feed mass transfer equation
        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Mass transfer term",
        )
        def eq_mass_flux_equal_mass_transfer(b, t, x, p, j):
            return (
                b.flux_mass_phase_comp[t, x, p, j] * b.width
                == -b.feed_side.mass_transfer_term[t, x, p, j]
            )

        # ==========================================================================
        # Final permeate mass flow rate (of solvent and solute) --> Mp,j, final = sum(Mp,j)

        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Permeate mass flow rates exiting unit",
        )
        def eq_permeate_production(b, t, p, j):
            return b.mixed_permeate[t].get_material_flow_terms(p, j) == sum(
                b.permeate_side[t, x].get_material_flow_terms(p, j)
                for x in b.difference_elements
            )

    def calculate_scaling_factors(self):
        if iscale.get_scaling_factor(self.dens_solvent) is None:
            sf = iscale.get_scaling_factor(
                self.feed_side.properties[0, 0].dens_mass_phase["Liq"]
            )
            iscale.set_scaling_factor(self.dens_solvent, sf)

        super().calculate_scaling_factors()

        # these variables should have user input, if not there will be a warning
        if iscale.get_scaling_factor(self.width) is None:
            sf = iscale.get_scaling_factor(self.width, default=1, warning=True)
            iscale.set_scaling_factor(self.width, sf)

        if iscale.get_scaling_factor(self.length) is None:
            sf = iscale.get_scaling_factor(self.length, default=10, warning=True)
            iscale.set_scaling_factor(self.length, sf)

        for (t, x, p, j), v in self.mass_transfer_phase_comp.items():
            sf = (
                iscale.get_scaling_factor(
                    self.feed_side.properties[t, x].get_material_flow_terms(p, j)
                )
                / iscale.get_scaling_factor(self.feed_side.length)
            ) * value(self.nfe)
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, sf)
            v = self.feed_side.mass_transfer_term[t, x, p, j]
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, sf)

        if hasattr(self, "deltaP"):
            for v in self.deltaP.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1e-4)
