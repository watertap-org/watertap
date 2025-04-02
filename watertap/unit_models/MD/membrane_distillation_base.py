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


from pyomo.environ import (
    NonNegativeReals,
    Param,
    Var,
    check_optimal_termination,
    exp,
    units as pyunits,
    Constraint,
    log,
)
from enum import Enum, auto

from idaes.core import UnitModelBlockData
from watertap.core.solvers import get_solver
from idaes.core.util import scaling as iscale
from idaes.core.util.exceptions import InitializationError
from idaes.core.util.misc import add_object_reference
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.tables import create_stream_table_dataframe
import idaes.logger as idaeslog
from watertap.core import InitializationMixin

from .MD_channel_base import (
    ConcentrationPolarizationType,
    TemperaturePolarizationType,
    MassTransferCoefficient,
)
from watertap.costing.unit_models.membrane_distillation import (
    cost_membrane_distillation,
)
from watertap.core.util.initialization import interval_initializer


__author__ = "Elmira Shamlou"


class MDconfigurationType(Enum):
    """
    DCMD: Direct Contact Membrane Distillation
    VMD: Vacuum Membrane Distillation
    GMD: Permeate Gap or Conductive Gap Membrane distillation
    AGMD: Air Gap Membrane Distillation
    """

    DCMD = auto()
    VMD = auto()
    GMD = auto()
    AGMD = auto()


class MembraneDistillationBaseData(InitializationMixin, UnitModelBlockData):
    def build(self):
        """
        Common variables and constraints for a DCMD unit model

        """

        super().build()

        self._make_MD_channel_control_volume("hot_ch", self.config, self.config.hot_ch)

        self.hot_ch.add_state_blocks(
            has_phase_equilibrium=False,
        )

        self.hot_ch._add_interface_stateblock(has_phase_equilibrium=False)
        self.hot_ch._add_vapor_stateblock(
            property_package_vapor=self.config.hot_ch.property_package_vapor,
            property_package_args_vapor=self.config.hot_ch.property_package_args_vapor,
            has_phase_equilibrium=False,
        )

        # self.hot_ch.set_config(self.config.hot_ch)

        self.hot_ch.add_material_balances(
            balance_type=self.config.hot_ch.material_balance_type,
            has_mass_transfer=True,
        )

        self.hot_ch.add_momentum_balances(
            balance_type=self.config.hot_ch.momentum_balance_type,
            pressure_change_type=self.config.hot_ch.pressure_change_type,
            has_pressure_change=self.config.hot_ch.has_pressure_change,
            friction_factor=self.config.hot_ch.friction_factor,
        )

        self.hot_ch.add_energy_balances(
            balance_type=self.config.hot_ch.energy_balance_type,
            has_heat_transfer=True,
            has_enthalpy_transfer=True,
        )

        # Add constraint for volumetric flow equality between interface and bulk
        self.hot_ch.add_extensive_flow_to_interface()

        self.hot_ch.add_concentration_polarization(
            concentration_polarization_type=self.config.hot_ch.concentration_polarization_type,
            mass_transfer_coefficient=self.config.hot_ch.mass_transfer_coefficient,
        )

        self.hot_ch.add_temperature_polarization(
            temperature_polarization_type=self.config.hot_ch.temperature_polarization_type,
        )

        try:
            self.hot_ch.apply_transformation()
        except AttributeError:
            pass

        self._make_MD_channel_control_volume(
            "cold_ch", self.config, self.config.cold_ch
        )

        self.cold_ch.add_state_blocks(has_phase_equilibrium=False)

        # cold channel
        # vacuum channel is VMD is also named as cold_ch

        if self.config.MD_configuration_Type in [
            MDconfigurationType.DCMD,
            MDconfigurationType.GMD,
            MDconfigurationType.AGMD,
        ]:
            self.cold_ch._add_interface_stateblock(has_phase_equilibrium=False)

            # Add constraint for volumetric flow equality between interface and bulk
            self.cold_ch.add_extensive_flow_to_interface()

        if self.config.MD_configuration_Type == MDconfigurationType.DCMD:
            self.cold_ch._add_vapor_stateblock(
                property_package_vapor=self.config.cold_ch.property_package_vapor,
                property_package_args_vapor=self.config.cold_ch.property_package_args_vapor,
                has_phase_equilibrium=False,
            )

        if self.config.MD_configuration_Type in [
            MDconfigurationType.DCMD,
            MDconfigurationType.VMD,
        ]:
            self.cold_ch.add_material_balances(
                balance_type=self.config.cold_ch.material_balance_type,
                has_mass_transfer=True,
            )
        elif self.config.MD_configuration_Type in [
            MDconfigurationType.GMD,
            MDconfigurationType.AGMD,
        ]:
            self.cold_ch.add_material_balances(
                balance_type=self.config.cold_ch.material_balance_type,
                has_mass_transfer=False,
            )

        if self.config.MD_configuration_Type in [
            MDconfigurationType.DCMD,
            MDconfigurationType.GMD,
            MDconfigurationType.AGMD,
        ]:
            self.cold_ch.add_momentum_balances(
                balance_type=self.config.cold_ch.momentum_balance_type,
                pressure_change_type=self.config.cold_ch.pressure_change_type,
                has_pressure_change=self.config.cold_ch.has_pressure_change,
                friction_factor=self.config.cold_ch.friction_factor,
            )

        if self.config.MD_configuration_Type == MDconfigurationType.VMD:
            self.cold_ch.add_momentum_balances(
                balance_type=self.config.cold_ch.momentum_balance_type,
                pressure_change_type=self.config.cold_ch.pressure_change_type,
                has_pressure_change=False,
                friction_factor=self.config.cold_ch.friction_factor,
            )

        if self.config.MD_configuration_Type == MDconfigurationType.DCMD:
            self.cold_ch.add_energy_balances(
                balance_type=self.config.cold_ch.energy_balance_type,
                has_heat_transfer=True,
                has_enthalpy_transfer=True,
            )
        elif self.config.MD_configuration_Type in [
            MDconfigurationType.GMD,
            MDconfigurationType.AGMD,
        ]:
            self.cold_ch.add_energy_balances(
                balance_type=self.config.cold_ch.energy_balance_type,
                has_heat_transfer=True,
                has_enthalpy_transfer=False,
            )
        elif self.config.MD_configuration_Type == MDconfigurationType.VMD:
            self.cold_ch.add_energy_balances(
                balance_type=self.config.cold_ch.energy_balance_type,
                has_enthalpy_transfer=True,
            )

        # Concentration polarization constraint is not accounted for in the below method; it is
        # written later in the base model (eq_concentration_polarization)
        if hasattr(self.cold_ch.config.property_package, "solute_set"):
            # If solute_set is defined, add concentration polarization configuration
            self.cold_ch.add_concentration_polarization(
                concentration_polarization_type=ConcentrationPolarizationType.none,
                mass_transfer_coefficient=MassTransferCoefficient.none,
            )
        if self.config.MD_configuration_Type in [
            MDconfigurationType.DCMD,
            MDconfigurationType.GMD,
            MDconfigurationType.AGMD,
        ]:
            self.cold_ch.add_temperature_polarization(
                temperature_polarization_type=self.config.cold_ch.temperature_polarization_type,
            )

        try:
            self.cold_ch.apply_transformation()
        except AttributeError:
            pass

        if self.config.MD_configuration_Type in [
            MDconfigurationType.DCMD,
            MDconfigurationType.GMD,
            MDconfigurationType.AGMD,
        ]:
            for t in self.flowsheet().config.time:
                for x in self.cold_ch.length_domain:

                    if (
                        "Vap" in self.cold_ch.config.property_package.phase_list
                        and "H2O" in self.cold_ch.config.property_package.component_list
                    ):

                        self.cold_ch.properties[t, x].flow_mass_phase_comp[
                            "Vap", "H2O"
                        ].fix(0)
                        self.cold_ch.properties_interface[t, x].flow_mass_phase_comp[
                            "Vap", "H2O"
                        ].fix(0)

        if self.config.MD_configuration_Type == MDconfigurationType.VMD:
            for t in self.flowsheet().config.time:
                for x in self.cold_ch.length_domain:
                    if (
                        "Vap" in self.cold_ch.config.property_package.phase_list
                        and "H2O" in self.cold_ch.config.property_package.component_list
                    ):

                        self.cold_ch.properties[t, x].flow_mass_phase_comp[
                            "Liq", "H2O"
                        ].fix(0)
                        self.cold_ch.properties[
                            t, self.cold_ch.length_domain.first()
                        ].flow_mass_phase_comp["Vap", "H2O"].fix(0)

        # gap channel
        # modification required for AGMD
        if self.config.MD_configuration_Type in [
            MDconfigurationType.GMD,
            MDconfigurationType.AGMD,
        ]:
            self._make_MD_channel_control_volume(
                "gap_ch", self.config, self.config.gap_ch
            )

            self.gap_ch.add_state_blocks(
                has_phase_equilibrium=False,
            )

            self.gap_ch._add_interface_stateblock(has_phase_equilibrium=False)

            self.gap_ch.add_material_balances(
                balance_type=self.config.gap_ch.material_balance_type,
                has_mass_transfer=True,
            )

            self.gap_ch.add_momentum_balances(
                balance_type=self.config.gap_ch.momentum_balance_type,
                pressure_change_type=self.config.gap_ch.pressure_change_type,
                has_pressure_change=self.config.gap_ch.has_pressure_change,
                friction_factor=self.config.gap_ch.friction_factor,
            )

            self.gap_ch.add_energy_balances(
                balance_type=self.config.gap_ch.energy_balance_type,
                has_heat_transfer=False,
                has_enthalpy_transfer=True,
            )

            self.gap_ch.add_extensive_flow_to_interface()

            try:
                self.gap_ch.apply_transformation()
            except AttributeError:
                pass

        if self.config.MD_configuration_Type in [
            MDconfigurationType.GMD,
        ]:

            for t in self.flowsheet().config.time:
                for x in self.gap_ch.length_domain:
                    self.gap_ch.properties[
                        t, self.gap_ch.length_domain.first()
                    ].flow_mass_phase_comp["Liq", "H2O"].fix(0)

                    if (
                        "Vap" in self.gap_ch.config.property_package.phase_list
                        and "H2O" in self.gap_ch.config.property_package.component_list
                    ):

                        self.gap_ch.properties[t, x].flow_mass_phase_comp[
                            "Vap", "H2O"
                        ].fix(0)

                        self.gap_ch.properties_interface[t, x].flow_mass_phase_comp[
                            "Vap", "H2O"
                        ].fix(0)

        add_object_reference(self, "length_domain", self.hot_ch.length_domain)
        add_object_reference(
            self, "difference_elements", self.hot_ch.difference_elements
        )
        add_object_reference(self, "first_element", self.hot_ch.first_element)
        add_object_reference(self, "nfe", self.hot_ch.nfe)

        # Add Ports
        self.add_inlet_port(name="hot_ch_inlet", block=self.hot_ch)
        self.add_outlet_port(name="hot_ch_outlet", block=self.hot_ch)
        self.add_inlet_port(name="cold_ch_inlet", block=self.cold_ch)
        self.add_outlet_port(name="cold_ch_outlet", block=self.cold_ch)
        if self.config.MD_configuration_Type == MDconfigurationType.GMD:
            self.add_inlet_port(name="gap_ch_inlet", block=self.gap_ch)
            self.add_outlet_port(name="gap_ch_outlet", block=self.gap_ch)

        self._add_heat_flux()
        self._add_mass_flux()

        if self.config.cold_ch.has_pressure_change:
            self.cold_ch._add_deltaP(
                pressure_change_type=self.config.cold_ch.pressure_change_type
            )
        if self.config.hot_ch.has_pressure_change:
            self.hot_ch._add_deltaP(
                pressure_change_type=self.config.hot_ch.pressure_change_type
            )

        self.recovery_mass = Var(
            self.flowsheet().config.time,
            initialize=0.1,
            bounds=(1e-11, 0.99),
            units=pyunits.dimensionless,
            doc="Mass-based recovery",
        )

        @self.Constraint(
            self.flowsheet().config.time,
            # self.config.hot_ch.property_package.component_list,
        )
        def eq_recovery_mass(b, t):
            return (
                b.recovery_mass[t]
                * b.hot_ch.properties[t, b.first_element].flow_mass_phase_comp[
                    "Liq", "H2O"
                ]
                == b.hot_ch.properties[t, b.first_element].flow_mass_phase_comp[
                    "Liq", "H2O"
                ]
                - b.hot_ch.properties[
                    t, b.hot_ch.length_domain.last()
                ].flow_mass_phase_comp["Liq", "H2O"]
            )

        self._add_mass_transfer()
        self._add_heat_transfer()

        @self.Expression(
            self.flowsheet().config.time,
            doc="Average thermal efficiency",
        )
        def thermal_efficiency(b, t):
            if self.config.MD_configuration_Type in [
                MDconfigurationType.DCMD,
                MDconfigurationType.GMD,
                MDconfigurationType.AGMD,
            ]:

                total_enth_flux = sum(
                    b.flux_enth_hot[t, x] for x in self.difference_elements
                )
                total_cond_heat_flux = sum(
                    b.flux_conduction_heat[t, x] for x in self.difference_elements
                )

                return total_enth_flux / (total_enth_flux + total_cond_heat_flux)
            else:
                return Constraint.Skip

        if self.config.MD_configuration_Type in [
            MDconfigurationType.DCMD,
            MDconfigurationType.GMD,
            MDconfigurationType.AGMD,
        ]:
            # to do: define effectiveness at the flowsheet level, particularly for VMD and DCMD configs
            @self.Expression(
                self.flowsheet().config.time,
                doc="module heat recovery: ratio of the actual heat recovered by cold side to the maximum ideal heat recovery",
            )
            def effectiveness(b, t):

                return (
                    b.cold_ch.properties[t, b.first_element].temperature
                    - b.cold_ch.properties[
                        t, b.cold_ch.length_domain.last()
                    ].temperature
                ) / (
                    b.hot_ch.properties[t, b.first_element].temperature
                    - b.cold_ch.properties[
                        t, b.cold_ch.length_domain.last()
                    ].temperature
                )

    def _add_mass_flux(self):

        solute_set = self.config.hot_ch.property_package.solute_set

        units_meta = (
            self.config.hot_ch.property_package.get_metadata().get_derived_units
        )

        # todo: add calculation method for permeability coefficient
        self.permeability_coef = Var(
            self.flowsheet().config.time,
            initialize=1e-10,
            bounds=(1e-11, 1e-9),
            units=units_meta("mass")
            * units_meta("time") ** -1
            * units_meta("length") ** -1
            * units_meta("pressure") ** -1,
            doc="membrane permeabiity coeffcient",
        )

        self.flux_mass = Var(
            self.flowsheet().config.time,
            self.difference_elements,
            initialize=1e-10,
            bounds=(1e-12, 1e-1),
            units=units_meta("mass")
            * units_meta("time") ** -1
            * units_meta("length") ** -2,
            doc="solvent vapor mass flux across the membrane",
        )

        self.flux_enth_hot = Var(
            self.flowsheet().config.time,
            self.difference_elements,
            initialize=1e3,
            bounds=(1e-10, 1e10),
            units=pyunits.J * pyunits.s**-1 * pyunits.m**-2,
            doc="hot side evaporation enthalpy flux",
        )

        self.flux_enth_cold = Var(
            self.flowsheet().config.time,
            self.difference_elements,
            initialize=1e3,
            bounds=(1e-10, 1e10),
            units=pyunits.J * pyunits.s**-1 * pyunits.m**-2,
            doc="cold side condensation enthalpy flux",
        )

        if self.config.MD_configuration_Type == MDconfigurationType.VMD:

            self.flux_expansion_heat = Var(
                self.flowsheet().config.time,
                self.difference_elements,
                initialize=1e3,
                bounds=(1e-10, 1e10),
                units=pyunits.J * pyunits.s**-1 * pyunits.m**-2,
                doc="heat of vapor expansion in VMD configuration",
            )

        self.dens_solvent = Param(
            initialize=1000,
            units=units_meta("mass") * units_meta("length") ** -3,
            doc="Pure water density",
        )

        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            doc="Solvent mass flux",
        )
        def eq_flux_mass(b, t, x):

            if self.config.MD_configuration_Type == MDconfigurationType.DCMD:
                return b.flux_mass[t, x] == b.permeability_coef[
                    t
                ] / b.membrane_thickness * (
                    b.hot_ch.properties_interface[t, x].pressure_sat
                    - b.cold_ch.properties_interface[t, x].pressure_sat
                )
            elif self.config.MD_configuration_Type in [
                MDconfigurationType.GMD,
                MDconfigurationType.AGMD,
            ]:
                return b.flux_mass[t, x] == b.permeability_coef[
                    t
                ] / b.membrane_thickness * (
                    b.hot_ch.properties_interface[t, x].pressure_sat
                    - b.gap_ch.properties_interface[t, x].pressure_sat
                )
            elif self.config.MD_configuration_Type == MDconfigurationType.VMD:
                return b.flux_mass[t, x] == b.permeability_coef[
                    t
                ] / b.membrane_thickness * (
                    b.hot_ch.properties_interface[t, x].pressure_sat
                    - b.cold_ch.properties[t, x].pressure
                )

        @self.Expression(
            self.flowsheet().config.time,
            doc="Average mass flux expression",
        )
        def flux_mass_avg(b, t):
            return sum(b.flux_mass[t, x] for x in self.difference_elements) / self.nfe

        if (
            self.config.hot_ch.concentration_polarization_type
            == ConcentrationPolarizationType.calculated
        ):

            @self.Constraint(
                self.flowsheet().config.time,
                self.difference_elements,
                solute_set,
                doc="External Concentration polarization in hot side",
            )
            def eq_concentration_polarization_hot_ch(b, t, x, j):
                jw = b.flux_mass[t, x] / self.dens_solvent
                exponent = jw / self.hot_ch.K[t, x, j]
                return b.hot_ch.properties_interface[t, x].conc_mass_phase_comp[
                    "Liq", j
                ] == b.hot_ch.properties[t, x].conc_mass_phase_comp["Liq", j] * exp(
                    exponent
                )

        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            doc="hot side evaporation enthalpy flux",
        )
        def eq_flux_hot_enth(b, t, x):
            return (
                b.flux_enth_hot[t, x]
                == b.flux_mass[t, x]
                * b.hot_ch.properties_vapor[t, x].enth_mass_phase["Vap"]
            )

        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            doc="cold side evaporation enthalpy flux",
        )
        def eq_flux_cold_enth(b, t, x):
            if self.config.MD_configuration_Type == MDconfigurationType.DCMD:
                return (
                    b.flux_enth_cold[t, x]
                    == b.flux_mass[t, x]
                    * b.cold_ch.properties_vapor[t, x].enth_mass_phase["Vap"]
                )
            elif self.config.MD_configuration_Type == MDconfigurationType.GMD:
                return (
                    b.flux_enth_cold[t, x]
                    == b.flux_mass[t, x]
                    * b.gap_ch.properties[t, x].enth_mass_phase["Liq"]
                )
            elif self.config.MD_configuration_Type == MDconfigurationType.VMD:
                return (
                    b.flux_enth_cold[t, x]
                    == b.flux_mass[t, x]
                    * b.cold_ch.properties[t, x].enth_mass_phase["Vap"]
                )

        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            doc="vapor expansion heat flux",
        )
        def eq_flux_expansion_heat(b, t, x):
            if self.config.MD_configuration_Type == MDconfigurationType.VMD:

                R = 8.3146261 * pyunits.J / pyunits.mol / pyunits.K
                molar_density = 18.01528e-3 * pyunits.kg / pyunits.mol
                T = b.hot_ch.properties_vapor[t, x].temperature
                P_f = b.hot_ch.properties_vapor[t, x].pressure
                P_p = b.cold_ch.properties[t, x].pressure

                return b.flux_expansion_heat[t, x] == b.flux_mass[
                    t, x
                ] / molar_density * R * T * log(P_f / P_p)
            else:
                return Constraint.Skip

        @self.Expression(
            self.flowsheet().config.time,
            doc="Average hot side enthalpy flux expression",
        )
        def flux_enth_hot_avg(b, t):
            return (
                sum(b.flux_enth_hot[t, x] for x in self.difference_elements) / self.nfe
            )

        @self.Expression(
            self.flowsheet().config.time,
            doc="Average cold side enthalpy flux expression",
        )
        def flux_enth_cold_avg(b, t):
            return (
                sum(b.flux_enth_cold[t, x] for x in self.difference_elements) / self.nfe
            )

        @self.Expression(
            self.flowsheet().config.time,
            doc="Average expansion heat",
        )
        def flux_expansion_heat_avg(b, t):
            if self.config.MD_configuration_Type == MDconfigurationType.VMD:
                return (
                    sum(b.flux_expansion_heat[t, x] for x in self.difference_elements)
                    / self.nfe
                )
            else:
                return Constraint.Skip

        @self.Constraint(
            self.flowsheet().time,
            self.difference_elements,
            doc="hot side Vapor temperature",
        )
        def eq_vapor_temperature_hot(b, t, x):
            return (
                b.hot_ch.properties_vapor[t, x].temperature
                == b.hot_ch.properties_interface[t, x].temperature
            )

        @self.Constraint(
            self.flowsheet().time,
            self.difference_elements,
            doc="hot  side Vapor pressure",
        )
        def eq_vapor_pressure_hot(b, t, x):
            return (
                b.hot_ch.properties_vapor[t, x].pressure
                == b.hot_ch.properties_interface[t, x].pressure_sat
            )

        @self.Constraint(
            self.flowsheet().time,
            self.difference_elements,
            doc="cold side Vapor temperature",
        )
        def eq_vapor_temperature_cold(b, t, x):
            if self.config.MD_configuration_Type == MDconfigurationType.DCMD:
                return (
                    b.cold_ch.properties_vapor[t, x].temperature
                    == b.cold_ch.properties_interface[t, x].temperature
                )
            elif self.config.MD_configuration_Type == MDconfigurationType.VMD:
                return (
                    b.cold_ch.properties[t, x].temperature
                    == b.hot_ch.properties_interface[t, x].temperature
                )
            else:
                return Constraint.Skip

        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            doc="gap bulk temperature in GMD",
        )
        def gap_bulk_temperature(b, t, x):
            # assuming linear temperature change across the gap

            if self.config.MD_configuration_Type == MDconfigurationType.GMD:
                if (
                    self.config.gap_ch.temperature_polarization_type
                    == TemperaturePolarizationType.fixed
                ):
                    return (
                        b.gap_ch.properties_interface[t, x].temperature
                        + b.cold_ch.properties_interface[t, x].temperature
                        == 2 * b.gap_ch.properties[t, x].temperature
                    )

            else:
                return Constraint.Skip

        @self.Constraint(
            self.flowsheet().time,
            self.difference_elements,
            doc="cold  side Vapor pressure",
        )
        def eq_vapor_pressure_cold(b, t, x):
            if self.config.MD_configuration_Type == MDconfigurationType.DCMD:
                return (
                    b.cold_ch.properties_vapor[t, x].pressure
                    == b.cold_ch.properties_interface[t, x].pressure_sat
                )
            else:
                return Constraint.Skip

        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            doc="Vapor flow rate and no liquid in vapor state block",
        )
        def eq_vapor_flow(b, t, x):
            lb = b.hot_ch.properties_vapor[t, x].flow_mass_phase_comp["Liq", "H2O"].lb
            b.hot_ch.properties_vapor[t, x].flow_mass_phase_comp["Liq", "H2O"].fix(lb)

            return (
                b.hot_ch.properties_vapor[t, x].flow_mass_phase_comp["Vap", "H2O"]
                == b.flux_mass[t, x] * b.area / self.nfe
            )

        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            doc="Equal vapor flow rate between hot_ch and cold_ch",
        )
        def eq_vapor_flow_equal(b, t, x):
            lb = b.hot_ch.properties_vapor[t, x].flow_mass_phase_comp["Liq", "H2O"].lb
            b.hot_ch.properties_vapor[t, x].flow_mass_phase_comp["Liq", "H2O"].fix(lb)
            if self.config.MD_configuration_Type == MDconfigurationType.DCMD:
                return (
                    b.hot_ch.properties_vapor[t, x].flow_mass_phase_comp["Vap", "H2O"]
                    == b.cold_ch.properties_vapor[t, x].flow_mass_phase_comp[
                        "Vap", "H2O"
                    ]
                )
            else:
                return Constraint.Skip

        if (
            self.config.hot_ch.temperature_polarization_type
            != TemperaturePolarizationType.none
        ):

            @self.Constraint(
                self.flowsheet().config.time,
                self.difference_elements,
                doc="Temperature polarization in hot channel",
            )
            def eq_temperature_polarization_hot(b, t, x):
                if self.config.MD_configuration_Type == MDconfigurationType.VMD:
                    return (
                        b.hot_ch.h_conv[t, x]
                        * (
                            b.hot_ch.properties[t, x].temperature
                            - b.hot_ch.properties_interface[t, x].temperature
                        )
                        == b.flux_expansion_heat[t, x]
                        + b.flux_enth_hot[t, x]
                        - b.flux_mass[t, x]
                        * b.hot_ch.properties[t, x].enth_mass_phase["Liq"]
                    )

                else:

                    return (
                        b.hot_ch.h_conv[t, x]
                        * (
                            b.hot_ch.properties[t, x].temperature
                            - b.hot_ch.properties_interface[t, x].temperature
                        )
                        == b.flux_conduction_heat[t, x]
                        + b.flux_enth_hot[t, x]
                        - b.flux_mass[t, x]
                        * b.hot_ch.properties[t, x].enth_mass_phase["Liq"]
                    )

        if (
            self.config.cold_ch.temperature_polarization_type
            != TemperaturePolarizationType.none
        ):

            @self.Constraint(
                self.flowsheet().config.time,
                self.difference_elements,
                doc="Temperature polarization in cold channel",
            )
            def eq_temperature_polarization_cold(b, t, x):
                if self.config.MD_configuration_Type == MDconfigurationType.GMD:
                    return (
                        b.cold_ch.h_conv[t, x]
                        * (
                            -b.cold_ch.properties[t, x].temperature
                            + b.cold_ch.properties_interface[t, x].temperature
                        )
                        == b.flux_conduction_heat_gap[t, x]
                    )
                elif self.config.MD_configuration_Type == MDconfigurationType.DCMD:
                    return (
                        b.cold_ch.h_conv[t, x]
                        * (
                            -b.cold_ch.properties[t, x].temperature
                            + b.cold_ch.properties_interface[t, x].temperature
                        )
                        == b.flux_conduction_heat[t, x]
                        + b.flux_enth_cold[t, x]
                        - b.flux_mass[t, x]
                        * b.cold_ch.properties[t, x].enth_mass_phase["Liq"]
                    )
                else:
                    return Constraint.Skip

            return self.eq_flux_mass

    def _add_heat_flux(self):

        # todo: add calculation method for permeability coefficient
        self.membrane_thickness = Var(
            initialize=1e-4,
            bounds=(1e-5, 1e-2),
            doc="membrane thickness",
            units=pyunits.m,
        )

        if self.config.MD_configuration_Type != MDconfigurationType.VMD:

            self.flux_conduction_heat = Var(
                self.flowsheet().config.time,
                self.difference_elements,
                initialize=1e3,
                bounds=(0, 1e10),
                units=pyunits.J * pyunits.s**-1 * pyunits.m**-2,
                doc="conduction heat flux",
            )
            self.membrane_thermal_conductivity = Var(
                initialize=0.2,
                bounds=(0, 1),
                units=pyunits.J * pyunits.s**-1 * pyunits.K**-1 * pyunits.m**-1,
                doc="Thermal conductivity coefficient of the membrane",
            )

        if self.config.MD_configuration_Type == MDconfigurationType.GMD:
            self.gap_thermal_conductivity = Var(
                initialize=0.06,
                bounds=(0, 1),
                units=pyunits.J * pyunits.s**-1 * pyunits.K**-1 * pyunits.m**-1,
                doc="Thermal conductivity coefficient of the gap",
            )

            self.gap_thickness = Var(
                initialize=1e-4,
                bounds=(1e-5, 1e-2),
                doc="gap thickness",
                units=pyunits.m,
            )

            self.flux_conduction_heat_gap = Var(
                self.flowsheet().config.time,
                self.difference_elements,
                initialize=10e3,
                bounds=(0, 1e20),
                units=pyunits.J * pyunits.s**-1 * pyunits.m**-2,
                doc="conduction heat across the gap",
            )

        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            doc="conduction heat flux",
        )
        def eq_flux_heat(b, t, x):
            if self.config.MD_configuration_Type == MDconfigurationType.DCMD:
                return b.flux_conduction_heat[
                    t, x
                ] == b.membrane_thermal_conductivity / b.membrane_thickness * (
                    b.hot_ch.properties_interface[t, x].temperature
                    - b.cold_ch.properties_interface[t, x].temperature
                )
            elif self.config.MD_configuration_Type == MDconfigurationType.GMD:
                return b.flux_conduction_heat[
                    t, x
                ] == b.membrane_thermal_conductivity / b.membrane_thickness * (
                    b.hot_ch.properties_interface[t, x].temperature
                    - b.gap_ch.properties_interface[t, x].temperature
                )
            elif self.config.MD_configuration_Type == MDconfigurationType.VMD:
                return Constraint.Skip

        @self.Expression(
            self.flowsheet().config.time,
            doc="Average conduction heat flux expression",
        )
        def flux_conduction_heat_avg(b, t):
            if self.config.MD_configuration_Type != MDconfigurationType.VMD:
                return (
                    sum(b.flux_conduction_heat[t, x] for x in self.difference_elements)
                    / self.nfe
                )
            else:
                return Constraint.Skip

        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            doc="heat conduction across the gap in GMD configuration",
        )
        def eq_flux_conduction_gap(b, t, x):
            if self.config.MD_configuration_Type == MDconfigurationType.GMD:

                return b.flux_conduction_heat_gap[
                    t, x
                ] == b.gap_thermal_conductivity / b.gap_thickness * (
                    b.gap_ch.properties_interface[t, x].temperature
                    - b.cold_ch.properties_interface[t, x].temperature
                )
            else:
                return Constraint.Skip

        @self.Expression(
            self.flowsheet().config.time,
            doc="Average conduction heat flux across the gap",
        )
        def flux_conduction_heat_gap_avg(b, t):
            if self.config.MD_configuration_Type == MDconfigurationType.GMD:
                return (
                    sum(
                        b.flux_conduction_heat_gap[t, x]
                        for x in self.difference_elements
                    )
                    / self.nfe
                )
            else:
                return Constraint.Skip

        return self.eq_flux_heat

    def _add_length_and_width(self):
        units_meta = (
            self.config.hot_ch.property_package.get_metadata().get_derived_units
        )
        self.length = Var(
            initialize=10,
            bounds=(0.1, 5e2),
            domain=NonNegativeReals,
            units=units_meta("length"),
            doc="Effective membrane length",
        )
        self.width = Var(
            initialize=1,
            bounds=(1e-1, 1e3),
            domain=NonNegativeReals,
            units=units_meta("length"),
            doc="Membrane width",
        )

    def _add_area(self, include_constraint=True):

        if not hasattr(self, "area"):
            self.area = Var(
                initialize=10,
                bounds=(1e-10, 1e10),
                domain=NonNegativeReals,
                units=pyunits.m**2,
                doc="Total Membrane area",
            )

        if include_constraint:
            if not hasattr(self, "eq_area"):

                @self.Constraint(doc="Total Membrane area")
                def eq_area(b):
                    return b.area == b.length * b.width

            else:
                raise ValueError(
                    "include_constraint was set to True inside of _add_area(), but area constraint already exists."
                )

    def initialize_build(
        self,
        initialize_guess=None,
        state_args_hot_ch=None,
        state_args_cold_ch=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):

        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        hot_ch_flags = self.hot_ch.initialize(
            state_args=state_args_hot_ch,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            initialize_guess=initialize_guess,
            type="hot_ch",
        )

        init_log.info("Initialization Step 1a (hot channel) Complete")

        cold_ch_flags = self.cold_ch.initialize(
            state_args=state_args_cold_ch,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            initialize_guess=initialize_guess,
            type="cold_ch",
        )

        init_log.info("Initialization Step 1b (cold channel) Complete")

        if self.config.MD_configuration_Type == MDconfigurationType.GMD:
            gap_ch_flags = self.gap_ch.initialize(
                state_args=state_args_cold_ch,
                outlvl=outlvl,
                optarg=optarg,
                solver=solver,
                initialize_guess=initialize_guess,
                type="cold_ch",
            )

            init_log.info("Initialization Step 1c (gap channel) Complete")

        if self.config.MD_configuration_Type == MDconfigurationType.VMD:
            for x in [self.cold_ch.length_domain.first()]:
                self.cold_ch_inlet.temperature[0].unfix()

        if self.config.MD_configuration_Type == MDconfigurationType.GMD:
            for x in [self.gap_ch.length_domain.first()]:
                self.gap_ch_inlet.temperature[0].unfix()

        if degrees_of_freedom(self) != 0:
            raise Exception(
                f"{self.name} degrees of freedom were not 0 at the beginning "
                f"of initialization. DoF = {degrees_of_freedom(self)}"
            )

        # pre-solve using interval arithmetic
        interval_initializer(self)

        # solver
        opt = get_solver(solver, optarg)

        # Solve unit *without* any flux equations
        self.eq_flux_mass.deactivate()
        self.eq_flux_heat.deactivate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info(f"Initialization Step 2 (No Flux) {idaeslog.condition(res)}")

        # Activate only the heat flux equations
        self.eq_flux_heat.activate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info(
            f"Initialization Step 2 (heat Flux only) {idaeslog.condition(res)}"
        )

        # Activate mass flux equations as well
        self.eq_flux_mass.activate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info(
            f"Initialization Step 4 (Heat and Mass Flux) {idaeslog.condition(res)}"
        )

        # Release inlet state
        self.cold_ch.release_state(cold_ch_flags, outlvl)
        self.hot_ch.release_state(hot_ch_flags, outlvl)
        if self.config.MD_configuration_Type == MDconfigurationType.GMD:
            self.gap_ch.release_state(gap_ch_flags, outlvl)

        init_log.info(f"Initialization Complete: {idaeslog.condition(res)}")

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {self.name} failed to initialize")

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Hot channel Inlet": self.hot_ch_inlet,
                "Hot channel Outlet": self.hot_ch_outlet,
                "Cold channel Inlet": self.cold_ch_inlet,
                "Cold channel Outlet": self.cold_ch_outlet,
            },
            time_point=time_point,
        )

    def _get_performance_contents(self, time_point=0):
        var_dict = {}
        expr_dict = {}
        var_dict["Mass Recovery Rate"] = self.recovery_mass[time_point]
        var_dict["Membrane Area"] = self.area
        if hasattr(self, "length"):
            var_dict["Membrane Length"] = self.length
        if hasattr(self, "width"):
            var_dict["Membrane Width"] = self.width

        expr_dict["Average Solute Flux"] = self.flux_mass_avg[time_point]

        return {"vars": var_dict, "exprs": expr_dict}

    def calculate_scaling_factors(self):

        if iscale.get_scaling_factor(self.dens_solvent) is None:
            sf = iscale.get_scaling_factor(
                self.hot_ch.properties[0, self.first_element].dens_mass_phase["Liq"]
            )
            iscale.set_scaling_factor(self.dens_solvent, sf)

        super().calculate_scaling_factors()
        for t, v in self.recovery_mass.items():
            sf = 1
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, sf)

        if iscale.get_scaling_factor(self.area) is None:
            sf = iscale.get_scaling_factor(self.area, default=10, warning=True)
            iscale.set_scaling_factor(self.area, sf)

        if iscale.get_scaling_factor(self.permeability_coef) is None:
            iscale.set_scaling_factor(self.permeability_coef, 1e10)

        if iscale.get_scaling_factor(self.membrane_thickness) is None:
            iscale.set_scaling_factor(self.membrane_thickness, 1e4)

        if hasattr(self, "membrane_thermal_conductivity"):
            if iscale.get_scaling_factor(self.membrane_thermal_conductivity) is None:
                iscale.set_scaling_factor(self.membrane_thermal_conductivity, 10)

        if hasattr(self, "gap_thermal_conductivity"):
            if iscale.get_scaling_factor(self.gap_thermal_conductivity) is None:
                iscale.set_scaling_factor(self.gap_thermal_conductivity, 10)

        if hasattr(self, "gap_thermal_conductivity"):
            if iscale.get_scaling_factor(self.gap_thickness) is None:
                iscale.set_scaling_factor(self.gap_thickness, 1e4)

        for (t, x), v in self.flux_mass.items():
            if iscale.get_scaling_factor(v) is None:

                sf_pressure_hot = iscale.get_scaling_factor(
                    self.hot_ch.properties_interface[t, x].pressure_sat
                )

                sf_permeability = iscale.get_scaling_factor(self.permeability_coef)
                sf_thickness = iscale.get_scaling_factor(self.membrane_thickness)

                sf_flux = sf_permeability * (sf_pressure_hot) / sf_thickness

                iscale.set_scaling_factor(v, sf_flux)

                iscale.constraint_scaling_transform(self.eq_flux_mass[t, x], sf_flux)

        for (t, x), v in self.flux_enth_hot.items():
            if iscale.get_scaling_factor(v) is None:
                sf_flux_enth = sf_flux * iscale.get_scaling_factor(
                    self.hot_ch.properties_vapor[t, x].enth_mass_phase["Vap"]
                )
                iscale.set_scaling_factor(v, sf_flux_enth)
            sf = iscale.get_scaling_factor(
                self.hot_ch.properties_vapor[t, x].flow_mass_phase_comp["Vap", "H2O"]
            )
            iscale.set_scaling_factor(
                self.hot_ch.properties_vapor[t, x].flow_mass_phase_comp["Vap", "H2O"],
                sf * 1000,
            )

        for (t, x), v in self.flux_enth_cold.items():
            if iscale.get_scaling_factor(v) is None:
                if self.config.MD_configuration_Type == MDconfigurationType.VMD:
                    sf_flux_enth = sf_flux * iscale.get_scaling_factor(
                        self.cold_ch.properties[t, x].enth_mass_phase["Vap"]
                    )
                    iscale.set_scaling_factor(v, sf_flux_enth)
                    sf = iscale.get_scaling_factor(
                        self.cold_ch.properties[t, x].flow_mass_phase_comp["Vap", "H2O"]
                    )
                    iscale.set_scaling_factor(
                        self.cold_ch.properties[t, x].flow_mass_phase_comp[
                            "Vap", "H2O"
                        ],
                        sf * 50,
                    )

                elif self.config.MD_configuration_Type == MDconfigurationType.DCMD:
                    sf_flux_enth = sf_flux * iscale.get_scaling_factor(
                        self.cold_ch.properties_vapor[t, x].enth_mass_phase["Vap"]
                    )
                    iscale.set_scaling_factor(v, sf_flux_enth)
                    sf = iscale.get_scaling_factor(
                        self.cold_ch.properties_vapor[t, x].flow_mass_phase_comp[
                            "Vap", "H2O"
                        ]
                    )
                    iscale.set_scaling_factor(
                        self.cold_ch.properties_vapor[t, x].flow_mass_phase_comp[
                            "Vap", "H2O"
                        ],
                        sf * 1000,
                    )

                    iscale.set_scaling_factor(
                        self.cold_ch.properties[t, x].pressure,
                        1e-4,
                    )

        if hasattr(self, "flux_conduction_heat"):
            for (t, x), v in self.flux_conduction_heat.items():
                if iscale.get_scaling_factor(v) is None:
                    sf_flux_cond = (
                        iscale.get_scaling_factor(self.membrane_thermal_conductivity)
                        / iscale.get_scaling_factor(self.membrane_thickness)
                        * iscale.get_scaling_factor(
                            self.hot_ch.properties_interface[t, x].temperature
                        )
                    )
                    iscale.set_scaling_factor(v, sf_flux_cond)

        if hasattr(self, "flux_conduction_heat_gap"):
            for (t, x), v in self.flux_conduction_heat_gap.items():
                if iscale.get_scaling_factor(v) is None:
                    sf_flux_cond = (
                        iscale.get_scaling_factor(self.membrane_thermal_conductivity)
                        / iscale.get_scaling_factor(self.membrane_thickness)
                        * iscale.get_scaling_factor(
                            self.hot_ch.properties_interface[t, x].temperature
                        )
                    )
                    iscale.set_scaling_factor(v, sf_flux_cond)

        if hasattr(self, "flux_expansion_heat"):
            for (t, x), v in self.flux_expansion_heat.items():
                if iscale.get_scaling_factor(v) is None:
                    sf_flux_expansion = 1e-3
                    iscale.set_scaling_factor(v, sf_flux_expansion)

        if hasattr(self, "length"):
            if iscale.get_scaling_factor(self.length) is None:
                iscale.set_scaling_factor(self.length, 1)

        if hasattr(self, "width"):
            if iscale.get_scaling_factor(self.width) is None:
                iscale.set_scaling_factor(self.width, 1)

    @property
    def default_costing_method(self):
        return cost_membrane_distillation
