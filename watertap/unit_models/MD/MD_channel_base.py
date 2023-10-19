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

""""
Base constraints and methods for membrane distillation channel
__author__ = "Elmira Shamlou"

"""
from copy import deepcopy
from enum import Enum, auto
from pyomo.environ import (
    Constraint,
    Param,
    NonNegativeReals,
    Var,
    units as pyunits,
)
from idaes.core import (
    FlowDirection,
)
from idaes.core.util import scaling as iscale
from idaes.core.util.exceptions import ConfigurationError


class TemperaturePolarizationType(Enum):
    """
    none: no temperature polarization
    fixed: convective heat transfer coefficient is a user specified value
    calculated: calculate convective heat transfer coefficient
    """

    none = auto()
    fixed = auto()
    calculated = auto()


class ConcentrationPolarizationType(Enum):
    """
    none: no concentration polarization
    fixed: concentration polarization modulus is a user specified value
    calculated: calculate concentration polarization (concentration at membrane interface)
    """

    none = auto()
    fixed = auto()
    calculated = auto()


class MassTransferCoefficient(Enum):
    """
    none: mass transfer coefficient not utilized for concentration polarization effect
    fixed: mass transfer coefficient is a user specified value
    calculated: mass transfer coefficient is calculated
    """

    none = auto()
    fixed = auto()
    calculated = auto()
    # TODO: add option for users to define their own relationship?


class PressureChangeType(Enum):
    """
    fixed_per_stage: pressure drop across membrane channel is a user-specified value
    fixed_per_unit_length: pressure drop per unit length across membrane channel is a user-specified value
    calculated: pressure drop across membrane channel is calculated
    """

    fixed_per_stage = auto()
    fixed_per_unit_length = auto()
    calculated = auto()


class FrictionFactor(Enum):
    """
    flat_sheet: Darcy's friction factor correlation by Guillen & Hoek
    spiral_wound: Darcy's friction factor correlation by Schock & Miquel, 1987
    """

    flat_sheet = auto()
    spiral_wound = auto()


class MDChannelMixin:
    def _skip_element(self, x):
        raise NotImplementedError()

    def _set_nfe(self):
        self.first_element = self.length_domain.first()
        self.last_element = self.length_domain.last()

        self.nfe = Param(
            initialize=(len(self.difference_elements)),
            units=pyunits.dimensionless,
            doc="Number of finite elements",
        )

    def add_total_pressure_balances(
        self,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        custom_term=None,
        friction_factor=FrictionFactor.flat_sheet,
    ):
        super().add_total_pressure_balances(
            has_pressure_change=has_pressure_change, custom_term=custom_term
        )

        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            doc="Pressure at interface",
        )
        def eq_equal_pressure_interface(b, t, x):
            if b._skip_element(x):
                return Constraint.Skip
            return b.properties_interface[t, x].pressure == b.properties[t, x].pressure

        if has_pressure_change:
            self._add_pressure_change(pressure_change_type=pressure_change_type)

        if pressure_change_type == PressureChangeType.calculated:
            self._add_calculated_pressure_change(friction_factor=friction_factor)

    def add_temperature_polarization(
        self,
        temperature_polarization_type=TemperaturePolarizationType.none,
    ):

        units_meta = self.config.property_package.get_metadata().get_derived_units

        if temperature_polarization_type == TemperaturePolarizationType.none:

            @self.Constraint(
                self.flowsheet().config.time,
                self.length_domain,
                doc="No temperature polarization",
            )
            def eq_no_temp_pol(b, t, x):
                if b._skip_element(x):
                    return Constraint.Skip
                return (
                    b.properties_interface[t, x].temperature
                    == b.properties[t, x].temperature
                )

            return self.eq_no_temp_pol

        elif temperature_polarization_type not in (
            TemperaturePolarizationType.fixed,
            TemperaturePolarizationType.calculated,
        ):
            raise ConfigurationError(
                f"Unrecognized temperature_polarization_type {temperature_polarization_type}"
            )

        self.h_conv = Var(
            self.flowsheet().config.time,
            self.length_domain,
            initialize=2000,
            bounds=(1, 5e4),
            domain=NonNegativeReals,
            units=pyunits.J * pyunits.s**-1 * pyunits.K**-1 * pyunits.m**-2,
            doc="Convective heat transfer coefficient",
        )

        if temperature_polarization_type == TemperaturePolarizationType.calculated:
            self._add_calculated_convective_heat_transfer_coefficient()

        return self.h_conv

    def _add_calculated_convective_heat_transfer_coefficient(self):
        self._add_calculated_pressure_change_mass_heat_transfer_components()

        self.N_Pr = Var(
            self.flowsheet().config.time,
            self.length_domain,
            initialize=1,
            bounds=(0.1, 50),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Prandtl number in membrane channel",
        )
        self.N_Nu = Var(
            self.flowsheet().config.time,
            self.length_domain,
            initialize=1e2,
            bounds=(1, 3e2),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Nusselt number in membrane channel",
        )

        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            doc="Convective heat transfer coefficient",
        )
        def eq_h_conv(b, t, x):
            if b._skip_element(x):
                return Constraint.Skip
            return (
                b.h_conv[t, x]
                == b.properties[t, x].therm_cond_phase["Liq"] * b.N_Nu[t, x] / b.dh
            )

        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            doc="Nusselt number",
        )
        def eq_N_Nu(b, t, x):
            if b._skip_element(x):
                return Constraint.Skip
            return b.N_Nu[t, x] == 0.162 * b.N_Re[t, x] ** 0.656 * b.N_Pr[t, x] ** 0.333

        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            doc="Prandtl number",
        )
        def eq_N_Pr(b, t, x):
            if b._skip_element(x):
                return Constraint.Skip
            return (
                b.N_Pr[t, x] * b.properties[t, x].therm_cond_phase["Liq"]
                == b.properties[t, x].visc_d_phase["Liq"]
                * b.properties[t, x].cp_mass_phase["Liq"]
            )

        return self.eq_h_conv

    def add_concentration_polarization(
        self,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
    ):

        solute_set = self.config.property_package.solute_set
        units_meta = self.config.property_package.get_metadata().get_derived_units

        if concentration_polarization_type == ConcentrationPolarizationType.none:

            @self.Constraint(
                self.flowsheet().config.time,
                self.length_domain,
                solute_set,
                doc="Unit concentration polarization modulus",
            )
            def eq_cp_modulus(b, t, x, j):
                if b._skip_element(x):
                    return Constraint.Skip
                return (
                    b.properties_interface[t, x].conc_mass_phase_comp["Liq", j]
                    == b.properties[t, x].conc_mass_phase_comp["Liq", j]
                )

            return self.eq_cp_modulus

        if concentration_polarization_type not in (
            ConcentrationPolarizationType.fixed,
            ConcentrationPolarizationType.calculated,
        ):
            raise ConfigurationError(
                f"Unrecognized concentration_polarization_type {concentration_polarization_type}"
            )

        self.cp_modulus = Var(
            self.flowsheet().config.time,
            self.length_domain,
            solute_set,
            initialize=1.1,
            bounds=(0.1, 3),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Concentration polarization modulus",
        )

        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            solute_set,
            doc="Concentration polarization modulus",
        )
        def eq_cp_modulus(b, t, x, j):  # pylint: disable=function-redefined
            if b._skip_element(x):
                return Constraint.Skip
            return (
                self.properties_interface[t, x].conc_mass_phase_comp["Liq", j]
                == self.properties[t, x].conc_mass_phase_comp["Liq", j]
                * self.cp_modulus[t, x, j]
            )

        if concentration_polarization_type == ConcentrationPolarizationType.calculated:
            if mass_transfer_coefficient == MassTransferCoefficient.none:
                raise ConfigurationError()

            # mass_transfer_coefficient is either calculated or fixed
            self.K = Var(
                self.flowsheet().config.time,
                self.length_domain,
                solute_set,
                initialize=5e-5,
                bounds=(1e-6, 1e-3),
                domain=NonNegativeReals,
                units=units_meta("length") * units_meta("time") ** -1,
                doc="Membrane channel mass transfer coefficient",
            )

            if mass_transfer_coefficient == MassTransferCoefficient.calculated:
                self._add_calculated_mass_transfer_coefficient()

        return self.eq_cp_modulus

    ## should be called by add concentration polarization
    def _add_calculated_mass_transfer_coefficient(self):
        self._add_calculated_pressure_change_mass_heat_transfer_components()

        solute_set = self.config.property_package.solute_set

        self.N_Sc_comp = Var(
            self.flowsheet().config.time,
            self.length_domain,
            solute_set,
            initialize=5e2,
            bounds=(1e2, 2e3),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Schmidt number in membrane channel",
        )
        self.N_Sh_comp = Var(
            self.flowsheet().config.time,
            self.length_domain,
            solute_set,
            initialize=1e2,
            bounds=(1, 3e2),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Sherwood number in membrane channel",
        )

        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            self.config.property_package.solute_set,
            doc="Mass transfer coefficient in membrane channel",
        )
        def eq_K(b, t, x, j):
            if b._skip_element(x):
                return Constraint.Skip
            return (
                b.K[t, x, j] * b.dh
                # TODO: add diff coefficient to SW prop and consider multi-components
                == b.properties[t, x].diffus_phase_comp["Liq", j] * b.N_Sh_comp[t, x, j]
            )

        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            self.config.property_package.solute_set,
            doc="Sherwood number",
        )
        def eq_N_Sh_comp(b, t, x, j):
            return (
                b.N_Sh_comp[t, x, j]
                == 0.2 * b.N_Re[t, x] ** 0.57 * b.N_Sc_comp[t, x, j] ** 0.4
            )

        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            self.config.property_package.solute_set,
            doc="Schmidt number",
        )
        def eq_N_Sc_comp(b, t, x, j):
            return (
                b.N_Sc_comp[t, x, j]
                * b.properties[t, x].dens_mass_phase["Liq"]
                * b.properties[t, x].diffus_phase_comp["Liq", j]
                == b.properties[t, x].visc_d_phase["Liq"]
            )

        return self.eq_K

    def _add_calculated_pressure_change_mass_heat_transfer_components(self):
        # NOTE: This function could be called by either
        # `_add_calculated_pressure_change` *and/or*
        # `_add_calculated_mass_transfer_coefficient`.
        # Therefore, we add this simple guard against it being called twice.
        if hasattr(self, "channel_height"):
            return

        if not hasattr(self, "width"):
            raise ConfigurationError(
                f"Due to either a calculated mass transfer coefficient or a calculated pressure change, a ``width`` variable needs to be supplied to `add_geometry` for this MembraneChannel"
            )

        units_meta = self.config.property_package.get_metadata().get_derived_units

        if not hasattr(self, "area"):
            # comes from ControlVolume1D for 1DMC
            self.area = Var(
                initialize=1e-3 * 1 * 0.95,
                bounds=(0, 1e3),
                domain=NonNegativeReals,
                units=units_meta("length") ** 2,
                doc="Cross sectional area",
            )

        self.channel_height = Var(
            initialize=1e-3,
            bounds=(1e-4, 5e-3),
            domain=NonNegativeReals,
            units=units_meta("length"),
            doc="membrane-channel height",
        )

        self.dh = Var(
            initialize=1e-3,
            bounds=(1e-4, 5e-3),
            domain=NonNegativeReals,
            units=units_meta("length"),
            doc="Hydraulic diameter of membrane channel",
        )

        self.spacer_porosity = Var(
            initialize=0.95,
            bounds=(0.1, 0.99),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="membrane-channel spacer porosity",
        )

        self.N_Re = Var(
            self.flowsheet().config.time,
            self.length_domain,
            initialize=5e2,
            bounds=(10, 5e3),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Reynolds number in membrane channel",
        )

        @self.Constraint(doc="Hydraulic diameter")  # eqn. 17 in Schock & Miquel, 1987
        def eq_dh(b):
            return b.dh == 4 * b.spacer_porosity / (
                2 / b.channel_height + (1 - b.spacer_porosity) * 8 / b.channel_height
            )

        @self.Constraint(doc="Cross-sectional area")
        def eq_area(b):
            return b.area == b.channel_height * b.width * b.spacer_porosity

        @self.Constraint(
            self.flowsheet().config.time, self.length_domain, doc="Reynolds number"
        )
        def eq_N_Re(b, t, x):
            return (
                b.N_Re[t, x] * b.area * b.properties[t, x].visc_d_phase["Liq"]
                == sum(
                    b.properties[t, x].flow_mass_phase_comp["Liq", j]
                    for j in b.config.property_package.component_list
                )
                * b.dh
            )

    def _add_calculated_pressure_change(
        self, friction_factor=FrictionFactor.flat_sheet
    ):
        self._add_calculated_pressure_change_mass_heat_transfer_components()

        units_meta = self.config.property_package.get_metadata().get_derived_units

        self.velocity = Var(
            self.flowsheet().config.time,
            self.length_domain,
            initialize=0.5,
            bounds=(1e-2, 5),
            domain=NonNegativeReals,
            units=units_meta("length") / units_meta("time"),
            doc="Crossflow velocity in feed channel",
        )
        self.friction_factor_darcy = Var(
            self.flowsheet().config.time,
            self.length_domain,
            initialize=0.5,
            bounds=(1e-2, 5),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Darcy friction factor in feed channel",
        )

        # Constraints active when MassTransferCoefficient.calculated
        # Mass transfer coefficient calculation

        ## ==========================================================================
        # Crossflow velocity
        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            doc="Crossflow velocity constraint",
        )
        def eq_velocity(b, t, x):
            return b.velocity[t, x] * b.area == b.properties[t, x].flow_vol_phase["Liq"]

        ## ==========================================================================
        # Darcy friction factor based on eq. S27 in SI for Cost Optimization of Osmotically Assisted Reverse Osmosis
        if friction_factor == FrictionFactor.flat_sheet:

            @self.Constraint(
                self.flowsheet().config.time,
                self.length_domain,
                doc="Darcy friction factor constraint for flat sheet membranes",
            )
            def eq_friction_factor(b, t, x):
                return (b.friction_factor_darcy[t, x] - 0.42) * b.N_Re[t, x] == 189.3

        # Darcy friction factor based on eq. 24 in Mass transfer and pressure loss in spiral wound modules (Schock & Miquel, 1987)
        elif friction_factor == FrictionFactor.spiral_wound:

            @self.Constraint(
                self.flowsheet().config.time,
                self.length_domain,
                doc="Darcy friction factor constraint for spiral-wound membranes",
            )
            def eq_friction_factor(b, t, x):
                return b.friction_factor_darcy[t, x] == 6.23 * b.N_Re[t, x] ** -0.3

        else:
            raise ConfigurationError(
                f"Unrecognized friction_factor type {friction_factor}"
            )

        ## ==========================================================================
        # Pressure change per unit length due to friction
        # -1/2*f/dh*density*velocity^2
        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            doc="pressure change per unit length due to friction",
        )
        def eq_dP_dx(b, t, x):
            return (
                b.dP_dx[t, x] * b.dh
                == -0.5
                * b.friction_factor_darcy[t, x]
                * b.properties[t, x].dens_mass_phase["Liq"]
                * b.velocity[t, x] ** 2
            )

    def _add_interface_stateblock(self, has_phase_equilibrium=None):
        """
        This method constructs the interface state blocks for the
        control volume.

        Args:
            has_phase_equilibrium: indicates whether equilibrium calculations
                                    will be required in state blocks
        Returns:
            None
        """
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = has_phase_equilibrium
        tmp_dict["defined_state"] = False  # these blocks are not inlets or outlets

        self.properties_interface = self.config.property_package.build_state_block(
            self.flowsheet().config.time,
            self.length_domain,
            doc="Material properties of feed-side membrane interface",
            **tmp_dict,
        )

    def _add_vapor_stateblock(
        self,
        property_package_vapor=None,
        property_package_args_vapor=None,
        has_phase_equilibrium=None,
    ):
        """
        This method constructs the vapor state blocks for the
        control volume.
        """
        if property_package_vapor is None or property_package_args_vapor is None:
            raise ValueError(
                "Vapor property package and its arguments must be provided."
            )

        tmp_dict = dict(**property_package_args_vapor)
        tmp_dict["has_phase_equilibrium"] = has_phase_equilibrium
        tmp_dict["defined_state"] = False  # these blocks are not inlets or outlets

        self.properties_vapor = property_package_vapor.build_state_block(
            self.flowsheet().config.time,
            self.length_domain,
            doc="Material properties of vapor phase",
            **tmp_dict,
        )

    def add_extensive_flow_to_interface(self):
        # VOLUMETRIC FLOWRATE
        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            doc="Volumetric flow at interface of inlet",
        )
        def eq_equal_flow_vol_interface(b, t, x):
            if b._skip_element(x):
                return Constraint.Skip
            return (
                b.properties_interface[t, x].flow_vol_phase["Liq"]
                == b.properties[t, x].flow_vol_phase["Liq"]
            )

    def _get_state_args(self, initialize_guess, state_args):
        """
        Arguments:
            initialize_guess : a dict of guesses for recovery, Temp_change, and deltaP.
            state_args : a dict of arguments to be passed to the property
                        package to provide an initial state for the inlet feed.
        """

        if initialize_guess is None:
            initialize_guess = {}

        if "recovery" not in initialize_guess:
            initialize_guess["recovery"] = 0.1

        if "deltaP" not in initialize_guess:
            initialize_guess["deltaP"] = -5e5

        Temp_change = initialize_guess.get("Temp_change", 40)

        # Getting source state
        if self.flow_direction == FlowDirection.forward:
            source_idx = self.length_domain.first()
        else:
            source_idx = self.length_domain.last()
        source = self.properties[self.flowsheet().config.time.first(), source_idx]

        if state_args is None:
            state_args = {}
            state_dict = source.define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        hot_outlet_args = deepcopy(state_args)
        cold_outlet_args = deepcopy(state_args)

        # Adjust flow based on recovery
        for j in self.config.property_package.component_list:
            hot_outlet_args["flow_mass_phase_comp"]["Liq", j] *= (
                1 - initialize_guess["recovery"]
            )
            cold_outlet_args["flow_mass_phase_comp"]["Liq", j] *= (
                1 + initialize_guess["recovery"]
            )

        # Adjust temperature based on Temp_change_ratio
        hot_outlet_args["temperature"] = state_args["temperature"] - Temp_change
        cold_outlet_args["temperature"] = state_args["temperature"] + Temp_change

        return {
            "inlet": state_args,
            "hot_outlet": hot_outlet_args,
            "cold_outlet": cold_outlet_args,
        }

    def _get_average_state(self, prop_in, prop_out):
        """Helper function to compute average of two property states"""
        x = 0.5
        state_args_tx = {}

        for k in prop_in:
            if isinstance(prop_in[k], dict):
                if k not in state_args_tx:
                    state_args_tx[k] = {}
                for index in prop_in[k]:
                    state_args_tx[k][index] = (1.0 - x) * prop_in[k][
                        index
                    ] + x * prop_out[k][index]
            else:
                state_args_tx[k] = (1.0 - x) * prop_in[k] + x * prop_out[k]

        return state_args_tx

    def _get_state_args_interface(self, prop_in, prop_out):
        return self._get_average_state(prop_in, prop_out)

    def _get_state_args_vapor(self, prop_in, prop_out):
        state_args = self._get_average_state(prop_in, prop_out)

        state_args_vapor = {}

        # state_args_vapor["pressure"] = state_args["pressure"]
        state_args_vapor["temperature"] = state_args["temperature"]

        state_args_vapor["flow_mass_phase_comp"] = {
            ("Liq", "H2O"): self.properties_vapor[0, 0]
            .flow_mass_phase_comp["Liq", "H2O"]
            .lb,
            ("Vap", "H2O"): state_args["flow_mass_phase_comp"][("Liq", "H2O")] / 20,
        }

        return state_args_vapor

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        if hasattr(self, "channel_height"):
            if iscale.get_scaling_factor(self.channel_height) is None:
                iscale.set_scaling_factor(self.channel_height, 1e3)

        if hasattr(self, "spacer_porosity"):
            if iscale.get_scaling_factor(self.spacer_porosity) is None:
                iscale.set_scaling_factor(self.spacer_porosity, 1)

        if hasattr(self, "dh"):
            if iscale.get_scaling_factor(self.dh) is None:
                iscale.set_scaling_factor(self.dh, 1e3)

        if hasattr(self, "K"):
            for v in self.K.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1e4)

        if hasattr(self, "N_Re"):
            for t, x in self.N_Re.keys():
                if iscale.get_scaling_factor(self.N_Re[t, x]) is None:
                    iscale.set_scaling_factor(self.N_Re[t, x], 1e-2)

        if hasattr(self, "N_Sc_comp"):
            for v in self.N_Sc_comp.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1e-2)

        if hasattr(self, "N_Sh_comp"):
            for v in self.N_Sh_comp.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1e-2)

        if hasattr(self, "N_Pr"):
            for v in self.N_Pr.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1)

        if hasattr(self, "N_Nu"):
            for v in self.N_Nu.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1)

        if hasattr(self, "h_conv"):
            for (t, x), v in self.h_conv.items():
                if iscale.get_scaling_factor(v) is None:
                    if hasattr(self, "N_Nu"):
                        sf_hconv = (
                            iscale.get_scaling_factor(
                                self.properties[t, x].therm_cond_phase["Liq"]
                            )
                            * iscale.get_scaling_factor(self.N_Nu[t, x])
                            / iscale.get_scaling_factor(self.dh)
                        )
                        iscale.set_scaling_factor(v, sf_hconv)
                    else:
                        iscale.set_scaling_factor(v, 1e3)

        if hasattr(self, "velocity"):
            for v in self.velocity.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1)

        if hasattr(self, "friction_factor_darcy"):
            for v in self.friction_factor_darcy.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1)

        if hasattr(self, "cp_modulus"):
            for v in self.cp_modulus.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1)
