
from pyomo.environ import (
    NonNegativeReals,
    NegativeReals,
    Set,
    Var,
    value,
)
from idaes.core import (declare_process_block_class,
    EnergyBalanceType,
    FlowDirection,
    )
from idaes.core.util.misc import add_object_reference
from idaes.core.util.exceptions import BalanceTypeNotSupportedError
from idaes.core.base.control_volume0d import ControlVolume0DBlockData
from watertap.core.membrane_channel_base import (
    MembraneChannelMixin, MassTransferCoefficient, PressureChangeType, ConcentrationPolarizationType) 


@declare_process_block_class("MembraneChannel0D")
class MembraneChannel0DBlockData(ControlVolume0DBlockData, MembraneChannelMixin):

    def build(self):
        super().build()

    # overwrite CV0D `add_geometry`
    def add_geometry(self):
        self._add_area_total()

    def add_state_blocks(
        self, information_flow=FlowDirection.forward, has_phase_equilibrium=None
    ):
        """
        This method constructs the state blocks for the
        control volume.

        Args:
            information_flow: a FlowDirection Enum indicating whether
                               information flows from inlet-to-outlet or
                               outlet-to-inlet
            has_phase_equilibrium: indicates whether equilibrium calculations
                                    will be required in state blocks
            package_arguments: dict-like object of arguments to be passed to
                                state blocks as construction arguments
        Returns:
            None
        """
        super().add_state_blocks(information_flow, has_phase_equilibrium)
        # quack like a 1D model
        self.length_domain = Set(ordered=True, initialize=(0.0, 1.0))
        add_object_reference(self, "difference_elements", self.length_domain)
        self.first_element = self.length_domain.first()
        add_object_reference(
            self,
            "properties",
            {
                **{
                    (t, 0.0): self.properties_in[t]
                    for t in self.flowsheet().config.time
                },
                **{
                    (t, 1.0): self.properties_out[t]
                    for t in self.flowsheet().config.time
                },
            },
        )

        self._add_interface_blocks(information_flow, has_phase_equilibrium)

    def add_mass_transfer(self):
        self._add_recovery_rejection()

        units_meta = self.config.property_package.get_metadata().get_derived_units

        # mass transfer
        def mass_transfer_phase_comp_initialize(b, t, p, j):
            return value(
                self.properties_in[t].get_material_flow_terms("Liq", j)
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
                == -self.mass_transfer_term[t, p, j]
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
                b.mixed_permeate[t].get_material_flow_terms(p, j)
                == -b.mass_transfer_term[t, p, j]
            )

        # Different expression in 1DRO
        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Permeate production",
        )
        def eq_permeate_production(b, t, p, j):
            return (
                b.mixed_permeate[t].get_material_flow_terms(p, j)
                == b.area_total * b.flux_mass_phase_comp_avg[t, p, j]
            )

        # Not in 1DRO
        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            self.config.property_package.solute_set,
            doc="Permeate mass fraction",
        )
        def eq_mass_frac_permeate(b, t, x, j):
            return (
                b.permeate_side[t, x].mass_frac_phase_comp["Liq", j]
                * sum(
                    self.flux_mass_phase_comp[t, x, "Liq", jj]
                    for jj in self.config.property_package.component_list
                )
                == self.flux_mass_phase_comp[t, x, "Liq", j]
            )


    def add_total_enthalpy_balances(self,**kwargs):
        eb = super().add_total_enthalpy_balances(**kwargs)

        if self._constructed_material_balance_type not in (EnergyBalanceType.enthalpyTotal, EnergyBalanceType.none):
            raise BalanceTypeNotSupportedError(
                "{self.name} OD membrane channels do not support {self._constructed_material_balance_type}")

        if kwargs.get("has_enthalpy_transfer", False) and eb is not None:
            # Non-existent in MembraneChannel1D 
            @self.Constraint(
                self.flowsheet().config.time, doc="Enthalpy transfer from feed to permeate"
            )
            def eq_connect_enthalpy_transfer(b, t):
                return (
                    b.mixed_permeate[t].get_enthalpy_flow_terms("Liq")
                    == -b.enthalpy_transfer[t]
                )

        return eb

    def add_volumetric_flowrate_balance(self,**kwargs):
        # not in MC1D
        @self.Constraint(
            self.flowsheet().config.time, self.length_domain, doc="Permeate flowrate"
        )
        def eq_flow_vol_permeate(b, t, x):
            return (
                b.permeate_side[t, x].flow_vol_phase["Liq"]
                == b.mixed_permeate[t].flow_vol_phase["Liq"]
            )

        super().add_volumetric_flowrate_balance(**kwargs)

    def add_expressions(self):
        super().add_expressions()

        @self.Expression(self.flowsheet().config.time, doc="Over pressure ratio")
        def over_pressure_ratio(b, t):
            return (
                b.properties_out[t].pressure_osm_phase["Liq"]
                - b.permeate_side[t, 1.0].pressure_osm_phase["Liq"]
            ) / b.properties_out[t].pressure

    def apply_transformation(self):
        pass

    def _add_calculated_pressure_change_mass_transfer_components(self):
        # NOTE: This function could be called by either
        # `_add_calculated_pressure_change` *and/or* 
        # `_add_calculated_mass_transfer_coefficient`.
        # Therefore, we add this simple gaurd against it being called twice.
        if hasattr(self, "channel_height"):
            return

        units_meta = self.config.property_package.get_metadata().get_derived_units

        # comes from ControlVolume1D in 1DRO
        self.length = Var(
            initialize=10,
            bounds=(0.1, 5e2),
            domain=NonNegativeReals,
            units=units_meta("length"),
            doc="Effective membrane length",
        )

        # not optional in 1DRO
        self.width = Var(
            initialize=1,
            bounds=(0.1, 5e2),
            domain=NonNegativeReals,
            units=units_meta("length"),
            doc="Effective feed-channel width",
        )

        self._add_area_total_equation()

        # comes from ControlVolume1D in 1DRO
        self.area = Var(
            initialize=1e-3 * 1 * 0.95,
            bounds=(0, 1e3),
            domain=NonNegativeReals,
            units=units_meta("length") ** 2,
            doc="Cross sectional area",
        )

        super()._add_calculated_pressure_change_mass_transfer_components()

    def _add_pressure_change(self, pressure_change_type=PressureChangeType.calculated):

        units_meta = self.config.property_package.get_metadata().get_derived_units

        if pressure_change_type == PressureChangeType.fixed_per_unit_length:
            # Pressure change equation when dP/dx = user-specified constant,
            self.dP_dx = Var(
                self.flowsheet().config.time,
                initialize=-5e4,
                bounds=(-2e5, -1e3),
                domain=NegativeReals,
                units=units_meta("pressure") * units_meta("length") ** -1,
                doc="pressure drop per unit length across channel",
            )
            @self.Constraint(
                self.flowsheet().config.time, doc="pressure change due to friction"
            )
            def eq_pressure_change(b, t):
                return b.pressure_change_total[t] == b.dP_dx[t] * b.length

        else:
            self.dP_dx= Var(
                self.flowsheet().config.time,
                self.length_domain,
                initialize=-5e4,
                bounds=(-2e5, -1e3),
                domain=NegativeReals,
                units=units_meta("pressure") * units_meta("length") ** -1,
                doc="Pressure drop per unit length of channel at inlet and outlet",
            )
            self._add_pressure_change_equation()
